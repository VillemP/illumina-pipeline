import argparse
import multiprocessing
import os
import shlex
import subprocess
import sys
import time
import traceback
from Queue import Empty
from tempfile import NamedTemporaryFile

import TrusightOne.gene
import TrusightOne.gene_panel
import TrusightOne.truesightoneconfig
import miniseq.configvalidator
import pipeline_utility.file_utility
import pipeline_utility.sample
from TrusightOne.gene_panel import CombinedPanels
from TrusightOne.jsonhandler import JsonHandler
from TrusightOne.tso_excel_filters import TruesightOneFilters
from TrusightOne.tso_excel_filters import TruesightOneFormats
from TrusightOne.tso_excel_filters import TruesightOnePostprocess
from pipeline_utility import vcf_manipulator, adsplit, annotate_by_pos, db_update, converthgnc
from pipeline_utility.converthgnc import HgncConverter
from pipeline_utility.txttoxlsx_filtered import create_excel

try:
    from cStringIO import StringIO
except ImportError:
    print("cStringIO is not installed, I/O is slower.")
    from StringIO import StringIO

workingDir = os.getcwd()
processes = []
tso_genes = list()

config = TrusightOne.truesightoneconfig.loadCfg(os.path.join(os.path.dirname(os.path.normpath(__file__)), "tso.yaml"))
handler = JsonHandler(config.json_dir, config)
combinedpanels = None
# Number of samples in the DB
# By default, the base db number is 0 but can be set to an increased number, if
# a combined VCF file is used (e.g. from a previous batch)
total_samples = config.db_base_n
db_name_samples = config.db_name  # Updated later with current samples name
project = os.path.basename(os.path.normpath(workingDir))
assert os.path.exists(config.annotator)


def async_log(stdout, queue, pid=None):
    if stdout is not None and not stdout.closed:
        for line in iter(stdout.readline, b''):  # b'\n'-separated lines
            queue.put((pid, line))


def logdata(stdout, pid=None):
    print("LOG: Connecting to subpipe (PID) {0}: {1}".format(pid, stdout))
    with open(config.logfile, "a+") as log:
        if stdout is not None:
            for line in iter(stdout.readline, b''):  # b'\n'-separated lines
                log.write(str(time.time()) + "\t" + line)
                print("{}".format(line.rstrip()))
        else:
            sys.stderr.write("\nPIPELINE ERROR: Log connected to empty subprocess!\n")


def sortlog(logfile):
    with open(logfile, "w+") as log:
        with open(logfile + ".sorted.txt", "wb+") as out:
            lines = log.readlines()
            lines.sort()
            for line in lines:
                text = line.split("\t")
                # Rejoin the text with any tabs that were split out
                # Removes the timestamp.
                out.write("\t".join(text[1:]) + "\n")


# deprecated, used with --old
def create_configs():
    prefixes = pipeline_utility.file_utility.write_prefixes_list(workingDir, "prefixes.list")
    vcflist = pipeline_utility.file_utility.write_vcfs_list(workingDir, "vcfs.list")
    bamlist = pipeline_utility.file_utility.write_bams_list(workingDir, "bams.list")

    if miniseq.configvalidator.validate_config("bams.list", "vcfs.list", "prefixes.list"):
        return prefixes, vcflist, bamlist
    else:
        raise IOError("Database updating failed due to incorrect files.")


# deprecated, used with --old
def create_arguments_file(dbnr):
    with open("arguments.txt", "wb+") as f:
        f.write("wd:\n")
        f.write(workingDir + "\n")
        f.write("dbname:\n")
        f.write(dbnr + "\n")
        f.write("project:\n")
        f.write(project)


def update_vcf_list(vcfs_list):
    global db_name_samples
    global total_samples

    # If the file doesn't exist (on first run, the length is automatically 0)
    if os.path.exists(config.db_vcf_dir):
        curr_len = pipeline_utility.file_utility.file_len(config.db_vcf_dir)
    else:
        curr_len = 0

    with open(config.db_vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()
        # db_vcfs.seek(0)

        i = 0
        skipped = 0
        for vcf in vcfs_list:
            name = os.path.basename(vcf).rsplit('.')[0]
            if not any(name in x.rstrip() for x in data):
                line = "V:{0} {1}".format(name, os.path.join(config.vcf_storage_location,
                                                             os.path.basename(vcf)))
                i += 1
                db_vcfs.write(line + "\n")
            else:
                skipped += 1
    total_samples += curr_len + i
    db_name_samples = db_name_samples + str(total_samples)
    create_arguments_file(str(total_samples))
    print ("Updated {0} with {1} unique samples. "
           "Skipped {2} preexisting samples. New database name is {3}.").format(
        os.path.join(config.db_directory, config.db_vcf_list_name), i, skipped, db_name_samples)
    return i, skipped


def combine_variants(out, vcflist=config.db_vcf_dir):
    # Combining variant files into a single reference to be used for statistical purposes
    args = shlex.split('java -Xmx10g -jar {0} '
                       '-T CombineVariants '
                       '-R {1} '
                       '-V {2} '
                       '-o {3} '
                       '-log {4} '
                       '--genotypemergeoption UNIQUIFY'.format(config.toolkit, config.reference, vcflist,
                                                               out, config.logfile))

    proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
    processes.append(proc)

    logdata(proc.stderr)
    proc.wait()

    args = shlex.split('java -Xmx10g -jar {0} '
                       '-T VariantsToTable '
                       '-R {1} '
                       '-V {2}{3}.vcf '
                       '-F CHROM -F POS -F REF -F ALT -F AC -F HET -F HOM-VAR '
                       '--splitMultiAllelic --showFiltered '
                       '-o {4}{5}.txt '.format(config.toolkit, config.reference, config.db_directory, db_name_samples,
                                               config.db_directory,
                                               db_name_samples, config.logfile))

    proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
    processes.append(proc)
    logdata(proc.stderr)
    proc.wait()
    # print proc.returncode
    # proc.communicate()


def update_database(samples, replace, testmode):
    if not testmode:
        assert len(samples) > 0, "List of samples to be updated into the database cannot be empty!"
        vcfslist = list()
        for sample in samples:
            vcfslist.append(sample.vcflocation)
        copied, skipped = pipeline_utility.file_utility.copy_vcf(vcfslist, config.vcf_storage_location, replace)
        updated, skip = update_vcf_list(vcfslist)
        if copied + updated > 0:
            # If there were any variant files updated (copied) or
            combine_variants(os.path.join(config.db_directory, db_name_samples + ".vcf"))
    pass


def genderCheck(args, samples):
    out = project + ".merged.vcf"
    vcflist = ""
    for sample in samples:
        vcflist = vcflist + "--variant:{0} {1}\\".format(sample.name, sample.vcflocation)
        # print vcflist
    if not args.testmode:
        args = shlex.split('java -Xmx10g -jar {0} '
                           '-T CombineVariants '
                           '-R {1} \\'
                           '{2} \\'
                           '-o {3} \\'
                           '-log {4} \\'
                           '--genotypemergeoption UNIQUIFY'.format(config.toolkit, config.reference, vcflist,
                                                                   out, config.logfile))

        proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        processes.append(proc)
        logdata(proc.stderr)
        logdata(proc.stdout)
        proc.wait()

        het_args = shlex.split("vcftools --vcf {0} "
                               "--remove-filtered-all --chr chrX "
                               "--from-bp 2699520 --to-bp 154931043 "
                               "--het --out {1}".format(out, project))
        proc = subprocess.Popen(het_args, shell=False, stderr=subprocess.PIPE)
        processes.append(proc)
        logdata(proc.stderr)
        proc.wait()


def convert_to_hgnc(converter, infile, outfile, args):
    print("Converting VCF gene names to HGNC names. Unknown gene symbols will be preserved.")
    converter.convert(infile, outfile, True)


def annotate(sample, args, converter):
    """

    :param converter:
    :type sample: pipeline_utility.sample as sample
    """
    print("Starting annotation for file: {0}".format(sample.name))
    global total_samples
    global db_name_samples

    # A hack to have VcfManipulator insert an annotation, but it takes in a file object if accessed from subprocess
    with NamedTemporaryFile(delete=False, prefix=sample.name + ".genes.") as tempfile:
        sample.genes_tempfile = tempfile
        for line in sample.final_order:
            tempfile.file.write(line + "\n")
        tempfile.file.seek(0)

        if not args.testmode:
            sample.reduced_variants_vcf = sample.vcflocation

            annotation = shlex.split('perl {0} "{1}" "{2}" -buildver hg19 '
                                     '-out "{3}" '
                                     "-remove -protocol "
                                     "refGene,avsnp147,1000g2015aug_all,1000g2015aug_eur,exac03,ljb26_all,clinvar_20170130 "
                                     "-argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' "
                                     "-operation g,f,f,f,f,f,f "
                                     "-nastring . "
                                     "-otherinfo "
                                     "-vcfinput".format(config.annotator,
                                                        sample.reduced_variants_vcf.replace(" ", "\\ "),
                                                        config.annotation_db,
                                                        sample.name))

            print(args)
            proc = subprocess.Popen(annotation, shell=False, stderr=subprocess.PIPE)
            with proc.stderr:
                # logdata(proc.stderr)
                pass
            processes.append(proc)

            proc.wait()

        # annotator output file --> add custom and external gene name based annotations
        # (if VARIANT.refGene == SYMBOL, insert custom_list[SYMBOL]=...) (vcfmanipulator)
        disease_name = os.path.join(config.custom_annotation_dir, "gene.omim_disease_name.synonyms.txt")
        disease_nr = os.path.join(config.custom_annotation_dir, "gene.disease.txt")
        hpo = os.path.join(config.custom_annotation_dir, "gene.hpoterm.txt")
        panels = os.path.join(config.custom_annotation_dir, "TSO_genepanels.txt")

        # Add extra annotations to the VCF
        args_0 = shlex.shlex("{0} {1} {2} {3} {4}".format("python " + os.path.abspath(converthgnc.__file__), "--hgnc "
                                                          + config.hgncPath, "--input -", "--output -", "--progress"))
        args_1 = shlex.shlex(
            "{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), disease_name, "Disease.name"))
        args_2 = shlex.shlex(
            "{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), disease_nr, "Disease.nr"))
        args_3 = shlex.shlex("{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), hpo, "HPO"))
        args_4 = shlex.shlex(
            "{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), panels, "Panel"))
        args_gene_panel = shlex.shlex("{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__),
                                                           sample.genes_tempfile.name,
                                                           "GeneReq"))
        # Extract the fields from the annotated VCF
        args_5 = shlex.shlex('java -Xmx4g -jar {0} '
                             'extractFields '
                             '- '
                             '-e . -s ";" CHROM POS avsnp147 REF ALT QUAL FILTER AC AF DP '
                             'Gene.refGene GeneReq Func.refGene GeneDetail.refGene ExonicFunc.refGene AAChange.refGene '
                             '1000g2015aug_all 1000g2015aug_eur ExAC_ALL ExAC_NFE ExAC_FIN SIFT_score SIFT_pred '
                             'Polyphen2_HVAR_score Polyphen2_HVAR_pred MutationTaster_score MutationTaster_pred '
                             'CADD_raw CADD_phred phyloP46way_placental phyloP100way_vertebrate CLINSIG CLNDBN CLNACC '
                             'CLNDSDB CLNDSDBID '
                             'Disease.name Disease.nr HPO Panel GEN[0].GT GEN[0].DP GEN[0].AD'
                             .format(config.snpsift))
        if args.testmode:
            total_samples = 3
            db_name_samples = db_name_samples + str(total_samples)
        # Splits the last column (allele 1 depth, allele 2 depth into 3 columns:
        # (reference reads, variant reads, percentage),
        args_6 = shlex.shlex('python {0}'.format(os.path.abspath(adsplit.__file__)))
        # Adds the custom local DB frequencies to the n'th column (default is the 23. column)
        args_7 = shlex.shlex('python {0} {1}.txt {2}'.format(os.path.abspath(annotate_by_pos.__file__),
                                                             os.path.join(config.db_directory,
                                                                          db_name_samples),
                                                             total_samples))
        # This approach retains quotation marks and complete whitespace delimited args
        slx = list([args_0, args_1, args_2, args_3, args_4, args_gene_panel, args_5, args_6, args_7])
        for arg in slx:
            arg.whitespace_split = True

        # try:
        with open(sample.name + ".hg19_multianno.vcf") as annotated:
            proc_0 = subprocess.Popen([a for a in args_0], shell=False, stdin=annotated, stdout=subprocess.PIPE)
            proc_1 = subprocess.Popen([a for a in args_1], shell=False, stdin=proc_0.stdout, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            proc_2 = subprocess.Popen([a for a in args_2], shell=False, stdin=proc_1.stdout, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            proc_3 = subprocess.Popen([a for a in args_3], shell=False, stdin=proc_2.stdout, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            proc_4 = subprocess.Popen([a for a in args_4], shell=False, stdin=proc_3.stdout, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)

            proc_gene_panel = subprocess.Popen([a for a in args_gene_panel],
                                               shell=False, stdin=proc_4.stdout,
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc_6 = subprocess.Popen([a for a in args_5], shell=False, stdin=proc_gene_panel.stdout,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # Proc 6 takes in a table
            proc_7 = subprocess.Popen([a for a in args_6], shell=False, stdin=proc_6.stdout, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            proc_8 = subprocess.Popen([a for a in args_7], shell=False, stdin=proc_7.stdout, stdout=subprocess.PIPE)

            with open(sample.name + ".annotated.table", "w+") as table:
                table.write(proc_8.communicate()[0])
                sample.table_files.append(os.path.abspath(table.name))
                if proc_8.returncode == 0:
                    sample.annotated = True
                    print("Finished annotating sample {0}".format(sample.name))
                    print(sample)
            jobs = multiprocessing.Pool(8)
            q = multiprocessing.Queue()
            # t = Thread(target=async_log, args=(proc_8.stdout, q))
            # t.daemon = True
            # t.start()

            try:
                pid, line = q.get_nowait()  # or q.get(timeout=.1)
                pass
            except Empty:
                pass
            else:  # got line
                print(pid, line)
            for proc in (proc_1, proc_2, proc_3, proc_4, proc_gene_panel, proc_6, proc_7, proc_8):
                processes.append(proc)
                # jobs.apply(logdata, proc.stderr,{'pid':proc.pid})


def calc_coverage(sample):
    # grep $'\t'${gene}$'\.' < ${targets} >> ${covdir}${prefix}.genes.bed
    # sort -k1,1V -k2,2n -k3,3n < ${covdir}${prefix}.genes.bed > ${covdir}${prefix}.genes.sorted.bed
    # TODO: replace with BEDtools?
    with NamedTemporaryFile(delete=False, prefix=sample.name + ".target.", suffix=".bed") as temptarget:
        pipeline_utility.file_utility.write_targetfile(sample.order_list, targetfile=config.targetfile,
                                                       out=temptarget.file)
    sample.temptargetfile = temptarget.name

    with NamedTemporaryFile(delete=False, prefix=sample.name + ".reference.", suffix=".refseq") as temprefseq:
        pipeline_utility.file_utility.write_refseq(sample.order_list, config.refseq, out=temprefseq.file)

    sample.temprefseq = temprefseq.name

    path = os.path.join(workingDir, sample.name + "_coverage")
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    depth_args = shlex.split('java -Xmx4g -jar {0} '
                             '-T DepthOfCoverage '
                             '-R {1} '
                             '-I {2} '
                             '-L {3} '
                             '-geneList {4} '
                             '-ct 1 -ct 10 -ct 20 -ct 50 '
                             '-o {5}'.format(config.toolkit, config.reference, sample.bamlocation,
                                             sample.temptargetfile,
                                             sample.temprefseq, os.path.join(path, sample.name + ".requested")))
    diagnose_args = shlex.split('java -Xmx4g -jar {0} '
                                '-T DiagnoseTargets '
                                '-R {1} '
                                '-I {2} '
                                '-L {3} '
                                '-min 20 '
                                '-o {4}'.format(config.toolkit, config.reference, sample.bamlocation,
                                                sample.temptargetfile,
                                                os.path.join(path, sample.name + ".diagnoseTargets")))
    table_args = shlex.split('java -Xmx4g -jar {0} '
                             '-T VariantsToTable '
                             '-R {1} '
                             '-V {2} '
                             '-F CHROM -F POS -F END -F FILTER -F IDP -F IGC '
                             '--showFiltered '.format(config.toolkit, config.reference,
                                                      os.path.join(path, sample.name + ".diagnoseTargets")))
    depth = subprocess.Popen(depth_args, shell=False)
    # depth  = None
    # diagnose = subprocess.Popen(diagnose_args, shell=False)
    diagnose = None
    for proc in (depth, diagnose):
        processes.append(proc)
    # table_name = os.path.join(path, sample.name + ".diagnoseTargets.table")
    # with open(table_name, "wb+") as f:
    #    f.write(to_table.communicate()[0])
    # diagnose.wait()
    depth.wait()
    # to_table = subprocess.Popen(table_args, shell=False, stdout=subprocess.PIPE)

    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_summary')))
    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_gene_summary')))


def create_excel_table(sample):
    # Creates an excel filterset that will be automatically applied to every output file
    """
    This is the function applied to every sample to create the output excel file. Also creates the excel filterset and
    applies the postprocessing rules (e.g. change a column into hyperlinks) and
    formats (make hyperlinks look like hyperlinks, add bold or italics, etc). Uses the txttoxlsx_filtered.py script but
    applies TrusightOne specific custom filters, postprocess rules and formats.
    :param sample: Sample object that contains the output files.
    """
    filters = TruesightOneFilters(total_samples, sample.table_files)
    postprocess = TruesightOnePostprocess(sample.table_files)
    formats = TruesightOneFormats(sample.table_files)

    if args.annotate and sample.annotated:
        create_excel(".".join([sample.name, "xlsx"]), sample.table_files, filters, postprocess, formats)
        sample.finished = True
    else:
        sys.stderr.write("PIPELINE ERROR: Cannot create excel file for {0} "
                         "due to incomplete annotations!\n".format(sample.name))


def getGeneOrder(samples, args):
    """
    This function opens the file from --panels (-p) and tries to compile a list of genes for every sample.
    Each variant corresponds to a gene (added by the annotator) regardless if it is an exonic, intronic,
    splicing or some other variant. The variants that match the gene within the GeneOrder
    have a GeneReq column annotated corresponding to the order. This enables a final filtering of variants that match
    the ordered genes (i.e. if only one Panel or one Gene is ordered, the final table will automatically be filtered to
    show variants matching genes from this order)
    :param samples: Samples to match to the order (indexed by sample_name)
    :param args: args the annotation command was started with
    """
    print("Getting orders from {0}...".format(args.panels.name))
    rows = {}
    try:
        for line in args.panels.readlines():
            # Don't read empty lines or comments
            if len(line) > 0 and not line.startswith("#"):
                (prefix, panels, genes) = line.split('\t')
                rows[prefix] = ([pan.strip().upper() for pan in panels.split(",")],
                                [g.strip().upper() for g in genes.split(",")])
        if len(rows) > len(samples):
            sys.stderr.write("WARNING: some samples represented in the --panels file "
                             "have no matching samples in the batch. "
                             "The following will be ignored: {0}\n"
                             .format([pref for pref in rows.iterkeys() if pref not in [s.name for s in samples]]))
            pass
        print("-" * 90)
        for sample in samples:
            if sample.name in rows.iterkeys():
                print("{0} {1}".format(sample.name, rows[sample.name]))
                # Orderkey is the panel identifier used (e.g. ID, METABO, name or panel.id)
                for orderkey in rows[sample.name][0]:
                    results = TrusightOne.gene_panel.match_order_to_panels(orderkey, combinedpanels, handler)
                    if len(results) > 0:
                        for panel in results:
                            sample.panels.append((panel, orderkey))
                # Orderkey is the gene name
                for orderkey in rows[sample.name][1]:
                    result = TrusightOne.gene.find_gene(orderkey)  # The gene name
                    if result is not None:
                        sample.genes.append((result, "ORDERED"))
            else:
                sys.stderr.write("WARNING: sample ({0}) represented in the batch "
                                 "has no selected panels or genes in the --panels file. "
                                 "This sample will be annotated with 'ALL'.\n".format(sample.name))
        print("-" * 90)
    except IOError as e:
        sys.stderr.write("PIPELINE ERROR: {0}\nDoes your batch contain the panels.txt file?\n"
                         "Creating a new panel.txt file in the working directory {1}.\n".
                         format(e.message, workingDir))
        if not os.path.exists(os.path.join(workingDir, "panels.txt")):
            with open(os.path.join(workingDir, "panels.txt"), "w+") as panelsfile:
                panelsfile.write("#SAMPLE\tGENES (comma-seperated)\tPANELS (comma-seperated)")
                for sample in samples:
                    print("\t".join((sample.name, "-", "-")))
                    panelsfile.write("\t".join((sample.name, "-", "-")))
        else:
            sys.stderr.write("ABORTING: panels.txt file already exists but might be corrupt.\n")
            sys.exit(1)


def run_samples(args, sample_list):
    if not handler.loaded:
        handler.get_all_panels(False)  # Get local data

    finished_samples = []
    unfinished_samples = []
    samples = list()

    for sample in sample_list:
        location = pipeline_utility.file_utility.find_file(workingDir, sample + ".vcf")[1]
        bamlocation = pipeline_utility.file_utility.find_file(workingDir, sample + ".bam")[1]
        try:
            samp = pipeline_utility.sample.Sample(sample, location, bamlocation)
            samples.append(samp)
        except IOError as e:
            sys.stderr.write("PIPELINE ERROR: {0}\nIs your sample name correct? "
                             "Do the BAM and VCF files have the same filename? "
                             "Do not include file endings!\n"
                             "After fixing the errors, rerun {1} --s {2}\n".
                             format(e.message, __file__, sample))
            pass
    if args.gendercheck:
        genderCheck(args, samples)
    getGeneOrder(samples, args)
    if args.annotate:
        converter = HgncConverter(config.hgncPath)
    # DB is already updated
    if not args.old:
        update_database(samples, args.no_replace, args.testmode)
    for sample in samples:
        try:
            if args.annotate:
                annotate(sample, args, converter)
                sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.vcf"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.txt"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".annotated.table"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".avinput"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".converted.vcf"))
                sample.trash.append(sample.genes_tempfile.name)
            if args.coverage:
                calc_coverage(sample)
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage", sample.name + ".requested"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".sample_cumulative_coverage_counts"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".sample_cumulative_coverage_proportions"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".sample_interval_statistics"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".sample_statistics"))
                sample.trash.append(sample.temprefseq)
                sample.trash.append(sample.temptargetfile)
            create_excel_table(sample)

            # Delete the intermediary files, unnecessary files and files that have already been inserted
            # into the output file.
            if not args.keep:
                for trashfile in sample.trash:
                    try:
                        os.remove(trashfile)
                    # The file might've been already deleted or did not exist in the first place.
                    except(OSError) as oserror:
                        sys.stderr.write("Couldn't remove file: {0}\n{1}\n".format(trashfile, oserror.message))

            if sample.finished:
                print("Finished sample {0}".format(sample.name))
                finished_samples.append(sample)
            else:
                print("Could not finish sample {0}".format(sample))
                unfinished_samples.append(sample)

        except Exception as error:
            sys.stderr.write("PIPELINE ERROR: {0}\nTrace: ".format(sample))
            traceback.print_exc(file=sys.stderr)
            raise error
            # TODO: Run conifer for the whole batch
    print("-" * 40)
    print("Annotated {0} samples of {1} ordered/found.".format(len(finished_samples), len(samples)))
    print("-" * 40)
    if len(unfinished_samples) > 0:
        for sample in unfinished_samples:
            print("{0} is unfinished. Check for errors.".format(sample))
        print("-" * 40)
        print("You can rerun with --samples {0}".format(" ".join([s.name for s in unfinished_samples])))
    sortlog(config.logfile)


def query_symbol(symbol_list):
    for inp in symbol_list:
        synonyms = TrusightOne.gene.find_synonyms(inp)
        print("'{0}' is HGNC declared symbol: {1}".format(inp, TrusightOne.gene.is_hgnc(inp)))
        if len(synonyms) > 0:
            print("HGNC synonyms for {0}: {1}".format(inp, synonyms))

        result = TrusightOne.gene.find_gene(inp)
        if result is not None:
            print("Gene is covered on TSO as {1}.\n"
                  "HGNC symbol: {0}".format(result, result._name))
        else:
            # Is it a panel?
            result = TrusightOne.gene_panel.match_order_to_panels(inp, combinedpanels, handler)
            if type(result) is list:
                if len(result) > 0:
                    print("Found a panel match for input '{0}'!".format(inp))
                    total_genes = {}
                    for match in result:
                        for g in match.tso_genes:
                            total_genes[g.name] = g
                    print("Genes in combined panel: {0}".format(len(total_genes)))
                    for match in result:
                        print(match)
                        for gene in match.tso_genes:
                            # print("{0}\t{1}".format(gene, gene.coverage))
                            pass

                else:
                    print ("No match found for '{0}'".format(inp))
            else:
                print("Found a panel match for input '{0}'!".format(inp))
                print(result)
                for gene in result.tso_genes:
                    print("{0}\t{1}".format(gene, gene.coverage))
        print("-" * 40)


def main(args):
    print("Running TSO pipeline tool with {0}".format(args))
    global combinedpanels
    samples = list()
    prefixes = list()
    yesChoice = ['yes', 'y']
    noChoice = ['no', 'n']

    # Will skip loading the panels from external (boolean) or local data
    if not args.update and not args.legacy:
        handler.get_all_panels(False)  # Get local data
        combinedpanels = CombinedPanels(handler)

    # Create list files for running legacy scripts for CONIFER
    if args.legacy:
        prefixes, vcfslist, bamlist = create_configs()

    # True if replace, False if --no_replace, passed as an argument in update_database()
    if not args.no_replace:
        sys.stderr.write(
            "WARNING: VCF files with file names that already exist in the database will not be updated! "
            "Old variant files still exist. (--no_replace is active)\n")

    if args.batch:
        # This means args.legacy is False and we already have our list of prefixes
        if len(prefixes) == 0 and not args.legacy:
            prefixes = pipeline_utility.file_utility.write_prefixes_list(workingDir, "prefixes.list")
        run_samples(args, prefixes)

    elif args.old:
        prefixes, vcfslist, bamlist = create_configs()
        for i in range(0, len(prefixes)):
            samples.append(pipeline_utility.sample.Sample(prefixes[i], vcfslist[i], bamlist[i]))
        update_database(samples, args.no_replace, args.testmode)
    elif args.samples:
        print("Input is {0} samples: {1}".format(len(args.samples), "\t".join(args.samples)))
        run_samples(args, args.samples)
    elif args.json:
        while True:
            try:

                input = raw_input("WARNING: You are about to download the newest gene panels."
                                  " Do you wish to continue {0} or load from local data and"
                                  " write human-readable gene tables {1} [Ctrl-c to quit]: \n"
                                  .format(yesChoice, noChoice)).lower()
                if input in yesChoice:
                    handler.write_version_one_panels(True)
                elif input in noChoice:
                    print("Loading and writing panel tables instead...")
                    combinedpanels.write_table()
            except KeyboardInterrupt:
                print("Exiting...")
                break

    elif args.update:
        input = raw_input("WARNING: You are about to download the OMIM and HPO terms."
                          " The updated files will be downloaded to {0}.\n"
                          "Do you wish to continue {1} or not {2}: "
                          .format(config.custom_annotation_dir, yesChoice, noChoice)).lower()
        if input in yesChoice:
            db_update.update_all(config.custom_annotation_dir, args.testmode)
        elif input in noChoice:
            print("Quitting...")

    elif args.search is not None:
        try:
            print("TSO PIPELINE GENE TOOL")
            print("-" * 40)
            print("Input the gene name or panel name/panel tag (MITO, etc)/panel id code (PanelApp) to find "
                  "if it is covered TSO or return the genes from the panel that are covered on TSO.")
            if len(args.search) > 0:
                query = args.search
                query_symbol([query])

            while True:
                query = raw_input("Input (case-sensitive) [Ctrl+C to exit]: \n")
                query = query.split(",")
                query_symbol(query)
        except KeyboardInterrupt:
            print("Exiting...")
    else:
        print("No valid command input. Use --help or --manual to view commands.")
        print("Printing cfg values.")
        for line in config:
            print(line)


class Logger(object):
    def __init__(self, output):
        self.terminal = sys.stdout
        self.log = open(output, "a+")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.log.flush()
        self.terminal.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="TruesightOne pipeline command-line tool")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-b", "--batch",
                       help="Find all unique .vcf files and their matching .bams. in the current directory tree. "
                            "Program will only run if each vcf has a matching .bam file.",
                       action="store_true")
    group.add_argument("-s", "--samples",
                       help="Followed by a list of unique sample identifiers e.g. "
                            "-s E0000001 E000002 E0000003 which are to be run through the pipeline. "
                            "Will only run if samples have a matching vcf and bam file.",
                       action="store", nargs='+', type=str)
    group.add_argument("-o", "--old", help="Creates text based config files for shell scripts. "
                                           "Automatically updates the variant database.",
                       action="store_true")
    group.add_argument("-j", "--json", action="store_true", default=False,
                       help="Updates the existing JSON panel database. "
                            "Can also write the tab delimited table of all panels "
                            "(version 1.0 and upper) and/or table of combined panels.")
    group.add_argument("-u", "--update", action="store_true", default=False,
                       help="Updates the custom annotations for HPO, OMIM terms.")
    group.add_argument("-m", "--search", action="store", nargs="?", type=str, const="",
                       help="The gene tool is a tool to help search for genes and panels covered by Truesight One "
                            "from custom panels, PanelApp panels and HGNC gene symbols.")
    parser.add_argument("-p", "--panels", type=argparse.FileType('r'),
                        help="The input file that determines the genes or panels that will be marked in the final excel"
                             " file (column GeneReq) with an extra annotation "
                             "(ORDERED for genes and PANEL_KEY/PANEL_ID/PANEL_NAME) for "
                             "each gene represented in the panel). "
                             "Genes can have numerous annotations depending on the complexity of the order. An example "
                             "--panels panels.txt file has one or more samples (matching the --samples or --batch arg"
                             " with the structure\n\n"
                             "\tSAMPLE_ID_1<tab>PANEL_NAME_1, PANEL_2, PANEL_ID, etc<tab>GENE1, DMD, FBN1 or "
                             "-<newline>\n"
                             "SAMPLE_ID_2<tab>PANEL_2<tab>-\\n\n"
                             "SAMPLE_ID_3<tab>-<tab>GENE2, DMD, FBN1\\n\n\n")
    parser.add_argument("-r", "--no_replace",
                        help="Will skip copying the VCF file to the VCF directory specified in the config file"
                             " (don't overwrite VCF in the DB)."
                             "Useful in a test scenario when running large batches repeatedly or if the VCF is "
                             "not supposed to be inserted into the DB.",
                        action="store_false", default=True)
    parser.add_argument("-a", "--annotate", default=False, action="store_true",
                        help="The output will contain a table with the annotated variants. "
                             "Run this explicitly for annotations")
    parser.add_argument("-c", "--coverage", help="Calculates the coverage tables for the sample's order and adds it to "
                                                 "the output excel file new pages.",
                        action="store_true", default=False)
    parser.add_argument("-t", "--testmode", help="Skips the variant DB recompilation (CombineVariants) step, "
                                                 "expects the DB to exist.",
                        action="store_true", default=False)
    parser.add_argument("-k", "--keep", help="Keep intermediary annotation and coverage files.",
                        action="store_true", default=False)
    parser.add_argument("-l", "--legacy", action="store_true", default=False, help="Create the list files required "
                                                                                   "for running legacy bash scripts.")
    parser.add_argument("-g", "--gendercheck", default=False, action="store_true",
                        help="Calculate the heterogenity for each sample and output it as a table.")

    args = parser.parse_args()
    if args.batch and args.panels is None:
        parser.error("--batch requires --panels\nCheck {0} --help panels for more info".format(__file__))
    if args.samples is not None and args.panels is None:
        parser.error("--samples requires --panels\nCheck {0} --help panels for more info".format(__file__))

    main(args)

    for p in processes:
        try:
            if p is not None:
                sys.stderr.write("Killing subprocesses {0}\n".format(p.pid))
                p.terminate()
        except OSError:
            pass  # process is already dead
        finally:
            sys.stderr.flush()
