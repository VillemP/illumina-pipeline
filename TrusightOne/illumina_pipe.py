import argparse
import os
import shlex
import subprocess
import sys
from tempfile import NamedTemporaryFile

import gene
import gene_panel
import miniseq.configvalidator
import pipeline_utility.file_utility
import pipeline_utility.sample
import truesightoneconfig
from TrusightOne.jsonhandler import JsonHandler
from gene_panel import CombinedPanels
from pipeline_utility import vcf_manipulator, adsplit, annotate_by_pos
from pipeline_utility.txttoxlsx_filtered import create_excel
from tso_excel_filters import TruesightOneFilters

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

workingDir = os.getcwd()
processes = []
tso_genes = list()

config = truesightoneconfig.loadCfg(os.path.join(os.path.dirname(os.path.normpath(__file__)), "tso.yaml"))
handler = JsonHandler(config.json_dir, config)
handler.get_all_panels(False)  # Get local data
combinedpanels = CombinedPanels(handler)
total_samples = 0
db_name_samples = config.db_name  # Updated later with current samples name
project = os.path.basename(os.path.normpath(workingDir))
assert os.path.exists(config.annotator)


def logdata(stdout):
    # TODO: Fix logging
    with open(config.logfile, "a+") as log:
        for line in iter(stdout.readline, b''):  # b'\n'-separated lines
            # log.write(time.time() + "\t" + line + "\n")
            log.write(line)
            print ("{}".format(line.rstrip()))


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
    total_samples = curr_len + i
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


def annotate(sample, testmode):
    """

    :type sample: pipeline_utility.sample as sample
    """
    print("Starting annotation for file: {0}".format(sample.name))
    global total_samples
    global db_name_samples

    # A hack to have VcfManipulator insert an annotation, but it takes in a file object if accessed from subprocess
    with NamedTemporaryFile(delete=False) as tempfile:
        sample.genes_tempfile = tempfile
        for line in sample.final_order:
            tempfile.file.write(line + "\n")
        tempfile.file.seek(0)

        if not testmode:
            sample.reduced_variants_vcf = sample.vcflocation

            args = shlex.split("perl {0} {1} {2} -buildver hg19 "
                               "-out {3} "
                               "-remove -protocol "
                               "refGene,avsnp147,1000g2015aug_all,1000g2015aug_eur,exac03,ljb26_all,clinvar_20150629 "
                               "-argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' "
                               "-operation g,f,f,f,f,f,f "
                               "-nastring . "
                               "-otherinfo "
                               "-vcfinput".format(config.annotator, sample.reduced_variants_vcf, config.annotation_db,
                                                  sample.name))

            proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
            with proc.stderr:
                logdata(proc.stderr)
            processes.append(proc)

            proc.wait()

        # annotator output file --> add custom and external gene name based
        disease_name = os.path.join(config.custom_annotation_dir, "gene.omim_disease_name.synonyms.txt")
        disease_nr = os.path.join(config.custom_annotation_dir, "gene.disease.txt")
        hpo = os.path.join(config.custom_annotation_dir, "gene.hpoterm.txt")
        panels = os.path.join(config.custom_annotation_dir, "TSO_genepanels.txt")

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
        args_5 = shlex.shlex('java -jar {0} '
                             'extractFields '
                             '- '
                             '-e . -s ";" CHROM POS avsnp147 REF ALT QUAL FILTER AC AF DP '
                             'Gene.refGene GeneReq Func.refGene GeneDetail.refGene ExonicFunc.refGene AAChange.refGene '
                             '1000g2015aug_all 1000g2015aug_eur ExAC_ALL ExAC_NFE ExAC_FIN SIFT_score SIFT_pred '
                             'Polyphen2_HVAR_score Polyphen2_HVAR_pred MutationTaster_score MutationTaster_pred '
                             'CADD_raw CADD_phred phyloP46way_placental phyloP100way_vertebrate clinvar_20150629 '
                             'Disease.name Disease.nr HPO Panel GEN[0].GT GEN[0].DP GEN[0].AD'
                             .format(config.snpsift))
        if testmode:
            total_samples = 1
        args_6 = shlex.shlex('python {0}'.format(os.path.abspath(adsplit.__file__)))
        args_7 = shlex.shlex('python {0} {1}.txt {2}'.format(os.path.abspath(annotate_by_pos.__file__),
                                                             os.path.join(config.db_directory,
                                                                          db_name_samples),
                                                             total_samples))
        # This approach retains quotation marks and complete whitespace delimited args
        slx = list([args_1, args_2, args_3, args_4, args_gene_panel, args_5, args_6, args_7])
        for arg in slx:
            arg.whitespace_split = True

        with open(sample.name + ".hg19_multianno.vcf") as annotated:
            proc_1 = subprocess.Popen([a for a in args_1], shell=False, stdin=annotated, stdout=subprocess.PIPE)
            proc_2 = subprocess.Popen([a for a in args_2], shell=False, stdin=proc_1.stdout, stdout=subprocess.PIPE)
            proc_3 = subprocess.Popen([a for a in args_3], shell=False, stdin=proc_2.stdout, stdout=subprocess.PIPE)
            proc_4 = subprocess.Popen([a for a in args_4], shell=False, stdin=proc_3.stdout, stdout=subprocess.PIPE)

            proc_gene_panel = subprocess.Popen([a for a in args_gene_panel],
                                               shell=False, stdin=proc_4.stdout,
                                               stdout=subprocess.PIPE)
            proc_5 = subprocess.Popen([a for a in args_5], shell=False, stdin=proc_gene_panel.stdout,
                                      stdout=subprocess.PIPE)
            # Proc 6 takes in a table
            proc_6 = subprocess.Popen([a for a in args_6], shell=False, stdin=proc_5.stdout, stdout=subprocess.PIPE)
            proc_7 = subprocess.Popen([a for a in args_7], shell=False, stdin=proc_6.stdout, stdout=subprocess.PIPE)
            # logdata(proc_7.stderr)
            with open(sample.name + ".annotated.table", "w+") as table:
                table.write(proc_7.communicate()[0])
                sample.table_files.append(os.path.abspath(table.name))
            for proc in (proc_1, proc_2, proc_3, proc_4, proc_gene_panel, proc_5, proc_6, proc_7):
                processes.append(proc)
        # TODO: Switch to multiprocessing?
        # jobs = multiprocessing.Pool(1)
        # annotator = multiprocessing.Process()
        proc_7.wait()
        if proc_7.returncode == 0:
            sample.annotated = True
            print("Finished annotating sample {0}".format(sample.name))
            print(sample)


def calc_coverage(sample):
    # TODO: Create a refSeq and target file for coverage analysis
    # grep $'\t'${gene}$'\.' < ${targets} >> ${covdir}${prefix}.genes.bed
    # sort -k1,1V -k2,2n -k3,3n < ${covdir}${prefix}.genes.bed > ${covdir}${prefix}.genes.sorted.bed
    sample.targetfile = None
    sample.refseq = None
    # TODO: replace with BEDtools?
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
                             '-o {5}'.format(config.toolkit, config.reference, sample.bamlocation, sample.targetfile,
                                             sample.refseq, os.path.join(path, sample.name + ".requested")))
    diagnose_args = shlex.split('java -Xmx4g -jar {0} '
                                '-T DiagnoseTargets '
                                '-R {1} '
                                '-I {2} '
                                '-L {3} '
                                '-min 20 '
                                '-o {4}'.format(config.toolkit, config.reference, sample.bamlocation, sample.targetfile,
                                                os.path.join(path, sample.name + ".diagnoseTargets")))
    table_args = shlex.split('java -Xmx4g -jar {0} '
                             '-T VariantsToTable '
                             '-R {1} '
                             '-V {2} '
                             '-F CHROM -F POS -F END -F FILTER -F IDP -F IGC '
                             '--showFiltered '.format(config.toolkit, config.reference,
                                                      os.path.join(path, sample.name + ".diagnoseTargets")))
    depth = subprocess.Popen(depth_args, shell=False)
    diagnose = subprocess.Popen(diagnose_args, shell=False)
    to_table = subprocess.Popen(table_args, shell=False, stdout=subprocess.PIPE)
    for proc in (depth, diagnose, to_table):
        processes.append(proc)
    table_name = os.path.join(path, sample.name + ".diagnoseTargets.table")
    with open(table_name, "wb+") as f:
        f.write(to_table.communicate()[0])
    diagnose.wait()
    depth.wait()

    # with to_table.stderr:
    #    logdata(to_table.stderr)

    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_summary')))
    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_summary')))


def create_excel_table(sample):
    filters = TruesightOneFilters(sample.table_files)
    if sample.annotated:
        create_excel(".".join([sample.name, "xlsx"]), filters, sample.table_files, total_samples)
    else:
        sys.stderr.write("PIPELINE ERROR: Cannot create excel file for {0} "
                         "due to incomplete annotations!\n".format(sample.name))


def getGeneOrder(samples, args):
    print("Getting the order from {0}...".format(args.panels.name))
    rows = {}
    try:
        for line in args.panels.readlines():
            (prefix, panels, genes) = line.split('\t')
            rows[prefix] = ([pan.strip().upper() for pan in panels.split(",")],
                            [g.strip().upper() for g in genes.split(",")])
        if len(rows) > len(samples):
            sys.stderr("WARNING: some samples represented in the --panels file "
                       "have no matching samples in the batch. These will be ignored.")
        for sample in samples:
            if sample.name in rows.iterkeys():
                print("{0} {1}".format(sample.name, rows[sample.name]))
                # Orderkey is the panel identifier used (e.g. ID, METABO, name or panel.id)
                for orderkey in rows[sample.name][0]:
                    results = gene_panel.match_order_to_panels(orderkey, combinedpanels, handler)
                    if len(results) > 0:
                        for panel in results:
                            sample.panels.append((panel, orderkey))
                # Orderkey is the gene name
                for orderkey in rows[sample.name][1]:
                    result = gene.find_gene(orderkey)
                    if result is not None:
                        sample.genes.append((result, "ORDERED"))
            else:
                sys.stderr.write("WARNING: sample ({0}) represented in the batch "
                                 "has no selected panels or genes in the --panels file. "
                                 "This sample will be annotated with 'ALL'.\n".format(sample.name))
    except IOError as e:
        sys.stderr.write("PIPELINE ERROR: {0}\nDoes your batch contain the panels.txt file?\n"
                         "Creating a new panel.txt file in the working directory {1}.\n".
                         format(e.message, workingDir))
        if not os.path.exists(os.path.join(workingDir, "panels.txt")):
            with open(os.path.join(workingDir, "panels.txt"), "w+") as panelsfile:
                for sample in samples:
                    print("\t".join((sample.name, "-", "-")))
                    panelsfile.write("\t".join((sample.name, "-", "-")))
        else:
            sys.stderr.write("ABORTING: panels.txt file already exists.\n")


def run_samples(args, sample_list):
    if not handler.loaded:
        handler.get_all_panels(False)  # Get local data

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
    # TODO: check for sexes
    getGeneOrder(samples, args)
    update_database(samples, args.no_replace, args.test)
    for sample in samples:
        try:
            annotate(sample, args.test)
            # calc_coverage(sample)
            create_excel_table(sample)

            if not args.keep and sample.annotated:
                os.remove(os.path.join(workingDir, sample.name + ".hg19_multianno.vcf"))
                os.remove(os.path.join(workingDir, sample.name + ".hg19_multianno.txt"))
                os.remove(os.path.join(workingDir, sample.name + ".annotated.table"))
                os.remove(os.path.join(workingDir, sample.name + ".avinput"))
                os.remove(os.path.join(workingDir, sample.name + "_coverage", sample.name + ".requested"))
                os.remove(os.path.join(workingDir, sample.name + "_coverage",
                                       sample.name + ".sample_cumulative_coverage_counts"))
                os.remove(os.path.join(workingDir, sample.name + "_coverage",
                                       sample.name + ".sample_cumulative_coverage_proportions"))
                os.remove(os.path.join(workingDir, sample.name + "_coverage",
                                       sample.name + ".sample_interval_statistics"))
                os.remove(os.path.join(workingDir, sample.name + "_coverage",
                                       sample.name + ".sample_statistics"))

            print("Finished sample {0}".format(sample.name))
            # TODO: Run conifer
        except Exception:
            pass  # continue running for other samples


def main(args):
    samples = list()

    # True if replace, False if --no_replace, passed as an argument in update_database()
    if not args.no_replace:
        sys.stderr.writelines(
            "WARNING: VCF files with file names that already exist in the database will not be updated! "
            "Old variant files still exist. (--no_replace is active)\n")

    if args.batch:
        prefixes = pipeline_utility.file_utility.write_prefixes_list(workingDir, "prefixes.list")
        run_samples(args, prefixes)

    elif args.old:
        prefixes, vcfslist, bamlist = create_configs()
        for i in range(0, len(prefixes)):
            samples.append(pipeline_utility.sample.Sample(prefixes[i], vcfslist[i], bamlist[i]))
        update_database(samples, args.no_replace, args.test)
    elif args.samples:
        print("Input is {0} samples: {1}".format(len(args.samples), "\t".join(args.samples)))
        run_samples(args, args.samples)
    elif args.json:
        yesChoice = ['yes', 'y']
        noChoice = ['no', 'n']
        input = raw_input("WARNING: You are about to download the newest gene panels."
                          " Do you wish to continue (y) or load from local data (n): ").lower()
        if input in yesChoice:
            handler.write_version_one_panels()
        elif input in noChoice:
            # print("Quitting.")
            print("Loading and writing panel tables instead...")
            combinedpanels.write_table()
    else:
        print "No valid command input."
        print "Printing cfg values."
        print config


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="TruesightOne pipeline command-line tool")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-b", "--batch",
                       help="Find all unique .vcf files and their matching .bams."
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
    parser.add_argument("-p", "--panels", type=argparse.FileType('r'))
    parser.add_argument("-r", "--no_replace",
                        help="Will skip copying the VCF file to the VCF directory specified in the config file"
                             " (don't overwrite VCF in the DB)."
                             "Useful in a test scenario when running large batches repeatedly or if the VCF is "
                             "not supposed to be inserted into the DB.",
                        action="store_false", default=True)
    parser.add_argument("-t", "--test", help="Skips the variant DB recompilation (CombineVariants) step, "
                                             "expects the DB to exist.",
                        action="store_true", default=False)
    parser.add_argument("-k", "--keep", help="Keep intermediary annotation and coverage files.",
                        action="store_true", default=False)

    args = parser.parse_args()
    if args.batch and args.panels is None:
        parser.error("--batch requires --panels\nCheck {0} --help panels for more info".format(__file__))
    if args.samples is not None and args.panels is None:
        parser.error("--samples requires --panels\nCheck {0} --help panels for more info".format(__file__))

    main(args)

    for p in processes:
        try:
            p.terminate()
        except OSError:
            pass  # process is already dead
