import argparse
import atexit
import datetime
import os
import shlex
import shutil
import subprocess
import sys
import thread
import traceback
from multiprocessing import Pool

import pipeline_utility.file_utility
import pipeline_utility.sample
# TODO: Run a set of commands from STDOUT -> STDIN
from miniseq.myeloid_excel_filters import MyeloidFilters, MyeloidFormats, MyeloidPostprocess
from miniseq.myeloidconfig import MyeloidConfig
from pipeline_utility import vcf_manipulator, adsplit, annotate_by_pos, local_db_tool
from pipeline_utility.txttoxlsx_filtered import create_excel

try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

processes = []
submitted = 0
cluster = None
# Location of the cfg file is in the same folder as this script
cfg_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'tsm.yaml')
config = MyeloidConfig(cfg_path)
config = config.load()
total_samples = 0
db_name_samples = config.db_name  # Updated later with current samples name
workingDir = os.getcwd()
project = os.path.basename(os.path.normpath(workingDir))

assert os.path.exists(config.annotator), "Did you fill out your config fields correctly?"


class JobStats():
    def __init__(self):
        self.lock = thread.allocate_lock()
        self.finished_samples = []
        self.unfinished_samples = []

    def change(self, sample, value):
        self.lock.acquire()
        if value is True:
            self.finished_samples.append(sample)
        else:
            self.unfinished_samples.append(sample)
        self.lock.release()


def logdata(stdout):
    # TODO: Fix logging
    with open(config.logfile, "a+") as log:
        for line in iter(stdout.readline, b''):  # b'\n'-separated lines
            # log.write(time.time() + "\t" + line + "\n")
            log.write(line)
            print ("{}".format(line.rstrip()))


# This function is used as a signature that is called later to update the global vars in this module.
def update_sample_stats(dns, ts):
    global db_name_samples
    global total_samples

    db_name_samples = dns
    total_samples = ts


def combine_variants():
    # Combining variant files into a single reference to be used for statistical purposes
    # ip = interval padding, required for inclusion of splicing variants (bed targets are too precise for exons)
    args = shlex.split('java -Xmx10g -jar {0} '
                       '-T CombineVariants '
                       '-R {1} '
                       '-V {2} '
                       '-L {3} '
                       '-o {4} '
                       '-ip {5} '
                       '-log {6} '
                       '--genotypemergeoption UNIQUIFY'.format(config.toolkit, config.reference, config.db_vcf_dir,
                                                               config.targetfile,
                                                               os.path.join(config.db_directory,
                                                                            db_name_samples + ".vcf"),
                                                               config.padding, config.logfile))

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


def annotate(sample, args, config):
    """
    Takes in a VCF file, reduces it to match coding areas only, annotate it using ANNOVAR
    Also creates a VCF file using ITDSeek if it exists
    Uses an ugly subprocess approach.
    :param sample:
    :param args:
    :param config:
    """
    print("Starting annotation for file: {0}".format(sample.name))
    global total_samples
    global db_name_samples

    # Reduce the amount of variants to work with
    if not args.testmode:
        args_gz = shlex.shlex("bgzip {0} -f".format(sample.vcflocation.strip(".gz")))
        args_idx = shlex.shlex("tabix {0}".format(sample.vcflocation))

        args_gz.whitespace_split = True
        args_idx.whitespace_split = True
        if args.recall:
            if not os.path.exists(sample.vcflocation + ".bak"):
                shutil.copy(sample.vcflocation, sample.vcflocation + ".bak")
            gz_proc = subprocess.Popen([a for a in args_gz], shell=False)
            gz_proc.wait()
            idx_proc = subprocess.Popen([a for a in args_idx], shell=False)
            idx_proc.wait()

        input_vcf = sample.vcflocation

        annotate_args = shlex.split("perl {0} {1} {2} -buildver hg19 "
                                    "-out {3} "
                                    "-remove -protocol "
                                    "refGene,cytoBand,avsnp150,gnomad_genome,gnomad_exome,ljb26_all,clinvar_20180603,cosmic84 "
                                    "-argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' "
                                    "-operation g,r,f,f,f,f,f,f "
                                    "-nastring . "
                                    "-otherinfo "
                                    "-vcfinput".format(config.annotator, input_vcf, config.annotation_db,
                                                       sample.name))
        if not sample.error:
            proc = subprocess.Popen(annotate_args, shell=False, stderr=subprocess.PIPE)
            with proc.stderr:
                logdata(proc.stderr)
            processes.append(proc)

            proc.wait()

    # annotator output file --> add custom and external gene name based
    disease_name = os.path.join(config.custom_annotation_dir, "gene.omim_disease_name.synonyms.txt")
    disease_nr = os.path.join(config.custom_annotation_dir, "gene.disease.txt")
    hpo = os.path.join(config.custom_annotation_dir, "gene.hpoterm.txt")
    panels = os.path.join(config.custom_annotation_dir, "TSC_genepanels.txt")

    args_1 = shlex.shlex(
        "{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), disease_name, "Disease.name"))
    args_2 = shlex.shlex(
        "{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), disease_nr, "Disease.nr"))
    args_3 = shlex.shlex("{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), hpo, "HPO"))
    args_4 = shlex.shlex("{0} {1} {2}".format("python " + os.path.abspath(vcf_manipulator.__file__), panels, "Panel"))
    args_5 = shlex.shlex('java -jar {0} '
                         'extractFields '
                         '- '
                         '-e . -s ";" CHROM POS cytoBand avsnp150 REF ALT QUAL FILTER DP GEN[0].AD[0] GEN[0].VF '
                         'Gene.refGene Func.refGene GeneDetail.refGene ExonicFunc.refGene AAChange.refGene '
                         'gnomAD_genome_ALL gnomAD_genome_NFE gnomAD_exome_ALL gnomAD_exome_NFE gnomAD_exome_FIN SIFT_score SIFT_pred '
                         'Polyphen2_HVAR_score Polyphen2_HVAR_pred MutationTaster_score MutationTaster_pred '
                         'CADD_raw CADD_phred phyloP46way_placental phyloP100way_vertebrate CLNSIG CLNALLELEID CLNDN '
                         'CLNREVSTAT CLNDISDB cosmic84 '
                         'Disease.name Disease.nr HPO Panel GEN[0].GT GEN[0].DP GEN[0].AD'
                         .format(config.snpsift))
    if args.testmode:
        total_samples = 1
    args_6 = shlex.shlex('python {0}'.format(os.path.abspath(adsplit.__file__)))

    # This adds the local database frequencies.
    args_7 = shlex.shlex(
        'python {0} db_freq --annotation {1}.txt --totalsamples {2} --positions {3} --insertpos {4}'.format(
            os.path.abspath(annotate_by_pos.__file__),
            os.path.join(config.db_directory,
                                                                              db_name_samples),
            total_samples, "0,1,4,5",
            str(21)))  # "0,1,4,5" are column positions
    args_8 = shlex.shlex('bash {0} {1} {2} {3}'.format(config.itdseek, sample.bamlocation, config.reference,
                                                       config.samtools))
    args_9 = shlex.shlex('java -jar {0} '
                         'extractFields '
                         '- '
                         '-e . -s ";" CHROM POS ID REF ALT QUAL FILTER DP2 LEN SEQ CLIPPING INSERTION VAF'.format(
        config.snpsift))
    # This approach retains quotation marks and complete whitespace delimited args
    slx = list([args_1, args_2, args_3, args_4, args_5, args_6, args_7, args_8, args_9])
    for arg in slx:
        arg.whitespace_split = True

    if not sample.error:
        with open(sample.name + ".hg19_multianno.vcf") as annotated:
            proc_1 = subprocess.Popen([a for a in args_1], shell=False, stdin=annotated, stdout=subprocess.PIPE)
            proc_2 = subprocess.Popen([a for a in args_2], shell=False, stdin=proc_1.stdout, stdout=subprocess.PIPE)
            proc_3 = subprocess.Popen([a for a in args_3], shell=False, stdin=proc_2.stdout, stdout=subprocess.PIPE)
            proc_4 = subprocess.Popen([a for a in args_4], shell=False, stdin=proc_3.stdout, stdout=subprocess.PIPE)
            proc_5 = subprocess.Popen([a for a in args_5], shell=False, stdin=proc_4.stdout, stdout=subprocess.PIPE)
            # Proc 6 takes in a table
            proc_6 = subprocess.Popen([a for a in args_6], shell=False, stdin=proc_5.stdout, stdout=subprocess.PIPE)
            proc_7 = subprocess.Popen([a for a in args_7], shell=False, stdin=proc_6.stdout, stdout=subprocess.PIPE)
            # with proc_7.stderr:
            #    logdata(proc_7.stderr)
            with open(sample.name + ".annotated.table", "w+") as table:
                table.write(proc_7.communicate()[0])
                sample.table_files.append(os.path.abspath(table.name))

        if proc_7.returncode == 0:
            sample.annotated = True
            print("Finished annotating sample {0}".format(sample.name))
            print(sample)

        # Only do this if ITDSeek has been installed
        # https://github.com/tommyau/itdseek
        if config.itdseek is not None:
            with open(sample.name + ".itdseek.table", "w+") as itdseek:
                proc_8 = subprocess.Popen([a for a in args_8], shell=False, stdout=subprocess.PIPE)
                proc_9 = subprocess.Popen([a for a in args_9], shell=False, stdin=proc_8.stdout, stdout=subprocess.PIPE)
                itdseek.write(proc_9.communicate()[0])
                sample.table_files.append(os.path.abspath(itdseek.name))

    else:
        # There was a problem with SelectVariants
        sample.annotated = False


def calc_coverage(sample):
    path = os.path.join(workingDir, sample.name + "_coverage")
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    depth_args = shlex.split('java -Xmx8g -jar {0} '
                             '-T DepthOfCoverage '
                             '-R {1} '
                             '-I {2} '
                             '-L {3} '
                             '-geneList {4} '
                             '-ct 20 -ct 100 -ct 500 -ct 1000 '
                             '-o {5}'.format(config.toolkit, config.reference, sample.bamlocation, config.targetfile,
                                             config.refseq, os.path.join(path, sample.name + ".requested")))
    diagnose_args = shlex.split('java -Xmx4g -jar {0} '
                                '-T DiagnoseTargets '
                                '-R {1} '
                                '-I {2} '
                                '-L {3} '
                                '-min 20 '
                                '-o {4}'.format(config.toolkit, config.reference, sample.bamlocation, config.targetfile,
                                                os.path.join(path, sample.name + ".diagnoseTargets")))
    table_args = shlex.split('java -Xmx4g -jar {0} '
                             '-T VariantsToTable '
                             '-R {1} '
                             '-V {2} '
                             '-F CHROM -F POS -F END -F FILTER -F IDP -F IGC '
                             '--showFiltered '.format(config.toolkit, config.reference,
                                                      os.path.join(path, sample.name + ".diagnoseTargets")))
    depth = subprocess.Popen(depth_args, shell=False)
    # diagnose = subprocess.Popen(diagnose_args, shell=False)
    # table_name = os.path.join(path, sample.name + ".diagnoseTargets.table")
    # diagnose.wait()
    depth.wait()
    # to_table = subprocess.Popen(table_args, shell=False, stdout=subprocess.PIPE)
    # with open(table_name, "wb+") as f:
    # f.write(to_table.communicate()[0])
    #    pass

    # with to_table.stderr:
    #    logdata(to_table.stderr)

    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_summary')))
    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_gene_summary')))

    # Annotate the interval summary file based on BED file targets
    pipeline_utility.file_utility.add_symbols_to_interval_summary(config.targetfile, os.path.abspath(
        os.path.join(workingDir, sample.name + "_coverage", sample.name + '.requested.sample_interval_summary')),
                                                                  os.path.abspath(
                                                                      os.path.join(workingDir,
                                                                                   sample.name + "_coverage",
                                                                                   sample.name +
                                                                                   ".annotated.sample_interval_summary")))

    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage", sample.name +
                                                           ".annotated.sample_interval_summary")))
    sample.error = depth.returncode


def create_excel_table(sample):
    filters = MyeloidFilters(sample.table_files)
    post = MyeloidPostprocess(sample.table_files)
    formats = MyeloidFormats(sample.table_files)
    if args.annotate and sample.annotated:
        if args.recall:
            create_excel(".".join([sample.name, "BQ10", str(config.padding), "xlsx"]), sample.table_files,
                         filters,
                         post, formats)
        else:
            create_excel(".".join([sample.name, str(config.padding), "xlsx"]), sample.table_files, filters, post,
                         formats)

    elif args.testmode:
        sample.table_files.append(sample.name + ".annotated.table")
        create_excel(".".join([sample.name, str(config.padding), "xlsx"]), sample.table_files, filters, post, formats)

    else:
        sys.stderr.write("PIPELINE ERROR: Cannot create excel file for {0} "
                         "due to incomplete annotations!\n".format(sample.name))
    if args.annotate and sample.annotated:
        sample.finished = True
    elif args.annotate and not sample.annotated:
        sample.finished = False
        sys.stderr.write("Coverage was ordered but was not finished for sample {0}!\n".format(sample.name))


def single_sample(sample):
    try:
        if args.annotate:
            annotate(sample, args, config)
            sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.vcf"))
            sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.txt"))
            sample.trash.append(os.path.join(workingDir, sample.name + ".annotated.table"))
            sample.trash.append(os.path.join(workingDir, sample.name + ".itdseek.table"))
            sample.trash.append(os.path.join(workingDir, sample.name + ".avinput"))
            # sample.trash.append(os.path.join(workingDir, sample.name + ".converted.vcf"))
        if args.coverage:
            calc_coverage(sample)
            sample.trash.append(os.path.join(workingDir, sample.name + "_coverage", sample.name + ".requested"))
            # These are used with DiagnoseTargets
            # sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
            #                                 sample.name + ".sample_cumulative_coverage_counts"))
            # sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
            #                                 sample.name + ".sample_cumulative_coverage_proportions"))
            # sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
            #                                 sample.name + ".sample_interval_statistics"))
            # sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
            #                                 sample.name + ".sample_statistics"))

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
            # finished_samples.append(sample)
        else:
            print("Could not finish sample {0}".format(sample))
            # unfinished_samples.append(sample)
        return sample
    except Exception as error:
        sys.stderr.write("WORKER ERROR for sample: {0}\nTrace: ".format(sample))
        traceback.print_exc(file=sys.stderr)
        return sample
        # raise error


def cleanup(processpool):
    processpool.terminate()
    print("Exiting...")


def run_samples(args, sample_list):
    global cluster

    timeStart = datetime.datetime.now()
    samples = list()
    finished = JobStats()
    for sample in sample_list:
        location = pipeline_utility.file_utility.find_file(workingDir, sample + ".vcf.gz")[1]
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
    print("Found {0} samples:".format(len(samples)))
    for sample in samples:
        print(sample)
    local_db_tool.update_database(samples, args, config, combine_variants, update_sample_stats,
                                  db_name_samples, total_samples, idx=".tbi")

    pool = Pool(processes=args.ncpus)
    atexit.register(cleanup, pool)
    return_samples = pool.map(single_sample, samples)

    timeEnd = datetime.datetime.now()

    for sample in return_samples:
        print(sample)
        if sample.finished:
            finished.finished_samples.append(sample)
        else:
            finished.unfinished_samples.append(sample)
    print("-" * 40)
    print("Finish time: {0}".format(timeEnd))
    print("Finished {0}/{1} of ordered/found samples in {2}:".format(len(finished.finished_samples), len(samples),
                                                                     str(timeEnd - timeStart)))
    print("-" * 40)
    if len(finished.unfinished_samples) > 0:
        for sample in finished.unfinished_samples:
            print("{0} is unfinished. Check for errors. ".format(sample))
        print("You can rerun with --samples {0}".format(
            str(list(s.name for s in finished.unfinished_samples)).strip("[]'")))


def main(args):
    samples = list()
    vcfslist = list()

    if args.overwrite:
        sys.stderr.writelines(
            "WARNING: VCF files with file names that already exist in the database will be overwritten! "
            "Old variant files will not exist. (--overwrite is active)\n")

    if args.batch:
        # prefixes, vcfslist, bamlist = create_configs()
        # for i in range(0, len(prefixes)):
        #    samp = Sample(prefixes[i], vcfslist[i], bamlist[i])
        prefixes = pipeline_utility.file_utility.find_prefixes(workingDir, "vcf.gz")
        if prefixes is not None:
            run_samples(args, prefixes)

    elif args.samples:
        print("Input is {0} samples: {1}".format(len(args.samples), "\t".join(args.samples)))
        run_samples(args, args.samples)
    else:
        print "No valid command input."
        print "Printing config values (including hidden variables)."
        print config


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(prog="Trusight Myeloid pipeline command-line tool.")
        group = parser.add_mutually_exclusive_group()
        group.add_argument("-b", "--batch",
                           help="Find all unique .vcf files and their matching .bams. "
                                "Program will only run if each vcf has a matching .bam file.",
                           action="store_true")
        group.add_argument("-s", "--samples",
                           help="Followed by a list of unique sample identifiers e.g. "
                                "-s E0000001 E000002 E0000003 which are to be run through the pipeline. "
                                "Will only run if samples have a matching vcf and bam file.",
                           action="store", nargs='+', type=str)

        parser.add_argument("-r", "--overwrite",
                            help="Overwrites the VCF file in the db with the current sample if it already exists "
                                 "in the DB."
                                 "Useful in a test scenario when running large batches repeatedly or if the VCF is "
                                 "not supposed to be inserted into the DB for each run.",
                            action="store_true", default=False)
        parser.add_argument("-t", "--testmode", help="Skips the variant DB recompilation (CombineVariants) step, "
                                                     "expects the DB to exist.",
                            action="store_true", default=False)
        parser.add_argument("-k", "--keep", help="Keep intermediary annotation files.",
                            action="store_true", default=False)
        parser.add_argument("-a", "--annotate", default=False, action="store_true",
                            help="The output will contain a table with the annotated variants. "
                                 "Run this explicitly for annotations")
        parser.add_argument("-c", "--coverage",
                            help="Calculates the coverage tables for the sample's order and adds it to "
                                 "the output excel file new pages.",
                            action="store_true", default=False)
        parser.add_argument("--ncpus", help="The number of CPUs dedicated to the annotation process.", default=1,
                            action="store", type=int)
        parser.add_argument("--recall", action="store_true", help="Takes in the .vcf file, " \
                                                                  "bgzips and uses tabix to index it.", default=False)

        args = parser.parse_args()
        main(args)

    except Exception:
        raise
