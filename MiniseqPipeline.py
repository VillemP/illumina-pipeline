import argparse
import os
import shlex
import subprocess
import sys

import miniseq.configvalidator
import pipeline_utility.file_utility
import pipeline_utility.sample
# TODO: Run a set of commands from STDOUT -> STDIN
from miniseq.miniseq_excel_filters import MiniseqFilters, MiniseqPostprocess, MiniseqFormats
from miniseq.miniseqconfig import MiniseqConfig
from pipeline_utility import vcf_manipulator, adsplit, annotate_by_pos
from pipeline_utility.txttoxlsx_filtered import create_excel

try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

processes = []
# Location of the cfg file is in the same folder as this script
cfg_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'miniseq.yaml')
config = MiniseqConfig(cfg_path)
config = config.load()
total_samples = 0
db_name_samples = config.db_name  # Updated later with current samples name
workingDir = os.getcwd()
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

    try:
        os.makedirs(config.db_directory)
    except OSError:
        if not os.path.isdir(config.db_directory):
            raise

    with open(config.db_vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()

        curr_len = pipeline_utility.file_utility.file_len(config.db_vcf_dir)

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


def update_database(samples, replace, testmode):
    if not testmode:
        assert len(samples) > 0, "List of samples to be updated into the database cannot be empty!"
        vcfslist = list()
        for sample in samples:
            vcfslist.append(sample.vcflocation)
        copied, renamed, skipped = pipeline_utility.file_utility.copy_vcf(vcfslist, config.vcf_storage_location,
                                                                          replace)
        updated, skip = update_vcf_list(vcfslist)
        if len(copied + renamed) + updated > 0:
            # If there were any variant files updated (copied) or new variants added to vcf list
            combine_variants()
    pass


def annotate(sample, testmode):
    print("Starting annotation for file: {0}".format(sample.name))
    global total_samples
    global db_name_samples
    # Reduce the amount of variants to work with
    if not testmode:
        outfile = "{0}.targeted.padding{1}bp.vcf".format(sample.name, config.padding)
        args = shlex.split("java -Xmx10g -jar {0} "
                           "-T SelectVariants "
                           "-R {1} "
                           "-V {2} "
                           "-L {3} "
                           "-o {4} "
                           "-ip {5}".format(config.toolkit, config.reference,
                                            sample.vcflocation, config.targetfile,
                                            outfile, config.padding))

        proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
        processes.append(proc)

        with proc.stderr:
            logdata(proc.stderr)
        processes.append(proc)

        proc.wait()
        if proc.returncode == 0:
            sample.reduced_variants_vcf = outfile
        else:
            sample.error = True


        args = shlex.split("perl {0} {1} {2} -buildver hg19 "
                           "-out {3} "
                           "-remove -protocol "
                           "refGene,avsnp147,1000g2015aug_all,1000g2015aug_eur,exac03,ljb26_all,clinvar_20170130 "
                           "-argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' "
                           "-operation g,f,f,f,f,f,f "
                           "-nastring . "
                           "-otherinfo "
                           "-vcfinput".format(config.annotator, outfile, config.annotation_db, sample.name))
        if not sample.error:
            proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
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
                         '-e . -s ";" CHROM POS avsnp147 REF ALT QUAL FILTER AC AF DP '
                         'Gene.refGene Func.refGene GeneDetail.refGene ExonicFunc.refGene AAChange.refGene '
                         '1000g2015aug_all 1000g2015aug_eur ExAC_ALL ExAC_NFE ExAC_FIN SIFT_score SIFT_pred '
                         'Polyphen2_HVAR_score Polyphen2_HVAR_pred MutationTaster_score MutationTaster_pred '
                         'CADD_raw CADD_phred phyloP46way_placental phyloP100way_vertebrate CLINSIG CLNDBN CLNACC '
                         'CLNDSDB CLNDSDBID '
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
    slx = list([args_1, args_2, args_3, args_4, args_5, args_6, args_7])
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

        # TODO: Switch to multiprocessing?
        # jobs = multiprocessing.Pool(1)
        # annotator = multiprocessing.Process()

        if proc_7.returncode == 0:
            sample.annotated = True
            print("Finished annotating sample {0}".format(sample.name))
            print(sample)
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
    depth_args = shlex.split('java -Xmx4g -jar {0} '
                             '-T DepthOfCoverage '
                             '-R {1} '
                             '-I {2} '
                             '-L {3} '
                             '-geneList {4} '
                             '-ct 1 -ct 10 -ct 20 -ct 50 '
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
    diagnose = subprocess.Popen(diagnose_args, shell=False)
    table_name = os.path.join(path, sample.name + ".diagnoseTargets.table")
    diagnose.wait()
    depth.wait()
    # to_table = subprocess.Popen(table_args, shell=False, stdout=subprocess.PIPE)
    with open(table_name, "wb+") as f:
        # f.write(to_table.communicate()[0])
        pass

    # with to_table.stderr:
    #    logdata(to_table.stderr)

    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_summary')))
    sample.table_files.append(os.path.abspath(os.path.join(workingDir, sample.name + "_coverage",
                                                           sample.name + '.requested.sample_gene_summary')))


def create_excel_table(sample):
    filters = MiniseqFilters(sample.table_files)
    post = MiniseqPostprocess(sample.table_files)
    formats = MiniseqFormats(sample.table_files)
    if sample.annotated:
        create_excel(".".join([sample.name, str(config.padding), "xlsx"]), sample.table_files, filters, post, formats)
    else:
        sys.stderr.write("PIPELINE ERROR: Cannot create excel file for {0} "
                         "due to incomplete annotations!\n".format(sample.name))


def run_samples(args, sample_list):
    samples = list()
    finished_samples = []
    unfinished_samples = []
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
    print("Found {0} samples:".format(len(samples)))
    for sample in samples:
        print(sample)
    update_database(samples, args.no_replace, args.test)
    for sample in samples:
        if not args.test:
            annotate(sample, args.test)
            if not sample.error:
                calc_coverage(sample)
                create_excel_table(sample)
            if not args.keep and sample.annotated:
                sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.vcf"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".hg19_multianno.txt"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".annotated.table"))
                sample.trash.append(os.path.join(workingDir, sample.name + ".avinput"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage", sample.name + ".requested"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".requested.sample_cumulative_coverage_counts"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".requested.sample_cumulative_coverage_proportions"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".requested.sample_interval_statistics"))
                sample.trash.append(os.path.join(workingDir, sample.name + "_coverage",
                                                 sample.name + ".requested.sample_statistics"))

                for trashfile in sample.trash:
                    try:
                        os.remove(trashfile)
                    # The file might've been already deleted or did not exist in the first place.
                    except(OSError) as oserror:
                        sys.stderr.write("Couldn't remove file: {0}\n{1}\n".format(trashfile, oserror))
                    finally:
                        if sample.annotated:
                            sample.finished = True

            if sample.finished and not sample.error:
                print("Finished sample {0}".format(sample.name))
                finished_samples.append(sample)
            else:
                print("Could not finish sample {0}".format(sample))
                unfinished_samples.append(sample)
    print("-" * 40)
    print("Annotated {0}/{1} of ordered/found samples:".format(len(finished_samples), len(samples)))
    print("-" * 40)
    for sample in samples:
        print(sample)
    if len(unfinished_samples) > 0:
        for sample in unfinished_samples:
            print("{0} is unfinished. Check for errors. ".format(sample))
        print("You can rerun with --samples {0} (strip the brackets)".format(s.name for s in unfinished_samples))


def main(args):
    samples = list()
    vcfslist = list()

    # True if replace, False if --no_replace, passed as an argument in update_database()
    if not args.no_replace:
        sys.stderr.writelines(
            "WARNING: VCF files with file names that already exist in the database will not be updated! "
            "Old variant files still exist. (--no_replace is active)\n")

    if args.batch:
        # prefixes, vcfslist, bamlist = create_configs()
        # for i in range(0, len(prefixes)):
        #    samp = Sample(prefixes[i], vcfslist[i], bamlist[i])
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
    else:
        print "No valid command input."
        print "Printing cfg values."
        print config


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(prog="MiniSeq pipeline command-line tool.")
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
        group.add_argument("-o", "--old", help="Creates text based config files for shell scripts. "
                                               "Automatically updates the variant database.",
                           action="store_true")
        parser.add_argument("-r", "--no_replace",
                            help="Will skip copying the VCF file to the VCF directory specified in the config file"
                                 " (don't overwrite VCF in the DB)."
                                 "Useful in a test scenario when running large batches repeatedly or if the VCF is "
                                 "not supposed to be inserted into the DB.",
                            action="store_false", default=True)
        parser.add_argument("-t", "--test", help="Skips the variant DB recompilation (CombineVariants) step, "
                                                 "expects the DB to exist..",
                            action="store_true", default=False)
        parser.add_argument("-k", "--keep", help="Keep intermediary annotation files.",
                            action="store_true", default=False)

        args = parser.parse_args()
        main(args)
    except Exception:
        raise
    finally:
        for p in processes:
            try:
                p.terminate()
            except OSError:
                pass  # process is already dead
