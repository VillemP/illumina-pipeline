import argparse
import os
import shlex
import subprocess

import yaml

import miniseq.configvalidator
import pipeline_utility.file_utility
import pipeline_utility.sample
# TODO: Run a set of commands from STDOUT -> STDIN
from miniseq.miniseqconfig import MiniseqConfig
from pipeline_utility import vcf_manipulator

processes = []
yaml.add_constructor(MiniseqConfig.yaml_tag, MiniseqConfig.cfg_constructor)
cfg_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'miniseq.yaml')
config = MiniseqConfig(cfg_path)
config = config.load()
db_name_samples = config.db_name
workingDir = os.getcwd()
project = os.path.basename(os.path.normpath(workingDir))

assert os.path.exists(config.annotator)


def logdata(stdout):
    with open(config.logfile, "a+") as log:
        for line in iter(stdout.readline, b''):  # b'\n'-separated lines
            # log.write(time.time() + "\t" + line + "\n")
            log.write(line)
            print ("{}".format(line.rstrip()))


# soon to be deprecated
def create_configs():
    prefixes = pipeline_utility.file_utility.write_prefixes_list(workingDir, "prefixes.list")
    vcflist = pipeline_utility.file_utility.write_vcfs_list(workingDir, "vcfs.list")
    bamlist = pipeline_utility.file_utility.write_bams_list(workingDir, "bams.list")

    if miniseq.configvalidator.validate_config("bams.list", "vcfs.list", "prefixes.list"):
        return prefixes, vcflist, bamlist
    else:
        raise IOError("Database updating failed due to incorrect files.")


# soon to be deprecated
def create_arguments_file(dbnr):
    with open("arguments.txt", "wb+") as f:
        f.write("wd:\n")
        f.write(workingDir + "\n")
        f.write("dbname:\n")
        f.write(dbnr + "\n")
        f.write("project:\n")
        f.write(project)


def update_vcf_list(vcfs_list, overwrite=False):
    global db_name_samples
    data = ""

    with open(config.db_vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()
        # db_vcfs.seek(0)
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
    create_arguments_file(str(curr_len + i))
    db_name_samples = db_name_samples + str((curr_len + i))
    print ("Updated {0} with {1} unique samples. "
           "Skipped {2} preexisting samples. New database name is {3}.").format(
        os.path.join(config.db_directory, config.db_vcf_list_name), i, skipped, db_name_samples)


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

    with proc.stderr:
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
    with proc.stderr:
        logdata(proc.stderr)
    proc.wait()
    # print proc.returncode
    # proc.communicate()


def update_database(samples, replace):
    assert len(samples) > 0, "List of samples to be updated into the database cannot be empty!"
    vcfslist = list()
    for sample in samples:
        vcfslist.append(sample.vcflocation)
    pipeline_utility.file_utility.copy_vcf(vcfslist, config.vcf_storage_location, replace)
    update_vcf_list(vcfslist, True)
    combine_variants()


def annotate(sample):
    # Reduce the amount of variants to work with
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
    proc.wait()
    if proc.returncode == 0:
        sample.reduced_variants_vcf = outfile

    args = shlex.split("perl {0} {1} {2} -buildver hg19 "
                       "-out {3} "
                       "-remove -protocol "
                       "refGene,avsnp147,1000g2015aug_all,1000g2015aug_eur,exac03,ljb26_all,clinvar_20150629 "
                       "-argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' "
                       "-operation g,f,f,f,f,f,f "
                       "-nastring . "
                       "-otherinfo "
                       "-vcfinput".format(config.annotator, outfile, config.annotation_db, sample.name))

    proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
    processes.append(proc)

    with proc.stderr:
        logdata(proc.stderr)
    proc.wait()
    with open(sample.name + "hg19_multianno.vcf") as annotated_file:
        proc = subprocess.Popen(vcf_manipulator.annotate(annotated_file,
                                                         os.path.join(config.custom_annotation_dir,
                                                                      "gene.omim_disease_name.synonyms.txt",
                                                                      ), "Disease.name"),
                                shell=False)
        proc.communicate()



def main(args):
    samples = list()
    vcfslist = list()

    # True if replace, False if --no_replace, passed as an argument in update_database()
    if not args.no_replace:
        print("WARNING: VCF files that have overlapping names were not copied into the database! "
              "Old variant files still exist.")

    if args.batch:
        # prefixes, vcfslist, bamlist = create_configs()
        # for i in range(0, len(prefixes)):
        #    samp = Sample(prefixes[i], vcfslist[i], bamlist[i])
        prefixes = pipeline_utility.file_utility.write_prefixes_list(workingDir, "prefixes.list")
        for prefix in prefixes:
            vcf = pipeline_utility.file_utility.find_file(workingDir, prefix + ".vcf")[1]
            bam = pipeline_utility.file_utility.find_file(workingDir, prefix + ".bam")[1]
            samples.append(pipeline_utility.sample.Sample(prefix, vcf, bam))
        update_database(samples, args.no_replace)
        for sample in samples:
            # annotate(sample)
            # calc_coverage(sample)
            # create_excel_table(sample)
            pass

    elif args.old:
        prefixes, vcfslist, bamlist = create_configs()
        for i in range(0, len(prefixes)):
            samples.append(pipeline_utility.sample.Sample(prefixes[i], vcfslist[i], bamlist[i]))
        update_database(samples, args.no_replace)
    elif args.samples:
        for sample in args.samples:
            vcfname, location = pipeline_utility.file_utility.find_file(workingDir, sample + ".vcf")
            bamname, bamlocation = pipeline_utility.file_utility.find_file(workingDir, sample + ".bam")
            samp = pipeline_utility.sample.Sample(sample, location, bamlocation)
            samples.append(samp)
            update_database(samples, args.no_replace)
        for sample in samples:
            annotate(sample)
            # calc_coverage(sample)
            # create_excel_table(sample)
            pass
    else:
        print "No valid command input."
        print "Printing cfg values."
        print config


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(prog="MiniSeq pipeline command-line tool.")
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
        parser.add_argument("-r", "--no_replace", action="store_false", default=True)

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
