import argparse
import os
import shlex
import subprocess

import miniseq.configvalidator
import miniseq.file_utility
# TODO: store variables in a config file
import pipeline_utility.sample

# TODO: Create a config file system
# TODO: Run a set of commands from STDOUT -> STDIN
vcf_storage_location = "/media/kasutaja/data/TSC_temp/miniseq_pipe/vcfs/"
db_vcf_list_name = "vcfs-sample-path.list"
db_location = "/media/kasutaja/data/NGS_data/var_db_miniseq/"
db_vcf_dir = os.path.join(db_location, db_vcf_list_name)
db_name = "miniseq-named-targeted-merged-n"
db_dir = os.path.join(db_location, db_name)
workingDir = os.getcwd()
project = os.path.basename(os.path.normpath(workingDir))
logfile = "Miniseq-log-{0}.txt".format(project)

processes = []


def logdata(stdout):
    with open(logfile, "a+") as log:
        for line in iter(stdout.readline, b''):  # b'\n'-separated lines
            # log.write(time.time() + "\t" + line + "\n")
            log.write(line)
            print ("{}".format(line.rstrip()))


# soon to be deprecated
def create_configs():
    prefixes = miniseq.file_utility.write_prefixes_list(workingDir, "prefixes.list")
    vcflist = miniseq.file_utility.write_vcfs_list(workingDir, "vcfs.list")
    bamlist = miniseq.file_utility.write_bams_list(workingDir, "bams.list")

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
    global db_name
    data = ""

    with open(db_vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()
        # db_vcfs.seek(0)
        curr_len = miniseq.file_utility.file_len(db_vcf_dir)

        i = 0
        skipped = 0
        for vcf in vcfs_list:
            name = os.path.basename(vcf).rsplit('.')[0]
            if not any(name in x.rstrip() for x in data):
                line = "V:{0} {1}".format(name, os.path.join(vcf_storage_location,
                                                             os.path.basename(vcf)))
                i += 1
                db_vcfs.write(line + "\n")
            else:
                skipped += 1
    create_arguments_file(str(curr_len + i))
    db_name = db_name + str((curr_len + i))
    print ("Updated {0} with {1} unique samples. "
           "Skipped {2} preexisting samples. New database name is {3}.").format(
        os.path.join(db_location, db_vcf_list_name), i, skipped, db_name)


def combine_variants():
    args = shlex.split('java -Xmx10g -jar /home/sander/NGS_programs/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar '
                       '-T CombineVariants '
                       '-R /media/kasutaja/data/NGS_data/hg19/ucsc.hg19.fasta '
                       '-V {0}{1} '
                       '-L /media/kasutaja/data/TSC_temp/miniseq_pipe/coverage/trusight_cancer_manifest_aUsed.bed '
                       '-o {2}{3}.vcf '
                       '-ip 10 '
                       '--genotypemergeoption UNIQUIFY '.format(db_location, db_vcf_list_name, db_location, db_name,
                                                                logfile))

    proc = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
    processes.append(proc)

    with proc.stderr:
        logdata(proc.stderr)
    proc.wait()

    args = shlex.split('java -Xmx10g -jar /home/sander/NGS_programs/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar '
                       '-T VariantsToTable -R '
                       '/media/kasutaja/data/NGS_data/hg19/ucsc.hg19.fasta '
                       '-V {0}{1}.vcf '
                       '-F CHROM -F POS -F REF -F ALT -F AC -F HET -F HOM-VAR '
                       '--splitMultiAllelic --showFiltered '
                       '-o {2}{3}.txt '.format(db_location, db_name, db_location, db_name, logfile))

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
    miniseq.file_utility.copy_vcf(vcfslist, vcf_storage_location, replace)
    update_vcf_list(vcfslist, True)
    # create_arguments_file()
    combine_variants()


def main(args):
    # update_database(single sample)
    # annotate(sample)
    # calc_coverage(sample)
    samples = list()
    vcfslist = list()

    # True if replace, False if --no_replace
    if not args.no_replace:
        print("WARNING: VCF files that have overlapping names were not copied into the database! "
              "Old variant files still exist.")

    if args.batch:
        # prefixes, vcfslist, bamlist = create_configs()
        # for i in range(0, len(prefixes)):
        #    samp = Sample(prefixes[i], vcfslist[i], bamlist[i])
        prefixes = miniseq.file_utility.write_prefixes_list(workingDir, "prefixes.list")
        for prefix in prefixes:
            vcf = miniseq.file_utility.find_file(workingDir, prefix + ".vcf")[1]
            bam = miniseq.file_utility.find_file(workingDir, prefix + ".bam")[1]
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
    else:
        for sample in args.samples:
            vcfname, location = miniseq.file_utility.find_file(workingDir, sample + ".vcf")
            bamname, bamlocation = miniseq.file_utility.find_file(workingDir, sample + ".bam")
            samp = pipeline_utility.sample.Sample(sample, location, bamlocation)
            samples.append(samp)
            update_database(samples, args.no_replace)


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
        parser.add_argument("-u", "--no_replace", action="store_false", default=True)

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