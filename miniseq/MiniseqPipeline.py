import os
import shlex
import subprocess

import configvalidator
import file_utility

vcf_storage_location = "/media/kasutaja/data/TSC_temp/illumina_pipe/vcfs/"
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
            log.write(line + "\n")


def create_configs():
    prefixes = file_utility.write_prefixes_list(workingDir, "prefixes.list")
    vcflist = file_utility.write_vcfs_list(workingDir, "vcfs.list")
    bamlist = file_utility.write_bams_list(workingDir, "bams.list")

    if configvalidator.validate_config("bams.list", "vcfs.list", "prefixes.list"):
        return prefixes, vcflist, bamlist
    else:
        raise IOError("Database updating failed due to incorrect files.")


def create_arguments_file(dbnr):
    with open("arguments.txt", "wb+") as f:
        f.write("wd:\n")
        f.write(workingDir + "\n")
        f.write("dbname:\n")
        f.write(dbnr + "\n")
        f.write("project:\n")
        f.write(project)


def update_vcf_list(samples, overwrite=False):
    global db_name
    data = ""

    with open(db_vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()
        # db_vcfs.seek(0)
        curr_len = file_utility.file_len(db_vcf_dir)

        i = 0
        skipped = 0
        for vcf in samples:
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
                       '-L /media/kasutaja/data/TSC_temp/illumina_pipe/coverage/trusight_cancer_manifest_aUsed.bed '
                       '-o {2}{3}.vcf '
                       '-ip 10 '
                       '--genotypemergeoption UNIQUIFY '.format(db_location, db_vcf_list_name, db_location, db_name))

    print args
    proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    processes.append(proc)

    with proc.stdout:
        logdata(proc.stdout)
    proc.wait()

    args = shlex.split('java -Xmx10g -jar /home/sander/NGS_programs/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar '
                       '-T VariantsToTable -R '
                       '/media/kasutaja/data/NGS_data/hg19/ucsc.hg19.fasta '
                       '-V {0}{1}.vcf '
                       '-F CHROM -F POS -F REF -F ALT -F AC -F HET -F HOM-VAR '
                       '--splitMultiAllelic --showFiltered '
                       '-o {2}{3}.txt '.format(db_location, db_name, db_location, db_name))

    proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    with proc.stdout:
        logdata(proc.stdout)
    proc.wait()
    # print proc.returncode


def update_database():
    prefixes, vcfslist, bamlist = create_configs()

    file_utility.copy_vcf(vcfslist, vcf_storage_location, True)
    update_vcf_list(vcfslist, True)
    # create_arguments_file()
    combine_variants()


if __name__ == "__main__":
    try:
        update_database()
    except Exception:
        raise
    finally:
        for p in processes:
            p.kill()
# MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq
