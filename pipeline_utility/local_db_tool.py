# TODO: Replace the DB handling in MiniseqPipeline.py and illumina_pipe.py so that both use the same commands
import os
import shlex

import miniseq
import miniseq.configvalidator
import pipeline_utility
import pipeline_utility.file_utility

workingDir = ""
project = ""


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


def update_vcf_list(vcfs_list, config):
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


def combine_variants(out, processes, logdata, config, vcflist):
    # Combining variant files into a single reference to be used for statistical purposes
    args = shlex.split('java -Xmx10g -jar {0} '
                       '-T CombineVariants '
                       '-R {1} '
                       '-V {2} '
                       '-o {3} '
                       '-log {4} '
                       '--genotypemergeoption UNIQUIFY'.format(config.toolkit, config.reference, vcflist,
                                                               out, config.logfile))

    proc = os.subprocess.Popen(args, shell=False, stderr=os.subprocess.PIPE)
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

    proc = os.subprocess.Popen(args, shell=False, stderr=os.subprocess.PIPE)
    processes.append(proc)
    logdata(proc.stderr)
    proc.wait()
    # print proc.returncode
    # proc.communicate()


def update_database(samples, replace, testmode, config):
    if not testmode:
        assert len(samples) > 0, "List of samples to be updated into the database cannot be empty!"
        vcfslist = list()
        for sample in samples:
            vcfslist.append(sample.vcflocation)
        copied, skipped = pipeline_utility.file_utility.copy_vcf(vcfslist, config.vcf_storage_location, replace)
        updated, skip = update_vcf_list(vcfslist, config)
        if copied + updated > 0:
            # If there were any variant files updated (copied) or
            # TODO: Fill the parameters
            combine_variants(os.path.join(config.db_directory, db_name_samples + ".vcf"))
    pass