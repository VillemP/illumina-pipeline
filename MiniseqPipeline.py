import os

import ConfigValidation
import DatabaseCrawler

vcf_storage_location = "/media/kasutaja/data/TSC_temp/miniseq_pipe/vcfs/"
vcf_list_name = "vcfs-sample-path.list"
db_location = "/media/kasutaja/data/NGS_data/var_db_miniseq/"
vcf_dir = os.path.join(db_location, vcf_list_name)
db_name = "miniseq-named-targeted-merged-n"
db_dir = os.path.join(db_location, db_name)
workingDir = os.getcwd()


def create_configs():
    prefixes = DatabaseCrawler.write_prefixes_list(workingDir, "prefixes.list")
    vcflist = DatabaseCrawler.write_vcfs_list(workingDir, "vcfs.list")
    bamlist = DatabaseCrawler.write_bams_list(workingDir, "bams.list")

    if ConfigValidation.validate_config("bams.list", "vcfs.list", "prefixes.list"):
        return prefixes, vcflist
    else:
        raise IOError("Database updating failed due to incorrect files.")


def update_vcf_list(samples, overwrite=False):
    assert os.path.exists(vcf_dir)
    global db_name
    data = ""
    curr_len = DatabaseCrawler.file_len(vcf_dir)

    with open(vcf_dir, "a+") as db_vcfs:
        data = db_vcfs.readlines()
        # db_vcfs.seek(0)

        i = 0
        skipped = 0
        for vcf in samples:
            name = os.path.basename(vcf).rsplit('.')[0]
            if not any(name in x.rstrip() for x in data):
                line = "V:{0} {1}".format(name, vcf)
                i += 1
                db_vcfs.write(line + "\n")
            else:
                skipped += 1
    db_name = db_name + str((curr_len + i))
    print ("Updated {0} with {1} unique samples. "
           "Skipped {2} preexisting samples. New database name is {3}.").format(
        os.path.realpath(vcf_list_name), i, skipped, db_name)


def combine_variants(vcfs_list):
    pass


def update_database():
    prefixes, vcfslist = create_configs()

    DatabaseCrawler.copy_vcf(vcfslist, vcf_storage_location, True)
    update_vcf_list(vcfslist, True)
    #combine_variants(vcf_dir)


if __name__ == "__main__":
    try:
        update_database()
    except Exception:
        raise
# MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq
