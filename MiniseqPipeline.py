import DatabaseCrawler
import ConfigValidation
import os

vcf_location = "/media/kasutaja/data/TSC_temp/miniseq_pipe/vcfs-test/"
vcf_list_name = "vcfs-sample-path-test.list"
vcf_dir = os.path.join(vcf_location, vcf_list_name)
db_location = "/media/kasutaja/data/NGS_data/var_db_miniseq/"
db_name = "miniseq-named-targeted-merged-n"
db_dir = os.path.join(db_location, db_name)
workingDir = os.getcwd()


def create_configs():
    samples = DatabaseCrawler.find_vcfs(workingDir)

    with open("prefixes.list", "wb") as prefixes:
        for prefix in samples.iterkeys():
            prefixes.write(prefix.rsplit('.')[0] + "\n")
    return samples


def update_vcf_list(samples):
    assert os.path.exists(vcf_dir)
    global db_name
    data = ""

    with open(vcf_dir, "r") as db_vcfs:
        data = db_vcfs.read()

    with open(vcf_dir, "a") as db_vcfs:
        i = 0
        s = 0

        curr_len = len(data.splitlines(True))
        for sample in samples.iteritems():
            if sample[0] not in data:
                fullpath = os.path.join(sample[1], sample[0])
                line = "V:{0} {1}".format(sample[0].rsplit('.')[0], fullpath)
                i += 1
                db_vcfs.write(line + "\n")
            else:
                s += 1
    db_name = db_name + str((curr_len + i))
    print ("Updated {0} with {1} unique samples. "
           "Skipped {2} samples. New database name is {3}.").format(
        vcf_list_name, i, s, db_name)


def combine_variants(vcfs_list):
    pass


def update_database():
    global db_name
    samples = create_configs()
    batch_size = len(samples)
    i = 0


    db_name = db_name + str(i)
    DatabaseCrawler.copy_vcf(samples, vcf_location)
    update_vcf_list(samples)
    #combine_variants(vcf_dir)


if __name__ == "__main__":
    try:
        update_database()
    except (Exception):
        raise



# MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq
