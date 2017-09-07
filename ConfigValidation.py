import os, sys, numpy


def validate_config(bams, vcfs, prefixes):
    assert os.path.exists(bams)
    assert os.path.exists(vcfs)
    assert os.path.exists(prefixes)

    bamfile = open(bams, "r")
    vcffile = open(vcfs, "r")
    prefixfile = open(prefixes, "r")
    config_files = (bamfile, vcffile, prefixfile)

    bamdata = numpy.asanyarray(bamfile.readlines())
    vcfdata = numpy.asanyarray(vcffile.readlines())
    prefixdata = numpy.asanyarray(prefixfile.readlines())

    bam_vcf_equal = len(bamdata) == len(vcfdata)
    bam_prefix_equal = len(bamdata) == len(prefixdata)
    if bam_vcf_equal and bam_prefix_equal:
        lineerror = False
        for config in config_files:
            count = 0
            if not validate_file(config):
                raise ValueError("Scripts will not run with the current configs! "
                                 "Check your bam, vcf, prefix lists for file length and line endings.")

        for i in range(0, len(prefixdata)):
            if not (prefixdata[i].rstrip() in bamdata[i].rstrip() and prefixdata[i].rstrip() in vcfdata[i].rstrip()):
                print("Config files on line {0} do not match!".format(i+1))
                print("Prefix: {0}\nVCF: {1}\nBAM: {2}".format(
                    prefixdata[i], vcfdata[i], bamdata[i]))
                lineerror = True
                # raise ValueError(prefixdata[i], bamdata[i], vcfdata[i])
        if not lineerror:
            print("Config files have no discernible errors and are ready to start analysis.")
            return True

    elif bam_vcf_equal or bam_prefix_equal:
        print "Files {0} and {1} do not match in length!".format(
            (bamfile.name, len(bamdata)), (vcffile.name, len(vcfdata)) if not bam_vcf_equal
            else (prefixfile.name, len(prefixdata))
            if not bam_vcf_equal or not bam_prefix_equal else None)
    else:
        print("Files {0} length:{1} lines\n{2} length: {3}\n{4} length:{5} are of inequal length!".format(
            bamfile.name, len(bamdata), vcffile.name,
            len(vcfdata), prefixfile.name, len(prefixdata)))

    for config in config_files:
        config.close()

    return False


# Check if file exists, doesn't have any trailing whitespace, ASCII only
def validate_file(f):
    assert f is not None
    f.seek(0)
    text = f.read()
    if len(text) <= 1:
        print("File {0} is empty! Length of file: {1}".format(f.name, len(text)))
        return False
    if '\r\n' in text:
        print("File {0} contains DOS line endings!".format(f.name))
        return False
    return True


if __name__ == "__main__":
    args = sys.argv
    validate_config(args[1], args[2], args[3])
