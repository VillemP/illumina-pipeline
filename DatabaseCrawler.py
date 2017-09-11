import os
import shutil
import sys

# MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq

commands = ['WriteVcfs', 'CopyVcfs']
duplicates = 0


def find_filetype(dir, filetype):
    assert os.path.exists(dir)
    duplicates = 0
    # unique_files = dict([])
    unique_files = list(())
    for (dirpath, dirnames, files) in os.walk(dir):
        for name in files:
            if name.endswith(filetype):
                if name not in unique_files:
                    unique_files.append((name, os.path.join(dirpath, name)))
                else:
                    duplicates += 1

    return unique_files


def write_filelist(dir, output, file_paths):
    try:
        with open(os.path.join(dir, output), "wb+") as f:
            all_paths = list()
            file_paths.sort(key=lambda pair: pair[0])
            for file_path_pair in file_paths:
                full_path = file_path_pair[1]
                all_paths.append(full_path)
                # print (full_path)
                f.write(full_path + "\n")
            print("Total files: {}".format(len(file_paths)))
            print("Duplicates: {}".format(duplicates))

        return all_paths

    except (NameError, KeyboardInterrupt, SyntaxError, AssertionError):
        if NameError or AssertionError:
            print("Unable to find {}".format(dir))
            raise
        if KeyboardInterrupt:
            print("Quitting.")


def find_vcfs(dir):
    return find_filetype(dir, '.vcf')


def find_bams(dir):
    return find_filetype(dir, '.bam')


def write_prefixes_list(dir, output):
    samples = find_vcfs(dir)
    samples.sort(key=lambda pair: pair[0])
    clean_prefixes = list()

    with open(output, "wb+") as prefixes:
        for prefix in samples:
            prefix_clean = prefix[0].rsplit('.')[0]
            prefixes.write(prefix_clean + "\n")
            clean_prefixes.append(prefix_clean)
    return clean_prefixes


def write_vcfs_list(dir, output):
    return write_filelist(dir, output, find_vcfs(dir))


def write_bams_list(dir, output):
    return write_filelist(dir, output, find_bams(dir))


def copy_vcf(files, dest, overwrite=False):
    assert os.path.exists(dest)
    nr = len(files)
    i = 0
    skipped = 0

    print("Starting copying of {0} files to destination folder {1}".format(nr, dest))

    for path in files:
        path = path.rstrip()
        if not os.path.isfile(os.path.join(dest, os.path.basename(path))) or overwrite:
            if os.path.isfile(path):
                shutil.copy(path, dest)
                i += 1
            else:
                "File does not exist: " + str(path)
                skipped += 1
        else:
            skipped += 1
    print("Finished copying of {0} files. Skipped {1} files. ".format(i, skipped))


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        cmd = str(sys.argv[1]).lower()
        if cmd == "writevcfs":
            directory = sys.argv[2]
            dest = str(sys.argv[3])
            write_vcfs_list(directory, dest)
        elif cmd == "copyvcfs":
            assert len(sys.argv) >= 4
            src = sys.argv[2]
            dest = sys.argv[3]
            f = open(src)  # dir
            copy_vcf(f.readlines(), dest)
            f.close()
    else:
        print("Incorrect command syntax. Insufficient args for command '{}'".format(sys.argv[-1]))
