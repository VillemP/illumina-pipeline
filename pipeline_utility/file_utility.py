import os

try:
    from scandir import walk
except ImportError:
    from os import walk

import shutil
import sys

# MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq

commands = ['WriteVcfs', 'CopyVcfs']
duplicates = 0


def find_file(main_dir, filename):
    """
    Will find the first file with the exact name match in a dir tree. Otherwise returns None.
    :param main_dir: The upmost directory to start the search in, will walk through subdirectories
    :param filename: The filename to be searched.
    :return: tuple(filename, directory)
    """
    assert os.path.exists(main_dir), "Path {} does not exist.".format(main_dir)

    for (dirpath, dirnames, files) in walk(main_dir):
        for name in files:
            if name == filename:
                return name, os.path.join(dirpath, name)
    return None, None


def find_filetype(dir, filetype):
    """
    Will find all files of a certain type (e.g. .vcf or .bam files) in a directory.
    Method will enter every subdirectory. Can look for only a single filetype at a time.
    :param dir: String of directory to walk.
    :param filetype: String of filetype to search for (e.g. .vcf or .bam)
    :return: list of tuples of file name and file directory
    """
    assert os.path.exists(dir), "Path {} does not exist.".format(dir)
    duplicates = 0
    # unique_files = dict([])
    unique_files = list(())
    for (dirpath, dirnames, files) in walk(dir):
        for name in files:
            if name.endswith(filetype):
                if name not in unique_files:
                    unique_files.append((name, os.path.join(dirpath, name)))
                else:
                    duplicates += 1

    return unique_files


def write_filelist(dir, file_name, file_paths):
    try:
        with open(os.path.join(dir, file_name), "wb+") as f:
            all_paths = list()
            file_paths.sort(key=lambda pair: pair[0])
            for file_path_pair in file_paths:
                full_path = file_path_pair[1]
                all_paths.append(full_path)
                # print (full_path)
                f.write(full_path + "\n")
            print("Writing filelist to {}".format(os.path.join(dir, file_name)))
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
            prefix_clean = prefix[0].rsplit('.')
            if len(prefix_clean) == 2:
                prefix_clean = prefix_clean[0]
                prefixes.write(prefix_clean + "\n")
                clean_prefixes.append(prefix_clean)
            else:
                pass
                # This file might be an intermediary vcf file

    return clean_prefixes


def write_vcfs_list(dir, output):
    return write_filelist(dir, output, find_vcfs(dir))


def write_bams_list(dir, output):
    return write_filelist(dir, output, find_bams(dir))


def copy_vcf(files, dest, overwrite=False):
    assert os.path.exists(dest), "Destination {} does not exist!".format(dest)
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
                print "File does not exist: " + str(path)
                skipped += 1
        else:
            skipped += 1
    print("Finished copying of {0} files. Skipped {1} files. ".format(i, skipped))
    return i, skipped


def file_len(fname):
    if os.path.exists(fname):
        with open(fname) as f:
            i, l = 0, 0
            for i, l in enumerate(f):
                pass
        return i + 1
    return None


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
