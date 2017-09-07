import os, sys, shutil

#MiniSeq: /media/kasutaja/ge_ngs/smb-share:server=srvfail,share=ge_ngs/MiniSeq

commands = ['WriteVcfs', 'CopyVcfs']
duplicates = 0


def find_vcfs(directory):
    assert os.path.exists(directory)
    duplicates = 0
    unique_files = dict([])
    for (dirpath, dirnames, files) in os.walk(directory):
        for f in files:
            if f.endswith('.vcf'):
                if f not in unique_files:
                    unique_files[f] = dirpath
                else:
                    duplicates += 1

    return unique_files


def write_vcfs_list(directory, output):
    try:

        file_paths = find_vcfs(str(directory))

        f = open(str(directory + output), "w+")
        for vcf_addr in file_paths.iteritems():
            full_path = os.path.join(vcf_addr[1] + vcf_addr[0])
            print (full_path)
            f.write(full_path + "\n")
        print("Total files: {}".format(len(file_paths)))
        print("Duplicates: {}".format(duplicates))

        f.close()

        return file_paths
    except (NameError, KeyboardInterrupt, SyntaxError, AssertionError):
        if NameError or AssertionError:
            print("Unable to find {}".format(directory))
            raise
        if KeyboardInterrupt:
            print("Quitting.")


def copy_vcf(files, dest):
    assert os.path.exists(dest)
    nr = len(files)
    i = 0
    print("Starting copying of {0} files to destination folder {1}".format(nr, dest))

    for path in files:
        path = path.rstrip()
        assert os.path.isfile(path), "File does not exist: " + str(path)
        shutil.copy(path, dest)
        print("Finished copying file {0}".format(path))
        i+= 1
        print("Files remaining: {}".format(nr - i))

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
            f = open(src) #dir
            copy_vcf(f.readlines(), dest)
            f.close()
    else:
        print("Incorrect command syntax. Insufficient args for command {}".format(sys.argv[1]))
        raise ValueError