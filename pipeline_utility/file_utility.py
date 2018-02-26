import argparse
import os
from distutils.version import LooseVersion as Version

import TrusightOne.gene as gene

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk

import shutil
import sys

commands = ['WriteVcfs', 'CopyVcfs']
duplicates = 0


def find_file(main_dir, filename):
    """
    Will find the first file with the exact name match in a dir tree. Otherwise returns None.
    :param main_dir: The upmost directory to start the search in, will walk through subdirectories
    :param filename: The filename to be searched.
    :return: tuple(filename, full path)
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


def find_prefixes(dir, extension):
    """
    Finds a list of prefixes with a given file extension e.g. E00001.vcf.gz --> E00001 in a directory

    :param dir: Directory to search in.
    :param extension: File extension string, must not include the seperator '.' ("vcf" or "vcf.gz" not ".vcf")
    """
    samples = find_filetype(dir, extension)
    clean_prefixes = list()
    for prefix in samples:
        prefix_clean = prefix[0].rsplit(".")
        # We only want exact matches, so E00001.genome.vcf.gz would be skipped
        # This tiny math ensures that long extensions don't confuse us.
        if len(prefix_clean) < 2 + len(extension.split(".")):
            prefix_clean = prefix_clean[0]
            clean_prefixes.append(prefix_clean)
    return clean_prefixes


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


def rename_file_idx(name_path, idx):
    basename = os.path.basename(name_path).split(".")
    basename[0] = basename[0] + "_re" + str(idx)
    new_name = ".".join(basename)
    return new_name

def copy_vcf(files, dest, overwrite=False):
    try:
        os.makedirs(dest)
    except OSError:
        if not os.path.isdir(dest):
            raise
    nr = len(files)
    unique_copied = []
    renamed = []
    not_found = []

    print("Starting copying of {0} files to destination folder {1}".format(nr, dest))

    for path in files:
        path = path.rstrip()
        if not os.path.isfile(os.path.join(dest, os.path.basename(path))) or overwrite:
            if os.path.isfile(path):
                shutil.copy(path, dest)
                fname = os.path.join(dest, os.path.basename(path))
                unique_copied.append(fname)
            else:
                print "File does not exist: " + str(path)
                not_found.append(str(path))
        else:
            re_idx = 1

            while os.path.isfile(os.path.join(dest, rename_file_idx(path, re_idx))):
                re_idx += 1
            fname = rename_file_idx(path, re_idx)
            shutil.copy(path, os.path.join(dest, fname))

            renamed.append(os.path.join(dest, fname))
    print("Finished copying of {0} files. Renamed {1} files. ".format(len(unique_copied), len(renamed)))
    return unique_copied, renamed, not_found


def file_len(fname):
    if os.path.exists(fname):
        with open(fname) as f:
            i, l = 0, 0
            for i, l in enumerate(f):
                pass
        return i + 1
    return None


def count_unique_names(infile, col, seperator="\t"):
    names = []
    if os.path.exists(infile):
        with open(infile) as f:
            for i, l in enumerate(f):
                cols = l.split(seperator)
                if not col > len(cols):
                    names.append(cols[col])
                else:
                    raise IndexError("Is your seperator correct? "
                                     "There weren't enough columns "
                                     "after splitting the line #{0}\n{1}!".format(i, l))
    unique = set(names)
    return len(unique)


# TODO: convert sort args to **nargs and sort as many positions as required
def filter_targetfile(geneslist, targetfile, genecolumn=3, sort0=0, sort1=1, sort2=2):
    """
    Creates a custom sorted targetfile (default) or refseq from a sample-specific order of gene symbols.
    Default is equivalent to:
    grep $'\t'${gene}$'\.' < ${targets.bed} >> ${sample}.genes.bed
    sort -k1,1V -k2,2n -k3,3n < ${sample}.genes.bed > ${prefix}.genes.sorted.bed

    Also tries to convert the gene name in the target file into a HGNC symbol since internally
    all gene names in the application are HGNC symbols if possible (see gene.py)
    :param geneslist: list of Gene.names (HGNC gene symbols)
    """
    g = None
    final_targets = []
    # Convert the targetfile to a list with stripped values
    with open(targetfile) as targetfile:
        targets = targetfile.readlines()
    split_targets = []  # Final edited list
    for i, line in enumerate(targets):
        target_line = []
        for element in line.split('\t'):
            target_line.append(element.strip())
        # The last element of the target line is a string GENE_name.EXON_START_int.EXON_END_int
        # Get the gene name from the target
        gene_name = target_line[genecolumn].split(".")[0]  # The gene name is the first element in both BED and refseq
        if i == 0:
            g = gene.find_gene(gene_name)
        else:
            if g is not None:
                if g.name != gene_name:  # Don't go looking for the gene again, it's the same gene but a different exon
                    g = gene.find_gene(gene_name)  # The target is for a new gene
            else:
                g = gene.find_gene(gene_name)  # g was set to None (no gene found), therefore look for it
        if g is not None:
            gene_name = g.name  # Ensure the HGNC name for this target, if none found, gene_name is the same as before
        target_line.append(gene_name)  # The final line of the target is now a converted HGNC symbol or
        # the string grepped from the column
        split_targets.append(target_line)  # List of all lines

    for ge in geneslist:
        # The final element is the converted HGNC symbol, also convert the query symbol (ge) to a HGNC name
        g = gene.find_gene(ge)  # Returns a gene object
        if g is not None:
            gene_name = g.name  # Found a gene
        else:
            gene_name = ge  # Try to find the gene from the targets using the ordered symbol
        gene_targets = [line for line in split_targets if line[-1] == gene_name]
        for t in gene_targets:
            final_targets.append(t)
    s = sorted(final_targets,
               key=lambda final_target: Version(final_target[sort0]))  # sort by chromosome (numeric not lexicographic)
    s = sorted(s, key=lambda final_target: final_target[sort1])  # sort by start index
    s = sorted(s, key=lambda final_target: final_target[sort2])  # sort by end
    return s


def write_targetfile(geneslist, targetfile, out=sys.stdout):
    for line in filter_targetfile(geneslist, targetfile):
        out.write("\t".join(line) + "\n")


def write_refseq(geneslist, refseq, out=sys.stdout):
    for line in filter_targetfile(geneslist, refseq, genecolumn=12, sort0=2, sort1=4, sort2=5):
        out.write("\t".join(line) + "\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="Annotation pipeline command-line file tool.")
    subparsers = parser.add_subparsers(title="commands", dest="command")

    writevcfs = subparsers.add_parser("writevcfs", help="Find all the VCFs in the source directory "
                                                        "and write their paths to a destination file")
    writevcfs.add_argument("-s", "--source",
                           help="Find all the VCFs in the source directory and write their paths to a destination file",
                           action="store", type=str)
    writevcfs.add_argument("-o", "--output",
                           help="Find all the VCFs in a directory and copy them to a destination directory",
                           action="store", type=str)

    copyvcfs = subparsers.add_parser("copyvcfs", help="This tool copies all the VCFs listed in a file into a directory")
    copyvcfs.add_argument("-s", "--source",
                          help="Source file of paths to be copied",
                          action="store", type=str)
    copyvcfs.add_argument("-d", "--destination",
                          help="Destination directory for the VCFS to be copied into.",
                          action="store", type=str)

    targetfile = subparsers.add_parser("targetfile")
    targetfile.add_argument("-g", "--genes", help="List of genes to be grepped and sorted into a new targetfile",
                            action="store", type=str, nargs="+")
    targetfile.add_argument("-t", "--targetfile", help="The input targetfile.", type=argparse.FileType('r'))
    targetfile.add_argument("-s", "--hgnc", help="The HGNC gene symbols reference file.", type=str, action="store")
    targetfile.add_argument("-r", "--tsogenes", help="The genes reference file for "
                                                     "genes covered on the sequencing panel", type=str, action="store")
    targetfile.add_argument("-o", "--output", help="The output file to be created", action="store", type=str)
    args = parser.parse_args()

    if args.command == "writevcfs":
        write_vcfs_list(args.source, args.output)
    if args.command == "copyvcfs":
        with open(args.source) as f:
            copy_vcf(f.readlines(), args.destination)
    if args.command == "targetfile":
        gene.load_hgnc_genes(args.hgnc)
        gene.load_tso_genes(args.tsogenes)
        # write_targetfile(args.genes, args.targetfile)
        write_refseq(args.genes, args.targetfile.name)
