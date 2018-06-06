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

reflines = None


class Interval():
    def __init__(self, chrom, start, stop, symbol):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.symbol = symbol

    def between(self, range):
        if self.chrom == range[0]:
            if self.start <= range[1] and self.stop >= range[2]:
                return True
        return False

    def __str__(self):
        return "\t".join((self.chrom, str(self.start), str(self.stop), str(self.symbol)))

    def compare(self, other):
        return self.symbol == other.symbol and self.chrom == other.chrom


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
# TODO: if the area covered is outside the regular chromosomes then the lexicographic sort will add
# contigs inbetween the chromosomes e.g. chr6 chr6_hapt2 chr7, which is incompatible with the regular UCSC reference
# file order (which has the extra contigs after chrY) and therefore GATK.
def filter_targetfile(hgncHandler, geneslist, targets, genecolumn=3, sort0=0, sort1=1, sort2=2):
    """
    Creates a custom sorted targetfile (default) or refseq from a sample-specific order of gene symbols.
    Default is equivalent to:
    grep $'\t'${gene}$'\.' < ${targets.bed} >> ${sample}.genes.bed
    sort -k3,3V -k4,4n -k5,5n < ${sample}.genes.bed > ${prefix}.genes.sorted.bed

    Also tries to convert the gene name in the target file into a HGNC symbol since internally
    all gene names in the application are HGNC symbols if possible (see gene.py)
    :param geneslist: list of Gene.names (HGNC gene symbols)
    """
    g = None
    final_targets = []

    split_targets = []  # Final edited list
    for i, line in enumerate(targets):
        target_line = []
        if not type(line) is list:
            for element in line.split('\t'):
                target_line.append(element.strip())
        else:
            for element in line:
                target_line.append(element.strip())
        assert len(
            target_line) >= genecolumn, "You input bedfile is malformed! " \
                                        "The bedfile should contain a column with gene symbols " \
                                        "(CHR START STOP GENE)"
        # The last element of the target line is a string GENE_name.EXON_START_int.EXON_END_int
        # Get the gene name from the target
        gene_name = target_line[genecolumn].split(".")[0]  # The gene name is the first element in the column
        if i == 0:
            # If working with a large "unfiltered" refseq,
            # don't report not finding genes that exists in the refseq but not in the covered genes
            g = hgncHandler.find_gene(gene_name, verbose=False)
        else:
            if g is not None:
                if g.name != gene_name:  # Don't go looking for the gene again, it's the same gene but a different exon
                    g = hgncHandler.find_gene(gene_name)  # The target is for a new gene
            else:
                g = hgncHandler.find_gene(gene_name)  # g was set to None (no gene found), therefore look for it
        if g is not None:
            gene_name = g.name  # Ensure the HGNC name for this target, if none found, gene_name is the same as before
        target_line.append(gene_name)  # The final line of the target is now a converted HGNC symbol or
        # the string grepped from the column
        split_targets.append(target_line)  # List of all lines

    for ge in geneslist:
        # The final element is the converted HGNC symbol, also convert the query symbol (ge) to a HGNC name
        g = hgncHandler.find_gene(ge)  # Returns a gene object
        if g is not None:
            gene_name = g.name  # Found a gene
        else:
            gene_name = ge  # Try to find the gene from the targets using the ordered symbol
        gene_targets = [line for line in split_targets if line[-1] == gene_name]
        for t in gene_targets:
            final_targets.append(t)
    s = sorted(final_targets,
               key=lambda final_target: (Version(final_target[sort0]), int(final_target[sort1]),
                                         int(final_target[sort2])))  # sort by chromosome (numeric not lexicographic)
    return s


def readdict(reference_dictionary, verbose=False):
    """
    Reads in the reference dict file and returns the contigs in order.
    :param verbose: Write the output to sys.stderr
    :param reference_dictionary: The dictionary file of the reference fasta (e.g. ucsc.hg19.dict)
    :return: list
    """
    contigs = []
    with open(reference_dictionary) as ref:
        for i, line in enumerate(ref.readlines()):
            if line.startswith("@SQ"):  # contig line
                # contig is "SN:chrM" (2. column)
                contig = line.split("\t")[1].strip().split(":")[1]
                contigs.append(contig)  # create a list of all contigs and the positional order they're in
    if verbose:
        sys.stderr.write("Reference dict contig order: \n")
        sys.stderr.write(str(contigs) + "\n")  # print the reference contig order
    return contigs


def readrefseq(refSeq):
    """
    Converts the refseq lines to a list. Returns it as a 2-dimensional table of lists)
    :param refSeq: list of lines in the refseq file
    :return: list
    """
    allCols = []
    with open(refSeq, "r") as reads:
        for i, line in enumerate(reads):
            if i != 0:  # skip the header
                allCols.append(
                    line.split("\t"))  # adds a column-split line to the list (the final column contains a LF)

    return allCols


def sort_refseq(ref_dict, reflines):
    """
    Sort the refseq file according to chromosome, startpos, endpos based on the input reference sequence dictionary file
    Returns a 2-dimensional sorted table of lists.
    :param ref_dict: Path to the reference dict (e.g. ucsc.hg19.dict)
    :param refseq: List of lines in the refseq file
    :return: list
    """
    contigs = readdict(ref_dict)
    reflines.sort(key=lambda line: (contigs.index(line[2]), int(line[4]), int(line[5])))
    return reflines


def write_targetfile(hgncHandler, geneslist, targetfile, out=sys.stdout):
    with open(targetfile) as targetf:
        lines = targetf.readlines()
        for line in filter_targetfile(hgncHandler, geneslist, lines):
            out.write("\t".join(line) + "\n")


def write_refseq(hgncHandler, geneslist, refseq, refdict=None, out=sys.stdout):
    # global reflines
    reflines = readrefseq(refseq)
    gene_filtered_refseq = filter_targetfile(hgncHandler, geneslist, reflines, genecolumn=12, sort0=2, sort1=4, sort2=5)
    # Definitely make sure our sort catches all sorting orders of refseq.
    # This will sort the refseq according to the reference dictionary file.
    final_sort = sort_refseq(refdict, gene_filtered_refseq)
    for line in final_sort:
        out.write("\t".join(line) + "\n")


def add_symbols_to_interval_summary(targetfile, intervalsummary, out=None, verbose=False):
    interval_list = []
    with open(targetfile) as targetfile:
        lines = targetfile.readlines()
        for line in lines:
            cols = line.split("\t")
            if len(cols) >= 4:
                interval = Interval(cols[0], int(cols[1]) - 1, int(cols[2]) - 1,
                                    cols[3].strip())  # GATK reports intervals as +1
                match = [item for item in interval_list if item.compare(interval)]
                if len(match) >= 1:
                    if match[0].start > interval.start:
                        match[0].start = interval.start
                    if match[0].stop < interval.stop:
                        match[0].stop = interval.stop
                else:
                    interval_list.append(interval)
            else:
                raise IndexError("The interval file contains too few columns or your file is not tab-delimited! "
                                 "The 4th column should contain the symbol for the interval e.g. 'CFTR.exon.3'.")
    if verbose:
        for i in interval_list:
            sys.stderr.write(str(i) + "\n")

    with open(intervalsummary) as interv:
        summary = interv.readlines()
        output = []
        header = ""
        for i, line in enumerate(summary):
            if i != 0:  # skip the header
                cols = line.split("\t")
                chr, range = cols[0].split(":")
                range = range.split("-")  # chr2:25475032-25475256
                range = (chr, int(range[0]) - 1, int(range[1]) - 1)  # convert to tuple (chr1, 25475031, 25475255)
                for interval in interval_list:
                    if interval.between(range):
                        cols.insert(0, interval.symbol)
                output.append("\t".join(cols))
            else:
                header = line
                header = header.split("\t")
                header.insert(0, "Location")
    if out is not None:
        with open(out, "w+") as outfile:
            # Insert the header
            outfile.write("\t".join(header))
            for line in output:
                outfile.write(line)
    return output


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
    targetfile.add_argument("-t", "--targetfile", help="The input targetfile.", type=str, action="store")
    targetfile.add_argument("--dict", help="The reference sequence dictionary file.", type=argparse.FileType('r'))
    targetfile.add_argument("-s", "--hgnc", help="The HGNC gene symbols reference file.", type=str, action="store")
    targetfile.add_argument("-r", "--tsogenes", help="The genes reference file for "
                                                     "genes covered on the sequencing panel", type=str, action="store")
    targetfile.add_argument("-o", "--output", help="The output file to be created", action="store", type=str)
    targetfile.add_argument("--bed",
                            help="The interval file (should also contain the symbols corresponding to each"
                                 " interval.", type=str, action="store")

    intervalsummary = subparsers.add_parser("intervalsummary")
    intervalsummary.add_argument("--interval",
                                 help="Add gene/exon symbols to an interval based summary file created by "
                                      "GATK DepthOfCoverage.")
    intervalsummary.add_argument("--bed",
                                 help="The interval file (should also contain the symbols corresponding to each"
                                      " interval.", type=argparse.FileType('r'))
    intervalsummary.add_argument("--summary", help="The interval_summary file created by GATK that is to be annotated.")
    intervalsummary.add_argument("--out", help="The output file, default is stdout", default="stdout")
    refseq = subparsers.add_parser("refseq")
    refseq.add_argument("-g", "--genes", help="List of genes to be grepped and sorted into a new targetfile",
                        action="store", type=str, nargs="+")
    refseq.add_argument("-t", "--targetfile", help="The input targetfile.", type=str, action="store")
    refseq.add_argument("--dict", help="The reference sequence dictionary file.", type=argparse.FileType('r'))
    refseq.add_argument("-s", "--hgnc", help="The HGNC gene symbols reference file.", type=str, action="store")
    refseq.add_argument("-r", "--tsogenes", help="The genes reference file for "
                                                 "genes covered on the sequencing panel", type=str, action="store")
    refseq.add_argument("-o", "--output", help="The output file to be created", action="store", type=str)
    refseq.add_argument("--bed",
                        help="The interval file (should also contain the symbols corresponding to each"
                             " interval.", type=str, action="store")


    args = parser.parse_args()

    if args.command == "writevcfs":
        write_vcfs_list(args.source, args.output)
    if args.command == "copyvcfs":
        with open(args.source) as f:
            copy_vcf(f.readlines(), args.destination)
    if args.command == "refseq":
        handler = gene.HgncHandler(args.hgnc, args.tsogenes)
        write_refseq(handler, args.genes, args.targetfile, refdict=args.dict.name)
    if args.command == "targetfile":
        handler = gene.HgncHandler(args.hgnc, args.tsogenes)
        write_targetfile(handler, args.genes, args.bed)
    if args.command == "intervalsummary":
        product = add_symbols_to_interval_summary(args.bed, args.summary)
        if args.out == "stdout":
            for line in product:
                sys.stdout.write(line)
        else:
            with open(args.out, "w+") as out:
                for line in product:
                    out.write(line)
