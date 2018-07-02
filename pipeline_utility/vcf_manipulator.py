#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
from collections import defaultdict


def load_annotation_file(annotation_file):
    with open(annotation_file) as f:
        d = create_gene_dict(f)
    return d


# Useful in a future setting skipping tempfiles and subprocessing
def create_gene_dict(lines):
    d = defaultdict(set)
    for line in lines:
        # Skip comment lines
        if not line.startswith('#'):
            try:
                (key, value) = line.split("\t")
                if len(value) > 1:
                    val = value[:-1]  # Get the last column GENE\tANNOTATION
                else:
                    val = "no_annotation"
                v = val.replace(" ", "_").replace(";", "|").replace("=",
                                                                    "_").strip()  # Escape newlines and replace illegal chars
                if key in d:
                    # Only add unique values
                    if v not in d[key]:
                        d[key].append(v)
                    else:
                        # Some debugging
                        # sys.stderr.write(key + str(d[key]) + "\n")
                        pass
                else:
                    d[key] = [v]
            except ValueError as e:
                sys.stderr.write("ANNOTATOR error: {0}\n"
                                 "Your annotation file format might be incorrect (spaces instead of tabs)\n", e.message)
    return d



def annotate(gene_dict, annotation, vcf=sys.stdin, output=sys.stdout, progress=True):
    sys.stderr.write("Starting VCF manipulator with '{0}' annotation...\n".format(annotation))
    # read in the vcf
    headers = 0
    variants = 0

    for i, row in enumerate(vcf):
        if row[0] == "#":
            headers += 1

            if row[0:6] == "#CHROM":
                # final header line, write our own header
                headerline = "##INFO=<ID={0},Number=.,Type=String,Description=\"Custom gene-based annotation\">\n".format(
                    annotation.strip())
                output.write(headerline)
            output.write(row)
        else:
            variant = row.split()
            # The eight column is the INFO tag, containing the annotations
            try:
                info = variant[7].split(';')
            except IndexError as e:
                sys.stderr.write("\t".join(variant) + "\n")
                raise e
            refgene = filter(lambda x: 'Gene.refGene=' in x, info)
            genelist = refgene[0].split('=')
            gene = genelist[1]
            if gene in gene_dict:
                info.append(annotation + "=" + ",".join(gene_dict[gene]))  # escape any newlines
            else:
                info.append(annotation + "=no_annotation")
            variant[7] = ';'.join(info)
            output.write("\t".join(variant) + "\n")
            variants += 1
        if progress:
            if i % 4000 == 0:
                if i != 0:
                    sys.stderr.write("VCF manipulator done with {0} lines ({1})...\n".format(i, annotation))
    sys.stderr.write("Finished VCF manipulator with '{0}' annotation...\n".format(annotation))


if __name__ == "__main__":
    args = sys.argv
    if len(args) >= 5:
        print ""
        print "Annotation file: " + args[1]
        print "Annotation name: " + args[2]
        print "Input vcf: " + args[3]
        print "Output is written to: " + args[4]
        print ""

        annotation_file = args[1]
        gene_dict = load_annotation_file(annotation_file)

        with open(args[4], "w+") as out:
            with open(args[3]) as infile:
                annotate(gene_dict, args[2], infile, out)
    elif len(args) >= 3:
        annotation_file = args[1]
        gene_dict = load_annotation_file(annotation_file)
        annotate(gene_dict, args[2])
    else:
        print("Too few arguments! Missing annotation file.")
