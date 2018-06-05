#!/usr/bin/python
# -*- coding: UTF-8 -*-
import argparse
import errno
import os
import sys

import TrusightOne.gene as gene


class HgncConverterTool:
    def __init__(self, hgnc_path, tso, hgnc_genes=None, hgncHandler=None):
        if hgncHandler is not None:
            self.hgncHandler = hgncHandler
        else:
            self.hgncHandler = gene.HgncHandler(hgnc_path, tsoPath=tso, verbosity=gene.EXCEPTIONS)
        # Load our HGNC terms if not already loaded in the module
        if hgnc_genes is not None:
            self.hgnc_genes = hgnc_genes
        elif len(self.hgncHandler.hgnc_genes) == 0:
            self.hgnc_genes = self.load(hgnc_path)
        else:
            self.hgnc_genes = self.hgncHandler.hgnc_genes

    def load(self, hgnc_path):
        if self.hgncHandler.hgnc_genes is not None:
            if len(self.hgncHandler.hgnc_genes) == 0:
                self.hgncHandler.load_hgnc_genes(hgnc_path)

        return self.hgncHandler.hgnc_genes

    def convert(self, vcf=sys.stdin, output=sys.stdout, progress=False):
        sys.stderr.write("Starting gene symbol conversion to HGNC terms for input {0}...\n".format(str(vcf)))
        # read in the vcf
        headers = 0
        variants = 0
        edited = {}
        edited_variants = 0
        last_changed = 0
        for i, row in enumerate(vcf):
            if row[0] == "#":
                headers += 1
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
                # Some genes are paired with their synonyms: TSIX\x3bXIST
                # genelist[1] is the raw representation
                symbol = genelist[1].split("\\x3b")
                # symbol = [genelist[1], None]
                hgnc = self.hgncHandler.find_gene(symbol[0], verbose=False)
                if hgnc is not None:
                    if hgnc.name != genelist[1]:
                        for field_idx, field in enumerate(info):
                            new_field = field
                            if symbol[0] in field:
                                edited[genelist[1]] = hgnc.name  # Log
                                if i != last_changed:
                                    edited_variants += 1
                                last_changed = i
                                new_field = field.replace(genelist[1], hgnc.name)
                            info[field_idx] = new_field

                variant[7] = ';'.join(info)
                try:
                    output.write("\t".join(variant) + "\n")
                except IOError as e:
                    if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
                        # Stop loop on "Invalid pipe" or "Invalid argument".
                        # No sense in continuing with broken pipe.
                        break
                    else:
                        # Raise any other error.
                        raise
                variants += 1
            if progress:
                if i % 400 == 0:
                    if i != 0:
                        sys.stderr.write("HGNC converter done with {0} lines...\n".format(i))
        sys.stderr.write("Finished the HGNC conversion with {0} genes replaced with a HGNC term in {1} variants."
                         "\nOriginal:HGNC {2}".format(len(edited), edited_variants, edited))
        #sys.stderr.flush()

    def convert_gene(self, symbol):
        if not self.hgncHandler.is_hgnc(symbol):
            hgnc = self.hgncHandler.synonyms_to_hgnc(symbol)
            if hgnc is not None:
                return hgnc
            else:
                return symbol  # Don't replace anything

    def __str__(self):
        return "Converter with {0} genes.".format(len(self.hgnc_genes))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Annotation pipeline VCF-variant gene symbol to HGNC converter. "
                                          "Requires a local representation of HGNC.")
    parser.add_argument("-l", "--hgnc", help="The path to the HGNC local representation.", action="store", type=str)
    parser.add_argument("-t", "--tso", help="The path to the TSO covered genes local representation.", action="store",
                        type=str)
    parser.add_argument("-i", "--input", help="The input file, default stdin", action="store", type=str,
                        default="-")
    parser.add_argument("-o", "--output", help="The output file, default stdout", action="store", type=str,
                        default="-")
    parser.add_argument("-p", "--progress", help="Reports the progress of conversion every 1000 lines.",
                        action="store_true", default=False)

    args = parser.parse_args()
    if not os.path.exists(args.hgnc):
        parser.error("Could not find the local representation of HGNC at {0}".format(args.hgnc))
    if args.input == "-":
        inp = sys.stdin
    if args.output == "-":
        out = sys.stdout
    sys.stderr.write(str(args) + "\n")
    handler = HgncConverterTool(args.hgnc, tso=args.tso)

    if args.input != "-" and args.output != "-":
        with open(args.input) as infile:
            with open(args.output, "w+") as outfile:
                handler.convert(vcf=infile, output=outfile, progress=args.progress)
    elif args.input == "-" and args.output != "-":
        with open(args.output) as outfile:
            handler.convert(output=outfile, progress=args.progress)
    elif args.input != "-" and args.output == "-":
        with open(args.input) as infile:
            handler.convert(vcf=infile, progress=args.progress)
    else:
        handler.convert(vcf=sys.stdin, output=sys.stdout, progress=args.progress)
