#!/usr/bin/python
# -*- coding: UTF-8 -*-
import argparse
import os
import sys

import TrusightOne.gene as gene


class HgncConverterTool:
    def __init__(self, hgnc_path, hgnc_genes=None, hgncHandler=None):
        if hgncHandler is not None:
            self.hgncHandler = hgncHandler
        else:
            self.hgncHandler = gene.HgncHandler(hgnc_path, None)
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
                symbol = genelist[1]

                hgnc = self.hgncHandler.find_gene(symbol)
                if hgnc is not None:
                    if hgnc.name != symbol:
                        for field in info:
                            if symbol in field:
                                edited[symbol] = hgnc.name  # Log
                                edited_variants += 1
                                field.replace(symbol, hgnc)
                    else:
                        pass  # Don't replace anything
                variant[7] = ';'.join(info)
                # print("\t".join(variant) + "\n")
                output.write("\t".join(variant) + "\n")
                variants += 1
            if progress:
                if i % 400 == 0:
                    if i != 0:
                        sys.stderr.write("HGNC converter done with {0} lines...\n".format(i))
        sys.stderr.write("Finished the HGNC conversion with {0} genes replaced with a HGNC term in {1} variants."
                         "\nOriginal:HGNC {2}".format(len(edited), edited_variants, edited))

        sys.stderr.flush()

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
    handler = HgncConverterTool(args.hgnc)

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
        handler.convert(progress=args.progress)
