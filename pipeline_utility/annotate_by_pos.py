#!/usr/bin/python
# -*- coding: UTF-8 -*-
import argparse
import sys


def annotate_pos(anno_file, ac, het, hom, chrpos=0, g_pos=1, ref_pos=3, alt_pos=4,
                 pos=21, empty=".", table=sys.stdin, output=sys.stdout):
    # sys.stderr.write("I'm here!\n")
    ann_list = []
    ann_dict = {}
    with open(anno_file, "r") as annotationfile:
        for line in annotationfile:
            split = line.split()
            key = "-".join(split[0:4])
            ann_dict[key] = split
            # ann_list.append(split)
    annotation_hash = set([v for v in ann_dict.keys()])  # Create a hashmap of the list

    for index, line in enumerate(table):
        row = line.split()
        if index == 0:
            row.insert((int(pos)), ac)
            row.insert((int(pos) + 1), het)
            row.insert((int(pos) + 2), hom)
            output.write("\t".join(row) + "\n")
        else:

            current_variant = "-".join([row[chrpos], row[g_pos], row[ref_pos], row[alt_pos]])
            # var = next((v for v in ann_set if v[0:4] == current_variant), None)
            if current_variant in annotation_hash:

                var = ann_dict[current_variant]
                ac_val = var[4]
                het_val = var[5]
                hom_val = var[6]
                row.insert((int(pos)), str(ac_val))
                row.insert((int(pos) + 1), str(het_val))
                row.insert((int(pos) + 2), str(hom_val))
            else:
                row.insert((int(pos)), empty)
                row.insert((int(pos) + 1), empty)
                row.insert((int(pos) + 2), empty)

            output.write("\t".join(row) + "\n")


if __name__ == '__main__':
    # positional arguments e.g. 0,1,2,3 correlates to CHR,POS,REF,ALT in the indexes 0,1,2,3
    # Insertion position is 21 by default.
    parser = argparse.ArgumentParser(prog="Annotation pipeline command-line file tool for local DB annotation.")
    subparsers = parser.add_subparsers(title="commands", dest="command")

    db_freq = subparsers.add_parser("db_freq", help="Write the frequency of variants in a local DB")
    db_freq.add_argument("-a", "--annotation", action="store", type=str)
    db_freq.add_argument("-i", "--input", help="The input file, default stdin", action="store", type=str,
                         default="-")
    db_freq.add_argument("-o", "--output", help="The output file, default stdout", action="store", type=str,
                         default="-")
    db_freq.add_argument("--totalsamples", type=int, action="store")
    db_freq.add_argument("--positions",
                         help="The representation of CHROM,POS,REF,ALT index positions in the DB file (default is '0,1,3,4')",
                         action="store", type=str, default="0,1,3,4")
    db_freq.add_argument("--insertpos", help="The index of the column where the DB columns are inserted. Default is 21",
                         type=int, action="store", default=21)

    args = parser.parse_args()

    positions = args.positions.split(",")
    positions = map(int, positions)
    assert len(positions) == 4, "The position string was not 4 elements long. (Check your seperator)" \
                                "All positions need to be described in the positions string (e.g. 0,1,3,4)"

    if args.command == "db_freq":
        if args.input is not "-" and args.output is not "-":
            with open(args.input) as infile:
                with open(args.output, "w+") as outfile:
                    annotate_pos(args.annotation, ac=str(args.totalsamples) + "AC", het=str(args.totalsamples) + "HET",
                                 hom=str(args.totalsamples) + "HOM", chrpos=positions[0], g_pos=positions[1],
                                 ref_pos=positions[2], alt_pos=positions[3], table=infile, output=outfile,
                                 pos=args.insertpos)
        elif args.input is not "-":
            with open(args.input) as infile:
                annotate_pos(args.annotation, ac=str(args.totalsamples) + "AC", het=str(args.totalsamples) + "HET",
                             hom=str(args.totalsamples) + "HOM", chrpos=positions[0], g_pos=positions[1],
                             ref_pos=positions[2], alt_pos=positions[3], table=infile)
        elif args.output is not "-":
            with open(args.output) as outfile:
                annotate_pos(args.annotation, ac=str(args.totalsamples) + "AC", het=str(args.totalsamples) + "HET",
                             hom=str(args.totalsamples) + "HOM", chrpos=positions[0], g_pos=positions[1],
                             ref_pos=positions[2], alt_pos=positions[3], output=outfile)
        else:
            annotate_pos(args.annotation, ac=str(args.totalsamples) + "AC", het=str(args.totalsamples) + "HET",
                         hom=str(args.totalsamples) + "HOM", chrpos=positions[0], g_pos=positions[1],
                         ref_pos=positions[2], alt_pos=positions[3])

    else:
        sys.stderr.write("Unknown command {0}\n".format(args.command))
