#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys


def annotate_pos(anno_file, ac, het, hom, chrpos=0, g_pos=1, ref_pos=3, alt_pos=4,
                 pos=21, empty=".", table=sys.stdin, output=sys.stdout):
    ann_set = []
    with open(anno_file, "r") as annotationfile:
        sys.stderr.write(anno_file + "\n")
        for line in annotationfile:
            ann_set.append(line.split())

    for index, line in enumerate(table):
        row = line.split()
        if index == 0:
            row.insert((int(pos)), ac)
            row.insert((int(pos) + 1), het)
            row.insert((int(pos) + 2), hom)
            output.write("\t".join(row) + "\n")
        else:
            var = [variant for variant in ann_set if variant[0:4] == [row[chrpos], row[g_pos], row[ref_pos],
                                                                      row[alt_pos]]]
            sys.stderr.write(str(var) + "\n")
            if len(var) > 0:
                ac_val = var[0][4]
                het_val = var[0][5]
                hom_val = var[0][6]
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
    args = sys.argv
    if len(args) >= 6:
        print ""
        print "Input file: " + args[6]
        print "Output is written to: " + args[7]
        print ""
        with open(args[3], "w+") as output:
            annotate_pos(args[1], args[2], args[3], args[4], args[5], output=output)

        print ""
        print "Finished!"
    else:
        if len(args) == 5:
            positions = args[3].split(",")
            positions = map(int, positions)
            if len(positions) == 4:
                # We need all the positions to map the variant (create a unique representation of it)
                sys.stderr.write(str(positions) + "\n")
                annotate_pos(args[1], args[2] + ".AC", args[2] + ".HET", args[2] + ".HOM", chrpos=positions[0],
                             g_pos=positions[1], ref_pos=positions[2], alt_pos=positions[3],
                             pos=args[4])
        else:
            annotate_pos(args[1], args[2] + ".AC", args[2] + ".HET", args[2] + ".HOM")
