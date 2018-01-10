#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys


def annotate_pos(anno_file, ac, het, hom, pos=21, empty=".", table=sys.stdin, output=sys.stdout):
    ann_set = []
    with open(anno_file, "r") as annotationfile:
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
            var = [variant for variant in ann_set if variant[0:4] == [row[0], row[1], row[3], row[4]]]
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
    args = sys.argv
    if len(args) >= 5:
        print ""
        print "Input file: " + args[6]
        print "Output is written to: " + args[7]
        print ""
        with open(args[3], "w+") as output:
            annotate_pos(args[1], args[2], args[3], args[4], args[5], output=output)

        print ""
        print "Finished!"
    else:
        annotate_pos(args[1], args[2] + ".AC", args[2] + ".HET", args[2] + ".HOM")
