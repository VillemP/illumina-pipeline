import sys


def annotate(annotation_file, annotation, vcf=sys.stdin, output=sys.stdout):
    # make dictionary from annotation file

    d = {}
    with open(annotation_file) as f:
        for line in f:
            (key, value) = line.split("\t")
            if len(value) > 1:
                val = value[:-1]
            else:
                val = "no_annotation"
            v = val.replace(" ", "_")
            if key in d:
                d[key].append(v)
            else:
                d[key] = [v]

    # read in the vcf
    headers = 0
    variants = 0

    for row in vcf:
        if row[0] == "#":
            headers += 1

            if row[0:6] == "#CHROM":
                # final header line, write our own header
                headerline = "##INFO=<ID={0},Number=.,Type=String,Description=\"Custom gene-based annotation\">\n".format(
                    annotation)
                output.write(headerline)
            output.write(row)
        else:
            variant = row.split()
            info = variant[7].split(';')
            refgene = filter(lambda x: 'Gene.refGene=' in x, info)
            genelist = refgene[0].split('=')
            gene = genelist[1]
            if gene in d:
                info.append(annotation + "=" + ",".join(d[gene]))
            else:
                info.append(annotation + "=no_annotation")
            variant[7] = ';'.join(info)
            output.write("\t".join(variant) + "\n")
            variants += 1


if __name__ == "__main__":
    args = sys.argv
    if len(args) >= 5:
        print ""
        print "Annotation file: " + args[1]
        print "Annotation name: " + args[2]
        print "Input vcf: " + args[3]
        print "Output is written to: " + args[4]
        print ""

        with open(args[4], "w+") as out:
            with open(args[3]) as infile:
                annotate(args[1], args[2], infile, out)
    else:
        annotate(args[1], args[2])
