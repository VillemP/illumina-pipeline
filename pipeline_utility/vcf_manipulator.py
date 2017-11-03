import sys


def load_annotation_file(annotation_file):
    with open(annotation_file) as f:
        d = create_gene_dict(f)
    return d


# Useful in a future setting skipping tempfiles and subprocessing
def create_gene_dict(lines):
    d = {}
    for line in lines:
        try:
            (key, value) = line.split("\t")
        except ValueError as e:
            sys.stderr.write("ANNOTATOR error: {0}\n"
                             "Your annotation file format might be incorrect (spaces instead of tabs)\n", e.message)
        if len(value) > 1:
            val = value[:-1]
        else:
            val = "no_annotation"
        v = val.replace(" ", "_")
        if key in d:
            d[key].append(v)
        else:
            d[key] = [v]
    return d


def annotate(gene_dict, annotation, vcf=sys.stdin, output=sys.stdout):
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
            if gene in gene_dict:
                info.append(annotation + "=" + ",".join(gene_dict[gene]))
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
