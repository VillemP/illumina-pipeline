"""
Simple gene-based annotation script, that adds info to vfc INFO field.

Usage: python Vcf-manipulation.py input.vcf annotation.file annotation.name output.vcf

Input vcf: previously annovar annotatated vcf file containing Gene.refGene= in INFO field.

Annotation file:	GENE	INFO 	(2 cols, without header)

You can prepare input with oneliner similar to this depending on source file:

cat genes_to_diseases.txt | cut -f2,3 | tail -n +2 > gene.disease.txt
cat geneMap2.txt | cut -f6,12 -d"|" | awk -F'|' '{ split($1,a,", "); for (i in a) print a[i],"\t",$2}' | sed 's/ \t /\t/' | tr ";" "," > gene.omim_disease_name.synonyms.txt

Input file can have same gene on many rows, all info is added

You can use this annotator many times in a row with different annotation files

Caution! This annotation relies on matching gene names in files, be aware of synonyms.
"""

import sys


def annotate(vcf, annotation_file, annotation, output=sys.stdout):
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

    for line in vcf:
        if line[0] == "#":
            headers += 1

            if line[0:6] == "#CHROM":
                # final header line, write our own header
                headerline = "##INFO=<ID={0},Number=.,Type=String,Description=\"Custom gene-based annotation\">\n".format(
                    annotation)
                output.write(headerline)
            output.write(line)
        else:
            variant = line.split()
            info = variant[7].split(';')
            refgene = "".join(filter(lambda x: 'Gene.refGene=' in x, info))
            genelist = refgene.split('=')
            gene = genelist[1]
            if gene in d:
                info.append(annotation + "=" + ",".join(d[gene]))
            else:
                info.append(annotation + "=no_annotation")
            variant[7] = ';'.join(info)
            output.write(" ".join(variant))
            variants += 1


if __name__ == "__main__":
    args = sys.argv
    print ""
    print "Input vcf: " + sys.argv[1]
    print "Annotation file: " + sys.argv[2]
    print "Annotation name: " + sys.argv[3]
    print "Output is written to: " + sys.argv[4]
    print ""

    with open(sys.argv[4], "w+") as out:
        with open(sys.argv[1]) as infile:
            annotate(infile, args[2], args[3], out)
