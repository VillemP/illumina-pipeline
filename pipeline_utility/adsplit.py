import sys


def split_ad(col, empty_field):
    """
    Returns the variance percent
    :param col: Unedited column to work on (e.g. "127,124")
    :param empty_field: String to denote empty field
    :return: tuple(int, int, int) reference reads, variant reads, percentage
    """
    ref_reads = 0
    var_reads = 0
    variance = 0.0
    a = col.split(",")
    if len(a) == 1:
        a[0] = empty_field
        a.append(empty_field)
    if a[0] == empty_field or a[1] == empty_field:
        c = 0
    elif int(a[0]) == 0 and int(a[1]) == 0:
        c = 0
    else:
        c = round((float(a[1]) / (float(a[0]) + float(a[1])) * 100), 1)
    ref_reads = a[0]
    var_reads = a[1]
    variance = c
    return ref_reads, var_reads, variance


def add_variance_pct(empty_field=".", pos=-1, table=sys.stdin, output=sys.stdout):
    """
    Takes in a table in which the specified column is AD and splits it into five columns:
    Genotype    Depth   Ref.reads   Var.reads   %var
    :param empty_field: Denotes the string to represent empty fields (default '.')
    :param pos: Which column contains the AD; 0-indexed (default is last or -1)
    Expects The column upstream to be Depth, and two columns upstream to be Genotype
    (will rename the annotators header)
    :type table: Input stream (e.g. a file object or stdin [default])
    :param output: Output stream (e.g. a file object or stdout [default])
    """

    for index, row in enumerate(table):
        seperated = row.split()
        if index == 0:
            # Edit the header
            seperated[pos - 2] = "Genotype"
            seperated[pos - 1] = "Depth"
            seperated[pos] = "Ref.reads"
            seperated.append("Var.reads")
            seperated.append("%var")
            output.write('\t'.join(seperated) + "\n")
        else:
            col = seperated[pos]
            ref_reads, var_reads, variance = split_ad(col, empty_field)
            seperated[pos] = str(ref_reads)
            seperated.append(str(var_reads))
            seperated.append(str(variance))
            output.write('\t'.join(seperated) + "\n")


if __name__ == '__main__':
    args = sys.argv
    if len(args) >= 3:
        print ""
        print "Input file: " + args[1]
        print "Output is written to: " + args[2]
        print ""
        with open(args[2], "w+") as out:
            add_variance_pct(output=out)

        print ""
        print "Finished!"
    else:
        add_variance_pct()
