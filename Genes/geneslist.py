import os

# import numpy

geneslist_dir = ""


def read_genes():
    genes = list()
    assert os.path.exists(geneslist_dir)
    with open(geneslist_dir, "r") as geneslist:
        lines = geneslist.readlines()
        for line in lines:
            genes.append(line.rstrip())
        return genes


if __name__ == "__main__":
    genes = read_genes()
