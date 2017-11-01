import sys

GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"

# TODO: Find these strings from the config. Problematic circular dependancy.
tso_genes = list()
tso_genes_tuple = list()


def load_tso_genes(tsoPath):
    try:
        with open(tsoPath) as tso:
            data = tso.readlines()
            for line in data:
                # print data[i].rstrip().split("\t")
                gene_cov = line.rstrip().split("\t")
                gene = Gene(panel=None, json=None, on_TSO=True)
                gene.name = gene_cov[0]
                gene.__coverage = gene_cov[1]
                tso_genes.append(gene)
    except IOError as error:
        sys.stderr.write("PIPELINE ERROR: {0}\n"
                         "You are missing the coverage file in the location {1}\n"
                         "The file should contain the gene name and coverage, tab seperated (GENE_NAME\\tCOVERAGE\\n)"
                         .format(error.message, tsoPath))
    except IndexError as error:
        sys.stderr.write("PIPELINE ERROR: {0}\n Is your file format correct?"
                         "The file should contain the gene name and coverage, tab seperated (GENE_NAME\\tCOVERAGE\\n)"
                         .format(error.message))


def get_tso_status(gene):
    # type: (gene) -> Gene
    if len(tso_genes) == 0:
        # Load TSO genes if currently not loaded.
        load_tso_genes(gene.panel.config.tsoGenes)
    # get the coverage from the tuple (name, coverage)
    match_coverage = next((gene.coverage for gene in tso_genes if gene.name == gene.name), None)
    if match_coverage is not None:
        gene._Gene__coverage = float(match_coverage)
        return True
    return False


def find_gene(name):
    genes = [g for g in tso_genes if g.name == name]
    if len(genes) > 0:
        if len(genes) == 1:
            return genes[0]
        else:
            sys.stderr.write("PIPELINE ERROR: {0}\nThe gene {1} had several matches among the covered genes!\n"
                             .format(name))
            # TODO: What happens if there are several matches?
    else:
        sys.stderr.write("PIPELINE ERROR: \nThe gene {1} was not found among the covered genes.\n".format(name))


def find_synonyms(gene):
    # TODO: get gene synonyms
    pass


class Gene(object):
    def __init__(self, panel, json=None, on_TSO=None):
        # type: (GenePanel, dict, bool) -> Gene

        if json is not None:
            self.ensemblegeneids = json['EnsembleGeneIds']
            self.name = json['GeneSymbol']
            self.synonyms = []
            self.panel = panel
            self.panel_name = panel.name
            self.__confidence = json['LevelOfConfidence']
            self.modeofinheritance = json['ModeOfInheritance']
            self.modeofpathogenicity = json['ModeOfPathogenicity']
            self.penetrance = json['Penetrance']
            self.phenotypes = json['Phenotypes']
            self.raw_json = json
            self.__coverage = 0.0
        else:
            self.ensemblegeneids = None
            self.name = None
            self.synonyms = []
            self.panel = None
            self.__confidence = 1
            self.modeofinheritance = None
            self.modeofpathogenicity = None
            self.penetrance = None
            self.phenotypes = None
            self.raw_json = None
            self.__coverage = 0
        if on_TSO is None:
            self.on_TSO = get_tso_status(self)
        else:
            self.on_TSO = on_TSO

    @property
    def evidence_level(self):
        return self.__confidence

    @property
    def coverage(self):
        return self.__coverage

    def __str__(self):
        return self.name
