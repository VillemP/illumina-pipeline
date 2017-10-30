GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"

# TODO: Find these strings from the config. Problematic circular dependancy.
tso_genes = list()


def load_tso_genes(tsoPath):
    with open(tsoPath) as tso:
        data = tso.readlines()
        for line in data:
            # print data[i].rstrip().split("\t")
            tso_genes.append((line.rstrip().split("\t")))


def get_tso_status(gene):
    # type: (gene) -> Gene
    if len(tso_genes) == 0:
        load_tso_genes(gene.panel.config.tsoGenes)
    # get the coverage from the tuple (name, coverage)
    match_coverage = next((gene_cov[1] for gene_cov in tso_genes if gene_cov[0] == gene.name), None)
    if match_coverage is not None:
        gene._Gene__coverage = float(match_coverage)
        return True
    return False


def find_synonyms(gene):
    # TODO: get gene synonyms
    pass


class Gene(object):
    def __init__(self, panel, json):
        """

        :type json: dict
        """
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
        self.on_TSO = get_tso_status(self)

    @property
    def evidence_level(self):
        return self.__confidence

    @property
    def coverage(self):
        return self.__coverage

    def __str__(self):
        pass
