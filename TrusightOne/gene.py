import sys

GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"

tso_genes = list()
hgnc_genes = {}
synonym_hgnc = {}


def load_hgnc_genes(hgncPath):
    try:
        with open(hgncPath) as hgnc:
            data = hgnc.readlines()
            for line in data:
                symbols = line.split("\t")
                # Match every HGNC symbol to a list of synonyms
                hgnc_genes[symbols[0]] = symbols[1].strip().split(',') + symbols[2].strip().split(',')
                # Match every synonym to it's HGNC symbol to enable quicker searching
        for keypair in hgnc_genes.iteritems():
            for symbol in keypair[1]:  # values
                if symbol is not None and symbol != "":
                    synonym_hgnc[symbol] = keypair[0]  # Match every synonym to its HGNC approved symbol

    except IOError as error:
        sys.stderr.write("PIPELINE ERROR: {0}\n"
                         "You are missing the HGNC gene symbols file! Unable to continue. "
                         "\nTry creating it with the command argument --update and ensure it is "
                         "correctly referenced in the config under hgncPath!\n"
                         "Exiting...\n".format(error))
        sys.exit(1)


def load_tso_genes(tsoPath):
    try:
        with open(tsoPath) as tso:
            data = tso.readlines()
            for line in data:
                # print data[i].rstrip().split("\t")
                gene_cov = line.rstrip().split("\t")
                gene = Gene(panel=None, json=None, on_TSO=True)
                gene._name = gene_cov[0]
                gene.coverage = gene_cov[1]
                tso_genes.append(gene)
    except IOError as error:
        sys.stderr.write("PIPELINE ERROR: {0}\n"
                         "You are missing the coverage file in the location {1}\n"
                         "The file should contain the gene name and coverage, tab seperated (GENE_NAME\\tCOVERAGE\\n)"
                         " The coverage can be 0 (GENE\t0)\nExiting..."
                         .format(error.message, tsoPath))
        sys.exit(1)  # Unable to continue
    except IndexError as error:
        sys.stderr.write("PIPELINE ERROR: {0}\n Is your file format correct?"
                         "The file should contain the gene name and coverage, tab seperated (GENE_NAME\\tCOVERAGE\\n)"
                         .format(error.message))


def get_tso_status(gene):
    # type: (gene) -> Gene
    if len(tso_genes) == 0:
        # Load TSO genes if currently not loaded.
        load_tso_genes(gene.panel.config.tsoGenes)
    # get the coverage from TSO genes (name based matching, if match found --> gene is covered)
    match_coverage = next((g.coverage for g in tso_genes if g.name == gene.name), None)
    if match_coverage is not None:
        gene.coverage = float(match_coverage)
        return True
    return False


def find_gene(name):
    genes = [g for g in tso_genes if g.name == name]
    if len(genes) > 0:
        if len(genes) == 1:
            return genes[0]
        else:
            sys.stderr.write("PIPELINE ERROR: The gene {0} had several matches among the covered genes!\n"
                             .format(name))
            # TODO: What happens if there are several matches?
    else:
        sys.stderr.write("PIPELINE ERROR: The gene {0} was not found among the covered genes.\n".format(name))


def synonyms_to_hgnc(symbol):
    # for synonyms in hgnc_genes.values():
    #   if symbol in synonyms:
    #        return synonyms
    return synonym_hgnc.get(symbol, None)
    # return None


def find_synonyms(symbol):
    """
    Get the synonyms for a HGNC symbol, if it isn't a HGNC Approved Name, try to find it among
    HGNC synonyms, in this case
    :param symbol: The gene symbol to be searched against.
    :return: Return a list of every possible synonym except self
    """
    result = ()
    result = hgnc_genes.get(symbol, None)
    if result is None:
        # Perhaps a synonym was used
        hgnc = synonyms_to_hgnc(symbol)
        if hgnc is not None:
            result = hgnc_genes[hgnc][:]
            result.append(hgnc)
    return result


def find_hgnc(gene):
    """
    This function returns the correct HGNC declared gene symbol regardless if
    the gene was created with a synonymous symbol. If the symbol is not a HGNC Approved Name, HGNC previous symbol or
    HGNC synonym --> return None
    :param symbol: Symbol the Gene object was created with
    :return:
    """
    symbol = gene._name
    if len(hgnc_genes) > 0:
        if hgnc_genes.has_key(symbol):
            return symbol
        else:
            return synonyms_to_hgnc(symbol)
    else:
        load_hgnc_genes(gene.panel.config.hgncPath)


class Gene(object):
    def __init__(self, panel, json=None, on_TSO=None):
        # type: (GenePanel, dict, bool) -> Gene

        if json is not None:
            self.ensemblegeneids = json['EnsembleGeneIds']
            self._name = str(json['GeneSymbol'])
            self.panel = panel
            self.panel_name = str(panel.name)
            self.__confidence = json['LevelOfConfidence']
            self.modeofinheritance = json['ModeOfInheritance']
            self.modeofpathogenicity = json['ModeOfPathogenicity']
            self.penetrance = json['Penetrance']
            self.phenotypes = json['Phenotypes']
            self.raw_json = json
            self.coverage = -1.0
        else:
            self.ensemblegeneids = None
            self._name = None
            self.panel = None
            self.__confidence = 1
            self.modeofinheritance = None
            self.modeofpathogenicity = None
            self.penetrance = None
            self.phenotypes = None
            self.raw_json = None
            self.coverage = -1.0
        if on_TSO is None:
            self.on_TSO = get_tso_status(self)
        else:
            self.on_TSO = on_TSO
            #self.synonyms = find_synonyms(self._name)

    @property
    def evidence_level(self):
        return self.__confidence

    @property
    def name(self):
        # Every time the gene name is used, the name property will either return the HGNC Approved name
        # or throw an AssertionError, see also the find_synonyms() docstring.
        result = find_hgnc(self)
        assert result is not None, "Tried to use a gene object with a nonexistant symbol: {0}" \
            .format(self._name)
        return result

    def __str__(self):
        return self.name
