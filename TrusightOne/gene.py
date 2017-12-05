import sys

GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"

tso_genes = list()
hgnc_genes = {}
synonym_hgnc = {}


def load_hgnc_genes(hgncPath):
    """
    This function loads the HGNC genes from local data into a dictionary, with unique symbols corresponding to
    1) lists of symbol synonyms (hgnc_genes dict)
    2) unique synonyms to their HGNC symbol (synonym_hgnc dict)
    :param hgncPath: Direct path to the HGNC containing text file (which can be created with the -update tool)
    """
    try:
        with open(hgncPath) as hgnc:
            data = hgnc.readlines()
            for line in data:
                clean_list = []
                symbols = line.split("\t")
                # Match every HGNC symbol to a list of synonyms
                synonym_list = symbols[1].strip().split(',') + symbols[2].strip().split(',')
                for synonym in synonym_list:
                    if synonym != "" and synonym is not None:
                        clean_list.append(synonym.strip())
                hgnc_genes[symbols[0]] = clean_list
                # Match every synonym to it's HGNC symbol to enable quicker searching
        for key, value in hgnc_genes.items():
            for symbol in value:  # values
                synonym_hgnc[symbol] = key  # Match every synonym to its HGNC approved symbol

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
                gene_symbol, coverage = line.rstrip().split("\t")
                gene = Gene(panel=None, json=None, on_TSO=True)
                gene._name = gene_symbol
                gene.coverage = coverage
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
            result = hgnc
    return result


def find_hgnc(gene):
    # type: (Gene) -> str
    """
    This function returns the correct HGNC declared gene symbol regardless if
    the gene was created with a synonymous symbol. If the symbol is not a HGNC Approved Name, HGNC previous symbol or
    HGNC synonym --> return None
    :param symbol: string symbol the Gene object was created with
    :return:
    """
    # If the correct HGNC term was found in an earlier call, start using it from here on
    if gene._hgnc is not None:
        return gene._hgnc
    symbol = gene._name
    if len(hgnc_genes) > 0:
        if symbol in hgnc_genes:
            gene._hgnc = symbol
            return symbol
        else:
            result = synonyms_to_hgnc(symbol)
            gene._hgnc = result
            return result
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
            self._hgnc = None
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
            self._hgnc = None
        # This is used for creating Genes from JSONs
        if on_TSO is None:
            self.on_TSO = get_tso_status(self)
        # Genes from local text file
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
        if result is None:
            # sys.stderr.write("Tried to use a gene object with a HGNC-nonexistant symbol: {0}\n".format(self._name))
            return self._name
        return result

    def __str__(self):
        return self.name
