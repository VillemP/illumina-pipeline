import sys
GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"
WARNINGS = 2  # Report all
EXCEPTIONS = 1  # Report potentially impactful events

class HgncHandler:

    def __init__(self, hgncPath, tsoPath, verbosity=WARNINGS):
        self.verbosity = verbosity
        self.tso_genes = list()
        self.hgnc_genes = {}
        self.synonym_hgnc = {}
        self.load_hgnc_genes(hgncPath)
        if tsoPath is not None:  # the HgncConverterTool only requires HGNC genes, thus tsoGenes are not loaded
            self.load_tso_genes(tsoPath)

    def load_hgnc_genes(self, hgncPath):
        """
        This function loads the HGNC genes from local data into a dictionary, with unique symbols corresponding to
        1) lists of symbol synonyms (hgnc_genes dict)
        2) unique synonyms to their HGNC symbol (synonym_hgnc dict)
        TODO: Make it case-insensitive
        :param hgncPath: Direct path to the HGNC containing text file (which can be created with the -update tool)
        """
        try:
            sys.stderr.write("Loading HGNC genes from {0}:\n".format(hgncPath))
            with open(hgncPath) as hgnc:
                data = hgnc.readlines()
                all_synonyms = set()
                for line in data:
                    clean_set = set()
                    symbols = line.split("\t")
                    # Match every HGNC symbol to a list of synonyms

                    synonym_list = symbols[1].strip().split(',') + symbols[2].strip().split(',')
                    for synonym in synonym_list:
                        if synonym != "" and synonym is not None:
                            synonym = synonym.strip()
                            # Only use unique synonyms
                            if synonym not in all_synonyms:
                                all_synonyms.add(synonym)
                                clean_set.add(synonym)
                                clean_set.add(
                                    synonym.upper())  # The synonym contains lower-chars, add these to the posible synonym bank
                            else:
                                # sys.stderr.write("Found a non-unique synonym: {0}\n".format(synonym))
                                pass

                    self.hgnc_genes[symbols[0]] = clean_set
                    # Match every synonym to it's HGNC symbol to enable quicker searching
            for key, value in self.hgnc_genes.items():
                for symbol in value:  # values
                    self.synonym_hgnc[symbol] = key  # Match every synonym to its HGNC approved symbol
            sys.stderr.write(
                "Loaded {0} HGNC genes and {1} HGNC synonyms.\n".format(len(self.hgnc_genes), len(self.synonym_hgnc)))

        except IOError as error:
            sys.stderr.write("PIPELINE ERROR: {0}\n"
                             "You are missing the HGNC gene symbols file! Unable to continue. "
                             "\nTry creating it with the command argument --update and ensure it is "
                             "correctly referenced in the config under hgncPath!\n"
                             "Exiting...\n".format(error))
            sys.exit(1)

    def load_tso_genes(self, tsoPath):
        try:
            with open(tsoPath) as tso:
                data = tso.readlines()
                for line in data:
                    # print data[i].rstrip().split("\t")
                    gene_symbol, coverage = line.rstrip().split("\t")
                    gene = Gene(panel=None, hgncHandler=self, json=None, on_TSO=True)
                    gene._name = gene_symbol
                    gene.coverage = coverage
                    self._set_hgnc(gene)
                    self.tso_genes.append(gene)
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

    def get_tso_status(self, gene):
        # type: (gene) -> bool
        if len(self.tso_genes) == 0:
            # Load TSO genes if currently not loaded.
            self.load_tso_genes(gene.panel.config.tsoGenes)
        # get the coverage from TSO genes (name based matching, if match found --> gene is covered)
        name = gene.name
        match_coverage = next((g.coverage for g in self.tso_genes if g._hgnc == name.upper()), None)
        if match_coverage is not None:
            gene.coverage = float(match_coverage)
            return True
        return False

    def find_gene(self, name, verbose=False):
        """
        Find the Gene object from the genes created from config.tsoGenes.
        Compares the gene name (which is converted to HGNC if possible) to the input,
        which is also converted to HGNC if possible. Thus returns a gene regardless of the original
        name and the synonym it was searched with.
        :param name: Gene symbol
        :return: gene.Gene | None
        """
        # hgnc = self.hgnc_genes.get(name, None)
        try:
            hgnc = self.hgnc_genes[name]
        except KeyError:
            hgnc = None
        # HGNC is None -> might be a synonym
        if hgnc is None:
            hgnc = self.synonyms_to_hgnc(name)
            if hgnc is not None:
                name = hgnc  # We found a match!
        else:
            name = name  # for clarity, it was already a HGNC symbol

        if hgnc is not None:
            genes = [g for g in self.tso_genes if g._hgnc == name]
            if len(genes) > 0:
                if len(genes) == 1:
                    return genes[0]
                else:
                    if self.verbosity > EXCEPTIONS or verbose:
                        sys.stderr.write("PIPELINE ERROR: The gene {0} had several matches among the covered genes!\n"
                                         .format(name))
                    # TODO: What happens if there are several matches?
            else:
                if self.verbosity >= WARNINGS or verbose:
                    sys.stderr.write(
                        "PIPELINE ERROR: The gene {0} was not found among the covered genes.\n".format(name))
                return None
        else:
            if self.verbosity >= WARNINGS or verbose:
                sys.stderr.write("The gene {0} was not found among HGNC names and synonyms.\n".format(name))

    def synonyms_to_hgnc(self, symbol):
        """
        Convert a synonym to HGNC. HGNC symbols have to be preloaded (load_hgnc_genes)
        :param symbol: gene name
        :return: HGNC symbol or None
        """
        assert len(self.synonym_hgnc) > 0, "HGNC symbols have not been loaded!"
        return self.synonym_hgnc.get(symbol, None)
        # return None

    def find_synonyms(self, symbol):
        """
        Get the synonyms for a HGNC symbol, if it isn't a HGNC Approved Name, try to find it among
        HGNC synonyms, in this case return a list of every possible synonym except self
        :param symbol: str | The gene symbol to be searched against.
        :return: list | Return a list of every possible synonym except self
        """
        result = []
        result = self.hgnc_genes.get(symbol, [])

        if result is None:
            # Perhaps a synonym was used
            hgnc = self.synonyms_to_hgnc(symbol)
            if hgnc is not None:
                result = self.hgnc_genes.get(hgnc, [])
                if len(result) > 0:
                    result.append(hgnc)  # return a full list of synonyms, including HGNC symbol
        return result

    def is_hgnc(self, name):
        if self.hgnc_genes.get(name, None) is not None:
            return True
        return False

    def _set_hgnc(self, gene):
        # type: (Gene) -> str
        """
        This function returns the correct HGNC declared gene symbol regardless if
        the gene was created with a synonymous symbol. If the symbol is not a HGNC Approved Name, HGNC previous symbol or
        HGNC synonym --> return None
        Some synonymous symbols exist match up to several HGNC names
        :param symbol: string symbol the Gene object was created with
        :return:
        """
        # If the correct HGNC term was found in an earlier call, start using it from here on
        if gene._hgnc is not None:
            return gene._hgnc
        symbol = gene._name
        if len(self.hgnc_genes) > 0:
            if symbol in self.hgnc_genes:
                gene._hgnc = symbol
                return symbol
            else:
                result = self.synonyms_to_hgnc(symbol)
                gene._hgnc = result
                return result
        else:
            self.load_hgnc_genes(gene.panel.config.hgncPath)
            # _set_hgnc(gene)


class Gene(object):
    def __init__(self, panel, hgncHandler=None, json=None, on_TSO=None):
        # type: (GenePanel, HgncHandler, dict, bool) -> self
        if hgncHandler is not None:
            self.hgncHandler = hgncHandler
        else:
            self.hgncHandler = panel.hgncHandler
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
        self.hgncHandler._set_hgnc(self)
        if on_TSO is None:
            self.on_TSO = self.hgncHandler.get_tso_status(self)
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
        # or set its own name
        result = self.hgncHandler._set_hgnc(self)
        if result is None:
            # sys.stderr.write("Tried to use a gene object with a HGNC-nonexistant symbol: {0}\n".format(self._name))
            return self._name
        return result

    def __str__(self):
        if self._hgnc == self._name:
            return self._hgnc
        return "{0} (HGNC: {1})".format(self._name, self._hgnc)
