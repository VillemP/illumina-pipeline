import itertools
import os
import sys
from distutils.version import LooseVersion

import numpy

import gene
from gene import Gene
from pipeline_utility.txttoxlsx_filtered import create_excel

genesdict = {}

table_format = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}" \
    .format("Name", "id", "version", "diseasegroup", "diseasesubgroup",
            "Total_genes",
            "green_genes",
            "avg_coverage_on_tso",
            "green_genes",
            "Covered_green_genes",
            "Covered_green_genes",
            "amber_genes",
            "avg_coverage_amber",
            "amber_genes",
            "Covered_amber_genes",
            "Covered_amber_genes")


class GenePanel(object):
    def __init__(self, hgnchandler, panel_json, config=None):
        self.name = replace_illegal_chars(panel_json['Name'])
        self.id = panel_json['Panel_Id']
        self.diseasegroup = panel_json['DiseaseGroup']
        self.diseasesubgroup = panel_json['DiseaseSubGroup']
        # Actually an "unpacked" JSON, so it's a dict
        self.json = panel_json
        self.genes = list()
        self.version = LooseVersion(panel_json['CurrentVersion'])
        self.genes_json = None
        # the config file for the Jsonhandler and pipeline
        self.config = config
        self.hgncHandler = hgnchandler

    def add_genes(self, unpacked_json):
        self.genes_json = unpacked_json
        for gene in self.genes_json['result']['Genes']:
            # Use preexisting gene objects to save time on searching for name matches and so forth
            # if str(gene['GeneSymbol']) in genesdict:
            #    self.genes.append(genesdict[gene['GeneSymbol']])
            # else:
            g = Gene(self, json=gene, hgncHandler=self.hgncHandler)
            self.genes.append(g)
            # genesdict[gene['GeneSymbol']] = g

    def __str__(self):
        return "Name={0}, Id={1}, Version={2}, Total genes={3}, Covered genes={4}" \
            .format(self.name, self.id, self.version, len(self.genes), len(self.tso_genes))

    def contains_gene(self, gene):
        return [g for g in self.genes if g.name == gene]

    @property
    def tso_genes(self):
        return [g for g in self.genes if g.on_TSO]

    @property
    def avg_coverage_on_tso(self):
        genes = [float(g.coverage) for g in self.tso_genes if float(g.coverage) != -1]
        try:
            return sum(genes) / float(len(genes))
        except TypeError as e:
            raise e
            # return None

    @property
    def green_genes(self):
        return [g for g in self.genes if g.evidence_level == gene.GREEN]

    @property
    def amber_genes(self):
        return [g for g in self.genes if g.evidence_level == gene.AMBER]

    @property
    def avg_coverage_GREEN(self):
        genes_cov = [g.coverage for g in self.tso_genes if g.evidence_level == gene.GREEN]
        if len(genes_cov) > 0:
            return sum(genes_cov) / float(len(genes_cov))
        return 0.0

    @property
    def avg_coverage_AMBER(self):
        genes_cov = [g.coverage for g in self.tso_genes if g.evidence_level == gene.AMBER]
        if len(genes_cov) > 0:
            return sum(genes_cov) / float(len(genes_cov))
        return 0.0

    @property
    def as_table(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}" \
            .format(self.name, self.id, self.version, self.diseasegroup, self.diseasesubgroup,
                    len(self.genes),
                    len(self.green_genes),
                    round(self.avg_coverage_on_tso, 2),
                    [g.name.encode("ascii") for g in self.green_genes],
                    len([g for g in self.tso_genes if g.evidence_level == gene.GREEN]),
                    [g.name.encode("ascii") for g in self.tso_genes if g.evidence_level == gene.GREEN],
                    len(self.amber_genes),
                    round(self.avg_coverage_AMBER, 2),
                    [g.name.encode("ascii") for g in self.amber_genes],
                    len([g for g in self.tso_genes if g.evidence_level == gene.AMBER]),
                    [g.name.encode("ascii") for g in self.tso_genes if g.evidence_level == gene.AMBER])


def replace_illegal_chars(text):
    chars = "\\`*{}[]()<>#+.!$/"
    for c in chars:
        text = text.replace(c, "_")
    return text


def compare_versions(panel1, panel2):
    """
    Returns True if panel1 version is greater or equal to the version of panel2.
    :rtype: Boolean
    """
    if panel1.version > panel2.version:
        return True
    return False


# Compile panels into a dict with the key matching the panel annotation and values being panels.
class CombinedPanels(dict):
    def __init__(self, handler):
        super(CombinedPanels, self).__init__([])
        self.handler = handler
        self.hgncHandler = handler.hgncHandler
        if len(handler.panels) > 0:
            all_genes_panel = GenePanel(self.hgncHandler,
                                        {'Name': 'All covered genes', 'Panel_Id': '0000', 'DiseaseGroup': 'None',
                                         'DiseaseSubGroup': 'None', 'CurrentVersion': '1.3'})
            if len(self.hgncHandler.hgnc_genes) == 0:
                self.hgncHandler.load_hgnc_genes(handler.config.hgncPath)
            if len(self.hgncHandler.tso_genes) == 0:
                self.hgncHandler.load_tso_genes(handler.config.tsoGenes)
            all_genes_panel.genes = self.hgncHandler.tso_genes
            handler.panels.append(all_genes_panel)
            self[('ALL', 'KOGU')] = [panel for panel in handler.panels if panel.id == '0000']
            self[('VAM', 'ID', 'Intellectual disability')] = [panel for panel in handler.panels
                                                              if panel.name == "Intellectual disability"
                                                              and panel.id == '558aa423bb5a16630e15b63c']
            self[('EP', 'Epilepsy')] = [panel for panel in handler.panels
                                        if panel.diseasesubgroup == 'Inherited Epilepsy Syndromes']
            self[('MITO', 'Mitochondrial')] = [panel for panel in handler.panels
                                               if panel.id == '55928cf522c1fc4f7d26e960']
            self[('METABO', 'AV', 'Undiagnosed metabolic disorders')] = [panel for panel in handler.panels
                                                                         if panel.id == '5763f1518f620350a22bccdb']
            self[('DEAF', 'Kuulmislangus')] = [panel for panel in handler.panels
                                               if panel.diseasegroup == 'Hearing and ear disorders']
            self[('EYE')] = [panel for panel in handler.panels
                             if panel.diseasegroup == 'Ophthalmological disorders']
            self[('CMT', "Charcot-Marie-Tooth")] = [panel for panel in handler.panels
                                                    if panel.id == '55ad205422c1fc7041340234']
            self[('HSP', 'Hereditary spastic paraplegy')] = [panel for panel in handler.panels
                                                             if panel.id == '55ad019f22c1fc7042059038']
            self[('MUSC', 'Muscle', 'Lihas')] = [panel for panel in handler.panels
                                                 if panel.diseasesubgroup == 'Neuromuscular disorders']
            self[('ARRHYTHMIA', 'ARR')] = [panel for panel in handler.panels
                                           if panel.diseasesubgroup == 'Cardiac arrhythmia'
                                           or panel.diseasesubgroup == 'Arrhythmogenic Right Ventricular Cardiomyopathy']
            self[('CARDIOMYOPATHY')] = [panel for panel in handler.panels
                                        if panel.diseasesubgroup == 'Cardiomyopathy']
            self[('MARFAN', 'MARFAN-LIKE')] = [panel for panel in handler.panels
                                               if panel.id == '5596735822c1fc4f7d26e96d']
            self[('SKELETAL')] = [panel for panel in handler.panels
                                  if panel.diseasesubgroup == 'Skeletal dysplasias']
            self[('MODY', 'Diabetes')] = [panel for panel in handler.panels
                                          if panel.id == '55a9238422c1fc6711b0c6c3'
                                          or panel.id == '553f9745bb5a1616e5ed45e9'
                                          or panel.id == '55a9041e22c1fc6711b0c6c0']
            self[('DEMENTIA')] = [panel for panel in handler.panels
                                  if panel.id == '55b6173522c1fc05fc7a1855']
            self[('PARKINSON')] = [panel for panel in handler.panels
                                   if panel.id == '58078e6e8f62030e233a8157']
            self[('LEUKODYSTROPHY')] = [panel for panel in handler.panels
                                        if panel.id == '568f920822c1fc1c79ca177a']
            self[('ATAXIA')] = [panel for panel in handler.panels
                                if panel.id == '559a7d1022c1fc58ad67fc97']
            self[('DYSTONIA')] = [panel for panel in handler.panels
                                  if panel.id == '553f95c9bb5a1616e5ed45bf']
            self[('ANAEMIA')] = [panel for panel in handler.panels
                                 if panel.id == '58a70e858f62037e8779b2e8']
            self[('BLEEDING', 'CLOTTING')] = [panel for panel in handler.panels
                                              if panel.id == '5763f32a8f620350a22bccde']
            self[('SEX')] = [panel for panel in handler.panels
                             if panel.id == '569380ac22c1fc251660faf8']
            self[('EHLERS-DANLOS')] = [panel for panel in handler.panels
                                       if panel.id == '588728f38f62030cf7152165']

    def __getitem__(self, item):
        return super(CombinedPanels, self).__getitem__(
            [key for key in self.iterkeys() if [True for element in key if item == element or item == key]][0])

    def __delitem__(self, item):
        return super(CombinedPanels, self).__delitem__(
            [key for key in self.iterkeys() if [True for element in key if item == element or item == key]][0])

    def get(self, item, default=None):
        return super(CombinedPanels, self).get(
            [key for key in self.iterkeys() if [True for element in key if item == element or item == key]][0], default)

    def __contains__(self, item):
        return super(CombinedPanels, self).__contains__(
            [key for key in self.iterkeys() if [True for element in key if item == element or item == key]][0])

    def table(self):
        lines = list()
        for key, panels in self.iteritems():
            genes = list()

            for panel in panels:
                genes.extend([g.name for g in panel.tso_genes if g.name not in genes])
            # In order for the izip to work properly, there can't be any spaces

            if type(key) is tuple:
                name = ",".join(str(symbol).replace(" ", "_") for symbol in key)
            else:
                name = key
            # First three lines are name, length of genes, mean panel-based coverage
            line = [name, len(genes), numpy.mean([panel.avg_coverage_on_tso for panel in panels])]
            line.extend(genes)
            lines.append(line)
        #
        lines.sort(key=lambda col: col[0])
        table = itertools.izip_longest(*lines, fillvalue="None")
        e = []
        for l in table:
            e.append(l)
        return e

    def write_table(self, out="combined_panels_summary.txt"):
        if not self.handler.loaded:
            self.handler.get_all_panels(False)
        print("Writing the panel combinations to {0}".format(os.path.join(self.handler.config.json_dir, out)))
        with open(os.path.join(self.handler.config.json_dir, out), "w+") as f:
            f.write(table_format + '\n')
            for panelcombo in self.iteritems():
                for panel in panelcombo[1]:
                    f.write(panel.as_table + '\t{}'.format(panelcombo[0]) + '\n')
        print("Writing the combined panels gene table to {0}".format(self.handler.config.gene_table))
        with open(self.handler.config.gene_table, "wb+") as f:
            tbl = self.table()
            # print type(tbl)
            for i, row in enumerate(tbl):
                if i == 0:
                    f.write("\t".join([str(tupl) for tupl in row]) + "\n")
                else:
                    f.write("\t".join([str(col) for col in row]).replace("None", ".") + "\n")
        # Currently you have to use excel's find and replace function to get rid of the empty cel fillvalue
        # TODO: Create a postprocess for every column (i--> can be a 1-indexed column, just need the total columns)
        create_excel(self.handler.config.gene_table + ".xlsx", files=[self.handler.config.gene_table])

        # Writes a tab-delimited list of GENE - PANELNAME associations for annotations.
        with open(os.path.join(self.handler.config.json_dir, "TSO_genepanels.txt"), "w+") as f:
            for panelcombo in self.iteritems():
                for panel in panelcombo[1]:
                    for gene in panel.tso_genes:
                        # Choose the whole string if it's a string, only the first element if it is a tuple
                        f.write("\t".join(
                            [gene.name, "".join(panelcombo[0][0] if type(panelcombo[0]) is tuple else [name for name in
                                                                                                       panelcombo[
                                                                                                           0]])]) + "\n")


def match_order_to_panels(key, combinedpanels, handler):
    try:
        # TODO: currently combinedpanels['-'] will return Ehlers-Danlos, fix this!
        if key != "" and key != '-':
            panel = combinedpanels[key]
            return panel
        else:
            IndexError("Empty key.")
    except IndexError:
        # This key didn't yield a combined panel result, perhaps it is not a custom compiled panel but a single panel
        pass
    if len(handler.panels) > 0:
        match = [panel for panel in handler.panels if panel.name == key]
        if len(match) == 0:
            # Maybe it is an ID based query
            match = [panel for panel in handler.panels if
                     panel.id.upper() == key]  # The keys are stored in a lower case, but our search tool converts all input to upper
            if len(match) > 0:
                return match
            else:
                sys.stderr.write("PIPELINE ERROR: Unknown panel key {0}\n".format(key))
                return []
        elif len(match) > 1:
            # Found at least two matches with this name, this shouldn't happen
            raise ValueError("The panel with the key {0} in the order yielded multiple results! \nPanels {1}"
                             .format(key, match))
        else:
            return match[0]
