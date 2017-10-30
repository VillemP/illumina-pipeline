from distutils.version import LooseVersion

import gene
from gene import Gene


class GenePanel(object):
    def __init__(self, panel_json, config=None):
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

    def add_genes(self, unpacked_json):
        self.genes_json = unpacked_json
        for gene in self.genes_json['result']['Genes']:
            g = Gene(self, gene)
            self.genes.append(g)

    def __str__(self):
        return "Name={0}, Id={1}, Version={2}, Total genes={3}" \
            .format(self.name, self.id, self.version, len(self.genes))

    @property
    def tso_genes(self):
        return [g for g in self.genes if g.on_TSO]

    @property
    def avg_coverage_on_tso(self):
        genes = [g.coverage for g in self.tso_genes]
        return sum(genes) / float(len(genes))

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
                    round(self.avg_coverage_GREEN, 2),
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


class CombinedPanels(dict):
    def __init__(self, handler):
        super(CombinedPanels, self).__init__([])
        if len(handler.panels) > 0:
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
