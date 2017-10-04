from distutils.version import LooseVersion

import gene
from gene import Gene


class GenePanel(object):
    def __init__(self, panel_json):
        self.name = replace_illegal_chars(panel_json['Name'])
        self.panel_id = panel_json['Panel_Id']
        self.diseasegroup = panel_json['DiseaseGroup']
        self.diseasesubgroup = panel_json['DiseaseSubGroup']
        # Actually an "unpacked" JSON, so it's a dict
        self.json = panel_json
        self.genes = list()
        self.version = LooseVersion(panel_json['CurrentVersion'])
        self.genes_json = None

    def add_genes(self, unpacked_json):
        self.genes_json = unpacked_json
        for gene in self.genes_json['result']['Genes']:
            g = Gene(self.name, gene)
            self.genes.append(g)

    def __str__(self):
        return "Name={0}, Id={1}, Version={2}, Total genes={3}" \
            .format(self.name, self.panel_id, self.version, len(self.genes))

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
            .format(self.name, self.panel_id, self.version, self.diseasegroup, self.diseasesubgroup,
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
