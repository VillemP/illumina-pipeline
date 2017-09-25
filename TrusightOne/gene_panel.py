from distutils.version import LooseVersion

import gene
from gene import Gene


class GenePanel:
    def __init__(self, panel_json):
        self.name = replace_illegal_chars(panel_json['Name'])
        self.panel_id = panel_json['Panel_Id']
        self.diseasegroup = panel_json['DiseaseGroup']
        self.diseasesubgroup = panel_json['DiseaseSubGroup']
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
        return "Name={0}, Id={1}, JSON={2}, Total genes={3}" \
            .format(self.name, self.panel_id, self.json, len(self.genes))

    @property
    def as_table(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}" \
            .format(self.name, self.panel_id, self.version, self.diseasegroup, self.diseasesubgroup,
                    len(self.genes),
                    len([lambda: g for g in self.genes if g.confidence == gene.GREEN]),
                    len([lambda: g for g in self.genes if g.confidence == gene.AMBER]),
                    [g.name.encode("ascii") for g in self.genes if g.confidence == gene.GREEN])


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
    if panel1.version >= panel2.version:
        return True
    return False
