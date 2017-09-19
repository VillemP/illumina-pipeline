from distutils.version import LooseVersion


class GenePanel:
    def __init__(self, panel_json):
        self.name = replace_illegal_chars(panel_json['Name'])
        self.panel_id = panel_json['Panel_Id']
        self.json = panel_json
        self.genes = list()
        self.version = LooseVersion(panel_json['CurrentVersion'])
        self.genes_json = None

    def add_genes(self, unpacked_json):
        self.genes_json = unpacked_json
        for gene in self.genes_json['result']['Genes']:
            self.genes.append(gene)

    def __str__(self):
        return "Name={0}, Id={1}, JSON={2}, Total genes={3}".format(self.name, self.panel_id, self.json, len(self.genes))


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
