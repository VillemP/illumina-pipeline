class GenePanel:
    def __init__(self, name, panel_id, json, version):
        self.name = replace_illegal_chars(name)
        self.panel_id = panel_id
        self.json = json
        self.genes = list()
        self.version = version
        self.genes_json = None

    def add_genes(self, unpacked_json):
        self.genes_json = unpacked_json
        for gene in self.genes_json:
            self.genes.append(gene)

    def __str__(self):
        return "Name={0}, Id={1}, JSON={2}".format(self.name, self.panel_id, self.json)


def replace_illegal_chars(text):
    chars = "\\`*{}[]()<>#+.!$/"
    for c in chars:
        text = text.replace(c, "_")
    return text