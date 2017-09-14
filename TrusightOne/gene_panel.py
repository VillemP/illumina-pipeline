class genepanel:
    def __init__(self, name, panel_id, json, version):
        self.name = name
        self.panel_id = panel_id
        self.json = json
        self.genes = list()
        self.version = version
        self.genes_json = None

    def add_genes(self, json):
        self.genes_json = json
        for gene in json:
            self.genes.append(gene)

    def __str__(self):
        return "Name={0}, Id={1}, JSON={2}".format(self.name, self.panel_id, self.json)