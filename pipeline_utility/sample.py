import os


class Sample(object):
    def __init__(self, name, vcflocation, bamlocation):
        super(Sample, self).__init__()
        if vcflocation is not None and bamlocation is not None:
            self.error = False
            self.name = name.rstrip()
            self.vcflocation = vcflocation.rstrip()
            self.bamlocation = bamlocation.rstrip()
            self.finished = False
            self.annotated = False
            self.reduced_variants_vcf = None
            self.genes_tempfile = None
            self.panels = list()
            self.genes = list()
            self.targetfile = None
            self.table_files = list()
            assert os.path.exists(self.vcflocation)
            assert os.path.exists(self.bamlocation)
            assert self.name in self.bamlocation and self.name in self.vcflocation
        else:
            self.error = True
            raise IOError("VCF and BAM files for sample {0} not found!".format(name))

    def __str__(self):
        return "{0} Annotated:{1} Finished:{2}".format(self.name, str(self.annotated), str(self.finished))

    @property
    def final_order(self):
        order = []
        for panel_order in self.panels:
            for gene in panel_order[0].tso_genes:
                order.append("\t".join((gene.name, panel_order[1])))
        for g_order in self.genes:
            order.append("\t".join((g_order[0].name, g_order[1])))

        return order
