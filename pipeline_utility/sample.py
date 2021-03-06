import os


class Sample(object):
    def __init__(self, name, vcflocation, bamlocation):
        super(Sample, self).__init__()
        if vcflocation is not None and bamlocation is not None:
            self.error = False
            self.exceptions = []
            self.dns = None  # Samples in the database
            self.ts = None  # Total samples
            self.name = name.rstrip()
            self.vcflocation = vcflocation.rstrip()
            self.bamlocation = bamlocation.rstrip()
            self.finished = False
            self.annotated = False
            self.reduced_variants_vcf = None  # Currently not in use but can be used to reduce variants before annotation
            self.genes_tempfile = None  # The tempfile containing the self.final_order output (for vcf_manipulator)
            self.panels = list()  # The panels ordered to be annotated (panel, annotation_string)
            self.genes = list()  # The genes ordered to be annotated (Gene, annotation_string) tuple
            self.trash = list()  # Files to be deleted that are byproducts of annotators/toolkits
            self.temptargetfile = None  # The target.bed file specific for this sample's panel/gene order
            self.temprefseq = None
            self.table_files = list()  # The output files that are converted into excel files
            assert os.path.exists(self.vcflocation)  # The input VCF
            assert os.path.exists(self.bamlocation)  # The input BAM
            assert self.name in self.bamlocation and self.name in self.vcflocation  # Assure all the sample specific
            # input names are the same (there shouldn't be mixups)
        else:
            self.error = True
            raise IOError("VCF and BAM files for sample {0} not found!".format(name))

    def __str__(self):
        if not self.finished:
            result = "{0} Annotated:{1} Finished:{2} Files {3} TEMP: {4}".format(self.name, str(self.annotated),
                                                                                 str(self.finished),
                                                                                 str(self.table_files),
                                                                                 str([self.temptargetfile,
                                                                                      self.temprefseq,
                                                                                      self.genes_tempfile]))
        else:
            result = "{0} Annotated:{1} Finished:{2} Files {3}".format(self.name, str(self.annotated),
                                                                       str(self.finished),
                                                                       str(self.table_files))
        return result

    @property
    def final_order(self):
        """
        Custom genes annotation ordered for this sample.
        :return: List of tab separated GENE   ANNOTATION lines to be used with vcf_manipulator
        that adds a new column to the output excel table where the variant matches the gene.
        """
        order = []
        for panel_order in self.panels:
            for gene in panel_order[0].tso_genes:
                order.append("\t".join((gene.name, panel_order[1])))
        for g_order in self.genes:
            order.append("\t".join((g_order[0].name, g_order[1])))

        return order

    @property
    def order_list(self):
        """
        Unique gene symbols (already converted to HGNC) that are to be annotated
        :return: Gene list
        """
        order = []
        for panel_order in self.panels:
            # panel_order is a (panel, panel annotation) tuple
            for gene in panel_order[0].tso_genes:
                order.append(gene.name)
        for g_order in self.genes:
            # gene order is a (gene, annotation string) tuple
            order.append(g_order[0].name)
        return order
