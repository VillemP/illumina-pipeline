import os


class Sample:
    def __init__(self, name, vcflocation, bamlocation):
        self.name = name.rstrip()
        self.vcflocation = vcflocation.rstrip()
        self.bamlocation = bamlocation.rstrip()
        self.finished = False
        self.annotated = False
        self.reduced_variants_vcf = None
        assert os.path.exists(self.vcflocation)
        assert os.path.exists(self.bamlocation)
        assert self.name in self.bamlocation and self.name in self.vcflocation

    def __str__(self):
        return "{0} Annotated:{1} Finished:{2}".format(self.name, str(self.annotated), str(self.finished))
