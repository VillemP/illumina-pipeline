import os


class Sample:
    def __init__(self, name, vcflocation, bamlocation):
        self.name = name.rstrip()
        self.vcflocation = vcflocation.rstrip()
        self.bamlocation = bamlocation.rstrip()
        assert os.path.exists(self.vcflocation)
        assert os.path.exists(self.bamlocation)
        assert self.name in self.bamlocation and self.name in self.vcflocation