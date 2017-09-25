GREEN = "HighEvidence"
AMBER = "ModerateEvidence"
RED = "LowEvidence"


class Gene:
    def __init__(self, panel, json):
        """

        :type json: dict
        """
        self.ensemblegeneids = json['EnsembleGeneIds']
        self.name = json['GeneSymbol']
        self.panel_name = panel
        self.confidence = json['LevelOfConfidence']
        self.modeofinheritance = json['ModeOfInheritance']
        self.modeofpathogenicity = json['ModeOfPathogenicity']
        self.penetrance = json['Penetrance']
        self.phenotypes = json['Phenotypes']
        self.raw_json = json

    def get_evidence_level(self):
        return self.confidence

    def __str__(self):
        pass
