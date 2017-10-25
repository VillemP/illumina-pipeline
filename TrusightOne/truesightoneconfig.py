from pipeline_utility.baseconfig import BaseConfig


class TruesightOneConfig(BaseConfig):
    yaml_tag = u"!TruesightOneConfig"

    def __init__(self, filepath, json_dir='', name='', ):
        super(TruesightOneConfig, self).__init__(filepath)
        self.json_dir = json_dir
        self.name = name
        self.vcf_storage_location = vcf_storage_location
        self.db_vcf_list_name = db_vcf_list_name
        self.db_directory = db_directory
        self.db_vcf_dir = os.path.join(self.db_directory, self.db_vcf_list_name)
        self.db_name = db_name
        self.logfile = logfile
        self.toolkit = toolkit
        self.reference = reference
        self.targetfile = targetfile
        self.refseq = refseq
        self.annotator = annotator
        self.annotation_db = annotation_db
        self.custom_annotation_dir = custom_annotation_dir
        self.snpsift = snpsift
