import os

from pipeline_utility.baseconfig import BaseConfig


class MiniseqConfig(BaseConfig):
    yaml_tag = u"!MiniseqConfig"
    BaseConfig.hidden_fields.append('db_vcf_dir')

    def __init__(self, filepath, vcf_storage_location='', db_vcf_list_name='', db_directory='', db_name='',
                 logfile="Miniseq-log.txt"):
        super(MiniseqConfig, self).__init__(filepath)
        self.vcf_storage_location = vcf_storage_location
        self.db_vcf_list_name = db_vcf_list_name
        self.db_directory = db_directory
        self.db_vcf_dir = os.path.join(self.db_directory, self.db_vcf_list_name)
        self.db_name = db_name
        self.logfile = "Miniseq-log.txt"
