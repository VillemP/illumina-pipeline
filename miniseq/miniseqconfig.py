import os

from pipeline_utility.baseconfig import BaseConfig


class MiniseqConfig(BaseConfig):
    yaml_tag = u"!MiniseqConfig"
    BaseConfig.hidden_fields.append('db_vcf_dir')

    def __init__(self, filepath, vcf_storage_location='', db_vcf_list_name='', db_directory='', db_name='',
                 logfile="Miniseq-log.txt", toolkit='directory to GATK', reference='directory to ucsc.hg19.fasta',
                 targetfile=".bed target file", padding=10, annotator=None, annotation_db=None,
                 custom_annotation_dir=None):
        super(MiniseqConfig, self).__init__(filepath)
        self.vcf_storage_location = vcf_storage_location
        self.db_vcf_list_name = db_vcf_list_name
        self.db_directory = db_directory
        self.db_vcf_dir = os.path.join(self.db_directory, self.db_vcf_list_name)
        self.db_name = db_name
        self.logfile = logfile
        self.toolkit = toolkit
        self.reference = reference
        self.targetfile = targetfile
        self.padding = padding
        self.annotator = annotator
        self.annotation_db = annotation_db
        self.custom_annotation_dir = custom_annotation_dir
