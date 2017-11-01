import os

import yaml

from pipeline_utility.baseconfig import BaseConfig


class TruesightOneConfig(BaseConfig):
    yaml_tag = u"!TruesightOneConfig"
    BaseConfig.hidden_fields.append('db_vcf_dir')

    def __init__(self, filepath, json_dir='', name='', vcf_storage_location='', db_vcf_list_name='', db_directory='',
                 db_name='',
                 logfile="Miniseq-log.txt", toolkit='directory to GenomeAnalysisTK.jar',
                 reference='directory to ucsc.hg19.fasta',
                 targetfile=".bed target file", refseq=".refSeq file", annotator=None, annotation_db=None,
                 custom_annotation_dir=None, snpsift="directory to Snpsift.jar", tsoGenes="TSO_coverage.txt",
                 gene_table=None):
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
        self.tsoGenes = tsoGenes  # Location of the gene\tcoverage list
        self.gene_table = gene_table


def loadCfg(cfg_path):
    yaml.add_constructor(TruesightOneConfig.yaml_tag, TruesightOneConfig.cfg_constructor)

    cfg = TruesightOneConfig(cfg_path)
    cfg = cfg.load()

    return cfg
