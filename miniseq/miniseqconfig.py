import os

from pipeline_utility.baseconfig import BaseConfig


class MiniseqConfig(BaseConfig):
    def __init__(self, filepath):
        super(MiniseqConfig, self).__init__(filepath)
        self.vcf_storage_location = "/media/kasutaja/data/TSC_temp/miniseq_pipe/vcfs/"
        self.db_vcf_list_name = "vcfs-sample-path.list"
        self.db_location = "/media/kasutaja/data/NGS_data/var_db_miniseq/"
        self.db_vcf_dir = os.path.join(self.db_location, self.db_vcf_list_name)
        self.db_name = "miniseq-named-targeted-merged-n"
        self.db_dir = os.path.join(self.db_location, self.db_name)
        self.workingDir = os.getcwd()
        self.project = os.path.basename(os.path.normpath(self.workingDir))
        self.logfile = "Miniseq-log-{0}.txt".format(self.project)
