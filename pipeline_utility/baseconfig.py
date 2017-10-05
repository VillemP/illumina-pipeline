import os
from copy import deepcopy

import yaml


class BaseConfig(yaml.YAMLObject):
    __hidden_fields = ["_filepath"]
    yaml_tag = u"!BaseConfig"
    yaml_loader = yaml.SafeLoader
    yaml_dumper = yaml.SafeDumper
    yaml_flow_style = False

    def __init__(self, filepath):
        super(BaseConfig, self).__init__()
        # self.yaml_flow_style =
        self._filepath = filepath

    @property
    def filepath(self):
        return self._filepath

    @classmethod
    def to_yaml(cls, dumper, data):
        new_data = deepcopy(data)
        for item in cls.__hidden_fields:
            del new_data.__dict__[item]
        return dumper.represent_yaml_object(cls.yaml_tag, new_data, cls,
                                            flow_style=cls.yaml_flow_style)

    def save(self):
        with open(self._filepath, "wb+") as f:
            yaml.safe_dump(self, f)
            # cfg = yaml.safe_dump(self, f)
        return self

    def load(self):
        # Load variables from yaml into the config class
        if not os.path.exists(self._filepath):
            self.save()
            print("Created a new empty config file!")
        with open(self._filepath, "r") as f:
            yaml.safe_load(f)
        return self


vcf_storage_location = "/media/kasutaja/data/TSC_temp/miniseq_pipe/vcfs/"
db_vcf_list_name = "vcfs-sample-path.list"
db_location = "/media/kasutaja/data/NGS_data/var_db_miniseq/"
db_vcf_dir = os.path.join(db_location, db_vcf_list_name)
db_name = "miniseq-named-targeted-merged-n"
db_dir = os.path.join(db_location, db_name)
workingDir = os.getcwd()
project = os.path.basename(os.path.normpath(workingDir))
logfile = "Miniseq-log-{0}.txt".format(project)
