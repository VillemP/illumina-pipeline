import os
from distutils.version import LooseVersion as Version

import yaml

from TrusightOne.jsonhandler import JsonHandler
from pipeline_utility.baseconfig import BaseConfig


class TruesightOneConfig(BaseConfig):
    yaml_tag = u"!TruesightOneConfig"

    def __init__(self, filepath, json_dir='', name='', important_data=''):
        super(TruesightOneConfig, self).__init__(filepath)
        self.json_dir = json_dir
        self.name = name
        self.important_data = important_data


json_dir = "/media/kasutaja/data/NGS_data/panels_json_main/"


def loadCfg():
    workingdir = os.getcwd()
    cfg_path = os.path.join(workingdir, "tso.yaml")

    yaml.add_constructor(TruesightOneConfig.yaml_tag, TruesightOneConfig.cfg_constructor)

    cfg = TruesightOneConfig(cfg_path)
    cfg = cfg.load()

    return cfg


def loadPanels():
    handler = JsonHandler(json_dir)
    panels = handler.get_all_panels(external=True)

    # Selecting and sorting for tables above version 1.0, ordered by diseasegroup and subgroup
    currentlist = [p for p in panels if p.version >= Version("1.0")]
    currentlist.sort(key=lambda panel: (panel.diseasegroup, panel.diseasesubgroup))

    with open(os.path.join(json_dir, "query.txt"), "w+") as f:
        for panel in currentlist:
            f.write(panel.as_table + '\n')


def main():
    loadCfg()
    #loadPanels()


if __name__ == "__main__":
    main()