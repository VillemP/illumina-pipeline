import os
from distutils.version import LooseVersion as Version

import yaml

from TrusightOne.jsonhandler import JsonHandler
from pipeline_utility.baseconfig import BaseConfig


class TruesightOneConfig(BaseConfig):
    yaml_tag = u"!TruesightOneConfig"

    def __init__(self, filepath):
        super(TruesightOneConfig, self).__init__(filepath)
        self.json_dir = ""
        self.name = "asdasd"


json_dir = "/media/kasutaja/data/NGS_data/panels_json_main/"


def createCfg():
    workingdir = os.getcwd()
    # cfg = TruesightOneConfig(os.path.join(workingdir, "tso.yaml"))
    with open(os.path.join(workingdir, "tso.yaml"), "r") as f:
        cfg = yaml.safe_load(f)
    print("JsonDir = " + cfg.json_dir)
    print cfg.name


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
    createCfg()
    #loadPanels()


if __name__ == "__main__":
    main()