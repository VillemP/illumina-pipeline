import os
from distutils.version import LooseVersion as Version

from TrusightOne.jsonhandler import JsonHandler

json_dir = "/media/kasutaja/data/NGS_data/panels_json_main/"


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
    loadPanels()


if __name__ == "__main__":
    main()
