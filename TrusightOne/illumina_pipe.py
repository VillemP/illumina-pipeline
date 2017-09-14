import os
import pickle

from TrusightOne.jsonhandler import jsonhandler

json_dir = "/media/kasutaja/data/NSG_data/panels_json/"


def main():
    handler = jsonhandler()
    panels = handler.get_all_panels()
    for panel in panels:
        print panel
        if not os.path.exists(os.path.join(json_dir, panel.name)):
            os.makedirs(os.path.join(json_dir, panel.name))
        with open(os.path.join(os.path.join(json_dir, panel.name),
                               panel.name + ".json"), "wb+") as json_file:
            pickle.dump(panel.json, json_file)
        with open(os.path.join(os.path.join(json_dir, panel.name),
                               panel.name + "genes.json", "wb+")) as genes_file:
            pickle.dump(panel.genes_json, genes_file)

if __name__ == "__main__":
    main()