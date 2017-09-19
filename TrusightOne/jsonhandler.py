import json
import os
import time

import requests

from TrusightOne.gene_panel import GenePanel
from miniseq import file_utility


class JsonHandler:
    def __init__(self, db_dir):
        self.panels = list()
        self.json_dir = db_dir

    def open_page(self, url, args=None):
        r = requests.get(url, args)
        # return r.json()
        return r

    def query(self, url, args=None):
        response = self.open_page(url, args)
        return response, response.json()

    def get_genes_per_panel(self, panel):
        q, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/get_panel/"+panel.panel_id)
        print("Got {0} genes for panel '{1}'".format(len(data['result']['Genes']), panel.name))
        return q, data

    def save_panel(self, panel, location):
        # Create directory for panel (JSON_dir/panel/panel.json
        # TODO: save into MySQL?
        if not os.path.exists(os.path.join(self.json_dir, panel.name)):
            os.makedirs(os.path.join(self.json_dir, panel.name))
        with open(os.path.join(os.path.join(self.json_dir, panel.name),
                               panel.name + ".json"), "wb+") as json_file:
            json.dump(panel.json, json_file, sort_keys=True, indent=4)

        # Create directory and file for panel's gene (JSON_dir/panel/panel.genes
        # Panel.genes is also a json
        # TODO: Join the panel json and genes json or load directly into MySQL
        with open(os.path.join(os.path.join(self.json_dir, panel.name),
                               panel.name + ".genes"), "wb+") as genes_file:
            # pickle.dump(panel.genes_json, genes_file)
            json.dump(panel.genes_json, genes_file, sort_keys=True, indent=4)

    def load_panels(self, location):
        panels = file_utility.find_filetype(location, ".json")
        if len(panels) > 1:
            for panel in panels:
                file_location = panel[1]
                with open(file_location, "r") as f:
                    newpanel = GenePanel(json.load(f))
                    genelist = file_utility.find_filetype(os.path.dirname(file_location), ".genes")
                    for genes in genelist:
                        with open(genes[1], "r") as gene_json:
                            newpanel.add_genes(json.load(gene_json))
                            print "Loaded {0} genes into panel {1}".format(len(newpanel.genes), newpanel.name)
                    self.panels.append(newpanel)
            print "Loaded {0} panels from local data.".format(len(self.panels))
            return self.panels
        else:
            return None

    def update_panel(self):
        pass

    def get_all_panels(self, external=True):
        latest_panels = list()
        json_response, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/list_panels",
                           [{'format','json'}])
        print("Got all panels ({0}).".format(len(data['result'])))
        if external:
            for panel in data['result']:
                g_panel = GenePanel(panel)
                genes_json, genes = self.get_genes_per_panel(g_panel)
                g_panel.add_genes(genes)
                latest_panels.append(g_panel)
                self.panels.append(g_panel)  # append only if new
                time.sleep(1)  # sleep for 1 second to avoid DDoS safeguards
                self.save_panel(g_panel, self.json_dir)
            # TODO: Check for version differences in panels
            # if latest_panels and load_panels (self.panels) don't match, update the .json and .genes file/MySQL DB
        else:
            self.load_panels(self.json_dir)
        return self.panels
