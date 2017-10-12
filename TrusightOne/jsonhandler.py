import json
import os

import requests

import TrusightOne.gene_panel
from TrusightOne.gene_panel import GenePanel
from pipeline_utility import file_utility


# TODO: Currently presumes data is static and files will not go missing. Check for validity of JSONs
class JsonHandler:
    def __init__(self, db_dir):
        self.panels = list()
        self.panel_index = self.load_index(db_dir)
        self.json_dir = db_dir

    def open_page(self, url, args=None):
        r = requests.get(url, args)
        # return r.json()
        return r

    def query(self, url, args=None):
        response = self.open_page(url, args)
        return response, response.json()

    def get_genes_per_panel(self, panel):
        q, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/get_panel/" + panel.panel_id)
        print("Downloaded {0} genes for panel '{1}'".format(len(data['result']['Genes']), panel.name))
        return q, data

    @staticmethod
    def save_index(handler, index_json, location):
        indexfile = os.path.join(location, "panels.index")
        handler.panel_index = index_json
        with open(indexfile, "wb+") as f:
            json.dump(index_json, f, sort_keys=True, indent=4)
        print("Updated index file panels.index")

    @staticmethod
    def load_index(location):
        indexfile = file_utility.find_file(location, "panels.index")
        if indexfile is not None:
            with open(indexfile[1]) as f:
                index_json = json.load(f)
            return index_json
        else:
            return None

    def save_panel(self, panel, location):
        # Create directory for panel (JSON_dir/panel/panel.json
        # TODO: save into MySQL?
        if not os.path.exists(os.path.join(location, panel.name)):
            os.makedirs(os.path.join(location, panel.name))
        with open(os.path.join(os.path.join(location, panel.name),
                               panel.name + ".json"), "wb+") as json_file:
            json.dump(panel.json, json_file, sort_keys=True, indent=4)

        # Creates directory and file for panel's gene (JSON_dir/panel/panel.genes)
        # Panel.genes is also a json
        # TODO: Join the panel json and genes json or load directly into MySQL
        with open(os.path.join(os.path.join(location, panel.name),
                               panel.name + ".genes"), "wb+") as genes_file:
            json.dump(panel.genes_json, genes_file, sort_keys=True, indent=4)

    def load_panels(self, location):
        panels = file_utility.find_filetype(location, ".json")
        if len(self.panels) > 0:
            # TODO: Enable post activation loading of panels.
            # Currently loading panels can only run at the start of the app.
            print("Warning! There are preexisting panels in the handler. Emptying the panels queue.")
            self.panels[:] = []
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

    def get_all_panels(self, external=True):
        update_index = False
        if external:
            json_response, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/list_panels",
                                             [{'format', 'json'}])
            print("Got {0} panels from PanelApp API.".format(len(data['result'])))
            if self.panel_index is None:
                print("Creating new index file panel.index")
                self.save_index(self, data, self.json_dir)
            for panel in data['result']:
                g_panel = GenePanel(panel)

                # TODO: Too complicated? Replace with load all panels every time, only save ones which differ?
                # time.sleep(1)  # sleep for 1 second to avoid DDoS safeguards
                match = [_oldpanel for _oldpanel in self.panel_index['result']
                         if g_panel.panel_id == _oldpanel['Panel_Id']]
                # Should return only 1 panel with matching ID
                if len(match) == 1:
                    # Creates a temporary panel object from json for comparison with the new panel object
                    # Serves the same purpose as self.load_panels(), latest_panels comparing which is slower
                    old_panel = GenePanel(match[0])
                    # Download data for only new panels based on panel version
                    if TrusightOne.gene_panel.compare_versions(g_panel, old_panel):
                        print("Found a panel {0} with new version ({1} vs {2})."
                              .format(g_panel.name, g_panel.version, old_panel.version))
                        genes_json, genes = self.get_genes_per_panel(g_panel)
                        g_panel.add_genes(genes)
                        # latest_panels.append(g_panel)
                        self.save_panel(g_panel, self.json_dir)
                        update_index = True
                    else:
                        # print("{} {} {}".format(g_panel.name, g_panel.version, old_panel.version))
                        pass
                elif len(match) < 1:
                    raise LookupError("No matching panel found for {}".format(g_panel))
                else:
                    raise LookupError("There are panels with overlapping IDs (panels: {})! Aborting!"
                                      .format([matched_panel.name.encode("ascii") for matched_panel in match]))
            # Finally load panels with new updates included
            if update_index:
                self.save_index(self, data, self.json_dir)
            self.load_panels(self.json_dir)

        else:
            self.load_panels(self.json_dir)
        return self.panels
