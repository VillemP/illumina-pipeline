import json
import os
import sys
import time
from distutils.version import LooseVersion as Version

from conda._vendor.progressbar import ProgressBar

import TrusightOne.gene_panel
import gene
from TrusightOne.gene_panel import GenePanel
from pipeline_utility import file_utility
# TODO: Currently presumes data is static and files will not go missing. Check for validity of JSONs
from pipeline_utility.jsonhandlerbase import JsonHandlerBase


class JsonHandler(JsonHandlerBase):
    def __init__(self, db_dir, config):
        super(JsonHandler, self).__init__()
        self.config = config
        self.panels = list()
        self.panel_index = self.load_index(db_dir)
        self.json_dir = db_dir
        self.pbar = None
        self.hgncHandler = gene.HgncHandler(config.hgncPath, config.tsoGenes, verbosity=gene.EXCEPTIONS)

    @property
    def loaded(self):
        return len(self.panels) > 0

    def get_genes_per_panel(self, panel):
        q, data = self.query("https://panelapp.genomicsengland.co.uk/WebServices/get_panel/" + panel.id)
        print("Downloaded {0} genes for panel '{1}'".format(len(data['result']['Genes']), panel.name))
        return q, data

    def save_index(self, index_json, location, name="panels.index"):
        return super(JsonHandler, self).save_index(index_json, location, name)

    def load_index(self, location, name="panels.index"):
        return super(JsonHandler, self).load_index(location, name)

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

    def load_panels(self, location, callback=None):
        print("Initializing panels...")
        start_time = time.time()
        panels = file_utility.find_filetype(location, ".json")
        len_panels = len(panels)
        if len(self.panels) > 0:
            # TODO: Enable post activation loading of panels.
            # Currently loading panels can only run at the start of the app.
            print("Warning! There are preexisting panels in the handler. Emptying the panels list.")
            self.panels[:] = []
        if len_panels > 1:
            for i, panel in enumerate(panels):
                file_location = panel[1]
                with open(file_location, "r") as f:
                    newpanel = GenePanel(self.hgncHandler, json.load(f), self.config)
                    genelist = file_utility.find_filetype(os.path.dirname(file_location), ".genes")
                    for genes in genelist:
                        with open(genes[1], "r") as gene_json:
                            newpanel.add_genes(json.load(gene_json))
                            # print "Loaded {0} genes into panel {1}".format(len(newpanel.genes), newpanel.name)
                    self.panels.append(newpanel)
                    if callback is not None:
                        callback(i + 1, len_panels)
            duration = time.time() - start_time
            print("Loaded {0} panels from local data in {1}s".format(len(self.panels), duration))
            return self.panels
        else:
            self.get_all_panels()

    def download_panel(self, g_panel):
        genes_json, genes = self.get_genes_per_panel(g_panel)
        g_panel.add_genes(genes)
        self.save_panel(g_panel, self.json_dir)

    def progressbar(self, current, total):
        if self.pbar is None:
            self.pbar = ProgressBar(maxval=total, fd=sys.stderr).start()
            self.pbar._time_sensitive = True

        if current < total:
            self.pbar.update(current)
            # sys.stdout.flush()
        else:
            self.pbar.finish()
            self.pbar = None

    def get_all_panels(self, external=True):
        if len(self.hgncHandler.hgnc_genes) == 0:
            self.hgncHandler.load_hgnc_genes(self.config.hgncPath)
        if len(self.hgncHandler.tso_genes) == 0:
            self.hgncHandler.load_tso_genes(self.config.tsoGenes)
        update_index = False
        if external:
            json_response, data = self.query("https://panelapp.genomicsengland.co.uk/WebServices/list_panels",
                                             [{'format', 'json'}])
            print("Found {0} panels from PanelApp API.".format(len(data['result'])))
            local_panels = file_utility.find_filetype(self.json_dir, ".json")
            if self.panel_index is None:
                print("Creating new index file panel.index")
                self.save_index(data, self.json_dir)
            for panel in data['result']:
                g_panel = GenePanel(self.hgncHandler, panel, self.config)

                local = [l for l in local_panels if g_panel.name == l[0].split('.')[0]]
                if len(local) == 0:
                    # Download panel if it doesn't exist. For example on first run.
                    self.download_panel(g_panel)
                else:
                    # Look for differences in local and external files
                    # time.sleep(1)  # sleep for 1 second to avoid DDoS safeguards
                    match = [_oldpanel for _oldpanel in self.panel_index['result']
                             if g_panel.id == _oldpanel['Panel_Id']]
                    # Should return only 1 panel with matching ID
                    if len(match) == 1:
                        # Creates a temporary panel object from json for comparison with the new panel object
                        # Serves the same purpose as self.load_panels(), latest_panels comparing which is slower
                        old_panel = GenePanel(self.hgncHandler, match[0])
                        # Download data for only new panels based on panel version
                        if TrusightOne.gene_panel.compare_versions(g_panel, old_panel):
                            print("Found a panel {0} with new version ({1} vs {2})."
                                  .format(g_panel.name, g_panel.version, old_panel.version))
                            self.download_panel(g_panel)
                            update_index = True

                    elif len(match) < 1:
                        # It must be a new panel, it's not in the current index
                        print("Found a new panel: {0}".format(g_panel.name))
                        self.download_panel(g_panel)
                        self.save_panel(g_panel, self.json_dir)
                    else:
                        raise LookupError("There are panels with overlapping IDs (panels: {})! "
                                          "There can't be more than one panel with a unique Id!"
                                          .format([matched_panel.name.encode("ascii") for matched_panel in match]))
            # Finally load panels with new updates included
            if update_index:
                self.save_index(data, self.json_dir)
            self.load_panels(self.json_dir, self.progressbar)

        else:
            self.load_panels(self.json_dir, self.progressbar)
        return self.panels

    def write_version_one_panels(self, update=False):
        if not self.loaded or update:
            self.get_all_panels(True)
        # Selecting and sorting for tables above version 1.0, ordered by diseasegroup and subgroup
        currentlist = [p for p in self.panels if p.version >= Version("1.0")]
        currentlist.sort(key=lambda panel: (panel.diseasegroup, panel.diseasesubgroup))

        print("Writing the table for panels with versions >=1.0 to {}".format(
            os.path.join(self.config.json_dir, "query.txt")))
        with open(os.path.join(self.config.json_dir, "query.txt"), "w+") as f:
            f.write(TrusightOne.gene_panel.table_format + '\n')
            for panel in currentlist:
                f.write(panel.as_table + '\n')
