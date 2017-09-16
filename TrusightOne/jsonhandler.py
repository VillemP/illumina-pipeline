import time

import requests

from TrusightOne.gene_panel import GenePanel


class JsonHandler:
    def __init__(self):
        pass

    def open_page(self, url, args=None):
        r = requests.get(url, args)
        #return r.json()
        return r

    def query(self, url, args=None):
        response = self.open_page(url, args)
        return response, response.json()

    def get_genes_per_panel(self, panel):
        q, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/get_panel/"+panel.panel_id)
        print("Got {0} genes for panel '{1}'".format(len(data['result']['Genes']), panel.name))
        return q, data

    def get_all_panels(self, external=True):
        panels = list()
        if external:
            #TODO: Check for version differences in panels
            json, data = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/list_panels",
                           [{'format','json'}])
            print("Got all panels.")
            for panel in data['result']:
                g_panel = GenePanel(panel['Name'], panel['Panel_Id'], data, panel['CurrentVersion'])
                genes_json, genes = self.get_genes_per_panel(g_panel)
                g_panel.add_genes(genes)
                panels.append(g_panel)
                time.sleep(1) #sleep for 10 seconds to avoid DDoS safeguards
        else:
            #TODO:Open locally stored panel data
            pass
        return panels
