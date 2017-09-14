import time

import requests

from TrusightOne.gene_panel import genepanel


class jsonhandler:
    def __init__(self):
        pass

    def open_page(self, url, args=None):
        r = requests.get(url, args)
        return r.json()

    def query(self, url, args=None):
        response = self.open_page(url, args)
        return response

    def get_genes_per_panel(self, panel_id):
        q = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/get_panel/"+panel_id)
        print("Got genes for panel {0}".format(panel_id))
        return q

    def get_all_panels(self, external=True):
        panels = list()
        if external:
            #TODO: Check for version differences in panels
            q = self.query("https://panelapp.extge.co.uk/crowdsourcing/WebServices/list_panels",
                           [{'format','json'}])
            print("Got all panels.")
            for panel in q['result']:
                g_panel = genepanel(panel['Name'], panel['Panel_Id'], panel, panel['CurrentVersion'])
                g_panel.add_genes(self.get_genes_per_panel(g_panel.panel_id))
                panels.append(g_panel)
                time.sleep(1) #sleep for 1 second to avoid DDoS safeguards
        else:
            #TODO:Open locally stored panel data
            pass
        return panels
