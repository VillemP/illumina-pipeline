import json
import os

import requests

from pipeline_utility import file_utility


class JsonHandlerBase(object):
    def __init__(self):
        self.indices = list()

    def open_page(self, url, args=None):
        r = requests.get(url, args)
        # return r.json()
        return r

    def query(self, url, args=None):
        response = self.open_page(url, args)
        return response, response.json()

    def save_index(self, index_json, location, name):
        indexfile = os.path.join(location, name)
        self.indices.append(index_json)
        with open(indexfile, "wb+") as f:
            json.dump(index_json, f, sort_keys=True, indent=4)
        print("Updated index file {0}".format(name))

    def load_index(self, location, name):
        indexfile = file_utility.find_file(location, name)
        if indexfile[1] is not None:
            with open(indexfile[1]) as f:
                index_json = json.load(f)
            self.indices.append(index_json)
            return index_json
        else:
            return None
