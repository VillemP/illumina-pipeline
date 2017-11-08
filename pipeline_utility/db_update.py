import os
import sys
import urllib

from conda._vendor.progressbar import ProgressBar, Bar, Percentage
from jenkinsapi.jenkins import Jenkins

from pipeline_utility import file_utility
from pipeline_utility.baseconfig import BaseConfig
from pipeline_utility.jsonhandlerbase import JsonHandlerBase


class HPOURLopener(urllib.FancyURLopener):
    version = "Chrome/56.0.2924.87"


class AnnotationDbCfg(BaseConfig):
    yaml_tag = u"!AnnotationDbCfg"

    def __init__(self, filepath, omim_auth_key="", omim_host="api.omim.org",
                 hpo_host="http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/api/python",
                 directory=""):
        super(AnnotationDbCfg, self).__init__(filepath)
        self.omim_auth_key = omim_auth_key
        self.omim_host = omim_host
        self.hpo_host = hpo_host
        self.directory = directory
        self.hpo_jenkins = None


class AnnotationHandler(JsonHandlerBase):
    def __init__(self, config):
        super(AnnotationHandler, self).__init__()
        self.config = config
        if not os.path.exists(self.config.directory):
            os.makedirs(self.config.directory)
        self.hpo_handler = Jenkins(self.config.hpo_host)

        self.pbar = None

    def load_hpo_index(self):
        data = self.hpo_handler.get_data(self.config.hpo_host + "api/python")
        return data

    def show_progress(self, block_num, block_size, total_size):
        if self.pbar is None:
            self.pbar = ProgressBar(maxval=total_size, widgets=[Bar('=', '[', ']'), ' ', Percentage()])
            self.pbar._time_sensitive = True
            self.pbar.start()

        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
            sys.stdout.flush()
        else:
            self.pbar.finish()
            self.pbar = None

    def download_hpo_terms(self, data):
        filelist = list()
        # urllib._urlopener = HPOURLopener()
        if 'artifacts' in data:
            for artifact in data['artifacts']:
                filename, headers = urllib.urlretrieve(self.config.hpo_host + 'artifact/' + artifact['relativePath'],
                                                       os.path.join(self.config.directory, artifact['fileName']),
                                                       self.show_progress)
                filelist.append(filename)
        print("Finished download: {0}".format(filelist))
        return data['id'], filelist

    def load_omim_index(self):
        name = "omim_index.txt"
        index = self.load_index(self.config.directory, name)
        if index:
            return index
        else:
            return self.save_index(self.query(self.config.omim_host)[1], self.config.directory, name)

    def splitHPOterms(self, infile, outfile):
        with open(infile) as f:
            with open(outfile, "wb+") as out:
                line = f.readline().split("\t")
                if len(line) >= 3:
                    out.write("\t".join((line[1], line[2])) + "\n")


def loadCfg(cfg_path):
    cfg = AnnotationDbCfg(cfg_path)
    cfg = cfg.load()

    return cfg


def update_all(annotationFolder):
    cfg = loadCfg(os.path.join(os.path.dirname(os.path.normpath(__file__)), "annotations.yaml"))
    handler = AnnotationHandler(cfg)
    version, files = handler.download_hpo_terms(handler.load_hpo_index())
    # Check if the gene_phenotypes file was among the downloaded files
    match = [f for f in files if f.strip().lower() ==
             os.path.join(cfg.directory, "ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt".lower())]
    if len(match) > 0:
        newfile = os.path.join(annotationFolder, "gene.hpoterm.txt")
        oldfile = os.path.join(annotationFolder, "gene.hpoterm.OLD.txt")
        # Make a backup of the old phenotypes file
        if os.path.exists(newfile):
            os.rename(newfile, oldfile)
        genes_to_phenotypes = file_utility.find_file(cfg.directory, match[0])[1]
        handler.splitHPOterms(genes_to_phenotypes, newfile)
        newlen = file_utility.file_len(newfile)
        newgenes = file_utility.count_unique_names(newfile, 1)

        oldlen = file_utility.file_len(oldfile)
        oldgenes = file_utility.count_unique_names(newfile, 1)
        print("New {0} file contains {1} lines for {2} unique gene names.\n"
              "The old file is saved as a backup in {3} and contained {4} lines for {5} unique genes.".
              format(newfile, newlen, newgenes, oldfile, oldlen, oldgenes))
    print version
    print files
