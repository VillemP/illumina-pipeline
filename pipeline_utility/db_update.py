import os
import urllib
import urlparse

from jenkinsapi.jenkins import Jenkins
from progressbar import ProgressBar, Bar, Percentage

from pipeline_utility import file_utility
from pipeline_utility.baseconfig import BaseConfig
from pipeline_utility.jsonhandlerbase import JsonHandlerBase


class AnnotationDbCfg(BaseConfig):
    yaml_tag = u"!AnnotationDbCfg"

    def __init__(self, filepath, omim_api_key="", omim_host="api.omim.org",
                 hpo_host="http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/api/python",
                 directory="", hgncPath="~/hgnc_symbols.txt", hgnc_api="https://rest.genenames.org/fetch/",
                 hgnc_url="https://www.genenames.org/cgi-bin/"):
        super(AnnotationDbCfg, self).__init__(filepath)
        self.omim_api_key = omim_api_key
        self.omim_host = omim_host
        self.hpo_host = hpo_host
        self.directory = directory
        self.hpo_jenkins = None
        self.hgncPath = hgncPath
        self.hgnc_api = hgnc_api
        self.hgnc_url = hgnc_url


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
            # sys.stdout.flush()
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

    def downloadOMIMterms(self, fileName):
        filename, headers = urllib.urlretrieve(urlparse.urljoin(self.config.omim_host,
                                                                os.path.join(self.config.omim_api_key,
                                                                             fileName)),
                                               os.path.join(self.config.directory, fileName),
                                               self.show_progress)
        return filename

    def downloadHGNCsymbols(self):
        args = "download?col=gd_app_sym&col=gd_prev_sym&" \
               "col=gd_aliases&status=Approved&status_opt=2&" \
               "where=&order_by=gd_app_sym_sort&format=text&" \
               "limit=&hgnc_dbtag=on&submit=submit"

        response = urllib.urlopen(urlparse.urljoin(self.config.hgnc_url, args))

        with open(self.config.hgncPath, "wb+") as hgnc:
            for i, line in enumerate(response.readlines()):
                # Skip the header
                if i != 0:
                    hgnc.write(line)

    def splitOMIMterms(self, infile, outfile):
        with open(infile) as f:
            with open(outfile, "wb+") as out:
                for line in f.readlines():
                    if not line.startswith("#"):
                        line = line.split("\t")
                        if len(line) >= 13:
                            # OMIM gene symbols
                            genes = []
                            genes.extend([col.strip() for col in line[6].split(",")])

                            # HGNC gene symbol
                            if line[8] not in genes:
                                genes.append(line[8])
                            # Clear up any empty lines
                            genes = [g for g in genes if g != "" and g is not None]
                            for gene in genes:
                                # GENE<tab>Phenotype/disease
                                out.write("\t".join((gene, line[12])) + "\n")

    def splitHPOterms(self, infile, outfile):
        with open(infile) as f:
            with open(outfile, "wb+") as out:
                for line in f.readlines():
                    if not line.startswith("#"):
                        line = line.split("\t")
                        if len(line) >= 3:
                            out.write("\t".join((line[1], line[2])) + "\n")


def loadCfg(cfg_path):
    cfg = AnnotationDbCfg(cfg_path)
    cfg = cfg.load()

    return cfg


cfg = loadCfg(os.path.join(os.path.dirname(os.path.normpath(__file__)), "annotations.yaml"))
handler = AnnotationHandler(cfg)


def _update_db(annotationFolder, externalFilename, outfile, downloadedFiles, splithook):
    match = [f for f in downloadedFiles if f.strip() ==
             os.path.join(cfg.directory, externalFilename)]
    if len(match) > 0:
        newfile = os.path.join(annotationFolder, outfile)
        oldfile = os.path.join(annotationFolder, outfile + ".bak")
        # Make a backup of the old annotation file
        if os.path.exists(newfile):
            print "Creating backup: {0}".format(oldfile)
            os.rename(newfile, oldfile)
        if os.path.exists(match[0]):
            genes_to_phenotypes = match[0]
            splithook(genes_to_phenotypes, newfile)
            newlen = file_utility.file_len(newfile)
            newgenes = file_utility.count_unique_names(newfile, 0)

            oldlen = file_utility.file_len(oldfile)
            oldgenes = file_utility.count_unique_names(oldfile, 0)
            print("New {0} file contains {1} lines for {2} unique gene symbols.\n"
                  "The old file is saved as a backup in {3} and contained {4} lines for {5} unique genes.".
                  format(newfile, newlen, newgenes, oldfile, oldlen, oldgenes))
    else:
        print("ERROR: The update could not be started for {0}".format(outfile))


def update_omim(annotationFolder, test):
    if test:
        files = file_utility.find_filetype(annotationFolder, "txt")
        files = [f[1] for f in files]
    else:
        files = handler.downloadOMIMterms("genemap2.txt")
    _update_db(annotationFolder, "genemap2.txt", "gene.omim_disease_name.synonyms.txt", files,
               handler.splitOMIMterms)


def update_hpo(annotationFolder, test):
    if test:
        files = file_utility.find_filetype(annotationFolder, "txt")
        files = [f[1] for f in files]
    else:
        version, files = handler.download_hpo_terms(handler.load_hpo_index())
    _update_db(annotationFolder, "ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", "gene.hpoterm.txt", files,
               handler.splitHPOterms)


def updated_hgnc():
    print("Downloading HGNC symbols to {0}".format(handler.config.hgncPath))
    handler.downloadHGNCsymbols()
    print("Finished!")


def update_all(annotationFolder, test):
    updated_hgnc()
    update_hpo(annotationFolder, test)
    update_omim(annotationFolder, test)
