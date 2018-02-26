import pipeline_utility.ruleset as filt
from pipeline_utility.ruleset import Ruleset


class MyeloidFilters(list):
    def __init__(self, sheets):
        super(MyeloidFilters, self).__init__([dict([])])
        for index, sheet in enumerate(sheets):
            self.append(dict([]))
            # expand the list of filters
        # self.filters = [dict([])]
        # Format is: sheet, column, python logic ('.' is unfiltered automatically), excel logic,
        # is it a list or a single filter

        self[0]['G'] = Ruleset('G', '{} >= 30',
                               'QUAL >= 30', filt.FILTER_COLUMN)  # QUAL > 30

        # noinspection PyTypeChecker
        self[0]['Q'] = Ruleset('Q', '{} <= 0.05',
                               '1000gAll <= 0.05 or 1000gAll == .',
                               filt.FILTER_COLUMN)

        self[0]['S'] = Ruleset('S',
                               '{} <= 0.05',
                               'ExacAll <= 0.05 or ExacAll == .',
                               filt.FILTER_COLUMN)

        self[0]['O'] = Ruleset('O', '"synonymous_SNV" != "{0}"',
                               ['.', 'Blanks',
                                'unknown',
                                'stopgain',
                                'frameshift_deletion',
                                'frameshift_insertion',
                                'nonframeshift_deletion',
                                'nonframeshift_insertion',
                                'nonsynonymous_SNV',
                                'stoploss'],
                               filt.FILTER_COLUMN_LIST)


class MyeloidPostprocess(list):
    def __init__(self, sheets):
        super(MyeloidPostprocess, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())

        # Prettify the column to be more human-readable.
        # Triple-escaped quotes ''' ''' ensure all characters in the string are escaped
        self[0]['AJ'] = Ruleset('AJ', '\'\'\'{0}\'\'\'.replace("_", " ")', None,
                                filt.EDITVALUE)  # Disease name from ANNO
        #        self[0]['AN'] = Ruleset('AN', '\'\'\'{0}\'\'\'.decode("ascii").replace("_", " ").replace(",", ", ")',
        #                                '\'=HYPERLINK("http://cancer.sanger.ac.uk/cosmic/mutation/overview?id={2}", "{0}")\'',
        #                                filt.FORMULA, '"".join(re.findall("ID=COSM(-?\d+(?:\d+)?)", "{0}"))') # Beautify cosmic data and convert to link
        self[0]['AN'] = Ruleset('AN', '\'\'\'{0}\'\'\'.decode("ascii").replace("_", " ").replace(",", ", ")',
                                None,
                                filt.EDITVALUE)  # Beautify cosmic data
        self[0]['AO'] = Ruleset('AO', '\'\'\'{0}\'\'\'.replace("_", " ")', None, filt.EDITVALUE)  # Disease.name
        self[0]['AQ'] = Ruleset('AQ', '\'\'\'{0}\'\'\'.replace("_", " ").replace(",", ", ")', None,
                                filt.EDITVALUE)  # HPO
        # Convert rs-values into links, does not check if links are valid. Empty values are skipped. {0} is cell.value,
        # {1} is cell.rawdata
        self[0]['D'] = Ruleset('D', None, '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{0}", "{1}")\'', filt.FORMULA)
        # Regex to get the number only. Old style dbSNP page
        # self[0]['C'] = Ruleset('C', '"".join(re.findall("\\d+", "{0}"))',
        #                       '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={0}", "{1}")\'',
        #                       filt.FORMULA)


class MyeloidFormats(list):
    def __init__(self, sheets):
        super(MyeloidFormats, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())

        # The hyperlink format for column C
        self[0]['D'] = {'font_color': 'blue', 'underline': 1}
        # self[0]['AN'] = {'font_color': 'blue'}
