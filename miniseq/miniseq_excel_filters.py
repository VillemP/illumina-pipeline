import pipeline_utility.ruleset as filt
from pipeline_utility.ruleset import Ruleset


class MiniseqFilters(list):
    def __init__(self, sheets):
        super(MiniseqFilters, self).__init__([dict([])])
        for index, sheet in enumerate(sheets):
            self.append(dict([]))
            # expand the list of filters
        # self.filters = [dict([])]
        # Format is: sheet, column, python logic ('.' is unfiltered automatically), excel logic,
        # is it a list or a single filter

        self[0]['F'] = Ruleset('F', '{} >= 30',
                               'QUAL >= 30', filt.FILTER_COLUMN)  # QUAL > 30

        # noinspection PyTypeChecker
        self[0]['P'] = Ruleset('P', '{} <= 0.05',
                              '1000gAll <= 0.05 or 1000gAll == .',
                               filt.FILTER_COLUMN)

        self[0]['R'] = Ruleset('R',
                              '{} <= 0.05',
                              'ExacAll <= 0.05 or ExacAll == .',
                               filt.FILTER_COLUMN)

        self[0]['N'] = Ruleset('N', '"synonymous_SNV" != "{0}"',
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


class MiniseqPostprocess(list):
    def __init__(self, sheets):
        super(MiniseqPostprocess, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())

        # Prettify the column to be more human-readable.
        # Triple-escaped quotes ''' ''' ensure all characters in the string are escaped
        self[0]['AI'] = Ruleset('AI', '\'\'\'{0}\'\'\'.replace("_", " ")', None,
                                filt.EDITVALUE)  # Disease name from ANNO
        self[0]['AM'] = Ruleset('AM', '\'\'\'{0}\'\'\'.replace("_", " ")', None, filt.EDITVALUE)  # Disease.name
        self[0]['AO'] = Ruleset('AO', '\'\'\'{0}\'\'\'.replace("_", " ").replace(",", ", ")', None,
                                filt.EDITVALUE)  # HPO
        # Convert rs-values into links, does not check if links are valid. Empty values are skipped.
        self[0]['C'] = Ruleset('C', None, '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{0}", "{1}")\'', filt.FORMULA)
        # Regex to get the number only. Old style dbSNP page
        # self[0]['C'] = Ruleset('C', '"".join(re.findall("\\d+", "{0}"))',
        #                       '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={0}", "{1}")\'',
        #                       filt.FORMULA)


class MiniseqFormats(list):
    def __init__(self, sheets):
        super(MiniseqFormats, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())

        # The hyperlink format for column C
        self[0]['C'] = {'font_color': 'blue', 'underline': 1}
