import pipeline_utility.ruleset as filt
from pipeline_utility.ruleset import Ruleset


# TODO: Switch eval(str) with lambda arguments
class TruesightOneFilters(list):
    def __init__(self, ac_5_percent, sheets):
        super(TruesightOneFilters, self).__init__([dict([])])
        self._ac_5_percent = ac_5_percent

        # Format is: filters[sheet.index]['COL']
        # Filter (COL, python logic ('.' is unfiltered automatically), excel logic, is it listtype (several filtered))
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict([]))
        self[0]['F'] = Ruleset('F', '{} >= 30', 'QUAL >= 30', filt.FILTER_COLUMN)  # QUAL > 30
        self[0]['V'] = Ruleset('V', '{{0}} <= {}'.format(self._ac_5_percent),
                              'db_AC <= {} or db_AC == .'.format(self._ac_5_percent),
                               filt.FILTER_COLUMN)  # n(db_value)_AC <= 5% of alleles
        self[0]['S'] = Ruleset('S', '{} <= 0.01', 'ExacAll <= 0.01 or ExacAll == .', filt.FILTER_COLUMN)
        self[0]['Q'] = Ruleset('Q', '{} <= 0.01', '1000gAll <= 0.01 or 1000gAll == .', filt.FILTER_COLUMN)
        self[0]['M'] = Ruleset('M', '"{0}" == "exonic" or "{0}" == "splicing"', ['exonic', 'splicing'],
                               filt.FILTER_COLUMN_LIST)
        self[0]['O'] = Ruleset('O', '"synonymous_SNV" != "{0}"', ['.', 'Blanks',
                                                                 'unknown',
                                                                 'stopgain',
                                                                 'frameshift_deletion',
                                                                 'frameshift_insertion',
                                                                 'nonframeshift_deletion',
                                                                 'nonframeshift_insertion',
                                                                 'nonsynonymous_SNV',
                                                                  'stoploss'], filt.FILTER_COLUMN_LIST)

        self[0]['L'] = Ruleset('L', '"{0}" != "no_annotation"', "GeneReq != no_annotation", filt.FILTER_COLUMN)

    @property
    def ac_5_percent(self):
        return self._ac_5_percent


class TruesightOnePostprocess(list):
    def __init__(self, sheets):
        super(TruesightOnePostprocess, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())
        # Prettify the column to be more human-readable.
        # Triple-escaped quotes ''' ''' ensure all characters in the string are escaped
        self[0]['P'] = Ruleset('P', '\'\'\'{0}\'\'\'.replace(",", ", ")', None, filt.EDITVALUE)
        self[0]['AN'] = Ruleset('AN', '\'\'\'{0}\'\'\'.replace("_", " ")', None, filt.EDITVALUE)
        self[0]['AP'] = Ruleset('AP', '\'\'\'{0}\'\'\'.replace("_", " ").replace(",", ", ")', None, filt.EDITVALUE)
        # Convert rs-values into links, does not check if links are valid. Empty values are skipped.
        self[0]['C'] = Ruleset('C', None, '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{0}", "{1}")\'', filt.FORMULA)
        # Regex to get the number only. This can be commented out to use the old-style dbSNP page
        # self[0]['C'] = Ruleset('C', '"".join(re.findall("\\d+", "{0}"))',
        #                      '\'=HYPERLINK("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={0}", "{1}")\'',
        #                       filt.FORMULA)


class TruesightOneFormats(list):
    def __init__(self, sheets):
        super(TruesightOneFormats, self).__init__([dict()])
        # Accomodate the filters container to have room for each sheet
        for index, sheet in enumerate(sheets):
            self.append(dict())

        # The hyperlink format for column C
        self[0]['C'] = {'font_color': 'blue', 'underline': 1}
