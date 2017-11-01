from pipeline_utility.filter import Filter


class TruesightOneFilters(list):
    def __init__(self, ac_5_percent):
        super(TruesightOneFilters, self).__init__([dict([])])
        # Format is: filters[sheet.index]['COL']
        # Filter (COL, python logic ('.' is unfiltered automatically), excel logic, is it listtype (several filtered))

        self[0]['F'] = Filter('F', '{} >= 30', 'QUAL >= 30', False)  # QUAL > 30
        self[0]['V'] = Filter('V', '{{0}} <= {}'.format(ac_5_percent),
                              'db_AC <= {} or db_AC == .'.format(ac_5_percent),
                              False)  # n(db_value)_AC <= 5% of alleles
        self[0]['S'] = Filter('S', '{} <= 0.01', 'ExacAll <= 0.01 or ExacAll == .', False)
        self[0]['Q'] = Filter('Q', '{} <= 0.01', '1000gAll <= 0.01 or 1000gAll == .', False)
        self[0]['M'] = Filter('M', '"{0}" == "exonic" or "{0}" == "splicing"', ['exonic', 'splicing'], True)
        self[0]['O'] = Filter('O', '"synonymous_SNV" != "{0}"', ['.', 'Blanks',
                                                                 'unknown',
                                                                 'stopgain',
                                                                 'frameshift_deletion',
                                                                 'frameshift_insertion',
                                                                 'nonframeshift_deletion',
                                                                 'nonframeshift_insertion',
                                                                 'nonsynonymous_SNV',
                                                                 'stoploss'], True)

        self[0]['L'] = Filter('L', '"{0}" != "no_annotation"', "GeneReq != no_annotation", False)