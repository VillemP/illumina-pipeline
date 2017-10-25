from pipeline_utility.filter import Filter


class MiniseqFilters(list):
    def __init__(self, sheets):
        super(MiniseqFilters, self).__init__([dict([])])
        for index, sheet in enumerate(sheets):
            self.append(dict([]))
            # expand the list of filters
        # self.filters = [dict([])]
        # Format is: sheet, column, python logic ('.' is unfiltered automatically), excel logic,
        # is it a list or a single filter

        self[0]['F'] = Filter('F', '{} >= 30',
                              'QUAL >= 30', False)  # QUAL > 30

        # noinspection PyTypeChecker
        self[0]['P'] = Filter('P', '{} <= 0.05',
                              '1000gAll <= 0.05 or 1000gAll == .',
                              False)

        self[0]['R'] = Filter('R',
                              '{} <= 0.05',
                              'ExacAll <= 0.05 or ExacAll == .',
                              False)

        self[0]['N'] = Filter('N', '"synonymous_SNV" != "{0}"',
                              ['.', 'Blanks',
                               'unknown',
                               'stopgain',
                               'frameshift_deletion',
                               'frameshift_insertion',
                               'nonframeshift_deletion',
                               'nonframeshift_insertion',
                               'nonsynonymous_SNV',
                               'stoploss'],
                              True)
