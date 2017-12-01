FILTER_COLUMN_LIST = True
FILTER_COLUMN = False
FORMULA = '='
EDITVALUE = "edit"


class Ruleset:
    def __init__(self, column, python_notation, excel_notation, filtertype):
        # type: (str, str or None, list or str, bool) -> Ruleset

        self.column = column
        assert type(column) is str
        self.python_notation = python_notation
        self.excel_notation = excel_notation

        if filtertype == FILTER_COLUMN_LIST:
            assert type(excel_notation) is list

        # The filter type is used to determine the XlsxWriter
        self.filter_type = filtertype

    def __str__(self):
        return 'Column: {0}. Excel notation: {1}'.format(self.column, self.excel_notation)
