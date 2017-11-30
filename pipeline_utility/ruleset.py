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

        if filtertype == FILTER_COLUMN_LIST:
            if type(excel_notation) is not list:
                self.excel_notation = [excel_notation]  # convert to list
        else:
            self.excel_notation = excel_notation
        # The filter type is used to determine the XlsxWriter
        self.filter_type = filtertype

    def __str__(self):
        return 'Column: {0}. Excel notation: {1}'.format(self.column, self.excel_notation)
