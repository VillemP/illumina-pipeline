class Filter:
    def __init__(self, column, python_notation, excel_notation, is_list):
        # type: (str, str, list or str, bool) -> Filter

        self.column = column
        assert type(column) is str
        self.python_notation = python_notation

        if type(excel_notation) is not list:
            self.excel_notation = [excel_notation]  # convert to list

        else:
            self.excel_notation = excel_notation
        self.filter_type = is_list
        # assert self.filter_type is str or self.filter_type is int
        assert type(self.filter_type) is bool

    def __str__(self):
        return 'Column: {0}. Excel notation: {1}'.format(self.column, self.excel_notation)
