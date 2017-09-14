from xlsxwriter.utility import xl_col_to_name


class Cell:
    def __init__(self, row_i, col_i, value, sheet):
        self.row_index = row_i
        self.column_index = col_i
        self.value = value
        self.sheet = sheet
        self.col_name = xl_col_to_name(self.column_index)
        self.is_toprow = self.row_index == 0
        self.filtered = False
        assert value is not None
        assert self.row_index >= 0 and self.column_index >= 0, \
            "Cell can only exist in positive range."

    @property
    def __str__(self):
        return "Cell: row={0}, column={1}, value={2}".format(self.row_index,
                                                             self.column_index, self.value)
