"""
Simple script, that takes n text files (tab delimited) and outputs excel workbook with n sheets. 

Usage: python TSO_txttoxlsx.py input.txt input2.txt ... inputn.txt output.xlsx

"""

import re
import sys

import xlsxwriter
from xlsxwriter.utility import xl_col_to_name

from miniseq.miniseq_excel_filters import MiniseqFilters


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
        assert self.row_index >= 0 and self.column_index >= 0

    @property
    def __str__(self):
        return "Cell: row={0}, column={1}, value={2}".format(self.row_index,
                                                             self.column_index, self.value)


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


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def add_filters(sheet, filters):
    for key, filter in filters[sheet.index].items():
        assert key == filter.column
        if not filter.filter_type:
            sheet.filter_column(filter.column, filter.excel_notation[0])
        elif filter.filter_type:
            sheet.filter_column_list(filter.column, filter.excel_notation)
        else:
            raise ValueError("Filter type not valid. Expected boolean. \n{0} type: {1}"
                             .format(filter, type(filter.filter_type)))


# noinspection PyTypeChecker
def current_db(cell):
    # type: (Cell) -> int
    assert cell.value is not None
    assert type(cell.value) is str
    ac_5_percent = int(''.join(re.findall('\d+', cell.value)))  # regex to find the current db number
    # noinspection PyTypeChecker
    return int(round(ac_5_percent * 0.1))


    # Returns the boolean on whether the cell's value should be hidden
    # True for hidden, false for shown


def check_filter(cell):
    # type: (Cell) -> bool
    filter = filters[cell.sheet].get(cell.col_name, None)
    if cell.is_toprow:
        return False
        # print cell

    if filter is not None and cell.value != '.':
        # check whether filter exists for a column
        # value {100} >= 30 --> true --> must be filtered? return False.
        return not eval(filter.python_notation.format(cell.value))
    return False


def create_excel(name, filters, files):
    wbook = xlsxwriter.Workbook()
    wbook.filename = name
    try:
        ws = []
        for i, infile in enumerate(files):
            ws.append(("worksheet%d", i))
            print "\nStarting input file " + str(i + 1) + ": " + infile
            ws[i] = wbook.add_worksheet()
            data = []
            filters.append(dict([]))

            table_width = 0
            f = open(sys.argv[i + 1], "r")
            for line in f.readlines():
                l = line.split()
                table_width = len(l)
                data.append(l)
            f.close()
            table_length = len(data)

            if i == 0:
                assert table_length > 0
                assert table_width > 0
                main_sheet = ws[0]  # First worksheet containing variants
                main_sheet.autofilter(0, 0, table_length, table_width)

                ac_5_percent = current_db(Cell(0, 21, data[0][21], 0))

                # Format is: sheet, column, python logic ('.' is unfiltered automatically), excel logic,
                # is it a list or a single filter

                # noinspection PyTypeChecker
                filters[0]['F'] = Filter('F', '{} >= 30', 'QUAL >= 30', False)  # QUAL > 30
                # noinspection PyTypeChecker
                filters[0]['V'] = Filter('V', '{{0}} <= {}'.format(ac_5_percent),
                                         'db_AC <= {} or db_AC == .'.format(ac_5_percent),
                                         False)  # n(db_value)_AC <= 5% of alleles
                filters[0]['S'] = Filter('S', '{} <= 0.01', 'ExacAll <= 0.01 or ExacAll == .', False)
                # noinspection PyTypeChecker
                filters[0]['Q'] = Filter('Q', '{} <= 0.01', '1000gAll <= 0.01 or 1000gAll == .', False)
                filters[0]['M'] = Filter('M', '"{0}" == "exonic" or "{0}" == "splicing"', ['exonic', 'splicing'], True)
                filters[0]['O'] = Filter('O', '"synonymous_SNV" != "{0}"', ['.', 'Blanks',
                                                                            'unknown',
                                                                            'stopgain',
                                                                            'frameshift_deletion',
                                                                            'frameshift_insertion',
                                                                            'nonframeshift_deletion',
                                                                            'nonframeshift_insertion',
                                                                            'nonsynonymous_SNV',
                                                                            'stoploss'], True)

                filters[0]['L'] = Filter('L', '"{0}" != "no_annotation"', "GeneReq != no_annotation", False)
                add_filters(main_sheet, filters)
                # filters['O'] = Filter(main_sheet, 'O', '"{0}" is not "synonymous_snv"','funcReg != *synonymous', int)
                # main_sheet.freeze_panes(0, table_width)

            row_index = 0
            # noinspection PyRedeclaration
            col_index = 0
            filtered = 0
            # noinspection PyRedeclaration
            row_filtered = False

            if i != 1:
                # freeze top row for sheets 1, 3, 4
                ws[i].freeze_panes(1, 0)

            for row in data:
                for col in row:
                    cell = Cell(row_index, col_index, col, i)
                    if is_number(cell.value):
                        ws[i].write_number(cell.row_index, cell.column_index, num(cell.value))
                    else:
                        ws[i].write(cell.row_index, cell.column_index, cell.value)
                    # Hide False values (filter conditions not met)
                    if check_filter(cell):
                        ws[i].set_row(cell.row_index, options={'hidden': True})
                        row_filtered = True
                    col_index += 1

                if row_filtered:
                    filtered += 1  # count filtered rows

                row_index += 1
                col_index = 0  # start from the beginning of the row
                row_filtered = False

            print "Total variants: {0}\nFiltered variants: {1}".format(table_length, filtered)
            # print "Filters active {0}.\n{1}".format(len(filters[i]),
            #                                       '\n'.join(str(v) for v in filters[i].itervalues()))
    except (Exception) as error:
        print "Error in creating excel file {0}".format(wbook.filename)
        print error
    finally:
        wbook.close()


if __name__ == "__main__":
    reload(sys)
    sys.setdefaultencoding('utf8')
    inputnr = len(sys.argv) - 2

    print "\nNr of inputfiles: " + str(inputnr)

    for i in range(0, inputnr):
        print "Input file " + str(i + 1) + ": " + sys.argv[i + 1]

    print "\nOutput is written to: " + sys.argv[-1]

    # TODO: Replace current_db requirement
    filters = MiniseqFilters(5)
    create_excel(sys.argv[-1], filters)
    print "Done with excel file {}".format(sys.argv[-1])
