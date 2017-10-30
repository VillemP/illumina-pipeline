"""
Simple script, that takes n text files (tab delimited) and outputs excel workbook with n sheets. 

Usage: python TSO_txttoxlsx.py input.txt input2.txt ... inputn.txt output.xlsx

"""

import re
import sys

import xlsxwriter
from xlsxwriter.utility import xl_col_to_name

from pipeline_utility.filter import Filter


class Cell:
    def __init__(self, row_i, col_i, value, sheet):
        self.row_index = row_i
        self.column_index = col_i
        self.value = value
        self.sheet = sheet  # 0-indexed
        self.col_name = xl_col_to_name(self.column_index)
        self.is_toprow = self.row_index == 0
        self.filtered = False
        assert value is not None
        assert self.row_index >= 0 and self.column_index >= 0

    @property
    def __str__(self):
        return "Cell: row={0}, column={1}, value={2}".format(self.row_index,
                                                             self.column_index, self.value)

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


def check_filter(cell, filter):
    # type: (Cell) -> bool
    if cell.is_toprow:
        return False
        # print cell

    if filter is not None and cell.value != '.':
        # check whether filter exists for a column
        # value {100} >= 30 --> true --> must be filtered? return False.
        return not eval(filter.python_notation.format(cell.value))
    return False


def create_excel(name, filterset, files, ac_5_percent=None):
    wbook = xlsxwriter.Workbook()
    wbook.filename = name
    print("Starting excel file {0} with {1} input files.".format(name, len(files)))

    try:
        ws = []
        for i, infile in enumerate(files):
            ws.append(("worksheet%d", i))
            print "\nStarting input file " + str(i + 1) + ": " + infile
            ws[i] = wbook.add_worksheet()
            data = []
            # filters.append(dict([])) # expand the filters according to sheets, each sheet has a dict, filters[sheet]

            table_width = 0
            with open(infile, "r") as f:
                for line in f.readlines():
                    l = line.split()
                    table_width = len(l)
                    data.append(l)
                table_length = len(data)

                assert table_length > 0
                assert table_width > 0

                # Activate autofilter (mandatory) for any sheets containing filters
                if len(filterset) >= i:
                    if len(filterset[i]) > 0:
                        ws[i].autofilter(0, 0, table_length, table_width)
                if i == 0:

                    # Terrible workaround for local DB sample count, if the number isn't provided in args.
                    # Uses regex to get the count from the header
                    if ac_5_percent is None:
                        ac_5_percent = current_db(Cell(0, 21, data[0][21], 0))
                        filterset[0] = filterset[0]['V'] = Filter('V', '{{0}} <= {}'.format(ac_5_percent),
                                                                  'db_AC <= {} or db_AC == .'.format(ac_5_percent),
                                                                  False)  # n(db_value)_AC <= 5% of alleles

                    add_filters(ws[0], filterset)

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

                        filter = filterset[cell.sheet].get(cell.col_name, None)
                        # Hide False values (filter conditions not met)
                        if check_filter(cell, filter):
                            ws[i].set_row(cell.row_index, options={'hidden': True})
                            row_filtered = True
                        col_index += 1

                    if row_filtered:
                        filtered += 1  # count filtered rows

                    row_index += 1
                    col_index = 0  # start from the beginning of the row
                    row_filtered = False

                print "Total rows: {0}\nFiltered rows: {1}".format(table_length, filtered)
                # print "Filters active {0}.\n{1}".format(len(filters[i]),
                #                                       '\n'.join(str(v) for v in filters[i].itervalues()))
    except (Exception) as error:
        print "Error in creating excel file {0}".format(wbook.filename)
        raise error
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

    # TODO: Enable filtering from script mode
    # Currently will not activate any filters from script, filters must be created manually e.g.
    # filters=TruesightOneFilters(current_db)
    filters = [dict([])]
    create_excel(sys.argv[-1], filters, sys.argv[1:-1], 0)
    print "Done with excel file {}".format(sys.argv[-1])
