"""
A script that takes n+1 text files (tab delimited) and outputs excel workbook with n sheets.
Can be used to post-process data before writing and activate custom filters on the data.
The functions embedded rely on the pipeline_utility package and thus can't be used as a single script.

This package can be used to post-process data before writing and activate custom filters on the data.
Use the create_excel() function with the input/output args and a list of dictionaries [dict([])] containing Filters.
Each dictionary in the list equates to a sheet, each value equates to a filter or postprocess (pipelineutility.ruleset)
or format (dict), KEY is the column name (A, B, C..).
Currently supports one format, filter and postprocessing per cell per column.

Structure of the formats, filter and postprocessing lists:

FormatsList[dict()]
     <index> 0 (sheet 1): Dict{'column_name<A>':Dict{'bold':True, 'font_color':'blue'...}}
     # dicts correspond to Xlsxwriter Format workbook.addformat(dict<>)
     # or wb.addformat(formats_list[sheet][KEY]

FilterList[dict()]
    [0]=Dict{'column_name<A>':Filter('A', args), 'column_name<B>':Filter('B', args, filter.SINGLE_FILTER), n+1}
    [0]=Dict{'column_name<B>':Filter('B', args), 'column_name<C>':Filter('C', args, filter.SINGLE_FILTER), n+1}

PostprocessingList[dict()]
    <index> 0 (sheet 1): Dict{'column_name<A>':Filter('A', args), 'column_name<B>':Filter('B', args), n+1}
    <index> 1 (sheet 2): Dict{'column_name<A>':Filter('A', args), 'column_name<B>':Filter('B', args), n+1}
    <index> 2 (sheet 3): Dict{'column_name<A>':Filter('A', args), 'column_name<B>':Filter('B', args), n+1}
    [3]={'AN':Filter('AN', '"{0}.replace(" ", "_")"', None, filter.EDIT)}
    [3]={'AN':Filter('AN', None, '"=SQRT({0})"', filter.FORMULA)}
    [3]={'CN':Filter('CN', '"{0}.replace(" ", "_")"', '\'=HYPERLINK("http://localhost:10000/request={0}","{1}")\'',
        filter.FORMULA)} # in excel_notation, two args can be formatted {0} for cell.value and {1} for cell.rawdata

#Usage: python txttoxlsx_filtered.py input.txt input2.txt ... inputn.txt output.xlsx
"""

import re
import sys
import traceback

import xlsxwriter
from xlsxwriter.utility import xl_col_to_name

import pipeline_utility.ruleset as filt
from TrusightOne.tso_excel_filters import TruesightOnePostprocess, TruesightOneFormats

EMPTY_CELL = (None, "Blank", ".", "")


class Cell:
    def __init__(self, row_i, col_i, value, sheet):
        self.row_index = row_i
        self.column_index = col_i
        self._rawdata = value  # The raw data is the data read in from the input file and won't be processed
        self.value = value  # The value can be post-processed before writing to the table
        self.displaydata = None  # If this data is set to other than None, this data will be displayed over cell.value
        self.sheet = sheet  # 0-indexed
        self.col_name = xl_col_to_name(self.column_index)
        self.is_toprow = self.row_index == 0
        self.filtered = False
        self.type = get_type(value)
        self.formula = None
        self.visuals = None
        self.is_blank = self.value in EMPTY_CELL
        self.format = None
        assert self.row_index >= 0 and self.column_index >= 0

    @property
    def rawdata(self):
        return self._rawdata
    @property
    def __str__(self):
        return "Cell: row={0}, column={1}, value={2}".format(self.row_index,
                                                             self.column_index, self.value)


class Formula(object):
    def __init__(self, value):
        super(Formula, self).__init__()
        self.str = value
        assert self.str.startswith("="), "Tried to create a formula with incorrect first character: {0}".format(self)


def get_type(s):
    if s is not None:
        if is_number(s):
            return int
        else:
            return str


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
    # TODO: Check for AttributeError in Rulesets?
    for key, filter in filters[sheet.index].items():
        assert key == filter.column
        if filter.filter_type == filt.FILTER_COLUMN:
            sheet.filter_column(filter.column, filter.excel_notation)
        elif filter.filter_type == filt.FILTER_COLUMN_LIST:
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


def check_format(cell, formats_list, formats, wbook):
    currentformat = None
    # skip empty cells
    # If you want formats for empty cells, post-process the value or change the EMPTY_CELLS tuple
    if cell.value not in EMPTY_CELL:
        # Try and get a format for this column, otherwise return None
        if cell.format is None:
            cell.format = formats_list[cell.sheet].get(cell.col_name, None)
            if cell.format is not None:
                # Add a new format object to be used as a format
                currentformat = formats[cell.col_name] = wbook.add_format(cell.format)
            else:
                currentformat = None
    return currentformat

def check_filter(cell, filter):
    # type: (Cell) -> bool
    if cell.is_toprow:
        return False
        # print cell

    # Empty cells will not be used for filtering
    if filter is not None and cell.value not in EMPTY_CELL:
        # check whether filter exists for a column
        # value {100} >= 30 --> true --> must be filtered? return False.
        return not eval(filter.python_notation.format(cell.value))
    return False


def post_process(cell, ruleset):
    """
    Uses the Filter class to post-process values using python eval (e.g. {0}.replace('.', '_')
    where {0} is the current cell.value), create a formula independent of the cell value with python_notation as None
    or both with the filter containing both python_notation and excel_notation.

    Skips EMPTY_CELL values for formulas so they have to be changed  beforehand.

    Examples:
    Filter('A', '{0}.replace('.', "DEFAULT")', '=HYPERLINK("http://localhost:10000/show?request={0}")', filter.FORMULA)
    result is '=HYPERLINK("http://localhost:10000/show?request={0}"'.format(cell.value)) # cell value is processed and changed
    Filter('A', None, '=HYPERLINK("http://localhost:10000/show?request={0}")', filter.FORMULA)
    result is '=HYPERLINK("http://localhost:10000/show?request={1}"'.format(cell._rawdata)) # (None == no processing)

    Runs the Filter.python_notation only, if the filter exists, otherwise return the unchanged cell.

    The Filter class sets the post processing ruleset to a whole column.

    :param cell: (txttoxlsx_filtered.Cell)
    :param ruleset: The ruleset object to be run against the current cell
    :return: Returns the same cell object with a changed or unchanged cell.value. The cell value can be confidently
    written to the excel table. However, if the cell.formula has now been set,
    the XlsxWriter.worksheet.formula_writer() will be used instead.
    """
    if ruleset is not None:
        if not cell.is_toprow:
            try:
                # First run any post processing on the value of the cell (e.g. data selection or cleanup,
                # removing unwanted characters)
                if ruleset.python_notation is not None:
                    cell.value = eval(ruleset.python_notation.format(cell.value))
                # Secondly check whether any special display data is used added (e.g. cell.value has been converted but
                # something other should be displayed
                if ruleset.display_python_notation is not None:
                    cell.displaydata = eval(ruleset.display_python_notation.format(cell.value))
                # Thirdly try to create a valid excel formula using the new data
                # if the ruleset is of type filter.FORMULA, other filter types' excel notation is used in add_filter()
                # Skip cells that are empty. If empty cells need to be converted to formulas the value has to be changed
                if ruleset.filter_type == filt.FORMULA:
                    if ruleset.excel_notation is not None and len(ruleset.excel_notation) > 0:
                        if cell.value not in EMPTY_CELL:
                            cell.formula = eval(
                                ruleset.excel_notation.format(cell.value, cell.rawdata, cell.displaydata))
                    else:
                        raise ValueError("The filter {0} had an empty formula in the excel_notation argument!"
                                         .format(ruleset))
            except(SyntaxError) as error:
                sys.stderr.write("Your filter contains a syntax error in the python or excel notation!\n"
                                 "Filter:{0}".format(ruleset))
                traceback.print_exc(file=sys.stderr)
                raise error
    return cell


def create_excel(name, files, filterset=None, postprocess_list=None, formats_list=None):
    if postprocess_list is None:
        postprocess_list = [dict()]
    if formats_list == None:
        formats_list = [dict()]
    if filterset is None:
        filterset = [dict([])]
    formats = dict()

    reload(sys)
    sys.setdefaultencoding('utf8')
    wbook = xlsxwriter.Workbook()
    wbook.filename = name
    print("Starting excel file {0} with {1} input files.".format(name, len(files)))

    try:
        ws = []
        for i, infile in enumerate(files):
            ws.append(("worksheet%d" % i))
            print("Starting input file " + str(i + 1) + ": " + infile)
            ws[i] = wbook.add_worksheet()
            data = []

            table_width = 0
            with open(infile, "r") as f:
                for line in f.readlines():
                    l = line.split()
                    table_width = len(l)
                    data.append(l)
                table_length = len(data)

                assert table_length > 0, "Table cannot be empty."
                assert table_width > 0

                # Activate autofilter (mandatory) for any sheets containing filters
                if len(filterset) >= i:
                    if len(filterset[i]) > 0:
                        ws[i].autofilter(0, 0, table_length, table_width)
                if i == 0:
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
                        # Try and get a ruleset for this column, otherwise return None
                        post = postprocess_list[cell.sheet].get(cell.col_name, None)
                        cell = post_process(cell, post)
                        # Try and get a format for this column, otherwise return None
                        currentformat = check_format(cell, formats_list, formats, wbook)
                        if cell.formula is not None:
                            ws[i].write_formula(cell.row_index, cell.column_index, str(cell.formula),
                                                currentformat)
                        elif cell.type is int:
                            # Ensure a number format
                            ws[i].write_number(cell.row_index, cell.column_index, num(cell.value),
                                               currentformat)
                        elif cell.type is str:
                            # Write a pretty string format
                            ws[i].write_string(cell.row_index, cell.column_index, cell.value, currentformat)
                        else:
                            # Try to handle all other formats (datetime, etc)
                            ws[i].write(cell.row_index, cell.column_index, cell.value, currentformat)

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

                print("Total rows: {0}\nFiltered rows: {1}".format(table_length, filtered))
                # print "Filters active {0}.\n{1}".format(len(filters[i]),
                #                                       '\n'.join(str(v) for v in filters[i].itervalues()))
    except (Exception) as error:
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write("Error in creating excel file {0}: {1}\n".format(wbook.filename, error))
        raise error
    finally:
        print ("Closing excel file...")
        wbook.close()


if __name__ == "__main__":
    inputnr = len(sys.argv) - 2

    print("Nr of inputfiles: " + str(inputnr))

    for i in range(0, inputnr):
        print "Input file " + str(i + 1) + ": " + sys.argv[i + 1]

    print("Output is written to: " + sys.argv[-1])

    # TODO: Enable filtering from script mode
    # Currently will not activate any filters from script, filters must be created manually e.g.
    # filters=TruesightOneFilters(current_db)
    filters = [dict([])]
    post = TruesightOnePostprocess(sys.argv[1:-1])
    formats = TruesightOneFormats(sys.argv[1:-1])
    create_excel(sys.argv[-1], sys.argv[1:-1], filters)
    print("Done with excel file {}".format(sys.argv[-1]))
