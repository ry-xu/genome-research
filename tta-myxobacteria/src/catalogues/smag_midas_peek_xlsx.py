"""Print sheet names, header rows and shape of xlsx files (inspection helper for SMAG/MiDAS)."""
import sys
import openpyxl
for path in sys.argv[1:]:
    wb = openpyxl.load_workbook(path, read_only=True)
    print("=====", path)
    for ws in wb.worksheets:
        print("--- sheet", ws.title, ws.max_row, ws.max_column)
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            print([c for c in row][:25])
            if i >= 3:
                break
