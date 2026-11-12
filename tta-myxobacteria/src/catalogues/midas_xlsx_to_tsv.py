"""Dump one xlsx sheet to TSV (skipping leading '#' title rows). Usage: in.xlsx out.tsv"""
import sys, csv
import openpyxl
wb = openpyxl.load_workbook(sys.argv[1], read_only=True)
ws = wb.worksheets[0]
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    for row in ws.iter_rows(values_only=True):
        if row and isinstance(row[0], str) and row[0].startswith("#"):
            continue
        if all(c is None for c in row):
            continue
        w.writerow(["" if c is None else str(c).strip() for c in row])
