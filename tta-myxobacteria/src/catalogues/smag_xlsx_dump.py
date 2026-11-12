"""Dump every sheet of an xlsx to TSV via python-calamine. Usage: in.xlsx outprefix"""
import sys, csv
from python_calamine import CalamineWorkbook
wb = CalamineWorkbook.from_path(sys.argv[1])
for name in wb.sheet_names:
    rows = wb.get_sheet_by_name(name).to_python()
    out = f"{sys.argv[2]}.{name.replace(' ', '_').replace('/', '_')}.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for r in rows:
            w.writerow(["" if c is None else str(c).strip() for c in r])
    print(name, len(rows), out)
