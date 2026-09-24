#!/usr/bin/env python3
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
#         Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Check the complete manuscript log, bibliography, and LaTeX house style."""
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[1]
main = sys.argv[1] if len(sys.argv) == 2 else "DD-LV-Effective-Markov-Traces"
failures = []
log = (root / (main + ".log")).read_text()
blg = (root / (main + ".blg")).read_text()
warning = re.compile(r"LaTeX (?:Font )?Warning:|Package .* Warning:|Class .* Warning:|"
                     r"(?:Underfull|Overfull) \\[hv]box|Missing character:|multiply defined|pdfTeX warning")
allowed_log = "pdfTeX warning (font expansion): font should be expanded before its first use"
for line in log.splitlines():
    if warning.search(line) and allowed_log not in line:
        failures.append(line)
allowed_bib = {"Warning--missing pages in " + key for key in
               ("AltassanLuca2021", "LagisquetEtAl2021", "LeeEtAl2023")}
for line in blg.splitlines():
    if line.startswith("Warning--") and line not in allowed_bib:
        failures.append(line)
for path in [root / (main + ".tex"), *sorted((root / "sections").glob("*.tex")),
             *sorted((root / "figures").glob("*.tex"))]:
    for number, line in enumerate(path.read_text().splitlines(), 1):
        active = re.split(r"(?<!\\)%", line, maxsplit=1)[0]
        if re.search(r"\\(?:leq?|geq?)(?!slant|[A-Za-z])|(?<!\\)\\[\[\]()]", active):
            failures.append(f"{path.relative_to(root)}:{number}: disallowed math delimiter or inequality")
        if re.match(r"[A-Za-z]", active) and re.search(r"\S {2,}\S", active):
            failures.append(f"{path.relative_to(root)}:{number}: repeated prose space")
if failures:
    print("\n".join(failures))
    raise SystemExit(1)
print("Complete LaTeX log, bibliography, and source-style checks passed.")
