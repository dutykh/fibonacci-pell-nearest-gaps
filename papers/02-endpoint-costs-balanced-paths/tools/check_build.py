#!/usr/bin/env python3
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
#         Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Check the complete manuscript log, bibliography and LaTeX house style."""
from pathlib import Path
import re
import sys

if not __debug__ or len(sys.argv) > 2:
    raise SystemExit("Run without -O; usage: check_build.py [manuscript-stem]")

ROOT = Path(__file__).resolve().parents[1]
MAIN = sys.argv[1] if len(sys.argv) == 2 else "DD-LV-Endpoint-Costs-Balanced-Paths"
if not re.fullmatch(r"[A-Za-z0-9_-]+", MAIN):
    raise SystemExit("The manuscript stem must be a filename without an extension.")


def main():
    required = [ROOT / (MAIN + suffix) for suffix in (".tex", ".log", ".blg", ".pdf")]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Incomplete manuscript build; run make rebuild. Missing: "
                         + ", ".join(missing))
    failures = []
    log = (ROOT / (MAIN + ".log")).read_text()
    blg = (ROOT / (MAIN + ".blg")).read_text()
    warning = re.compile(
        r"LaTeX (?:Font )?Warning:|Package .* Warning:|Class .* Warning:|"
        r"(?:Underfull|Overfull) \\[hv]box|Missing character:|multiply defined|pdfTeX warning"
    )
    allowed_log = "pdfTeX warning (font expansion): font should be expanded before its first use"
    for line in log.splitlines():
        if line.startswith("!") or (warning.search(line) and allowed_log not in line):
            failures.append(line)
    if "Output written on " + MAIN + ".pdf" not in log:
        failures.append("The LaTeX log has no successful PDF-output marker.")
    for line in blg.splitlines():
        if line.startswith("Warning--"):
            failures.append(line)
        if re.search(r"\(There (?:was|were) .*error message|I couldn't open|"
                     r"I found no|Illegal, another|Repeated entry|You're missing", line):
            failures.append(line)
    sources = [ROOT / (MAIN + ".tex"), *sorted((ROOT / "sections").glob("*.tex")),
               *sorted((ROOT / "figures").glob("*.tex"))]
    for path in sources:
        for number, line in enumerate(path.read_text().splitlines(), 1):
            active = re.split(r"(?<!\\)%", line, maxsplit=1)[0]
            if re.search(r"\\(?:leq?|geq?)(?!slant|[A-Za-z])|(?<!\\)\\[\[\]()]", active):
                failures.append(f"{path.relative_to(ROOT)}:{number}: "
                                "disallowed math delimiter or inequality")
            if re.match(r"[A-Za-z]", active) and re.search(r"\S {2,}\S", active):
                failures.append(f"{path.relative_to(ROOT)}:{number}: repeated prose space")
    if failures:
        raise SystemExit("\n".join(failures))
    print("Complete LaTeX log, bibliography and source-style checks passed.")


if __name__ == "__main__":
    main()
