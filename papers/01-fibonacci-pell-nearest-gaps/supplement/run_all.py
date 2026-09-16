#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Run the exact standard-library checks used by the article."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


MANUSCRIPT = Path(__file__).resolve().parents[1]
CHECKS = Path(__file__).resolve().parent / "checks"

COMMANDS: list[tuple[str, ...]] = [
    ("check_two_arms.py",),
    ("check_logarithmic_h_bound.py",),
    ("check_effective_advances.py",),
    ("check_adjusted_logarithmic_form.py",),
    ("check_continued_fraction_closure.py",),
    ("check_terminal_reduction_second.py",),
    ("check_continued_fraction_closure_second.py",),
    ("check_nearest_gap_classification.py", "--q-max", "2000"),
    (
        "check_nearest_gap_classification_second.py",
        "--q-max",
        "2000",
    ),
    ("check_window_densities.py", "--q-max", "5000"),
    ("check_two_sign_odd_pell_orbit.py",),
    ("check_even_pell_anchor_orbits.py", "--q-max", "5000"),
    (
        "check_even_pell_anchor_orbits_second.py",
        "--q-max",
        "5000",
    ),
    ("check_general_seed_orbit.py",),
    ("check_general_seed_orbit_second.py",),
    ("check_uniform_local_clocks.py",),
    ("check_even_unit_gap_criterion.py", "--q-max", "2000"),
    ("check_two_clock_obstruction.py",),
    ("check_two_clock_obstruction_second.py",),
    ("check_plus_branch_quartic.py",),
    ("check_plus_branch_quartic_second.py",),
]


def main() -> int:
    for specification in COMMANDS:
        script = CHECKS / specification[0]
        if not script.is_file():
            print(f"MISSING: supplement/checks/{script.name}", file=sys.stderr)
            return 2
        command = [sys.executable, "-B", str(script), *specification[1:]]
        display = [
            "python3",
            "-B",
            f"supplement/checks/{script.name}",
            *specification[1:],
        ]
        print(f"\n$ {' '.join(display)}", flush=True)
        result = subprocess.run(command, cwd=MANUSCRIPT, check=False)
        if result.returncode:
            print(f"FAILED with exit status {result.returncode}", file=sys.stderr)
            return result.returncode
    print("\nPASS: every check completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
