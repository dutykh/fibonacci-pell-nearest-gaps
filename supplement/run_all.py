#!/usr/bin/env python3
#
# Fibonacci-Pell nearest gaps: an all-exponent classification and
# quadratic-unit orbit rigidity
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
"""Run the exact standard-library certificates used by the manuscript."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


MANUSCRIPT = Path(__file__).resolve().parents[1]
CERTIFICATES = Path(__file__).resolve().parent / "certificates"

COMMANDS: list[tuple[str, ...]] = [
    ("check_c0082_logarithmic_h_bound_curator.py",),
    ("check_c0082_effective_advances_independent_audit_agent_ea.py",),
    ("check_c0082_full8_adjusted_form_agent_c.py",),
    ("check_c0082_full8_dp_closure_agent_c.py",),
    ("check_c0082_full8_complete_exclusion_independent_audit_agent_fra.py",),
    ("check_full8_adjusted_dp_independent_agent_b.py",),
    ("check_all_exponent_nearest_gap_agent_aeng.py", "--q-max", "2000"),
    (
        "check_all_exponent_nearest_gap_independent_audit_agent_aei.py",
        "--q-max",
        "2000",
    ),
    ("check_candidate_window_density_curator.py", "--q-max", "5000"),
    ("check_two_sign_odd_pell_orbit_independent_audit_agent_tsa.py",),
    ("check_even_pell_anchor_orbits_agent_epa.py", "--q-max", "5000"),
    (
        "check_even_pell_anchor_orbits_independent_audit_agent_eaa.py",
        "--q-max",
        "5000",
    ),
    ("check_general_seed_orbit_nearest_gap_curator.py",),
    ("check_general_seed_orbit_nearest_gap_independent_audit_agent_gsa.py",),
    ("check_uniform_local_clock_generalisation_agent_lcg.py",),
    ("check_c0068_nearest_even_unit_gap_curator.py", "--q-max", "2000"),
    ("check_c0077_imprimitive_two_clock_recurrence_agent_cte.py",),
    ("check_c0077_imprimitive_two_clock_recurrence_independent_audit_agent_cia.py",),
    ("check_c0068_plus_j1_quartic_curator.py",),
    ("check_j1_plus_branch_quartic_independent_audit_agent_j1.py",),
]


def main() -> int:
    for specification in COMMANDS:
        script = CERTIFICATES / specification[0]
        if not script.is_file():
            print(f"MISSING: supplement/certificates/{script.name}", file=sys.stderr)
            return 2
        command = [sys.executable, "-B", str(script), *specification[1:]]
        display = [
            "python3",
            "-B",
            f"supplement/certificates/{script.name}",
            *specification[1:],
        ]
        print(f"\n$ {' '.join(display)}", flush=True)
        result = subprocess.run(command, cwd=MANUSCRIPT, check=False)
        if result.returncode:
            print(f"FAILED with exit status {result.returncode}", file=sys.stderr)
            return result.returncode
    print("\nPASS: every manuscript certificate completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
