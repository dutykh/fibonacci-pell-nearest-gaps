#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Reproduce all exact checks, with per-process time limits."""
from pathlib import Path
import subprocess
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed runner accepts no arguments.")

try:
    import resource
except ImportError:
    resource = None


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (128*1024*1024, 128*1024*1024))
    resource.setrlimit(resource.RLIMIT_CPU, (30, 30))


here = Path(__file__).resolve().parent
commands = []
for script in ("check_decoder.py", "check_examples.py", "check_displaced.py", "check_centres.py",
               "check_nested_and_lengths.py", "check_centre_decoder.py"):
    commands.append([sys.executable, "-B", str(here / "checks" / script)])
for command in commands:
    result = subprocess.run(command, cwd=here, text=True, capture_output=True,
                            timeout=30, preexec_fn=limits if resource else None)
    if result.returncode:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
    print(" ".join(command[2:]) + ": PASS")
    print(result.stdout.strip().splitlines()[-1])
print("All checks reproduced successfully.")
