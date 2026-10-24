#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Optional bounded preparation of the prescribed Gaussian data, one fibre at a time."""
import argparse
from math import gcd
from pathlib import Path
import resource
import subprocess
import sys

if not __debug__:
    raise SystemExit("Assertions are required; omit -O.")

HERE = Path(__file__).resolve().parent
WINDOWS = {8: 28, 9: 34, 10: 38, 11: 75, 12: 100, 13: 145}


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (128*1024*1024,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (30, 30))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", required=True,
                        help="explicitly prepare the entire fixed finite window")
    parser.add_argument("--resume", action="store_true", help="skip already saved source records")
    args = parser.parse_args()
    for p, bound in WINDOWS.items():
        for q in range(p+1, bound+1):
            if gcd(p, q) != 1:
                continue
            target = HERE.parent / "data" / "gaussian" / f"{p:02d}-{q:03d}.json"
            if args.resume and target.exists():
                continue
            command = [sys.executable, "-B", str(HERE/"export_gaussian_data.py"), str(p), str(q)]
            try:
                result = subprocess.run(command, text=True, capture_output=True,
                                        timeout=30, preexec_fn=limits)
            except subprocess.TimeoutExpired as error:
                print(f"CAP REACHED P={p} Q={q}", flush=True)
                print(error.stdout.decode() if isinstance(error.stdout, bytes) else error.stdout)
                raise SystemExit(1)
            if result.returncode:
                print(result.stdout)
                print(result.stderr)
                raise SystemExit(result.returncode)
            print(result.stdout.strip().splitlines()[-1], flush=True)
    print("Prepared the fixed Gaussian window; independent verification is still required.")


if __name__ == "__main__":
    main()
