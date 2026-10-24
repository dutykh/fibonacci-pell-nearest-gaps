#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Run the thirteen fixed checks under per-process limits and compare each result with the record.

The checks, in order: the first seven rows; the universal identities; the finite Gaussian
windows (first verification) and their independent second verification; the algebra of the
tail bounds; the covers of the unbounded tails (first verification, replaying the recorded
trees) and their independent second verification; the six junctions and the named constants;
the labels of a row and the source comparisons up to the second junctions (first verification)
and their independent second verification; the two ladders below 3/5 (first verification); the
census of margins and the ladders again (second, independent verification); the finite facts on
directives, compared with the census of code/directives_census.jl (the second verification), with
controls of the Gaussian recursion along the Farey tree and of the square roots of -1.
The final record of each check is compared with data/verification-reference.json, so that any
change of result fails loudly.

Usage, from anywhere:
    python3 -B supplement/run_all.py                    # run and compare
    python3 -B supplement/run_all.py --write-reference  # run and record the results
The limits are generous on purpose (600 s of wall time, 300 s of processor time and 1 GiB of
address space per check); on an idle machine the thirteen checks take about a minute.
"""
from pathlib import Path
import json
import subprocess
import sys

if not __debug__ or sys.argv[1:] not in ([], ["--write-reference"]):
    raise SystemExit("Run without -O; the only accepted argument is --write-reference.")

try:
    import resource
except ImportError:
    resource = None

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "data" / "verification-reference.json"
WALL, CPU, MEMORY = 600, 300, 1 << 30
CHECKS = (
    ("early/check_early_rows.py", ()),
    ("checks/check_foundation_identities.py", ()),
    ("checks/verify_gaussian.py", ()),
    ("checks/independent_gaussian.py", ()),
    ("checks/check_tail_bounds.py", ()),
    ("tail/verify_all.py", ("--run",)),
    ("checks/independent_tail_cover.py", ()),
    ("checks/check_assembly.py", ()),
    ("checks/check_row_bounds.py", ()),
    ("checks/independent_row_bounds.py", ()),
    ("checks/check_margin_family.py", ()),
    ("checks/independent_margins.py", ()),
    ("checks/check_directives.py", ()),
)


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    resource.setrlimit(resource.RLIMIT_CPU, (CPU, CPU))


def output(value, stream):
    if value:
        if isinstance(value, bytes):
            value = value.decode(errors="replace")
        print(value.rstrip(), file=stream, flush=True)


def final_record(stdout):
    """The whole output if it is one JSON value, otherwise its last line read as JSON."""
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return json.loads(stdout.strip().splitlines()[-1])


def main():
    write = sys.argv[1:] == ["--write-reference"]
    missing = [name for name, _ in CHECKS if not (HERE / name).is_file()]
    if missing:
        raise SystemExit("Incomplete package; missing: " + ", ".join(missing))
    reference = {} if write else json.loads(REFERENCE.read_text())
    if not write and sorted(reference) != sorted(name for name, _ in CHECKS):
        raise SystemExit("The reference record does not list the thirteen checks; run with --write-reference.")
    if resource is None:
        print(f"Resource module unavailable: only the {WALL}-second wall limit applies.", flush=True)
    results = {}
    for name, arguments in CHECKS:
        command = [sys.executable, "-B", str(HERE / name), *arguments]
        print("Checking " + name + "...", flush=True)
        try:
            result = subprocess.run(
                command, cwd=HERE, text=True, capture_output=True, timeout=WALL,
                preexec_fn=limits if resource is not None else None,
            )
        except subprocess.TimeoutExpired as error:
            output(error.stdout, sys.stdout)
            output(error.stderr, sys.stderr)
            raise SystemExit(f"{name}: FAILED ({WALL}-second wall limit exceeded)") from None
        except (OSError, subprocess.SubprocessError) as error:
            raise SystemExit(name + ": could not start: " + str(error)) from None
        if result.returncode:
            output(result.stdout, sys.stdout)
            output(result.stderr, sys.stderr)
            detail = ("terminated by signal " + str(-result.returncode)
                      if result.returncode < 0 else "exit status " + str(result.returncode))
            raise SystemExit(name + ": FAILED (" + detail + ")")
        if not result.stdout.strip():
            raise SystemExit(name + ": FAILED (no check output)")
        output(result.stderr, sys.stderr)
        try:
            record = final_record(result.stdout)
        except json.JSONDecodeError:
            raise SystemExit(name + ": FAILED (the final record is not JSON)") from None
        if record.get("status") != "PASS":
            raise SystemExit(name + ": FAILED (status " + str(record.get("status")) + ")")
        if not write and record != reference[name]:
            print(json.dumps(record, sort_keys=True), flush=True)
            raise SystemExit(name + ": FAILED (the result differs from the reference record)")
        results[name] = record
        print(name + ": PASS" + ("" if write else ", identical to the reference record"), flush=True)
        print(json.dumps(record, sort_keys=True), flush=True)
    if write:
        REFERENCE.write_text(json.dumps(results, indent=1, sort_keys=True) + "\n")
        print("Reference record written for the thirteen checks.")
    else:
        print("All thirteen checks passed, each identical to its reference record.")


if __name__ == "__main__":
    main()
