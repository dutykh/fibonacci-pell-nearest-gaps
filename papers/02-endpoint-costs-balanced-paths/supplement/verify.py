#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Reconstruct the norm-89 reduction graph and check the analytic identities."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import signal
import sys

if not __debug__:
    raise SystemExit("Optimized execution is forbidden: the assertions are part of the check.")

HERE = Path(__file__).resolve().parent


def resource_limits():
    for kind, requested in ((resource.RLIMIT_AS, 128 * 1024 * 1024),
                            (resource.RLIMIT_CPU, 30)):
        soft, hard = resource.getrlimit(kind)
        limit = requested
        if soft != resource.RLIM_INFINITY:
            limit = min(limit, soft)
        if hard != resource.RLIM_INFINITY:
            limit = min(limit, hard)
        resource.setrlimit(kind, (limit, limit))

    def expired(_signum, _frame):
        raise RuntimeError("INCONCLUSIVE: the fixed 30-second wall-time cap was reached")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-reference", action="store_true",
                        help="regenerate the full potential table and reference JSON after all checks pass")
    args = parser.parse_args()
    resource_limits()
    sys.path.insert(0, str(HERE / "code"))
    import endpoint
    import critical
    import witnesses

    endpoint_record, graph = endpoint.check()
    critical_record = critical.check(endpoint_record, graph)
    assert critical_record["L"] == "58" and critical_record["K_graph"] == "59"
    witness_record = witnesses.check()
    record = {"format_version": 1, "status": "PASS", "endpoint": endpoint_record,
              "critical": critical_record, "fixed_identities": witness_record,
              "balanced_count_bound": "min_theta |m-theta*n| < 59*(Delta+1)+2",
              "genuine_count_bound": "min_theta |P-theta*Q| <= 118*Delta+121"}
    # JSON normalization prevents tuple/list differences from masquerading as
    # a failed mathematical comparison on an ordinary read-only run.
    rendered = json.dumps(record, sort_keys=True, indent=2) + "\n"
    normalized = json.loads(rendered)
    reference = HERE / "data" / "reference-output.json"
    states, _index, _edges, _seeds, distance, _accepting = graph
    potential_table = [[list(states[i]), distance[i]]
                       for i in sorted(range(len(states)), key=states.__getitem__)]
    potential_file = HERE / "data" / "axis89-potential.json"
    if args.write_reference:
        reference.parent.mkdir(parents=True, exist_ok=True)
        reference.write_text(rendered, encoding="utf-8")
        # One complete state and integer potential per line is easier to inspect
        # than a fingerprint or a multi-megabyte indented expansion.
        table_text = "[\n" + ",\n".join("  " + json.dumps(row) for row in potential_table) + "\n]\n"
        potential_file.write_text(table_text, encoding="utf-8")
    else:
        if not reference.is_file():
            raise SystemExit("Missing data/reference-output.json; the package is incomplete.")
        if not potential_file.is_file():
            raise SystemExit("Missing data/axis89-potential.json; the integer table is incomplete.")
        if json.loads(potential_file.read_text(encoding="utf-8")) != potential_table:
            raise SystemExit("The stored state and potential table differs from the reconstructed one.")
        if json.loads(reference.read_text(encoding="utf-8")) != normalized:
            raise SystemExit("The recorded values differ from the reference file.")
    print(json.dumps(normalized, sort_keys=True))
    signal.alarm(0)


if __name__ == "__main__":
    main()
