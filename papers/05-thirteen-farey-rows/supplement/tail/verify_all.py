#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Recompute all 28 unbounded tail covers and compare their actual records.

Default: size estimate only. --run: full exact verification, with no file writes.
The manuscript's six rows, cost thresholds and defect ranges are fixed
below; a modified or incomplete manifest cannot change the covered domain.
"""

import sys
if not __debug__:
    raise SystemExit("assertions are required; do not run with -O")

import argparse
from hashlib import sha256
import json
from pathlib import Path

from corner_cover import cover
from exact import fixed_checks

CASES = ((8, 34, (2, 3)),
         (9, 41, (2, 3, 4)),
         (10, 46, (2, 3, 4, 5)),
         (11, 84, (2, 3, 4, 5, 6)),
         (12, 110, (2, 3, 4, 5, 6, 7)),
         (13, 156, (2, 3, 4, 5, 6, 7, 8)))
DATA = Path(__file__).resolve().parent/"data"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="store_true")
    args = p.parse_args()
    if not args.run:
        print(json.dumps({"status": "PREFLIGHT_ONLY", "groups": 28,
                          "sectors": 1976, "known_nodes": 1023768,
                          "per_sector_node_cap": 100000,
                          "per_group_node_cap": 1000000,
                          "external_limits": "128 MiB, 30 seconds"}, sort_keys=True))
        return
    fixed_checks()
    all_totals = {key: 0 for key in ("sectors", "nodes", "roots", "below", "above", "ineligible")}
    manifest = sha256()
    for row, cost, counts in CASES:
        for k in counts:
            group_nodes = 0
            path = DATA/f"row{row}-k{k}.jsonl"
            with path.open() as source:
                for sector in range(2**(k+1)):
                    remaining = 1000000-group_nodes
                    if remaining < 1:
                        raise RuntimeError("per-group node cap exceeded")
                    record = cover(row, k, cost, sector, min(100000, remaining), 4096)
                    if k == row-5 and record["above"]:
                        raise AssertionError("high defect layer has an above leaf")
                    saved = source.readline()
                    if not saved or json.loads(saved) != record:
                        raise AssertionError("manifest mismatch: %s sector %d" % (path.name, sector))
                    encoded = json.dumps(record, sort_keys=True, separators=(",", ":"))+"\n"
                    manifest.update(encoded.encode())
                    group_nodes += record["nodes"]
                    all_totals["sectors"] += 1
                    for key in all_totals:
                        if key != "sectors":
                            all_totals[key] += record[key]
                if source.read():
                    raise AssertionError("extra manifest records: %s" % path.name)
            print(json.dumps({"row": row, "defects": k, "nodes": group_nodes,
                              "status": "PASS"}), flush=True)
    with (DATA/"lower-cost-two-defects.jsonl").open() as source:
        for sector in range(8):
            record = cover(13, 2, 155, sector, 100000, 4096,
                           below_threshold=((71, 0), 100))
            saved = source.readline()
            if not saved or json.loads(saved) != record:
                raise AssertionError("lower-cost manifest mismatch")
            encoded = json.dumps(record, sort_keys=True, separators=(",", ":"))+"\n"
            manifest.update(encoded.encode())
            all_totals["sectors"] += 1
            for key in all_totals:
                if key != "sectors":
                    all_totals[key] += record[key]
        if source.read():
            raise AssertionError("extra lower-cost records")
    assert all_totals == {"sectors": 1976, "nodes": 1023768, "roots": 14280,
                          "below": 342676, "above": 176348, "ineligible": 0}
    print(json.dumps({"status": "PASS", "groups": 28, **all_totals,
                      "combined_manifest_sha256": manifest.hexdigest()}, sort_keys=True))


if __name__ == "__main__":
    main()
