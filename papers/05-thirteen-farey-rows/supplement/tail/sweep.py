#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Stream every parity sector of one exact tail check.

With no --run/--verify, report a size estimate only. --run writes a new JSONL
manifest after every sector has closed. --verify recomputes every covering
tree and compares every record, rather than trusting a saved digest.
Run with external time and memory limits; each group is capped and an exceeded cap is a failure, never a partial result.
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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--row", type=int, choices=(8, 9, 10, 11, 12, 13), required=True)
    p.add_argument("--defects", type=int, choices=range(1, 9), required=True)
    p.add_argument("--minimum-cost", type=int, required=True)
    p.add_argument("--node-cap", type=int, default=100000)
    p.add_argument("--total-node-cap", type=int, default=1000000)
    p.add_argument("--endpoint-cap", type=int, default=4096)
    p.add_argument("--require-below", action="store_true")
    action = p.add_mutually_exclusive_group()
    action.add_argument("--run", type=Path, metavar="OUTPUT_JSONL")
    action.add_argument("--verify", type=Path, metavar="SAVED_JSONL")
    args = p.parse_args()
    if min(args.node_cap, args.total_node_cap, args.endpoint_cap) <= 0:
        p.error("positive caps required")
    sectors = 2**(args.defects+1)
    if not (args.run or args.verify):
        print(json.dumps({"status": "PREFLIGHT_ONLY", "sectors": sectors,
                          "row": args.row, "defects": args.defects,
                          "minimum_cost": args.minimum_cost,
                          "per_sector_node_cap": args.node_cap,
                          "total_node_cap": args.total_node_cap}, sort_keys=True))
        return
    fixed_checks()
    destination = args.run
    pending = destination.with_suffix(destination.suffix+".partial") if destination else None
    if destination:
        destination.parent.mkdir(parents=True, exist_ok=True)
    source = args.verify.open() if args.verify else None
    output = pending.open("w") if pending else None
    totals = {name: 0 for name in ("nodes", "below", "above", "ineligible", "roots")}
    manifest = sha256()
    try:
        for sector in range(sectors):
            remaining = args.total_node_cap-totals["nodes"]
            if remaining < 1:
                raise RuntimeError("aggregate node cap exhausted")
            record = cover(args.row, args.defects, args.minimum_cost, sector,
                           min(args.node_cap, remaining), args.endpoint_cap)
            if args.require_below and record["above"]:
                raise AssertionError("a below-only group has an above leaf")
            if source:
                line = source.readline()
                if not line or json.loads(line) != record:
                    raise AssertionError("saved sector mismatch: %d" % sector)
            encoded = json.dumps(record, sort_keys=True, separators=(",", ":"))+"\n"
            if output:
                output.write(encoded)
                output.flush()
            manifest.update(encoded.encode())
            for name in totals:
                totals[name] += record[name]
            if (sector+1) % 16 == 0:
                print(json.dumps({"checkpoint_sectors": sector+1, "nodes": totals["nodes"]}), flush=True)
        if source and source.read():
            raise AssertionError("unexpected extra manifest records")
    finally:
        if source:
            source.close()
        if output:
            output.close()
    assert totals["nodes"] == 2*(totals["below"]+totals["above"]+totals["ineligible"])-totals["roots"]
    if pending:
        pending.replace(destination)
    print(json.dumps({"status": "PASS", "row": args.row, "defects": args.defects,
                      "minimum_cost": args.minimum_cost, "sectors": sectors,
                      **totals, "manifest_sha256": manifest.hexdigest()}, sort_keys=True))


if __name__ == "__main__":
    main()
