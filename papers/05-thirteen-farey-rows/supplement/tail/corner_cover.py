#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Exact matched-corner cover of one prescribed parity sector.

Default mode reports the size estimate and performs no cover. Explicit
--run requests the stated sector only. A cap, unresolved singleton, or
unexhausted stack is a failure, never a partial check of the domain.
"""

import sys

if not __debug__:
    raise SystemExit("assertions are required; do not run with -O")

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
import time

from exact import SOURCE_CEILINGS, compare, fixed_checks, ratio, source_limit


@dataclass(frozen=True)
class Domain:
    lo: int
    hi: int | None

    def __post_init__(self):
        assert self.lo >= 0
        assert self.hi is None or (self.hi >= self.lo and (self.hi-self.lo) % 2 == 0)


def roots(k, minimum_sum, sector):
    dimension = k+1
    assert 0 <= sector < 2**dimension
    parity = tuple((sector >> (dimension-1-i)) & 1 for i in range(dimension))
    threshold = (minimum_sum+dimension-1)//dimension
    assert threshold >= 2
    for first in range(dimension):
        box = []
        for index, bit in enumerate(parity):
            if index < first:
                hi = threshold-1-(threshold-1-bit) % 2
                assert hi >= bit
                box.append(Domain(bit, hi))
            elif index == first:
                box.append(Domain(threshold+(threshold-bit) % 2, None))
            else:
                box.append(Domain(bit, None))
        yield tuple(box)


def corner(box, maximize):
    result = []
    for index, domain in enumerate(box):
        increasing = index+1 < len(box) and domain.lo % 2 == 0
        use_upper = increasing == maximize
        result.append(domain.hi if use_upper else domain.lo)
    return tuple(result)


def split(box):
    choices = [(domain.lo, index) for index, domain in enumerate(box)
               if domain.hi is None or domain.hi > domain.lo]
    if not choices:
        raise AssertionError("unresolved eligible singleton: %r" % (box,))
    _, index = min(choices)
    domain = box[index]
    if domain.hi is None:
        cut = domain.lo+2
    else:
        count = (domain.hi-domain.lo)//2+1
        cut = domain.lo+2*(count//2)
    assert cut > domain.lo and (domain.hi is None or cut <= domain.hi)
    left, right = list(box), list(box)
    left[index] = Domain(domain.lo, cut-2)
    right[index] = Domain(cut, domain.hi)
    return tuple(left), tuple(right)


def cover(row, defects, minimum_cost, sector, node_cap, endpoint_cap,
          below_threshold=None):
    minimum_sum = minimum_cost-defects-1
    stack = list(roots(defects, minimum_sum, sector))
    root_count = len(stack)
    lower, upper = source_limit(row), SOURCE_CEILINGS[row]
    if below_threshold is not None:
        lower = upper = below_threshold
    counts = {"nodes": 0, "below": 0, "above": 0, "ineligible": 0,
              "max_endpoint": 0, "max_stack": len(stack)}
    events = sha256()
    while stack:
        if counts["nodes"] >= node_cap:
            raise RuntimeError("node cap reached with %d pending boxes" % len(stack))
        box = stack.pop()
        counts["nodes"] += 1
        maximum = max(max(d.lo, d.hi or 0) for d in box)
        if maximum > endpoint_cap:
            raise RuntimeError("gap endpoint cap reached")
        counts["max_endpoint"] = max(counts["max_endpoint"], maximum)
        key = [(d.lo, d.hi) for d in box]
        if all(d.hi is not None for d in box) and sum(d.hi for d in box) < minimum_sum:
            category = "ineligible"
        elif compare(ratio(corner(box, True)), lower) < 0:
            category = "below"
        elif compare(ratio(corner(box, False)), upper) > 0:
            category = "above"
        else:
            left, right = split(box)
            stack.extend((left, right))
            counts["max_stack"] = max(counts["max_stack"], len(stack))
            continue
        counts[category] += 1
        events.update(json.dumps([category, key], separators=(",", ":")).encode()+b"\n")
    leaves = counts["below"]+counts["above"]+counts["ineligible"]
    assert counts["nodes"] == 2*leaves-root_count
    result = {"status": "PASS", "row": row, "defects": defects,
              "minimum_cost": minimum_cost, "minimum_sum": minimum_sum,
              "sector": sector, "roots": root_count, **counts,
              "terminal_sha256": events.hexdigest()}
    if below_threshold is not None:
        if counts["above"]:
            raise AssertionError("a below-only threshold has an above leaf")
        result["below_threshold"] = {"numerator": list(below_threshold[0]),
                                     "denominator": below_threshold[1]}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row", type=int, choices=tuple(SOURCE_CEILINGS), required=True)
    parser.add_argument("--defects", type=int, choices=range(1, 9), required=True)
    parser.add_argument("--minimum-cost", type=int, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--node-cap", type=int, default=100000)
    parser.add_argument("--endpoint-cap", type=int, default=4096)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if args.node_cap < 1 or args.endpoint_cap < 1:
        parser.error("positive resource caps required")
    minimum_sum = args.minimum_cost-args.defects-1
    planned_roots = tuple(roots(args.defects, minimum_sum, args.sector))
    if not args.run:
        print(json.dumps({"status": "PREFLIGHT_ONLY", "row": args.row,
                          "defects": args.defects, "sector": args.sector,
                          "minimum_sum": minimum_sum, "roots": len(planned_roots),
                          "node_cap": args.node_cap, "endpoint_cap": args.endpoint_cap}, sort_keys=True))
        return
    started = time.monotonic()
    checks = fixed_checks()
    result = cover(args.row, args.defects, args.minimum_cost, args.sector,
                   args.node_cap, args.endpoint_cap)
    result["fixed_normalization_checks"] = checks
    result["elapsed_seconds"] = round(time.monotonic()-started, 6)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
