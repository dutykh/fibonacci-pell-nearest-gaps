#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Reproduce one fixed small-occurrence check with exact matrices.

Python 3.10+, standard library. Run from any directory, for example:
  python3 -B checks/classify_small_occurrences.py 3
Add --write to regenerate its small reference JSON. The full finite domain is
justified in the paper. This program
does not identify tree occurrences by their numerical trace.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Assertions must be enabled; do not use python -O.")
SUPPLEMENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SUPPLEMENT / "code"))
import trace_decision as td

BIT_CAP = 8192
OUTPUT_CAP = 16 * 1024 * 1024
MODULUS = 1 << 256


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest_integer(value):
    return int.from_bytes(hashlib.sha256(encode(value).encode()).digest(), "big")


def classify(h):
    bound = td.bounds(h)
    assert bound["comparison_bound"] <= 871_000
    matches = []
    result = {"h": h, "bound": bound, "branches": {}}
    patterns = [t for t in range(1, h) if gcd(h, t) == 1]
    for branch, height, indices, qhi in (
        ("positive", bound["positive_height"], range(1, 3*h), 3*h//2),
        ("minimum", bound["minimum_height"], (0,), 5*h//2),
    ):
        seeds = comparisons = words = max_bits = seed_sum = comparison_sum = 0
        branch_matches = 0
        for seed in td.iter_seeds(height):
            seeds += 1
            assert seed.p <= height
            assert seed.x <= seed.y and seed.x <= 3*seed.p*seed.x-seed.y
            assert seed.p**2+seed.x**2+seed.y**2 == 3*seed.p*seed.x*seed.y
            assert td.tr(seed.A) == 3*seed.x and td.tr(seed.B) == 3*seed.p
            assert td.tr(td.mul(seed.A, seed.B)) == 3*seed.y
            seed_sum += digest_integer([branch, seed.path, seed.orientation,
                                        seed.p, seed.x, seed.y])
            ray_cache = {}
            for l in indices:
                for t in patterns:
                    target_matrix = td.word_matrix(seed, h, t, l)
                    target = td.tr(target_matrix)
                    assert target > 2 and target % 3 == 0
                    assert target_matrix[0]*target_matrix[3]-target_matrix[1]*target_matrix[2] == 1
                    words += 1
                    for q in range(-((h+1)//2), qhi+1):
                        m = h*l+t+q
                        if m < 1:
                            continue
                        for sign in (-1, 1):
                            n = sign*m
                            if n not in ray_cache:
                                ray_cache[n] = td.ray_trace(seed, n)
                            single = ray_cache[n]
                            assert single > 2 and single % 3 == 0
                            max_bits = max(max_bits, abs(target).bit_length(), abs(single).bit_length())
                            assert max_bits <= BIT_CAP
                            comparisons += 1
                            comparison_sum += digest_integer([branch, seed.path,
                                seed.orientation, l, t, q, n, target, single])
                            if target != single:
                                continue
                            v = td.linear(seed.alpha, 1, seed.beta, n)
                            w = td.linear(seed.alpha, h, seed.beta, h*l+t)
                            difference = l-n
                            k = difference if difference >= 1 else -difference-1 if difference <= -2 else None
                            original_t = t if difference >= 1 else h-t if difference <= -2 else None
                            step = 1 if difference >= 1 else -1
                            witness = td.isometry_witness(v, w)
                            assert gcd(*v) == gcd(*w) == 1
                            assert witness is not None, ("nonisometric equality", h, seed, l, t, n)
                            row = {"branch": branch, "h": h, "path": seed.path,
                                "orientation": seed.orientation, "seed": [seed.p, seed.x, seed.y],
                                "l": l, "t": t, "q": q, "n": n, "trace": target,
                                "source": list(v), "target": list(w),
                                "intersection": abs(td.det(v, w)), "original_k": k,
                                "original_t": original_t, "isometry": witness,
                                "character": [single, 3*seed.p, td.ray_trace(seed, n+step)]}
                            if k is not None:
                                assert row["intersection"] == h*k+original_t
                            matches.append(row)
                            branch_matches += 1
            assert seeds <= bound[branch+"_seed_bound"]
        assert comparisons <= bound[branch+"_comparison_bound"]
        result["branches"][branch] = {"seeds": seeds, "literal_words": words,
            "comparisons": comparisons, "equality_rows": branch_matches,
            "max_trace_bits": max_bits,
            "seed_digest_sum": f"{seed_sum % MODULUS:064x}",
            "comparison_digest_sum": f"{comparison_sum % MODULUS:064x}"}
    result["equalities"] = sorted(matches, key=lambda row: encode(row))
    result["scoped_equality_rows"] = sum(row["original_k"] is not None for row in matches)
    result["original_types"] = sorted({(row["original_k"], row["original_t"])
                                      for row in matches if row["original_k"] is not None})
    result["status"] = "PASS: complete fixed domain; every equality has an actual isometry witness"
    # Convert tuples before comparing against JSON decoded from the release.
    return json.loads(encode(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h", type=int, choices=range(2, 7))
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = classify(args.h)
    output = encode(result)+"\n"
    assert len(output.encode()) <= OUTPUT_CAP
    path = SUPPLEMENT / "data" / f"occurrences-{args.h}.json"
    if args.write:
        path.write_text(output)
    else:
        assert json.loads(path.read_text()) == result, "Reference check mismatch"
    print(encode({key: value for key, value in result.items() if key != "equalities"}))
    print("reference_sha256="+hashlib.sha256(output.encode()).hexdigest())


if __name__ == "__main__":
    main()
