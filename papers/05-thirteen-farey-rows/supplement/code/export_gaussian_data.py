#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Prepare one prescribed Gaussian fibre and its exact prime checks.

Factor discovery is optional preparation, never used by the public verifier.
Run P Q for a member of the fixed window, under the documented resource cap.
"""
import argparse
from itertools import product
import json
from math import gcd, isqrt, prod
from pathlib import Path
import sys
from sympy import factorint
from source_words import source, euclidean_data

if not __debug__:
    raise SystemExit("Assertions are required; do not use -O.")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WINDOWS = {8: 28, 9: 34, 10: 38, 11: 75, 12: 100, 13: 145}


def read_json(name, default):
    p = DATA / name
    return json.loads(p.read_text()) if p.exists() else default


def atomic_json(name, value):
    p = DATA / name
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_suffix(p.suffix + ".tmp")
    temp.write_text(json.dumps(value, sort_keys=True, indent=2)+"\n")
    temp.replace(p)


hints = read_json("factor_hints.json", {})
checks = read_json("lucas_primes.json", {})
splits = read_json("gaussian_splits.json", {})


def factors_of(n):
    fixed = hints["minus_one_factors"].get(str(n+1))
    if fixed is None:
        print("factor", n, flush=True)
        fixed = sorted((int(p), int(e)) for p, e in factorint(n).items())
    assert prod(p**e for p, e in fixed) == n
    assert all(p >= 2 and e >= 1 for p, e in fixed)
    return fixed


def prove_prime(p):
    if str(p) in checks:
        return
    if p == 2:
        checks["2"] = {"witness": 1, "factors": []}
        return
    assert p > 2 and p % 2
    factors = factors_of(p-1)
    for q, _ in factors:
        assert q < p
        prove_prime(q)
    for a in range(2, 100):
        if pow(a, p-1, p) == 1 and all(gcd(pow(a, (p-1)//q, p)-1, p) == 1 for q, _ in factors):
            checks[str(p)] = {"witness": a, "factors": factors}
            return
    raise AssertionError(("No full Lucas witness below 100", p))


def split_prime(p):
    if str(p) in splits:
        return tuple(splits[str(p)])
    prove_prime(p)
    assert p % 4 == 1
    a = 2
    while pow(a, (p-1)//2, p) != p-1:
        a += 1
    root = pow(a, (p-1)//4, p)
    assert root*root % p == p-1
    x, y = p, root
    while y*y > p:
        x, y = y, x % y
    z = isqrt(p-y*y)
    assert y*y+z*z == p
    pair = tuple(sorted((y, z)))
    assert pair[0] > 0 and gcd(*pair) == 1
    splits[str(p)] = pair
    return pair


def gaussian_product(x, y):
    return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def gaussian_power(z, exponent):
    result = (1, 0)
    for _ in range(exponent):
        result = gaussian_product(result, z)
    return result


def export(p, q):
    m, occurrence = source(p, q)
    factors = hints["source_factors"].get(f"{p}/{q}")
    if factors is None:
        factors = factors_of(m)
    assert prod(r**e for r, e in factors) == m
    for r, _ in factors:
        prove_prime(r)
    two = dict(factors).get(2, 0)
    assert two in (0, 1)
    odd = [(r, e) for r, e in factors if r != 2]
    assert odd and all(r % 4 == 1 for r, _ in odd)
    blocks = [gaussian_power(split_prime(r), e) for r, e in odd]
    pairs = set()
    for signs in product((1, -1), repeat=len(blocks)-1):
        z = (1, 1) if two else (1, 0)
        for block, sign in zip(blocks, (1,)+signs):
            z = gaussian_product(z, (block[0], sign*block[1]))
        pair = tuple(sorted(map(abs, z)))
        assert 0 < pair[0] < pair[1] and gcd(*pair) == 1
        assert pair[0]*pair[0]+pair[1]*pair[1] == m
        pairs.add(pair)
    assert len(pairs) == 2**(len(odd)-1)
    rows = []
    for pair in pairs:
        quotients, cost = euclidean_data(pair)
        rows.append({"pair": pair, "quotients": quotients, "cost": cost})
    rows.sort(key=lambda row: (row["cost"], row["pair"]))
    assert occurrence in pairs
    assert rows[0]["cost"] == p+q-2
    assert rows[0]["pair"] == occurrence
    assert len(rows) == 1 or rows[1]["cost"] > p+q-2
    payload = {"P": p, "Q": q, "M": m, "factors": factors,
               "occurrence": occurrence, "occurrence_cost": p+q-2,
               "representations": rows}
    atomic_json(f"gaussian/{p:02d}-{q:03d}.json", payload)
    atomic_json("lucas_primes.json", checks)
    atomic_json("gaussian_splits.json", splits)
    print(f"PASS P={p} Q={q} pairs={len(rows)} prime_nodes={len(checks)}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("P", type=int)
    parser.add_argument("Q", type=int)
    args = parser.parse_args()
    p, q = args.P, args.Q
    if not (p in WINDOWS and p < q <= WINDOWS[p] and gcd(p, q) == 1):
        parser.error("P,Q must belong to the fixed manuscript window")
    export(p, q)


if __name__ == "__main__":
    main()
