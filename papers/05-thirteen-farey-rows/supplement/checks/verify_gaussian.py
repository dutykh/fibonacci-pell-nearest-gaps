#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Verify the entire finite Gaussian package using only exact integers.

Each label is computed twice and compared: by Stern--Brocot matrix mutation and from the
floor increments of the lower Christoffel word. The occurrence pair is recomputed from the
central Gram factor H of the occurrence (mu(z) = H^T H) and compared with the record.
Independent fibre construction: all signs before quotienting, with twofold multiplicity.
No factorization, probable-prime predicate or external module is used.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import gcd, prod
from pathlib import Path
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

DATA = Path(__file__).resolve().parents[1] / "data"
WINDOWS = {8: 28, 9: 34, 10: 38, 11: 75, 12: 100, 13: 145}


def read(name):
    return json.loads((DATA/name).read_text())


def ordered_factors(raw):
    factors = tuple((int(p), int(e)) for p, e in raw)
    assert factors == tuple(sorted(set(factors)))
    assert len({p for p, _ in factors}) == len(factors)
    assert all(p >= 2 and e >= 1 for p, e in factors)
    return factors


prime_data = {int(p): record for p, record in read("lucas_primes.json").items()}
assert 2 in prime_data
proved = set()
for p in sorted(prime_data):
    check = prime_data[p]
    factors = ordered_factors(check["factors"])
    witness = check["witness"]
    if p == 2:
        assert witness == 1 and factors == ()
    else:
        assert p > 2 and p % 2 == 1 and 1 < witness < p
        assert prod(r**e for r, e in factors) == p-1
        assert all(r in proved and r < p for r, _ in factors)
        assert pow(witness, p-1, p) == 1
        assert all(gcd(pow(witness, (p-1)//r, p)-1, p) == 1 for r, _ in factors)
    proved.add(p)

splits = {int(p): tuple(pair) for p, pair in read("gaussian_splits.json").items()}
for p, (a, b) in splits.items():
    assert p in proved and p % 4 == 1
    assert 0 < a < b and a*a+b*b == p and gcd(a, b) == 1


def matrix_product(x, y):
    return [[sum(x[i][k]*y[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]


def source_by_mediants(p, q):
    lp, lq, rp, rq = 0, 1, 1, 1
    left = [[2, 1], [1, 1]]
    right = [[5, 2], [2, 1]]
    for _ in range(q):
        a, b = lp+rp, lq+rq
        middle = matrix_product(left, right)
        difference = p*b-a*q
        if difference == 0:
            assert a == p and b == q
            assert middle[0][0]+middle[1][1] == 3*middle[0][1]
            assert middle[0][0]*middle[1][1]-middle[0][1]*middle[1][0] == 1
            return middle[0][1]
        if difference < 0:
            rp, rq, right = a, b, middle
        else:
            lp, lq, left = a, b, middle
    raise AssertionError("Stern--Brocot construction did not terminate")


def source_by_floor_word(p, q):
    word = ["b" if (i*p)//q - ((i-1)*p)//q else "a" for i in range(1, q+1)]
    value = [[1, 0], [0, 1]]
    for letter in word:
        value = matrix_product(value, [[2, 1], [1, 1]] if letter == "a" else [[5, 2], [2, 1]])
    assert value[0][0]+value[1][1] == 3*value[0][1]
    return value[0][1], word


def occurrence_by_gram(p, q):
    _, word = source_by_floor_word(p, q)
    z = word[1:-1]
    assert z == z[::-1]
    n = len(z)
    h = [[1, 0], [0, 1]]
    for letter in z[(n+1)//2:]:
        h = matrix_product(h, [[2, 1], [1, 1]] if letter == "a" else [[5, 2], [2, 1]])
    if n % 2:
        h = matrix_product([[1, 1], [1, 0]] if z[n//2] == "a" else [[2, 1], [1, 0]], h)
    return tuple(sorted((2*h[0][0]+h[0][1], 2*h[1][0]+h[1][1])))


def gaussian_product(x, y):
    # Deliberately use a real 2x2 multiplication instead of the exporter formula.
    value = matrix_product([[x[0], -x[1]], [x[1], x[0]]], [[y[0], 0], [y[1], 0]])
    return value[0][0], value[1][0]


def power(z, e):
    result = (1, 0)
    while e:
        if e % 2:
            result = gaussian_product(result, z)
        z = gaussian_product(z, z)
        e //= 2
    return result


def euclid(pair):
    a, b = max(pair), min(pair)
    quotients = []
    while b:
        quotients.append(a//b)
        a, b = b, a % b
    assert a == 1 and quotients[-1] >= 2
    return quotients, sum(quotients)-1


expected = {(p, q) for p, bound in WINDOWS.items() for q in range(p+1, bound+1) if gcd(p, q) == 1}
paths = sorted((DATA/"gaussian").glob("*.json"))
assert len(paths) == len(expected) == 248
seen = set()
summary = {p: {"sources": 0, "primitive_pairs": 0, "largest_fibre": 0} for p in WINDOWS}
full_records = []
for path in paths:
    row = json.loads(path.read_text())
    p, q, m = row["P"], row["Q"], row["M"]
    assert (p, q) in expected and (p, q) not in seen
    seen.add((p, q))
    assert path.name == f"{p:02d}-{q:03d}.json"
    assert source_by_mediants(p, q) == m
    assert source_by_floor_word(p, q)[0] == m
    assert occurrence_by_gram(p, q) == tuple(row["occurrence"])
    factors = ordered_factors(row["factors"])
    assert all(r in proved for r, _ in factors)
    assert prod(r**e for r, e in factors) == m
    exponent_two = dict(factors).get(2, 0)
    assert exponent_two in (0, 1)
    odd = [(r, e) for r, e in factors if r != 2]
    assert odd and all(r in splits for r, _ in odd)
    blocks = [power(splits[r], e) for r, e in odd]
    multiplicities = Counter()
    for signs in product((-1, 1), repeat=len(blocks)):
        value = (1, 1) if exponent_two else (1, 0)
        for z, sign in zip(blocks, signs):
            value = gaussian_product(value, (z[0], sign*z[1]))
        pair = tuple(sorted(abs(c) for c in value))
        assert 0 < pair[0] < pair[1] and gcd(*pair) == 1
        assert pair[0]**2+pair[1]**2 == m
        multiplicities[pair] += 1
    assert len(multiplicities) == 2**(len(odd)-1)
    assert set(multiplicities.values()) == {2}
    reconstructed = []
    for pair in multiplicities:
        quotients, cost = euclid(pair)
        reconstructed.append({"pair": list(pair), "quotients": quotients, "cost": cost})
    reconstructed.sort(key=lambda r: (r["cost"], r["pair"]))
    assert row["representations"] == reconstructed
    assert row["occurrence_cost"] == p+q-2
    assert row["occurrence"] == reconstructed[0]["pair"]
    assert reconstructed[0]["cost"] == p+q-2
    assert len(reconstructed) == 1 or reconstructed[1]["cost"] > p+q-2
    summary[p]["sources"] += 1
    summary[p]["primitive_pairs"] += len(reconstructed)
    summary[p]["largest_fibre"] = max(summary[p]["largest_fibre"], len(reconstructed))
    full_records.append(row)
assert seen == expected
assert [summary[p]["primitive_pairs"] for p in (8, 9, 10, 11, 12, 13)] == \
    [22, 56, 33, 324, 211, 992]
assert [summary[p]["sources"] for p in (8, 9, 10, 11, 12, 13)] == [10, 17, 11, 59, 29, 122]
payload = {"sources": full_records, "prime_checks": prime_data, "gaussian_splits": splits}
digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
print(json.dumps({"status": "PASS", "sources": len(seen), "prime_checks": len(proved),
                  "gaussian_splits": len(splits), "by_row": summary,
                  "reference_sha256": digest}, sort_keys=True))
