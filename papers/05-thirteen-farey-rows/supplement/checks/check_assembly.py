#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Exact six source junctions and constant signs used in Section 6.

The independent arithmetic here uses rational coefficients of sqrt(5),
rather than the tail engine's integer coefficients of phi. Source labels
are reconstructed both from floor words and Stern--Brocot products.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from source_words import source, fibonacci

def add(x, y):
    return x[0]+y[0], x[1]+y[1]

def neg(x):
    return -x[0], -x[1]

def mul(x, y):
    return x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0]

def power(x, n):
    out = (F(1), F(0))
    for _ in range(n):
        out = mul(out, x)
    return out

def sign(x):
    a, b = x
    if not a:
        return (b > 0)-(b < 0)
    if not b or (a > 0) == (b > 0):
        return (a > 0)-(a < 0)
    delta = a*a-5*b*b
    assert delta
    return ((a > 0)-(a < 0))*((delta > 0)-(delta < 0))

def rational(n, d=1):
    return F(n, d), F(0)

def positive_difference(x, y):
    assert sign(add(x, neg(y))) > 0

phi = (F(1, 2), F(1, 2))
phi_inverse = (F(-1, 2), F(1, 2))
row_ratio = (F(3, 2), F(-3, 10))
constants = {p: mul(power(phi, 2), power(row_ratio, p-1))
             for p in range(2, 14)}
ceilings = {8: (707, 1000), 9: (2929, 5000), 10: (4857, 10000),
            11: (20111, 50000), 12: (333511, 10**6),
            13: (2765404267, 10**10)}
starts = {8: 29, 9: 35, 10: 39, 11: 76, 12: 101, 13: 146}
costs = {8: 34, 9: 41, 10: 46, 11: 84, 12: 110, 13: 156}

checks = {
    "C11 > 2/5": (constants[11], rational(2, 5)),
    "C12 > 1/3": (constants[12], rational(1, 3)),
    "C13 > 71/260": (constants[13], rational(71, 260)),
    "C13 > 3/20": (constants[13], rational(3, 20)),
    "5 phi^-4 > 29/40": (mul(rational(5), power(phi_inverse, 4)), rational(29, 40)),
    "29/40 > U8": (rational(29, 40), rational(*ceilings[8])),
    "phi^-2 > 3/8": (power(phi_inverse, 2), rational(3, 8)),
    "3/8 > U12": (rational(3, 8), rational(*ceilings[12])),
    "5 phi^-6 > 5/18": (mul(rational(5), power(phi_inverse, 6)), rational(5, 18)),
    "5/18 > U13": (rational(5, 18), rational(*ceilings[13])),
}
for name, values in checks.items():
    assert sign(add(values[0], neg(values[1]))) > 0, name
for p in range(2, 13):
    positive_difference(constants[p], constants[p+1])
for p in range(8, 13):
    positive_difference(rational(*ceilings[p]), rational(*ceilings[p+1]))
for p in starts:
    positive_difference(rational(*ceilings[p]), constants[p])
    positive_difference(rational(1), rational(*ceilings[p]))

def matrix_product(x, y):
    return tuple(sum(x[2*i+k]*y[2*k+j] for k in (0, 1))
                 for i in (0, 1) for j in (0, 1))

def mediant_source(p, q):
    left, right = (0, 1), (1, 1)
    a, b = (2, 1, 1, 1), (5, 2, 2, 1)
    for _ in range(q):
        middle = (left[0]+right[0], left[1]+right[1])
        product = matrix_product(a, b)
        delta = p*middle[1]-q*middle[0]
        if delta == 0:
            assert middle == (p, q)
            assert product[0]+product[3] == 3*product[1]
            return product[1]
        if delta < 0:
            right, b = middle, product
        else:
            left, a = middle, product
    raise AssertionError("Mediant construction did not terminate")

records = []
for p, q in starts.items():
    assert gcd(p, q) == 1 and p+q-3 == costs[p]
    m, pair = source(p, q)
    assert m == mediant_source(p, q)
    index = 2*(p+q-2)+1
    f = fibonacci(index)
    numerator, denominator = ceilings[p]
    margin = numerator*f-denominator*m
    assert margin > 0
    records.append({"P": p, "Q": q, "candidate_cost": costs[p],
                    "M": m, "occurrence_pair": list(pair),
                    "fibonacci_index": index, "fibonacci_value": f,
                    "ceiling": [numerator, denominator], "positive_margin": margin})
expected = json.loads((ROOT / "data/source-junctions.json").read_text())
assert expected == records
digest = sha256(json.dumps({"sources": records, "sign_checks": sorted(checks)},
                          sort_keys=True, separators=(",", ":")).encode()).hexdigest()
print(json.dumps({"status": "PASS", "source_junctions": len(records),
                  "named_constant_signs": len(checks), "reference_sha256": digest},
                 sort_keys=True))
