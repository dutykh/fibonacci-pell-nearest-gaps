#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Independent left-morphism reconstruction of the five table examples."""
import hashlib
import json
import math
import sys
from pathlib import Path
from fractions import Fraction as F

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed checker accepts no arguments.")


def multiply(x, y):
    return [[sum(x[i][k]*y[k][j] for k in (0, 1)) for j in (0, 1)] for i in (0, 1)]


def transpose(x):
    return [list(row) for row in zip(*x)]


def matrix(word):
    out = [[1, 0], [0, 1]]
    for c in word:
        out = multiply(out, {"a": [[2, 1], [1, 1]], "b": [[5, 2], [2, 1]]}[c])
    return out


def left_image(directive, word):
    for c in directive[::-1]:
        rules = {"a": "a", "b": "ab"} if c == "a" else {"a": "ba", "b": "b"}
        word = "".join(rules[x] for x in word)+c
    return word


path = Path(__file__).resolve().parents[1] / "data" / "midpoint_examples.json"
raw = path.read_bytes()
rows = json.loads(raw)
assert len(rows) == 5
expected = [("", 5, 2, 1130, 437, 467), ("a", 13, 5, 19786, 7571, 7649),
            ("b", 29, 12, 219530, 90753, 90927),
            ("ab", 194, 75, 65712650, 25403793, 25404957),
            ("ba", 433, 179, 730645066, 302043659, 302046257)]
p = [[2, 1], [1, 0]]
for row, values in zip(rows, expected):
    d, m0, r0, m, lo, hi = values
    assert row["directive"] == d
    ancestor = left_image(d, "")
    assert row["ancestor"] == ancestor
    g0 = multiply(multiply(p, matrix(ancestor)), p)
    assert g0[0] == [m0, r0]
    assert row["ancestor_gram"] == {"M": m0, "rho": r0, "H": g0[1][1]}
    grams = []
    for saved, preimage, rho in zip(row["images"], ("baab", "abba"), (hi, lo)):
        word = left_image(d, preimage)
        assert word == word[::-1] == saved["word"]
        g = multiply(multiply(p, matrix(word)), p)
        assert g[0] == [m, rho] and g[0][1] == g[1][0]
        assert g[0][0]*g[1][1]-g[0][1]**2 == 1
        assert saved["gram"] == {"M": m, "rho": rho, "H": g[1][1]}
        counts = [word.count("a")+1, word.count("b")+1]
        assert counts == saved["completed_counts"]
        assert counts == [3*(ancestor.count(c)+1) for c in "ab"]
        assert math.gcd(*counts) == 3
        grams.append(g)
    assert (lo+hi)*m0 == 2*m*r0 and hi-lo == 6*m0
    j = [[1, F(2*r0, m0)], [0, -1]]
    assert multiply(multiply(transpose(j), grams[0]), j) == grams[1]
print("PASS: five independently reconstructed Gram/count/reflection examples")
print("SHA256 " + hashlib.sha256(raw).hexdigest())
