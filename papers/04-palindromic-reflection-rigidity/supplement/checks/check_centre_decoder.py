#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact controls for the decoder at every palindromic centre and its corollary.

Fixed scopes, chosen to run within the limits of ``run_all.py``:

1. for a general symmetric matrix G = ((g, h), (h, k)) with symbolic entries, and with
   M = e^T G e, rho = (G e)_1, K = 2 e e^T G - M I: the closed form of K and of its transpose,
   K A B = B A K^T, K B A = A B K^T, B^{-1} K A = M ((0, 1), (1, -1 - tau)) with tau = 2 rho / M,
   K^2 = M^2 I, K e = M e, G K = K^T G, K^T G e = M G e, K^T (1, -2) = -M (1, -2), and the two
   slope formulas; no determinant condition is used;
2. every rational constant of the proof: the first-letter intervals, the window of tau, the
   bounds of the first two letters, the two quotients 131/181 and 71/55, the bounds
   11098/7421 < 21/13 and 158/55 > 169/70, the bound 29/168, and the images of the negative
   half-line under A^{-1} and B^{-1};
3. the theorem, by exhaustion: for every palindromic centre of length at most four and all words
   x, y of length at most nine, the lines [R_p mu(x) e] and [mu(y) e] coincide exactly for the
   block words x in {ab p, ba p}* and their exchanges, with exact equality of columns;
4. the corollary, on all palindromes of label at most 10^9: every pair of twins whose root
   midpoint is the ratio of a palindrome is an exchange around that palindrome.

All arithmetic is exact.  The last line prints a digest of the scope summary.
"""
import hashlib
import itertools
import sys
from fractions import Fraction as F

import sympy as sp

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed checker accepts no arguments.")

summary = []

# 1: symbolic identities for a general symmetric centre matrix
g, h, k, t = sp.symbols("g h k t")
A = sp.Matrix([[2, 1], [1, 1]])
B = sp.Matrix([[5, 2], [2, 1]])
e = sp.Matrix([2, 1])
G = sp.Matrix([[g, h], [h, k]])
M = (e.T * G * e)[0]
rho = (G * e)[0]
tau = 2 * rho / M
K = 2 * e * e.T * G - M * sp.eye(2)
f = sp.Matrix([1, -2])
identities = {
    "closed form of K": K - sp.Matrix([[4 * rho - M, 4 * M - 8 * rho], [2 * rho, M - 4 * rho]]),
    "K A B = B A K^T": K * A * B - B * A * K.T,
    "K B A = A B K^T": K * B * A - A * B * K.T,
    "B^-1 K A": B.inv() * K * A - M * sp.Matrix([[0, 1], [1, -1 - tau]]),
    "K^2 = M^2 I": K * K - M ** 2 * sp.eye(2),
    "K e = M e": K * e - M * e,
    "G K = K^T G": G * K - K.T * G,
    "K^T G e = M G e": K.T * G * e - M * G * e,
    "K^T (1,-2) = -M (1,-2)": K.T * f + M * f,
    "slope action of K": K * sp.Matrix([t, 1]) - M * sp.Matrix([2 * (1 + tau * (t - 2)) + 2 - t, 1 + tau * (t - 2)]),
    "slope action of K^T": K.T * sp.Matrix([t, 1]) - M * sp.Matrix([(2 * tau - 1) * t + tau, (4 - 4 * tau) * t + 1 - 2 * tau]),
}
for name, expr in identities.items():
    assert sp.simplify(expr) == sp.zeros(*expr.shape), name
summary.append(f"identities {len(identities)}")


# 2: rational constants
def act(m, x):
    return F(m[0][0] * x + m[0][1], m[1][0] * x + m[1][1])


A_ = ((2, 1), (1, 1))
B_ = ((5, 2), (2, 1))
J = (F(8, 5), F(29, 12))
assert (act(A_, J[0]), act(A_, J[1])) == (F(21, 13), F(70, 41))
assert (act(B_, J[0]), act(B_, J[1])) == (F(50, 21), F(169, 70))
assert J[0] < F(21, 13) and F(70, 41) < 2 < F(50, 21) and F(169, 70) < J[1]
window = tuple(2 * x / (2 * x + 1) for x in J)
assert window == (F(16, 21), F(29, 35))
assert 1 - F(2, 5) * window[1] == F(117, 175) > 0
assert 1 / (1 - window[0]) == F(21, 5) > F(169, 70)
assert F(70, 41) < 1 + window[0] == F(37, 21)
assert (37 - 42 * window[1]) / 5 == F(11, 25) > 0
qa = (F(21, 13) + F(2, 5)) / (4 - F(21, 13) + F(2, 5))
qb = (F(50, 21) + 1) / (4 - F(50, 21) + 1)
assert (qa, qb) == (F(131, 181), F(71, 55))
upper_a = F(70, 41) - F(12, 41) * qa
lower_b = F(50, 21) + F(8, 21) * qb
assert upper_a == F(11098, 7421) < F(21, 13) and lower_b == F(158, 55) > F(169, 70)
assert (F(169, 70) - 2) / (2 + F(2, 5)) == F(29, 168) < 1
for x in (F(-1, 10 ** 6), F(-1), F(-10 ** 6)):
    ai = F(x - 1, 2 - x)
    bi = F(x - 2, 5 - 2 * x)
    assert -1 < ai < F(-1, 2) and F(-1, 2) < bi < F(-2, 5)
assert sp.simplify(sp.Rational(1) * ((t - 1) / (2 - t) + 1) - 1 / (2 - t)) == 0
assert sp.simplify((t - 2) / (5 - 2 * t) + sp.Rational(2, 5) - t / (5 * (5 - 2 * t))) == 0
summary.append("constants 18")


# 3: the theorem by exhaustion
def mul(x, y):
    return tuple(tuple(sum(x[i][l] * y[l][j] for l in (0, 1)) for j in (0, 1)) for i in (0, 1))


def mu(word):
    out = ((1, 0), (0, 1))
    for c in word:
        out = mul(out, A_ if c == "a" else B_)
    return out


def column(word):
    m = mu(word)
    return (2 * m[0][0] + m[0][1], 2 * m[1][0] + m[1][1])


def k_matrix(p):
    r, s2 = column(p)
    m = 2 * r + s2
    return ((4 * r - m, 4 * m - 8 * r), (2 * r, m - 4 * r)), m


words = ["".join(w) for n in range(10) for w in itertools.product("ab", repeat=n)]
slope_of = {}
for w in words:
    c = column(w)
    slope_of.setdefault(F(c[0], c[1]), []).append(w)
centres = [p for p in ("".join(w) for n in range(5) for w in itertools.product("ab", repeat=n)) if p == p[::-1]]
EXCH = {"ab": "ba", "ba": "ab"}
pairs = 0
for p in centres:
    kp, mp = k_matrix(p)
    expected = set()
    for n in range(0, 10):
        for blocks in itertools.product(("ab", "ba"), repeat=n):
            x = "".join(b + p for b in blocks)
            if len(x) > 9:
                break
            expected.add((x, "".join(EXCH[b] + p for b in blocks)))
    found = set()
    for x in words:
        c = column(x)
        v = (kp[0][0] * c[0] + kp[0][1] * c[1], kp[1][0] * c[0] + kp[1][1] * c[1])
        if v[1] == 0:
            continue
        for y in slope_of.get(F(v[0], v[1]), []):
            found.add((x, y))
            cy = column(y)
            assert v == (mp * cy[0], mp * cy[1]), (p, x, y)
    assert found == {(x, y) for x, y in expected if len(y) <= 9}, p
    pairs += len(found)
summary.append(f"decoder {len(centres)} {len(words)} {pairs}")

# 4: the corollary on all palindromes of label at most 10^9
pal = {}
stack = [("", 2, 1)]
while stack:
    w, x, y = stack.pop()
    for z in (w[::-1] + w, w[::-1] + "a" + w, w[::-1] + "b" + w):
        c = column(z)
        m = 2 * c[0] + c[1]
        if m <= 10 ** 9:
            pal[z] = (m, c[0])
    if x * x + y * y <= 10 ** 9:
        stack.append(("a" + w, 2 * x + y, x + y))
        stack.append(("b" + w, 5 * x + 2 * y, 2 * x + y))
by_ratio = {F(r, m): z for z, (m, r) in pal.items()}
by_label = {}
for z, (m, r) in pal.items():
    by_label.setdefault(m, []).append(z)
axis = 0
for m, group in by_label.items():
    for z1, z2 in itertools.combinations(sorted(group), 2):
        mid = F(pal[z1][1] + pal[z2][1], 2 * m)
        p = by_ratio.get(mid)
        if p is None:
            continue
        b1 = z1[len(p):] if z1.startswith(p) else None
        b2 = z2[len(p):] if z2.startswith(p) else None
        assert b1 is not None and b2 is not None, (z1, z2, p)
        size = len(p) + 2
        blocks = [b1[i:i + size] for i in range(0, len(b1), size)]
        assert all(b[:2] in EXCH and b[2:] == p for b in blocks), (z1, p)
        assert "".join(EXCH[b[:2]] + p for b in blocks) == b2, (z1, z2, p)
        axis += 1
summary.append(f"axis {len(pal)} {axis}")

text = "; ".join(summary)
print("PASS: decoder at every palindromic centre, its identities and constants, and the classification at every centre (" + text + ")")
print("SHA256 " + hashlib.sha256(text.encode()).hexdigest())
