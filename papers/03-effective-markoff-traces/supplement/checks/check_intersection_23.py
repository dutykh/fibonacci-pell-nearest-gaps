#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Intersection twenty-three: exact trace identities and polynomial sign checks.

Both compared traces are rebuilt from the universal character algebra of the
basis (X, Y), with no numerical root search and no parameter census. The
elimination, the fifteen positive coefficient lists, the two power-trace
identities, the four grouped sign estimates and the two rational margins are
all confirmed. Initial external limits are 128 MiB / 30 seconds. Requires
Python 3.10+ and SymPy.
Run: python3 -B checks/check_intersection_23.py
Add --write to regenerate the accompanying complete coefficient data.
"""
import hashlib
import json
from pathlib import Path
import sys

if not __debug__ or sys.argv[1:] not in ([], ["--write"]):
    raise SystemExit("Run without -O; the only optional argument is --write.")

SUPPLEMENT = Path(__file__).resolve().parents[1]
DATA = SUPPLEMENT / "data" / "intersection-23-polynomials.json"
WRITE = sys.argv[1:] == ["--write"]

import sympy as S

p, r, s, u, z = S.symbols("p r s u z")
t = 3*r*s - p                                   # so that tr(XY) = 3t, tr(XY^-1) = 3p


def trace_word(word):
    """Right multiplication in (I,X,Y,XY); lowercase means the inverse letter."""
    tx, ty, txy = 3*r, 3*s, 3*t
    vector = (S.Integer(1), S.Integer(0), S.Integer(0), S.Integer(0))
    for letter in word:
        a0, b0, c0, d0 = vector
        if letter.upper() == "X":
            product = (-b0 + (txy - tx*ty)*c0 - ty*d0,
                       a0 + tx*b0 + ty*c0 + txy*d0, tx*c0 + d0, -c0)
            trace = tx
        else:
            assert letter.upper() == "Y"
            product = (-c0, -d0, a0 + ty*c0, b0 + ty*d0)
            trace = ty
        if letter.islower():
            product = tuple(trace*v - w for v, w in zip(vector, product))
        vector = tuple(S.expand(v) for v in product)
    a0, b0, c0, d0 = vector
    return S.expand(2*a0 + tx*b0 + ty*c0 + txy*d0)


# The mechanical gap pattern of homology (7,23), and its marking orbit.
assert [23*(j+1)//7 - 23*j//7 for j in range(7)] == [3, 3, 3, 4, 3, 3, 4]
assert {7, -7 % 23, pow(7, -1, 23), -pow(7, -1, 23) % 23} == {7, 10, 13, 16}

# B = X^-1 Y and A = X B^-3, so the original class is the word X (Y^-1 X)^3.
assert S.expand(trace_word("xY") - 3*p) == 0
original = (27*p**3 - 6*p)*r - (9*p*p - 1)*s
assert S.expand(trace_word("XyXyXyX")/3 - original) == 0

# The competing class X^3 Y X^2 Y = X W^2 with W = X^2 Y.
a = (9*r*r - 1)*s - 3*p*r
b = (27*r**3 - 6*r)*s - (9*r*r - 1)*p
assert S.expand(trace_word("XXY")/3 - a) == 0
assert S.expand(trace_word("XXXY")/3 - b) == 0
target = 3*a*b - r
assert S.expand(trace_word("XXXYXXY")/3 - target) == 0

markov = p*p + r*r + s*s - 3*p*r*s
difference = S.expand(original - target)
quotient, remainder = S.div(difference, markov, p)
linear = S.Poly(remainder, p)
assert linear.degree() == 1
L, N = linear.nth(1), -linear.nth(0)
assert S.expand(difference - quotient*markov - remainder) == 0
Q = S.expand(N*N - 3*r*s*N*L + (r*r + s*s)*L*L)
assert S.expand(S.resultant(difference, markov, p) - Q) == 0

# Fifteen coefficient lists: every shifted coefficient is positive on the half-line.
rows = {}
for name, expression, parameter in (("s_ge_r", Q.subs(s, r + u), r),
                                    ("r_ge_s", Q.subs(r, s + u), s)):
    polynomial = S.Poly(expression, u)
    entries = []
    for degree in range(polynomial.degree(), -1, -1):
        coefficient = S.Poly(polynomial.coeff_monomial(u**degree), parameter)
        exponents = [e[0] for e, _ in coefficient.terms()]
        parity = exponents[0] % 2
        assert all(e % 2 == parity for e in exponents)
        shifted = S.Poly(sum(c*(z + 1)**((e[0] - parity)//2)
                             for e, c in coefficient.terms()), z)
        values = [int(v) for v in shifted.all_coeffs()]
        top = degree == polynomial.degree()
        # The leading coefficient is positive; every lower one is negative, so
        # each b_j or c_j is a polynomial in z with positive coefficients.
        assert all((v > 0) if top else (v < 0) for v in values)
        entries.append([degree, parity, values if top else [-v for v in values]])
    # Exactly one sign change, hence exactly one positive root.
    assert polynomial.coeff_monomial(u**polynomial.degree()) > 0
    rows[name] = entries
assert [entry[0] for entry in rows["s_ge_r"]] == list(range(6, -1, -1))
assert [entry[0] for entry in rows["r_ge_s"]] == list(range(10, -1, -1))
# The two constant terms are the same polynomial, both equal to Q(a,a), so the
# eighteen entries above carry fifteen distinct nonleading coefficient lists.
assert rows["s_ge_r"][-1] == rows["r_ge_s"][-1]
assert rows["s_ge_r"][0][2] == [81] and rows["r_ge_s"][0][2] == [6561]

# The two power traces, and the factorisations that place them above the diagonal.
R9 = 6561*r**9 - 6561*r**7 + 2187*r**5 - 270*r**3 + 9*r
R5 = 81*s**5 - 45*s**3 + 5*s
power = [S.Integer(2), 3*z]
for j in range(2, 10):
    power.append(S.expand(3*z*power[-1] - power[-2]))
assert S.expand(power[9]/3 - R9.subs(r, z)) == 0
assert S.expand(power[5]/3 - R5.subs(s, z)) == 0
assert S.expand(R9 - 1926*r**9 - r*(r*r - 1)*(4635*r**6 - 1926*r**4 + 261*r*r - 9)) == 0
assert S.expand(R5 - 41*s**5 - 5*s*(s*s - 1)*(8*s*s - 1)) == 0

# The two groupings.
grouped_s = (81*s**5*(s - R9) + (324*r*r + 27)*s**4
             + (295245*r**9 - 164025*r**7 + 19683*r**5 + 486*r**3 - 63*r)*s**3
             + (972*r**4 + 1)*s*s
             + (-32805*r**9 + 3645*r**7 + 729*r**5 - 54*r**3 + 2*r)*s
             + 6561*r**10 - 729*r**8 + 567*r**6 + 18*r**4 + r*r)
grouped_r = (6561*r**9*(r - R5) - 729*r**8
             + (531441*s**5 - 164025*s**3 + 3645*s)*r**7 + 567*r**6
             + (-177147*s**5 + 19683*s**3 + 729*s)*r**5
             + (972*s*s + 18)*r**4 + (21870*s**5 + 486*s**3 - 54*s)*r**3
             + (324*s**4 + 1)*r*r + (-729*s**5 - 63*s**3 + 2*s)*r
             + 81*s**6 + 27*s**4 + s*s)
assert S.expand(Q - grouped_s) == 0
assert S.expand(Q - grouped_r) == 0

# The four grouped sign estimates.
assert 295245 - 164025 - 63 == 131157
assert 32805 + 54 == 32859
assert 131157 > 32859
assert S.expand(531441*s**5 - 164025*s**3 + 3645*s - 729*R5
                - (472392*s**5 - 131220*s**3)) == 0
assert 472392 - 131220 == 341172 > 177147
assert 21870 - 54 == 21816 > 729 + 63 == 792

# The two rational margins.
margin_s = (S.Rational(351, 1925) + S.Rational(315414, 1925**2)
            + S.Rational(973, 1925**3) + S.Rational(4376, 1925**4)
            + S.Rational(7147, 1925**5))
margin_r = (S.Rational(535086, 40**2) + S.Rational(567, 40**3)
            + S.Rational(20412, 40**4) + S.Rational(990, 40**5)
            + S.Rational(22356, 40**6) + S.Rational(325, 40**7)
            + S.Rational(2, 40**8) + S.Rational(109, 40**9))
assert margin_s < 1 and 81 - margin_s > 80
assert margin_r < 335 and 6561 - margin_r > 6226

# Both endpoint signs, by a shift to the half-line.
endpoints = {}
for name, expression, parameter, sign in (
        ("s_upper", Q.subs(s, R9), r, 1), ("s_lower", Q.subs(s, R9 - 1), r, -1),
        ("r_upper", Q.subs(r, R5), s, 1), ("r_lower", Q.subs(r, R5 - 1), s, -1)):
    polynomial = S.Poly(expression, parameter)
    shifted = S.Poly(polynomial.as_expr().subs(parameter, z + 1), z)
    assert all(sign*c > 0 for c in shifted.all_coeffs())
    endpoints[name] = {"degree": polynomial.degree(),
                       "at_one": int(polynomial.eval(1)),
                       "shifted_coefficient_sign": sign}

payload = {"status": "PASS",
           "scope": "Exact intersection-23 elimination and both ordered integer gaps",
           "shift_rows": rows,
           "endpoint_polynomials": endpoints,
           "s_lower_positive_remainder_bound": str(margin_s),
           "r_lower_positive_remainder_bound": str(margin_r)}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
if WRITE:
    DATA.write_text(encoded)
else:
    assert json.loads(DATA.read_text()) == payload
print(json.dumps({"status": "PASS",
                  "distinct_coefficient_lists":
                      sum(len(v) - 1 for v in rows.values()) - 1,
                  "endpoint_polynomials": sorted(endpoints)}, sort_keys=True))
print("reference_sha256=" + hashlib.sha256(encoded.encode()).hexdigest())
