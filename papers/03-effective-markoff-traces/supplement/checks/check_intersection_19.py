#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Intersection nineteen: exact trace identities and polynomial sign checks.

No source-code imports, parameter census or numerical root search. Initial
external limits are 128 MiB / 30 seconds. Requires Python 3.10+ and SymPy.
Run: python3 -B checks/check_intersection_19.py
Add --write to regenerate the accompanying complete coefficient data.
"""
import hashlib
import json
from pathlib import Path
import sys

if not __debug__ or sys.argv[1:] not in ([], ["--write"]):
    raise SystemExit("Run without -O; the only optional argument is --write.")

SUPPLEMENT = Path(__file__).resolve().parents[1]
DATA = SUPPLEMENT / "data" / "intersection-19-polynomials.json"
WRITE = sys.argv[1:] == ["--write"]

import sympy as S

r, s, t, p, u, z = S.symbols("r s t p u z")


def trace_word(word):
    """Right multiplication in (I,X,Y,XY); lowercase means inverse."""
    tx, ty, txy = 3*r, 3*s, 3*t
    vector = (S.Integer(1), S.Integer(0), S.Integer(0), S.Integer(0))
    for letter in word:
        a, b, c, d = vector
        if letter.upper() == "X":
            product = (-b+(txy-tx*ty)*c-ty*d,
                       a+tx*b+ty*c+txy*d, tx*c+d, -c)
            trace = tx
        else:
            assert letter.upper() == "Y"
            product = (-c, -d, a+ty*c, b+ty*d)
            trace = ty
        if letter.islower():
            product = tuple(trace*v-w for v, w in zip(vector, product))
        vector = tuple(S.expand(v) for v in product)
    a, b, c, d = vector
    return S.expand(2*a+tx*b+ty*c+txy*d)


def P(r, s, t):
    pp = 3*r*s-t
    return (9*pp*pp-1)*r-3*pp*s


M = trace_word("XyXyX") / 3
N = trace_word("XYYXYYY") / 3
assert S.expand(M-P(r, s, t)) == 0
assert S.expand(N-P(s, t, r)) == 0
markov = r*r+s*s+t*t-3*r*s*t
difference = S.expand(M-N)
# Eliminate t, whereas the larger-first source eliminates p=3rs-t.
Q = S.expand(S.resultant(markov, difference, t))
R = 729*s**7-567*s**5+126*s**3-7*s
AA = 32805*s**7-13851*s**5+1134*s**3+9*s
BB = 324*s**4+1
CC = 3645*s**7-567*s**5+27*s**3+2*s
DD = 6561*s**10-1458*s**8+243*s**6-9*s**4+s*s
assert S.expand(Q-(81*r**5*(r-R)+27*r**4+AA*r**3+BB*r*r-CC*r+DD)) == 0
# Independently reconstruct an identity that includes the linear-degeneracy case.
quotient, remainder = S.div(difference, markov, t)
linear = S.Poly(remainder, t)
assert linear.degree() == 1
ll, nn = linear.nth(1), linear.nth(0)
assert S.expand(Q-(nn*nn+3*r*s*ll*nn+(r*r+s*s)*ll*ll)) == 0
assert S.expand(difference-quotient*markov-remainder) == 0


def positive_shift(poly, strict=True):
    """Independent positivity coordinate s=1+z, not source s^2=1+v."""
    shifted = S.Poly(S.expand(poly.subs(s, z+1)), z)
    assert all(c >= 0 for c in shifted.all_coeffs())
    if strict:
        assert shifted.nth(0) > 0
    return [int(v) for v in reversed(shifted.all_coeffs())]


shifted_q = S.Poly(S.expand(Q.subs(r, s+u)), u)
assert shifted_q.nth(6) == 81
shift_signs = [positive_shift(-shifted_q.nth(j)) for j in range(1, 6)]
assert S.expand(shifted_q.nth(0)+729*s**6*(s*s-1)*(9*s*s-2)**2) == 0
assert shifted_q.nth(0).subs(s, 1) == 0
assert shifted_q.nth(1).subs(s, 1) < 0

# Direct endpoint margins give the same human strict-sign proof with different
# coefficient checks; there is no expansion of Q(R-1,s).
endpoint_polynomials = [CC, AA-CC, BB, DD, R-s-1,
                        81*(R-1)-AA, 2*(R-1)-BB, (R-1)**2-DD]
endpoint_signs = [positive_shift(poly) for poly in endpoint_polynomials]
assert R.subs(s, 1) == 281
assert -81*2+108+S.Rational(2, 2)+S.Rational(1, 4) == -S.Rational(211, 4)
assert S.expand(Q.subs(r, s)+729*s**6*(s*s-1)*(9*s*s-2)**2) == 0
assert S.factor((p*p+r*r+s*s-3*p*r*s).subs({r: 1, s: 1})) == (p-1)*(p-2)
diagonal = [[int(M.subs({r: 1, s: 1, t: 3-pp})),
             int(N.subs({r: 1, s: 1, t: 3-pp}))] for pp in (1, 2)]
assert diagonal == [[5, 194], [29, 29]]

# Every preliminary comparison is an exact identity or the displayed rational
# margin. Their sign hypotheses are proved in the article.
assert S.expand((3*s*t-r)-(3*r*s-t)-(3*s+1)*(t-r)) == 0
assert S.expand(t*(3*r*s-t)-r*(3*s*t-r)-r*r+t*t) == 0
assert S.expand((3*t-1)**2-4*t-(t-1)*(9*t-1)) == 0
assert S.Rational(42, 5)-1 > 0 and S.Rational(32, 27) > 1
F = lambda a, b, c: a*a+b*b+c*c-3*a*b*c
assert S.expand(F(r, s, 3*r*s-t)-F(r, s, t)) == 0


class Poly:
    """Sparse integer polynomial in exactly three formal variables."""
    def __init__(self, data):
        if isinstance(data, int):
            data = {(0, 0, 0): data} if data else {}
        self.data = {k: v for k, v in data.items() if v}

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = Poly(other)
        out = self.data.copy()
        for exponent, coefficient in other.data.items():
            out[exponent] = out.get(exponent, 0)+coefficient
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.data.items()})

    def __sub__(self, other):
        return self+(-other if isinstance(other, Poly) else -int(other))

    def __rsub__(self, other):
        return (-self)+other

    def __mul__(self, other):
        if not isinstance(other, Poly):
            other = Poly(other)
        out = {}
        for a, ca in self.data.items():
            for b, cb in other.data.items():
                exponent = (a[0]+b[0], a[1]+b[1], a[2]+b[2])
                out[exponent] = out.get(exponent, 0)+ca*cb
        return Poly(out)

    __rmul__ = __mul__


A = Poly({(1, 0, 0): 1})
B = Poly({(0, 1, 0): 1})
C = Poly({(0, 0, 1): 1})


def descended(a, w, v):
    tt = 3*w*v-a
    rr = 3*tt*w-v
    ss = 3*tt*rr-w
    return P(rr, ss, tt)-P(ss, tt, rr)


base = descended(A, B, C)
assert len(base.data) == 151 and max(map(sum, base.data)) == 19
del base
cones = [
    (A+1, A+B+1, A+B+C+2, -1, 1234, -7969111823,
     "3f73f989110affc77b587dea9e4fb5980edd996b6b2d291f8e15f50d432203c9"),
    (A+B+2, A+1, A+B+C+2, -1, 659, -1010269301,
     "f3f3e5ec3909b1970d86f7c20ea30dcfce61c04493937be880e9908697ca5cbb"),
    (A+1, A+B+C+2, A+B+1, 1, 1436, 207308780444,
     "d9c984a7a87757dc2e0cc67c57c50037b551e7be4e8a74de02c74cad211eaa26"),
    (A+B+2, A+B+C+2, A+1, 1, 1091, 47654507165,
     "c801845cdbfcb28ec931fa9a081246c5980a455453e97d7ab4cd9d1c3bc58072"),
]
cone_results = []
cone_terms = []
for aa, ww, vv, sign, count, constant, digest in cones:
    poly = descended(aa, ww, vv)
    assert len(poly.data) == count <= 1540
    assert max(map(sum, poly.data)) == 19
    assert all(sign*value > 0 for value in poly.data.values())
    assert poly.data[(0, 0, 0)] == constant
    terms = [[list(k), poly.data[k]] for k in sorted(poly.data, reverse=True)]
    encoding = json.dumps(terms, separators=(",", ":"))
    assert hashlib.sha256(encoding.encode()).hexdigest() == digest
    cone_results.append([count, constant, digest])
    cone_terms.append(terms)
    del poly, terms, encoding

assert P(5, 29, 2)-P(29, 2, 5) == 945951 and F(5, 29, 2) == 0

# The actual modular normalizer, plus exact based-word and homology dictionaries.
E = S.Matrix([[2, 1], [1, 1]])
FF = S.Matrix([[2, -1], [-1, 1]])
G = S.Matrix([[0, 1], [-1, 1]])
MG = S.Matrix([[0, -1], [1, -1]])
assert G*E*G.inv() == FF and G*FF*G.inv() == (E*FF).inv()
assert G**3 == -S.eye(2)
assert S.trace((E*FF.inv())**2*E) == S.trace((E*FF**2)**2*FF) == 87
original_a = (E*FF.inv())**2*E
original_b = E.inv()*FF
assert tuple(S.trace(matrix) for matrix in (original_a, original_b, original_a*original_b)) == (87, 6, 15)
assert MG*S.Matrix([3, -2]) == S.Matrix([2, 5])
assert S.det(S.Matrix([[3, 2], [-2, 5]])) == 19
assert [19*i//7-19*(i-1)//7 for i in range(1, 8)] == [2, 3, 3, 2, 3, 3, 3]
original_basis = S.Matrix([[3, -1], [-2, 1]])
assert original_basis.det() == 1
assert original_basis*S.Matrix([1, 2]) == S.Matrix([1, 0])
assert original_basis*S.Matrix([1, 3]) == S.Matrix([0, 1])
assert original_basis*S.Matrix([7, 19]) == S.Matrix([2, 5])

expected_orbits = [{1, 18}, {2, 9, 10, 17}, {3, 6, 13, 16},
                   {4, 5, 14, 15}, {7, 8, 11, 12}]
orbits = []
for h in (1, 2, 3, 4, 7):
    inverse = pow(h, -1, 19)
    orbit = {h, -h % 19, inverse, -inverse % 19}
    orbits.append(orbit)
assert orbits == expected_orbits
assert set.union(*orbits) == set(range(1, 19))
assert sum(map(len, orbits)) == 18
assert [(19//h, 19 % h) for h in (2, 3, 4)] == [(9, 1), (6, 1), (4, 3)]

payload = {"status": "PASS", "scope": "both I19 trace orderings and all18 unit residues",
           "larger_first_shifted_coefficients": shift_signs,
           "independent_endpoint_margins": endpoint_signs,
           "sparse_cones": cone_results, "diagonal": diagonal,
           "forced_endpoint_difference": 945951,
           "residue_orbits": [sorted(orbit) for orbit in orbits]}
full_data = {
    "description": "Exact coefficients for the intersection-nineteen proof",
    "variable_order": ["A", "B", "C"],
    "cone_substitutions": [
        ["A+1", "A+B+1", "A+B+C+2"],
        ["A+B+2", "A+1", "A+B+C+2"],
        ["A+1", "A+B+C+2", "A+B+1"],
        ["A+B+2", "A+B+C+2", "A+1"],
    ],
    "coefficient_convention": "Each entry is [[degree_A,degree_B,degree_C],integer_coefficient]",
    "all_cone_coefficients": cone_terms,
    "summary": payload,
}
if WRITE:
    DATA.write_text(json.dumps(full_data, sort_keys=True, separators=(",", ":"))+"\n")
else:
    assert json.loads(DATA.read_text()) == full_data, "Coefficient data differ from exact reconstruction"
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
print(encoded)
print("reference_sha256="+hashlib.sha256(encoded.encode()).hexdigest())
