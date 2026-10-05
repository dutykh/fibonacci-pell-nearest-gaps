#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact identities and rational bounds for the finite-word reflection decoder.

Python 3.10+ and SymPy; no research-code imports and no word or directive range.
The only literal mixed controls are d=ab,ba, with one four-letter token word.
Run from the manuscript directory:
  prlimit --as=134217728 --cpu=30 -- timeout 30s python3 -B supplement/checks/check_decoder.py
The proof and the scope of this fixed check are given in Section 3.
"""

import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without optimization and without arguments.")

import hashlib
import json
import sympy as S


identities = []
margins = {}


def same(left, right, name):
    difference = left-right
    if isinstance(difference, S.MatrixBase):
        valid = all(S.cancel(value) == 0 for value in difference)
    else:
        valid = S.cancel(difference) == 0
    if not valid:
        raise AssertionError(name)
    identities.append(name)


def positive(value, name):
    value = S.Rational(value)
    if value <= 0:
        raise AssertionError(name)
    margins[name] = str(value)


def act(matrix, value):
    return S.cancel((matrix[0, 0]*value+matrix[0, 1]) /
                    (matrix[1, 0]*value+matrix[1, 1]))


q = S.Rational
A = S.Matrix([[2, 1], [1, 1]])
B = S.Matrix([[5, 2], [2, 1]])
P = S.Matrix([[2, 1], [1, 0]])
e = S.Matrix([2, 1])
eye = S.eye(2)
g11, g12, g22 = S.symbols("g11 g12 g22")
G = S.Matrix([[g11, g12], [g12, g22]])
M0 = (e.T*G*e)[0]
rho0 = (G*e)[0]
K = 2*e*e.T*G-M0*eye
same(A*B+B*A, 6*e*e.T, "constant block sum")
same(S.trace(A*B*G), 3*M0, "symmetric trace ABG")
same(S.trace(B*A*G), 3*M0, "symmetric trace BAG")
same((P.inv()*A*B*G*P)[1, 0], -M0, "negative lower-left block")
same((P.inv()*B*A*G*P)[1, 0], M0, "positive lower-left block")
same(P.inv()*K*P, S.Matrix([[M0, 2*rho0], [0, -M0]]),
     "ancestor cusp conjugation")
same(K*K, M0*M0*eye, "reflection square")
same(K*e, M0*e, "terminal eigenvector")

# A general cusp basis, with the sole exact parameter relation substituted.
u, ell, alpha, z = S.symbols("u L alpha z")
v = ell+3-u
beta = ell*(3-ell)-alpha
delta = (alpha-beta)/ell
tau = u-v+delta
cu = S.Matrix([[u, u*u-3*u+alpha], [-1, 3-u]])
cv = S.Matrix([[v, 3*v-v*v-beta], [1, 3-v]])
kr = S.Matrix([[1, tau], [0, -1]])
same(cu.det(), alpha, "first normalized determinant")
same(cv.det(), beta, "second normalized determinant")
same((cv*cu)[1, 0], ell, "M0 equals pqL")
same((cu*cv+cv*cu)[1, 1], 0, "parameter-curve entry")
same(S.trace(cu*cv), 3*ell, "product trace in cusp coordinates")
same(cu*cv+cv*cu-3*ell*eye, 3*ell*kr,
     "centered reflection coefficient")
same(kr*cu*cv, cv*cu*kr, "exchange ab to ba")
same(kr*cv*cu, cu*cv*kr, "exchange ba to ab")
same(act(cu, z), -u-alpha/(z+u-3), "first cusp function")
same(act(cv, z), v-beta/(z+3-v), "second cusp function")
same(act(kr, z), -z-tau, "affine reflection")
ef = beta/(delta-alpha/(z+u-3))+v-3
df = alpha/(delta-beta/(z+3-v))+3-u
same(act(cv.inv()*kr*cu, z), ef, "first residual formula")
same(act(cu.inv()*kr*cv, z), df, "second residual formula")
same((cv.inv()*kr*cu)*cv, cu*kr, "first residual return")
same((cu.inv()*kr*cv)*cu, cv*kr, "second residual return")

# The four rearrangements use independent parameters and no curve equation.
aa, bb, ll, ss = S.symbols("aa bb ll ss")
dd = (aa-bb)/ll
eta = -dd
same(aa/dd-ll, ll*bb/(aa-bb), "a terminal positive gap")
same(aa/(dd-bb/ss)-ll, bb*(1+ll/ss)/(dd-bb/ss),
     "a finite positive gap")
same(ll-bb/eta, -ll*aa/(bb-aa), "b terminal negative gap")
same(ll-bb/(eta-aa/ss), -aa*(1+ll/ss)/(eta-aa/ss),
     "b finite negative gap")

# These rational margins are exactly the uniform denominator and gap estimates.
positive(q(12, 5), "first positive denominator")
same(q(5, 2)+3-q(41, 12), q(25, 12),
     "first negative denominator magnitude")
same(q(12, 5)+3-q(21, 8), q(111, 40),
     "second positive denominator magnitude")
same(q(5, 2)-(3-q(5, 2)), 2, "second negative denominator magnitude")
same((1-q(1, 25))/3, q(8, 25), "a delta lower coefficient")
same(q(1, 25)/q(8, 25), q(1, 8), "a terminal ratio")
same(q(8, 25)-q(1, 50), q(3, 10), "a residual denominator margin")
same((1-q(4, 25))/3, q(7, 25), "b eta lower coefficient")
same(q(4, 25)/q(7, 25), q(4, 7), "b terminal ratio")
same(q(7, 25)-q(2, 25), q(1, 5), "b residual denominator margin")
same(q(5, 12)+q(4, 7), q(83, 84), "b central forbidden interval")
positive(q(5, 2)-q(1, 2), "positive cylinder gap coefficient")
positive(q(5, 2)/q(1, 4)-q(5, 12), "negative cylinder gap coefficient")
positive(3-q(40, 111)*q(1, 25)-2, "a finite s exceeds two")
positive(3-q(12, 25)*q(1, 25)-2, "b finite s exceeds two")
positive(-q(5, 2)-(-3+q(12, 25)), "negative cylinder in H")
positive(q(5, 2)-q(40, 111)*q(1, 4)-q(12, 5), "positive cylinder in H")
positive(q(5, 2)-q(83, 84), "b forbidden interval outside H")
same(1-q(2, 5), q(3, 5), "reflection denominator infimum")
same(2-q(5, 12), q(19, 12), "reflection positive lower bound")

lo, hi = q(8, 5), q(29, 12)
for generator, endpoints in (
    (A, (q(21, 13), q(70, 41))),
    (B, (q(50, 21), q(169, 70))),
):
    for endpoint, expected in zip((lo, hi), endpoints):
        same(act(generator, endpoint), expected, "generator interval endpoint")
        if not lo < expected < hi:
            raise AssertionError("interval invariance")
same(2-act(A, z), 1/(z+1), "A sends positive input below two")
same(act(B, z)-2, z/(2*z+1), "B sends positive input above two")
tau1 = S.symbols("tau1")
rt = S.Matrix([[2*tau1-1, 4*(1-tau1)], [tau1, 1-2*tau1]])
same(act(rt, z), 2-(z-2)/(1+tau1*(z-2)),
     "original reflection coordinate")

# Complete fixed empty-directive check, with poles visible in the matrices.
k0 = S.Matrix([[3, 4], [4, -3]])
ee, de = B.inv()*k0*A, A.inv()*k0*B
same(ee, S.Matrix([[0, 5], [5, -9]]), "empty first residual")
same(de, S.Matrix([[9, 5], [5, 0]]), "empty second residual")
same(ee*A, S.Matrix([[5, 5], [1, -4]]), "empty wrong aa map")
same(de*B, S.Matrix([[55, 23], [25, 10]]), "empty wrong bb map")
same(act(ee, 2), 5, "empty first terminal")
same(act(de, 2), q(23, 10), "empty second terminal")
positive(5-hi, "empty first terminal outside I")
positive(4-hi, "empty aa pole beyond I")
positive(q(50, 21)-q(23, 10), "empty second terminal below B")
positive(q(23, 10)-2, "empty second terminal above two")
positive(q(50, 21)-q(111, 50), "empty bb image below B")
same(act(de*B, lo), q(111, 50), "empty bb upper endpoint")
for endpoint in (lo, hi):
    if not act(ee*A, endpoint) < 0:
        raise AssertionError("empty wrong aa negative")
    if not 2 < act(de*B, endpoint) <= q(111, 50):
        raise AssertionError("empty wrong bb gap")
if (de*B).det() >= 0:
    raise AssertionError("empty bb map must be decreasing")


def word_matrix(word):
    result = S.eye(2)
    for letter in word:
        result *= A if letter == "a" else B
    return result


controls = []
for directive, center in (("ab", "aba"), ("ba", "bab")):
    uw, vw = "a", "b"
    for letter in directive:
        if letter == "a":
            vw = vw+uw
        else:
            uw = uw+vw
    uu, vv, gg = word_matrix(uw), word_matrix(vw), word_matrix(center)
    mm = (e.T*gg*e)[0]
    kk = 2*e*e.T*gg-mm*eye
    same(uu*vv, A*B*gg, "mixed literal ab identity")
    same(vv*uu, B*A*gg, "mixed literal ba identity")
    pp, qq = S.trace(uu)/3, S.trace(vv)/3
    ux = (P.inv()*uu*P)[0, 0]/pp
    vx = (P.inv()*vv*P)[0, 0]/qq
    if not 3 <= ux <= q(41, 12) or not q(5, 2) <= vx <= q(21, 8):
        raise AssertionError("mixed coordinate bounds")
    if directive[-1] == "a" and not qq >= 5*pp:
        raise AssertionError("mixed a dominance")
    if directive[-1] == "b" and not pp >= q(5, 2)*qq:
        raise AssertionError("mixed b dominance")
    same(kk*uu*vv*vv*uu*e, mm*vv*uu*uu*vv*e,
         "mixed prescribed abba token vector")
    controls.append({"directive": directive, "center": center,
                     "p": int(pp), "q": int(qq), "M0": int(mm)})

payload = {
    "status": "PASS",
    "scope": "fixed exact decoder identities and rational exclusions; no enumeration",
    "symbolic_and_exact_identities": len(identities),
    "strict_rational_margins": margins,
    "literal_controls": controls,
}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
print(encoded)
print("reference_sha256="+hashlib.sha256(encoded.encode()).hexdigest())
