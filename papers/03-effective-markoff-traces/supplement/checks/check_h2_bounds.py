#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Fixed algebra and rational endpoint controls; no ray or phase search."""

import sys

if not __debug__:
    raise SystemExit("optimization is forbidden")
if len(sys.argv) != 1:
    raise SystemExit("unexpected arguments")

from fractions import Fraction as F
import hashlib
import json

import sympy as sp

checks = []


def identity(name, expression):
    if sp.cancel(expression) != 0:
        raise AssertionError(name)
    checks.append(name)


c, d, b, lam, shift = sp.symbols("c d b lam shift", nonzero=True)
a0 = sp.Matrix([[3 * c, b], [b, 3 * d]])
word = a0 * sp.diag(shift, 1 / shift) * a0 * sp.diag(
    shift * lam, 1 / (shift * lam)
)
identity(
    "literal_h2_four_paths",
    sp.trace(word) / 3
    - 3 * c**2 * shift**2 * lam
    - 3 * d**2 / (shift**2 * lam)
    - b**2 * (lam + 1 / lam) / 3,
)

u, v = sp.symbols("u v", nonzero=True)
ba, bb = (u + 1 / u) / 3, (v + 1 / v) / 3
bc, bd = (u * v + 1 / (u * v)) / 3, (u / v + v / u) / 3
identity("baseline_subtraction", 3 * ba * bb - bc - bd)
identity("baseline_quotient_defect", ba**2 + bb**2 - bc * bd - sp.Rational(4, 9))
ea, eb, ec = sp.symbols("ea eb ec")
actual_new = ((ba + ea) ** 2 + (bb + eb) ** 2) / (bc + ec)
identity(
    "full_quadratic_error_identity",
    (bc + ec) * (actual_new - bd)
    - (2 * ba * ea + ea**2 + 2 * bb * eb + eb**2 - bd * ec + sp.Rational(4, 9)),
)

s, alpha, delta = sp.symbols("s alpha delta", nonzero=True)
identity(
    "h2_forward_coefficient",
    ((s * alpha / 3) * s - (1 + delta) * 3 * (s * alpha / 3) ** 2)
    .subs(delta, 1 / alpha - 1),
)
identity(
    "h2_backward_coefficient",
    ((1 + delta) * alpha**3 / (3 * s * alpha)) * s**3
    - (1 + delta) * 3 * (s * alpha / 3) ** 2,
)

h = 2
scale = 100 * h**2
margins = {
    "positive_delta": F(h) - F(8 * h, 27) - F(16, 625),
    "p64_exceeds_3h": F(4 ** (h + 1) - 3 * h),
    "positive_quotient": 5 * (F(1, 3) - F(5, 256))
    - (F(50, 48) + F(50, 4096) + F(4, 9)),
    "initial_factor": F(5) - F(16384, 3825),
    "positive_terminal": F(1) - F(2, 3) - F(5, 16),
    "coarse_activity": F(1, 16) - F(4 * h**2, scale**2),
    "coarse_mixed": F(1, 16) - F(64 * h**2, 15 * scale**2),
    "opposite_constant": F(1, 16) - F(81, 8192),
    "total_coefficient": F(15 * h**2) - F(8) - F(64 * h**2, 5),
    "initial_label_segment": F(1, 12) - F(20 * h**2, scale**2),
    "bootstrap": F(1, 16) - F(2, 125) - F(256, 10**8),
    "spectral_terminal": F(1) - F(271, 360),
}
if not all(value > 0 for value in margins.values()):
    raise AssertionError("nonpositive endpoint margin")
if margins["positive_quotient"] != F(1303, 18432):
    raise AssertionError("quotient margin mismatch")
if F(320 * h**3, scale**2) != F(2, 125):
    raise AssertionError("bootstrap endpoint mismatch")
if 3 * h**2 + 3 * h - 2 + 3 * h // 2 != 19:
    raise AssertionError("gap endpoint mismatch")
if scale ** (h + 1) != 64000000:
    raise AssertionError("minimum threshold mismatch")

payload = {"identities": checks, "positive_margins": {k: str(v) for k, v in margins.items()}}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
print(encoded)
print("SHA256 " + hashlib.sha256(encoded.encode()).hexdigest())
