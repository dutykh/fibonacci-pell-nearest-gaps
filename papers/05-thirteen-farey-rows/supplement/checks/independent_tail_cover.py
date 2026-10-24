#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Second, independent verification of parts (ii) and (iii) of the finite verifications.

Written from the text of the article alone, with a different formula for the normalized
ratio: numbers a + b sqrt5 with integer a, b, and
    R = sqrt5 phi^(-2k-5) |W|^2 / (1 - q^(2k+5+2S)),
    W = (Pi+ + x_k Pi-) A ... A (Pi+ + x_0 Pi-) z0,   x_i = q^g_i,
where Pi+ and Pi- are the spectral projectors of B. The root boxes, matched corners and
splitting rule are those of the section on covers. Every tree is rebuilt, every leaf is
required to close, and the counts are compared with tail/data/catalog.json.
Standard library only. Exit status 0 when every check passes.
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

from math import gcd

def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def mul(x, y): return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def smul(n, x): return (n*x[0], n*x[1])

def sign(x):
    a, b = x
    if b == 0: return (a > 0) - (a < 0)
    if a == 0: return (b > 0) - (b < 0)
    if (a > 0) == (b > 0): return 1 if a > 0 else -1
    return (1 if a > 0 else -1) * (1 if a*a - 5*b*b > 0 else -1)

ZERO, ONE, R5 = (0, 0), (1, 0), (0, 1)

@lru_cache(maxsize=None)
def qnum(n):                      # (-3 + sqrt5)^n, so q^n = qnum(n) / 2^n
    r = ONE
    for _ in range(n): r = mul(r, (-3, 1))
    return r

@lru_cache(maxsize=None)
def one_plus_r5_pow(n):
    r = ONE
    for _ in range(n): r = mul(r, (1, 1))
    return r

A = ((ONE, ZERO), (ONE, ONE))
MP = (((-1, 1), (2, 0)), ((2, 0), (1, 1)))     # 2 sqrt5 Pi+ = 2B - (1 - sqrt5) I
MM = (((1, 1), (-2, 0)), ((-2, 0), (-1, 1)))   # 2 sqrt5 Pi- = (1 + sqrt5) I - 2B

def mat_vec(M, v):
    return (add(mul(M[0][0], v[0]), mul(M[0][1], v[1])),
            add(mul(M[1][0], v[0]), mul(M[1][1], v[1])))

def factor_apply(g, v):           # (2^g M+ + Qg M-) v ; M+ v for an infinite gap
    if g is None: return mat_vec(MP, v)
    a, b, Qg = mat_vec(MP, v), mat_vec(MM, v), qnum(g)
    return (add(smul(2**g, a[0]), mul(Qg, b[0])), add(smul(2**g, a[1]), mul(Qg, b[1])))

@lru_cache(maxsize=1 << 16)
def ratio_parts(gaps):
    k = len(gaps) - 1
    v = ((1, 0), (2, 0))
    for i, g in enumerate(gaps):
        v = factor_apply(g, v)
        if i < k: v = mat_vec(A, v)
    w2 = add(mul(v[0], v[0]), mul(v[1], v[1]))
    sfin = sum(g for g in gaps if g is not None)
    num = mul(R5, smul(2**(2*k + 5), w2))
    den = smul(20**(k + 1) * 4**sfin, one_plus_r5_pow(2*k + 5))
    if all(g is not None for g in gaps):
        n = 2*k + 5 + 2*sfin
        num, den = smul(2**n, num), mul(den, sub((2**n, 0), qnum(n)))
    assert sign(den) > 0 and sign(num) > 0
    return num, den

def compare(gaps, t):             # sign(R - tnum/tden)
    num, den = ratio_parts(gaps)
    return sign(sub(smul(t[1], num), mul(t[0], den)))

def C_threshold(P):               # phi^2 = (3+sqrt5)/2, 3/(sqrt5 phi) = (15-3sqrt5)/10
    num, den = (3, 1), 2
    for _ in range(P - 1): num, den = mul(num, (15, -3)), den * 10
    g = gcd(gcd(abs(num[0]), abs(num[1])), den)
    return (num[0]//g, num[1]//g), den//g

U = {8: ((707, 0), 1000), 9: ((2929, 0), 5000), 10: ((4857, 0), 10000),
     11: ((20111, 0), 50000), 12: ((333511, 0), 10**6), 13: ((2765404267, 0), 10**10)}
CP = {8: 34, 9: 41, 10: 46, 11: 84, 12: 110, 13: 156}

def corner(box, maximize):
    k = len(box) - 1
    return tuple(hi if ((i < k and lo % 2 == 0) == maximize) else lo
                 for i, (lo, hi) in enumerate(box))

def roots(k, S0, parity):
    H = -(-S0 // (k + 1)); assert H >= 2
    for f in range(k + 1):
        box = []
        for i, b in enumerate(parity):
            if i < f: box.append((b, H - 1 if (H - 1) % 2 == b else H - 2))
            elif i == f: box.append((H if H % 2 == b else H + 1, None))
            else: box.append((b, None))
        yield tuple(box)

def split(box):
    lo, i = min((lo, i) for i, (lo, hi) in enumerate(box) if hi is None or hi > lo)
    hi = box[i][1]
    if hi is None: left, right = (lo, lo), (lo + 2, None)
    else:
        cut = lo + 2*(((hi - lo)//2 + 1)//2)
        left, right = (lo, cut - 2), (cut, hi)
    b1, b2 = list(box), list(box); b1[i], b2[i] = left, right
    return tuple(b1), tuple(b2)

def cover(k, S0, lower, upper, below_only=False):
    st = dict(nodes=0, below=0, above=0, ineligible=0, roots=0)
    for sector in range(2**(k + 1)):
        stack = list(roots(k, S0, tuple((sector >> (k - i)) & 1 for i in range(k + 1))))
        st["roots"] += len(stack)
        while stack:
            box = stack.pop(); st["nodes"] += 1
            if all(hi is not None for _, hi in box) and sum(hi for _, hi in box) < S0:
                st["ineligible"] += 1
            elif compare(corner(box, True), lower) < 0:
                st["below"] += 1
            elif compare(corner(box, False), upper) > 0:
                assert not below_only; st["above"] += 1
            else:
                stack.extend(split(box))      # raises on an unresolved singleton
    return st

groups = {}
total = dict(nodes=0, below=0, above=0, ineligible=0, roots=0)
for P in range(8, 14):
    for k in range(2, P - 4):                 # middle 2..P-6, final P-5
        st = cover(k, CP[P] - k - 1, C_threshold(P), U[P], below_only=(k == P - 5))
        groups[f"row{P}-k{k}"] = st
        for key in total: total[key] += st[key]
st = cover(2, 152, ((71, 0), 100), ((71, 0), 100), below_only=True)
groups["lower-cost-two-defects"] = st
for key in total: total[key] += st[key]
assert total["ineligible"] == 0
catalog = json.loads((Path(__file__).resolve().parents[1] / "tail" / "data" / "catalog.json").read_text())
recorded = {}
for g in catalog["groups"]:
    name = "lower-cost-two-defects" if g["file"].startswith("lower-cost") else f"row{g['row']}-k{g['defects']}"
    recorded[name] = {key: g[key] for key in ("nodes", "below", "above", "ineligible", "roots")}
assert recorded == groups, "independent cover differs from the recorded catalog"
print(json.dumps({"status": "PASS", "groups": len(groups), "total": total,
                  "agrees_with_catalog": True}, sort_keys=True))
