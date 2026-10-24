#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Second, independent verification of the margins of Section 7.

1. The recorded census (data/margins-census.txt, written by code/margins_census.gp from the Gaussian sign
   cube) is recomputed with other formulas: each label as the continuant of [2, sigma(z), 2], the occurrence
   pair as the numerator and denominator of the right half of that continued fraction, the product of the
   recorded prime powers, the BPSW test (SymPy) for each prime, every primitive representation from the
   square roots of -1 modulo the label (Hensel lifting, Chinese remainders, Cornacchia descent), the
   Euclidean costs, and the margin; the domain must be every reduced P/Q with Q <= 100 or P + Q <= 160,
   and every margin positive.
2. On the two ladders below 3/5 the competitor is obtained by Gaussian division of the occurrence pair by
   1 + 2i or 1 - 2i and multiplication by the conjugate, with no matrix; its cost exceeds the occurrence
   cost by exactly k at the k-th rung where 5 divides the label, k <= 40.
3. The thresholds: among the recorded occurrences with P + Q <= 160 the largest P + Q with margin at most k
   is 24k + 3 for k = 1, ..., 6, attained only at (9k+1)/(15k+2), which also gives the largest P and Q.
Exit status 0 when every check passes.
"""
import json
from math import gcd, isqrt
from pathlib import Path
import sys

from sympy import isprime
from sympy.ntheory.modular import crt

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

DATA = Path(__file__).resolve().parents[1] / "data" / "margins-census.txt"

def word(P, Q):
    return "".join("b" if (i*P)//Q - ((i-1)*P)//Q else "a" for i in range(1, Q+1))

def cf_value(digits):                 # numerator and denominator of [d0; d1, ..., dn]
    h, h1, k, k1 = 1, 0, 0, 1
    for d in digits:
        h, h1 = d*h + h1, h
        k, k1 = d*k + k1, k
    return h, k

def label_and_pair(P, Q):
    z = word(P, Q)[1:-1]
    full = [2] + [d for c in z for d in ((1, 1) if c == "a" else (2, 2))] + [2]
    m = cf_value(full)[0]
    n = len(z); w = z[(n + 1)//2:]
    half = [d for c in w for d in ((1, 1) if c == "a" else (2, 2))] + [2]
    if n % 2:
        half = [1 if z[n//2] == "a" else 2] + half
    s, r = cf_value(half)
    assert r*r + s*s == m and gcd(r, s) == 1
    return m, (min(r, s), max(r, s))

def cost(x, y):
    a, b, total = max(x, y), min(x, y), 0
    while b:
        total += a//b; a, b = b, a % b
    assert a == 1
    return total - 1

def root_minus_one(p, e):
    c = 2
    while pow(c, (p - 1)//2, p) != p - 1:
        c += 1
    r, pk = pow(c, (p - 1)//4, p), p
    for _ in range(1, e):
        pk *= p
        r = (r - (r*r + 1)*pow(2*r, -1, pk)) % pk
    assert (r*r + 1) % pk == 0
    return r, pk

def representations(m, factors):
    mods, roots = [], []
    for p, e in factors:
        if p == 2:
            assert e == 1; mods.append(2); roots.append([1])
        else:
            assert p % 4 == 1; r, pk = root_minus_one(p, e); mods.append(pk); roots.append([r, pk - r])
    reps = set()
    def combine(i, chosen):
        if i == len(roots):
            u = int(crt(mods, chosen)[0]) % m
            a, b, s = m, u, isqrt(m)
            while b > s:
                a, b = b, a % b
            y = isqrt(m - b*b)
            assert y*y == m - b*b and gcd(b, y) == 1
            reps.add((min(b, y), max(b, y)))
            return
        for r in roots[i]:
            combine(i + 1, chosen + [r])
    combine(0, [])
    return reps

# 1. the census
lines = [ln.split() for ln in DATA.read_text().splitlines() if ln.strip()]
expected = {(P, Q) for Q in range(2, 160) for P in range(1, Q)
            if gcd(P, Q) == 1 and (Q <= 100 or P + Q <= 160)}
seen, margins = set(), {}
for P, Q, margin, fact in lines:
    P, Q = int(P), int(Q)
    assert (P, Q) in expected and (P, Q) not in seen
    seen.add((P, Q))
    m, occ = label_and_pair(P, Q)
    factors = []
    for item in fact.split("*"):
        p, _, e = item.partition("^")
        factors.append((int(p), int(e) if e else 1))
    product = 1
    for p, e in factors:
        assert isprime(p)
        product *= p**e
    assert product == m
    reps = representations(m, factors)
    assert len(reps) == 2**(sum(1 for p, _ in factors if p != 2) - 1) and occ in reps
    assert cost(*occ) == P + Q - 2
    others = [cost(*r) for r in reps if r != occ]
    computed = min(others) - (P + Q - 2) if others else None
    assert (margin == "-" and computed is None) or (margin != "-" and int(margin) == computed)
    assert computed is None or computed > 0
    margins[(P, Q)] = computed
assert seen == expected

# 2. the two ladders, by Gaussian division
def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gdiv_exact(x, y):
    n = y[0]**2 + y[1]**2
    num = gmul(x, (y[0], -y[1]))
    if num[0] % n or num[1] % n:
        return None
    return (num[0]//n, num[1]//n)
for P0, Q0, first in ((1, 2, 3), (2, 3, 2)):
    for k in range(1, 41):
        j = first + 3*(k - 1)
        P, Q = P0 + 3*j, Q0 + 5*j
        m, occ = label_and_pair(P, Q)
        assert m % 5 == 0
        z = (occ[1], occ[0])
        comps = []
        for g in ((1, 2), (1, -2)):
            q = gdiv_exact(z, g)
            if q is not None:
                comps.append(gmul(q, (g[0], -g[1])))
        assert len(comps) == 1
        c = (abs(comps[0][0]), abs(comps[0][1]))
        assert c[0]**2 + c[1]**2 == m and gcd(*c) == 1
        assert cost(*c) - cost(*occ) == k
        if (P, Q) in margins:
            assert margins[(P, Q)] <= k

# 3. the thresholds
for k in range(1, 7):
    small = [(P, Q) for (P, Q), mg in margins.items() if P + Q <= 160 and mg is not None and mg <= k]
    assert max(P + Q for P, Q in small) == 24*k + 3
    assert [pq for pq in small if pq[0] + pq[1] == 24*k + 3] == [(9*k + 1, 15*k + 2)]
    assert max(P for P, _ in small) == 9*k + 1 and max(Q for _, Q in small) == 15*k + 2
    assert margins[(9*k + 1, 15*k + 2)] == k
count = {c: sum(1 for mg in margins.values() if (mg is None and c == "single") or
                (mg is not None and {"1-2": mg <= 2, "3-9": 3 <= mg <= 9, "10-99": 10 <= mg <= 99,
                                     "100+": mg >= 100}.get(c, False)))
         for c in ("1-2", "3-9", "10-99", "100+", "single")}
print(json.dumps({"status": "PASS", "occurrences": len(margins), "by_margin": count,
                  "ladder_rungs": 80, "thresholds_checked": 6}, sort_keys=True))
