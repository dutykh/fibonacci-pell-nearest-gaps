#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Second, independent verification of part (i) of the finite verifications.

Labels from the floor increments of the lower Christoffel word; the occurrence pair from
the central Gram factor; primality of the recorded factors by the BPSW test (SymPy), not
by the Lucas chains of checks/verify_gaussian.py; every primitive representation by the
Cornacchia descent from each square root of -1 modulo the label, not by the Gaussian sign
cube. The occurrence pair must be the unique representation of minimum cost P + Q - 2.
It also checks, with its own code, the recorded Lucas witnesses of every prime factor of every label,
recursively down to the prime two, and the six junction comparisons of part (iv), with labels taken as
continuants of the palindromic continued fraction [2, sigma(z), 2] and Fibonacci numbers by doubling.
Exit status 0 when every check passes.
"""
import json
import sys

from itertools import product
from math import gcd, isqrt
from pathlib import Path
from sympy import isprime
from sympy.ntheory.modular import crt

def mm(X, Y):
    return ((X[0][0]*Y[0][0]+X[0][1]*Y[1][0], X[0][0]*Y[0][1]+X[0][1]*Y[1][1]),
            (X[1][0]*Y[0][0]+X[1][1]*Y[1][0], X[1][0]*Y[0][1]+X[1][1]*Y[1][1]))
MA, MB, P1, P2, I2 = ((2,1),(1,1)), ((5,2),(2,1)), ((1,1),(1,0)), ((2,1),(1,0)), ((1,0),(0,1))
def word(P, Q): return ''.join('b' if (i*P)//Q - ((i-1)*P)//Q else 'a' for i in range(1, Q+1))
def mu(w):
    M = I2
    for c in w: M = mm(M, MA if c == 'a' else MB)
    return M
def occurrence_pair(P, Q):
    z = word(P, Q)[1:-1]; n = len(z); H = mu(z[(n+1)//2:])
    if n % 2: H = mm(P1 if z[n//2] == 'a' else P2, H)
    return tuple(sorted((2*H[0][0] + H[0][1], 2*H[1][0] + H[1][1])))
def cost(x, y):
    a, b, s = max(x, y), min(x, y), 0
    while b: s += a//b; a, b = b, a % b
    return s - 1
def root_minus_one(p, e):
    c = 2
    while pow(c, (p-1)//2, p) != p-1: c += 1
    r, pk = pow(c, (p-1)//4, p), p
    for _ in range(1, e):
        pk *= p; r = (r - (r*r+1)*pow(2*r, -1, pk)) % pk
    return r, pk
def fibre(m, factors):
    mods, roots = [], []
    for p, e in factors:
        if p == 2: mods.append(2); roots.append([1])
        else:
            assert p % 4 == 1; r, pk = root_minus_one(p, e); mods.append(pk); roots.append([r, pk-r])
    reps = set()
    for ch in product(*roots):
        u = int(crt(mods, list(ch))[0]) % m
        a, b, s = m, u, isqrt(m)
        while b > s: a, b = b, a % b
        y = isqrt(m - b*b); assert y*y == m - b*b and gcd(b, y) == 1
        reps.add(tuple(sorted((b, y))))
    return reps
if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")
DATA = Path(__file__).resolve().parents[1] / "data" / "gaussian"; total = [0, 0]
for P, Qmax in {8: 28, 9: 34, 10: 38, 11: 75, 12: 100, 13: 145}.items():
    for Q in range(P+1, Qmax+1):
        if gcd(P, Q) != 1: continue
        rec = json.loads((DATA / f"{P:02d}-{Q:03d}.json").read_text())
        M = mu(word(P, Q)); m = (M[0][0] + M[1][1]) // 3; assert rec["M"] == m
        factors = [(int(p), int(e)) for p, e in rec["factors"]]
        assert all(isprime(p) for p, _ in factors)
        prod = 1
        for p, e in factors: prod *= p**e
        assert prod == m
        reps = fibre(m, factors)
        assert len(reps) == 2**(sum(1 for p, _ in factors if p != 2) - 1)
        costs = sorted((cost(*r), r) for r in reps); occ = occurrence_pair(P, Q)
        assert costs[0] == (P+Q-2, occ) and (len(costs) == 1 or costs[1][0] > P+Q-2)
        total[0] += 1; total[1] += len(reps)

# Second verification of the primality proofs: the Lucas criterion with the recorded witnesses, recomputed here
# with its own code, recursively from the prime factors of every label down to the prime two.
LUCAS = {int(p): rec for p, rec in json.loads((DATA.parent / "lucas_primes.json").read_text()).items()}
proved_here = {2}
def lucas_proved(n, depth=0):
    if n in proved_here:
        return True
    assert depth < 200 and n in LUCAS and n > 2 and n % 2 == 1
    rec = LUCAS[n]
    a = int(rec["witness"])
    factors = [(int(r), int(e)) for r, e in rec["factors"]]
    product = 1
    for r, e in factors:
        product *= r ** e
    assert product == n - 1 and 1 < a < n
    assert pow(a, n - 1, n) == 1
    for r, _ in factors:
        assert gcd(pow(a, (n - 1) // r, n) - 1, n) == 1
        assert lucas_proved(r, depth + 1)
    proved_here.add(n)
    return True
label_primes = set()
for path in sorted(DATA.glob("*.json")):
    for p, _ in json.loads(path.read_text())["factors"]:
        label_primes.add(int(p))
assert all(lucas_proved(p) for p in label_primes)

# Second verification of the six junction comparisons, part (iv): labels as continuants of the palindromic
# continued fraction [2, sigma(z), 2] (Canakci--Schiffler), Fibonacci numbers by doubling.
def continuant_label(P, Q):
    w = word(P, Q)
    digits = [2] + [d for c in w[1:-1] for d in ((1, 1) if c == "a" else (2, 2))] + [2]
    h, h_prev = 1, 0                  # numerator of the continued fraction, by the forward recurrence
    for d in digits:
        h, h_prev = d * h + h_prev, h
    return h
def fib_doubling(n):
    def pair(m):                      # (F(m), F(m+1))
        if m == 0:
            return 0, 1
        f, g = pair(m // 2)
        c, d = f * (2 * g - f), f * f + g * g
        return (d, c + d) if m % 2 else (c, d)
    return pair(n)[0]
CEILINGS = {8: (707, 1000), 9: (2929, 5000), 10: (4857, 10000), 11: (20111, 50000),
            12: (333511, 10**6), 13: (2765404267, 10**10)}
JUNCTIONS = {8: 29, 9: 35, 10: 39, 11: 76, 12: 101, 13: 146}
margins = {}
for P, Qp in JUNCTIONS.items():
    m = continuant_label(P, Qp)
    M = mu(word(P, Qp)); assert m == (M[0][0] + M[1][1]) // 3
    num, den = CEILINGS[P]
    margin = num * fib_doubling(2 * (P + Qp - 2) + 1) - den * m
    assert margin > 0
    margins[P] = margin
assert margins[8] == 25054682255203
print(json.dumps({"status": "PASS", "sources": total[0], "primitive_pairs": total[1],
                  "primes_proved": len(proved_here), "junction_comparisons": len(margins)}, sort_keys=True))
