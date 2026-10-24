#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""First verification of the labels of a row and of the source comparisons that replace the import.

Exact arithmetic in Q(sqrt5), numbers written a + b*sqrt5 with rational a, b; signs by the integer test.
1. The eigenbasis identities of the proposition on the labels of a row: with V = [[phi, 1], [1, -phi]],
   V^T mu(a) V = (phi^2+1) diag(phi^2, phi^-2) and V^T mu(b) V = (phi^2+1) Bhat, with phi^2 + 1 = sqrt5 phi.
2. Controls of the trace formula m = (1/3) tr prod diag(phi^(2u_k), phi^(-2u_k)) Bhat, for P <= 13, Q <= 80,
   and of the two-sided bound C_P < rho < C_P (1 + phi^(-4n-6))^P, or rho < C_1 when P = 1, for P <= 13, Q <= 120.
   These controls illustrate statements proved in the text; they are not part of the proof.
3. The comparisons used in the proof: C_P (1 + phi^(-4 n_P - 6))^P < U_P for P = 8..13 with
   n_P = 3, 4, 4, 7, 8, 11 (and the failure at n_P - 1, which shows each n_P is the least possible),
   C_7 (1 + phi^-10)^7 < 1, the 38 integer comparisons m_{P/Q} < U_P F_{2(P+Q-2)+1} for the reduced Q with
   Q_P <= Q < P(n_P + 1), and the six comparisons m_{7/Q} < F_{2Q+11}, 8 <= Q <= 13.
Standard library only. Exit status 0 when every check passes.
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
from source_words import source, fibonacci, cohn_word

def add(x, y): return x[0] + y[0], x[1] + y[1]
def sub(x, y): return x[0] - y[0], x[1] - y[1]
def mul(x, y): return x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0]
def scal(c, x): return c*x[0], c*x[1]
def inv(x):
    n = x[0]*x[0] - 5*x[1]*x[1]
    assert n != 0
    return x[0]/n, -x[1]/n
def power(x, k):
    out = (F(1), F(0))
    base = x if k >= 0 else inv(x)
    for _ in range(abs(k)):
        out = mul(out, base)
    return out
def sign(x):
    a, b = x
    if b == 0: return (a > 0) - (a < 0)
    if a == 0: return (b > 0) - (b < 0)
    if (a > 0) == (b > 0): return 1 if a > 0 else -1
    d = a*a - 5*b*b
    assert d != 0
    return (1 if a > 0 else -1) * (1 if d > 0 else -1)
def less(x, y): return sign(sub(y, x)) > 0
def rat(n, d=1): return F(n, d), F(0)

PHI = (F(1, 2), F(1, 2))
ONE, ZERO = rat(1), rat(0)
SQRT5 = (F(0), F(1))
LAMBDA = mul(rat(3), inv(mul(SQRT5, PHI)))            # 3/(sqrt5 phi)
def C(P): return mul(power(PHI, 2), power(LAMBDA, P - 1))

def mat_mul(X, Y):
    return [[add(mul(X[i][0], Y[0][j]), mul(X[i][1], Y[1][j])) for j in range(2)] for i in range(2)]
def transpose(X): return [[X[0][0], X[1][0]], [X[0][1], X[1][1]]]
def int_mat(rows): return [[rat(v) for v in row] for row in rows]

# 1. Eigenbasis identities
V = [[PHI, ONE], [ONE, scal(-1, PHI)]]
norm2 = add(power(PHI, 2), ONE)
assert norm2 == mul(SQRT5, PHI)
lhs_a = mat_mul(mat_mul(transpose(V), int_mat([[2, 1], [1, 1]])), V)
assert lhs_a == [[mul(norm2, power(PHI, 2)), ZERO], [ZERO, mul(norm2, power(PHI, -2))]]
BHAT = [[scal(F(1, 5), mul(SQRT5, scal(3, power(PHI, 3)))), scal(F(1, 5), scal(2, SQRT5))],
        [scal(F(1, 5), scal(2, SQRT5)), scal(F(1, 5), mul(SQRT5, scal(3, power(PHI, -3))))]]
lhs_b = mat_mul(mat_mul(transpose(V), int_mat([[5, 2], [2, 1]])), V)
assert lhs_b == [[mul(norm2, BHAT[i][j]) for j in range(2)] for i in range(2)]
assert add(BHAT[0][0], BHAT[1][1]) == rat(6)
assert sub(mul(BHAT[0][0], BHAT[1][1]), mul(BHAT[0][1], BHAT[1][0])) == ONE
assert all(sign(BHAT[i][j]) > 0 for i in range(2) for j in range(2))

def runs(P, Q):
    word, out, count = cohn_word(P, Q), [], 0
    for letter in word:
        if letter == "a": count += 1
        else: out.append(count); count = 0
    assert len(out) == P and sum(out) == Q - P
    n = (Q - P) // P
    assert all(n <= u <= n + 1 for u in out)
    return out, n

# 2. Controls of the trace formula and of the two-sided bound
expansion_controls = 0
for Q in range(2, 81):
    for P in range(1, min(Q, 14)):
        if gcd(P, Q) != 1: continue
        u, _ = runs(P, Q)
        T = [[ONE, ZERO], [ZERO, ONE]]
        for uk in u:
            T = mat_mul(mat_mul(T, [[power(PHI, 2*uk), ZERO], [ZERO, power(PHI, -2*uk)]]), BHAT)
        assert scal(F(1, 3), add(T[0][0], T[1][1])) == rat(source(P, Q)[0])
        expansion_controls += 1
bound_controls = 0
for Q in range(2, 121):
    for P in range(1, min(Q, 14)):
        if gcd(P, Q) != 1: continue
        _, n = runs(P, Q)
        rho = rat(source(P, Q)[0], fibonacci(2*(P + Q - 2) + 1))
        if P == 1:
            assert less(rho, C(1))
        else:
            assert less(C(P), rho)
            assert less(rho, mul(C(P), power(add(ONE, power(PHI, -4*n - 6)), P)))
        bound_controls += 1

# 3. The comparisons used in the proof
N_P = {8: 3, 9: 4, 10: 4, 11: 7, 12: 8, 13: 11}
Q_P = {8: 29, 9: 35, 10: 39, 11: 76, 12: 101, 13: 146}
U_P = {8: (707, 1000), 9: (2929, 5000), 10: (4857, 10000), 11: (20111, 50000),
       12: (333511, 10**6), 13: (2765404267, 10**10)}
corner = {}
for P, n in N_P.items():
    bound = lambda k: mul(C(P), power(add(ONE, power(PHI, -4*k - 6)), P))
    assert less(bound(n), rat(*U_P[P]))
    assert not less(bound(n - 1), rat(*U_P[P]))
    corner[P] = n
assert less(mul(C(7), power(add(ONE, power(PHI, -10)), 7)), ONE)
assert not less(mul(C(7), power(add(ONE, power(PHI, -6)), 7)), ONE)
junctions = []
for P in range(8, 14):
    for Q in range(Q_P[P], P*(N_P[P] + 1)):
        if gcd(P, Q) != 1: continue
        m, pair = source(P, Q)
        f = fibonacci(2*(P + Q - 2) + 1)
        num, den = U_P[P]
        margin = num*f - den*m
        assert margin > 0
        junctions.append({"P": P, "Q": Q, "M": m, "positive_margin": margin})
assert len(junctions) == 38
assert [r["positive_margin"] for r in junctions if (r["P"], r["Q"]) in Q_P.items()] == \
    [r["positive_margin"] for r in json.loads((ROOT / "data/source-junctions.json").read_text())]
assert junctions[0]["positive_margin"] == 25054682255203
row_seven = []
for Q in range(8, 14):
    m, _ = source(7, Q)
    f = fibonacci(2*Q + 11)
    assert m < f
    row_seven.append([Q, m, f])
assert row_seven[0][1:] == [195025, 196418]
payload = {"junctions": junctions, "row_seven": row_seven, "corner_n": corner}
digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
print(json.dumps({"status": "PASS", "expansion_controls": expansion_controls,
                  "bound_controls": bound_controls, "corner_comparisons": len(corner) + 1,
                  "junction_comparisons": len(junctions), "row_seven_comparisons": len(row_seven),
                  "reference_sha256": digest}, sort_keys=True))
