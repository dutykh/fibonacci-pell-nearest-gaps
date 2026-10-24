#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""First verification of the finite part of the proposition on the two ladders below 3/5.

The two ladders are the Farey fractions (1+3j)/(2+5j) and (2+3j)/(3+5j); their lower Christoffel words
are ab (ababb)^j and abb (ababb)^j.  Checked, in exact integer arithmetic:
1. the central words satisfy z_{j+2} = rev(u) z_j u, with u = babab on the first ladder and u = abbab on
   the second (for every j up to 40, a control of the identity proved in the text), so that
   H_{j+2} = H_j mu(u) for the central Gram factors;
2. for each ladder and each parity of j the transfer T = H_j mu(u) H_j^(-1), which therefore depends on
   the parity only, and the occurrence vectors satisfy v_{j+2} = T v_j;
3. on the rungs where 5 divides the label, the rotation R (multiplication by (-3-4i)/5 or (-3+4i)/5,
   whichever gives an integral vector) and a sign matrix S give a positive competitor S R v_j, and the
   conjugated six-rung transfer N = S R T^3 R^(-1) S^(-1) is a nonnegative integer matrix of determinant 1
   whose word in L_a = [[1,1],[0,1]], L_b = [[1,0],[1,1]] has 50 letters, while T^3 has 48;
4. the base cases: the competitor costs exactly 1 and 2 more than the occurrence at the first two rungs
   of each family, and, as a control, exactly k more at the k-th rung for k <= 40.
Standard library only. Exit status 0 when every check passes.
"""
from fractions import Fraction as Fr
import json
from math import gcd
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

A, B, P1, P2, I2 = ((2, 1), (1, 1)), ((5, 2), (2, 1)), ((1, 1), (1, 0)), ((2, 1), (1, 0)), ((1, 0), (0, 1))
def mm(X, Y):
    return tuple(tuple(sum(X[i][k]*Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def mv(X, v): return (X[0][0]*v[0] + X[0][1]*v[1], X[1][0]*v[0] + X[1][1]*v[1])
def det(X): return X[0][0]*X[1][1] - X[0][1]*X[1][0]
def inv(X):
    d = det(X); assert d in (1, -1)
    return ((X[1][1]*d, -X[0][1]*d), (-X[1][0]*d, X[0][0]*d))
def mu(w):
    M = I2
    for c in w: M = mm(M, A if c == "a" else B)
    return M
def christoffel(P, Q): return "".join("b" if (i*P)//Q - ((i-1)*P)//Q else "a" for i in range(1, Q+1))
def gram(P, Q):
    z = christoffel(P, Q)[1:-1]; n = len(z); H = mu(z[(n+1)//2:])
    if n % 2: H = mm(P1 if z[n//2] == "a" else P2, H)
    return H, z
def occurrence(P, Q):
    H, _ = gram(P, Q)
    return mv(H, (2, 1))
def cost(x, y):
    a, b, s = max(abs(x), abs(y)), min(abs(x), abs(y)), 0
    while b: s += a//b; a, b = b, a % b
    assert a == 1
    return s - 1
def monoid_length(N):
    """Number of letters of N in the free monoid on L_a, L_b; None if N is not in it."""
    if det(N) != 1 or min(min(r) for r in N) < 0: return None
    count = 0
    while N != I2:
        (a, b), (c, d) = N
        if a >= c and b >= d: N = ((a - c, b - d), (c, d))
        elif c >= a and d >= b: N = ((a, b), (c - a, d - b))
        else: return None
        if min(min(r) for r in N) < 0: return None
        count += 1
    return count
LADDERS = {"1/2": (1, 2, 0, "babab"), "2/3": (2, 3, 2, "abbab")}   # (P0, Q0, residue of j with 5 | label, u)
ROTATIONS = (((-3, 4), (-4, -3)), ((-3, -4), (4, -3)))   # 5 times multiplication by (-3-4i)/5, (-3+4i)/5
SIGNS = (((1, 0), (0, 1)), ((-1, 0), (0, -1)), ((1, 0), (0, -1)), ((-1, 0), (0, 1)))
records = {}
for name, (P0, Q0, res, U) in LADDERS.items():
    rung = lambda j: (P0 + 3*j, Q0 + 5*j)
    for j in range(1, 41):                           # 1. the word identity
        z_j, z_j2 = gram(*rung(j))[1], gram(*rung(j + 2))[1]
        assert z_j2 == U[::-1] + z_j + U
        assert gram(*rung(j + 2))[0] == mm(gram(*rung(j))[0], mu(U))
    transfers = {}
    for parity in (0, 1):                            # 2. one transfer per parity
        j0 = 2 if parity == 0 else 1
        H = gram(*rung(j0))[0]
        T = mm(mm(H, mu(U)), inv(H))
        for j in range(j0, 41, 2):
            assert mv(T, occurrence(*rung(j))) == occurrence(*rung(j + 2))
        transfers[parity] = T
    lengths = []
    for j in (res + 3, res + 6):                     # 3. the six-rung identity on the rungs of each parity
        if j < 1: continue
        T = transfers[j % 2]; T3 = mm(mm(T, T), T)
        v0, v6 = occurrence(*rung(j)), occurrence(*rung(j + 6))
        found = False
        for R0 in ROTATIONS:
            if any(x % 5 for x in mv(R0, v0)): continue
            for S in SIGNS:
                c0 = tuple(x // 5 for x in mv(mm(S, R0), v0))
                c6 = tuple(x // 5 for x in mv(mm(S, R0), v6))
                if min(c0) <= 0 or min(c6) <= 0: continue
                # N = S R T^3 R^-1 S^-1 with R = R0/5; R0^-1 = R0^T / 25
                R0T = ((R0[0][0], R0[1][0]), (R0[0][1], R0[1][1]))
                N25 = mm(mm(mm(mm(S, R0), T3), R0T), S)
                assert all(x % 25 == 0 for r in N25 for x in r)
                N = tuple(tuple(x // 25 for x in r) for r in N25)
                assert mv(N, c0) == c6
                assert monoid_length(N) == 50 and monoid_length(T3) == 48
                assert cost(*c6) - cost(*c0) == 50 and cost(*v6) - cost(*v0) == 48
                lengths.append([j, j % 2, 50, 48]); found = True
        assert found, (name, j)
    assert len(lengths) == 2 and {x[1] for x in lengths} == {0, 1}
    excess = []                                      # 4. base cases and control
    for k in range(1, 41):
        j = res + 3*k if res else 3*k
        if name == "2/3": j = 3*(k - 1) + 2
        P, Q = rung(j); v = occurrence(P, Q); m = v[0]**2 + v[1]**2
        assert m % 5 == 0 and cost(*v) == P + Q - 2
        comps = [tuple(x // 5 for x in mv(R0, v)) for R0 in ROTATIONS if not any(x % 5 for x in mv(R0, v))]
        assert len(comps) == 1
        c = comps[0]; assert c[0]**2 + c[1]**2 == m and gcd(*c) == 1
        excess.append(cost(*c) - cost(*v))
    assert excess == list(range(1, 41))
    records[name] = {"six_rung_identities": lengths, "excess_first_two": excess[:2],
                     "transfers": {str(p): [list(r) for r in transfers[p]] for p in (0, 1)}}
assert records["1/2"]["transfers"]["0"] == [[1105, 463], [463, 194]] or records["1/2"]["transfers"]["1"] == [[1105, 463], [463, 194]]
print(json.dumps({"status": "PASS", "ladders": 2, "six_rung_identities": 4,
                  "excess_checked_up_to": 40, "records": records}, sort_keys=True))
