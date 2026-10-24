#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Second, independent verification of the source comparisons beyond the six junctions.

Written from the text of the article, with formulas different from those of checks/check_row_bounds.py:
labels are continuants of the palindromic continued fraction [2, sigma(z), 2] of Canakci and Schiffler,
Fibonacci numbers come from the doubling formulas, and the constants are bounded with rigorous rational
enclosures of sqrt5 (sixty decimal digits) instead of the exact sign test in Q(sqrt5).
Checked: m_{P/Q} < U_P F_{2(P+Q-2)+1} for the 38 reduced Q with Q_P <= Q < P(n_P + 1), P = 8..13;
m_{7/Q} < F_{2Q+11} for 8 <= Q <= 13; C_P (1 + phi^(-4 n_P - 6))^P < U_P for P = 8..13 and
C_7 (1 + phi^-10)^7 < 1, from upper enclosures; and, from lower enclosures, that n_P - 1 would not do.
Standard library only. Exit status 0 when every check passes.
"""
from fractions import Fraction as Fr
import json
from math import gcd, isqrt
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

def word(P, Q):
    return "".join("b" if (i*P)//Q - ((i-1)*P)//Q else "a" for i in range(1, Q+1))

def continuant_label(P, Q):
    w = word(P, Q)
    digits = [2] + [d for c in w[1:-1] for d in ((1, 1) if c == "a" else (2, 2))] + [2]
    h, h_prev = 1, 0
    for d in digits:
        h, h_prev = d*h + h_prev, h
    return h

def fib(n):
    def pair(k):                      # (F(k), F(k+1))
        if k == 0:
            return 0, 1
        f, g = pair(k // 2)
        c, d = f*(2*g - f), f*f + g*g
        return (d, c + d) if k % 2 else (c, d)
    return pair(n)[0]

DIGITS = 60
s = isqrt(5 * 10**(2*DIGITS))
SQRT5_LO, SQRT5_HI = Fr(s, 10**DIGITS), Fr(s + 1, 10**DIGITS)
assert SQRT5_LO**2 < 5 < SQRT5_HI**2
PHI_LO, PHI_HI = (1 + SQRT5_LO)/2, (1 + SQRT5_HI)/2

def enclosure(P, n):
    """Lower and upper rational bounds for C_P (1 + phi^(-4n-6))^P, C_P = phi^(3-P) (3/sqrt5)^(P-1)."""
    def value(phi, root5):
        return phi**(3 - P) * (Fr(3)/root5)**(P - 1) * (1 + phi**(-4*n - 6))**P
    # For P >= 3 every factor decreases when phi or sqrt5 increases.
    assert P >= 3
    return value(PHI_HI, SQRT5_HI), value(PHI_LO, SQRT5_LO)

N_P = {8: 3, 9: 4, 10: 4, 11: 7, 12: 8, 13: 11}
Q_P = {8: 29, 9: 35, 10: 39, 11: 76, 12: 101, 13: 146}
U_P = {8: Fr(707, 1000), 9: Fr(2929, 5000), 10: Fr(4857, 10000), 11: Fr(20111, 50000),
       12: Fr(333511, 10**6), 13: Fr(2765404267, 10**10)}

corners = 0
for P, n in N_P.items():
    lo, hi = enclosure(P, n)
    assert hi < U_P[P]
    lo1, _ = enclosure(P, n - 1)
    assert lo1 > U_P[P]
    corners += 1
lo7, hi7 = enclosure(7, 1)
assert hi7 < 1
assert enclosure(7, 0)[0] > 1
corners += 1

junctions = 0
for P in range(8, 14):
    for Q in range(Q_P[P], P*(N_P[P] + 1)):
        if gcd(P, Q) != 1:
            continue
        m = continuant_label(P, Q)
        f = fib(2*(P + Q - 2) + 1)
        assert U_P[P].numerator*f - U_P[P].denominator*m > 0
        junctions += 1
assert junctions == 38
assert U_P[8].numerator*fib(71) - U_P[8].denominator*continuant_label(8, 29) == 25054682255203

row_seven = 0
for Q in range(8, 14):
    assert continuant_label(7, Q) < fib(2*Q + 11)
    row_seven += 1
assert (continuant_label(7, 8), fib(27)) == (195025, 196418)

print(json.dumps({"status": "PASS", "corner_comparisons": corners, "junction_comparisons": junctions,
                  "row_seven_comparisons": row_seven, "sqrt5_digits": DIGITS}, sort_keys=True))
