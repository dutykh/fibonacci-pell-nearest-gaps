#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Independent fixed algebra schemas for Sections 2 and 4; no path scan."""

import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("assertions required; no arguments accepted")

import hashlib
import json
import sympy as s

x, y, a, b, d, u, t, v, eps = s.symbols("x y a b d u t v eps")
phi = (1+s.sqrt(5))/2
A = s.Matrix([[2, 1], [1, 1]])
B = s.Matrix([[5, 2], [2, 1]])
G = s.Matrix([[a, b], [b, d]])
M = A*G*B
assert s.expand(M[0, 1]-(4*a+4*b+d)) == 0
assert s.expand(s.trace(M)-3*M[0, 1]) == 0

norm = x*x+y*y
perron = x+phi*y
assert s.simplify(y*y+(x+y)**2-phi**-2*norm-phi**-1*perron**2) == 0
assert s.simplify(y+phi*(x+y)-phi*perron) == 0

# a=F_(c-2), b=F_(c-1). No Cassini relation is needed here.
early = s.Matrix([2*a+b, a+3*b])
last = s.Matrix([a+b, 2*a+3*b])
assert early-last == s.Matrix([a, -a])
assert s.expand((last.T*last)[0]-(early.T*early)[0]-4*a*b) == 0
assert s.expand((early.T*early)[0]-5*(b*b+(a+b)**2)) == 0

# x=F_j, y=F_(j-1), t=F_r, v=F_(r+1), eps=(-1)^j.
Fc = x*v+y*t
Fc2 = (2*x+y)*v+(x+y)*t
X, Y = Fc+t*x, Fc2-t*y
deficit = s.expand(Fc**2+Fc2**2-X**2-Y**2)
printed = t*(2*(x*y+eps)*v-(x-y)*(x+y)*t)
cassini = x*x-x*y-y*y+eps
assert s.expand(deficit-printed+2*t*v*cassini) == 0
assert s.expand(2*x*y-(x-y)*(x+y)-(2*y*y-(x-y)**2)) == 0

normalized_E = u*u+(1+u)**2
normalized_F = 1+(1+u)**2
assert s.expand(6*normalized_E-5*normalized_F-(7*u*u+2*u-4)) == 0
assert (7*u*u+2*u-4).subs(u, s.Rational(5, 8)) == -s.Rational(1, 64)
assert s.Rational(3, 5) < 1/(1+s.Rational(5, 8))
assert 1/(1+s.Rational(3, 5)) == s.Rational(5, 8)

theta = 3/(s.sqrt(5)*phi)
C7 = s.simplify(phi**2*theta**6)
C6 = s.simplify(phi**2*theta**5)
assert s.simplify(C7-729/(125*phi**4)) == 0
assert s.simplify(C6-243/(25*s.sqrt(5)*phi**3)) == 0
assert s.simplify(C7-s.Rational(5, 6)).is_positive
assert s.simplify(C6-1).is_positive
assert 4373**2-5*1875**2 == 1545004
assert 59**2-5*25**2 == 356

register = {"C6": str(C6), "C7": str(C7), "envelope": str(printed),
            "scope": "universal matrix, Fibonacci and radical identities"}
encoded = json.dumps(register, sort_keys=True, separators=(",", ":")).encode()
print(json.dumps({"status": "PASS", "sha256": hashlib.sha256(encoded).hexdigest(),
                  **register}, sort_keys=True))
