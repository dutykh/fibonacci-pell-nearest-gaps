#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Independent fixed algebra checks for the all-defect parity theorem.

No gap lattice is enumerated. Universal identities and four prescribed
integer paths are checked; the all-parameter proof is in the article.
Run under the documented external 128 MiB / 30 second limits.
"""

import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("assertions required; this checker accepts no arguments")

import hashlib
import json
import sympy as s

root = s.sqrt(5)
phi = (1 + root) / 2
r = (3 - root) / 2
q = -r
eye = s.eye(2)
A = s.Matrix([[1, 0], [1, 1]])
B = s.Matrix([[0, 1], [1, 1]])
T = s.Matrix([[1, phi], [phi, -1]])
z = s.Matrix([1, 2])


def zero(expr):
    return s.simplify(s.expand(expr)) == 0


def matrix_equal(left, right):
    assert left.shape == right.shape
    assert all(zero(entry) for entry in left-right)


matrix_equal(T*T, (1+phi**2)*eye)
matrix_equal(T*B, phi*s.diag(1, q)*T)
matrix_equal(root/phi*T*A, s.Matrix([[2, 1], [q, -2*q]])*T)
matrix_equal(T*z/phi, s.Matrix([phi**2, 2-root]))
matrix_equal(T*s.Matrix([phi**-2, 1]), s.Matrix([2, q]))
matrix_equal(T*s.Matrix([phi**-1, phi**-3]), s.Matrix([1, -2*q]))

u = s.Matrix([phi**-2, 1])
v = s.Matrix([phi**-1, phi**-3])
assert all(s.simplify(a).is_nonnegative for a in v-phi**-3*u)
assert all(s.simplify(a).is_nonnegative for a in phi*u-v)
rho = s.symbols("rho", positive=True)
h = (rho-phi)/(1+phi*rho)
assert zero(h.subs(rho, 2)-phi**-5)
assert zero(s.diff(h, rho)-(1+phi**2)/(1+phi*rho)**2)
assert zero(1/phi-h-(1+phi**-2)/(rho+phi**-1))
assert zero(phi**8-r**-4)
assert zero(phi**-10-r**5)
assert s.simplify(1-r**4-r**9).is_positive
assert s.simplify(r**5-r**7).is_positive

C, a, b, d, t = s.symbols("C a b d t", positive=True)
f = (C+2*a*t+b*t*t)/(1+d*t*t)
assert zero(s.diff(f, t)*(1+d*t*t)**2-2*(a+b*t-C*d*t-a*d*t*t))
X, Y, zz = s.symbols("X Y zz", real=True)
f = (X*X+Y*Y*zz)/(1+d*zz)
assert zero(s.diff(f, zz)*(1+d*zz)**2-(Y*Y-d*X*X))

controls = ((4, 0, 3), (1, 3, 0, 2), (8, 2, 1, 0, 2), (0, 0, 0))
register = []
for gaps in controls:
    k = len(gaps)-1
    physical = z
    channel = s.Matrix([phi**2, 2-root])
    for j, g in enumerate(gaps):
        physical = B**g*physical
        channel = s.diag(1, q**g)*channel
        if j < k:
            physical = A*physical
            channel = s.Matrix([[2, 1], [q, -2*q]])*channel
        channel = channel.applyfunc(s.simplify)
    cost = sum(gaps)+k+1
    norm = int((physical.T*physical)[0])
    P = s.simplify((channel.T*channel)[0])
    K = 5**k*phi**4*(1-q**(2*cost+3))
    assert zero(P*s.fibonacci(2*cost+3)-K*norm)
    register.append([list(gaps), cost, norm, int(s.fibonacci(2*cost+3))])

payload = json.dumps(register, separators=(",", ":")).encode()
print(json.dumps({"status": "PASS", "scope": "fixed all-defect algebra schemas",
                  "controls": register,
                  "sha256": hashlib.sha256(payload).hexdigest()}, sort_keys=True))
