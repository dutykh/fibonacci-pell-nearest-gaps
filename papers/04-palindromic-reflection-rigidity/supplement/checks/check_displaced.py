#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Fixed exact identities for Section 5; no word or parameter enumeration.

Requires Python 3 and SymPy. Run without arguments, with assertions enabled.
The rational Gram controls check arithmetic only: no word realization is
claimed. The universal inequalities and central-word hypotheses are proved
in the manuscript; this program is a separate algebraic verification.
"""

import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("run with assertions enabled and without arguments")

from fractions import Fraction as F
from hashlib import sha256
from math import gcd
import json

import sympy as S


checked = []


def zero(name, expression):
    entries = list(expression) if isinstance(expression, S.MatrixBase) else [expression]
    assert all(S.cancel(entry) == 0 for entry in entries), name
    checked.append(name)


def truth(name, statement):
    assert statement, name
    checked.append(name)


x, y, M, rx, ry, delta = S.symbols("x y M rx ry delta", nonzero=True)
G = S.Matrix([[x, y], [y, (y*y + 1)/x]])
P = S.Matrix([[2, 1], [1, 0]])
A = S.Matrix([[2, 1], [1, 1]])
B = S.Matrix([[5, 2], [2, 1]])
I = S.eye(2)
e = S.Matrix([2, 1])
f = S.Matrix([1, -2])
omega = G*e
M0 = S.cancel((e.T*omega)[0])
rho0 = omega[0]
r0 = S.cancel(rho0/M0)
a0 = (G.inv()*f).applyfunc(S.cancel)
R0 = ((2*e*e.T*G - M0*I)/M0).applyfunc(S.cancel)

zero("generic Gram determinant", G.det() - 1)
zero("primitive displacement column", a0 - S.Matrix([M0-2*rho0, -rho0]))
zero("orthogonal covectors", (omega.T*a0)[0])
zero("terminal displacement weight", (f.T*a0)[0] - M0)
zero("ancestor involution", R0*R0 - I)
zero("ancestor terminal", R0*e-e)
zero("ancestor adjoint", R0.T*G-G*R0)
zero("odd covector", f.T*R0+f.T)

tau = 2*r0 + delta
J = S.Matrix([[1, tau], [0, -1]])
Rt = (P*J*P.inv()).applyfunc(S.cancel)
Lt = (G.inv()*Rt.T*G).applyfunc(S.cancel)
zero("right rank one", Rt-R0-delta*e*f.T)
zero("left rank one", Lt-R0-delta*a0*omega.T)
zero("right terminal fixed", Rt*e-e)
zero("left terminal moving", Lt*e-e-delta*M0*a0)
zero("weighted scale fixed", omega.T*Lt-omega.T)
zero("parabolic product", Lt*R0-I-delta*a0*omega.T)
zero("nilpotent defect", (a0*omega.T)**2)


def gram(label, root):
    return S.Matrix([[label, root], [root, (root*root+1)/label]])


HX, HY = gram(M, rx), gram(M, ry)
endpoint_delta = (rx+ry)/M-2*r0
Je = J.subs(delta, endpoint_delta).applyfunc(S.cancel)
Re = Rt.subs(delta, endpoint_delta).applyfunc(S.cancel)
Le = Lt.subs(delta, endpoint_delta).applyfunc(S.cancel)
NX = (G.inv()*P.inv()*HX*P.inv()).applyfunc(S.cancel)
NY = (G.inv()*P.inv()*HY*P.inv()).applyfunc(S.cancel)
zero("full Gram congruence", Je.T*HX*Je-HY)
zero("full two-sided word equation", NY-Le*NX*Re)
zero("actual common weighted scale", (omega.T*NX*e)[0]-M)
zero("actual moving column", NY*e-R0*NX*e-endpoint_delta*M*a0)

v1, v2 = S.symbols("v1 v2")
v = S.Matrix([v1, v2])


def eta_column(column):
    numerator = S.cancel(M0*(f.T*column)[0])
    denominator = S.cancel((omega.T*column)[0])
    return S.cancel(numerator/denominator)


zero("homogeneous cusp formula", eta_column(v)-(v1-2*v2)/(v2+r0*(v1-2*v2)))
zero("exact moving terminal coordinate", eta_column(Lt*v)+eta_column(v)-delta*M0*M0)
zero("actual cusp sum", eta_column(NX*e)+eta_column(NY*e)-endpoint_delta*M0*M0)

expected_trace = M/M0+M0/M+M0/M*(rx-M*rho0/M0)**2
zero("relative trace completed square", S.trace(NX)-expected_trace)
zero("relative trace difference", S.trace(NY)-S.trace(NX)-M0*(ry-rx)*endpoint_delta)

frame = P*G
u, vv, s, t = frame[0, 0], frame[0, 1], frame[1, 0], frame[1, 1]
zero("prefix determinant", u*t-vv*s+1)
zero("prefix special column", u-2*s-t)
zero("prefix ancestor coordinates", M0-2*u-vv)
zero("Farey interval width", s/u-t/vv-1/(u*vv))

c0, Sigma, alpha, beta = S.symbols("c0 Sigma alpha beta", nonzero=True)
ac = (vv*Sigma-2*M*t)/c0
bc = (2*M*s-u*Sigma)/c0
zero("primitive direction first coordinate", ac*u+bc*vv-2*M/c0)
zero("primitive direction second coordinate", ac*s+bc*t-Sigma/c0)
zero("integral midpoint defect", M0*Sigma-2*rho0*M-c0*(ac-2*bc))
Q = alpha*u+beta*vv
Sigma_reduced = alpha*s+beta*t
rr = alpha/beta
xi = M0*M0*(2*Sigma_reduced/Q-2*r0)
zero("primitive direction nonlinear mean", xi-2*(rr-2)/(1+r0*(rr-2)))
kap, eps = S.symbols("kap eps", nonzero=True)
zero("all parity trace defect", (c0*(alpha-2*beta))*(kap*Q/eps)/(c0*Q/2)-2*kap*(alpha-2*beta)/eps)

controls = []
for label, rlow, rhigh in ((65, 8, 18), (2725, 1057, 1232),
                          (130, 47, 57), (610, 133, 233)):
    total, difference = rlow+rhigh, rhigh-rlow
    common = gcd(total, 2*label)
    denominator = 2*label//common
    parity = gcd(denominator, 2)
    kappa = F(parity*difference, denominator)
    # Use the empty-prefix unimodular frame, without claiming literal words.
    aa, bb = total//common, (2*label-2*total)//common
    truth("root control %d" % label,
          0 < rlow < rhigh < F(label, 2)
          and (rlow*rlow+1) % label == 0
          and (rhigh*rhigh+1) % label == 0
          and kappa.denominator == 1 and kappa > 0
          and gcd(aa, bb) == 1 and aa > 0 and bb > 0
          and 2*aa+bb == denominator
          and 5*total-4*label == common*(aa-2*bb)
          and F(difference, label) == F(2*kappa, parity*common))
    controls.append([label, rlow, rhigh, common, denominator, parity, int(kappa)])
truth("both label and denominator parities", {(r[0] % 2, r[4] % 2) for r in controls} == {(0, 0), (0, 1), (1, 0), (1, 1)})

for name, short in (("empty", I), ("a", A), ("b", B)):
    difference = A*A-short
    truth("short-core domination " + name,
          all(entry >= 0 for entry in difference) and any(entry > 0 for entry in difference))
truth("strict short-core ordering", all(entry >= 0 for entry in A-I)
      and any(entry > 0 for entry in A-I)
      and all(entry >= 0 for entry in B-A) and any(entry > 0 for entry in B-A))

us, vs, astar, bstar, eu, ev = S.symbols("us vs astar bstar eu ev")
L = us+vs-3
dstar = (astar-bstar)/L
rr0 = (us-vs+dstar)/2
yx, yy = -us+eu+rr0, vs+ev+rr0
zero("left centered parent cylinder", yx+(L+3)/2-dstar/2-eu)
zero("right centered parent cylinder", yy-(L+3)/2-dstar/2-ev)
zero("parent-center cancellation", 1/yx+1/yy-(dstar+eu+ev)/(yx*yy))
truth("cylinder lower a coefficient", F(1, 3)-F(5, 12) == -F(1, 12))
truth("cylinder lower b coefficient", -F(2, 5)-F(40, 111) == -F(422, 555))
truth("cylinder upper a coefficient", F(2, 5)+F(12, 25) == F(22, 25))
truth("cylinder upper b coefficient", -F(1, 3)+F(1, 2) == F(1, 6))
truth("absolute numerator coefficient", max(F(1, 12), F(422, 555), F(22, 25), F(1, 6)) == F(22, 25))
truth("cusp denominator bound", F(22, 25)/F(24, 5) == F(11, 60))
truth("global cusp bound", F(11, 60)*F(5, 4) == F(11, 48))
truth("positive inversion denominator", 2-F(1, 2)*F(11, 48) == F(181, 96))
truth("parent-sensitive cone", F(11, 60)/F(181, 96) == F(88, 905) < F(1, 10))
truth("gamma bound", 1-F(1, 2)*F(1, 8) == F(15, 16))
truth("uniform quartic coefficient", F(3, 1)/(10*F(15, 16)) == F(8, 25))

p, q = S.symbols("p q", positive=True)
zero("parent imbalance conversion", (p/q+q/p)/(p*q)-(p**-2+q**-2))
z = S.symbols("z", positive=True)
zero("prefix derivative", S.diff(z/(2*z+1)**2, z)-(1-2*z)/(2*z+1)**3)
zero("prefix radical endpoint", ((1+S.sqrt(2))/(3+2*S.sqrt(2))**2)-(5*S.sqrt(2)-7))
zero("prefix reciprocal constant", (5*S.sqrt(2)-7)*(7+5*S.sqrt(2))-1)
phi, upper = (1+S.sqrt(5))/2, 1+S.sqrt(2)
map_a = lambda value: (2*value+1)/(value+1)
map_b = lambda value: (5*value+2)/(2*value+1)
zero("left cone fixed point", S.radsimp(map_a(phi)-phi))
zero("right cone fixed point", S.radsimp(map_b(upper)-upper))
truth("invariant cone inner endpoints", S.simplify(map_a(upper) < upper) is S.true
      and S.simplify(map_b(phi) > phi) is S.true and phi < 2 < upper)

payload = {"checks": checked, "arithmetic_controls": controls,
           "scope": "fixed exact algebra, not a word or collision census"}
digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
print(json.dumps({"status": "PASS", "checks": len(checked),
                  "arithmetic_controls": len(controls), "reference_sha256": digest}, sort_keys=True))
