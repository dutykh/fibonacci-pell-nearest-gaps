# Independent audit: complete quartic exclusion of the plus $j=1$ branch

- **Date:** 2026-08-23
- **Author/agent:** agent j1 (`/root/j1_quartic_attack`)
- **Objective:** Classify the exact C0068 plus-branch boundary equation
  $$
  P_M=3F_q^2+2(-1)^q
  $$
  for positive odd $M$, and in particular exclude it for $q\geqslant7$.
- **Status:** complete computer-assisted branch theorem, independently
  reproduced. The only pairs with $q\geqslant0$ and positive odd $M$ are
  $(q,M)=(1,1),(2,3),(4,5)$. Thus the plus normalized-gap exponent $j=1$
  never occurs in the unresolved range $q\geqslant7$. This does not exclude
  the other normalized exponents or solve Markoff--Frobenius uniqueness.
- **Inputs:** C0068, the exact boundary identity derived in
  `2026-08-23-nearest-even-gap-recurrence-agent-neg.md`, the negative Pell
  identity, and Magma V2.29-9's documented complete integral-quartic
  routine.

## Exact reduction

Write

$$
y=F_q,
\qquad
\varepsilon=(-1)^q.
\tag{1}
$$

If C0068's plus gap has normalized unit exponent $j=1$, its exact
reconstruction has $M=n_+-1$ and the radical coefficient identity is

$$
P_M=F_{q-1}^2+F_{q+1}^2
=3F_q^2+2(-1)^q.
\tag{2}
$$

The index $M$ is positive and odd. Hence

$$
C_M^2-2P_M^2=-1.
\tag{3}
$$

Substitution of (1)--(2) into (3) gives the necessary quartic equation

$$
\boxed{
C_M^2=18y^4+24\varepsilon y^2+7.
}
\tag{4}
$$

In fact (4) is equivalent to the assertion that
$3y^2+2\varepsilon$ is the radical coordinate of some positive odd Pell
unit: the positive solution of $C^2-2P^2=-1$ lies on the unique positive
odd orbit. Only necessity is needed below.

The constant coefficient $7$ in (4) is not a square, so the one-argument
form of Magma's `IntegralQuarticPoints` cannot be applied directly. Use the
integral translation

$$
z=y-1,
\qquad
y=z+1.
\tag{5}
$$

For $\varepsilon=1$, equation (4) becomes

$$
\boxed{
C^2
=18z^4+72z^3+132z^2+120z+49.
}
\tag{6}
$$

For $\varepsilon=-1$, it becomes

$$
\boxed{
C^2
=18z^4+72z^3+84z^2+24z+1.
}
\tag{7}
$$

Thus both translated constant coefficients are squares. The translation
(5) is a bijection on the integral abscissae, so it neither loses nor
invents an integral solution. Both quartics are nonsingular; their
polynomial discriminants are

$$
10{,}450{,}944\ne0.
\tag{8}
$$

## Complete integral-point lists

The official Magma handbook specifies that, when the constant coefficient
is a square, `IntegralQuarticPoints([a,b,c,d,e])` returns all integral
points modulo ordinate negation on
$C^2=az^4+bz^3+cz^2+dz+e$. See the
[integral-quartic specification](https://magma.maths.usyd.edu.au/magma/handbook/text/1567).

I submitted
`../scripts/check_j1_plus_branch_quartic_independent_audit_agent_j1.m`
to the University of Sydney's public Magma calculator. This was a run
independent of the curator's first calculation. It used the distinct fixed
seeds $42424243$ and $98765431$ for the two quartics. Magma V2.29-9
completed in $0.870$ seconds with $85.16$ MB of memory. After normalizing
the ordinate sign, it returned

```text
Qplus_normalized [
    [ -4, 41 ],
    [ -2, 7 ],
    [ 0, 7 ],
    [ 2, 41 ]
]
Qminus_normalized [
    [ -2, 1 ],
    [ 0, 1 ]
]
```

The full transcript of the certificate-relevant output was

```text
VERSION V2.29-9
PLUS_SEED 42424243 0 3571
Qplus_normalized [
    [ -4, 41 ],
    [ -2, 7 ],
    [ 0, 7 ],
    [ 2, 41 ]
]
MINUS_SEED 98765431 0 1812
Qminus_normalized [
    [ -2, 1 ],
    [ 0, 1 ]
]

Using model [ 0, 0, 0, -696, -7040 ]
Torsion Subgroup = Z/2
Analytic rank = 2
The 2-Selmer group has rank 3
New point of infinite order (x = 272)
New point of infinite order (x = -712/49)
After 2-descent:
    2 <= Rank(E) <= 2
    Sha(E)[2] is trivial
(Searched up to height 100 on the 2-coverings.)

Using model [ 0, 0, 0, -696, 7040 ]
Torsion Subgroup = Z/2
Analytic rank = 1
     ==> Rank(E) = 1
The 2-Selmer group has rank 2
New point of infinite order (x = 17)
After 2-descent:
    1 <= Rank(E) <= 1
    Sha(E)[2] is trivial
(Searched up to height 100 on the 2-coverings.)

DISCRIMINANTS 10450944 10450944
ASSOCIATED_PLUS Elliptic Curve defined by y^2 = x^3 - 696*x - 7040 over
Rational Field
PLUS_RANK_INFO [ 2, 2 ]
PLUS_MW_STATUS true true
ASSOCIATED_MINUS Elliptic Curve defined by y^2 = x^3 - 696*x + 7040 over
Rational Field
MINUS_RANK_INFO [ 1, 1 ]
MINUS_MW_STATUS true true
PASS: both complete shifted-quartic lists and Mordell-Weil statuses were
independently reproduced.
```

The random-state step counts $3571$ and $1812$ confirm that both separately
seeded integral-point calls executed. The exact $2$-descent rank intervals
are $[2,2]$ and $[1,1]$, respectively. Calls to `MordellWeilGroup` returned
both status flags `true` for both associated curves. The
[Mordell--Weil interface](https://magma.maths.usyd.edu.au/magma/handbook/text/1570)
states that these flags certify, respectively, the rank and the full
Mordell--Weil group. Thus the usual generic subgroup-completeness caveat is
absent in these two calculations. The analytic-rank lines do not carry the
proof; the exact descent bounds and full-group flags do.

It follows from (5) and the complete lists that all integral solutions of
(4), modulo the sign of $C$, have

$$
\begin{aligned}
\varepsilon=1:&\quad y\in\{-3,-1,1,3\},\\
\varepsilon=-1:&\quad y\in\{-1,1\}.
\end{aligned}
\tag{9}
$$

## Fibonacci and Pell filtering

Now impose $y=F_q\geqslant0$ and $\varepsilon=(-1)^q$. The Fibonacci
sequence is strictly increasing from index $2$ onward. Equation (9) leaves

$$
q=1,
\qquad
q=2,
\qquad
q=4.
\tag{10}
$$

The corresponding right sides of (2) are

$$
1=P_1,
\qquad
5=P_3,
\qquad
29=P_5.
\tag{11}
$$

The positive Pell coordinates $P_M$ are strictly increasing, so (11)
recovers the unique odd indices. Therefore

$$
\boxed{
P_M=3F_q^2+2(-1)^q,\quad M>0\text{ odd}
\quad\Longleftrightarrow\quad
(q,M)\in\{(1,1),(2,3),(4,5)\}.
}
\tag{12}
$$

In particular, (2) has no solution for $q\geqslant7$, whether or not the
additional condition $3\nmid q$ is imposed. This completely excludes the
C0068 plus branch with normalized gap exponent $j=1$ in the unresolved
range.

## Independent exact arithmetic checker and scope

The standard-library checker
`../scripts/check_j1_plus_branch_quartic_independent_audit_agent_j1.py`
expands the shifted quartics, verifies every returned integral point,
undoes the shift, reconstructs the Pell square, filters the Fibonacci
indices, and reproduces (10)--(12). It also performs a separate bounded
falsification scan for $|y|\leqslant10{,}000$. That bounded scan is only a
regression check. Completeness for all integral $y$ comes from the documented
Magma integral-point calculation and the proved full Mordell--Weil groups,
not from the scan.

This result closes exactly one normalized-unit boundary, $j=1$ on the plus
branch. It says nothing by itself about the other odd values of $j$, the
minus branch, or the remaining half-gap scalar balance. It is therefore a
genuine unbounded branch exclusion but neither a proof nor a counterexample
to Markoff--Frobenius uniqueness.
