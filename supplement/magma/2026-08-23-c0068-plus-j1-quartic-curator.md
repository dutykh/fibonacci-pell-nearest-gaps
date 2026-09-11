# Attempt: Complete plus-branch $j=1$ quartic exclusion

- **Date:** 2026-08-23
- **Author/curator:** Codex `/root`
- **Objective:** classify C0068 plus-branch hits whose normalized nearest-gap
  exponent is $j=1$.
- **Status:** complete computer-assisted classification, independently rerun
  with distinct seeds and exact rank certificates. The only solutions of
  the resulting recurrence equation are the three low pairs displayed
  below. In C0068's admissible domain the coefficient equation retains
  $q=2,4$, but its nearest-gap filter keeps only the known plus hit at $q=4$.
  There is no such hit for $q\geqslant7$.

## Exact reduction

Write

$$
\lambda^M=C_M+P_M\sqrt2,
\qquad
\lambda=1+\sqrt2,
$$

with $M>0$ odd. The shifted-root identity in
../attempts/2026-08-23-nearest-even-gap-recurrence-agent-neg.md shows that a
plus-branch collision with normalized gap exponent $j=1$ must satisfy

$$
P_M
=F_{q-1}^2+F_{q+1}^2
=3F_q^2+2(-1)^q.
\tag{1}
$$

Conversely, this note only classifies equation (1); C0068 supplies the
selector interpretation when the nearest-gap hypotheses also hold.

Put $y=F_q$ and $s=(-1)^q$. Since $M$ is odd,

$$
C_M^2-2P_M^2=-1.
\tag{2}
$$

Substitution of (1) into (2) gives the necessary quartic

$$
\boxed{
C_M^2=18y^4+24sy^2+7.
}
\tag{3}
$$

Set $y=z+1$. For $s=1$ and $s=-1$, respectively, equation (3) becomes

$$
\begin{aligned}
C_M^2
&=18z^4+72z^3+132z^2+120z+49,\\
C_M^2
&=18z^4+72z^3+84z^2+24z+1.
\end{aligned}
\tag{4}
$$

Both quartics are nonsingular; their polynomial discriminant is
$10{,}450{,}944$. Their constant terms are the squares $7^2$ and $1^2$, so
Magma's documented `IntegralQuarticPoints` routine applies and returns every
integral point, modulo ordinate negation.

## Complete official Magma run

The exact source is
../scripts/check_c0068_plus_j1_quartic_curator.m. It was submitted to the
University of Sydney public calculator at
`https://magma.maths.usyd.edu.au/xml/calculator.xml`, with fixed seed
$8{,}675{,}309$. Magma V2.29-9 completed in $0.660$ seconds and returned

```text
Qeven:
[ -2, -7 ]
[ 0, -7 ]
[ 2, -41 ]
[ -4, 41 ]
Qodd:
[ -2, -1 ]
[ 0, -1 ]
```

The coordinates here are $(z,C_M)$. Undoing $y=z+1$ shows

$$
|y|\in\{1,3\}\quad(s=1),
\qquad
|y|=1\quad(s=-1).
\tag{5}
$$

For nonnegative Fibonacci indices, monotonicity from $q\geqslant2$ and the
parity condition $s=(-1)^q$ leave exactly

$$
\boxed{
(q,M)=(1,1),(2,3),(4,5).
}
\tag{6}
$$

Direct substitution gives $P_1=1$, $P_3=5$, and $P_5=29$, so every pair in
(6) really solves (1). Therefore the admissible $q\geqslant2$, $3\nmid q$
solutions of the coefficient equation are just $(2,3)$ and $(4,5)$. At
$q=2$ the plus nearest-gap window is absent; at $q=4$ one recovers the known
plus gap. Thus no $j=1$ plus-branch collision can occur for
$q\geqslant7$.

## Certificate and scope

The standard-library verifier
../scripts/check_c0068_plus_j1_quartic_curator.py checks the shifts,
nonsingularity, every returned point, the Fibonacci/Pell filtering, and a
bounded recurrence regression. Completeness comes from Magma V2.29-9's
documented elliptic-logarithm routine, not from that bounded loop.

The independent rerun, exact Jacobian-rank audit, and separate verifier are
in
../attempts/2026-08-23-j1-plus-branch-quartic-independent-audit-agent-j1.md,
../scripts/check_j1_plus_branch_quartic_independent_audit_agent_j1.m, and
../scripts/check_j1_plus_branch_quartic_independent_audit_agent_j1.py.

This closes exactly the plus $j=1$ stratum. Other odd normalized exponents
and the entire possible single-branch collision remain open. It is not a
proof or counterexample to Markoff--Frobenius uniqueness.
