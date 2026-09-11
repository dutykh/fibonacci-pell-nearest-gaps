#!/usr/bin/env python3
"""Exact checks for the adjusted full-eight-term logarithmic form.

Mathematical target: reconstruct C0082's eight terms in all four sign
sectors, verify their Galois action, and check the exact regrouping

    lambda^n = alpha_{sigma,h} phi^(2q) + R

conditional on the trace equation. The script also verifies the full-norm
formula proving nonvanishing and the rational arithmetic in the proposed
q < 10^32 reduction.

This is a standard-library checker. It performs no exponent search and makes
no bounded no-hit claim.

Deterministic command from math-sandbox/:
    python3 -B scripts/check_c0082_full8_adjusted_form_agent_c.py
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


Q = Fraction


@dataclass(frozen=True)
class MQ:
    """Element a + b*sqrt(2) + c*sqrt(5) + d*sqrt(10)."""

    a: Q = Q(0)
    b: Q = Q(0)
    c: Q = Q(0)
    d: Q = Q(0)

    def __add__(self, other: object) -> "MQ":
        if not isinstance(other, MQ):
            other = MQ(Q(other))
        return MQ(
            self.a + other.a,
            self.b + other.b,
            self.c + other.c,
            self.d + other.d,
        )

    __radd__ = __add__

    def __neg__(self) -> "MQ":
        return MQ(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, other: object) -> "MQ":
        return self + (-as_mq(other))

    def __rsub__(self, other: object) -> "MQ":
        return as_mq(other) - self

    def __mul__(self, other: object) -> "MQ":
        if not isinstance(other, MQ):
            other = MQ(Q(other))
        a, b, c, d = self.a, self.b, self.c, self.d
        e, f, g, h = other.a, other.b, other.c, other.d
        return MQ(
            a * e + 2 * b * f + 5 * c * g + 10 * d * h,
            a * f + b * e + 5 * c * h + 5 * d * g,
            a * g + c * e + 2 * b * h + 2 * d * f,
            a * h + d * e + b * g + c * f,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "MQ":
        if isinstance(other, MQ):
            raise TypeError("only rational scalar division is needed")
        return self * Q(1, other)

    def __pow__(self, exponent: int) -> "MQ":
        if exponent < 0:
            raise ValueError("use the explicit unit inverse for negative powers")
        answer = MQ(Q(1))
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power >>= 1
        return answer

    def auto(self, root_two_sign: int, root_five_sign: int) -> "MQ":
        return MQ(
            self.a,
            root_two_sign * self.b,
            root_five_sign * self.c,
            root_two_sign * root_five_sign * self.d,
        )


def as_mq(value: object) -> MQ:
    return value if isinstance(value, MQ) else MQ(Q(value))


ONE = MQ(Q(1))
R2 = MQ(b=Q(1))
R5 = MQ(c=Q(1))
PHI = (ONE + R5) / 2
PHI_INV = PHI - 1
LAMBDA = ONE + R2
LAMBDA_INV = R2 - 1


def unit_power(unit: MQ, inverse: MQ, exponent: int) -> MQ:
    if exponent >= 0:
        return unit**exponent
    return inverse ** (-exponent)


def phi_power(exponent: int) -> MQ:
    return unit_power(PHI, PHI_INV, exponent)


def lambda_power(exponent: int) -> MQ:
    return unit_power(LAMBDA, LAMBDA_INV, exponent)


def fibonacci(index: int) -> int:
    a, b = 0, 1
    for _ in range(index):
        a, b = b, a + b
    return a


def pell_coordinates(index: int) -> tuple[int, int]:
    """Return C_index, P_index from lambda^index=C+P*sqrt(2)."""
    c_value, p_value = 1, 0
    for _ in range(index):
        c_value, p_value = c_value + 2 * p_value, c_value + p_value
    return c_value, p_value


def field_norm(value: MQ) -> Q:
    product = ONE
    for root_two_sign in (1, -1):
        for root_five_sign in (1, -1):
            product *= value.auto(root_two_sign, root_five_sign)
    assert product.b == product.c == product.d == 0
    return product.a


def eight_terms(q: int, n: int, h: int, sigma: int, epsilon: int) -> list[MQ]:
    """Return C0082's T_0,...,T_7 for M=n+h."""
    u_plus = sigma + R2 * PHI
    u_minus = -sigma + R2 * PHI_INV
    v_plus = sigma - R2 * PHI
    v_minus = -sigma - R2 * PHI_INV
    m_big = n + h
    return [
        5 * lambda_power(m_big),
        5 * lambda_power(-m_big),
        -(u_plus**2) * lambda_power(h) * phi_power(2 * q),
        -(u_minus**2) * lambda_power(h) * phi_power(-2 * q),
        -2 * epsilon * (1 - sigma * R2) * lambda_power(h),
        -(v_plus**2) * lambda_power(-h) * phi_power(2 * q),
        -(v_minus**2) * lambda_power(-h) * phi_power(-2 * q),
        -2 * epsilon * (1 + sigma * R2) * lambda_power(-h),
    ]


def adjusted_data(
    q: int, n: int, h: int, sigma: int, epsilon: int
) -> tuple[MQ, MQ, MQ, MQ]:
    """Return D, Dbar, alpha, and the bounded remainder R."""
    f_value = fibonacci(q)
    x_value = fibonacci(q + 1)
    d_value = sigma * f_value + x_value * R2
    d_bar = sigma * f_value - x_value * R2

    u_value = sigma + R2 * PHI
    v_value = -sigma + R2 * PHI_INV
    w_value = sigma - R2 * PHI
    z_value = -sigma - R2 * PHI_INV
    alpha = (u_value**2 + (w_value**2) * lambda_power(-2 * h)) / 5
    remainder = (
        Q(2 * epsilon, 5)
        * ((1 - sigma * R2) + (1 + sigma * R2) * lambda_power(-2 * h))
        + Q(1, 5)
        * (v_value**2 + (z_value**2) * lambda_power(-2 * h))
        * phi_power(-2 * q)
        - lambda_power(-n - 2 * h)
    )
    return d_value, d_bar, alpha, remainder


Interval = tuple[Q, Q]


def interval_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def interval_multiply(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def interval_scale(interval: Interval, scalar: Q) -> Interval:
    if scalar >= 0:
        return scalar * interval[0], scalar * interval[1]
    return scalar * interval[1], scalar * interval[0]


def square_root_interval(radicand: int, digits: int = 70) -> Interval:
    denominator = 10**digits
    numerator = isqrt(radicand * denominator * denominator)
    return Q(numerator, denominator), Q(numerator + 1, denominator)


SQRT_TWO_INTERVAL = square_root_interval(2)
SQRT_FIVE_INTERVAL = square_root_interval(5)
SQRT_TEN_INTERVAL = interval_multiply(SQRT_TWO_INTERVAL, SQRT_FIVE_INTERVAL)


def real_interval(value: MQ) -> Interval:
    answer = (value.a, value.a)
    for coefficient, basis_interval in (
        (value.b, SQRT_TWO_INTERVAL),
        (value.c, SQRT_FIVE_INTERVAL),
        (value.d, SQRT_TEN_INTERVAL),
    ):
        answer = interval_add(answer, interval_scale(basis_interval, coefficient))
    return answer


def abs_upper(interval: Interval) -> Q:
    return max(abs(interval[0]), abs(interval[1]))


def check_galois_and_support() -> None:
    p_two = (1, 0, 5, 6, 7, 2, 3, 4)
    p_five = (0, 1, 3, 2, 4, 6, 5, 7)
    supports = (
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 2),
        (0, 1, -2),
        (0, 1, 0),
        (0, -1, 2),
        (0, -1, -2),
        (0, -1, 0),
    )

    dominant_pairs = {
        "identity": (0, 2),
        "tau_5": (0, 3),
        "tau_2": (1, 5),
        "tau_2_tau_5": (1, 6),
    }
    differences = []
    for pair in dominant_pairs.values():
        differences.append(
            tuple(supports[pair[0]][j] - supports[pair[1]][j] for j in range(3))
        )
    assert differences == [
        (1, -1, -2),
        (1, -1, 2),
        (-1, 1, -2),
        (-1, 1, 2),
    ]
    assert differences[2] == tuple(-entry for entry in differences[1])
    assert differences[3] == tuple(-entry for entry in differences[0])

    for q in (7, 8, 13):
        for n in (2, 8, 14):
            for h in (3, 5):
                for sigma in (1, -1):
                    for epsilon in (1, -1):
                        terms = eight_terms(q, n, h, sigma, epsilon)
                        for index, term in enumerate(terms):
                            assert term.auto(-1, 1) == -terms[p_two[index]]
                            assert term.auto(1, -1) == terms[p_five[index]]

                        signs = []
                        for term in terms:
                            interval = real_interval(term)
                            assert interval[0] > 0 or interval[1] < 0
                            signs.append(1 if interval[0] > 0 else -1)
                        kappa = sigma * epsilon
                        assert tuple(signs) == (1, 1, -1, -1, kappa, -1, -1, -kappa)


def check_exact_regrouping() -> None:
    for q in (7, 8, 13, 20):
        for n in (2, 8, 14):
            for h in (3, 5, 9):
                for sigma in (1, -1):
                    epsilon = -1 if q % 2 else 1
                    d_value, d_bar, alpha, remainder = adjusted_data(
                        q, n, h, sigma, epsilon
                    )
                    terms = eight_terms(q, n, h, sigma, epsilon)
                    m_big = n + h
                    trace_residual = (
                        lambda_power(m_big)
                        + lambda_power(-m_big)
                        - d_value**2 * lambda_power(h)
                        - d_bar**2 * lambda_power(-h)
                    )
                    assert sum(terms, MQ()) == 5 * trace_residual

                    regrouped = (
                        d_value**2
                        + d_bar**2 * lambda_power(-2 * h)
                        - lambda_power(-n - 2 * h)
                    )
                    assert alpha * phi_power(2 * q) + remainder == regrouped
                    assert lambda_power(n) - regrouped == (
                        trace_residual * lambda_power(-h)
                    )

                    alpha_interval = real_interval(alpha)
                    remainder_interval = real_interval(remainder)
                    assert alpha_interval[0] > Q(1, 5)
                    assert alpha_interval[1] < 7
                    assert abs_upper(remainder_interval) < 4

    # phi^14=377*phi+233 gives the exact q>=7 smallness threshold.
    assert PHI**14 == 377 * PHI + 233
    assert Q(40, 1597) < Q(1, 39)
    assert Q(20 * 39, 38) < 21


def check_norm_nonvanishing() -> None:
    for sigma in (1, -1):
        for h in range(1, 32, 2):
            c_value, p_value = pell_coordinates(h)
            assert c_value * c_value - 2 * p_value * p_value == -1
            u_value = sigma + R2 * PHI
            w_value = sigma - R2 * PHI
            alpha = (u_value**2 + (w_value**2) * lambda_power(-2 * h)) / 5
            norm_factor = 4 + 3 * p_value * p_value - 2 * sigma * c_value * p_value
            assert field_norm(alpha) == Q(64, 625) * norm_factor * norm_factor
            assert field_norm(alpha) != 1


def check_effective_cutoff() -> None:
    # The already audited Matveev prefactor from C0085, excluding A_3.
    matveev_prefix = 9_275_000_000_000
    q_zero = 10**32

    # At q_zero, H=1+log(2q)<79 follows from log(2)<1 and log(10)<12/5.
    h_upper = 79
    rhs_at_cutoff = Q(5, 4) * (
        matveev_prefix * (960_000_000_000_000 * h_upper + 44) * h_upper + 4
    )
    assert rhs_at_cutoff == 69_462_330_000_000_040_299_875_000_000_005
    assert rhs_at_cutoff < q_zero

    # Once q<10^32, C0085 gives this explicit h-box.
    assert 240_000_000_000_000 * 79 == 18_960_000_000_000_000


def main() -> None:
    check_galois_and_support()
    check_exact_regrouping()
    check_norm_nonvanishing()
    check_effective_cutoff()
    print(
        "PASS: C0082 full-eight adjusted form; four sign sectors and Galois "
        "dominant-pair supports reconstructed; exact alpha_h regrouping and "
        "norm nonvanishing verified; rational Matveev cutoff gives q<10^32 "
        "conditional on C0085 and the archived Matveev specialization."
    )


if __name__ == "__main__":
    main()
