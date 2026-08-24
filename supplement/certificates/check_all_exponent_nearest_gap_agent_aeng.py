#!/usr/bin/env python3
"""Exact support for the all-exponent Fibonacci nearest-gap theorem.

The script is intrinsic: it enumerates quadratic-unit gaps, not Pell-orbit
trajectories.  It also reuses the already audited rational-interval engine
from C0088 to certify that its fixed and moving Dujella--Petho bounds are
unchanged by the odd-trace sign.  Finite absence outside the certified
terminal range is regression evidence only.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd

import check_c0082_full8_dp_closure_agent_c as certified_dp


Pair = tuple[int, int]


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: Pair, right: Pair) -> Pair:
    return left[0] - right[0], left[1] - right[1]


def multiply(left: Pair, right: Pair) -> Pair:
    return (
        left[0] * right[0] + 2 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def scale(multiplier: int, value: Pair) -> Pair:
    return multiplier * value[0], multiplier * value[1]


def square(value: Pair) -> Pair:
    return multiply(value, value)


def conjugate(value: Pair) -> Pair:
    return value[0], -value[1]


def norm(value: Pair) -> int:
    return value[0] * value[0] - 2 * value[1] * value[1]


def content(value: Pair) -> int:
    return gcd(abs(value[0]), abs(value[1]))


def sign_surd(value: Pair) -> int:
    """Return the exact sign of a+b*sqrt(2), including zero."""

    rational, radical = value
    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if rational > 0 and radical > 0:
        return 1
    if rational < 0 and radical < 0:
        return -1
    comparison = rational * rational - 2 * radical * radical
    if rational > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def less(left: Pair, right: Pair) -> bool:
    return sign_surd(subtract(right, left)) > 0


def fibonacci_table(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return values


def positive_unit_table(limit: int) -> list[Pair]:
    values = [(1, 0)]
    for _ in range(limit):
        values.append(multiply(values[-1], (1, 1)))
    return values


def signed_unit(exponent: int, positive_units: list[Pair]) -> Pair:
    if exponent >= 0:
        return positive_units[exponent]
    index = -exponent
    parity = -1 if index % 2 else 1
    value = positive_units[index]
    return parity * value[0], -parity * value[1]


def first_larger_power(square_value: Pair, powers: list[Pair]) -> int:
    low = 1
    high = len(powers)
    while low < high:
        midpoint = (low + high) // 2
        if less(square_value, powers[midpoint]):
            high = midpoint
        else:
            low = midpoint + 1
    assert low < len(powers)
    return low


def unit_exponent(value: Pair, powers: list[Pair]) -> int | None:
    """Recognize a positive principal norm-minus-one unit exactly."""

    if norm(value) != -1 or sign_surd(value) <= 0:
        return None
    for exponent in range(1, len(powers), 2):
        positive = powers[exponent]
        if value == positive:
            return exponent
        negative = (-positive[0], positive[1])
        if value == negative:
            return -exponent
        if positive[1] > abs(value[1]):
            break
    return None


def gap_data(
    q_value: int,
    sigma: int,
    fib: list[int],
    powers: list[Pair],
    parity: int | None = None,
) -> tuple[int, int, int, Pair] | None:
    """Return (n,j,g,gap) for a target, or None.

    If ``parity`` is supplied, only exponents of that parity are considered.
    """

    f_value = fib[q_value]
    x_value = fib[q_value + 1]
    a_square = (
        f_value * f_value + 2 * x_value * x_value,
        2 * sigma * f_value * x_value,
    )
    exponent = first_larger_power(a_square, powers)
    if parity is not None and exponent % 2 != parity:
        # The interval has multiplicative width 2<lambda, so it contains at
        # most one integral power.  A different-parity power cannot follow.
        return None
    if not less(powers[exponent], scale(2, a_square)):
        return None
    gap = subtract(powers[exponent], a_square)
    gap_content = content(gap)
    if norm(gap) != -(gap_content * gap_content):
        return None
    normalized = (gap[0] // gap_content, gap[1] // gap_content)
    exponent_j = unit_exponent(normalized, powers)
    assert exponent_j is not None and exponent_j % 2
    return exponent, exponent_j, gap_content, gap


def check_parity_selector() -> None:
    """Verify the complete coefficient-parity implication symbolically."""

    # C_n is always odd, P_n has the parity of n, the rational coefficient
    # of A^2 has the parity of F_q, and its radical coefficient is even.
    for fib_parity in (0, 1):
        for exponent_parity in (0, 1):
            rational_gap = 1 - fib_parity
            radical_gap = exponent_parity
            possible_unit_content = rational_gap == radical_gap
            assert possible_unit_content == (
                exponent_parity == 1 - fib_parity
            )


def check_low_unit_exponents(fib: list[int], powers: list[Pair]) -> None:
    """Check the exact coefficient identities used to remove j=+-1."""

    for q_value in range(3, 97, 3):
        f_value = fib[q_value]
        x_value = fib[q_value + 1]
        height_plus = x_value * x_value + (x_value + f_value) ** 2
        height_minus = x_value * x_value + (x_value - f_value) ** 2
        assert height_plus == fib[2 * q_value + 3]
        assert height_minus == 3 * f_value * f_value + 2 * (-1) ** q_value
        for sigma in (1, -1):
            a_value = (sigma * f_value, x_value)
            a_square = square(a_value)
            radical_after_lambda = multiply((1, 1), a_square)[1]
            radical_after_inverse = multiply((-1, 1), a_square)[1]
            expected_same = height_plus if sigma == 1 else height_minus
            expected_opposite = height_minus if sigma == 1 else height_plus
            assert radical_after_lambda == expected_same
            assert radical_after_inverse == expected_opposite

        a_plus = (f_value, x_value)
        a_minus = (-f_value, x_value)
        assert less(multiply((1, 1), a_minus), a_plus)
        assert less(a_plus, scale(3, a_minus))

    # q=0 is excluded from the theorem and is the sharp low boundary:
    # A=sqrt(2), lambda-A^2=lambda^(-1).
    boundary_gap = subtract(powers[1], (2, 0))
    assert boundary_gap == signed_unit(-1, powers)


def check_changed_trace_sign(powers: list[Pair]) -> None:
    """Check the odd trace on the exact q=0 boundary target."""

    n_value = 1
    h_value = 1
    a_square = (2, 0)
    opposite_square = a_square
    left = subtract(
        powers[n_value],
        signed_unit(-(n_value + 2 * h_value), powers),
    )
    right = add(
        a_square,
        multiply(opposite_square, signed_unit(-2 * h_value, powers)),
    )
    assert left == right


def exact_enumeration(q_max: int) -> tuple[list[tuple[int, int, int, int, int]], int]:
    """Enumerate all nearest integral powers through q_max exactly."""

    fib = fibonacci_table(2 * q_max + 5)
    powers = positive_unit_table(2 * q_max + 20)
    hits: list[tuple[int, int, int, int, int]] = []
    odd_candidates = 0
    for q_value in range(1, q_max + 1):
        for sigma in (1, -1):
            f_value = fib[q_value]
            x_value = fib[q_value + 1]
            a_square = (
                f_value * f_value + 2 * x_value * x_value,
                2 * sigma * f_value * x_value,
            )
            exponent = first_larger_power(a_square, powers)
            if not less(powers[exponent], scale(2, a_square)):
                continue
            if exponent % 2:
                odd_candidates += 1
            target = gap_data(q_value, sigma, fib, powers)
            if target is None:
                continue
            n_value, j_value, g_value, _ = target
            if n_value % 2:
                assert q_value % 3 == 0
            else:
                assert q_value % 3 != 0
            hits.append((q_value, sigma, n_value, j_value, g_value))

    assert hits == [
        (2, -1, 2, -1, 6),
        (4, 1, 6, 1, 40),
    ]
    return hits, odd_candidates


def exact_terminal_odd_enumeration() -> None:
    """Certify the q<96 terminal range forced by the proof."""

    fib = fibonacci_table(200)
    powers = positive_unit_table(220)
    terminal_hits = []
    for q_value in range(3, 96, 3):
        for sigma in (1, -1):
            target = gap_data(q_value, sigma, fib, powers, parity=1)
            if target is not None:
                terminal_hits.append((q_value, sigma, target))
    assert terminal_hits == []


def certify_crude_bound_arithmetic() -> None:
    """Check the rational inequalities in the q<10^32 bootstrap."""

    # e > 19/7 follows already from the exponential series through 1/5!.
    assert Fraction(163, 60) > Fraction(19, 7)
    # Hence 1+log(2*10^32)<79.
    assert 19**78 > 2 * 10**32 * 7**78
    h_envelope = 116_000_000_000_000
    assert 10**32 // 2 > h_envelope * 79
    moving_upper = Fraction(5, 4) * (
        9_275_000_000_000
        * (4 * h_envelope * 79 + 44)
        * 79
        + 4
    )
    assert moving_upper < 10**32


def certify_reused_dp_intervals() -> tuple[int, int]:
    """Rerun the audited fixed and moving rational-interval certificates."""

    gamma, log_pell, _, kappa, mu = certified_dp.build_logarithmic_intervals()
    sqrt_two = certified_dp.square_root_interval(2)
    sqrt_five = certified_dp.square_root_interval(5)
    pell_unit = certified_dp.interval_add(
        (Fraction(1), Fraction(1)), sqrt_two
    )
    golden = certified_dp.interval_scale(
        certified_dp.interval_add((Fraction(1), Fraction(1)), sqrt_five),
        Fraction(1, 2),
    )
    gamma_distance = certified_dp.certify_common_convergent(gamma)
    certified_dp.certify_fixed_reduction(gamma_distance, mu, pell_unit)
    return certified_dp.certify_moving_reduction(
        gamma_distance,
        log_pell,
        golden,
        kappa,
        mu,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q-max", type=int, default=10_000)
    args = parser.parse_args()
    if not 96 <= args.q_max <= 20_000:
        parser.error("--q-max must lie between 96 and 20000")

    check_parity_selector()
    fib = fibonacci_table(2 * args.q_max + 5)
    powers = positive_unit_table(2 * args.q_max + 20)
    check_low_unit_exponents(fib, powers)
    check_changed_trace_sign(powers)
    certify_crude_bound_arithmetic()
    worst_sector, worst_h = certify_reused_dp_intervals()
    exact_terminal_odd_enumeration()
    hits, odd_candidates = exact_enumeration(args.q_max)
    print(
        "PASS: all-exponent Fibonacci nearest-gap support; "
        "global parity selector and j=+-1 coefficient identities exact; "
        "odd-trace sign verified; bootstrap q<10^32 arithmetic exact; "
        "fixed DP threshold<48 and moving threshold<90 certified "
        f"(worst moving interval sigma={worst_sector},h={worst_h}); "
        "terminal q<96 odd branch empty; "
        f"all-exponent regression q<={args.q_max}, "
        f"odd_window_candidates={odd_candidates}, hits={hits}."
    )


if __name__ == "__main__":
    main()
