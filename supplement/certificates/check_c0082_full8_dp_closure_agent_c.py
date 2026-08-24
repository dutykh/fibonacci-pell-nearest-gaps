#!/usr/bin/env python3
"""Certified continued-fraction closure of the C0082 full block.

The checker uses rational intervals throughout. It certifies one convergent
of log(phi)/log(lambda), the two fixed-coefficient Dujella--Petho epsilon
bounds, all moving-coefficient epsilon bounds for odd 3 <= h < 192, and an
exact C0068 enumeration for q < 90.

It is a checker for the proof recorded in
attempts/2026-08-24-full8-dominance-recon-agent-c.md. It performs no
unbounded search and uses only the Python standard library.

Deterministic command from math-sandbox/:
    python3 -B scripts/check_c0082_full8_dp_closure_agent_c.py
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt


Q = Fraction
Interval = tuple[Q, Q]

CONVERGENT_NUMERATOR = 2_585_762_774_943_540_084_566_744_409_566_833
CONVERGENT_DENOMINATOR = 4_736_007_914_708_481_810_221_186_278_309_584
EXPONENT_CAP = 200_000_000_000_000_000_000_000_000_000_000


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


def interval_divide(left: Interval, right: Interval) -> Interval:
    assert right[0] > 0
    reciprocals = (Q(1, right[1]), Q(1, right[0]))
    return interval_multiply(left, reciprocals)


def interval_power(interval: Interval, exponent: int) -> Interval:
    assert interval[0] >= 0 and exponent >= 0
    answer = (Q(1), Q(1))
    base = interval
    power = exponent
    while power:
        if power & 1:
            answer = interval_multiply(answer, base)
        base = interval_multiply(base, base)
        power >>= 1
    return answer


def round_interval_outward(interval: Interval, digits: int = 75) -> Interval:
    """Replace large exact denominators by a certified fixed decimal grid."""
    scale = 10**digits
    lower_numerator = (interval[0].numerator * scale) // interval[0].denominator
    upper_numerator = -(
        (-interval[1].numerator * scale) // interval[1].denominator
    )
    return Q(lower_numerator, scale), Q(upper_numerator, scale)


def square_root_interval(radicand: int, digits: int = 90) -> Interval:
    denominator = 10**digits
    numerator = isqrt(radicand * denominator * denominator)
    return Q(numerator, denominator), Q(numerator + 1, denominator)


def log_point_interval(value: Q, terms: int) -> Interval:
    """Certified log interval from the atanh series at a rational point."""
    assert value > 0
    z_value = (value - 1) / (value + 1)
    z_square = z_value * z_value
    term = z_value
    partial = Q(0)
    for index in range(terms):
        partial += term / (2 * index + 1)
        term *= z_square
    approximation = 2 * partial
    remainder = 2 * abs(term) / ((2 * terms + 1) * (1 - z_square))
    return approximation - remainder, approximation + remainder


def log_interval(interval: Interval, terms: int = 60) -> Interval:
    """Use monotonicity of log and certified endpoint series."""
    lower = log_point_interval(interval[0], terms)[0]
    upper = log_point_interval(interval[1], terms)[1]
    assert lower < upper
    return lower, upper


def nearest_integer(value: Q) -> int:
    return (value.numerator * 2 + value.denominator) // (2 * value.denominator)


def distance_to_integer_interval(interval: Interval) -> tuple[Interval, int]:
    """Certify ||x|| on an interval contained in one nearest-integer cell."""
    midpoint = (interval[0] + interval[1]) / 2
    integer = nearest_integer(midpoint)
    assert interval[0] > integer - Q(1, 2)
    assert interval[1] < integer + Q(1, 2)
    if interval[1] <= integer:
        return (integer - interval[1], integer - interval[0]), integer
    if interval[0] >= integer:
        return (interval[0] - integer, interval[1] - integer), integer
    return (Q(0), max(integer - interval[0], interval[1] - integer)), integer


def scaled_interval(interval: Interval, integer: int) -> Interval:
    return interval_scale(interval, Q(integer))


def build_logarithmic_intervals() -> tuple[
    Interval,
    Interval,
    dict[int, Interval],
    dict[int, Interval],
    dict[int, Interval],
]:
    sqrt_two = square_root_interval(2)
    sqrt_five = square_root_interval(5)
    one = (Q(1), Q(1))
    two = (Q(2), Q(2))
    five = (Q(5), Q(5))

    pell_unit = interval_add(one, sqrt_two)
    golden = interval_scale(interval_add(one, sqrt_five), Q(1, 2))
    root_product = interval_multiply(sqrt_two, golden)
    # Scaling keeps every atanh parameter small. In particular, it avoids
    # applying the series directly at lambda (z about 0.414) or c_- (z about
    # -0.502), which would make an exact-Fraction certificate unnecessarily
    # expensive.
    log_two = log_interval(two, terms=82)
    log_pell = interval_add(
        log_interval(interval_scale(pell_unit, Q(1, 2)), terms=42),
        log_two,
    )
    log_golden = log_interval(golden, terms=64)
    gamma = interval_divide(log_golden, log_pell)

    c_intervals: dict[int, Interval] = {}
    kappa_intervals: dict[int, Interval] = {}
    mu_intervals: dict[int, Interval] = {}
    for sigma in (1, -1):
        if sigma == 1:
            u_value = interval_add(one, root_product)
            abs_w = (root_product[0] - 1, root_product[1] - 1)
        else:
            u_value = (root_product[0] - 1, root_product[1] - 1)
            abs_w = interval_add(one, root_product)
        c_value = interval_divide(interval_power(u_value, 2), five)
        ratio = interval_divide(abs_w, u_value)
        kappa = interval_power(ratio, 2)
        c_intervals[sigma] = c_value
        kappa_intervals[sigma] = kappa
        if sigma == 1:
            log_c = interval_add(
                log_interval(interval_scale(c_value, Q(1, 2)), terms=42),
                log_two,
            )
        else:
            log_c = interval_add(
                log_interval(interval_scale(c_value, Q(4)), terms=48),
                interval_scale(log_two, Q(-2)),
            )
        mu_intervals[sigma] = interval_divide(log_c, log_pell)

    # These elementary bounds are used in the written first reduction.
    assert c_intervals[1][0] > Q(1, 5)
    assert c_intervals[-1][0] > Q(1, 5)
    assert kappa_intervals[1][1] < 1
    assert kappa_intervals[-1][1] < 9
    assert gamma[0] > Q(1, 2)
    return gamma, log_pell, c_intervals, kappa_intervals, mu_intervals


def certify_common_convergent(gamma: Interval) -> Interval:
    numerator = CONVERGENT_NUMERATOR
    denominator = CONVERGENT_DENOMINATOR
    assert gcd(numerator, denominator) == 1
    assert denominator > 6 * EXPONENT_CAP

    rational = Q(numerator, denominator)
    legendre_radius = Q(1, 2 * denominator * denominator)
    assert gamma[0] > rational - legendre_radius
    assert gamma[1] < rational + legendre_radius

    gamma_scaled = scaled_interval(gamma, denominator)
    gamma_distance, nearest = distance_to_integer_interval(gamma_scaled)
    assert nearest == numerator
    return gamma_distance


def certify_fixed_reduction(
    gamma_distance: Interval,
    mu_intervals: dict[int, Interval],
    pell_unit: Interval,
) -> None:
    denominator = CONVERGENT_DENOMINATOR
    epsilon_lowers = {}
    for sigma in (1, -1):
        mu_scaled = scaled_interval(mu_intervals[sigma], denominator)
        mu_distance, _ = distance_to_integer_interval(mu_scaled)
        epsilon_lower = mu_distance[0] - EXPONENT_CAP * gamma_distance[1]
        epsilon_lowers[sigma] = epsilon_lower
        assert epsilon_lower > Q(7, 100)

        # This is the logarithm-free certificate that the DP threshold is <48.
        left = Q(42 * denominator, 1) / epsilon_lower
        right_lower = interval_power(pell_unit, 96)[0]
        assert left < right_lower

    assert epsilon_lowers[1] < epsilon_lowers[-1]


def certify_moving_reduction(
    gamma_distance: Interval,
    log_pell: Interval,
    golden: Interval,
    kappa_intervals: dict[int, Interval],
    mu_intervals: dict[int, Interval],
) -> tuple[int, int]:
    denominator = CONVERGENT_DENOMINATOR
    sqrt_two = square_root_interval(2)
    pell_unit = interval_add((Q(1), Q(1)), sqrt_two)
    pell_inverse = interval_divide((Q(1), Q(1)), pell_unit)
    pell_inverse_square = interval_power(pell_inverse, 2)
    step = interval_power(pell_inverse_square, 2)
    current_power = interval_power(pell_inverse_square, 3)
    right_lower = interval_power(golden, 180)[0]

    minimum: tuple[Q, int, int] | None = None
    for h in range(3, 192, 2):
        for sigma in (1, -1):
            if h < 55:
                t_value = interval_multiply(kappa_intervals[sigma], current_power)
                t_value = round_interval_outward(t_value)
                one_plus_t = interval_add((Q(1), Q(1)), t_value)
                log_correction = log_interval(one_plus_t, terms=18)
                correction_quotient = interval_divide(log_correction, log_pell)
            else:
                # kappa<9, lambda^2>5, log(1+t)<t, and log(lambda)>4/5.
                # This coarse rational interval is already narrower than
                # 0.0002/Q at h=55 and rapidly improves thereafter.
                correction_quotient = (Q(0), Q(45, 4 * 5**h))
            nu_interval = interval_add(mu_intervals[sigma], correction_quotient)
            nu_scaled = scaled_interval(nu_interval, denominator)
            nu_distance, _ = distance_to_integer_interval(nu_scaled)
            epsilon_lower = nu_distance[0] - EXPONENT_CAP * gamma_distance[1]
            assert epsilon_lower > Q(9, 1000)

            left = Q(33 * denominator, 1) / epsilon_lower
            assert left < right_lower

            if minimum is None or epsilon_lower < minimum[0]:
                minimum = (epsilon_lower, sigma, h)
        if h < 55:
            current_power = interval_multiply(current_power, step)

    assert minimum is not None
    assert (minimum[1], minimum[2]) == (1, 31)
    return minimum[1], minimum[2]


def sign_quadratic(pair: tuple[int, int]) -> int:
    """Exact sign of a+b*sqrt(2)."""
    rational, radical = pair
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if rational == 0:
        return (radical > 0) - (radical < 0)
    if rational > 0 and radical > 0:
        return 1
    if rational < 0 and radical < 0:
        return -1
    comparison = rational * rational - 2 * radical * radical
    if rational > 0:
        return (comparison > 0) - (comparison < 0)
    return (comparison < 0) - (comparison > 0)


def subtract_pair(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] - right[0], left[1] - right[1]


def scale_pair(pair: tuple[int, int], scalar: int) -> tuple[int, int]:
    return scalar * pair[0], scalar * pair[1]


def exact_nearest_gap_enumeration() -> list[tuple[int, int, int]]:
    """Enumerate C0068 targets exactly for 2 <= q < 90 and 3 not dividing q."""
    hits: list[tuple[int, int, int]] = []
    for q_index in range(2, 90):
        if q_index % 3 == 0:
            continue
        f_value = fibonacci(q_index)
        x_value = fibonacci(q_index + 1)
        for sigma in (1, -1):
            d_square = (
                f_value * f_value + 2 * x_value * x_value,
                2 * sigma * f_value * x_value,
            )

            # Advance lambda^(2r) monotonically until it first exceeds D^2.
            pell_even = (1, 0)
            exponent = 0
            while sign_quadratic(subtract_pair(pell_even, d_square)) <= 0:
                c_value, p_value = pell_even
                pell_even = (3 * c_value + 4 * p_value, 2 * c_value + 3 * p_value)
                exponent += 2

            if sign_quadratic(subtract_pair(pell_even, scale_pair(d_square, 2))) >= 0:
                continue

            delta = subtract_pair(pell_even, d_square)
            content = gcd(abs(delta[0]), abs(delta[1]))
            norm = delta[0] * delta[0] - 2 * delta[1] * delta[1]
            if norm == -(content * content):
                hits.append((q_index, sigma, exponent))

    assert hits == [(2, -1, 2), (4, 1, 6)]
    assert all(q_index < 7 for q_index, _, _ in hits)
    return hits


def fibonacci(index: int) -> int:
    first, second = 0, 1
    for _ in range(index):
        first, second = second, first + second
    return first


def main() -> None:
    gamma, log_pell, _, kappa_intervals, mu_intervals = (
        build_logarithmic_intervals()
    )
    sqrt_two = square_root_interval(2)
    sqrt_five = square_root_interval(5)
    pell_unit = interval_add((Q(1), Q(1)), sqrt_two)
    golden = interval_scale(interval_add((Q(1), Q(1)), sqrt_five), Q(1, 2))
    gamma_distance = certify_common_convergent(gamma)
    certify_fixed_reduction(gamma_distance, mu_intervals, pell_unit)
    minimum_sector, minimum_h = certify_moving_reduction(
        gamma_distance,
        log_pell,
        golden,
        kappa_intervals,
        mu_intervals,
    )
    hits = exact_nearest_gap_enumeration()
    print(
        "PASS: certified full-eight DP closure; fixed threshold<48; "
        f"moving threshold<90 with worst interval at (sigma,h)="
        f"({minimum_sector},{minimum_h}); exact q<90 C0068 hits={hits}; "
        "no residual q>=7 target remains."
    )


if __name__ == "__main__":
    main()
