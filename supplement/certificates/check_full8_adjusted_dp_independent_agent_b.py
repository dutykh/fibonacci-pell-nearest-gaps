#!/usr/bin/env python3
"""Independent fixed-point audit of the adjusted full-block closure.

This checker uses outward-rounded base-10 fixed-point intervals. Its interval
engine is structurally different from the source checker's ``Fraction``
engine. It certifies the shared convergent, both fixed-shift epsilon bounds,
all 190 moving-shift epsilon bounds, both strict exponent thresholds, and the
exact terminal C0068 enumeration.

Run from ``math-sandbox/`` with

    python3 -B scripts/check_full8_adjusted_dp_independent_agent_b.py
"""

from __future__ import annotations

from math import gcd, isqrt


DIGITS = 100
SCALE = 10**DIGITS
ONE = (SCALE, SCALE)

CONVERGENT_NUMERATOR = 2_585_762_774_943_540_084_566_744_409_566_833
CONVERGENT_DENOMINATOR = 4_736_007_914_708_481_810_221_186_278_309_584
COEFFICIENT_CAP = 2 * 10**32

Interval = tuple[int, int]


def ceil_div(numerator: int, denominator: int) -> int:
    """Return the mathematical ceiling of a rational integer quotient."""
    assert denominator > 0
    return -((-numerator) // denominator)


def add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def negate(value: Interval) -> Interval:
    return -value[1], -value[0]


def subtract(left: Interval, right: Interval) -> Interval:
    return add(left, negate(right))


def multiply(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products) // SCALE, ceil_div(max(products), SCALE)


def multiply_integer(value: Interval, integer: int) -> Interval:
    if integer >= 0:
        return value[0] * integer, value[1] * integer
    return value[1] * integer, value[0] * integer


def divide_integer(value: Interval, integer: int) -> Interval:
    assert integer > 0
    return value[0] // integer, ceil_div(value[1], integer)


def reciprocal(value: Interval) -> Interval:
    assert value[0] > 0
    return SCALE * SCALE // value[1], ceil_div(SCALE * SCALE, value[0])


def divide(left: Interval, right: Interval) -> Interval:
    return multiply(left, reciprocal(right))


def power(value: Interval, exponent: int) -> Interval:
    assert exponent >= 0
    answer = ONE
    base = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        remaining >>= 1
    return answer


def square_root_integer(value: int) -> Interval:
    lower = isqrt(value * SCALE * SCALE)
    return lower, lower + 1


def atanh_log_on_one_to_two(value: Interval, terms: int = 70) -> Interval:
    """Certify log(value) when 1 <= value <= 2 by the atanh series."""
    assert SCALE <= value[0] <= value[1] <= 2 * SCALE
    z_value = divide(subtract(value, ONE), add(value, ONE))
    assert z_value[0] >= 0 and z_value[1] <= SCALE // 3 + 2
    z_square = multiply(z_value, z_value)
    term = z_value
    partial = (0, 0)
    for index in range(terms):
        partial = add(partial, divide_integer(term, 2 * index + 1))
        term = multiply(term, z_square)
    partial = multiply_integer(partial, 2)
    # Since z <= 1/3, twice the omitted tail is less than
    # 3*z^(2N+1)/(2N+1). Outward rounding in the partial sum is retained.
    tail = ceil_div(3 * term[1], 2 * terms + 1)
    return partial[0] - tail, partial[1] + tail


LOG_TWO = atanh_log_on_one_to_two((2 * SCALE, 2 * SCALE))


def log_interval(value: Interval, terms: int = 70) -> Interval:
    """Certified natural logarithm of a positive narrow interval."""
    assert value[0] > 0
    reduced = value
    exponent = 0
    while reduced[1] < SCALE:
        reduced = multiply_integer(reduced, 2)
        exponent -= 1
    while reduced[0] >= 2 * SCALE:
        reduced = divide_integer(reduced, 2)
        exponent += 1
    assert SCALE <= reduced[0] and reduced[1] <= 2 * SCALE
    answer = atanh_log_on_one_to_two(reduced, terms)
    return add(answer, multiply_integer(LOG_TWO, exponent))


def nearest_distance(value: Interval) -> tuple[Interval, int]:
    """Return a certified distance to the unique nearest integer cell."""
    nearest = (value[0] + value[1] + SCALE) // (2 * SCALE)
    center = nearest * SCALE
    assert value[0] > center - SCALE // 2
    assert value[1] < center + SCALE // 2
    if value[1] <= center:
        return (center - value[1], center - value[0]), nearest
    if value[0] >= center:
        return (value[0] - center, value[1] - center), nearest
    return (0, max(center - value[0], value[1] - center)), nearest


def build_constants() -> tuple[
    Interval,
    Interval,
    Interval,
    dict[int, Interval],
    dict[int, Interval],
]:
    sqrt_two = square_root_integer(2)
    sqrt_five = square_root_integer(5)
    pell_unit = add(ONE, sqrt_two)
    golden = divide_integer(add(ONE, sqrt_five), 2)
    log_pell = log_interval(pell_unit)
    gamma = divide(log_interval(golden), log_pell)
    root_product = multiply(sqrt_two, golden)

    mu: dict[int, Interval] = {}
    kappa: dict[int, Interval] = {}
    for sigma in (1, -1):
        if sigma == 1:
            u_value = add(ONE, root_product)
            abs_w = subtract(root_product, ONE)
        else:
            u_value = subtract(root_product, ONE)
            abs_w = add(ONE, root_product)
        c_value = divide_integer(power(u_value, 2), 5)
        kappa[sigma] = power(divide(abs_w, u_value), 2)
        mu[sigma] = divide(log_interval(c_value), log_pell)

    assert gamma[0] > SCALE // 2
    assert kappa[1][1] < SCALE
    assert kappa[-1][1] < 9 * SCALE
    return gamma, pell_unit, golden, mu, kappa


def certify_convergent(gamma: Interval) -> Interval:
    numerator = CONVERGENT_NUMERATOR
    denominator = CONVERGENT_DENOMINATOR
    assert gcd(numerator, denominator) == 1
    assert denominator > 6 * COEFFICIENT_CAP

    # Legendre's strict 1/(2Q^2) criterion, checked by cross multiplication.
    assert (gamma[0] * denominator - numerator * SCALE) * (2 * denominator) > -SCALE
    assert (gamma[1] * denominator - numerator * SCALE) * (2 * denominator) < SCALE

    distance, nearest = nearest_distance(multiply_integer(gamma, denominator))
    assert nearest == numerator
    return distance


def certify_reductions(
    gamma_distance: Interval,
    pell_unit: Interval,
    golden: Interval,
    mu: dict[int, Interval],
    kappa: dict[int, Interval],
) -> tuple[int, int]:
    denominator = CONVERGENT_DENOMINATOR
    fixed_eps: dict[int, int] = {}
    for sigma in (1, -1):
        mu_distance, _ = nearest_distance(multiply_integer(mu[sigma], denominator))
        epsilon = mu_distance[0] - COEFFICIENT_CAP * gamma_distance[1]
        fixed_eps[sigma] = epsilon
        assert epsilon * 100 > 7 * SCALE
        # 42Q/epsilon < lambda^96 gives a strict DP threshold below 48
        # with base lambda^2.
        assert 42 * denominator * SCALE * SCALE < epsilon * power(pell_unit, 96)[0]
    assert fixed_eps[1] < fixed_eps[-1]

    pell_inverse = reciprocal(pell_unit)
    current = power(pell_inverse, 6)
    step = power(pell_inverse, 4)
    minimum: tuple[int, int, int] | None = None
    for h_value in range(3, 192, 2):
        for sigma in (1, -1):
            correction = log_interval(add(ONE, multiply(kappa[sigma], current)), 45)
            moving_mu = add(mu[sigma], divide(correction, log_interval(pell_unit)))
            moving_distance, _ = nearest_distance(
                multiply_integer(moving_mu, denominator)
            )
            epsilon = moving_distance[0] - COEFFICIENT_CAP * gamma_distance[1]
            assert epsilon * 1000 > 9 * SCALE
            # 33Q/epsilon < phi^180 gives a strict DP threshold below 90
            # with base phi^2.
            assert 33 * denominator * SCALE * SCALE < epsilon * power(golden, 180)[0]
            if minimum is None or epsilon < minimum[0]:
                minimum = epsilon, sigma, h_value
        current = multiply(current, step)
    assert minimum is not None
    assert (minimum[1], minimum[2]) == (1, 31)
    return minimum[1], minimum[2]


def fibonacci(index: int) -> int:
    first, second = 0, 1
    for _ in range(index):
        first, second = second, first + second
    return first


def sign_quadratic(value: tuple[int, int]) -> int:
    """Exact sign of a+b*sqrt(2)."""
    rational, radical = value
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


def exact_gap_enumeration() -> list[tuple[int, int, int]]:
    """Independently enumerate the exact C0068 criterion for q < 90."""
    hits: list[tuple[int, int, int]] = []
    for q_value in range(2, 90):
        if q_value % 3 == 0:
            continue
        f_value = fibonacci(q_value)
        x_value = fibonacci(q_value + 1)
        for sigma in (1, -1):
            d_square = (
                f_value * f_value + 2 * x_value * x_value,
                2 * sigma * f_value * x_value,
            )
            unit = (1, 0)
            exponent = 0
            while sign_quadratic(
                (unit[0] - d_square[0], unit[1] - d_square[1])
            ) <= 0:
                unit = (3 * unit[0] + 4 * unit[1], 2 * unit[0] + 3 * unit[1])
                exponent += 2
            if sign_quadratic(
                (unit[0] - 2 * d_square[0], unit[1] - 2 * d_square[1])
            ) >= 0:
                continue
            delta = unit[0] - d_square[0], unit[1] - d_square[1]
            content = gcd(abs(delta[0]), abs(delta[1]))
            if delta[0] * delta[0] - 2 * delta[1] * delta[1] == -(content**2):
                hits.append((q_value, sigma, exponent))
    assert hits == [(2, -1, 2), (4, 1, 6)]
    return hits


def main() -> None:
    gamma, pell_unit, golden, mu, kappa = build_constants()
    gamma_distance = certify_convergent(gamma)
    worst = certify_reductions(
        gamma_distance, pell_unit, golden, mu, kappa
    )
    hits = exact_gap_enumeration()
    print(
        "PASS: independent outward-rounded fixed-point audit; common "
        "convergent certified by Legendre; both fixed epsilons give w<48; "
        f"all moving epsilons give q<90 with worst (sigma,h)={worst}; "
        f"exact C0068 hits below 90 are {hits}, so none has q>=7."
    )


if __name__ == "__main__":
    main()
