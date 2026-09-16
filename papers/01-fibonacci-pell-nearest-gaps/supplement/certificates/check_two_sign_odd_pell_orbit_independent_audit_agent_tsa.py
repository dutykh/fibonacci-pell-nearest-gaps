#!/usr/bin/env python3
"""Independent exact audit of the proposed two-sign odd-Pell theorem.

This checker deliberately imports no project module.  It reconstructs the
Fibonacci and Pell sequences, the orbit map, and arithmetic in Z[sqrt(2)]
from their defining recurrences.  The finite grid is falsification evidence;
the universal result rests on the proof audited in the companion attempt.
"""

from __future__ import annotations

from math import gcd


Q_MAX = 80
D_MAX = 119
T_MAX = 40


def sequence(first: int, second: int, coefficient: int, limit: int) -> list[int]:
    """Return a second-order recurrence through ``limit`` inclusive."""

    values = [first, second]
    while len(values) <= limit:
        values.append(coefficient * values[-1] + values[-2])
    return values


FIB = sequence(0, 1, 1, 404)
PELL = sequence(0, 1, 2, 404)


def companion(index: int) -> int:
    """Return C_index, including the convention C_0=1."""

    assert index >= 0
    previous = 1 if index == 0 else PELL[index - 1]
    return PELL[index] + previous


def unit_power(exponent: int) -> tuple[int, int]:
    """Return coefficients of (1+sqrt(2))**exponent."""

    if exponent >= 0:
        return companion(exponent), PELL[exponent]
    index = -exponent
    sign = -1 if index % 2 else 1
    return sign * companion(index), -sign * PELL[index]


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] - right[0], left[1] - right[1]


def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (
        left[0] * right[0] + 2 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def square(value: tuple[int, int]) -> tuple[int, int]:
    return multiply(value, value)


def scale(multiplier: int, value: tuple[int, int]) -> tuple[int, int]:
    return multiplier * value[0], multiplier * value[1]


def conjugate(value: tuple[int, int]) -> tuple[int, int]:
    return value[0], -value[1]


def norm(value: tuple[int, int]) -> int:
    return value[0] ** 2 - 2 * value[1] ** 2


def content(value: tuple[int, int]) -> int:
    return gcd(abs(value[0]), abs(value[1]))


def sign_surd(value: tuple[int, int]) -> int:
    """Return the exact sign of a+b*sqrt(2)."""

    rational, radical = value
    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if rational > 0 and radical > 0:
        return 1
    if rational < 0 and radical < 0:
        return -1
    if rational > 0:
        return 1 if rational * rational > 2 * radical * radical else -1
    return 1 if 2 * radical * radical > rational * rational else -1


def forms(point: tuple[int, int]) -> tuple[int, int, int]:
    """Return (H,B,Q) at an integer point."""

    u, v = point
    return (
        u * u + v * v,
        v * v - u * u + 2 * u * v,
        v * v - 2 * u * v - u * u,
    )


def step(point: tuple[int, int]) -> tuple[int, int]:
    u, v = point
    return v, u + 2 * v


def inverse_step(point: tuple[int, int]) -> tuple[int, int]:
    u, v = point
    return v - 2 * u, u


def seed(q_index: int, sigma: int) -> tuple[int, int]:
    return FIB[q_index + 1], FIB[q_index + 1] + sigma * FIB[q_index]


def positive_integral_quotient(
    dividend: tuple[int, int], divisor: tuple[int, int]
) -> int | None:
    """Return positive integer g if dividend=g*divisor, else None."""

    assert divisor[0] != 0 and divisor[1] != 0
    if dividend[0] % divisor[0] or dividend[1] % divisor[1]:
        return None
    first = dividend[0] // divisor[0]
    second = dividend[1] // divisor[1]
    return first if first == second and first > 0 else None


def check_static_algebra() -> None:
    """Check seed signs, form identities, and window inequalities."""

    assert unit_power(0) == (1, 0)
    assert unit_power(1) == (1, 1)
    assert unit_power(-1) == (-1, 1)
    for exponent in range(-200, 201):
        power = unit_power(exponent)
        expected_norm = -1 if exponent % 2 else 1
        assert norm(power) == expected_norm
        assert multiply(power, unit_power(-exponent)) == (1, 0)

    for q_index in range(1, Q_MAX + 1):
        f = FIB[q_index]
        x = FIB[q_index + 1]
        plus = seed(q_index, 1)
        minus = seed(q_index, -1)
        assert plus == (x, FIB[q_index + 2])
        assert minus == (x, FIB[q_index - 1])

        for sigma, point in ((1, plus), (-1, minus)):
            height, bilinear, invariant = forms(point)
            assert invariant == f * f - 2 * x * x
            assert bilinear * bilinear + invariant * invariant == 2 * height * height
            gamma = (-sigma * f, x)
            a_sigma = (sigma * f, x)
            assert conjugate(gamma) == scale(-1, a_sigma)
            assert multiply(unit_power(-1), square(gamma)) == (-bilinear, height)

        if q_index >= 2:
            invariant = forms(plus)[2]
            product = (
                FIB[q_index - 1]
                * FIB[q_index + 1]
                * FIB[q_index + 2]
                * FIB[q_index + 4]
            )
            assert invariant * invariant - 1 == product > 0

            a_plus = (f, x)
            a_minus = (-f, x)
            assert sign_surd(a_minus) > 0
            assert sign_surd(subtract(a_plus, a_minus)) > 0
            assert sign_surd(subtract(scale(3, a_minus), a_plus)) > 0
            for crossing_time in range(1, T_MAX + 1):
                exponent = 2 * crossing_time + 1
                assert sign_surd(
                    subtract(a_minus, multiply(a_plus, unit_power(-exponent)))
                ) > 0

    assert sign_surd(subtract(unit_power(3), (3, 0))) > 0


def check_evolution_and_enumeration() -> tuple[int, list[tuple[int, int, int, int]]]:
    """Check the bridge, its anchored converse, and the stated finite grid."""

    comparisons = 0
    hits: list[tuple[int, int, int, int]] = []
    for q_index in range(1, Q_MAX + 1):
        f = FIB[q_index]
        x = FIB[q_index + 1]
        for sigma in (1, -1):
            origin = seed(q_index, sigma)
            height_zero, bilinear_zero, invariant_zero = forms(origin)
            point = origin
            for crossing_time in range(T_MAX + 1):
                height_now, bilinear_now, invariant_now = forms(point)
                even_index = 2 * crossing_time
                c_even = companion(even_index)
                p_even = PELL[even_index]
                assert height_now == c_even * height_zero + p_even * bilinear_zero
                assert bilinear_now == 2 * p_even * height_zero + c_even * bilinear_zero
                assert invariant_now == (-1) ** crossing_time * invariant_zero

                for d_index in range(1, D_MAX + 1, 2):
                    comparisons += 1
                    target_index = d_index + 2 * crossing_time
                    hit = height_now == PELL[target_index]
                    if hit:
                        hits.append((sigma, q_index, d_index, crossing_time))

                    # The Pell pair evolves under exactly the same matrix.
                    assert PELL[target_index] == (
                        c_even * PELL[d_index]
                        + p_even * companion(d_index)
                    )
                    assert companion(target_index) == (
                        2 * p_even * PELL[d_index]
                        + c_even * companion(d_index)
                    )

                    if q_index == 1:
                        continue

                    difference_d = height_zero - PELL[d_index]
                    difference_e = companion(d_index) - bilinear_zero
                    r_index = (d_index - 1) // 2
                    gamma = (-sigma * f, x)
                    z_value = multiply(unit_power(r_index), gamma)
                    unconditional_left = subtract(square(z_value), (1, 0))
                    unconditional_right = multiply(
                        unit_power(d_index), (difference_e, difference_d)
                    )
                    assert unconditional_left == unconditional_right

                    a_sigma = (sigma * f, x)
                    gap = subtract(unit_power(d_index - 1), square(a_sigma))
                    odd_unit = unit_power(-(2 * crossing_time + 1))
                    gap_multiplier = positive_integral_quotient(gap, odd_unit)
                    # This is the exact fixed-anchor converse, not merely a
                    # test that every crossing maps to a gap.
                    assert (gap_multiplier is not None) == hit

                    if not hit:
                        continue

                    assert invariant_zero * invariant_zero > 1
                    assert abs(bilinear_now) < companion(target_index)
                    if crossing_time == 0:
                        assert difference_d == 0
                        assert difference_e > 0
                        bridge_multiplier = difference_e
                    else:
                        assert c_even * difference_d == p_even * difference_e
                        assert (
                            bilinear_now - companion(target_index)
                            == -difference_e // c_even
                        )
                        assert difference_d > 0 and difference_e > 0
                        bridge_multiplier = difference_d // p_even
                        assert difference_d % p_even == 0
                    assert difference_d == bridge_multiplier * p_even
                    assert difference_e == bridge_multiplier * c_even
                    assert bridge_multiplier == gcd(difference_d, difference_e)
                    assert gap_multiplier == bridge_multiplier
                    assert gap == scale(bridge_multiplier, odd_unit)
                    assert content(gap) == bridge_multiplier
                    assert norm(gap) == -(bridge_multiplier**2)

                    # Every plus hit and every positive-time minus hit must be
                    # in the strict nearest window.
                    if sigma == 1 or crossing_time >= 1:
                        a_square = square(a_sigma)
                        power = unit_power(d_index - 1)
                        assert sign_surd(subtract(power, a_square)) > 0
                        assert sign_surd(subtract(scale(2, a_square), power)) > 0

                point = step(point)

    expected = [(-1, 1, 1, 0), (-1, 2, 3, 0), (-1, 4, 5, 0)]
    expected.extend((1, 1, 3, crossing_time) for crossing_time in range(T_MAX + 1))
    assert sorted(hits) == sorted(expected)
    return comparisons, hits


def check_boundaries() -> tuple[int, list[tuple[int, int]]]:
    """Check q=1, all direct hits in a larger grid, and negative times."""

    for crossing_time in range(T_MAX + 41):
        plus_point = seed(1, 1)
        minus_point = seed(1, -1)
        for _ in range(crossing_time):
            plus_point = step(plus_point)
            minus_point = step(minus_point)
        assert forms(plus_point)[0] == PELL[3 + 2 * crossing_time]
        if crossing_time == 0:
            assert minus_point == (1, 0)
            assert forms(minus_point)[0] == PELL[1]
        else:
            assert minus_point == (PELL[crossing_time - 1], PELL[crossing_time])
            assert forms(minus_point)[0] == PELL[2 * crossing_time - 1]

    direct_hits: list[tuple[int, int]] = []
    direct_comparisons = 0
    for q_index in range(1, 201):
        direct_height = forms(seed(q_index, -1))[0]
        assert direct_height == 3 * FIB[q_index] ** 2 + 2 * (-1) ** q_index
        for d_index in range(1, 400, 2):
            direct_comparisons += 1
            if direct_height == PELL[d_index]:
                direct_hits.append((q_index, d_index))
    assert direct_hits == [(1, 1), (2, 3), (4, 5)]

    assert inverse_step(seed(2, 1)) == (-1, 2)
    assert forms(inverse_step(seed(2, 1)))[0] == PELL[3]
    assert 5 + 2 * (-1) == 3
    assert inverse_step(seed(4, 1)) == (-2, 5)
    assert forms(inverse_step(seed(4, 1)))[0] == PELL[5]
    assert 7 + 2 * (-1) == 5

    a_four_minus = (-FIB[4], FIB[5])
    gap = subtract(unit_power(4), square(a_four_minus))
    assert gap == scale(42, unit_power(-1))
    upper_residual = subtract(unit_power(4), scale(2, square(a_four_minus)))
    assert upper_residual == (-101, 72)
    assert sign_surd(upper_residual) > 0

    # The other nontrivial direct row is inside the nearest window and has
    # the normalized quotient specified by C0088.
    a_two_minus = (-FIB[2], FIB[3])
    assert subtract(unit_power(2), square(a_two_minus)) == scale(
        6, unit_power(-1)
    )

    return direct_comparisons, direct_hits


def main() -> None:
    check_static_algebra()
    comparisons, hits = check_evolution_and_enumeration()
    direct_comparisons, direct_hits = check_boundaries()
    print(
        "PASS: independent two-sign odd-Pell audit; "
        f"orbit_comparisons={comparisons}; "
        f"orbit_hits={hits}; direct_comparisons={direct_comparisons}; "
        f"direct_hits={direct_hits}; exact converse, q=1 boundaries, "
        "nearest-window inequalities, and negative-time witnesses verified."
    )


if __name__ == "__main__":
    main()
