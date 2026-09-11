#!/usr/bin/env python3
"""Independent integer audit of the general-seed orbit/gap theorem.

No project module or originating checker is imported.  The finite box is a
regression test only; the companion audit supplies the universal proof.
"""

from __future__ import annotations

from math import gcd


Pair = tuple[int, int]
Point = tuple[int, int]

MAX_COORDINATE = 60
MAX_D = 49
MAX_TIME = 12


def pell_values(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(2 * values[-1] + values[-2])
    return values


PELL = pell_values(200)


def companion(index: int) -> int:
    previous = 1 if index == 0 else PELL[index - 1]
    return PELL[index] + previous


def unit(exponent: int) -> Pair:
    """Return coefficients of lambda**exponent by Pell recurrence."""

    if exponent >= 0:
        return companion(exponent), PELL[exponent]
    index = -exponent
    parity = -1 if index % 2 else 1
    return parity * companion(index), -parity * PELL[index]


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
    return value[0] ** 2 - 2 * value[1] ** 2


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


def step(point: Point) -> Point:
    u, v = point
    return v, u + 2 * v


def forms(point: Point) -> tuple[int, int, int]:
    u, v = point
    return (
        u * u + v * v,
        v * v - u * u + 2 * u * v,
        v * v - 2 * u * v - u * u,
    )


def attached(point: Point) -> Pair:
    u, v = point
    return v - u, u


def inverse_attached(value: Pair) -> Point:
    a, b = value
    return b, a + b


def positive_integer_multiple(dividend: Pair, divisor: Pair) -> int | None:
    assert divisor[0] != 0 and divisor[1] != 0
    if dividend[0] % divisor[0] or dividend[1] % divisor[1]:
        return None
    first = dividend[0] // divisor[0]
    second = dividend[1] // divisor[1]
    return first if first == second and first > 0 else None


def check_unit_signs() -> None:
    for exponent in range(-150, 151):
        value = unit(exponent)
        expected_norm = -1 if exponent % 2 else 1
        assert norm(value) == expected_norm
        assert multiply(value, unit(-exponent)) == (1, 0)
        expected_conjugate = scale(expected_norm, unit(-exponent))
        assert conjugate(value) == expected_conjugate
        assert content(value) == 1


def check_box() -> tuple[int, list[tuple[int, int, int, int, int]], list[Point]]:
    comparisons = 0
    noncanonical_hits: list[tuple[int, int, int, int, int]] = []
    canonical_seen: list[Point] = []

    canonical_expected: dict[Point, int] = {}
    rank = 1
    while PELL[rank + 1] <= MAX_COORDINATE:
        canonical_expected[(PELL[rank], PELL[rank + 1])] = rank
        rank += 1

    for u in range(1, MAX_COORDINATE + 1):
        for v in range(u, MAX_COORDINATE + 1):
            origin = (u, v)
            a_value = attached(origin)
            height_zero, bilinear_zero, qform = forms(origin)

            # The map is a literal bijection between the ordered seed cone
            # and coefficient pairs (a,b) with a>=0 and b>0.
            assert a_value[0] >= 0 and a_value[1] > 0
            assert inverse_attached(a_value) == origin
            assert norm(a_value) == qform
            assert qform != 0
            assert content(a_value) == gcd(u, v)
            assert attached(step(origin)) == multiply(unit(1), a_value)

            gamma = (u - v, u)
            assert gamma == scale(-1, conjugate(a_value))
            assert multiply(unit(-1), square(gamma)) == (
                -bilinear_zero,
                height_zero,
            )
            # The strengthened proof uses this parity-free identity.
            assert multiply(unit(1), square(a_value)) == (
                bilinear_zero,
                height_zero,
            )

            # This checks the universal comparison independently of any hit.
            absolute_gamma = gamma if sign_surd(gamma) >= 0 else scale(-1, gamma)
            assert sign_surd(subtract(multiply(unit(1), a_value), absolute_gamma)) > 0
            assert sign_surd(subtract(a_value, (0, 1))) >= 0

            if u == v:
                assert a_value == (0, u)
                assert gamma == a_value
                assert qform == -2 * u * u
                assert qform * qform > 1

            is_canonical = qform * qform == 1
            assert is_canonical == (origin in canonical_expected)
            if is_canonical:
                canonical_seen.append(origin)
                r_index = canonical_expected[origin]
                assert a_value == unit(r_index)
                point = origin
                for crossing_time in range(MAX_TIME + 1):
                    assert point == (
                        PELL[r_index + crossing_time],
                        PELL[r_index + crossing_time + 1],
                    )
                    assert forms(point)[0] == PELL[2 * (r_index + crossing_time) + 1]
                    for d_index in range(1, MAX_D + 1):
                        hit = forms(point)[0] == PELL[d_index + 2 * crossing_time]
                        assert hit == (d_index == 2 * r_index + 1)
                        if hit:
                            assert subtract(unit(d_index - 1), square(a_value)) == (0, 0)
                    point = step(point)
                continue

            assert qform * qform > 1
            point = origin
            for crossing_time in range(MAX_TIME + 1):
                height_now, bilinear_now, qform_now = forms(point)
                assert attached(point) == multiply(unit(crossing_time), a_value)
                assert qform_now == (-1) ** crossing_time * qform
                assert qform_now * qform_now == qform * qform
                assert point[0] > 0 and point[0] <= point[1]
                assert bilinear_now > 0

                c_even = companion(2 * crossing_time)
                p_even = PELL[2 * crossing_time]
                assert height_now == c_even * height_zero + p_even * bilinear_zero
                assert bilinear_now == 2 * p_even * height_zero + c_even * bilinear_zero

                for d_index in range(1, MAX_D + 1):
                    comparisons += 1
                    target_index = d_index + 2 * crossing_time
                    hit = height_now == PELL[target_index]

                    difference_d = height_zero - PELL[d_index]
                    difference_e = companion(d_index) - bilinear_zero
                    # New all-parity unconditional identity (15).
                    assert subtract(
                        unit(d_index), multiply(unit(1), square(a_value))
                    ) == (
                        difference_e,
                        -difference_d,
                    )

                    gap = subtract(unit(d_index - 1), square(a_value))
                    odd_unit = unit(-(2 * crossing_time + 1))
                    multiplier = positive_integer_multiple(gap, odd_unit)
                    # Exact fixed-anchor converse across the entire box.
                    assert (multiplier is not None) == hit

                    if d_index == 1:
                        assert multiplier is None
                        assert sign_surd(gap) < 0

                    if not hit:
                        continue

                    assert multiplier is not None
                    g_value = multiplier
                    noncanonical_hits.append(
                        (u, v, d_index, crossing_time, qform)
                    )
                    parity_term = -1 if d_index % 2 else 1
                    assert bilinear_now * bilinear_now == (
                        companion(target_index) ** 2
                        - (qform * qform + parity_term)
                    )
                    assert 0 < bilinear_now < companion(target_index)
                    if crossing_time == 0:
                        assert difference_d == 0
                        assert difference_e > 0
                    else:
                        assert c_even * difference_d == p_even * difference_e
                        assert difference_d > 0 and difference_e > 0
                    assert difference_d == g_value * p_even
                    assert difference_e == g_value * c_even
                    assert gcd(difference_d, difference_e) == g_value
                    assert gap == scale(g_value, odd_unit)
                    assert content(gap) == g_value
                    assert norm(gap) == -g_value * g_value
                    assert sign_surd(gap) > 0
                    assert sign_surd(subtract(scale(2, square(a_value)), unit(d_index - 1))) > 0
                    assert d_index >= 2

                    # Direct-time form of the gap and the parity-uniform
                    # proof of its strict upper window.
                    a_now = attached(point)
                    direct_gap = subtract(unit(target_index - 1), square(a_now))
                    assert direct_gap == scale(g_value, unit(-1))
                    assert g_value == companion(target_index) - bilinear_now
                    direct_margin = subtract(
                        multiply(unit(1), square(a_now)), (g_value, 0)
                    )
                    correction_sign = 1 if target_index % 2 else -1
                    correction = add(
                        (2 * bilinear_now, 0),
                        scale(correction_sign, unit(-target_index)),
                    )
                    assert direct_margin == correction
                    assert sign_surd(direct_margin) > 0

                point = step(point)

    assert set(canonical_seen) == set(canonical_expected)
    return comparisons, noncanonical_hits, canonical_seen


def main() -> None:
    check_unit_signs()
    comparisons, hits, canonical = check_box()
    assert hits == [
        (1, 1, 2, 0, -2),
        (13, 47, 10, 0, 818),
        (16, 27, 9, 0, -391),
        (23, 43, 10, 0, -658),
    ]
    print(
        "PASS: independent general-seed orbit/gap audit; "
        f"coordinates<={MAX_COORDINATE}, all_d<={MAX_D}, "
        f"t<={MAX_TIME}, comparisons={comparisons}, "
        f"canonical_seeds={canonical}, noncanonical_hits={hits}; "
        "bijection, unit boundary, parity-uniform companion comparison, "
        "all-d exact converse, d=1 exclusion, u=v face, and direct-time "
        "strict factor-two bound verified."
    )


if __name__ == "__main__":
    main()
