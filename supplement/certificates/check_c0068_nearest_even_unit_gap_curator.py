#!/usr/bin/env python3
"""Exact checks for C0068's nearest-even-unit gap criterion.

The proof in results/C0068-nearest-even-unit-gap.md is symbolic. This
standard-library checker:

* constructs the at most two nearest even-unit gaps exactly;
* tests their coefficient content and normalized norm through q=2000;
* reconstructs every hit and verifies C0065's orbit and determinant phase;
* compares with a direct bounded orbit-coefficient search through q=80.

The finite absence beyond q=4 is regression evidence only.

Run with:
    python3 -B scripts/check_c0068_nearest_even_unit_gap_curator.py
"""

from __future__ import annotations

from math import gcd, isqrt


Pair = tuple[int, int]


def multiply(left: Pair, right: Pair) -> Pair:
    """Multiply coordinate pairs in Z[sqrt(2)]."""
    a_value, b_value = left
    c_value, d_value = right
    return (
        a_value * c_value + 2 * b_value * d_value,
        a_value * d_value + b_value * c_value,
    )


def power(base: Pair, exponent: int) -> Pair:
    """Raise a quadratic integer to a nonnegative power."""
    assert exponent >= 0
    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def lambda_power(exponent: int) -> Pair:
    """Return the coordinates of (1 + sqrt(2))^exponent."""
    if exponent >= 0:
        return power((1, 1), exponent)
    return power((-1, 1), -exponent)


def norm(value: Pair) -> int:
    """Return the quadratic norm."""
    rational, radical = value
    return rational * rational - 2 * radical * radical


def surd_sign(value: Pair) -> int:
    """Return the exact sign of a + b*sqrt(2)."""
    rational, radical = value
    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if rational > 0 and radical > 0:
        return 1
    if rational < 0 and radical < 0:
        return -1
    rational_square = rational * rational
    radical_square = 2 * radical * radical
    if rational > 0:
        return 1 if rational_square > radical_square else -1
    return 1 if radical_square > rational_square else -1


def less(left: Pair, right: Pair) -> bool:
    """Test the principal-embedding inequality left < right exactly."""
    difference = (right[0] - left[0], right[1] - left[1])
    return surd_sign(difference) > 0


def fibonacci_pair(index: int) -> Pair:
    """Return (F_index, F_(index + 1))."""
    first, second = 0, 1
    for _ in range(index):
        first, second = second, first + second
    return first, second


def selector_data(index: int) -> tuple[int, int, int, Pair, Pair]:
    """Return f, x, k, A^2, and B^2."""
    f_value, x_value = fibonacci_pair(index)
    candidate = 2 * x_value * x_value - f_value * f_value
    rational = f_value * f_value + 2 * x_value * x_value
    radical = 2 * f_value * x_value
    return (
        f_value,
        x_value,
        candidate,
        (rational, radical),
        (rational, -radical),
    )


def first_larger_even_power(
    square: Pair,
    even_powers: list[Pair],
) -> tuple[int, Pair]:
    """Return the least lambda^(2j) exceeding square."""
    while not less(square, even_powers[-1]):
        even_powers.append(multiply(even_powers[-1], (3, 2)))

    low = 0
    high = len(even_powers) - 1
    while low < high:
        middle = (low + high) // 2
        if less(square, even_powers[middle]):
            high = middle
        else:
            low = middle + 1
    return 2 * low, even_powers[low]


def nearest_gap(
    square: Pair,
    even_powers: list[Pair],
) -> tuple[int, Pair, int] | None:
    """Return n, lambda^n-square, and its content when the gap exists."""
    index, unit = first_larger_even_power(square, even_powers)
    twice_square = (2 * square[0], 2 * square[1])
    if not less(unit, twice_square):
        return None
    difference = (unit[0] - square[0], unit[1] - square[1])
    assert surd_sign(difference) > 0
    content = gcd(abs(difference[0]), abs(difference[1]))
    assert content > 0
    return index, difference, content


def find_unit_exponent(unit: Pair, index_bound: int) -> int:
    """Recover the unique odd exponent of a positive norm-minus-one unit."""
    assert norm(unit) == -1
    assert surd_sign(unit) > 0
    first_odd = -index_bound if index_bound % 2 else -index_bound + 1
    for exponent in range(first_odd, index_bound + 1, 2):
        if lambda_power(exponent) == unit:
            return exponent
    raise AssertionError("unit exponent outside proved finite recovery bound")


def rectangle_delta(
    candidate: int,
    companion: int,
    other: int,
) -> int:
    """Reconstruct the signed C0060 rectangle and return its determinant."""
    first = (companion - other) // 2
    second = (companion + other) // 2
    assert first > 0 and second > 0
    assert first * second == (candidate * candidate - 1) // 4

    valid: list[int] = []
    for delta in (1, -1):
        x_half = (candidate - delta) // 2
        y_half = (candidate + delta) // 2
        a_coord = gcd(first, x_half)
        d_coord = gcd(first, y_half)
        s_coord = x_half // a_coord
        t_coord = y_half // d_coord
        if (
            a_coord * d_coord == first
            and s_coord * t_coord == second
            and d_coord * t_coord - a_coord * s_coord == delta
        ):
            h_square = s_coord * s_coord + d_coord * d_coord
            j_square_numerator = t_coord * t_coord + a_coord * a_coord
            hypotenuse = isqrt(h_square)
            pell = isqrt(j_square_numerator // 2)
            if (
                hypotenuse * hypotenuse == h_square
                and 2 * pell * pell == j_square_numerator
                and gcd(a_coord * hypotenuse, d_coord * pell) == 1
            ):
                valid.append(delta)
    assert len(valid) == 1
    return valid[0]


def reconstruct_hit(
    q_value: int,
    branch: str,
    gap: tuple[int, Pair, int],
) -> tuple[int, str, int, int, int, int]:
    """Recover and verify M, e, W, and delta from one norm-minus-one gap."""
    index, difference, content = gap
    assert norm(difference) == -(content * content)
    quotient = (
        difference[0] // content,
        difference[1] // content,
    )
    exponent_j = find_unit_exponent(quotient, 2 * index + 3)
    pell_index = index - exponent_j
    assert pell_index > 0 and pell_index % 2
    companion, pell = lambda_power(pell_index)

    f_value, x_value, candidate, a_square, _ = selector_data(q_value)
    orbit_exponent = -exponent_j if branch == "plus" else exponent_j
    alpha = multiply(a_square, lambda_power(orbit_exponent))
    other, pell_coefficient = alpha
    assert pell_coefficient == pell
    assert other * other - 2 * pell * pell == -(candidate * candidate)
    assert gcd(candidate, other) == 1
    expected_other = (
        companion - content if branch == "plus" else content - companion
    )
    assert other == expected_other
    assert (other > 0) == (branch == "plus")

    delta = rectangle_delta(candidate, companion, other)
    expected_delta = (-1) ** (index // 2)
    if branch == "minus":
        expected_delta = -expected_delta
    assert delta == expected_delta
    assert delta == (-1) ** ((pell_index - orbit_exponent) // 2)

    scalar = x_value * other - f_value * pell
    assert scalar % candidate == 0
    return (
        q_value,
        branch,
        index,
        pell_index,
        orbit_exponent,
        delta,
    )


def direct_orbit_search(
    q_bound: int = 80,
    exponent_bound: int = 159,
) -> tuple[tuple[int, int, int], ...]:
    """Direct bounded C0065 coefficient search for comparison."""
    pell_to_index = {
        lambda_power(index)[1]: index
        for index in range(1, 401, 2)
    }
    hits: list[tuple[int, int, int]] = []
    for q_value in range(2, q_bound + 1):
        if q_value % 3 == 0:
            continue
        f_value, x_value, _, _, _ = selector_data(q_value)
        for exponent in range(-exponent_bound, exponent_bound + 1, 2):
            companion_e, pell_e = lambda_power(exponent)
            pell_coefficient = (
                2 * f_value * x_value * companion_e
                + (f_value * f_value + 2 * x_value * x_value) * pell_e
            )
            pell_index = pell_to_index.get(pell_coefficient)
            if pell_index is not None:
                hits.append((q_value, exponent, pell_index))
    return tuple(hits)


def main(q_bound: int = 2000) -> None:
    """Run the q-only criterion and all fixed regressions."""
    even_powers = [(1, 0)]
    candidate_counts = {"plus": 0, "minus": 0}
    hits: list[tuple[int, str, int, int, int, int]] = []

    for q_value in range(2, q_bound + 1):
        if q_value % 3 == 0:
            continue
        _, _, candidate, plus_square, minus_square = selector_data(q_value)
        assert candidate > 1 and candidate % 2
        assert norm(plus_square) == candidate * candidate
        assert norm(minus_square) == candidate * candidate

        for branch, square in (
            ("plus", plus_square),
            ("minus", minus_square),
        ):
            gap = nearest_gap(square, even_powers)
            if gap is None:
                continue
            candidate_counts[branch] += 1
            index, difference, content = gap
            assert index > 0 and index % 2 == 0
            assert norm(difference) < 0
            assert norm(difference) % (content * content) == 0
            if norm(difference) == -(content * content):
                hits.append(reconstruct_hit(q_value, branch, gap))

    assert hits == [
        (2, "minus", 2, 3, -1, 1),
        (4, "plus", 6, 5, -1, -1),
    ]

    direct_hits = direct_orbit_search()
    assert direct_hits == ((2, -1, 3), (4, -1, 5))

    q_two_minus = nearest_gap(selector_data(2)[4], even_powers)
    q_four_plus = nearest_gap(selector_data(4)[3], even_powers)
    assert q_two_minus == (2, (-6, 6), 6)
    assert q_four_plus == (6, (40, 40), 40)

    print(
        "PASS: C0068 nearest-even-unit gap; "
        f"q<={q_bound}; plus_candidates={candidate_counts['plus']}; "
        f"minus_candidates={candidate_counts['minus']}; "
        f"normalized_norm_minus_one_hits={hits}; "
        f"direct_orbit_q<=80_e<=159={direct_hits}. "
        "Finite absence beyond q=4 is evidence only."
    )


if __name__ == "__main__":
    main()
