#!/usr/bin/env python3
"""Independent exact audit of the even-Pell-anchor reduction.

This implementation deliberately does not import the curator's checker. It
uses direct Fibonacci and even-Pell recurrences, exhaustive orbit comparison
on a fixed box, and a separately organized crossing search. Large bounded
searches are falsification evidence, not a proof of unbounded nonexistence.
"""

from __future__ import annotations

import argparse
from math import gcd, isqrt


Pair = tuple[int, int]


def multiply(left: Pair, right: Pair) -> Pair:
    """Multiply two coefficient pairs in Z[sqrt(2)]."""

    a_value, b_value = left
    c_value, d_value = right
    return (
        a_value * c_value + 2 * b_value * d_value,
        a_value * d_value + b_value * c_value,
    )


def lambda_table(maximum: int) -> list[Pair]:
    """Return coefficient pairs of lambda^n for 0 <= n <= maximum."""

    values = [(1, 0)]
    for _ in range(maximum):
        values.append(multiply(values[-1], (1, 1)))
    return values


def lambda_pair(exponent: int, positive_table: list[Pair]) -> Pair:
    """Return the exact coefficient pair of lambda^exponent."""

    if exponent >= 0:
        return positive_table[exponent]
    index = -exponent
    companion, pell = positive_table[index]
    sign = -1 if index % 2 else 1
    return sign * companion, -sign * pell


def sign_of_surd(value: Pair) -> int:
    """Determine the sign of a+b*sqrt(2) using integer arithmetic."""

    rational, radical = value
    if rational == 0 and radical == 0:
        return 0
    if rational >= 0 and radical >= 0:
        return 1
    if rational <= 0 and radical <= 0:
        return -1
    comparison = rational * rational - 2 * radical * radical
    if rational > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def subtract(left: Pair, right: Pair) -> Pair:
    """Subtract coefficient pairs."""

    return left[0] - right[0], left[1] - right[1]


def forms(u_value: int, v_value: int) -> tuple[int, int, int]:
    """Return H, B, and Q."""

    return (
        u_value * u_value + v_value * v_value,
        v_value * v_value - u_value * u_value + 2 * u_value * v_value,
        v_value * v_value - 2 * u_value * v_value - u_value * u_value,
    )


def step(seed: Pair) -> Pair:
    """Apply S(u,v)=(v,u+2v)."""

    return seed[1], seed[0] + 2 * seed[1]


def fibonacci_pairs(maximum: int) -> list[Pair]:
    """Return (F_q,F_(q+1)) for 0 <= q <= maximum."""

    pairs = [(0, 1)]
    for _ in range(maximum):
        first, second = pairs[-1]
        pairs.append((second, first + second))
    return pairs


def verify_form_algebra() -> int:
    """Check the form identities and evolution independently on a box."""

    powers = lambda_table(24)
    comparisons = 0
    for u_value in range(-8, 9):
        for v_value in range(-8, 9):
            height, companion_form, invariant = forms(u_value, v_value)
            assert companion_form * companion_form + invariant * invariant == 2 * height * height
            current = (u_value, v_value)
            for time in range(9):
                height_now, companion_now, invariant_now = forms(*current)
                companion_even, pell_even = powers[2 * time]
                assert height_now == companion_even * height + pell_even * companion_form
                assert companion_now == 2 * pell_even * height + companion_even * companion_form
                assert invariant_now == (-1 if time % 2 else 1) * invariant
                current = step(current)
                comparisons += 1
    return comparisons


def gap_has_required_content(gap: Pair, time: int, powers: list[Pair]) -> tuple[bool, int]:
    """Test gap=g*lambda^(-(2t+1)) with an integer g>0."""

    unit = lambda_pair(-(2 * time + 1), powers)
    if gap[1] <= 0 or gap[1] % unit[1]:
        return False, 0
    content = gap[1] // unit[1]
    return gap[0] == content * unit[0] and content > 0, content


def verify_fibonacci_bridge_box(
    q_max: int = 120,
    d_max: int = 100,
    t_max: int = 30,
) -> int:
    """Compare direct hits with the asserted quadratic-unit equation."""

    powers = lambda_table(d_max + 2 * t_max + 2)
    fib = fibonacci_pairs(q_max)
    comparisons = 0
    hits: list[tuple[int, int, int, int]] = []
    for q_value in range(1, q_max + 1):
        f_value, x_value = fib[q_value]
        for sigma in (1, -1):
            seed = (x_value, x_value + sigma * f_value)
            height, companion_form, _ = forms(*seed)
            a_value = (sigma * f_value, x_value)
            a_square = multiply(a_value, a_value)
            current = seed
            for time in range(t_max + 1):
                current_height = forms(*current)[0]
                for d_value in range(2, d_max + 1, 2):
                    hit = current_height == powers[d_value + 2 * time][1]
                    gap = subtract(powers[d_value - 1], a_square)
                    gap_target, content = gap_has_required_content(gap, time, powers)
                    assert hit == gap_target
                    if hit:
                        hits.append((q_value, sigma, d_value, time))
                        difference_d = height - powers[d_value][1]
                        difference_e = powers[d_value][0] - companion_form
                        assert difference_d == content * powers[2 * time][1]
                        assert difference_e == content * powers[2 * time][0]
                        assert gcd(difference_d, difference_e) == content
                        assert q_value % 3 == 0
                        assert content % 2 == 1
                        assert (d_value + 2 * time) % 8 == 2
                    comparisons += 1
                current = step(current)
    assert hits == []
    return comparisons


def verify_local_filters() -> None:
    """Derive the parity, modulo-8, modulo-3, and Fibonacci filters."""

    powers = lambda_table(200)
    fib = fibonacci_pairs(120)
    for q_value in range(1, 121):
        f_value, x_value = fib[q_value]
        epsilon = -1 if q_value % 2 else 1
        for sigma in (1, -1):
            seed = (x_value, x_value + sigma * f_value)
            current = seed
            for time in range(17):
                height = forms(*current)[0]
                assert height % 2 == f_value % 2
                odd_pell = powers[2 * time + 1][1]
                assert (height - 2 * epsilon * odd_pell) % f_value == 0
                assert (height + epsilon * odd_pell) % x_value == 0
                if q_value % 3 == 0:
                    assert height % 8 == 2
                current = step(current)

    pell_even_residues = {index: powers[index][1] % 8 for index in range(0, 8, 2)}
    assert pell_even_residues == {0: 0, 2: 2, 4: 4, 6: 6}

    expected = {
        (0, 1): {0, 3},
        (0, -1): {0, 3},
        (1, 1): {1, 2},
        (1, -1): {1, 2},
        (2, 1): {2, 3},
        (2, -1): {0, 1},
        (3, 1): {0, 1},
        (3, -1): {2, 3},
    }
    representatives = {0: 12, 1: 3, 2: 6, 3: 9}
    for k_class, q_value in representatives.items():
        f_value, x_value = fib[q_value]
        for sigma in (1, -1):
            current = (x_value, x_value + sigma * f_value)
            allowed: set[int] = set()
            for time in range(4):
                if forms(*current)[0] % 3 == 2:
                    allowed.add(time)
                current = step(current)
            assert allowed == expected[(k_class, sigma)]


def append_even_pell(values: list[Pair]) -> None:
    """Append the next pair for lambda^(2r)."""

    companion, pell = values[-1]
    values.append((3 * companion + 4 * pell, 2 * companion + 3 * pell))


def first_even_power_above(values: list[Pair], target: Pair) -> int:
    """Find the first r>=1 for which lambda^(2r)>target."""

    while sign_of_surd(subtract(values[-1], target)) <= 0:
        append_even_pell(values)
    low = 1
    high = len(values) - 1
    while low < high:
        middle = (low + high) // 2
        if sign_of_surd(subtract(values[middle], target)) > 0:
            high = middle
        else:
            low = middle + 1
    return low


def verify_candidate_search(q_max: int) -> tuple[int, int, int, int]:
    """Run a separate exact one-crossing search through q_max."""

    fib = fibonacci_pairs(q_max)
    even_pell: list[Pair] = [(1, 0), (3, 2)]
    signed_seeds = 0
    positive_candidates = 0
    square_norms = 0
    exact_hits = 0

    for q_value in range(3, q_max + 1, 3):
        f_value, x_value = fib[q_value]
        for sigma in (1, -1):
            signed_seeds += 1
            u_value = x_value
            v_value = x_value + sigma * f_value
            height, companion_form, _ = forms(u_value, v_value)
            table_index = first_even_power_above(
                even_pell,
                (companion_form, height),
            )
            companion_d, pell_d = even_pell[table_index]
            difference_d = height - pell_d
            if difference_d <= 0:
                continue
            difference_e = companion_d - companion_form
            assert difference_e > 0
            norm = difference_e * difference_e - 2 * difference_d * difference_d
            assert norm > 0
            positive_candidates += 1

            root = isqrt(norm)
            if root * root != norm:
                continue
            square_norms += 1
            content = gcd(difference_d, difference_e)
            if root != content:
                continue
            exact_hits += 1
            normalized = (
                difference_e // content,
                difference_d // content,
            )
            normalized_index = first_even_power_above(
                even_pell,
                (normalized[0] - 1, normalized[1]),
            )
            assert even_pell[normalized_index] == normalized
            final_index = 2 * (table_index + normalized_index)
            evolved_height = (
                even_pell[normalized_index][0] * height
                + even_pell[normalized_index][1] * companion_form
            )
            while len(even_pell) <= table_index + normalized_index:
                append_even_pell(even_pell)
            assert evolved_height == even_pell[table_index + normalized_index][1]
            assert final_index % 8 == 2

    return signed_seeds, positive_candidates, square_norms, exact_hits


def verify_candidate_uniqueness_brutally(q_max: int = 600) -> int:
    """Enumerate every even anchor below H and count positive-norm rows."""

    fib = fibonacci_pairs(q_max)
    even_pell = [(1, 0), (3, 2)]
    maximum_height = 0
    for q_value in range(3, q_max + 1, 3):
        f_value, x_value = fib[q_value]
        for sigma in (1, -1):
            maximum_height = max(
                maximum_height,
                forms(x_value, x_value + sigma * f_value)[0],
            )
    while even_pell[-1][1] < maximum_height:
        append_even_pell(even_pell)

    checked = 0
    for q_value in range(3, q_max + 1, 3):
        f_value, x_value = fib[q_value]
        for sigma in (1, -1):
            height, companion_form, _ = forms(
                x_value,
                x_value + sigma * f_value,
            )
            candidates: list[int] = []
            for table_index, (companion_d, pell_d) in enumerate(even_pell[1:], start=1):
                if pell_d >= height:
                    break
                difference_d = height - pell_d
                difference_e = companion_d - companion_form
                norm = difference_e * difference_e - 2 * difference_d * difference_d
                if difference_e > 0 and norm > 0:
                    candidates.append(table_index)
                checked += 1
            assert len(candidates) <= 1
            if candidates:
                assert candidates[0] == first_even_power_above(
                    even_pell,
                    (companion_form, height),
                )
    return checked


def verify_specific_pseudohit() -> tuple[int, int, int, int, int]:
    """Independently reproduce the small-modulus survivor."""

    q_value = 357
    sigma = -1
    d_value = 390
    time = 70
    modulus = 2_042_040
    fib = fibonacci_pairs(q_value)
    f_value, x_value = fib[q_value]
    height, companion_form, _ = forms(x_value, x_value + sigma * f_value)
    powers = lambda_table(d_value + 2 * time)
    evolved_height = powers[2 * time][0] * height + powers[2 * time][1] * companion_form
    target = powers[d_value + 2 * time][1]
    assert evolved_height != target
    assert (evolved_height - target) % modulus == 0

    difference_d = height - powers[d_value][1]
    difference_e = powers[d_value][0] - companion_form
    norm = difference_e * difference_e - 2 * difference_d * difference_d
    assert difference_d > 0 and difference_e > 0 and norm > 0
    return q_value, sigma, d_value, time, modulus


def main() -> None:
    """Run all independent exact checks."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--q-max", type=int, default=5_000)
    args = parser.parse_args()

    form_comparisons = verify_form_algebra()
    bridge_comparisons = verify_fibonacci_bridge_box()
    verify_local_filters()
    anchor_comparisons = verify_candidate_uniqueness_brutally()
    counts = verify_candidate_search(args.q_max)
    pseudohit = verify_specific_pseudohit()
    print(
        "PASS: independent even-anchor audit; "
        f"form comparisons={form_comparisons}; "
        f"bridge comparisons={bridge_comparisons}; "
        f"anchor comparisons={anchor_comparisons}; "
        f"q<= {args.q_max}; signed seeds={counts[0]}; "
        f"positive candidates={counts[1]}; square norms={counts[2]}; "
        f"exact hits={counts[3]}; modular survivor={pseudohit}."
    )


if __name__ == "__main__":
    main()
