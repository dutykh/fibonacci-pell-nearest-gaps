#!/usr/bin/env python3
"""Exact support for the even-Pell-anchor orbit investigation.

The accompanying proof note gives the unbounded implications. This program
checks the algebra, the parity and small-modulus consequences, a fixed brute
box, and a larger exact candidate search. The bounded searches are
falsification evidence only and are not used as a proof of nonexistence.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left
from math import gcd, isqrt


Pair = tuple[int, int]


def multiply(left: Pair, right: Pair) -> Pair:
    """Multiply coefficient pairs in Z[sqrt(2)]."""

    a_value, b_value = left
    c_value, d_value = right
    return (
        a_value * c_value + 2 * b_value * d_value,
        a_value * d_value + b_value * c_value,
    )


def power(base: Pair, exponent: int) -> Pair:
    """Raise a quadratic integer to a nonnegative power."""

    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent >>= 1
    return result


def lambda_power(exponent: int) -> Pair:
    """Return the coefficients of (1+sqrt(2))**exponent."""

    if exponent >= 0:
        return power((1, 1), exponent)
    positive = power((1, 1), -exponent)
    sign = -1 if exponent % 2 else 1
    return sign * positive[0], -sign * positive[1]


def fibonacci_pair(index: int) -> Pair:
    """Return (F_index, F_(index+1)) by fast doubling."""

    if index == 0:
        return 0, 1
    first, second = fibonacci_pair(index // 2)
    even = first * (2 * second - first)
    odd = first * first + second * second
    if index % 2:
        return odd, even + odd
    return even, odd


def subtract(left: Pair, right: Pair) -> Pair:
    """Subtract coefficient pairs."""

    return left[0] - right[0], left[1] - right[1]


def sign_surd(value: Pair) -> int:
    """Return the exact sign of a+b*sqrt(2)."""

    rational, radical = value
    if rational == 0 and radical == 0:
        return 0
    if rational >= 0 and radical >= 0:
        return 1
    if rational <= 0 and radical <= 0:
        return -1
    norm = rational * rational - 2 * radical * radical
    if rational > 0:
        return 1 if norm > 0 else -1
    return -1 if norm > 0 else 1


def forms(u_value: int, v_value: int) -> tuple[int, int, int]:
    """Return the H, B, and Q forms."""

    return (
        u_value * u_value + v_value * v_value,
        v_value * v_value - u_value * u_value + 2 * u_value * v_value,
        v_value * v_value - 2 * u_value * v_value - u_value * u_value,
    )


def iterate_s(u_value: int, v_value: int, steps: int) -> Pair:
    """Apply S(u,v)=(v,u+2v) exactly."""

    for _ in range(steps):
        u_value, v_value = v_value, u_value + 2 * v_value
    return u_value, v_value


def seed_data(q_value: int, sigma: int) -> tuple[int, int, int, int, int]:
    """Return f, x, H, B, Q for one two-sign Fibonacci seed."""

    fib, fib_next = fibonacci_pair(q_value)
    height, companion_form, invariant = forms(
        fib_next,
        fib_next + sigma * fib,
    )
    return fib, fib_next, height, companion_form, invariant


def verify_unconditional_algebra(q_max: int = 50, d_max: int = 80) -> None:
    """Check all identities in the even-anchor bridge on a fixed box."""

    for q_value in range(1, q_max + 1):
        for sigma in (1, -1):
            fib, fib_next, height, companion_form, invariant = seed_data(
                q_value,
                sigma,
            )
            assert invariant == fib * fib - 2 * fib_next * fib_next
            assert companion_form * companion_form + invariant * invariant == 2 * height * height

            a_value = (sigma * fib, fib_next)
            gamma = (-sigma * fib, fib_next)
            assert multiply(lambda_power(-1), multiply(gamma, gamma)) == (
                -companion_form,
                height,
            )

            for d_value in range(2, d_max + 1, 2):
                companion_d, pell_d = lambda_power(d_value)
                difference_d = height - pell_d
                difference_e = companion_d - companion_form
                gap = subtract(lambda_power(d_value - 1), multiply(a_value, a_value))
                bridge = multiply(lambda_power(-1), (difference_e, -difference_d))
                assert gap == bridge

                left = (
                    multiply(lambda_power(d_value - 1), multiply(gamma, gamma))[0] + 1,
                    multiply(lambda_power(d_value - 1), multiply(gamma, gamma))[1],
                )
                right = multiply(
                    lambda_power(d_value),
                    (difference_e, difference_d),
                )
                assert left == right


def verify_parity_and_residue_package() -> None:
    """Verify the mod-2, mod-8, and mod-3 necessary residue package."""

    expected_mod_three = {
        (0, 1): {0, 3},
        (0, -1): {0, 3},
        (1, 1): {1, 2},
        (1, -1): {1, 2},
        (2, 1): {2, 3},
        (2, -1): {0, 1},
        (3, 1): {0, 1},
        (3, -1): {2, 3},
    }

    for q_value in range(1, 97):
        fib, fib_next, _, _, _ = seed_data(q_value, 1)
        for sigma in (1, -1):
            _, _, height, companion_form, _ = seed_data(q_value, sigma)
            assert height % 2 == fib % 2
            epsilon = -1 if q_value % 2 else 1
            for time in range(13):
                companion_even, pell_even = lambda_power(2 * time)
                evolved_height = (
                    companion_even * height + pell_even * companion_form
                )
                pell_odd = lambda_power(2 * time + 1)[1]
                assert (evolved_height - 2 * epsilon * pell_odd) % fib == 0
                assert (evolved_height + epsilon * pell_odd) % fib_next == 0

    for q_value in range(3, 49, 3):
        for sigma in (1, -1):
            _, _, height, companion_form, _ = seed_data(q_value, sigma)
            assert height % 8 == 2
            observed: set[int] = set()
            for time in range(4):
                companion_even, pell_even = lambda_power(2 * time)
                evolved_height = companion_even * height + pell_even * companion_form
                assert evolved_height % 8 == 2
                if evolved_height % 3 == 2:
                    observed.add(time)
            assert observed == expected_mod_three[(q_value // 3 % 4, sigma)]

    even_pell_mod_eight = {
        index: lambda_power(index)[1] % 8 for index in range(0, 8, 2)
    }
    assert even_pell_mod_eight == {0: 0, 2: 2, 4: 4, 6: 6}


def verify_candidate_uniqueness_box(q_max: int = 300) -> None:
    """Compare the one-candidate proof with exhaustive anchors in a box."""

    seeds = fibonacci_seeds(q_max)
    pell = even_pell_table(max(seed[2] for seed in seeds))
    for _, _, height, companion_form, _ in seeds:
        candidates: list[int] = []
        for table_index, (companion_d, pell_d) in enumerate(pell[1:], start=1):
            if pell_d >= height:
                break
            difference_d = height - pell_d
            difference_e = companion_d - companion_form
            norm = difference_e * difference_e - 2 * difference_d * difference_d
            if difference_e > 0 and norm > 0:
                candidates.append(table_index)
        assert len(candidates) <= 1
        if candidates:
            assert candidates[0] == first_power_above(
                pell,
                (companion_form, height),
            )


def verify_brute_box(
    q_max: int = 80,
    d_max: int = 120,
    t_max: int = 40,
) -> int:
    """Search a fixed three-dimensional box with direct orbit iteration."""

    hits: list[tuple[int, int, int, int]] = []
    pell_values = {
        index: lambda_power(index)[1]
        for index in range(2, d_max + 2 * t_max + 1, 2)
    }
    for q_value in range(1, q_max + 1):
        for sigma in (1, -1):
            fib, fib_next = fibonacci_pair(q_value)
            current = (fib_next, fib_next + sigma * fib)
            for time in range(t_max + 1):
                height = forms(*current)[0]
                for d_value in range(2, d_max + 1, 2):
                    if height == pell_values[d_value + 2 * time]:
                        hits.append((q_value, sigma, d_value, time))
                current = iterate_s(*current, 1)
    assert hits == []
    return q_max * 2 * (d_max // 2) * (t_max + 1)


def fibonacci_seeds(q_max: int) -> list[tuple[int, int, int, int, int]]:
    """Build all q=0 mod 3 two-sign seed data through q_max."""

    seeds: list[tuple[int, int, int, int, int]] = []
    fib, fib_next = 0, 1
    for q_value in range(1, q_max + 1):
        fib, fib_next = fib_next, fib + fib_next
        if q_value % 3:
            continue
        for sigma in (1, -1):
            height, companion_form, invariant = forms(
                fib_next,
                fib_next + sigma * fib,
            )
            seeds.append((q_value, sigma, height, companion_form, invariant))
    return seeds


def even_pell_table(max_height: int) -> list[Pair]:
    """Return lambda^(2r) until its Pell coordinate exceeds max_height."""

    values: list[Pair] = [(1, 0)]
    current = (1, 0)
    while current[1] <= max_height:
        current = multiply(current, (3, 2))
        values.append(current)
    return values


def first_power_above(values: list[Pair], target: Pair) -> int:
    """Find the first table index whose principal embedding exceeds target."""

    low = 1
    high = len(values) - 1
    while low < high:
        middle = (low + high) // 2
        if sign_surd(subtract(values[middle], target)) > 0:
            high = middle
        else:
            low = middle + 1
    return low


def verify_large_candidate_search(q_max: int) -> tuple[int, int, int, int]:
    """Search the unique positive-norm even anchor for each q and sign."""

    seeds = fibonacci_seeds(q_max)
    max_height = max(seed[2] for seed in seeds)
    pell = even_pell_table(max_height)
    pell_coordinates = [value[1] for value in pell]

    direct_hits: list[tuple[int, int, int]] = []
    positive_norm_candidates = 0
    square_norm_candidates: list[tuple[int, int, int]] = []
    exact_hits: list[tuple[int, int, int, int]] = []

    for q_value, sigma, height, companion_form, invariant in seeds:
        direct_index = bisect_left(pell_coordinates, height)
        if direct_index < len(pell) and pell[direct_index][1] == height:
            direct_hits.append((q_value, sigma, 2 * direct_index))

        table_index = first_power_above(pell, (companion_form, height))
        companion_d, pell_d = pell[table_index]
        if pell_d >= height:
            continue

        difference_d = height - pell_d
        difference_e = companion_d - companion_form
        norm = difference_e * difference_e - 2 * difference_d * difference_d
        if norm <= 0:
            continue
        positive_norm_candidates += 1

        root = isqrt(norm)
        if root * root != norm:
            continue
        square_norm_candidates.append((q_value, sigma, 2 * table_index))

        content = gcd(difference_d, difference_e)
        if root != content:
            continue
        normalized = (difference_e // content, difference_d // content)
        normalized_index = bisect_left(pell_coordinates, normalized[1])
        assert pell[normalized_index] == normalized
        time = normalized_index
        assert (2 * table_index + 2 * time) % 8 == 2
        exact_hits.append((q_value, sigma, 2 * table_index, time))

        height_at_time = (
            pell[time][0] * height + pell[time][1] * companion_form
        )
        assert height_at_time == lambda_power(2 * table_index + 2 * time)[1]
        assert invariant * invariant + 1 == content * (
            2 * lambda_power(2 * table_index + 2 * time)[0] - content
        )

    assert direct_hits == []
    assert square_norm_candidates == []
    assert exact_hits == []
    return len(seeds), positive_norm_candidates, len(square_norm_candidates), len(exact_hits)


def verify_modular_route_survivor() -> tuple[int, int, int, int, int]:
    """Record a positive-slope pseudohit surviving several small moduli."""

    q_value, sigma, d_value, time = 357, -1, 390, 70
    modulus = 2_042_040  # 8 * 3 * 5 * 7 * 11 * 13 * 17
    _, _, height, companion_form, _ = seed_data(q_value, sigma)
    companion_even, pell_even = lambda_power(2 * time)
    evolved_height = companion_even * height + pell_even * companion_form
    target = lambda_power(d_value + 2 * time)[1]
    assert evolved_height != target
    assert (evolved_height - target) % modulus == 0

    companion_d, pell_d = lambda_power(d_value)
    difference_d = height - pell_d
    difference_e = companion_d - companion_form
    assert difference_d > 0 and difference_e > 0
    assert difference_e * difference_e - 2 * difference_d * difference_d > 0
    return q_value, sigma, d_value, time, modulus


def main() -> None:
    """Run all exact regressions."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--q-max", type=int, default=5_000)
    args = parser.parse_args()

    verify_unconditional_algebra()
    verify_parity_and_residue_package()
    verify_candidate_uniqueness_box()
    comparisons = verify_brute_box()
    seeds, positive, squares, hits = verify_large_candidate_search(args.q_max)
    pseudohit = verify_modular_route_survivor()
    print(
        "PASS: even-Pell-anchor exact regressions; "
        f"brute comparisons={comparisons}; q<= {args.q_max}; "
        f"q=0 mod 3 signed seeds={seeds}; positive-norm candidates={positive}; "
        f"square norms={squares}; exact hits={hits}; "
        f"small-modulus survivor={pseudohit}."
    )


if __name__ == "__main__":
    main()
