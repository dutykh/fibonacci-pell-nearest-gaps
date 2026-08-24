#!/usr/bin/env python3
"""Fixed exact regressions for the uniform local-clock generalisation.

The infinitude statements are proved by equidistribution in the accompanying
attempt. This checker performs no prime search. It verifies the sharp rank
obstruction at p=241 and two fixed simultaneous-prime nearest-window rows.
"""

from __future__ import annotations

from math import gcd


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


def pell_pair(index: int) -> Pair:
    """Return (C_index, P_index)."""

    return power((1, 1), index)


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


def fibonacci_mod(index: int, modulus: int) -> int:
    """Return F_index modulo modulus."""

    def doubled(value: int) -> Pair:
        if value == 0:
            return 0, 1
        first, second = doubled(value // 2)
        even = first * (2 * second - first) % modulus
        odd = (first * first + second * second) % modulus
        if value % 2:
            return odd, (even + odd) % modulus
        return even, odd

    return doubled(index)[0]


def companion_mod(index: int, modulus: int) -> int:
    """Return C_index modulo modulus."""

    return power_mod((1, 1), index, modulus)[0]


def power_mod(base: Pair, exponent: int, modulus: int) -> Pair:
    """Raise one coefficient pair modulo an integer."""

    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = (
                (result[0] * base[0] + 2 * result[1] * base[1]) % modulus,
                (result[0] * base[1] + result[1] * base[0]) % modulus,
            )
        base = (
            (base[0] * base[0] + 2 * base[1] * base[1]) % modulus,
            2 * base[0] * base[1] % modulus,
        )
        exponent >>= 1
    return result


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


def exact_nearest_window(q_value: int, m_value: int) -> None:
    """Verify D_q < lambda^m < sqrt(2) D_q exactly."""

    fib, fib_next = fibonacci_pair(q_value)
    root = (fib, fib_next)
    unit = pell_pair(m_value)
    upper = (2 * fib_next, fib)
    assert sign_surd(subtract(unit, root)) > 0
    assert sign_surd(subtract(upper, unit)) > 0


def cross_gcd(q_value: int, m_value: int) -> tuple[int, int]:
    """Return the direct cross-gcd and its equal-phase diagonal product."""

    fib, fib_next = fibonacci_pair(q_value)
    companion, pell = pell_pair(m_value)
    norm = 2 * fib_next * fib_next - fib * fib
    sign = -1 if m_value % 2 else 1
    direct = gcd(norm - sign, 2 * fib_next * pell - fib * companion)
    diagonal = gcd(fib_next, companion) * gcd(
        fib + fib_next,
        pell_pair(m_value + 1)[0],
    )
    return direct, diagonal


def check_rank_obstruction() -> None:
    """Certify z_F(241)=120 and z_C(241)=20 without a range search."""

    prime = 241
    assert prime % 8 == 1
    assert fibonacci_mod(120, prime) == 0
    assert fibonacci_mod(60, prime) == 34
    assert fibonacci_mod(40, prime) == 12
    assert fibonacci_mod(24, prime) == 96
    assert companion_mod(20, prime) == 0
    assert companion_mod(4, prime) == 17


def check_simultaneous_rows() -> None:
    """Verify two fixed multi-prime rows with exact integer arithmetic."""

    rows = (
        (5191, 2835, 23681, 7 * 17),
        (11543, 6303, 161, 7 * 23),
    )
    for q_value, m_value, expected, forced_product in rows:
        assert q_value % 3 != 0
        assert q_value % 2 == m_value % 2
        exact_nearest_window(q_value, m_value)
        direct, diagonal = cross_gcd(q_value, m_value)
        assert direct == diagonal == expected
        assert direct % forced_product == 0

    q_value, m_value = 5191, 2835
    assert q_value % 8 == 7
    assert q_value % 9 == 7
    assert m_value % 6 == 3
    assert m_value % 8 == 3

    q_value, m_value = 11543, 6303
    assert q_value % 24 == 23
    assert m_value % 66 == 33


def main() -> None:
    """Run the fixed exact regression checks."""

    check_rank_obstruction()
    check_simultaneous_rows()
    print(
        "PASS: uniform local-clock regressions; "
        "sharp obstruction=(p,z_F,z_C)=(241,120,20); "
        "simultaneous rows=(5191,2835,23681),(11543,6303,161)."
    )


if __name__ == "__main__":
    main()
