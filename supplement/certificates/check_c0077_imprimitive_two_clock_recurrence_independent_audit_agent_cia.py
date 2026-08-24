#!/usr/bin/env python3
"""Independent exact checks for the imprimitive two-clock obstruction."""

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


def subtract(left: Pair, right: Pair) -> Pair:
    """Subtract two coefficient pairs."""

    return left[0] - right[0], left[1] - right[1]


def pell_pair(index: int) -> Pair:
    """Return (C_index, P_index), including negative indices."""

    if index < 0:
        companion, pell = pell_pair(-index)
        norm_sign = -1 if (-index) % 2 else 1
        return norm_sign * companion, -norm_sign * pell
    result = (1, 0)
    base = (1, 1)
    exponent = index
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent >>= 1
    return result


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


def sign_surd(value: Pair) -> int:
    """Return the exact sign of a+b*sqrt(2)."""

    rational, radical = value
    if rational >= 0 and radical >= 0:
        return 1
    if rational <= 0 and radical <= 0:
        return -1
    norm = rational * rational - 2 * radical * radical
    if rational > 0:
        return 1 if norm > 0 else -1
    return -1 if norm > 0 else 1


def valuation(value: int, prime: int) -> int:
    """Return the prime-adic valuation of a nonzero integer."""

    value = abs(value)
    assert value
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def exact_window(q_value: int, sigma: int, m_value: int) -> None:
    """Check that m is the least Pell exponent above D and is in its window."""

    fib, fib_next = fibonacci_pair(q_value)
    root = sigma * fib, fib_next
    unit = pell_pair(m_value)
    previous = pell_pair(m_value - 1)
    upper = 2 * fib_next, sigma * fib
    assert sign_surd(subtract(unit, root)) > 0
    assert sign_surd(subtract(root, previous)) >= 0
    assert sign_surd(subtract(upper, unit)) > 0


def aligned_bins(
    q_value: int,
    sigma: int,
    m_value: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Return C0077's two exact Fibonacci-companion index pairs."""

    epsilon = -1 if q_value % 2 else 1
    sign = -1 if m_value % 2 else 1
    if sign == epsilon:
        shift = 1 if sigma == 1 else -1
        return (q_value + 1, m_value), (q_value + 2, m_value + shift)
    if sigma == 1:
        return (q_value - 1, m_value - 1), (q_value + 4, m_value + 2)
    return (q_value - 1, m_value + 1), (q_value + 4, m_value - 2)


def exact_cross_gcd(q_value: int, sigma: int, m_value: int) -> tuple[int, int]:
    """Return the direct cross-gcd and independent diagonal product."""

    fib, fib_next = fibonacci_pair(q_value)
    companion, pell = pell_pair(m_value)
    norm = 2 * fib_next * fib_next - fib * fib
    sign = -1 if m_value % 2 else 1
    cross = 2 * fib_next * pell - sigma * fib * companion
    direct = gcd(norm - sign, cross)
    first, second = aligned_bins(q_value, sigma, m_value)
    diagonal = gcd(
        fibonacci_pair(first[0])[0],
        abs(pell_pair(first[1])[0]),
    ) * gcd(
        fibonacci_pair(second[0])[0],
        abs(pell_pair(second[1])[0]),
    )
    return direct, diagonal


def normalized_gap_norm(q_value: int, sigma: int, m_value: int) -> int:
    """Return N(lambda^(2m)-D_sigma^2) divided by content squared."""

    fib, fib_next = fibonacci_pair(q_value)
    root = sigma * fib, fib_next
    gap = subtract(
        multiply(pell_pair(m_value), pell_pair(m_value)),
        multiply(root, root),
    )
    content = gcd(abs(gap[0]), abs(gap[1]))
    norm = gap[0] * gap[0] - 2 * gap[1] * gap[1]
    assert norm % (content * content) == 0
    return norm // (content * content)


def check_rank_lifts() -> None:
    """Verify exact base ranks and several arbitrarily extendable lifts."""

    assert fibonacci_pair(8)[0] == 21
    assert pell_pair(3)[0] == 7
    assert fibonacci_pair(9)[0] == 34
    assert pell_pair(4)[0] == 17
    for prime, fib_rank, comp_rank, bound in (
        (7, 8, 3, 5),
        (17, 9, 4, 4),
    ):
        for exponent in range(1, bound + 1):
            multiplier = prime ** (exponent - 1)
            assert valuation(
                fibonacci_pair(fib_rank * multiplier)[0],
                prime,
            ) == exponent
            assert valuation(
                pell_pair(comp_rank * multiplier)[0],
                prime,
            ) == exponent


def check_progression_algebra() -> None:
    """Check parity and Markoff-index admissibility in both progressions."""

    for prime, fib_rank, comp_rank, q_factor in (
        (7, 8, 3, 3),
        (17, 9, 4, 2),
    ):
        for exponent in range(1, 6):
            fib_multiple = fib_rank * prime ** (exponent - 1)
            comp_multiple = comp_rank * prime ** (exponent - 1)
            for step in range(20):
                q_value = fib_multiple - 1 + q_factor * fib_multiple * step
                m_value = comp_multiple * (2 * step + 1)
                assert (q_value + 1) % fib_multiple == 0
                assert m_value % (2 * comp_multiple) == comp_multiple
                assert q_value % 3 != 0
                assert q_value % 2 == m_value % 2


def check_finite_rows() -> None:
    """Check the four progression witnesses and two unequal-phase rows."""

    rows = (
        (103, 1, 57, 7, 7, 1),
        (10807, 1, 5901, 49, 7, 2),
        (314, 1, 172, 17, 17, 1),
        (5354, 1, 2924, 289, 17, 2),
        (1118, -1, 610, 49, 7, 2),
        (1345, -1, 734, 49, 7, 2),
    )
    for q_value, sigma, m_value, expected, prime, exponent in rows:
        exact_window(q_value, sigma, m_value)
        direct, diagonal = exact_cross_gcd(q_value, sigma, m_value)
        assert direct == diagonal == expected
        assert direct % (prime**exponent) == 0
        assert normalized_gap_norm(q_value, sigma, m_value) < -1

    for q_value, m_value, prime, exponent in (
        (103, 57, 7, 1),
        (10807, 5901, 7, 2),
        (314, 172, 17, 1),
        (5354, 2924, 17, 2),
    ):
        fib_multiple = (8 if prime == 7 else 9) * prime ** (exponent - 1)
        comp_multiple = (3 if prime == 7 else 4) * prime ** (exponent - 1)
        q_modulus = (3 if prime == 7 else 2) * fib_multiple
        assert q_value % q_modulus == (fib_multiple - 1) % q_modulus
        assert m_value % (2 * comp_multiple) == comp_multiple


def main() -> None:
    """Run the independent exact audit checks."""

    check_rank_lifts()
    check_progression_algebra()
    check_finite_rows()
    print(
        "PASS: independent C0077 imprimitive recurrence audit; "
        "rank lifts, progression clocks, six exact nearest rows, "
        "diagonal gcds, and trace-conic failures verified."
    )


if __name__ == "__main__":
    main()
