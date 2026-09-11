#!/usr/bin/env python3
"""Checks for the imprimitive two-clock recurrence obstruction.

The accompanying proof uses equidistribution and is not computational. This
script verifies the rank and lifting data, finds exact nearest-window examples
in the prescribed progressions, and checks the two displayed 7-squared rows.
"""

from __future__ import annotations

from math import ceil, gcd, log


Pair = tuple[int, int]
LAMBDA = 1.0 + 2.0**0.5
PHI = (1.0 + 5.0**0.5) / 2.0
ALPHA = log(PHI) / log(LAMBDA)
BETA = log((1.0 + 2.0**0.5 * PHI) / 5.0**0.5) / log(LAMBDA)
THETA = log(2.0**0.5) / log(LAMBDA)


def multiply(left: Pair, right: Pair) -> Pair:
    """Multiply two elements of Z[sqrt(2)]."""

    a_value, b_value = left
    c_value, d_value = right
    return (
        a_value * c_value + 2 * b_value * d_value,
        a_value * d_value + b_value * c_value,
    )


def power(base: Pair, exponent: int) -> Pair:
    """Raise one quadratic integer to a nonnegative power."""

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

    def doubled(n_value: int) -> Pair:
        if n_value == 0:
            return 0, 1
        first, second = doubled(n_value // 2)
        even = first * (2 * second - first)
        odd = first * first + second * second
        if n_value % 2:
            return odd, even + odd
        return even, odd

    return doubled(index)


def sign_surd(value: Pair) -> int:
    """Return the exact sign of a+b*sqrt(2)."""

    rational, radical = value
    if rational >= 0 and radical >= 0:
        return 1
    if rational <= 0 and radical <= 0:
        return -1
    comparison = rational * rational - 2 * radical * radical
    if rational > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def subtract(left: Pair, right: Pair) -> Pair:
    """Subtract coordinate pairs."""

    return left[0] - right[0], left[1] - right[1]


def exact_nearest_window(
    q_value: int,
    m_value: int,
    sigma: int = 1,
) -> bool:
    """Check D_sigma(q)<lambda^m<sqrt(2)D_sigma(q) exactly."""

    fib, fib_next = fibonacci_pair(q_value)
    root = sigma * fib, fib_next
    unit = pell_pair(m_value)
    upper = 2 * fib_next, sigma * fib
    return (
        sign_surd(subtract(unit, root)) > 0
        and sign_surd(subtract(upper, unit)) > 0
    )


def modular_fibonacci(index: int, modulus: int) -> int:
    """Return F_index modulo modulus by fast doubling."""

    def doubled(n_value: int) -> Pair:
        if n_value == 0:
            return 0, 1
        first, second = doubled(n_value // 2)
        even = first * (2 * second - first) % modulus
        odd = (first * first + second * second) % modulus
        if n_value % 2:
            return odd, (even + odd) % modulus
        return even, odd

    return doubled(index)[0]


def modular_companion(index: int, modulus: int) -> int:
    """Return C_index modulo modulus."""

    result = (1, 0)
    base = (1, 1)
    exponent = index
    while exponent:
        if exponent & 1:
            result = (
                (result[0] * base[0] + 2 * result[1] * base[1])
                % modulus,
                (result[0] * base[1] + result[1] * base[0])
                % modulus,
            )
        base = (
            (base[0] * base[0] + 2 * base[1] * base[1]) % modulus,
            2 * base[0] * base[1] % modulus,
        )
        exponent >>= 1
    return result[0]


def first_progression_witness(
    prime: int,
    exponent: int,
    search_steps: int,
) -> tuple[int, int, int]:
    """Find one exact nearest row in the proof's prescribed progression."""

    if prime == 7:
        fib_multiple = 8 * prime ** (exponent - 1)
        q_start = fib_multiple - 1
        q_step = 3 * fib_multiple
        comp_multiple = 3 * prime ** (exponent - 1)
    else:
        assert prime == 17
        fib_multiple = 9 * prime ** (exponent - 1)
        q_start = fib_multiple - 1
        q_step = 2 * fib_multiple
        comp_multiple = 4 * prime ** (exponent - 1)

    m_modulus = 2 * comp_multiple
    modulus = prime**exponent
    for step in range(search_steps):
        q_value = q_start + q_step * step
        approximation = q_value * ALPHA + BETA
        m_value = ceil(approximation)
        if m_value % m_modulus != comp_multiple:
            continue
        if not 0.0 < m_value - approximation < THETA:
            continue
        if not exact_nearest_window(q_value, m_value):
            continue
        assert q_value % 3 != 0
        assert q_value % 2 == m_value % 2
        assert modular_fibonacci(q_value + 1, modulus) == 0
        assert modular_companion(m_value, modulus) == 0
        return q_value, m_value, step
    raise AssertionError("search bound did not contain a progression witness")


def cross_gcd(q_value: int, m_value: int, sigma: int) -> int:
    """Return gcd(R_q-s, 2F_(q+1)P_m-sigma F_q C_m)."""

    fib, fib_next = fibonacci_pair(q_value)
    companion, pell = pell_pair(m_value)
    norm = 2 * fib_next * fib_next - fib * fib
    sign = -1 if m_value % 2 else 1
    cross = 2 * fib_next * pell - sigma * fib * companion
    return gcd(norm - sign, cross)


def normalized_gap_norm(q_value: int, m_value: int, sigma: int) -> int:
    """Return the norm of the coefficient-primitive full nearest gap."""

    fib, fib_next = fibonacci_pair(q_value)
    root = sigma * fib, fib_next
    difference = subtract(multiply(pell_pair(m_value), pell_pair(m_value)), multiply(root, root))
    gap_content = gcd(abs(difference[0]), abs(difference[1]))
    gap_norm = difference[0] * difference[0] - 2 * difference[1] * difference[1]
    assert gap_norm % (gap_content * gap_content) == 0
    return gap_norm // (gap_content * gap_content)


def check_rank_data() -> None:
    """Verify the two base ranks and their first valuation lifts."""

    assert fibonacci_pair(8)[0] == 21
    assert pell_pair(3)[0] == 7
    assert fibonacci_pair(9)[0] == 34
    assert pell_pair(4)[0] == 17
    for prime, fib_rank, comp_rank in ((7, 8, 3), (17, 9, 4)):
        assert modular_fibonacci(fib_rank, prime) == 0
        assert modular_companion(comp_rank, prime) == 0
        assert modular_fibonacci(fib_rank, prime * prime) != 0
        assert modular_companion(comp_rank, prime * prime) != 0
        assert modular_fibonacci(fib_rank * prime, prime * prime) == 0
        assert modular_companion(comp_rank * prime, prime * prime) == 0


def main() -> None:
    """Run the exact finite regressions supporting the symbolic proof."""

    check_rank_data()
    witnesses = (
        (7, 1, 7, first_progression_witness(7, 1, 200)),
        (7, 2, 49, first_progression_witness(7, 2, 1000)),
        (17, 1, 17, first_progression_witness(17, 1, 300)),
        (17, 2, 289, first_progression_witness(17, 2, 3000)),
    )
    for prime, exponent, expected, (q_value, m_value, _) in witnesses:
        assert cross_gcd(q_value, m_value, 1) == expected
        assert expected % (prime**exponent) == 0
        assert normalized_gap_norm(q_value, m_value, 1) < -1

    for q_value, m_value in ((1118, 610), (1345, 734)):
        assert exact_nearest_window(q_value, m_value, -1)
        assert cross_gcd(q_value, m_value, -1) == 49
        assert normalized_gap_norm(q_value, m_value, -1) < -1

    print(
        "PASS: imprimitive two-clock recurrence; "
        f"progression_witnesses={witnesses}; "
        "exact exceptional rows=(1118,-1,610,49),(1345,-1,734,49)."
    )


if __name__ == "__main__":
    main()
