#!/usr/bin/env python3
"""Verify the fixed arithmetic in the independent j=1 quartic audit."""

from __future__ import annotations

from math import isqrt


PLUS_SHIFTED = [(-4, 41), (-2, 7), (0, 7), (2, 41)]
MINUS_SHIFTED = [(-2, 1), (0, 1)]


def shifted_quartic(x: int, sign: int) -> int:
    y = x + 1
    return 18 * y**4 + sign * 24 * y**2 + 7


def fibonacci(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def pell_coordinates(n: int) -> tuple[int, int]:
    c, p = 1, 0
    for _ in range(n):
        c, p = c + 2 * p, c + p
    return c, p


def main() -> None:
    for x, ordinate in PLUS_SHIFTED:
        assert ordinate * ordinate == shifted_quartic(x, 1)
    for x, ordinate in MINUS_SHIFTED:
        assert ordinate * ordinate == shifted_quartic(x, -1)

    plus_y = sorted(x + 1 for x, _ in PLUS_SHIFTED)
    minus_y = sorted(x + 1 for x, _ in MINUS_SHIFTED)
    assert plus_y == [-3, -1, 1, 3]
    assert minus_y == [-1, 1]

    # Recover the original Pell-coordinate equation exactly.
    for sign, values in ((1, plus_y), (-1, minus_y)):
        for y in values:
            pell_p = 3 * y * y + 2 * sign
            pell_c_squared = 2 * pell_p * pell_p - 1
            pell_c = isqrt(pell_c_squared)
            assert pell_c * pell_c == pell_c_squared

    # Filter the complete integer-y lists by y=F_q and sign=(-1)^q.
    fibonacci_hits: list[tuple[int, int, int, int]] = []
    for q in range(0, 20):
        y = fibonacci(q)
        sign = -1 if q % 2 else 1
        allowed = plus_y if sign == 1 else minus_y
        if y not in allowed:
            continue
        target = 3 * y * y + 2 * sign
        for m in range(1, 20, 2):
            _, pell_p = pell_coordinates(m)
            if pell_p == target:
                fibonacci_hits.append((q, m, y, sign))

    assert fibonacci_hits == [(1, 1, 1, -1), (2, 3, 1, 1), (4, 5, 3, 1)]
    assert fibonacci(7) == 13

    # A bounded independent falsification check beyond the complete lists.
    for sign, allowed in ((1, set(plus_y)), (-1, set(minus_y))):
        for y in range(-10_000, 10_001):
            value = 18 * y**4 + sign * 24 * y**2 + 7
            root = isqrt(value)
            assert (root * root == value) == (y in allowed)

    print(
        "PASS: independent j=1 shifted-quartic arithmetic; "
        f"plus_y={plus_y}; minus_y={minus_y}; "
        f"fibonacci_hits={fibonacci_hits}; bounded_abs_y=10000."
    )


if __name__ == "__main__":
    main()
