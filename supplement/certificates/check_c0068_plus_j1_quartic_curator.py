#!/usr/bin/env python3
"""Exact local verification of the C0068 plus-branch j=1 quartics.

Completeness of the hard-coded integral-point lists comes from the documented
Magma ``IntegralQuarticPoints`` run, not from the bounded recurrence loop.
"""

from __future__ import annotations

from fractions import Fraction


EVEN_POINTS = ((-2, -7), (0, -7), (2, -41), (-4, 41))
ODD_POINTS = ((-2, -1), (0, -1))


def evaluate(coefficients: tuple[int, ...], value: int) -> int:
    result = 0
    for coefficient in coefficients:
        result = result * value + coefficient
    return result


def trim(polynomial: list[Fraction]) -> list[Fraction]:
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def derivative(polynomial: list[Fraction]) -> list[Fraction]:
    return [index * coefficient for index, coefficient in enumerate(polynomial)][1:]


def remainder(
    dividend: list[Fraction], divisor: list[Fraction]
) -> list[Fraction]:
    value = trim(dividend[:])
    divisor = trim(divisor[:])
    while len(value) >= len(divisor):
        scale = value[-1] / divisor[-1]
        offset = len(value) - len(divisor)
        for index, coefficient in enumerate(divisor):
            value[index + offset] -= scale * coefficient
        trim(value)
    return value


def squarefree(coefficients: tuple[int, ...]) -> bool:
    polynomial = [Fraction(value) for value in reversed(coefficients)]
    divisor = derivative(polynomial)
    while divisor:
        polynomial, divisor = divisor, remainder(polynomial, divisor)
    return len(polynomial) == 1


def fibonacci(index: int) -> int:
    first, second = 0, 1
    for _ in range(index):
        first, second = second, first + second
    return first


def pell(index: int) -> int:
    first, second = 0, 1
    for _ in range(index):
        first, second = second, 2 * second + first
    return first


def main() -> None:
    q_even = (18, 72, 132, 120, 49)
    q_odd = (18, 72, 84, 24, 1)

    # These are exactly Q_s(z) for y=z+1 in
    # C^2=18y^4+24s*y^2+7, with s=+1 and s=-1.
    for z in range(-12, 13):
        y = z + 1
        assert evaluate(q_even, z) == 18 * y**4 + 24 * y**2 + 7
        assert evaluate(q_odd, z) == 18 * y**4 - 24 * y**2 + 7
    assert squarefree(q_even)
    assert squarefree(q_odd)

    for z, companion in EVEN_POINTS:
        assert companion * companion == evaluate(q_even, z)
    for z, companion in ODD_POINTS:
        assert companion * companion == evaluate(q_odd, z)

    even_fibonacci_values = {abs(z + 1) for z, _ in EVEN_POINTS}
    odd_fibonacci_values = {abs(z + 1) for z, _ in ODD_POINTS}
    assert even_fibonacci_values == {1, 3}
    assert odd_fibonacci_values == {1}

    recovered: list[tuple[int, int]] = []
    for q in range(0, 200):
        target = 3 * fibonacci(q) ** 2 + 2 * (-1) ** q
        for m in range(1, 400, 2):
            pell_value = pell(m)
            if pell_value == target:
                recovered.append((q, m))
            if pell_value > target and target >= 0:
                break
    assert recovered == [(1, 1), (2, 3), (4, 5)]

    print(
        "PASS: C0068 plus j=1 quartic; exact Magma point lists verified; "
        "complete recurrence pairs=(1,1),(2,3),(4,5); "
        "bounded q<200 loop is regression only"
    )


if __name__ == "__main__":
    main()
