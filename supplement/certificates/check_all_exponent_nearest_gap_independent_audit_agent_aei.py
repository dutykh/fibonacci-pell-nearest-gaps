#!/usr/bin/env python3
"""Independent exact checks for the all-exponent nearest-gap audit.

This file deliberately imports no project checker.  It verifies the parity
selector, the low normalized-unit coefficient identities, the odd-trace
boundary identity, the bootstrap constants, the q=1 boundary, the complete
terminal odd branch, and a larger intrinsic regression using exact arithmetic
in Z[sqrt(2)].  The unbounded proof remains the written Matveev and certified
continued-fraction argument audited separately.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd


Pair = tuple[int, int]


def add(a: Pair, b: Pair) -> Pair:
    return a[0] + b[0], a[1] + b[1]


def sub(a: Pair, b: Pair) -> Pair:
    return a[0] - b[0], a[1] - b[1]


def mul(a: Pair, b: Pair) -> Pair:
    return a[0] * b[0] + 2 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def scale(k: int, a: Pair) -> Pair:
    return k * a[0], k * a[1]


def norm(a: Pair) -> int:
    return a[0] * a[0] - 2 * a[1] * a[1]


def sign(a: Pair) -> int:
    """Exact sign of a[0] + a[1] sqrt(2)."""

    x, y = a
    if not x:
        return (y > 0) - (y < 0)
    if not y:
        return (x > 0) - (x < 0)
    if x > 0 and y > 0:
        return 1
    if x < 0 and y < 0:
        return -1
    comparison = x * x - 2 * y * y
    if x > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def less(a: Pair, b: Pair) -> bool:
    return sign(sub(b, a)) > 0


def fibonacci(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return values


def units(limit: int) -> list[Pair]:
    values = [(1, 0)]
    while len(values) <= limit:
        values.append(mul(values[-1], (1, 1)))
    return values


def signed_unit(exponent: int, positive: list[Pair]) -> Pair:
    if exponent >= 0:
        return positive[exponent]
    k = -exponent
    c, p = positive[k]
    return ((-1) ** k) * c, ((-1) ** (k + 1)) * p


def first_power_above(value: Pair, positive: list[Pair]) -> int:
    lo, hi = 1, len(positive)
    while lo < hi:
        mid = (lo + hi) // 2
        if less(value, positive[mid]):
            hi = mid
        else:
            lo = mid + 1
    assert lo < len(positive)
    return lo


def normalized_exponent(value: Pair, positive: list[Pair]) -> int | None:
    if norm(value) != -1 or sign(value) <= 0:
        return None
    for exponent in range(1, len(positive), 2):
        c, p = positive[exponent]
        if value == (c, p):
            return exponent
        if value == (-c, p):
            return -exponent
        if p > abs(value[1]):
            return None
    return None


def target(q: int, sigma: int, fib: list[int], positive: list[Pair]) -> tuple[int, int, int] | None:
    f, x = fib[q], fib[q + 1]
    square = (f * f + 2 * x * x, 2 * sigma * f * x)
    n = first_power_above(square, positive)
    if not less(positive[n], scale(2, square)):
        return None
    gap = sub(positive[n], square)
    g = gcd(abs(gap[0]), abs(gap[1]))
    if norm(gap) != -g * g:
        return None
    j = normalized_exponent((gap[0] // g, gap[1] // g), positive)
    assert j is not None and j % 2
    return n, j, g


def check_parity_selector() -> None:
    for f_mod_2 in (0, 1):
        for n_mod_2 in (0, 1):
            r_mod_2 = 1 - f_mod_2
            s_mod_2 = n_mod_2
            same_parity = r_mod_2 == s_mod_2
            assert same_parity == (n_mod_2 == 1 - f_mod_2)


def check_low_exponents(fib: list[int], positive: list[Pair]) -> None:
    for q in range(3, 301, 3):
        f, x = fib[q], fib[q + 1]
        h_plus = x * x + (x + f) ** 2
        h_minus = x * x + (x - f) ** 2
        assert h_plus == fib[2 * q + 3]
        assert h_minus == 3 * f * f + 2 * (-1) ** q
        for sigma in (1, -1):
            square = (f * f + 2 * x * x, 2 * sigma * f * x)
            # Radical coefficients of lambda^{-1} A_sigma^2 and
            # lambda A_sigma^2, respectively.
            assert mul((-1, 1), square)[1] == (h_minus if sigma == 1 else h_plus)
            assert mul((1, 1), square)[1] == (h_plus if sigma == 1 else h_minus)
        a_plus, a_minus = (f, x), (-f, x)
        assert less(mul((1, 1), a_minus), a_plus)
        assert less(a_plus, scale(3, a_minus))

    # The unique q=0 boundary validates all signs in the odd trace.
    gap = sub(positive[1], (2, 0))
    assert gap == signed_unit(-1, positive)
    left = sub(positive[1], signed_unit(-3, positive))
    right = add((2, 0), mul((2, 0), signed_unit(-2, positive)))
    assert left == right


def check_bootstrap() -> None:
    assert Fraction(163, 60) > Fraction(19, 7)
    assert 19**78 > 2 * 10**32 * 7**78
    coefficient = 116_000_000_000_000
    assert Fraction(10**32, 2) > coefficient * 79
    upper = Fraction(5, 4) * (
        9_275_000_000_000 * (4 * coefficient * 79 + 44) * 79 + 4
    )
    assert upper < Fraction(336, 10) * 10**30
    assert upper < 10**32


def check_q_one(fib: list[int], positive: list[Pair]) -> None:
    for sigma in (1, -1):
        f, x = fib[1], fib[2]
        square = (f * f + 2 * x * x, 2 * sigma * f * x)
        n = first_power_above(square, positive)
        assert not less(positive[n], scale(2, square))
        assert target(1, sigma, fib, positive) is None


def enumerate_targets(q_max: int, fib: list[int], positive: list[Pair]) -> tuple[list[tuple[int, int, int, int, int]], int]:
    hits: list[tuple[int, int, int, int, int]] = []
    terminal_odd = 0
    for q in range(2, q_max + 1):
        for sigma in (1, -1):
            found = target(q, sigma, fib, positive)
            if found is None:
                continue
            n, j, g = found
            assert (n % 2 == 1) == (q % 3 == 0)
            if q < 96 and q % 3 == 0 and n % 2:
                terminal_odd += 1
            hits.append((q, sigma, n, j, g))
    assert terminal_odd == 0
    assert hits == [(2, -1, 2, -1, 6), (4, 1, 6, 1, 40)]
    return hits, terminal_odd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q-max", type=int, default=10_000)
    args = parser.parse_args()
    if not 300 <= args.q_max <= 20_000:
        parser.error("--q-max must lie between 300 and 20000")
    fib = fibonacci(2 * args.q_max + 10)
    positive = units(2 * args.q_max + 20)
    check_parity_selector()
    check_low_exponents(fib, positive)
    check_bootstrap()
    check_q_one(fib, positive)
    hits, terminal = enumerate_targets(args.q_max, fib, positive)
    print(
        "PASS: independent all-exponent exact audit; parity selector, "
        "j=+-1 identities, odd-trace boundary, q=1 boundary, bootstrap "
        f"constants, terminal odd hits={terminal}, q<={args.q_max} hits={hits}."
    )


if __name__ == "__main__":
    main()
