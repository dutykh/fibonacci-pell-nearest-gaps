#!/usr/bin/env python3
"""Exact bounded regression for the general-seed orbit/gap bridge.

The proof, not this finite test, establishes the unbounded theorem.
"""

from __future__ import annotations

from math import gcd, isqrt


Pair = tuple[int, int]


def mul(x: Pair, y: Pair) -> Pair:
    """Multiply coefficient pairs representing a+b*sqrt(2)."""
    a, b = x
    c, d = y
    return a * c + 2 * b * d, a * d + b * c


def power(x: Pair, exponent: int) -> Pair:
    assert exponent >= 0
    out = (1, 0)
    base = x
    n = exponent
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n >>= 1
    return out


def lambda_power(exponent: int) -> Pair:
    return power((1, 1) if exponent >= 0 else (-1, 1), abs(exponent))


def pell_table(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(2 * values[-1] + values[-2])
    return values


def step(z: tuple[int, int]) -> tuple[int, int]:
    u, v = z
    return v, u + 2 * v


def iterate(z: tuple[int, int], t: int) -> tuple[int, int]:
    for _ in range(t):
        z = step(z)
    return z


def forms(z: tuple[int, int]) -> tuple[int, int, int]:
    u, v = z
    return u * u + v * v, v * v - u * u + 2 * u * v, v * v - 2 * u * v - u * u


def scalar_multiple_positive(x: Pair, unit: Pair) -> int | None:
    """Return positive g with x=g*unit, or None."""
    candidates: list[int] = []
    for value, divisor in zip(x, unit, strict=True):
        if divisor == 0:
            if value != 0:
                return None
        else:
            if value % divisor:
                return None
            candidates.append(value // divisor)
    if not candidates or any(candidate != candidates[0] for candidate in candidates):
        return None
    return candidates[0] if candidates[0] > 0 else None


def sign_surd(x: Pair) -> int:
    """Exact sign of a+b*sqrt(2)."""
    a, b = x
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    norm = a * a - 2 * b * b
    if a > 0:
        return 1 if norm > 0 else -1
    return -1 if norm > 0 else 1


def square(x: Pair) -> Pair:
    return mul(x, x)


def add(x: Pair, y: Pair) -> Pair:
    return x[0] + y[0], x[1] + y[1]


def subtract(x: Pair, y: Pair) -> Pair:
    return x[0] - y[0], x[1] - y[1]


def scale(k: int, x: Pair) -> Pair:
    return k * x[0], k * x[1]


def main() -> None:
    max_coordinate = 32
    max_d = 31
    max_t = 8
    pell = pell_table(max_d + 2 * max_t + 4)
    canonical = {(pell[r], pell[r + 1]): r for r in range(1, 20)}

    comparisons = 0
    hits = 0
    for u in range(1, max_coordinate + 1):
        for v in range(u, max_coordinate + 1):
            z = (u, v)
            h0, b0, qform = forms(z)
            a_seed = (v - u, u)

            assert a_seed[0] * a_seed[0] - 2 * a_seed[1] * a_seed[1] == qform
            assert (step(z)[1] - step(z)[0], step(z)[0]) == mul((1, 1), a_seed)

            if qform * qform == 1:
                assert z in canonical
                r = canonical[z]
                assert a_seed == lambda_power(r)
                for d in range(1, max_d + 1):
                    for t in range(max_t + 1):
                        zt = iterate(z, t)
                        hit = forms(zt)[0] == pell[d + 2 * t]
                        assert hit == (d == 2 * r + 1)
                continue

            assert qform != 0 and qform * qform > 1
            a_square = square(a_seed)
            for d in range(1, max_d + 1):
                for t in range(max_t + 1):
                    comparisons += 1
                    zt = iterate(z, t)
                    hit = forms(zt)[0] == pell[d + 2 * t]
                    gap = subtract(lambda_power(d - 1), a_square)
                    unit = lambda_power(-(2 * t + 1))
                    g = scalar_multiple_positive(gap, unit)
                    equation = g is not None
                    assert hit == equation, (z, d, t, hit, gap, unit, g)
                    if not hit:
                        continue
                    hits += 1
                    assert g is not None
                    d0 = h0 - pell[d]
                    c_d = pell[d] + pell[d - 1]
                    e0 = c_d - b0
                    c_2t = pell[2 * t] + (pell[2 * t - 1] if t else 1)
                    assert d0 == g * pell[2 * t]
                    assert e0 == g * c_2t
                    assert gcd(abs(d0), abs(e0)) == g
                    assert sign_surd(gap) > 0
                    assert sign_surd(subtract(a_square, lambda_power(d - 1))) < 0
                    assert sign_surd(subtract(scale(2, a_square), lambda_power(d - 1))) > 0
                    assert gap == scale(g, unit)
                    assert gcd(abs(gap[0]), abs(gap[1])) == g
                    assert gap[0] * gap[0] - 2 * gap[1] * gap[1] == -g * g

    # The box deliberately contains both direct and positive-time hits.
    assert hits > 0
    print(
        "PASS: general-seed orbit/nearest-gap regression; "
        f"coordinates<={max_coordinate}, all d<={max_d}, t<={max_t}; "
        f"comparisons={comparisons}, noncanonical hits={hits}."
    )


if __name__ == "__main__":
    main()
