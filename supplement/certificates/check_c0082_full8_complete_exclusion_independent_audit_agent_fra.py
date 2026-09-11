#!/usr/bin/env python3
"""Independent rational-interval audit of the C0082 terminal reduction.

The checker has three deliberately separate parts.

1. Fixed-point rational intervals certify one common convergent of
   log(phi)/log(lambda), both fixed-coefficient reduction epsilons, and all
   190 moving-coefficient epsilons for odd 3 <= h < 192.
2. Integer arithmetic checks the explicit Matveev cutoff q < 10^32 and the
   elementary constants used before the reduction.
3. A direct exact search in Q(sqrt(2)) reconstructs every q-only nearest gap
   for 2 <= q < 90.  It does not import any earlier nearest-gap checker.

All transcendental enclosures come from the rational series

    log(x) = 2 * sum_{j >= 0} y^(2j+1)/(2j+1),
    y = (x-1)/(x+1),

with an explicit geometric tail.  No floating-point or third-party package
is used.

Run from math-sandbox/:
    python3 -B scripts/check_c0082_full8_complete_exclusion_independent_audit_agent_fra.py
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt


DIGITS = 130
SCALE = 10**DIGITS
LOG_TERMS = 500

CONVERGENT_P = 2585762774943540084566744409566833
CONVERGENT_Q = 4736007914708481810221186278309584
EXPONENT_BOUND = 2 * 10**32


def ceil_div(numerator: int, denominator: int) -> int:
    """Return ceil(numerator / denominator) for a positive denominator."""
    assert denominator > 0
    return -((-numerator) // denominator)


@dataclass(frozen=True)
class Interval:
    """Closed rational interval [lower/SCALE, upper/SCALE]."""

    lower: int
    upper: int

    def __post_init__(self) -> None:
        assert self.lower <= self.upper

    @staticmethod
    def rational(numerator: int, denominator: int = 1) -> "Interval":
        assert denominator > 0
        return Interval(
            (numerator * SCALE) // denominator,
            ceil_div(numerator * SCALE, denominator),
        )

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lower + other.lower, self.upper + other.upper)

    def __neg__(self) -> "Interval":
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return Interval(
            min(products) // SCALE,
            ceil_div(max(products), SCALE),
        )

    def multiply_integer(self, value: int) -> "Interval":
        if value >= 0:
            return Interval(self.lower * value, self.upper * value)
        return Interval(self.upper * value, self.lower * value)

    def divide_integer(self, value: int) -> "Interval":
        assert value > 0
        return Interval(self.lower // value, ceil_div(self.upper, value))

    def reciprocal(self) -> "Interval":
        assert self.lower > 0
        return Interval(
            SCALE * SCALE // self.upper,
            ceil_div(SCALE * SCALE, self.lower),
        )

    def divide(self, other: "Interval") -> "Interval":
        return self * other.reciprocal()

    def power(self, exponent: int) -> "Interval":
        assert exponent >= 0
        answer = Interval.rational(1)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                answer = answer * base
            base = base * base
            remaining //= 2
        return answer


ZERO = Interval.rational(0)
ONE = Interval.rational(1)


def square_root_interval(radicand: int) -> Interval:
    lower = isqrt(radicand * SCALE * SCALE)
    return Interval(lower, lower + 1)


def positive_log_point(value: int) -> Interval:
    """Enclose log(value/SCALE) by a rational atanh series."""
    assert value > 0
    if value == SCALE:
        return ZERO

    if value > SCALE:
        numerator = value - SCALE
        denominator = value + SCALE
        sign = 1
    else:
        numerator = SCALE - value
        denominator = SCALE + value
        sign = -1

    y_value = Interval.rational(numerator, denominator)
    y_squared = y_value * y_value
    term = y_value
    partial = ZERO
    for index in range(LOG_TERMS):
        partial = partial + term.divide_integer(2 * index + 1)
        term = term * y_squared

    # The first omitted power is term.  Dropping all odd denominators gives
    # the rigorous geometric upper bound term/(1-y^2).
    tail = term.divide(ONE - y_squared)
    magnitude = Interval(
        2 * partial.lower,
        2 * partial.upper + 2 * tail.upper,
    )
    return magnitude if sign > 0 else -magnitude


def logarithm(value: Interval) -> Interval:
    """Use monotonicity to enclose log on a positive interval."""
    assert value.lower > 0
    lower_log = positive_log_point(value.lower)
    upper_log = positive_log_point(value.upper)
    return Interval(lower_log.lower, upper_log.upper)


def absolute_interval(value: Interval) -> Interval:
    if value.lower >= 0:
        return value
    if value.upper <= 0:
        return -value
    return Interval(0, max(-value.lower, value.upper))


def nearest_integer_distance(value: Interval, multiplier: int) -> tuple[Interval, int]:
    """Enclose ||multiplier*value|| and certify its nearest integer."""
    lower = value.lower * multiplier
    upper = value.upper * multiplier
    nearest = (lower + upper + SCALE) // (2 * SCALE)
    half = SCALE // 2
    assert lower > nearest * SCALE - half
    assert upper < nearest * SCALE + half
    displacement = Interval(
        lower - nearest * SCALE,
        upper - nearest * SCALE,
    )
    return absolute_interval(displacement), nearest


def epsilon_interval(mu_value: Interval, gamma_error: Interval) -> Interval:
    mu_distance, _ = nearest_integer_distance(mu_value, CONVERGENT_Q)
    penalty = gamma_error.multiply_integer(EXPONENT_BOUND)
    return mu_distance - penalty


def certify_threshold(
    epsilon: Interval,
    base: Interval,
    threshold: int,
    coefficient: int,
) -> None:
    """Certify log(coefficient*Q/epsilon)/log(base) < threshold."""
    assert epsilon.lower > 0
    product = epsilon * base.power(threshold)
    assert product.lower > coefficient * CONVERGENT_Q * SCALE


def analytic_reduction_audit() -> None:
    root_two = square_root_interval(2)
    root_three = square_root_interval(3)
    root_five = square_root_interval(5)
    pell_unit = ONE + root_two
    golden = (ONE + root_five).divide_integer(2)

    log_lambda = logarithm(pell_unit)
    log_phi = logarithm(golden)
    gamma = log_phi.divide(log_lambda)

    assert 5 * log_lambda.lower > 4 * SCALE
    assert 5 * log_phi.lower > 2 * SCALE
    assert 2 * log_phi.lower > log_lambda.upper  # gamma > 1/2

    gamma_distance, nearest = nearest_integer_distance(gamma, CONVERGENT_Q)
    assert nearest == CONVERGENT_P
    assert gcd(CONVERGENT_P, CONVERGENT_Q) == 1
    # Legendre's criterion: |gamma-p/Q| < 1/(2Q^2).
    assert 2 * CONVERGENT_Q * gamma_distance.upper < SCALE
    assert CONVERGENT_Q > 6 * EXPONENT_BOUND
    # Displayed sharper enclosure |Q*gamma - p| < 26/10^36, consumed by the margins.
    assert gamma_distance.upper * 10**36 < 26 * SCALE

    # Recompute the numerical factor inherited from the stated Matveev
    # specialization.  Here 3^(9/2)=81*sqrt(3).
    log_four = logarithm(Interval.rational(4))
    prefactor = Interval.rational(14, 10)
    prefactor = prefactor.multiply_integer(30**6)
    prefactor = prefactor * root_three.multiply_integer(81)
    prefactor = prefactor.multiply_integer(4**2)
    prefactor = prefactor * (ONE + log_four)
    prefactor = prefactor * log_lambda.multiply_integer(2)
    prefactor = prefactor * log_phi.multiply_integer(2)
    assert prefactor.upper < 9_275_000_000_000 * SCALE

    # The stronger adjusted-form cutoff at q_0=10^32.
    log_two = logarithm(Interval.rational(2))
    log_ten = logarithm(Interval.rational(10))
    h_at_cutoff = ONE + log_two + log_ten.multiply_integer(32)
    assert h_at_cutoff.upper < 79 * SCALE
    matveev_bracket = (
        9_275_000_000_000
        * (960_000_000_000_000 * 79 + 44)
        * 79
        + 4
    )
    assert 5 * matveev_bracket < 4 * 10**32
    assert 5 * matveev_bracket < 4 * 7 * 10**31
    assert 240_000_000_000_000 * 79 == 18_960_000_000_000_000

    # Elementary constants in the two normalized remainder estimates.
    assert 40 * 39 < 1597
    assert 39 * 20 < 38 * 21
    assert 201 * 2 + 500 < 1000
    assert 9 * 101**2 < 10 * 100**2
    assert 100 * (5**6 * 3**14 + 2**14 * 12**6) < 12**6 * 3**14

    golden_inverse = golden.reciprocal()
    for sigma in (1, -1):
        u_value = Interval.rational(sigma) + root_two * golden
        v_value = Interval.rational(sigma) - root_two * golden
        u_prime = Interval.rational(-sigma) + root_two * golden_inverse
        v_prime = Interval.rational(-sigma) - root_two * golden_inverse
        assert u_value.lower > SCALE
        assert absolute_interval(u_prime.divide(u_value)).upper < 2 * SCALE
        assert (
            absolute_interval(v_prime).divide(absolute_interval(v_value)).upper
            < 2 * SCALE
        )
        assert absolute_interval(v_value.divide(u_value)).upper < 3 * SCALE
    assert 200 * 2**14 < 3**14

    fixed_base = pell_unit * pell_unit
    fixed_epsilons: dict[int, Interval] = {}
    for sigma in (1, -1):
        u_value = Interval.rational(sigma) + root_two * golden
        u_square_over_five = (u_value * u_value).divide_integer(5)
        mu_value = logarithm(u_square_over_five).divide(log_lambda)
        epsilon = epsilon_interval(mu_value, gamma_distance)
        assert epsilon.lower > 0
        if sigma == 1:
            assert 100 * epsilon.lower > 7 * SCALE
        else:
            assert 5 * epsilon.lower > 2 * SCALE
        certify_threshold(epsilon, fixed_base, 48, 30)
        fixed_epsilons[sigma] = epsilon

    # Once min(h,q*gamma)<48, gamma>1/2, h<n, and n<2q give h<192.
    moving_base = golden * golden
    pell_inverse = pell_unit.reciprocal()
    moving_rows: list[tuple[int, int, Interval]] = []
    for sigma in (1, -1):
        u_value = Interval.rational(sigma) + root_two * golden
        v_value = Interval.rational(sigma) - root_two * golden
        u_squared = u_value * u_value
        v_squared = v_value * v_value
        for height in range(3, 192, 2):
            adjusted = (
                u_squared
                + v_squared * pell_inverse.power(2 * height)
            ).divide_integer(5)
            mu_value = logarithm(adjusted).divide(log_lambda)
            epsilon = epsilon_interval(mu_value, gamma_distance)
            assert epsilon.lower > 0
            assert 1000 * epsilon.lower > 9 * SCALE
            certify_threshold(epsilon, moving_base, 90, 27)
            moving_rows.append((sigma, height, epsilon))

    assert len(moving_rows) == 190
    weakest = min(moving_rows, key=lambda row: row[2].lower)
    assert weakest[0:2] == (1, 31)


Pair = tuple[int, int]


def multiply_sqrt_two(left: Pair, right: Pair) -> Pair:
    a_value, b_value = left
    c_value, d_value = right
    return a_value * c_value + 2 * b_value * d_value, a_value * d_value + b_value * c_value


def sign_sqrt_two(value: Pair) -> int:
    """Return the exact sign of a+b*sqrt(2)."""
    rational, radical = value
    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if (rational > 0) == (radical > 0):
        return 1 if rational > 0 else -1
    if rational * rational > 2 * radical * radical:
        return 1 if rational > 0 else -1
    return 1 if radical > 0 else -1


def subtract_sqrt_two(left: Pair, right: Pair) -> Pair:
    return left[0] - right[0], left[1] - right[1]


def fibonacci_pair(index: int) -> tuple[int, int]:
    first, second = 0, 1
    for _ in range(index):
        first, second = second, first + second
    return first, second


def exact_terminal_enumeration() -> None:
    """Reconstruct the q-only nearest gaps without logarithms or estimates."""
    targets: list[tuple[int, int, int, int, Pair]] = []
    pell_even: Pair = (1, 0)
    lambda_squared: Pair = (3, 2)

    # The largest needed even Pell power is generated once and reused.
    even_powers: list[Pair] = [pell_even]
    for _ in range(1, 2 * 90 + 1):
        even_powers.append(multiply_sqrt_two(even_powers[-1], lambda_squared))

    for q_index in range(2, 90):
        if q_index % 3 == 0:
            continue
        fibonacci, successor = fibonacci_pair(q_index)
        level = fibonacci * fibonacci + 2 * successor * successor
        for sigma in (1, -1):
            square = (level, 2 * sigma * fibonacci * successor)

            candidate: tuple[int, Pair] | None = None
            crossed = False
            for half_index, unit in enumerate(even_powers[1:], start=1):
                difference = subtract_sqrt_two(unit, square)
                if sign_sqrt_two(difference) <= 0:
                    continue
                crossed = True
                upper_difference = subtract_sqrt_two(
                    (2 * square[0], 2 * square[1]),
                    unit,
                )
                if sign_sqrt_two(upper_difference) > 0:
                    candidate = (2 * half_index, difference)
                break

            assert crossed

            if candidate is None:
                continue
            even_exponent, gap = candidate
            content = gcd(abs(gap[0]), abs(gap[1]))
            gap_norm = gap[0] * gap[0] - 2 * gap[1] * gap[1]
            if gap_norm != -(content * content):
                continue
            primitive = (gap[0] // content, gap[1] // content)
            targets.append((q_index, sigma, even_exponent, content, primitive))

    assert targets == [
        (2, -1, 2, 6, (-1, 1)),
        (4, 1, 6, 40, (1, 1)),
    ]
    assert not [row for row in targets if 7 <= row[0] < 90]


def main() -> None:
    analytic_reduction_audit()
    exact_terminal_enumeration()
    print(
        "PASS: independent C0082 terminal audit; q<10^32 cutoff checked; "
        "common convergent certified by rational log intervals; fixed form "
        "gives min(h,q*gamma)<48 and h<192; 190 moving forms give q<90; "
        "exact q-only targets below 90 are q=2 and q=4, so no residual "
        "q>=7 negative-gap target remains"
    )


if __name__ == "__main__":
    main()
