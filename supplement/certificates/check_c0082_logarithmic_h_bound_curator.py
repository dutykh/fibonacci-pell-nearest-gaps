#!/usr/bin/env python3
"""Exact constant checks for the C0082 logarithmic h-bound attempt.

Mathematical target: verify the rational inequalities and conservative
Matveev constant used in
attempts/2026-08-24-c0082-logarithmic-h-bound-curator.md.

Dependencies: Python 3 standard library only. Deterministic command:
    python3 -B scripts/check_c0082_logarithmic_h_bound_curator.py

This checks constants and sample algebraic sign bounds. It does not turn the
remaining unbounded q-range into a finite computation.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction


def main() -> None:
    """Verify every numerical enclosure used by the written proof."""
    assert 2 * 2**14 * 100 < 3**14

    delta_coefficient = Fraction(14) * Fraction(5, 12) ** 3
    assert delta_coefficient == Fraction(1750, 1728)

    binet_error = Fraction(201, 100) * delta_coefficient
    gap_error = 3 * Fraction(101, 100) ** 2
    total_error = binet_error + gap_error
    assert binet_error == Fraction(2345, 1152)
    assert gap_error == Fraction(30603, 10000)
    assert total_error == Fraction(3669041, 720000)
    assert total_error < Fraction(128, 25)
    assert Fraction(128, 25) * Fraction(5, 12) ** 3 == Fraction(10, 27)
    assert Fraction(27, 17) * Fraction(128, 25) < 9

    # The elementary n < 2q comparison only needs this fixed inequality.
    assert 7**14 > 32 * 5**14

    getcontext().prec = 80
    one = Decimal(1)
    two = Decimal(2)
    four = Decimal(4)
    sqrt_two = two.sqrt()
    sqrt_five = Decimal(5).sqrt()
    pell_unit = one + sqrt_two
    golden = (one + sqrt_five) / two

    for sigma in (one, -one):
        u_value = sigma + sqrt_two * golden
        v_value = -sigma + sqrt_two / golden
        assert u_value > 0
        assert abs(v_value / u_value) < 2
        alpha = u_value * u_value / Decimal(5)
        assert alpha > 0
        assert abs(alpha.ln()) < 20

    base_constant = (
        Decimal("1.4")
        * Decimal(30) ** 6
        * Decimal(3) ** Decimal("4.5")
        * four**2
        * (one + four.ln())
        * (2 * pell_unit.ln())
        * (2 * golden.ln())
    )
    assert base_constant < Decimal("9.275e12")
    matveev_constant = 20 * base_constant
    assert matveev_constant < Decimal("1.855e14")
    assert pell_unit.ln() > Decimal(4) / Decimal(5)
    assert (Decimal("1.855e14") + Decimal(9).ln()) / pell_unit.ln() < Decimal(
        "2.4e14"
    )

    print(
        "PASS: C0082 logarithmic h-bound constants; "
        "analytic_error<128/25*lambda^-h; |Lambda|<9*lambda^-h; "
        "Matveev_constant<1.855e14; rounded_h_coefficient<2.4e14. "
        "This does not bound q."
    )


if __name__ == "__main__":
    main()
