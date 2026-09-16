#!/usr/bin/env python3
"""Exact finite regression for the two candidate-window densities.

Window membership is decided with exact integer arithmetic. Floating-point
numbers are used only after counting, to compare the observed proportions
with the analytically proved limiting interval lengths.
"""

from __future__ import annotations

import argparse
from math import log, sqrt

from check_c0068_nearest_even_unit_gap_curator import nearest_gap, selector_data


def pattern_counts(q_max: int, admissible_only: bool) -> dict[str, int]:
    """Count the four exact window-availability patterns."""

    counts = {"both": 0, "plus": 0, "minus": 0, "neither": 0}
    even_powers = [(1, 0)]
    for q_value in range(2, q_max + 1):
        if admissible_only and q_value % 3 == 0:
            continue
        _, _, _, plus_square, minus_square = selector_data(q_value)
        plus = nearest_gap(plus_square, even_powers) is not None
        minus = nearest_gap(minus_square, even_powers) is not None
        if plus and minus:
            counts["both"] += 1
        elif plus:
            counts["plus"] += 1
        elif minus:
            counts["minus"] += 1
        else:
            counts["neither"] += 1
    return counts


def main(q_max: int = 10000) -> None:
    """Compare exact finite counts with the proved limiting densities."""

    assert q_max >= 5000
    lambda_value = 1.0 + sqrt(2.0)
    phi = (1.0 + sqrt(5.0)) / 2.0
    theta = log(sqrt(2.0)) / log(lambda_value)
    delta = log(
        (sqrt(2.0) * phi + 1.0)
        / (lambda_value * (sqrt(2.0) * phi - 1.0))
    ) / log(lambda_value)
    expected = {
        "both": theta - delta,
        "plus": delta,
        "minus": delta,
        "neither": 1.0 - theta - delta,
    }

    reports: dict[str, tuple[dict[str, int], dict[str, float]]] = {}
    for label, admissible_only in (("all", False), ("3nmidq", True)):
        counts = pattern_counts(q_max, admissible_only)
        total = sum(counts.values())
        proportions = {key: value / total for key, value in counts.items()}
        for key in expected:
            assert abs(proportions[key] - expected[key]) < 0.01
        reports[label] = counts, proportions

    print(
        "PASS: candidate-window density regression; "
        f"q<= {q_max}; theta={theta:.12f}; delta={delta:.12f}; "
        f"expected={expected}; reports={reports}. "
        "Limiting densities use the irrational-rotation proof."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--q-max", type=int, default=10000)
    arguments = parser.parse_args()
    main(arguments.q_max)

