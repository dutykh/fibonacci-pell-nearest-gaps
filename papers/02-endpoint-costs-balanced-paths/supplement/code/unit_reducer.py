# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Independent unit-by-unit reducer for the fixed norm-89 check.

This implementation performs each forced row subtraction separately. It has
no imports from the batched check engine. The finite input sign tree
records all direct terminals before cone entry and every stable-sign seed.
"""
from __future__ import annotations
from collections import deque
from dataclasses import dataclass
from math import gcd, isqrt

Matrix = tuple[int, int, int, int]
IDENTITY: Matrix = (1, 0, 0, 1)
P1: Matrix = (1, 1, 1, 0)
P2: Matrix = (2, 1, 1, 0)
BLOCKS = {1: (2, 1, 1, 1), 2: (5, 2, 2, 1)}
PREFIX_CAP = 64

@dataclass(frozen=True)
class Seed:
    reflection: int
    optional: int
    prefix: tuple[int, ...]
    state: Matrix
    weight: int


@dataclass(frozen=True)
class Direct:
    reflection: int
    optional: int
    prefix: tuple[int, ...]
    source: tuple[int, int]
    target_raw: tuple[int, int]
    change: int


def multiply(left: Matrix, right: Matrix) -> Matrix:
    """Multiply two row-major 2 by 2 matrices."""

    a, b, c, d = left
    e, f, g, h = right
    return a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h


def determinant(matrix: Matrix) -> int:
    """Return the determinant."""

    a, b, c, d = matrix
    return a * d - b * c


def vector_cost(x: int, y: int) -> int:
    """Count unit row subtractions until two positive coordinates agree."""

    x, y = abs(x), abs(y)
    if not x or not y:
        raise ValueError("vector cost needs two nonzero coordinates")
    cost = 0
    while x != y:
        if x > y:
            quotient = (x - 1) // y
            x -= quotient * y
        else:
            quotient = (y - 1) // x
            y -= quotient * x
        cost += quotient
    return cost


def reduce_rows(matrix: Matrix) -> tuple[Matrix, int]:
    """Emit every subtraction forced simultaneously by both columns."""

    if min(matrix) < 0:
        raise ValueError(f"row reduction needs nonnegative entries: {matrix}")
    a, b, c, d = matrix
    emitted = 0
    while True:
        if a >= c and b >= d and (a > c or b > d):
            a, b = a - c, b - d
            emitted += 1
        elif c >= a and d >= b and (c > a or d > b):
            c, d = c - a, d - b
            emitted += 1
        else:
            break
    return (a, b, c, d), emitted


def stable_row_sign(matrix: Matrix) -> int | None:
    """Return the sign making row two nonnegative, or None if it straddles 0."""

    c, d = matrix[2], matrix[3]
    if c >= 0 and d >= 0:
        return 1
    if c <= 0 and d <= 0:
        return -1
    return None


def orient_second_row(matrix: Matrix, sign: int) -> Matrix:
    """Apply the fixed sign to row two."""

    a, b, c, d = matrix
    return a, b, sign * c, sign * d


def append_terminal(matrix: Matrix) -> tuple[int, int]:
    """Return the first column after appending the terminal digit 2."""

    product = multiply(matrix, P2)
    return product[0], product[2]


def source_matrix(optional: int, prefix: tuple[int, ...]) -> Matrix:
    """Build P(optional) times the doubled input prefix."""

    product = IDENTITY if optional == 0 else (P1 if optional == 1 else P2)
    for digit in prefix:
        product = multiply(product, BLOCKS[digit])
    return product


def split_prefixes(reflection: Matrix, reflection_index: int) -> tuple[list[Seed], list[Direct], int, int]:
    """Partition all doubled words at the first stable second-row sign."""

    seeds: list[Seed] = []
    direct: list[Direct] = []
    visited = 0
    maximum_depth = 0
    for optional in (0, 1, 2):
        initial = IDENTITY if optional == 0 else (P1 if optional == 1 else P2)
        queue: deque[tuple[tuple[int, ...], Matrix]] = deque((((), multiply(reflection, initial)),))
        while queue:
            prefix, transformed = queue.popleft()
            visited += 1
            maximum_depth = max(maximum_depth, len(prefix))
            if visited > PREFIX_CAP:
                raise RuntimeError(f"prefix cap {PREFIX_CAP} reached")
            sign = stable_row_sign(transformed)
            if sign is not None:
                oriented = orient_second_row(transformed, sign)
                reduced, emitted = reduce_rows(oriented)
                input_cost = optional + 2 * sum(prefix)
                seeds.append(
                    Seed(reflection_index, optional, prefix, reduced, emitted - input_cost)
                )
                continue

            # If the doubled word stops here, its terminal 2 is a finite case.
            full_source = multiply(source_matrix(optional, prefix), P2)
            raw_target = multiply(reflection, full_source)
            tx, ty = raw_target[0], raw_target[2]
            block_norm = isqrt(abs(determinant(reflection)))
            if tx and ty and gcd(abs(tx), abs(ty)) == block_norm:
                sx, sy = full_source[0], full_source[2]
                source_cost_value = optional + 2 * sum(prefix) + 1
                direct.append(
                    Direct(
                        reflection_index,
                        optional,
                        prefix,
                        (sx, sy),
                        (tx, ty),
                        vector_cost(tx, ty) - source_cost_value,
                    )
                )
            for digit in (1, 2):
                queue.append((prefix + (digit,), multiply(transformed, BLOCKS[digit])))
    return seeds, direct, visited, maximum_depth


def consume(matrix: Matrix, digit: int) -> tuple[Matrix, int]:
    """Consume one doubled digit and return the normalized state and weight."""

    product = multiply(matrix, BLOCKS[digit])
    reduced, emitted = reduce_rows(product)
    return reduced, emitted - 2 * digit


def finish(matrix: Matrix) -> int:
    """Return target terminal cost minus the source terminal cost 1."""

    x, y = append_terminal(matrix)
    return vector_cost(x, y) - 1
