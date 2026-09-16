# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact primitive Gaussian reflections and componentwise Euclidean reduction.

All arithmetic is integral. Row subtractions are emitted in maximal batches;
their sum is the number of actual unit subtractions. No external package is
required. This module contains only the operations used by the fixed B=89
check and the printed fixed examples.
"""
from __future__ import annotations
from collections import deque
from dataclasses import dataclass
from math import gcd
from resource import RUSAGE_SELF, getrusage
from time import monotonic

Mat = tuple[int, int, int, int]
Word = tuple[int, ...]
Pair = tuple[int, int]
ONE: Mat = (1, 0, 0, 1)
J = {1: (1, 1, 1, 0), 2: (2, 1, 1, 0)}
DOUBLE = {1: (2, 1, 1, 1), 2: (5, 2, 2, 1)}
PREFIX_CAP = 64
WALL_CAP_SECONDS = 28.0
MEMORY_CAP_KB = 128 * 1024

@dataclass(frozen=True)
class Seed:
    branch: int
    optional: int | None
    prefix: Word
    state: Mat
    weight: int


@dataclass(frozen=True)
class Edge:
    target: Mat
    weight: int
    emitted: int


@dataclass(frozen=True)
class Direct:
    branch: int
    optional: int | None
    prefix: Word
    change: int


def mul(left: Mat, right: Mat) -> Mat:
    a, b, c, d = left
    e, f, g, h = right
    return a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h


def det(matrix: Mat) -> int:
    a, b, c, d = matrix
    return a * d - b * c


def content(values: tuple[int, ...]) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def inverse_unit(matrix: Mat) -> Mat:
    a, b, c, d = matrix
    determinant = det(matrix)
    if abs(determinant) != 1:
        raise AssertionError("requested inverse of a non-unimodular matrix")
    return d // determinant, -b // determinant, -c // determinant, a // determinant


def matrices(c: int, d: int, block_norm: int) -> tuple[Mat, Mat]:
    big_c = c * c - d * d
    big_f = 2 * c * d
    if c * c + d * d != block_norm or gcd(c, d) != 1:
        raise AssertionError("bad primitive Gaussian block")
    if gcd(big_c, big_f) != 1:
        raise AssertionError("reflection numerator has nontrivial content")
    result = (
        (big_c, big_f, big_f, -big_c),
        (big_f, big_c, big_c, -big_f),
    )
    for reflection in result:
        if det(reflection) != -(block_norm * block_norm):
            raise AssertionError("reflection determinant failed")
        if mul(reflection, reflection) != (
            block_norm * block_norm,
            0,
            0,
            block_norm * block_norm,
        ):
            raise AssertionError("reflection square failed")
    return result


def word_matrix(word: Word) -> Mat:
    answer = ONE
    for digit in word:
        answer = mul(answer, J[digit])
    return answer


def make_word(optional: int | None, blocks: Word) -> Word:
    initial = () if optional is None else (optional,)
    doubled = tuple(digit for block in blocks for digit in (block, block))
    return initial + doubled + (2,)


def source_prefix(optional: int | None, blocks: Word) -> Mat:
    answer = ONE if optional is None else J[optional]
    for block in blocks:
        answer = mul(answer, DOUBLE[block])
    return answer


def vector_cost(first: int, second: int) -> int:
    first, second = abs(first), abs(second)
    if not first or not second:
        raise AssertionError("zero coordinate in Euclidean cost")
    answer = 0
    while first != second:
        if first > second:
            step = (first - 1) // second
            first -= step * second
        else:
            step = (second - 1) // first
            second -= step * first
        answer += step
    return answer


def terminal_vector(matrix: Mat) -> Pair:
    product = mul(matrix, J[2])
    return product[0], product[2]


def terminal_weight(matrix: Mat) -> int:
    return vector_cost(*terminal_vector(matrix)) - 1


def row_two_sign(matrix: Mat) -> int | None:
    c, d = matrix[2], matrix[3]
    if c >= 0 and d >= 0:
        return 1
    if c <= 0 and d <= 0:
        return -1
    return None


def orient_row_two(matrix: Mat, sign: int) -> Mat:
    a, b, c, d = matrix
    oriented = a, b, sign * c, sign * d
    if min(oriented) < 0:
        raise AssertionError("prefix did not enter the nonnegative chamber")
    return oriented


def batch_reduce(matrix: Mat) -> tuple[Mat, int]:
    """Reach the same row-reduced state using maximal exact batches."""

    if min(matrix) < 0:
        raise AssertionError("negative entry in row reducer")
    a, b, c, d = matrix
    emitted = 0
    while True:
        if a >= c and b >= d and (a > c or b > d):
            ratios = []
            if c:
                ratios.append(a // c)
            if d:
                ratios.append(b // d)
            step = min(ratios)
            if step < 1:
                raise AssertionError("invalid first-row subtraction batch")
            a, b = a - step * c, b - step * d
            emitted += step
        elif c >= a and d >= b and (c > a or d > b):
            ratios = []
            if a:
                ratios.append(c // a)
            if b:
                ratios.append(d // b)
            step = min(ratios)
            if step < 1:
                raise AssertionError("invalid second-row subtraction batch")
            c, d = c - step * a, d - step * b
            emitted += step
        else:
            break
    state = a, b, c, d
    if content(state) != 1:
        raise AssertionError("state acquired hidden scalar content")
    return state, emitted


def audit_reduced_state(state: Mat, block_norm: int) -> None:
    """Check the determinant decomposition proving fixed-B finiteness."""

    if min(state) < 0 or abs(det(state)) != block_norm * block_norm:
        raise AssertionError("state left its fixed determinant shell")
    a, b, c, d = state
    if a >= c and b >= d and (a > c or b > d):
        raise AssertionError("first row still dominates")
    if c >= a and d >= b and (c > a or d > b):
        raise AssertionError("second row still dominates")
    if a < c:
        a, b, c, d = c, d, a, b
    if not (a > c and b < d):
        raise AssertionError("row-reduced state is not strictly crossed")
    u, v = a - c, d - b
    decomposition = c * v + b * u + u * v
    if decomposition != block_norm * block_norm:
        raise AssertionError("fixed-determinant finiteness identity failed")
    if max(b, c, u, v) > block_norm * block_norm:
        raise AssertionError("fixed-B coordinate bound failed")


def check_clock(started: float, stage: str) -> None:
    if monotonic() - started > WALL_CAP_SECONDS:
        raise RuntimeError(f"28-second cap reached during {stage}")
    if getrusage(RUSAGE_SELF).ru_maxrss > MEMORY_CAP_KB:
        raise RuntimeError(f"memory cap of 128 mebibytes reached during {stage}")


def prefix_partition(
    reflections: tuple[Mat, Mat],
    block_norm: int,
    started: float,
) -> tuple[list[Seed], list[Direct], int, int]:
    seeds: list[Seed] = []
    direct: list[Direct] = []
    visited = 0
    maximum_depth = 0
    for branch, reflection in enumerate(reflections):
        if min(reflection[0:2]) <= 0:
            raise AssertionError("the reflection's first row must be strictly positive")
        for optional in (None, 1, 2):
            transformed = mul(reflection, ONE if optional is None else J[optional])
            queue: deque[tuple[Word, Mat]] = deque((((), transformed),))
            while queue:
                prefix, current = queue.popleft()
                visited += 1
                maximum_depth = max(maximum_depth, len(prefix))
                if visited > PREFIX_CAP:
                    raise RuntimeError("prefix-partition cap reached")
                check_clock(started, "prefix partition")
                sign = row_two_sign(current)
                if sign is not None:
                    state, emitted = batch_reduce(orient_row_two(current, sign))
                    input_weight = (0 if optional is None else optional) + 2 * sum(prefix)
                    seeds.append(
                        Seed(branch, optional, prefix, state, emitted - input_weight)
                    )
                    continue

                # A straddling row has one positive rational zero. At most one
                # doubled child cylinder may continue to straddle that zero.
                unresolved_children = 0
                for digit in (1, 2):
                    child = mul(current, DOUBLE[digit])
                    if row_two_sign(child) is None:
                        unresolved_children += 1
                    queue.append((prefix + (digit,), child))
                if unresolved_children > 1:
                    raise AssertionError("rational zero entered two child cylinders")

                raw = terminal_vector(current)
                raw_content = gcd(abs(raw[0]), abs(raw[1]))
                if raw[0] and raw[1] and raw_content == block_norm:
                    source_cost = (0 if optional is None else optional) + 2 * sum(prefix) + 1
                    direct.append(
                        Direct(
                            branch,
                            optional,
                            prefix,
                            vector_cost(*raw) - source_cost,
                        )
                    )
    return seeds, direct, visited, maximum_depth


def step(state: Mat, digit: int) -> Edge:
    target, emitted = batch_reduce(mul(state, DOUBLE[digit]))
    return Edge(target, emitted - 2 * digit, emitted)
