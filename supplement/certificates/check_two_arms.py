#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#          and Technology, Abu Dhabi, UAE)
#          Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Certificate for the section "Why the Fibonacci and Pell arms".

Every predicate below is an exact integer or Fraction comparison; no
floating-point value enters a decision, and no third-party package is used.

The section separates proved statements from finite verifications, and so does
this checker.  Each block prints either PROVED (the block re-checks, on
examples, a statement that the manuscript proves in general) or VERIFIED (the
block establishes a statement that the manuscript asserts only for the finite
family scanned here).

Blocks:

1.  Cohn normalisation.  For every Christoffel word to depth 11, the Cohn
    matrix satisfies det = 1, trace = 3m with m the (1,2) entry, and
    1 <= D <= m, so the characteristic root is u = A - 2m and u^2 + 1 = 0
    mod m.  This re-checks the background normalisation quoted in the
    manuscript.

2.  Letter counts on the two arms.  The closed forms
    (m_a, m_b) = (F_2k, F_2k-1) on the Fibonacci arm and (C_2k, P_2k) on the
    Pell arm, together with the explicit Cohn matrices used in the proof.

3.  Branch recurrence.  Along a branch fixing m0, the letter-count vector
    obeys v_(k+1) = 3 m0 v_k - v_(k-1), checked on ten branches.

4.  Trace criterion.  The scan for 3*m0 - 2 a perfect square over the 4097
    Markoff numbers of depth 11.  This is the VERIFIED statement of the
    section's computational-scope remark.

5.  One displayed fixed-34 ray.  The explicit rational linear system showing
    that no integral half-step exists on the LLL|R^inf ray used in the
    manuscript.  This does not classify every ray fixing 34.

6.  Growth rates.  The exact trace comparison 322 > 198 > 194 over a common
    length of twelve letters, which refutes the tempting extremality claim
    that the section explicitly does not make.

Run from the manuscript directory:

    python3 -B supplement/certificates/check_two_arms.py
"""

from __future__ import annotations

import math
from fractions import Fraction

Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]

COHN_A: Matrix = ((2, 1), (1, 1))
COHN_B: Matrix = ((5, 2), (2, 1))
FIB_MATRIX: Matrix = ((1, 1), (1, 0))
PELL_MATRIX: Matrix = ((2, 1), (1, 0))

TREE_DEPTH = 11


def matmul(x: Matrix, y: Matrix) -> Matrix:
    return (
        (x[0][0] * y[0][0] + x[0][1] * y[1][0], x[0][0] * y[0][1] + x[0][1] * y[1][1]),
        (x[1][0] * y[0][0] + x[1][1] * y[1][0], x[1][0] * y[0][1] + x[1][1] * y[1][1]),
    )


def matpow(x: Matrix, n: int) -> Matrix:
    result: Matrix = ((1, 0), (0, 1))
    for _ in range(n):
        result = matmul(result, x)
    return result


def trace(x: Matrix) -> int:
    return x[0][0] + x[1][1]


def determinant(x: Matrix) -> int:
    return x[0][0] * x[1][1] - x[0][1] * x[1][0]


def cohn_matrix(word: str) -> Matrix:
    result: Matrix = ((1, 0), (0, 1))
    for letter in word:
        result = matmul(result, COHN_A if letter == "a" else COHN_B)
    return result


def christoffel_words(depth: int) -> list[str]:
    words: set[str] = set()

    def walk(left: str, right: str, remaining: int) -> None:
        words.add(left + right)
        if remaining:
            walk(left, left + right, remaining - 1)
            walk(left + right, right, remaining - 1)

    walk("a", "b", depth)
    return sorted(words, key=lambda w: (len(w), w))


def fibonacci(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def pell(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, 2 * b + a
    return a


def pell_companion(n: int) -> int:
    c, p = 1, 0
    for _ in range(n):
        c, p = c + 2 * p, c + p
    return c


def branch(prefix: str, tail: str, depth: int) -> list[tuple[int, int, int]]:
    left, right = "a", "b"
    for move in prefix:
        left, right = (left, left + right) if move == "L" else (left + right, right)
    rows: list[tuple[int, int, int]] = []
    for _ in range(depth):
        left, right = (left, left + right) if tail == "L" else (left + right, right)
        matrix = cohn_matrix(left + right)
        m = matrix[0][1]
        u = matrix[0][0] - 2 * m
        rows.append((m, m - u, u))
    return rows


BRANCHES = [
    ("", "L"), ("", "R"), ("L", "R"), ("R", "L"), ("LL", "R"),
    ("RR", "L"), ("LR", "L"), ("RL", "R"), ("LLL", "R"), ("LRL", "R"),
]


def block_cohn_normalisation() -> None:
    words = christoffel_words(TREE_DEPTH)
    for word in words:
        a, m = cohn_matrix(word)[0]
        c, d = cohn_matrix(word)[1]
        assert a * d - m * c == 1, word
        assert a + d == 3 * m, word
        assert 1 <= d <= m, (word, d, m)
        u = a - 2 * m
        assert 0 <= u < m, (word, u, m)
        assert (u * u + 1) % m == 0, (word, u, m)
    print(f"PROVED-CHECK  Cohn normalisation and u = A - 2m on all {len(words)} "
          f"Christoffel words of depth {TREE_DEPTH}")


def block_arm_letter_counts() -> None:
    for k in range(2, 25):
        matrix = matmul(matpow(COHN_A, k - 1), COHN_B)
        assert matpow(FIB_MATRIX, 2) == COHN_A
        assert matrix == matmul(matpow(FIB_MATRIX, 2 * k - 2), COHN_B)
        m = matrix[0][1]
        u = matrix[0][0] - 2 * m
        assert m == fibonacci(2 * k + 1) and u == fibonacci(2 * k - 1), k
        assert m - u == fibonacci(2 * k), k
        # determinant identity used in the proof
        assert fibonacci(2 * k + 1) * fibonacci(2 * k - 1) - fibonacci(2 * k) ** 2 == 1, k
    print("PROVED-CHECK  Fibonacci arm: (m_a, m_b) = (F_2k, F_2k-1), the coordinates "
          "of phi^(2k), for k = 2..24")

    for k in range(1, 25):
        matrix = matmul(COHN_A, matpow(COHN_B, k))
        assert matpow(PELL_MATRIX, 2) == COHN_B
        assert matrix == matmul(COHN_A, matpow(PELL_MATRIX, 2 * k))
        m = matrix[0][1]
        u = matrix[0][0] - 2 * m
        assert m == pell(2 * k + 1) and u == pell(2 * k), k
        assert m - u == pell_companion(2 * k), k
        assert pell(2 * k + 1) * pell(2 * k - 1) - pell(2 * k) ** 2 == 1, k
    print("PROVED-CHECK  Pell arm: (m_a, m_b) = (C_2k, P_2k), the coordinates "
          "of lambda^(2k), for k = 1..24")

    # The coincidence the manuscript records in its remark.
    matrix = matmul(matpow(COHN_A, 6), COHN_B)
    assert matrix[0][1] == 610 and matrix[0][0] - 2 * 610 == 233
    assert 610 - 233 == 377 == fibonacci(14) and 233 == fibonacci(13)
    print("PROVED-CHECK  the Markoff word of 610 has 377 letters a and 233 letters b, "
          "matching phi^14 = 377 phi + 233")


def block_branch_recurrence() -> None:
    for prefix, tail in BRANCHES:
        rows = branch(prefix, tail, 8)
        vectors = [(row[1], row[2]) for row in rows]
        # read m0 off the first valid instance of the recurrence
        m0_num = vectors[2][0] + vectors[0][0]
        m0_den = 3 * vectors[1][0]
        assert m0_num % m0_den == 0, (prefix, tail)
        m0 = m0_num // m0_den
        for k in range(1, len(vectors) - 1):
            expected = (
                3 * m0 * vectors[k][0] - vectors[k - 1][0],
                3 * m0 * vectors[k][1] - vectors[k - 1][1],
            )
            assert expected == vectors[k + 1], (prefix, tail, k)
        print(f"              branch {prefix or '(root)'}|{tail}^inf : m0 = {m0}, "
              f"trace 3*m0 = {3 * m0}")
    print("PROVED-CHECK  the letter-count vector obeys v_(k+1) = 3*m0*v_k - v_(k-1) "
          "on all ten branches")


def markoff_numbers_below(limit: int) -> list[int]:
    """Complete enumeration of the Markoff numbers below `limit` by Vieta descent."""
    triples = {(1, 1, 2)}
    stack = [(1, 1, 2)]
    seen: set[tuple[int, int, int]] = set()
    while stack:
        triple = tuple(sorted(stack.pop()))
        if triple in seen:
            continue
        seen.add(triple)
        a, b, c = triple
        for nxt in ((a, c, 3 * a * c - b), (b, c, 3 * b * c - a)):
            nxt = tuple(sorted(nxt))
            if nxt[2] <= limit and nxt not in seen:
                stack.append(nxt)
                triples.add(nxt)
    return sorted({x for t in triples for x in t})


def block_trace_criterion() -> None:
    complete = markoff_numbers_below(10**60)
    complete_hits = [m for m in complete if math.isqrt(3 * m - 2) ** 2 == 3 * m - 2]
    assert complete_hits == [1, 2, 34], complete_hits
    print(f"VERIFIED      complete enumeration below 10^60: {len(complete)} Markoff numbers, "
          f"3*m0 - 2 a perfect square exactly for {complete_hits}")

    numbers = sorted({1, 2} | {cohn_matrix(w)[0][1] for w in christoffel_words(TREE_DEPTH)})
    hits = []
    for m in numbers:
        # 3m + 2 is never a square: squares are 0 or 1 modulo 3.
        assert (3 * m + 2) % 3 == 2
        root = math.isqrt(3 * m - 2)
        if root * root == 3 * m - 2:
            hits.append((m, root))
    assert [m for m, _ in hits] == [1, 2, 34], hits
    print(f"VERIFIED      among the {len(numbers)} Markoff numbers of depth {TREE_DEPTH}, "
          f"3*m0 - 2 is a perfect square exactly for {[m for m, _ in hits]}")
    print("              (this is a finite scan, not a theorem; the manuscript says so)")


def block_no_half_step_on_displayed_34_ray() -> None:
    points = [(5, 8), (34, 89), (513, 814), (3468, 9077)]
    (x1, y1), (x2, y2), (x3, y3) = points[0], points[1], points[2]
    det = x1 * y2 - y1 * x2
    assert det != 0
    inverse = ((y2, -x2), (-y1, x1))
    image = ((x2, x3), (y2, y3))
    raw = matmul(image, inverse)
    entries = [[Fraction(raw[i][j], det) for j in (0, 1)] for i in (0, 1)]
    assert any(e.denominator != 1 for row in entries for e in row), entries
    assert entries[0][0] == Fraction(-1078, 173)
    assert entries[0][1] == Fraction(1409, 173)
    assert entries[1][0] == Fraction(1409, 173)
    assert entries[1][1] == Fraction(1044, 173)
    nxt = (
        entries[0][0] * x3 + entries[0][1] * y3,
        entries[1][0] * x3 + entries[1][1] * y3,
    )
    assert nxt == (Fraction(593912, 173), Fraction(1572633, 173)), nxt
    assert nxt != (Fraction(points[3][0]), Fraction(points[3][1]))
    print("PROVED-CHECK  the displayed LLL|R^inf fixed-34 ray carries no half-step: "
          "the unique rational matrix matching two steps has denominator 173 "
          "and fails at the third")


def block_growth_rates() -> None:
    left: Matrix = ((1, 1), (0, 1))
    right: Matrix = ((1, 0), (1, 1))
    ones = matmul(left, right)
    twos = matmul(matpow(left, 2), matpow(right, 2))
    mixed = matmul(left, matpow(right, 2))
    assert (trace(ones), determinant(ones)) == (3, 1)
    assert (trace(twos), determinant(twos)) == (6, 1)
    assert (trace(mixed), determinant(mixed)) == (4, 1)
    t_ones = trace(matpow(ones, 6))
    t_twos = trace(matpow(twos, 3))
    t_mixed = trace(matpow(mixed, 4))
    assert (t_ones, t_twos, t_mixed) == (322, 198, 194)
    assert t_ones > t_twos > t_mixed
    print(f"PROVED-CHECK  over twelve letters the Stern-Brocot traces are "
          f"{t_ones} (all ones), {t_twos} (all twos), {t_mixed} (alternating)")
    print("              so the Fibonacci arm is fastest and the Pell arm is not "
          "slowest; no extremality from below is claimed")


def main() -> None:
    block_cohn_normalisation()
    block_arm_letter_counts()
    block_branch_recurrence()
    block_trace_criterion()
    block_no_half_step_on_displayed_34_ray()
    block_growth_rates()
    print(
        "PASS: two-arms certificate; Cohn normalisation, arm letter counts, branch "
        "recurrence with trace 3*m0, trace criterion scan (finite), failure on "
        "the displayed fixed-34 ray, and the growth-rate trace comparison all "
        "check exactly."
    )


if __name__ == "__main__":
    main()
