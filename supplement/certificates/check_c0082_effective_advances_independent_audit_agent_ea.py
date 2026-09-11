#!/usr/bin/env python3
"""Independent audit of two proposed effective advances beyond C0082.

Questions checked
-----------------
1. For ``attempts/2026-08-24-c0082-effective-44-blocks-agent-eb.md``, derive
   the two Galois permutations from exact arithmetic in
   Q(sqrt(2), sqrt(5)); independently enumerate all unordered 4+4 support
   partitions; impose the real-sign and C0082 exceptional-partition filters;
   and verify the claimed row-span classes, trace obstruction, and elementary
   algebraic inequalities.
2. For ``attempts/2026-08-24-c0082-logarithmic-h-bound-curator.md``, verify
   the Binet signs, algebraic norms, rational endpoint constants, logarithmic
   height budget, Matveev prefactor, and final rounding.

Claim/attempt dependencies: C0068, C0070, C0073, C0080, C0082 and the two
attempts named above. This script does not promote a canonical claim and does
not turn a bounded regression into a universal proof.

Software: Python 3 standard library only. All structural and endpoint checks
use exact integers or ``fractions.Fraction``. Algebraic signs are certified by
rational isolating intervals obtained with ``math.isqrt``. Logarithms are
enclosed by the rational atanh series with an explicit geometric tail.

Reproduce from ``math-sandbox/`` with::

    python3 -B scripts/check_c0082_effective_advances_independent_audit_agent_ea.py

The computation is deterministic and has no random inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isqrt


Q = Fraction


def sqrt_interval(square: int, decimal_places: int = 45) -> tuple[Q, Q]:
    """Return a rigorously certified rational interval for sqrt(square)."""

    scale = 10**decimal_places
    floor = isqrt(square * scale * scale)
    assert floor * floor < square * scale * scale < (floor + 1) ** 2
    return Q(floor, scale), Q(floor + 1, scale)


SQRT2_INTERVAL = sqrt_interval(2)
SQRT3_INTERVAL = sqrt_interval(3)
SQRT5_INTERVAL = sqrt_interval(5)
SQRT10_INTERVAL = (
    SQRT2_INTERVAL[0] * SQRT5_INTERVAL[0],
    SQRT2_INTERVAL[1] * SQRT5_INTERVAL[1],
)


@dataclass(frozen=True)
class K:
    """An exact element a + b sqrt(2) + c sqrt(5) + d sqrt(10)."""

    coefficients: tuple[Q, Q, Q, Q]

    @staticmethod
    def scalar(value: int | Q) -> "K":
        return K((Q(value), Q(0), Q(0), Q(0)))

    def __add__(self, other: "K" | int) -> "K":
        other_k = other if isinstance(other, K) else K.scalar(other)
        return K(tuple(a + b for a, b in zip(self.coefficients, other_k.coefficients)))

    __radd__ = __add__

    def __neg__(self) -> "K":
        return K(tuple(-a for a in self.coefficients))

    def __sub__(self, other: "K" | int) -> "K":
        return self + (-other if isinstance(other, K) else -other)

    def __rsub__(self, other: int) -> "K":
        return K.scalar(other) - self

    def __mul__(self, other: "K" | int | Q) -> "K":
        if not isinstance(other, K):
            scalar = Q(other)
            return K(tuple(scalar * a for a in self.coefficients))
        out = [Q(0), Q(0), Q(0), Q(0)]
        for left_mask, left in enumerate(self.coefficients):
            for right_mask, right in enumerate(other.coefficients):
                overlap = left_mask & right_mask
                factor = (2 if overlap & 1 else 1) * (5 if overlap & 2 else 1)
                out[left_mask ^ right_mask] += factor * left * right
        return K(tuple(out))

    __rmul__ = __mul__

    def divide_scalar(self, value: int | Q) -> "K":
        return self * (Q(1) / Q(value))

    def __pow__(self, exponent: int) -> "K":
        assert exponent >= 0
        result = K.scalar(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def tau2(self) -> "K":
        a, b, c, d = self.coefficients
        return K((a, -b, c, -d))

    def tau5(self) -> "K":
        a, b, c, d = self.coefficients
        return K((a, b, -c, -d))

    def full_norm(self) -> "K":
        return self * self.tau2() * self.tau5() * self.tau2().tau5()

    def interval(self) -> tuple[Q, Q]:
        basis = (
            (Q(1), Q(1)),
            SQRT2_INTERVAL,
            SQRT5_INTERVAL,
            SQRT10_INTERVAL,
        )
        lower = Q(0)
        upper = Q(0)
        for coefficient, (root_lower, root_upper) in zip(self.coefficients, basis):
            if coefficient >= 0:
                lower += coefficient * root_lower
                upper += coefficient * root_upper
            else:
                lower += coefficient * root_upper
                upper += coefficient * root_lower
        return lower, upper

    def sign(self) -> int:
        lower, upper = self.interval()
        if lower > 0:
            return 1
        if upper < 0:
            return -1
        if self == K.scalar(0):
            return 0
        raise AssertionError(f"isolating interval did not determine sign: {self}")


ONE = K.scalar(1)
SQRT2 = K((Q(0), Q(1), Q(0), Q(0)))
SQRT5 = K((Q(0), Q(0), Q(1), Q(0)))
LAMBDA = ONE + SQRT2
LAMBDA_INV = SQRT2 - ONE
PHI = (ONE + SQRT5).divide_scalar(2)
PHI_INV = PHI - ONE
INV_SQRT5 = SQRT5.divide_scalar(5)


@dataclass(frozen=True)
class FormalTerm:
    """A fixed coefficient times lambda^(lm M + lh h) phi^(qexp q)."""

    coefficient: K
    lm: int
    lh: int
    qexp: int

    @property
    def support(self) -> tuple[int, int, int]:
        return self.lm, self.lh, self.qexp

    def tau2(self) -> "FormalTerm":
        # M and h are both odd, and tau2(lambda) = -lambda^(-1).
        parity_sign = -1 if (self.lm + self.lh) % 2 else 1
        return FormalTerm(
            parity_sign * self.coefficient.tau2(),
            -self.lm,
            -self.lh,
            self.qexp,
        )

    def tau5(self) -> "FormalTerm":
        # tau5(phi) = -phi^(-1); all q coefficients here are even.
        assert self.qexp % 2 == 0
        return FormalTerm(self.coefficient.tau5(), self.lm, self.lh, -self.qexp)


def formal_terms(sigma: int, epsilon: int) -> list[FormalTerm]:
    u_plus = sigma * ONE + SQRT2 * PHI
    u_minus = -sigma * ONE + SQRT2 * PHI_INV
    v_plus = sigma * ONE - SQRT2 * PHI
    v_minus = -sigma * ONE - SQRT2 * PHI_INV
    return [
        FormalTerm(K.scalar(5), 1, 0, 0),
        FormalTerm(K.scalar(5), -1, 0, 0),
        FormalTerm(-(u_plus**2), 0, 1, 2),
        FormalTerm(-(u_minus**2), 0, 1, -2),
        FormalTerm(-2 * epsilon * (ONE - sigma * SQRT2), 0, 1, 0),
        FormalTerm(-(v_plus**2), 0, -1, 2),
        FormalTerm(-(v_minus**2), 0, -1, -2),
        FormalTerm(-2 * epsilon * (ONE + sigma * SQRT2), 0, -1, 0),
    ]


def derive_permutation(
    terms: list[FormalTerm], automorphism: str
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Derive, rather than assume, the label permutation and scalar signs."""

    images = [getattr(term, automorphism)() for term in terms]
    permutation: list[int] = []
    signs: list[int] = []
    for transformed in images:
        matches = []
        for label, candidate in enumerate(terms):
            if transformed.support != candidate.support:
                continue
            if transformed.coefficient == candidate.coefficient:
                matches.append((label, 1))
            if transformed.coefficient == -candidate.coefficient:
                matches.append((label, -1))
        assert len(matches) == 1
        label, sign = matches[0]
        permutation.append(label)
        signs.append(sign)
    return tuple(permutation), tuple(signs)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(8))


def apply_permutation(block: frozenset[int], permutation: tuple[int, ...]) -> frozenset[int]:
    return frozenset(permutation[index] for index in block)


def complement(block: frozenset[int]) -> frozenset[int]:
    return frozenset(set(range(8)) - set(block))


def incidence(block: frozenset[int]) -> list[Q]:
    return [Q(int(index in block)) for index in range(8)]


def matrix_rank(rows: list[list[Q]]) -> int:
    matrix = [list(row) for row in rows]
    row = 0
    for column in range(8):
        pivot = next(
            (candidate for candidate in range(row, len(matrix)) if matrix[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        value = matrix[row][column]
        matrix[row] = [entry / value for entry in matrix[row]]
        for other in range(len(matrix)):
            if other == row or not matrix[other][column]:
                continue
            multiplier = matrix[other][column]
            matrix[other] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(matrix[other], matrix[row])
            ]
        row += 1
    return row


def row_span_contains(rows: list[list[Q]], target: list[Q]) -> bool:
    return matrix_rank(rows + [target]) == matrix_rank(rows)


def relation(**entries: int) -> list[Q]:
    row = [Q(0) for _ in range(8)]
    for label, coefficient in entries.items():
        row[int(label)] = Q(coefficient)
    return row


TARGET_F23 = relation(**{"2": 1, "3": -1})
TARGET_F56 = relation(**{"5": 1, "6": -1})
TARGET_A = relation(**{"0": 1, "1": 1, "4": -1, "7": -1})
TARGET_B04 = relation(**{"0": 1, "4": 1, "1": -1, "7": -1})
TARGET_B07 = relation(**{"0": 1, "7": 1, "1": -1, "4": -1})


def block_string(block: frozenset[int]) -> str:
    return "".join(str(index) for index in sorted(block))


def partition_orbit(
    block: frozenset[int], group: tuple[tuple[int, ...], ...]
) -> set[frozenset[int]]:
    images = {apply_permutation(block, permutation) for permutation in group}
    return images | {complement(image) for image in images}


def verify_four_plus_four_classification(
    p2: tuple[int, ...], p5: tuple[int, ...]
) -> None:
    identity = tuple(range(8))
    p25 = compose(p2, p5)
    assert p25 == compose(p5, p2)
    group = (identity, p2, p5, p25)
    full_relation = [Q(1) for _ in range(8)]
    exceptional = {
        frozenset({0, 2, 3, 4}),
        frozenset({0, 5, 6, 7}),
    }
    expected_representatives = {
        1: {
            "F": {"0127", "0157", "0235", "0256"},
            "I": {"0237", "0456"},
            "A": {"0123", "0125", "0126"},
            "B": {"0245", "0246", "0257", "0267"},
        },
        -1: {
            "F": {"0124", "0145", "0235", "0256"},
            "I": {"0237", "0456"},
            "A": {"0123", "0125", "0126"},
            "B": {"0245", "0246", "0257", "0267"},
        },
    }

    for kappa in (1, -1):
        signs = (1, 1, -1, -1, kappa, -1, -1, -kappa)
        candidates: set[frozenset[int]] = set()
        # Choosing label 0 selects exactly one support from every unordered
        # 4+4 partition.
        for other_labels in combinations(range(1, 8), 3):
            block = frozenset((0, *other_labels))
            other = complement(block)
            if len({signs[index] for index in block}) < 2:
                continue
            if len({signs[index] for index in other}) < 2:
                continue
            if block not in exceptional:
                candidates.add(block)
        assert len(candidates) == 28

        classifications: dict[str, set[frozenset[int]]] = {
            "F": set(),
            "I": set(),
            "A": set(),
            "B": set(),
        }
        for block in candidates:
            rows = [full_relation]
            rows.extend(
                incidence(apply_permutation(block, permutation))
                for permutation in group
            )
            has_f = row_span_contains(rows, TARGET_F23) and row_span_contains(
                rows, TARGET_F56
            )
            has_a = row_span_contains(rows, TARGET_A)
            has_b04 = row_span_contains(rows, TARGET_B04)
            has_b07 = row_span_contains(rows, TARGET_B07)
            flags = (has_f, has_a, has_b04, has_b07)
            if has_f:
                assert flags == (True, False, False, False)
                category = "F"
            elif has_a:
                assert flags == (False, True, False, False)
                category = "A"
            elif has_b04 or has_b07:
                assert has_b04 != has_b07
                category = "B"
            else:
                category = "I"
            classifications[category].add(block)

        assert {name: len(blocks) for name, blocks in classifications.items()} == {
            "F": 12,
            "I": 2,
            "A": 6,
            "B": 8,
        }

        representatives: dict[str, set[str]] = {name: set() for name in classifications}
        seen: set[frozenset[int]] = set()
        orbit_sizes: dict[str, list[int]] = {name: [] for name in classifications}
        for category, blocks in classifications.items():
            for block in sorted(blocks, key=block_string):
                if block in seen:
                    continue
                members = blocks & partition_orbit(block, group)
                seen.update(members)
                representative = min(members, key=block_string)
                representatives[category].add(block_string(representative))
                orbit_sizes[category].append(len(members))
        assert seen == candidates
        assert representatives == expected_representatives[kappa]
        assert sorted(orbit_sizes["F"]) == [2, 2, 4, 4]
        assert sorted(orbit_sizes["I"]) == [1, 1]
        assert sorted(orbit_sizes["A"]) == [2, 2, 2]
        assert sorted(orbit_sizes["B"]) == [2, 2, 2, 2]

    # The complete support is fixed by each label permutation, so its orbit
    # contributes only the original relation.
    full_support = frozenset(range(8))
    full_rows = [
        incidence(apply_permutation(full_support, permutation)) for permutation in group
    ]
    assert matrix_rank(full_rows) == 1


def fibonacci(index: int) -> int:
    left, right = 0, 1
    for _ in range(index):
        left, right = right, left + right
    return left


def term_values(q: int, sigma: int, h: int, m: int) -> list[K]:
    """Evaluate the eight C0082 terms exactly for a stress-test row."""

    epsilon = -1 if q % 2 else 1
    M = 2 * m + h
    formal = formal_terms(sigma, epsilon)
    values: list[K] = []
    for term in formal:
        lambda_factor = (
            LAMBDA ** (term.lm * M + term.lh * h)
            if term.lm * M + term.lh * h >= 0
            else LAMBDA_INV ** (-(term.lm * M + term.lh * h))
        )
        phi_exponent = term.qexp * q
        phi_factor = (
            PHI**phi_exponent
            if phi_exponent >= 0
            else PHI_INV ** (-phi_exponent)
        )
        values.append(term.coefficient * lambda_factor * phi_factor)
    return values


def verify_algebraic_identities() -> None:
    assert LAMBDA * LAMBDA_INV == ONE
    assert PHI * PHI_INV == ONE
    assert PHI**2 == PHI + ONE
    assert ONE - SQRT2 == -LAMBDA_INV
    assert ONE + SQRT2 == LAMBDA

    expected_p2 = (1, 0, 5, 6, 7, 2, 3, 4)
    expected_p5 = (0, 1, 3, 2, 4, 6, 5, 7)
    derived_pairs = set()
    for sigma in (1, -1):
        for epsilon in (1, -1):
            terms = formal_terms(sigma, epsilon)
            p2, signs2 = derive_permutation(terms, "tau2")
            p5, signs5 = derive_permutation(terms, "tau5")
            assert p2 == expected_p2
            assert p5 == expected_p5
            assert signs2 == (-1,) * 8
            assert signs5 == (1,) * 8
            kappa = sigma * epsilon
            assert tuple(term.coefficient.sign() for term in terms) == (
                1,
                1,
                -1,
                -1,
                kappa,
                -1,
                -1,
                -kappa,
            )
            derived_pairs.add((p2, p5))
    assert len(derived_pairs) == 1
    p2, p5 = derived_pairs.pop()
    verify_four_plus_four_classification(p2, p5)

    # Class F's exact comparison for sigma=-1.
    u_plus = SQRT2 * PHI - ONE
    u_minus = ONE + SQRT2 * PHI_INV
    assert PHI**2 * u_plus - u_minus == (SQRT2 - ONE) * (PHI + 2)
    assert ((SQRT2 - ONE) * (PHI + 2)).sign() == 1
    for sigma in (1, -1):
        u = sigma * ONE + SQRT2 * PHI
        v = -sigma * ONE + SQRT2 * PHI_INV
        vp = sigma * ONE - SQRT2 * PHI
        vm = -sigma * ONE - SQRT2 * PHI_INV
        for unit in (u, v, vp, vm):
            assert unit.full_norm() == K.scalar(-1)

    # The trace recurrence used in class I is genuinely periodic modulo 5.
    traces = [2, 6]
    for _ in range(12):
        traces.append(6 * traces[-1] - traces[-2])
    residues = [value % 5 for value in traces]
    assert residues[:6] == [2, 1, 4, 3, 4, 1]
    assert residues[6:12] == residues[:6]
    assert 0 not in residues[:6]

    # Check the exact T4-T7 trace formula and both Binet triples on several
    # parity rows. These are identity regressions, not a finite existence test.
    for sigma in (1, -1):
        for epsilon in (1, -1):
            for h in (3, 5, 7, 9):
                t4 = -2 * epsilon * (ONE - sigma * SQRT2) * (LAMBDA**h)
                t7 = -2 * epsilon * (ONE + sigma * SQRT2) * (LAMBDA_INV**h)
                if sigma == 1:
                    expected = 2 * epsilon * (
                        LAMBDA ** (h - 1) + LAMBDA_INV ** (h - 1)
                    )
                else:
                    expected = -2 * epsilon * (
                        LAMBDA ** (h + 1) + LAMBDA_INV ** (h + 1)
                    )
                assert t4 - t7 == expected
                assert (t4 - t7).coefficients[1:] == (Q(0), Q(0), Q(0))
                assert int((t4 - t7).coefficients[0]) % 5 != 0

    for q in range(7, 15):
        epsilon = -1 if q % 2 else 1
        f = fibonacci(q)
        x = fibonacci(q + 1)
        for sigma in (1, -1):
            D = sigma * f * ONE + x * SQRT2
            Dbar = D.tau2()
            values = term_values(q=q, sigma=sigma, h=3, m=1)
            assert sum(values[2:5], K.scalar(0)) == -5 * D**2 * LAMBDA**3
            assert sum(values[5:8], K.scalar(0)) == -5 * Dbar**2 * LAMBDA_INV**3


def log_rational_bounds(value: Q, terms: int = 120) -> tuple[Q, Q]:
    """Enclose log(value) for rational value > 1 by an exact atanh series."""

    assert value > 1
    z = (value - 1) / (value + 1)
    total = Q(0)
    z_power = z
    for index in range(terms):
        total += 2 * z_power / (2 * index + 1)
        z_power *= z * z
    first_omitted_denominator = 2 * terms + 1
    tail = 2 * z_power / (first_omitted_denominator * (1 - z * z))
    return total, total + tail


def log_interval_bounds(interval: tuple[Q, Q]) -> tuple[Q, Q]:
    lower_bound, _ = log_rational_bounds(interval[0])
    _, upper_bound = log_rational_bounds(interval[1])
    return lower_bound, upper_bound


def verify_logarithmic_bound() -> None:
    # Binet formula, signs, norms, and the factor k/D^2 used to control g.
    for sigma in (1, -1):
        u = sigma * ONE + SQRT2 * PHI
        v = -sigma * ONE + SQRT2 * PHI_INV
        assert (u - ONE).sign() == 1
        assert (K.scalar(4) - u).sign() == 1
        absolute_v = v if v.sign() > 0 else -v
        assert (2 * u - absolute_v).sign() == 1
        assert u.full_norm() == K.scalar(-1)
        alpha = (u**2).divide_scalar(5)
        assert (alpha - K.scalar(Q(1, 5))).sign() == 1
        assert (K.scalar(Q(16, 5)) - alpha).sign() == 1

        for q in range(7, 31):
            epsilon = -1 if q % 2 else 1
            f = fibonacci(q)
            x = fibonacci(q + 1)
            D = sigma * f * ONE + x * SQRT2
            binet = (u * PHI**q + epsilon * v * PHI_INV**q) * INV_SQRT5
            assert D == binet
            A = f * ONE + x * SQRT2
            B = x * SQRT2 - f * ONE
            k = 2 * x * x - f * f
            assert A * B == K.scalar(k)
            assert (A - LAMBDA * B).sign() == 1
            assert (3 * B - A).sign() == 1
            assert (3 * D**2 - K.scalar(k)).sign() == 1
            assert (4 * PHI**q - D).sign() == 1

    # Every endpoint appearing before Matveev is checked exactly.
    assert 2 * Q(2, 3) ** 14 < Q(1, 100)
    assert 2 * Q(16, 5) * Q(101, 100) ** 2 < 7
    assert 14 * Q(5, 12) ** 3 == Q(1750, 1728)
    decay_coefficient = Q(201, 100) * Q(1750, 1728) + 3 * Q(101, 100) ** 2
    assert decay_coefficient < Q(128, 25)
    assert Q(128, 25) * Q(5, 12) ** 3 == Q(10, 27)
    assert Q(27, 17) * Q(128, 25) < 9

    assert (LAMBDA - K.scalar(Q(12, 5))).sign() == 1
    assert (3 * ONE - SQRT2 * PHI).sign() == 1
    assert (5 * LAMBDA - 7 * PHI).sign() == 1
    assert Q(7, 5) ** 14 > 32

    # Independent rational enclosures for all logarithms in the height and
    # Matveev prefactor calculations.
    log_lambda = log_interval_bounds(
        (1 + SQRT2_INTERVAL[0], 1 + SQRT2_INTERVAL[1])
    )
    log_phi = log_interval_bounds(
        ((1 + SQRT5_INTERVAL[0]) / 2, (1 + SQRT5_INTERVAL[1]) / 2)
    )
    log2 = log_rational_bounds(Q(2))
    log4 = log_rational_bounds(Q(4))
    log5 = log_rational_bounds(Q(5))
    log9 = log_rational_bounds(Q(9))

    assert log_lambda[0] > Q(4, 5)
    assert 3 * log2[1] + log_phi[1] + log5[1] < 5
    assert log9[1] < 3

    # 3^(9/2) = 81 sqrt(3). Use upper endpoints throughout, so this is a
    # rigorous upper enclosure rather than a floating-point recomputation.
    matveev_prefactor_upper = (
        Q(7, 5)
        * 30**6
        * 81
        * SQRT3_INTERVAL[1]
        * 4**2
        * (1 + log4[1])
        * (2 * log_lambda[1])
        * (2 * log_phi[1])
    )
    assert matveev_prefactor_upper < 9_275_000_000_000
    assert 20 * matveev_prefactor_upper < 185_500_000_000_000

    # From h log(lambda) < 1.855e14 H + log(9), H >= 1 and
    # log(lambda) > 4/5, certify the advertised upward rounding.
    assert 185_500_000_000_000 + 3 < Q(4, 5) * 240_000_000_000_000


def main() -> None:
    verify_algebraic_identities()
    verify_logarithmic_bound()
    print(
        "PASS: independently derived both Galois actions; exhaustively split "
        "each 28-pattern 4+4 sector as 12+2+6+8 with the claimed row-span "
        "consequences; certified the mod-5 trace and size obstructions; and "
        "verified the Binet signs, norms, rational constants, height budget, "
        "Matveev prefactor, coefficient bound n<2q, and final h-bound rounding."
    )


if __name__ == "__main__":
    main()
