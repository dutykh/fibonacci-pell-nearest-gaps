# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Fixed examples and exact polynomial identities used by the analytic proofs."""
from math import gcd
import elementary as h


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def mul(*polynomials):
    result = {(0, 0, 0): 1}
    for polynomial in polynomials:
        product = {}
        for left, a in result.items():
            for right, b in polynomial.items():
                exponent = tuple(x + y for x, y in zip(left, right))
                product[exponent] = product.get(exponent, 0) + a * b
        result = {exponent: coefficient for exponent, coefficient in product.items() if coefficient}
    return result


def scale(polynomial, coefficient):
    return {exponent: value * coefficient for exponent, value in polynomial.items() if value * coefficient}


def markoff(a, b, c):
    return add(mul(a, a), mul(b, b), mul(c, c), scale(mul(a, b, c), -3))


def check():
    # A generic polynomial identity, not sampling of the recurrence parameter.
    a, b, c = ({(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1})
    next_c = add(scale(mul(b, c), 3), scale(a, -1))
    previous_a = add(scale(mul(a, b), 3), scale(c, -1))
    assert markoff(b, c, next_c) == markoff(a, b, c)
    assert add(scale(mul(b, c), 3), scale(next_c, -1)) == a
    assert add(scale(mul(a, b), 3), scale(previous_a, -1)) == c
    initial = (1, 34, 89)
    assert sum(x * x for x in initial) == 3 * initial[0] * initial[1] * initial[2]
    q4, q5 = h.word_matrix((1, 1, 2)), h.word_matrix((1, 1, 1, 2))
    assert (q4[0], q4[2]) == (5, 3) and q4[0] ** 2 + q4[2] ** 2 == 34
    assert (q5[0], q5[2]) == (8, 5) and q5[0] ** 2 + q5[2] ** 2 == 89
    assert 1 * 5 - 1 * 4 == 1
    # In Q(phi), these linear identities prove 1/5 < rho < 1/4.
    numerator, denominator = (1, 1), (4, 5)
    assert tuple(5 * x - y for x, y in zip(numerator, denominator)) == (1, 0)
    assert tuple(y - 4 * x for x, y in zip(numerator, denominator)) == (0, 1)

    # One literal pumped family; its all-length claim uses this exact loop.
    reflection89 = (39, 80, 80, -39)
    entry = tuple(map(int, "1221122112"))
    exit_word = tuple(map(int, "2122121"))
    state, emitted = h.batch_reduce(h.mul(reflection89, h.DOUBLE[2]))
    total = emitted - 4
    for letter in entry:
        edge = h.step(state, letter)
        state, total = edge.target, total + edge.weight
    prescribed = (101, 19, 19, 82)
    assert state in (prescribed, (prescribed[2], prescribed[3], prescribed[0], prescribed[1]))
    loop = h.step(state, 1)
    assert loop.target == state and loop.weight == 0 and loop.emitted == 2
    for letter in exit_word:
        edge = h.step(state, letter)
        state, total = edge.target, total + edge.weight
    assert gcd(*h.terminal_vector(state)) == 89
    assert total + h.terminal_weight(state) == 36
    prefix = (2,) + entry
    assert prefix.count(2) + exit_word.count(2) == 10
    assert len(prefix) + len(exit_word) == 18
    return {"recurrence": "F(a,b,c)=(b,c,3bc-a)", "polynomial_invariance": True,
            "polynomial_inverse": True, "initial_markoff_triple": initial,
            "modulus": 89**2, "finite_permutation_period_bound": 89**6,
            "limit_slope": "(1+phi)/(4+5phi)",
            "limit_interval": ["1/5", "1/4"],
            "zero_light_family": {"entry_prefix": prefix, "exit": exit_word,
                                  "state": prescribed, "gap": 36,
                                  "P": 21, "Q": "38+2j"}}
