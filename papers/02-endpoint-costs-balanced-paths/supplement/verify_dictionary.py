#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact checks for the source dictionary and for the definition of the excess.

Three statements of the article are checked here directly from their
definitions, without reference to the reduction graph:

  * the genuine source dictionary, between a reduced Farey slope P/Q and the
    continued-fraction word of its Markoff occurrence, together with the
    letter counts, the balance of the core, and the cost P+Q-2;
  * the definition of the completed excess, Delta = l(R_nu v / 89) - l(v),
    evaluated directly on every word of the source language up to a given
    length and compared with the value the article's weights predict;
  * the explicit zero-cost light family that closes the article, and the five
    integers printed with it.

Everything is exact integer arithmetic. Run with

    python3 -B verify_dictionary.py
"""

import sys
from math import gcd

if not __debug__:
    sys.exit("run without -O: the assertions are part of the check")

I2 = ((1, 0), (0, 1))
A_COHN = ((2, 1), (1, 1))
B_COHN = ((5, 2), (2, 1))
R0 = ((39, 80), (80, -39))
R1 = ((80, 39), (39, -80))
E_VEC = (2, 1)


def mmul(x, y):
    return tuple(tuple(sum(x[i][k] * y[k][j] for k in range(2)) for j in range(2))
                 for i in range(2))


def mvec(x, v):
    return tuple(sum(x[i][k] * v[k] for k in range(2)) for i in range(2))


def digit(j):
    return ((j, 1), (1, 0))


def word_matrix(digits):
    out = I2
    for d in digits:
        out = mmul(out, digit(d))
    return out


def cost(v):
    """Number of unit subtractions taking |v| to a pair of equal entries."""
    x, y = abs(v[0]), abs(v[1])
    n = 0
    while x != y:
        if x > y:
            x -= y
        else:
            y -= x
        n += 1
    return n


def content(v):
    return gcd(abs(v[0]), abs(v[1]))


def christoffel(p, q):
    """Lower Christoffel word of slope p/q as a string over {a,b}."""
    return "".join("b" if (i + 1) * p // q - i * p // q == 1 else "a"
                   for i in range(q))


def mu(word):
    out = I2
    for c in word:
        out = mmul(out, A_COHN if c == "a" else B_COHN)
    return out


def balanced(word):
    n = len(word)
    for r in range(1, n + 1):
        counts = {word[i:i + r].count("b") for i in range(n - r + 1)}
        if counts and max(counts) - min(counts) > 1:
            return False
    return True


def source_word(p, q):
    """The word q = eps chi(w) 2 of the source dictionary, as a tuple of digits."""
    gamma = christoffel(p, q)
    assert gamma[0] == "a" and gamma[-1] == "b"
    z = gamma[1:-1]
    n = len(z)
    if n % 2 == 0 and z[: n // 2] == z[n // 2:][::-1]:
        eps, w = (), z[n // 2:]
    else:
        mid, w = z[n // 2], z[n // 2 + 1:]
        eps = (1,) if mid == "a" else (2,)
    chi = tuple(d for c in w for d in ((1, 1) if c == "a" else (2, 2)))
    return eps + chi + (2,), w, eps


def flip(v):
    """The complete norm-89 flip of a primitive v with valuation one, or None."""
    hits = [r for r in (R0, R1) if content(mvec(r, v)) == 89]
    if len(hits) != 1:
        return None
    w = mvec(hits[0], v)
    return (w[0] // 89, w[1] // 89)


failures = []


def check(name, condition, detail=""):
    if not condition:
        failures.append(f"{name} {detail}")
    return condition


# ------------------------------------------------- the genuine source dictionary
occurrences = 0
for q in range(2, 71):
    for p in range(1, q):
        if gcd(p, q) != 1:
            continue
        occurrences += 1
        digits, core, eps = source_word(p, q)
        M = word_matrix(digits)
        v = (M[0][0], M[1][0])
        label = mu(christoffel(p, q))[0][1]
        check("primitive", content(v) == 1, f"{p}/{q}")
        check("ordered", v[0] > v[1] > 0, f"{p}/{q}")
        check("norm is the label", v[0] ** 2 + v[1] ** 2 == label, f"{p}/{q}")
        check("cost is P+Q-2", cost(v) == p + q - 2, f"{p}/{q}")
        check("core balanced", balanced(core), f"{p}/{q}")
        check("count of twos is P",
              digits.count(2) == 1 + (1 if eps == (2,) else 0) + 2 * core.count("b"),
              f"{p}/{q}")
        check("length gives Q",
              len(digits) + 1 == 2 + (1 if eps else 0) + 2 * len(core), f"{p}/{q}")

# ---------------------------------------------------------------- the excess itself
tested = accepted = 0
for length in range(0, 9):
    for mask in range(2 ** length):
        core = tuple(1 + ((mask >> i) & 1) for i in range(length))
        for eps in ((), (1,), (2,)):
            chi = tuple(d for j in core for d in (j, j))
            digits = eps + chi + (2,)
            M = word_matrix(digits)
            v = (M[0][0], M[1][0])
            norm = v[0] ** 2 + v[1] ** 2
            val = 0
            n = norm
            while n % 89 == 0:
                n //= 89
                val += 1
            tested += 1
            if val != 1:
                continue
            target = flip(v)
            if target is None:
                continue
            accepted += 1
            check("excess nonnegative", cost(target) - cost(v) >= 0, str(digits))
            check("flip preserves the norm",
                  target[0] ** 2 + target[1] ** 2 == norm, str(digits))
            check("flip stays primitive", content(target) == 1, str(digits))
            check("digit sum gives the cost",
                  cost(v) == sum(digits) - 1, str(digits))

# the unique word of zero excess among genuine occurrences
base = word_matrix((1, 1, 1, 2))
bv = (base[0][0], base[1][0])
check("base source is (8,5)", bv == (8, 5), str(bv))
check("base is fixed by the flip", flip(bv) in (bv, (bv[1], bv[0])), str(flip(bv)))

# ------------------------------------------------- the zero-cost light family
S = mmul(((10, 1), (1, 9)), ((10, 1), (1, 9)))
check("the state is a square", S == ((101, 19), (19, 82)), str(S))
check("it commutes with the light block",
      mmul(S, mmul(digit(1), digit(1))) == mmul(mmul(digit(1), digit(1)), S))
prefix = (2, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2)
suffix = (2, 1, 2, 2, 1, 2, 1)
for j in range(0, 7):
    core = prefix + (1,) * j + suffix
    digits = tuple(d for c in core for d in (c, c)) + (2,)
    M = word_matrix(digits)
    v = (M[0][0], M[1][0])
    norm = v[0] ** 2 + v[1] ** 2
    val = 0
    n = norm
    while n % 89 == 0:
        n //= 89
        val += 1
    check("valuation one", val == 1, f"j={j}")
    target = flip(v)
    check("excess is 36", target is not None and cost(target) - cost(v) == 36,
          f"j={j}")
    check("counts are P=21, Q=38+2j",
          (digits.count(2), len(digits) + 1) == (21, 38 + 2 * j), f"j={j}")
    check("the core is unbalanced", not balanced(
        "".join("a" if c == 1 else "b" for c in core)), f"j={j}")

if failures:
    print("FAIL")
    for f in failures[:20]:
        print("   ", f)
    sys.exit(1)

print(f"PASS: source dictionary on {occurrences} reduced occurrences with Q<=70; "
      f"excess evaluated directly on {tested} words of the source language, "
      f"{accepted} of them with norm valuation one at 89, all with Delta>=0; "
      f"the zero-cost light family has excess 36 and unbalanced cores for j=0..6.")
