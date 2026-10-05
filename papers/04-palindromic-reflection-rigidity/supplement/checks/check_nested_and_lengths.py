#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact controls for the cubes of nested exchanges and for the lengths of directing words.

Fixed scopes, chosen to run within the limits of ``run_all.py``:

1. the lemma on the length of a directing word, for every word z of length at most ten:
   |F_a(z)| = |z|_a + 2|z|_b + 1, the regular continued fraction of M/rho, with
   M = e^T mu(z) e and rho = (mu(z) e)_1, is exactly (2, sigma(z), 2), and its digit sum is
   2|F_a(z)| + 2;
2. part (i) of the proposition on Markov words, for the same words: the value of
   a F_a(z) b in the matrices ((1,1),(1,2)) for a and ((2,1),(1,1)) for b, read as the
   top-right entry, equals M(z);
3. part (ii), for every central word z with directive of length at most nine:
   a F_a(z) b is the lower Christoffel word with |z|_b + 1 letters b and |z| + 2 letters a,
   built from its lattice-path definition, M(z) equals the Markov number of that Farey
   fraction computed from the Markov tree of triples, and |F_a(z)| = p + q - 2;
4. the theorem on cubes of nested exchanges, for the centres empty, a, b, aa, aba and abba,
   base block sequences with two or four blocks, and depth at most three (depth three only for
   centres of length at most one): the 2^r words are distinct palindromes with one label and
   one pair of letter counts, the flip of the last coordinates from level i is the exchange
   around the palindrome of level i-1 with the root midpoint of that centre, and, with two
   blocks at every level, the roots are the 2^r patterns of signs modulo 9 M_i^2 + 1;
5. the corollary on the square for the empty centre and W = V = (ba, ab): the four words of
   label 12 986 074 130 and the midpoints 2/5 +- 3390/M of the two diagonals.

All arithmetic is exact.  The last line prints a digest of the scope summary.
"""
import hashlib
import itertools
import sys
from fractions import Fraction as F

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed checker accepts no arguments.")

A = ((2, 1), (1, 1))
B = ((5, 2), (2, 1))
AL = ((1, 1), (1, 2))
BL = ((2, 1), (1, 1))
I2 = ((1, 0), (0, 1))
EXCHANGE = {"ab": "ba", "ba": "ab"}


def mul(x, y):
    return tuple(tuple(sum(x[i][k] * y[k][j] for k in (0, 1)) for j in (0, 1)) for i in (0, 1))


def product(word, a, b):
    out = I2
    for c in word:
        out = mul(out, a if c == "a" else b)
    return out


def label_root(word):
    m = product(word, A, B)
    x, y = 2 * m[0][0] + m[0][1], 2 * m[1][0] + m[1][1]
    return 2 * x + y, x


def f_a(word):
    return "a" + "".join("a" if c == "a" else "ba" for c in word)


def regular_cf(n, d):
    digits = []
    while d:
        q, r = divmod(n, d)
        digits.append(q)
        n, d = d, r
    return digits


def pal(directive):
    word = ""
    for letter in directive:
        word += letter
        for i in range(len(word) + 1):
            if word[i:] == word[i:][::-1]:
                word += word[:i][::-1]
                break
    return word


def lower_christoffel(n, d):
    x = y = 0
    out = []
    while (x, y) != (d, n):
        if y < n and (y + 1) * d <= n * x:
            out.append("b")
            y += 1
        else:
            out.append("a")
            x += 1
    return "".join(out)


def markov_tree(depth):
    table = {(0, 1): 1, (1, 1): 2, (1, 2): 5}
    stack = [((0, 1), (1, 2), (1, 1), 1, 5, 2, 1)]
    while stack:
        left, mid, right, ml, mm, mr, level = stack.pop()
        if level > depth:
            continue
        c1, c2 = (left[0] + mid[0], left[1] + mid[1]), (mid[0] + right[0], mid[1] + right[1])
        m1, m2 = 3 * ml * mm - mr, 3 * mm * mr - ml
        table[c1], table[c2] = m1, m2
        stack.append((left, c1, mid, ml, m1, mm, level + 1))
        stack.append((mid, c2, right, mm, m2, mr, level + 1))
    return table


def angle(q, blocks):
    return q + "".join(b + q for b in blocks)


def parse(word, q):
    if not word.startswith(q):
        return None
    i, blocks = len(q), []
    while i < len(word):
        block = word[i:i + 2]
        if block not in EXCHANGE or word[i + 2:i + 2 + len(q)] != q:
            return None
        blocks.append(block)
        i += 2 + len(q)
    return blocks


summary = []

# 1 and 2: the length lemma and the bridge to the m-value
n1 = 0
for n in range(11):
    for letters in itertools.product("ab", repeat=n):
        z = "".join(letters)
        m, rho = label_root(z)
        sigma = [d for c in z for d in ((1, 1) if c == "a" else (2, 2))]
        digits = regular_cf(m, rho)
        length = len(f_a(z))
        assert length == z.count("a") + 2 * z.count("b") + 1, z
        assert digits == [2] + sigma + [2] and sum(digits) == 2 * length + 2, z
        assert product("a" + f_a(z) + "b", AL, BL)[0][1] == m, z
        n1 += 1
summary.append(f"lengths {n1}")

# 3: central words, Christoffel words and the Markov tree
tree = markov_tree(11)
n3 = 0
for n in range(10):
    for letters in itertools.product("ab", repeat=n):
        z = pal("".join(letters))
        m, _ = label_root(z)
        p, q = z.count("b") + 1, len(z) + 2
        assert "a" + f_a(z) + "b" == lower_christoffel(p, q), z
        assert tree[(p, q)] == m and len(f_a(z)) == p + q - 2, z
        n3 += 1
summary.append(f"central {n3}")

# 4: cubes of nested exchanges
bases = [["ba", "ab"], ["ab", "ba"], ["ab", "ab", "ba", "ba"]]
n4 = 0
for p in ("", "a", "b", "aa", "aba", "abba"):
    m_p, r_p = label_root(p)
    for depth in (1, 2, 3):
        if depth == 3 and len(p) > 1:
            continue
        for choice in itertools.product(range(len(bases)), repeat=depth):
            seqs = [bases[c] for c in choice]
            words, chains = {}, {}
            for s in itertools.product((0, 1), repeat=depth):
                q, chain = p, [p]
                for i, bit in enumerate(s):
                    q = angle(q, seqs[i] if bit == 0 else [EXCHANGE[x] for x in seqs[i]])
                    chain.append(q)
                words[s], chains[s] = q, chain
            values = list(words.values())
            assert len(set(values)) == 2 ** depth and all(w == w[::-1] for w in values)
            data = {s: label_root(w) for s, w in words.items()}
            assert len({d[0] for d in data.values()}) == 1
            assert len({(w.count("a"), w.count("b")) for w in values}) == 1
            m = data[(0,) * depth][0]
            for s in words:
                for i in range(depth):
                    t = tuple(bit ^ (1 if j >= i else 0) for j, bit in enumerate(s))
                    centre = chains[s][i]
                    blocks = parse(words[s], centre)
                    assert blocks is not None and len(blocks) >= 2
                    assert angle(centre, [EXCHANGE[x] for x in blocks]) == words[t]
                    m_c, r_c = label_root(centre)
                    assert m_c * (data[s][1] + data[t][1]) == 2 * m * r_c
            if all(len(x) == 2 for x in seqs):
                levels = [label_root(c)[0] for c in chains[(0,) * depth]]
                patterns = set()
                for s in words:
                    rho = data[s][1]
                    assert (rho - r_p) % m_p == 0
                    signs = []
                    for i in range(depth):
                        sign = (1 if seqs[i][0] == "ba" else -1) * (-1) ** s[i]
                        assert (rho - sign * 3 * levels[i]) % (9 * levels[i] ** 2 + 1) == 0
                        signs.append(sign)
                    patterns.add(tuple(signs))
                assert len(patterns) == 2 ** depth
            n4 += 1
summary.append(f"cubes {n4}")

# 5: the first square
w = ["ba", "ab"]
p1, p2 = angle("", w), angle("", ["ab", "ba"])
z1, z2, z3, z4 = angle(p1, w), angle(p1, ["ab", "ba"]), angle(p2, w), angle(p2, ["ab", "ba"])
assert (z1, z2, z3, z4) == ("baabbabaababbaab", "baababbaabbabaab", "abbabaabbaababba", "abbaababbabaabba")
(m, r1), (_, r2), (_, r3), (_, r4) = (label_root(x) for x in (z1, z2, z3, z4))
assert m == 12986074130 == 1130 * 11492101
assert F(r1 + r3, 2 * m) == F(2, 5) + F(3390, m) and F(r2 + r4, 2 * m) == F(2, 5) - F(3390, m)
assert r1 - r2 == 6 * 1130 and r4 - r3 == r2 - r1
summary.append(f"square {m} {r1} {r2} {r3} {r4}")

text = "; ".join(summary)
print("PASS: lengths of directing words, the bridge to the m-value, nested cubes and the first square (" + text + ")")
print("SHA256 " + hashlib.sha256(text.encode()).hexdigest())
