#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact controls for the remark on chains of exchanges in the closing section.

With u = ab and v = ba, the structured block sequences

    B  = uvuuuvvuuuvuuvvuvvvuuvvvuv,
    B' = uvvuuuvuuuvvuvuuvvvuvvvuuv

of 26 blocks give, around the empty word, two palindromes Z and Z' of 52 letters.  The check
verifies, in exact integer arithmetic:

1. B and B' are structured, Z and Z' are distinct palindromes, and the four palindromes built
   from B, its exchange, B' and its exchange are distinct and share the label
   M = 19 021 048 775 294 072 282 311 853 530 130, which exceeds 10^26;
2. each of the four words is structured around the empty word and around no other palindrome,
   found by testing every palindromic prefix as a centre;
3. the words reached from Z by chains of exchanges, each around a centre of the current word,
   are Z and its exchange only, and the same holds for Z', so no chain joins Z to Z';
4. all four words have 26 letters a and 26 letters b, and the root midpoint of Z and Z' differs
   from the ratio 2/5 of the empty word.

The last line prints a digest of the scope summary.
"""
import hashlib
import sys
from fractions import Fraction as F

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed checker accepts no arguments.")

A = ((2, 1), (1, 1))
B = ((5, 2), (2, 1))
P = ((2, 1), (1, 0))
EXCHANGE = {"ab": "ba", "ba": "ab"}


def mul(x, y):
    return tuple(tuple(sum(x[i][k] * y[k][j] for k in range(2)) for j in range(2)) for i in range(2))


def label_root(word):
    m = ((1, 0), (0, 1))
    for letter in word:
        m = mul(m, A if letter == "a" else B)
    g = mul(mul(P, m), P)
    assert g[0][1] == g[1][0] and g[0][0] * g[1][1] - g[0][1] ** 2 == 1
    return g[0][0], g[0][1]


def angle(q, blocks):
    return q + "".join(b + q for b in blocks)


def structured(blocks):
    n = len(blocks)
    return n >= 2 and n % 2 == 0 and all(blocks[n - 1 - i] == EXCHANGE[blocks[i]] for i in range(n))


def centres(word):
    """Every palindrome p with word = p<C> for a structured block sequence C."""
    found = []
    for length in range(len(word) + 1):
        p = word[:length]
        if p != p[::-1] or (len(word) - length) % (length + 2):
            continue
        blocks, i = [], length
        while i < len(word) and word[i:i + 2] in EXCHANGE and word[i + 2:i + 2 + length] == p:
            blocks.append(word[i:i + 2])
            i += 2 + length
        if i == len(word) and structured(blocks):
            found.append((p, blocks))
    return found


def chain_class(word):
    reached, todo = {word}, [word]
    while todo:
        w = todo.pop()
        for p, blocks in centres(w):
            partner = angle(p, [EXCHANGE[b] for b in blocks])
            if partner not in reached:
                reached.add(partner)
                todo.append(partner)
    return reached


summary = []
code = {"u": "ab", "v": "ba"}
seq = [code[c] for c in "uvuuuvvuuuvuuvvuvvvuuvvvuv"]
seq2 = [code[c] for c in "uvvuuuvuuuvvuvuuvvvuvvvuuv"]
assert len(seq) == len(seq2) == 26 and structured(seq) and structured(seq2)

# 1: four palindromes with one label
words = [angle("", s) for s in (seq, [EXCHANGE[b] for b in seq], seq2, [EXCHANGE[b] for b in seq2])]
z, z2 = words[0], words[2]
assert all(len(w) == 52 and w == w[::-1] for w in words) and len(set(words)) == 4
labels = {label_root(w)[0] for w in words}
assert labels == {19021048775294072282311853530130}
m = labels.pop()
assert m > 10 ** 26
summary.append(f"label {m}")

# 2: the empty word is the only centre of each word
for w in words:
    assert [p for p, _ in centres(w)] == [""]
summary.append("centres empty only")

# 3: the classes under chains of exchanges
assert chain_class(z) == {words[0], words[1]} and chain_class(z2) == {words[2], words[3]}
assert z2 not in chain_class(z)
summary.append("chain classes of size 2, disjoint")

# 4: letter counts and the displaced midpoint
assert all((w.count("a"), w.count("b")) == (26, 26) for w in words)
_, rho = label_root(z)
_, rho2 = label_root(z2)
assert F(rho + rho2, 2 * m) != F(2, 5)
summary.append(f"roots {rho} {rho2}")

text = "; ".join(summary)
print("PASS: two twins of 52 letters around the empty word, joined by no chain of exchanges (" + text + ")")
print("SHA256 " + hashlib.sha256(text.encode()).hexdigest())
