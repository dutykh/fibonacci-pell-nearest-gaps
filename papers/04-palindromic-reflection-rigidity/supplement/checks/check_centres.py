#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact controls for the Thue-Morse coding, the block exchange around an
arbitrary palindromic centre, its arithmetic, the four palindromes of label
12 986 074 130 and the ladder lemma.

Fixed scopes, chosen to run within the limits of ``run_all.py``:

1. the identity ``L_{theta(F_a(z))} = P1 mu(z) P1`` for every word of length at
   most ten, the length and the number of letters ``b`` of
   ``a Pal(theta(F_a(z))) b`` for every palindrome of length at most four (the
   word built by iterated palindromic closure), and the pair of the introduction
   in both codings;
2. ``F_a`` maps the palindromes onto the palindromes beginning with ``a`` and
   avoiding ``bb``, for words of length at most twelve;
3. the Thue-Morse directive equals the doubled directive (odd length) or its
   letter exchange (even length), for every palindrome of length at most ten;
4. for every palindrome ``p`` of length at most five and every block word of at
   most five blocks: equal labels; for palindromic block words the root midpoint,
   the completed counts, the closed form for one pair of blocks, the divisibility
   by ``M_p``, the coprime cofactor and the congruences of the roots;
5. the four palindromes of label 12 986 074 130 and their exchanges;
6. the ladder lemma for every pair of palindromes with equal labels at most
   ``10^9``.

All arithmetic is exact.  The last line prints a digest of the scope summary.
"""
import hashlib
import itertools
import sys
from fractions import Fraction as F
from math import gcd

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed checker accepts no arguments.")

A = ((2, 1), (1, 1))
B = ((5, 2), (2, 1))
P = ((2, 1), (1, 0))
P1 = ((1, 1), (1, 0))
LA = ((1, 1), (0, 1))
LB = ((1, 0), (1, 1))
I2 = ((1, 0), (0, 1))
EXCHANGE = {"ab": "ba", "ba": "ab"}


def mul(x, y):
    return tuple(tuple(sum(x[i][k] * y[k][j] for k in (0, 1)) for j in (0, 1)) for i in (0, 1))


def mu(word):
    out = I2
    for c in word:
        out = mul(out, A if c == "a" else B)
    return out


def incidence(directive):
    out = I2
    for c in directive:
        out = mul(out, LA if c == "a" else LB)
    return out


def gram(word):
    return mul(mul(P, mu(word)), P)


def theta(word):
    return "".join("ab" if c == "a" else "ba" for c in word)


def f_a(word):
    return "a" + "".join("a" if c == "a" else "ba" for c in word)


def exchange(word):
    return word.translate(str.maketrans("ab", "ba"))


def pal(directive):
    word = ""
    for c in directive:
        word += c
        for i in range(len(word) + 1):
            if word[i:] == word[i:][::-1]:
                word += word[:i][::-1]
                break
    return word


def words(n):
    for length in range(n + 1):
        for t in itertools.product("ab", repeat=length):
            yield "".join(t)


def palindromes(n):
    for length in range(n + 1):
        for half in itertools.product("ab", repeat=length // 2):
            h = "".join(half)
            for middle in ([""] if length % 2 == 0 else ["a", "b"]):
                yield h + middle + h[::-1]


def directive_of(q, s):
    out = []
    while (q, s) != (1, 1):
        if q > s:
            out.append("a")
            q -= s
        else:
            out.append("b")
            s -= q
    return "".join(out)


summary = []

# 1. The Thue-Morse coding.
n1 = 0
for z in words(10):
    assert incidence(theta(f_a(z))) == mul(mul(P1, mu(z)), P1), z
    n1 += 1
n1b = 0
for z in palindromes(4):
    g = gram(z)
    omega = "a" + pal(theta(f_a(z))) + "b"
    assert len(omega) == g[0][0] and omega.count("b") == g[1][0], z
    n1b += 1
for z, av, counts in (("baab", "abaaaba", (663, 467)), ("abba", "aababaa", (693, 437))):
    assert f_a(z) == av
    c = incidence(theta(av))
    assert (c[0][0] + c[0][1], c[1][0] + c[1][1]) == counts and sum(counts) == 1130
summary.append(f"coding {n1} {n1b}")

# 2. F_a onto the palindromes beginning with a and avoiding bb.
n2 = 0
for length in range(1, 13):
    targets = {u for u in ("".join(t) for t in itertools.product("ab", repeat=length))
               if u[0] == "a" and "bb" not in u and u == u[::-1]}
    images = {f_a(z) for z in palindromes(length - 1) if len(f_a(z)) == length}
    assert targets == images, length
    n2 += len(targets)
summary.append(f"bijection {n2}")

# 3. The Thue-Morse directive is a doubled directive.
n3 = 0
for z in palindromes(10):
    n = len(z)
    w = z[(n + 1) // 2:]
    t = mul(mu(w), P)
    if n % 2:
        t = mul(P1 if z[n // 2] == "a" else P, t)
    d = directive_of(t[0][0], t[1][0])
    doubled = exchange(d[::-1]) + d
    assert theta(f_a(z)) == (doubled if n % 2 else exchange(doubled)), z
    n3 += 1
summary.append(f"doubled {n3}")

# 4. Block exchange around every palindromic centre and its arithmetic.
n4 = n4p = 0
for p in palindromes(5):
    gp = gram(p)
    m_p, rho_p = gp[0][0], gp[1][0]
    for n in range(6):
        for blocks in itertools.product(("ab", "ba"), repeat=n):
            x = p + "".join(b + p for b in blocks)
            y = p + "".join(EXCHANGE[b] + p for b in blocks)
            gx, gy = gram(x), gram(y)
            m = gx[0][0]
            assert gy[0][0] == m and m % m_p == 0, (p, blocks)
            n4 += 1
            if x != x[::-1]:
                continue
            n4p += 1
            k = n // 2
            rx, ry = gx[1][0], gy[1][0]
            assert m_p * (rx + ry) == 2 * rho_p * m, (p, blocks)
            c = ((2 * k + 1) * (p.count("a") + 1), (2 * k + 1) * (p.count("b") + 1))
            assert (x.count("a") + 1, x.count("b") + 1) == c == (y.count("a") + 1, y.count("b") + 1)
            assert (rx - rho_p) % m_p == 0 and (ry - rho_p) % m_p == 0, (p, blocks)
            if k >= 1:
                assert x != y and gcd(*c) > 1 and gcd(m_p, m // m_p) == 1
                assert (rx + ry) % (m // m_p) == 0
            if k == 1:
                assert m == m_p * (9 * m_p ** 2 + 1)
                assert rx == (9 * m_p ** 2 + 1) * rho_p + (3 * m_p if blocks[0] == "ba" else -3 * m_p)
summary.append(f"exchange {n4} {n4p}")

# 5. The four palindromes of label 12 986 074 130.
Z = ("baab" + "ba" + "baab" + "ab" + "baab", "baab" + "ab" + "baab" + "ba" + "baab",
     "abba" + "ba" + "abba" + "ab" + "abba", "abba" + "ab" + "abba" + "ba" + "abba")
M = 12_986_074_130
assert M == 2 * 5 * 113 * 11_492_101 == 1130 * (9 * 1130 ** 2 + 1)
roots = []
for z in Z:
    g = gram(z)
    assert z == z[::-1] and g[0][0] == M
    roots.append(g[1][0])
assert roots == [5_366_814_557, 5_366_807_777, 5_022_051_527, 5_022_044_747]
blocks = [[z[i:i + 2] for i in range(0, len(z), 2)] for z in Z]
assert all(set(b) <= {"ab", "ba"} for b in blocks)
assert [EXCHANGE[x] for x in blocks[0]] == blocks[3] and [EXCHANGE[x] for x in blocks[1]] == blocks[2]
for i, j, sign in ((0, 2, 1), (1, 3, -1)):
    assert F(roots[i] + roots[j], 2 * M) == F(2, 5) + sign * F(3390, M)
summary.append("quadruple 4")


# 6. The ladder lemma for equal labels up to 10^9.
def all_palindromes(bound):
    stack = [("", 2, 1)]
    while stack:
        w, x, y = stack.pop()
        even = x * x + y * y
        if even > bound:
            continue
        rw = w[::-1]
        for label, z in ((even, rw + w), (2 * x * x + 2 * x * y + y * y, rw + "a" + w),
                         (5 * x * x + 4 * x * y + y * y, rw + "b" + w)):
            if label <= bound:
                yield z
        stack.append(("a" + w, 2 * x + y, x + y))
        stack.append(("b" + w, 5 * x + 2 * y, 2 * x + y))


points = []
for z in all_palindromes(10 ** 9):
    m = mu(z)
    col = (2 * m[0][0] + m[0][1], 2 * m[1][0] + m[1][1])
    points.append((F(col[0], col[1]), 2 * col[0] + col[1], z))
by_label = {}
for slope, label, z in points:
    by_label.setdefault(label, []).append((slope, z))
n6 = 0
for label, group in by_label.items():
    for (s1, z1), (s2, z2) in itertools.combinations(group, 2):
        q = ""
        while len(q) < min(len(z1), len(z2)) and z1[len(q)] == z2[len(q)]:
            q += z1[len(q)]
        lo, hi = min(s1, s2), max(s1, s2)
        inside = [(lab, z) for s, lab, z in points if lo < s < hi]
        assert all(z.startswith(q) for _, z in inside), (z1, z2)
        mq = mu(q)
        label_q = 2 * (2 * mq[0][0] + mq[0][1]) + (2 * mq[1][0] + mq[1][1])
        assert all(lab > label_q for lab, z in inside if z != q), (z1, z2)
        if q == q[::-1]:
            assert (label_q, q) in inside, (z1, z2)
        n6 += 1
summary.append(f"ladders {len(points)} {n6}")

text = "; ".join(summary)
print("PASS: Thue-Morse coding, block exchange at every centre, quadruple and ladders (" + text + ")")
print("SHA256 " + hashlib.sha256(text.encode()).hexdigest())
