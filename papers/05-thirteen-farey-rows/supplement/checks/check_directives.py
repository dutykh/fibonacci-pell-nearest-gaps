#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""First verification of the finite facts on directives, and controls of the Gaussian recursion.

The finite facts (the lemma on the directives of Christoffel words of a given length), in exact
integer arithmetic:
1. the Markov numbers m with 5 <= m <= 10^12 are generated from the Markov triples of the Farey
   tree, (m_L, m_R, m_D) -> m_M = 3 m_L m_R - m_D, starting from (1, 2, 1); each must arise once,
   and its root is the residue of -m_L/m_R modulo m reduced below m/2 (the characteristic number);
2. for each of them, every n < m/2 with m/n = [2; B_1, ..., B_j, 2], blocks (2) or (1, 1), is found
   by a search over the blocks from the left: a prefix Y of blocks gives the row (alpha, beta) =
   (2, 1) Y, a tail with column (p, q) gives m = alpha p + beta q, and once alpha (alpha + beta) > m
   at most one tail remains, which is solved for exactly and tested; the root must be the only n;
3. the integers 31 and 1130 have exactly the n listed in the article, by a complete scan of n;
4. the root of 9/10 and the continued fractions and costs quoted for m_{9/10} = 6 625 109;
5. the census written by code/directives_census.jl (the second, independent verification, which
   enumerates every continued fraction [2; B, 2] with numerator at most 10^12 and takes the roots
   from the Cohn matrices) agrees with items 1 to 4 line by line.
Controls of statements proved in the article, which the finite facts do not use:
6. the Gaussian recursion conj(alpha_M) = (m_L + i m_R)/alpha_D along the Farey tree, for every
   reduced P/Q with Q <= 100: exact divisions, norms equal to the labels of the Cohn matrices, and
   alpha_M equal to the occurrence pair of the central Gram factor up to a unit and a conjugation;
7. the square roots of -1: for every integer 5 <= m <= 20000 with a primitive representation, the
   canonical continued fraction of m/u is an even palindrome for each square root u of -1 in
   (1, m/2), its reversed first half is s/r with m = r^2 + s^2, ell(u, m) = 2 ell(r, s) + 1, and
   u -> (r, s) is a bijection onto the primitive representations;
8. the directives: for every 3 <= m <= 300 and every n < m/2 prime to m, the central word of the
   Christoffel word with n letters b and length m is Pal(d) for the word d whose runs are
   c_0 - 1, c_1, ..., c_{k-1}, c_k - 1, and d lies in {ab, ba}* exactly when m/n has the block form.
Standard library only. Exit status 0 when every check passes.
"""
import json
from math import gcd, isqrt
from pathlib import Path
import sys

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O and without arguments.")

HERE = Path(__file__).resolve().parents[1]
CENSUS = HERE / "data" / "directives-census.txt"
BOUND = 10 ** 12


def continued_fraction(p, q):
    out = []
    while q:
        out.append(p // q)
        p, q = q, p % q
    return out


def cost(x, y):
    return sum(continued_fraction(max(x, y), min(x, y))) - 1


def block_form(c):
    """True when c = [2; B_1, ..., B_j, 2] with every block B_i equal to (2) or (1, 1)."""
    if len(c) < 2 or c[0] != 2 or c[-1] != 2:
        return False
    middle, i = c[1:-1], 0
    while i < len(middle):
        if middle[i] == 2:
            i += 1
        elif middle[i] == 1 and i + 1 < len(middle) and middle[i + 1] == 1:
            i += 2
        else:
            return False
    return True


def block_tail(p, q):
    """True when p/q = [B_1, ..., B_j, 2] canonically, with j >= 0 blocks (2) or (1, 1)."""
    if gcd(p, q) != 1:
        return False
    c = continued_fraction(p, q)
    return block_form([2] + c)


def markov_roots(bound):
    """Labels m <= bound of the Farey tree with the residue of -m_L/m_R reduced below m/2."""
    out = {}
    stack = [(1, 2, 1)]
    while stack:
        ml, mr, md = stack.pop()
        m = 3 * ml * mr - md
        if m > bound:
            continue
        a = (-ml * pow(mr, -1, m)) % m
        out.setdefault(m, []).append(min(a, m - a))
        stack.append((ml, m, mr))
        stack.append((m, mr, ml))
    return out


T, A = (2, 1, 1, 0), (2, 1, 1, 1)    # the blocks (2) and (1, 1): P_2 and P_1 P_1 = mu(a)


def search(m):
    """Every n < m/2 with m/n of the block form, by a search over the blocks from the left."""
    found, nodes = [], 0
    stack = [(1, 0, 0, 1)]
    while stack:
        y11, y12, y21, y22 = stack.pop()
        nodes += 1
        al, be = 2 * y11 + y21, 2 * y12 + y22
        if al * (al + be) > m:
            # the tails (p, q) with al p + be q = m and p > q >= 1 have q < m/(al + be) < al
            q = (m * pow(be, -1, al)) % al or al
            while q * (al + be) < m:
                rest = m - be * q
                if rest % al == 0 and block_tail(rest // al, q):
                    found.append(y11 * (rest // al) + y12 * q)
                q += al
            continue
        if 2 * al + be == m:
            found.append(2 * y11 + y12)
        if 5 * al + 2 * be <= m:
            for x in (T, A):
                stack.append((y11 * x[0] + y12 * x[2], y11 * x[1] + y12 * x[3],
                              y21 * x[0] + y22 * x[2], y21 * x[1] + y22 * x[3]))
    return sorted(found), nodes


def scan(m):
    """Every n < m/2 prime to m with m/n of the block form, by a complete scan."""
    return [n for n in range(1, (m + 1) // 2)
            if 2 * n < m and gcd(n, m) == 1 and block_form(continued_fraction(m, n))]


# 1 and 2: the Markov numbers up to 10^12 and their unique integers
roots = markov_roots(BOUND)
labels = sorted(m for m in roots if m >= 5)
assert len(labels) == 150 and all(len(roots[m]) == 1 for m in labels)
nodes_total = 0
found = {}
for m in labels:
    f, nodes = search(m)
    nodes_total += nodes
    assert f == roots[m], (m, f, roots[m])
    found[m] = f

# the search agrees with the complete scan on every integer up to 3000, Markov or not
for m in range(3, 3001):
    assert search(m)[0] == scan(m), m

# 3: the two integers of the article that are not Markov numbers
assert 31 not in roots and 1130 not in roots
assert scan(31) == [12, 13] and continued_fraction(31, 12) == [2, 1, 1, 2, 2]
assert continued_fraction(31, 13) == [2, 2, 1, 1, 2]
sigma = {"a": [1, 1], "b": [2, 2]}
assert scan(1130) == [437, 467]
for n, z in ((437, "abba"), (467, "baab")):
    assert continued_fraction(1130, n) == [2] + [x for c in z for x in sigma[c]] + [2]

# 4: the label of 9/10
MU_A, MU_B = ((2, 1), (1, 1)), ((5, 2), (2, 1))


def mm(x, y):
    return ((x[0][0] * y[0][0] + x[0][1] * y[1][0], x[0][0] * y[0][1] + x[0][1] * y[1][1]),
            (x[1][0] * y[0][0] + x[1][1] * y[1][0], x[1][0] * y[0][1] + x[1][1] * y[1][1]))


def christoffel(p, q):
    return "".join("b" if (i * p) // q - ((i - 1) * p) // q else "a" for i in range(1, q + 1))


def cohn(p, q):
    x = ((1, 0), (0, 1))
    for c in christoffel(p, q):
        x = mm(x, MU_A if c == "a" else MU_B)
    return x


c910 = cohn(9, 10)
m910 = c910[0][1]
u910 = min(c910[0][0] % m910, m910 - c910[0][0] % m910)
assert (m910, u910) == (6625109, 2744210) and roots[m910] == [u910]
assert continued_fraction(m910, u910) == [2] * 18 and cost(u910, m910) == 35
n910 = 2561343
assert gcd(n910, m910) == 1 and 2 * n910 < m910 and cost(n910, m910) == 34
assert continued_fraction(m910, n910) == [2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2,
                                          1, 1, 1, 1, 2, 2]

# 5: the census of the second verification
header, rows = {}, {}
for line in CENSUS.read_text().splitlines():
    if line.startswith("#"):
        key, _, value = line[1:].partition(":")
        header[key.strip()] = value.strip()
        continue
    fields = line.split()
    rows[int(fields[0])] = (int(fields[1]), [int(x) for x in fields[2].split(",")])
assert int(header["bound"]) == BOUND and int(header["Markov numbers from 5"]) == 150
assert sorted(rows) == labels
assert all(rows[m] == (roots[m][0], found[m]) for m in labels)
assert header["31"] == "12,13" and header["1130"] == "437,467"
assert header["9/10"] == f"{m910} {u910} {cost(u910, m910)} {n910} {cost(n910, m910)}"
fractions_enumerated = int(header["continued fractions [2; B, 2]"])


# 6: control of the Gaussian recursion along the Farey tree, Q <= 100
def gauss_div(a, b):
    n = b[0] ** 2 + b[1] ** 2
    x, y = a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1]
    assert x % n == 0 and y % n == 0
    return (x // n, y // n)


def gram_pair(p, q):
    z = christoffel(p, q)[1:-1]
    n = len(z)
    h = ((1, 0), (0, 1))
    for c in z[(n + 1) // 2:]:
        h = mm(h, MU_A if c == "a" else MU_B)
    if n % 2:
        h = mm(((1, 1), (1, 0)) if z[n // 2] == "a" else ((2, 1), (1, 0)), h)
    return (2 * h[0][0] + h[0][1], 2 * h[1][0] + h[1][1])


tree = 0
stack = [((0, 1, 1, (1, 0)), (1, 1, 2, (1, 1)), (1, 0, 1, (1, 0)))]
while stack:
    left, right, diff = stack.pop()
    p, q = left[0] + right[0], left[1] + right[1]
    if q > 100:
        continue
    m = 3 * left[2] * right[2] - diff[2]
    x = cohn(p, q)
    assert x[0][0] + x[1][1] == 3 * m
    conj = gauss_div((left[2], right[2]), diff[3])
    alpha = (conj[0], -conj[1])
    assert alpha[0] ** 2 + alpha[1] ** 2 == m
    assert sorted(map(abs, alpha)) == sorted(gram_pair(p, q))
    mid = (p, q, m, alpha)
    tree += 1
    stack.append((left, mid, right))
    stack.append((mid, right, left))
assert tree == sum(1 for qq in range(2, 101) for pp in range(1, qq) if gcd(pp, qq) == 1)


# 7: control of the lemma on square roots of -1
def fiber(m):
    out = []
    for r in range(1, isqrt(m // 2) + 1):
        s2 = m - r * r
        s = isqrt(s2)
        if s * s == s2 and s > r and gcd(r, s) == 1:
            out.append((r, s))
    return out


def numerator(c):
    p, pp = 1, 0
    for x in c:
        p, pp = x * p + pp, p
    return p


labels_with_fibers = 0
for m in range(5, 20001):
    fib = fiber(m)
    if not fib:
        continue
    labels_with_fibers += 1
    us = [u for u in range(2, (m + 1) // 2) if 2 * u < m and (u * u + 1) % m == 0]
    images = []
    for u in us:
        c = continued_fraction(m, u)
        assert len(c) % 2 == 0 and c == c[::-1]
        half = c[:len(c) // 2][::-1]
        s, r = numerator(half), numerator(half[1:])
        assert s * s + r * r == m and 0 < r < s and gcd(r, s) == 1
        assert cost(u, m) == 2 * cost(r, s) + 1
        images.append((r, s))
    assert sorted(images) == sorted(fib)


# 8: control of the directives of Christoffel words
def closure(w):
    for i in range(len(w) + 1):
        if w[i:] == w[i:][::-1]:
            return w + w[:i][::-1]


def pal(d):
    p = ""
    for x in d:
        p = closure(p + x)
    return p


def runs_word(c):
    lengths = [c[0] - 2] if len(c) == 1 else [c[0] - 1] + c[1:-1] + [c[-1] - 1]
    out, letter = [], "a"
    for length in lengths:
        out.append(letter * length)
        letter = "b" if letter == "a" else "a"
    return "".join(out)


pairs = 0
for m in range(3, 301):
    for n in range(1, (m + 1) // 2):
        if 2 * n >= m or gcd(n, m) != 1:
            continue
        c = continued_fraction(m, n)
        d = runs_word(c)
        assert pal(d) == christoffel(n, m)[1:-1] and len(d) == cost(n, m) - 1
        in_image = len(d) % 2 == 0 and all(d[2 * i] != d[2 * i + 1] for i in range(len(d) // 2))
        assert in_image == block_form(c)
        pairs += 1

print(json.dumps({
    "markov_numbers_from_5_to_10^12": len(labels),
    "search_prefixes": nodes_total,
    "continued_fractions_in_census": fractions_enumerated,
    "non_markov_examples": {"31": [12, 13], "1130": [437, 467]},
    "label_9/10": [m910, u910, n910],
    "gaussian_recursion_fractions": tree,
    "square_root_labels": labels_with_fibers,
    "directive_pairs": pairs,
    "status": "PASS",
}, sort_keys=True))
