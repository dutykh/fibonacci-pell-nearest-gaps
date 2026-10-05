#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Enumeration of twins among palindromes, for the conjecture of the closing section.

Every palindrome ``z`` with Gram label ``M(z) = e^T mu(z) e <= --bound`` is
enumerated exactly.  Twins are distinct palindromes with the same label.  A
palindrome is structured around ``p`` when it reads
``p W1 p W2 ... p Wk p Wk^ ... p W1^ p`` with ``k >= 1``, ``p`` a palindrome,
blocks ``Wi`` in ``{ab, ba}`` and ``^`` the block exchange.  The run checks:

* part (i) of the conjecture: every palindrome with a twin is structured;
* part (ii): every class of twins is connected by block exchanges;
* every exchange pair has its root midpoint at its centre, ``M_p`` divides the
  label, the cofactor is prime to ``M_p`` and the roots agree modulo ``M_p``
  (the propositions on arbitrary centres and on the arithmetic of an exchange);
* it reports the labels at which two twins are not an exchange of each other
  (the remark on the midpoint hypothesis);
* the twins of every class have the same numbers of letters ``a`` and ``b``
  (the subsection on lengths of directing words);
* every class of four is a square of nested exchanges ``p<W><V>`` (the corollary
  on the square), and it reports the numbers of blocks of the squares.

With ``--ladder`` it also checks the ladder lemma on every exchange pair.

Run from the article folder, ideally under a memory limit:

    python3 -B supplement/code/twin_census.py [--bound 1e18] [--ladder]

Measured with the system Python 3: ``1e18``, 127,166 palindromes, under a
second; ``1e22``, 1,463,991 palindromes, about a second and 0.3 GB; ``1e26``,
16,743,538 palindromes, 8 s and 2.0 GB (the scope stated in the article).  The
run stops with an INCONCLUSIVE verdict and exit status 2 if more than
``--max-palindromes`` palindromes would be stored; a counterexample to either
part of the conjecture gives exit status 1.  Standard library only.
"""

import argparse
import bisect
import sys
from collections import defaultdict
from fractions import Fraction
from math import gcd

EXCHANGE = {"ab": "ba", "ba": "ab"}


def label_root(word):
    a11, a12, a21, a22 = 1, 0, 0, 1
    for letter in word:
        if letter == "a":
            a11, a12, a21, a22 = 2 * a11 + a12, a11 + a12, 2 * a21 + a22, a21 + a22
        else:
            a11, a12, a21, a22 = 5 * a11 + 2 * a12, 2 * a11 + a12, 5 * a21 + 2 * a22, 2 * a21 + a22
    x, y = 2 * a11 + a12, 2 * a21 + a22
    return 2 * x + y, x


def decode(code):
    kind = code & 3
    code >>= 2
    length = code & 255
    bits = code >> 8
    w = "".join("b" if bits >> i & 1 else "a" for i in range(length))
    return w[::-1] + ("", "a", "b")[kind] + w


def enumerate_labels(bound, cap):
    """Every palindrome w~ c w with label <= bound, as label -> packed codes.

    A palindrome is w~ w or w~ c w; its label is |v|^2 or v^T mu(c) v with
    v = mu(w) e, and v grows when a letter is prepended to w, which prunes.
    """
    first, twins = {}, defaultdict(list)
    stack = [(0, 0, 2, 1)]
    count = 0
    while stack:
        bits, length, x, y = stack.pop()
        even = x * x + y * y
        if even > bound:
            continue
        base = (bits << 8 | length) << 2
        for label, kind in ((even, 0), (2 * x * x + 2 * x * y + y * y, 1), (5 * x * x + 4 * x * y + y * y, 2)):
            if label <= bound:
                count += 1
                if count > cap:
                    return None, None, count
                code = base | kind
                seen = first.get(label)
                if seen is None:
                    first[label] = code
                else:
                    if not twins[label]:
                        twins[label].append(seen)
                    twins[label].append(code)
        stack.append((bits << 1, length + 1, 2 * x + y, x + y))
        stack.append((bits << 1 | 1, length + 1, 5 * x + 2 * y, 2 * x + y))
    return first, twins, count


def parse(word, p):
    if not word.startswith(p):
        return None
    i, blocks = len(p), []
    while i < len(word):
        block = word[i:i + 2]
        if block not in EXCHANGE:
            return None
        blocks.append(block)
        i += 2
        if word[i:i + len(p)] != p:
            return None
        i += len(p)
    return blocks


def structures(word):
    found = []
    for length in range(len(word) - 1):
        p = word[:length]
        if p == p[::-1]:
            blocks = parse(word, p)
            if blocks and len(blocks) >= 2:
                found.append((p, blocks))
    return found


def nested_square(words):
    """(inner, outer) numbers of blocks of a square p<W><V> equal to the four words, or None."""
    target = set(words)
    for w in words:
        for p1, outer in structures(w):
            for p, inner in (structures(p1) if len(p1) > 1 else []):
                p2 = p + "".join(EXCHANGE[b] + p for b in inner)
                hat_outer = [EXCHANGE[b] for b in outer]
                square = {p1 + "".join(b + p1 for b in outer), p1 + "".join(b + p1 for b in hat_outer),
                          p2 + "".join(b + p2 for b in outer), p2 + "".join(b + p2 for b in hat_outer)}
                if square == target:
                    return len(inner), len(outer)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bound", type=float, default=1e18, help="largest label")
    parser.add_argument("--ladder", action="store_true", help="also test the ladder property")
    parser.add_argument("--max-palindromes", type=int, default=20_000_000)
    arguments = parser.parse_args()
    bound = int(arguments.bound)

    first, twins, count = enumerate_labels(bound, arguments.max_palindromes)
    if first is None:
        print(f"INCONCLUSIVE  more than {arguments.max_palindromes} palindromes below {bound:.3e}; nothing checked")
        sys.exit(2)
    sizes = defaultdict(int)
    unstructured, disconnected, pairwise_failures = [], [], []
    exchange_pairs = []
    count_failures, square_failures, square_kinds = [], [], defaultdict(int)
    for label, codes in sorted(twins.items()):
        words = [decode(code) for code in codes]
        sizes[len(words)] += 1
        if len({(w.count("a"), w.count("b")) for w in words}) != 1:
            count_failures.append(label)
        if len(words) == 4:
            kind = nested_square(words)
            if kind is None:
                square_failures.append(label)
            else:
                square_kinds[kind] += 1
        members = set(words)
        neighbours = {w: set() for w in words}
        for w in words:
            found = structures(w)
            if not found:
                unstructured.append((label, w))
            for p, blocks in found:
                partner = p + "".join(EXCHANGE[b] + p for b in blocks)
                assert partner in members, (label, w, p, "the exchange is a twin")
                neighbours[w].add(partner)
                m_p, rho_p = label_root(p)
                m, rho = label_root(w)
                _, rho_partner = label_root(partner)
                assert m == label and m_p * (rho + rho_partner) == 2 * rho_p * m, (label, w, "midpoint")
                assert m % m_p == 0 and gcd(m_p, m // m_p) == 1, (label, w, "divisibility")
                assert (rho - rho_p) % m_p == 0, (label, w, "root congruence")
                if w < partner:
                    exchange_pairs.append((label, p, rho, rho_partner, m_p))
        reached, todo = {words[0]}, [words[0]]
        while todo:
            u = todo.pop()
            for v in neighbours[u] - reached:
                reached.add(v)
                todo.append(v)
        if reached != members:
            disconnected.append((label, words))
        if any(v not in neighbours[u] for u in words for v in words if u != v):
            pairwise_failures.append(label)

    print(f"bound {bound:.3e}: {count} palindromes, {len(first)} labels, {len(twins)} labels with twins,"
          f" group sizes {dict(sorted(sizes.items()))}")
    status = 0
    if unstructured:
        status = 1
        print(f"COUNTEREXAMPLE  part (i): {unstructured[:5]}")
    else:
        print("PASS  part (i): every palindrome with a twin is structured around a palindrome")
    if disconnected:
        status = 1
        print(f"COUNTEREXAMPLE  orbit form: {disconnected[:3]}")
    else:
        print("PASS  part (ii): every class of twins is connected by exchanges")
    print(f"PASS  {len(exchange_pairs)} exchange pairs have their midpoint at the centre,"
          " M_p | M with coprime cofactor, rho = rho_p mod M_p")
    if pairwise_failures:
        print(f"INFO  twins that are not exchanges of each other occur at {len(pairwise_failures)} labels,"
              f" the first being {pairwise_failures[0]}")
    if count_failures:
        status = 1
        print(f"COUNTEREXAMPLE  letter counts differ within a class at labels {count_failures[:5]}")
    else:
        print("PASS  the twins of every class have the same numbers of letters a and b")
    if square_failures:
        status = 1
        print(f"COUNTEREXAMPLE  classes of four that are not squares of nested exchanges: {square_failures[:5]}")
    else:
        print(f"PASS  every class of four is a square of nested exchanges; numbers of blocks (inner, outer)"
              f" {dict(sorted(square_kinds.items()))}")

    if arguments.ladder:
        points = []
        for label, code in first.items():
            m, rho = label_root(decode(code))
            points.append((Fraction(rho, m), m))
        for label, codes in twins.items():
            for code in codes[1:]:
                m, rho = label_root(decode(code))
                points.append((Fraction(rho, m), m))
        points.sort()
        slopes = [s for s, _ in points]
        holds = fails = 0
        for label, p, rho, rho_partner, m_p in exchange_pairs:
            lo, hi = sorted((Fraction(rho, label), Fraction(rho_partner, label)))
            inside = points[bisect.bisect_right(slopes, lo):bisect.bisect_left(slopes, hi)]
            lowest = min(m for _, m in inside)
            if lowest == m_p and sum(1 for _, m in inside if m == m_p) == 1:
                holds += 1
            else:
                fails += 1
                if fails <= 5:
                    print(f"INFO  ladder property fails at label {label}, centre {p!r}: lowest {lowest}, centre {m_p}")
        print(f"INFO  ladder property: holds for {holds} exchange pairs, fails for {fails}")
    sys.exit(status)


if __name__ == "__main__":
    main()
