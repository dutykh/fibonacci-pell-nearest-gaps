#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Decoding at an arbitrary palindromic centre, for the second conjecture of the closing section.

For a palindrome ``p`` put ``G = mu(p)``, ``M_p = e^T G e`` and
``K_p = 2 e e^T G - M_p I``, a multiple of the reflection ``R_p``.  The script
finds every pair of words ``x, y`` of length at most ``--max-length`` with
``[K_p mu(x) e] = [mu(y) e]`` and checks that ``x`` lies in ``{ab p, ba p}*``
and ``y`` is obtained by exchanging every block, the empty pair included.  Any
other solution is printed and gives exit status 1.

Run from the article folder:

    python3 -B supplement/code/centre_decoder.py [--max-length 12]

The default centres are the nine palindromes of length at most three and abba,
baab, aabaa, ababa, bbabb, abaaba, babbab.  Measured with the system Python 3:
length 12 in well under a second; length 17, the scope stated in the article,
in about 2 s with 0.1 GB.  Memory grows like ``2^(max-length + 1)``.
"""

import argparse
import sys
from math import gcd

EXCHANGE = {"ab": "ba", "ba": "ab"}
DEFAULT_CENTRES = ["", "a", "b", "aa", "bb", "aba", "bab", "aaa", "bbb",
                   "abba", "baab", "aabaa", "ababa", "bbabb", "abaaba", "babbab"]


def direction(u, v):
    g = gcd(u, v)
    u, v = u // g, v // g
    if v < 0 or (v == 0 and u < 0):
        u, v = -u, -v
    return u, v


def columns(max_length):
    """Yield (word, x1, x2) with (x1, x2) = mu(word) e, by prepending letters."""
    stack = [("", 2, 1)]
    while stack:
        word, x1, x2 = stack.pop()
        yield word, x1, x2
        if len(word) < max_length:
            stack.append(("a" + word, 2 * x1 + x2, x1 + x2))
            stack.append(("b" + word, 5 * x1 + 2 * x2, 2 * x1 + x2))


def expected_partner(x, p):
    i, blocks = 0, []
    while i < len(x):
        block = x[i:i + 2]
        if block not in EXCHANGE:
            return None
        blocks.append(block)
        i += 2
        if x[i:i + len(p)] != p:
            return None
        i += len(p)
    return "".join(EXCHANGE[b] + p for b in blocks)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-length", type=int, default=12)
    parser.add_argument("--centres", nargs="*", default=DEFAULT_CENTRES)
    arguments = parser.parse_args()

    table = {}
    for word, x1, x2 in columns(arguments.max_length):
        table.setdefault(direction(x1, x2), []).append(word)
    status = 0
    for p in arguments.centres:
        assert p == p[::-1], p
        g11, g12, g21, g22 = 1, 0, 0, 1
        for letter in p:
            if letter == "a":
                g11, g12, g21, g22 = 2 * g11 + g12, g11 + g12, 2 * g21 + g22, g21 + g22
            else:
                g11, g12, g21, g22 = 5 * g11 + 2 * g12, 2 * g11 + g12, 5 * g21 + 2 * g22, 2 * g21 + g22
        ge1, ge2 = 2 * g11 + g12, 2 * g21 + g22
        m_p = 2 * ge1 + ge2
        k = ((4 * ge1 - m_p, 4 * ge2), (2 * ge1, 2 * ge2 - m_p))
        found = unexpected = 0
        for x, x1, x2 in columns(arguments.max_length):
            u = k[0][0] * x1 + k[0][1] * x2
            v = k[1][0] * x1 + k[1][1] * x2
            if u <= 0 or v <= 0:
                continue
            for y in table.get(direction(u, v), []):
                if expected_partner(x, p) == y:
                    found += 1
                else:
                    unexpected += 1
                    status = 1
                    print(f"COUNTEREXAMPLE  centre {p!r}: x = {x}, y = {y}")
        print(f"{'PASS' if not unexpected else 'FAIL'}  centre {p or '(empty)':8} M_p = {m_p:6}:"
              f" {found} expected pairs, {unexpected} unexpected, words of length <= {arguments.max_length}")
    sys.exit(status)


if __name__ == "__main__":
    main()
