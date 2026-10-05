#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact words and Gram matrices for the palindromic reflection paper.

Python 3.10+, standard library. Words are literal finite strings over a,b.
Matrix tuples are row-major, and reflection entries are rational Fractions.
This is an example/reference implementation, not a finite proof of a language
classification. Run the fixed check suite described in ../README.md.
"""
from fractions import Fraction

A = (2, 1, 1, 1)
B = (5, 2, 2, 1)
P = (2, 1, 1, 0)
IDENTITY = (1, 0, 0, 1)
E = (2, 1)


def mul(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def transpose(a):
    return (a[0], a[2], a[1], a[3])


def det(a):
    return a[0]*a[3]-a[1]*a[2]


def tr(a):
    return a[0]+a[3]


def inv(a):
    d = det(a)
    if not d:
        raise ValueError("singular matrix")
    return tuple(Fraction(x, d) for x in (a[3], -a[1], -a[2], a[0]))


def act(a, v):
    return (a[0]*v[0]+a[1]*v[1], a[2]*v[0]+a[3]*v[1])


def check_word(word):
    if not isinstance(word, str) or any(c not in "ab" for c in word):
        raise ValueError("a word must be a string over a,b")


def word_matrix(word):
    check_word(word)
    out = IDENTITY
    for letter in word:
        out = mul(out, A if letter == "a" else B)
    return out


def palindrome_closure(word):
    """Shortest palindrome with the given word as prefix."""
    check_word(word)
    for i in range(len(word)+1):
        tail = word[i:]
        if tail == tail[::-1]:
            return word+word[:i][::-1]
    raise RuntimeError("the empty suffix is palindromic")


def pal(directive):
    check_word(directive)
    out = ""
    for letter in directive:
        out = palindrome_closure(out+letter)
    return out


def right_images(directive):
    check_word(directive)
    u, v = "a", "b"
    for letter in directive:
        if letter == "a":
            v = v+u
        else:
            u = u+v
    return u, v


def right_image(directive, word):
    check_word(word)
    u, v = right_images(directive)
    return "".join(u if c == "a" else v for c in word)


def standard_image(directive, palindrome):
    check_word(palindrome)
    if palindrome != palindrome[::-1]:
        raise ValueError("the input must be palindromic")
    return pal(directive)+right_image(directive, palindrome)


def gram(palindrome):
    check_word(palindrome)
    if palindrome != palindrome[::-1]:
        raise ValueError("the input must be palindromic")
    return mul(mul(P, word_matrix(palindrome)), P)


def gram_data(palindrome):
    g = gram(palindrome)
    return {"M": g[0], "rho": g[1], "H": g[3]}


def completed_counts(word):
    check_word(word)
    return (word.count("a")+1, word.count("b")+1)


def incidence(directive):
    u, v = right_images(directive)
    return (u.count("a"), v.count("a"), u.count("b"), v.count("b"))


def token_exchange(word):
    check_word(word)
    if len(word) % 2 or any(word[i:i+2] not in ("ab", "ba")
                            for i in range(0, len(word), 2)):
        raise ValueError("the word must be a concatenation of ab/ba tokens")
    return "".join("ba" if word[i:i+2] == "ab" else "ab"
                   for i in range(0, len(word), 2))


def paired_preimages(half_word):
    partner = token_exchange(half_word)
    return (half_word[::-1]+half_word, partner[::-1]+partner)


def ancestor_reflection(directive):
    g = word_matrix(pal(directive))
    w = act(g, E)
    m = E[0]*w[0]+E[1]*w[1]
    k = (2*E[0]*w[0]-m, 2*E[0]*w[1],
         2*E[1]*w[0], 2*E[1]*w[1]-m)
    return tuple(Fraction(x, m) for x in k)
