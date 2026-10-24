#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Exact Cohn sources, Gram occurrence pairs and Euclidean cost."""
from math import gcd

A = (2, 1, 1, 1)
B = (5, 2, 2, 1)
P1 = (1, 1, 1, 0)
P2 = (2, 1, 1, 0)
IDENTITY = (1, 0, 0, 1)


def multiply(x, y):
    return (x[0]*y[0]+x[1]*y[2], x[0]*y[1]+x[1]*y[3],
            x[2]*y[0]+x[3]*y[2], x[2]*y[1]+x[3]*y[3])


def cohn_word(p, q):
    if not (0 < p < q and gcd(p, q) == 1):
        raise ValueError("require 0 < P < Q and gcd(P,Q)=1")
    return "".join("b" if (j+1)*p//q > j*p//q else "a" for j in range(q))


def word_matrix(word):
    answer = IDENTITY
    for c in word:
        if c not in "ab":
            raise ValueError("binary word required")
        answer = multiply(answer, A if c == "a" else B)
    return answer


def source(p, q):
    word = cohn_word(p, q)
    matrix = word_matrix(word)
    assert word[0] == "a" and word[-1] == "b"
    center = word[1:-1]
    assert center == center[::-1]
    assert matrix[0]+matrix[3] == 3*matrix[1]
    half = center[(len(center)+1)//2:]
    factor = multiply(word_matrix(half), P2)
    if len(center) % 2:
        middle = center[len(center)//2]
        factor = multiply(P1 if middle == "a" else P2, factor)
    pair = tuple(sorted((factor[0], factor[2])))
    assert 0 < pair[0] < pair[1] and gcd(*pair) == 1
    assert pair[0]**2+pair[1]**2 == matrix[1]
    assert euclidean_data(pair)[1] == p+q-2
    return matrix[1], pair


def euclidean_data(pair):
    x, y = sorted(pair)
    if not (0 < x <= y and gcd(x, y) == 1):
        raise ValueError("positive primitive pair required")
    quotients = []
    while x:
        a, r = divmod(y, x)
        quotients.append(a)
        y, x = x, r
    return tuple(quotients), sum(quotients)-1


def fibonacci(n):
    if n < 0:
        raise ValueError("nonnegative index required")
    x, y = 0, 1
    for _ in range(n):
        x, y = y, x+y
    return x
