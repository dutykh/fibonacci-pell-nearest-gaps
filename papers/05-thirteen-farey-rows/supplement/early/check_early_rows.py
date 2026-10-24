#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Fixed source arithmetic and occurrence-cost controls; no parameter scan."""
import sys

if not __debug__:
    raise SystemExit("Python optimization is not supported.")
if len(sys.argv) != 1:
    raise SystemExit("This fixed check accepts no arguments.")

from hashlib import sha256
import json
from math import gcd

IDENTITY = (1, 0, 0, 1)
A = (2, 1, 1, 1)
B = (5, 2, 2, 1)


def mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h


def matrix(word):
    out = IDENTITY
    for letter in word:
        out = mul(out, A if letter == "a" else B)
    return out


def cost(x, y):
    total = 0
    while x != y:
        if x > y:
            q = (x - 1) // y
            x -= q * y
        else:
            q = (y - 1) // x
            y -= q * x
        total += q
    assert x == 1
    return total


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def source(p, q):
    assert 0 < p < q and gcd(p, q) == 1
    word = "".join("b" if i * p // q > (i - 1) * p // q else "a"
                   for i in range(1, q + 1))
    assert word[0] == "a" and word[-1] == "b"
    z = word[1:-1]
    assert z == z[::-1]
    length = len(z)
    w = z[(length + 1) // 2:]
    h = matrix(w)
    if length % 2:
        digit = 1 if z[length // 2] == "a" else 2
        h = mul((digit, 1, 1, 0), h)
    x, y = 2 * h[0] + h[1], 2 * h[2] + h[3]
    m = matrix(word)
    norm = x * x + y * y
    assert norm == m[1] == (m[0] + m[3]) // 3
    assert m[0] + m[3] == 3 * norm
    assert m[0] * m[3] - m[1] * m[2] == 1
    assert gcd(x, y) == 1 and cost(x, y) == p + q - 2
    return [p, q, norm, min(x, y), max(x, y), p + q - 2]


def sign5(a, b):
    """Exact sign of a+b*sqrt(5), for integer coefficients."""
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    delta = a * a - 5 * b * b
    assert delta != 0
    return ((delta > 0) - (delta < 0)) * (1 if a > 0 else -1)


SLOPES = (
    (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
    (2, 3), (2, 5), (3, 4), (7, 8),
    (8, 9), (8, 29), (9, 10), (9, 35),
)
rows = [source(p, q) for p, q in SLOPES]
labels = {(row[0], row[1]): row[2] for row in rows}
assert [row[2] for row in rows[:9]] == [
    5, 13, 34, 89, 233, 29, 194, 169, 195025
]
assert labels[8, 9] == 1136689
assert labels[9, 10] == 6625109
assert labels[9, 35] == 152098457523687709
assert 707 * fib(71) - 1000 * labels[8, 29] == 25054682255203
assert 2929 * fib(85) - 5000 * labels[9, 35] == 155822834239506465
assert 25 * labels[9, 10] < 18 * fib(35)
assert labels[8, 9] < fib(31) and labels[7, 8] < fib(27)
for row in rows[:8]:
    assert row[2] > fib(2 * row[5] + 1)

# Radical comparisons: C7>5/6, C6>1, C8>141/200, and
# 117/200<C9<2929/5000, then 5*phi^-4>29/40.
sign_controls = (
    (4373, -1875),
    (59, -25),
    (71725, -32076),
    (200 * 6561 * 9 - 117 * 625, -200 * 6561 * 4),
    (2929 * 625 - 5000 * 6561 * 9, 5000 * 6561 * 4),
    (671, -300),
)
assert all(sign5(a, b) == 1 for a, b in sign_controls)
payload = {"slopes": rows, "radical_sign_controls": list(sign_controls),
           "scope": "13 fixed slopes; 6 exact radical signs; no range scan"}
canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
print(json.dumps({"status": "PASS", "check": sha256(canonical.encode()).hexdigest(),
                  **payload}, sort_keys=True, indent=2))
