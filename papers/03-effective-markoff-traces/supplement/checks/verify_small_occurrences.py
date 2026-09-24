#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Independent integer-tree and scalar-character verification of one check.

Uses Python 3.10+ standard library only. Imports no primary implementation.
The tree begins at (1,1,1), is traversed by scalar Vieta mutations, and preserves
occurrences. Ray entries use the fundamental second-order recurrence. Literal
word traces use the four-dimensional universal character algebra. Every row of
the bounded domain contributes to an order-independent SHA-256 sum.

Run: python3 -B checks/verify_small_occurrences.py 3
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path

if not __debug__:
    raise SystemExit("Assertions must be enabled; do not use python -O.")


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def checksum(value):
    return int.from_bytes(hashlib.sha256(canonical(value).encode()).digest(), "big")


def integer_nodes(height):
    """Scalar largest-coordinate descent tree; words retain marked homology."""
    stack = [((1, 1, 1), "E", "F", "111")]
    while stack:
        (a, b, p), u, v, path = stack.pop()
        if p > height:
            continue
        assert a <= b <= p and a*a+b*b+p*p == 3*a*b*p
        assert gcd(a, b) == gcd(a, p) == gcd(b, p) == 1
        yield (a, b, p), u, v, path
        if (a, b, p) == (1, 1, 1):
            stack.append(((1, 1, 2), "E", "FE", "112"))
        elif (a, b, p) == (1, 1, 2):
            stack.append(((1, 2, 5), "E", "EFE", "125"))
        else:
            left = (a, p, 3*a*p-b)
            right = (b, p, 3*b*p-a)
            for child in (left, right):
                assert child[2] > p
                assert sorted((child[0], child[1], 3*child[0]*child[1]-child[2])) == [a, b, p]
            stack.append((right, v, u+v, path+"1"))
            stack.append((left, u, u+v, path+"0"))


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in (0, 1)) for j in (0, 1)) for i in (0, 1))


def mi(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    assert det in (-1, 1)
    return ((det*a[1][1], -det*a[0][1]), (-det*a[1][0], det*a[0][0]))


def neg(a):
    return tuple(tuple(-v for v in row) for row in a)


def action(a, v):
    return [sum(a[i][j]*v[j] for j in (0, 1)) for i in (0, 1)]


def normalizer_actions():
    ident = ((1, 0), (0, 1))
    e, f = ((2, 1), (1, 1)), ((2, -1), (-1, 1))
    rot, ref, hyp = ((0, 1), (-1, 1)), ((1, 0), (0, -1)), ((0, -1), (1, 0))
    assert mm(mm(rot, e), mi(rot)) == f
    assert mm(mm(rot, f), mi(rot)) == mi(mm(e, f))
    assert mm(mm(ref, e), mi(ref)) == f and mm(mm(ref, f), mi(ref)) == e
    assert mm(mm(hyp, e), mi(hyp)) == mi(e) and mm(mm(hyp, f), mi(hyp)) == mi(f)
    g, d = ((0, -1), (1, -1)), ((0, 1), (1, 0))
    current, out = ident, []
    for _ in range(3):
        for mat in (current, mm(current, d)):
            out.extend((mat, neg(mat)))
        current = mm(current, g)
    assert current == ident and len(set(out)) == 12
    return out


ACTIONS = normalizer_actions()


def ray(p, x, y, maximum):
    # Fundamental solution, not an iteration starting at the supplied ray.
    fundamental = [0, 1]
    for _ in range(maximum):
        fundamental.append(3*p*fundamental[-1]-fundamental[-2])
    out = {0: x}
    z = 3*p*x-y
    for j in range(1, maximum+1):
        out[j] = fundamental[j]*y-fundamental[j-1]*x
        out[-j] = fundamental[j]*z-fundamental[j-1]*x
        assert out[j] > 0 and out[-j] > 0
    return out


def literal_label(bits, low, high, p):
    tx, ty, txy = 3*low, 3*high, 3*(3*low*high-p)
    a, b, c, d = 1, 0, 0, 0
    for bit in bits:
        if bit:
            a, b, c, d = -c, -d, a+ty*c, b+ty*d
        else:
            a, b, c, d = (-b+(txy-tx*ty)*c-ty*d,
                           a+tx*b+ty*c+txy*d, tx*c+d, -c)
        assert max(abs(v).bit_length() for v in (a, b, c, d)) <= 8192
    raw = 2*a+tx*b+ty*c+txy*d
    assert raw > 2 and raw % 3 == 0
    return raw


def verify(h, reference):
    patterns = {t: tuple(int((i*t)%h >= h-t) for i in range(h))
                for t in range(1, h) if gcd(t, h) == 1}
    matches = []
    for branch, height, indices, qhi in (
        ("positive", 4**(h+1)-1, range(1, 3*h), 3*h//2),
        ("minimum", (100*h*h)**(h+1)-1, (0,), 5*h//2),
    ):
        seeds = comparisons = words = max_bits = seed_sum = comparison_sum = 0
        branch_matches = 0
        # Include the largest comparison index and its recorded neighbor.
        upper = h*max(indices)+(h-1)+qhi+1
        for (x, z, p), word_u, word_v, path in integer_nodes(height):
            alpha = (word_u.count("E"), word_u.count("F"))
            beta0 = (alpha[0]+word_v.count("E"), alpha[1]+word_v.count("F"))
            for orientation, y in ((-1, z), (1, 3*p*x-z)):
                seeds += 1
                assert x <= y and x <= 3*p*x-y
                seed_sum += checksum([branch, path, orientation, p, x, y])
                values = ray(p, x, y, upper)
                beta = (orientation*beta0[0], orientation*beta0[1])
                for l in indices:
                    for t, bits in patterns.items():
                        target = literal_label(bits, values[l], values[l+1], p)
                        words += 1
                        for q in range(-((h+1)//2), qhi+1):
                            m = h*l+t+q
                            if m < 1:
                                continue
                            for n in (-m, m):
                                single = 3*values[n]
                                max_bits = max(max_bits, target.bit_length(), single.bit_length())
                                assert max_bits <= 8192
                                comparisons += 1
                                comparison_sum += checksum([branch, path, orientation, l, t, q, n, target, single])
                                if target != single:
                                    continue
                                v = [alpha[i]+n*beta[i] for i in (0, 1)]
                                w = [h*alpha[i]+(h*l+t)*beta[i] for i in (0, 1)]
                                assert gcd(*v) == gcd(*w) == 1
                                witnesses = [i for i, mat in enumerate(ACTIONS) if action(mat, v) == w]
                                assert witnesses, ("Nonisometric equality", h, p, x, y, l, t, n)
                                original = ((l-n, t, 1) if l-n >= 1 else
                                            (n-l-1, h-t, -1) if n-l-1 >= 1 else (None, None, -1))
                                k, orig_t, step = original
                                matches.append({"branch": branch, "h": h, "path": path,
                                    "orientation": orientation, "seed": [p, x, y], "l": l,
                                    "t": t, "q": q, "n": n, "trace": target, "source": v,
                                    "target": w, "intersection": abs(v[0]*w[1]-v[1]*w[0]),
                                    "original_k": k, "original_t": orig_t, "isometry": witnesses[0],
                                    "character": [single, 3*p, 3*values[n+step]]})
                                branch_matches += 1
        reconstructed = {"seeds": seeds, "literal_words": words,
            "comparisons": comparisons, "equality_rows": branch_matches,
            "max_trace_bits": max_bits, "seed_digest_sum": f"{seed_sum % (1 << 256):064x}",
            "comparison_digest_sum": f"{comparison_sum % (1 << 256):064x}"}
        assert reconstructed == reference["branches"][branch], (branch, reconstructed)
    assert sorted(matches, key=canonical) == reference["equalities"]
    scoped = [row for row in matches if row["original_k"] is not None]
    assert len(scoped) == reference["scoped_equality_rows"]
    assert sorted({(row["original_k"], row["original_t"]) for row in scoped}) == [tuple(row) for row in reference["original_types"]]
    return {"h": h, "status": "PASS: independent complete seed, trace, and actual-isometry reconstruction",
            "comparisons": sum(v["comparisons"] for v in reference["branches"].values()),
            "equality_rows": len(matches), "scoped_equality_rows": len(scoped)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h", type=int, choices=range(2, 7))
    args = parser.parse_args()
    path = Path(__file__).resolve().parents[1] / "data" / f"occurrences-{args.h}.json"
    reference = json.loads(path.read_text())
    assert reference["h"] == args.h
    print(canonical(verify(args.h, reference)))


if __name__ == "__main__":
    main()
