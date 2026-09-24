#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact fixed-occurrence trace comparison on the modular torus.

The analytic bounds are proved in the accompanying article. This module
enumerates their finite domain without identifying occurrences by trace.
Only Python's standard library is required. All matrix arithmetic is integral.
"""

from __future__ import annotations

import argparse
import json
from math import gcd
from typing import NamedTuple

Matrix = tuple[int, int, int, int]
Vector = tuple[int, int]
IDENTITY: Matrix = (1, 0, 0, 1)
E: Matrix = (2, 1, 1, 1)
F: Matrix = (2, -1, -1, 1)
C: Matrix = (5, 2, 2, 1)


def mul(a: Matrix, b: Matrix) -> Matrix:
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def inv(a: Matrix) -> Matrix:
    if a[0]*a[3]-a[1]*a[2] != 1:
        raise ValueError("inverse requires determinant one")
    return (a[3], -a[1], -a[2], a[0])


def pow_mat(a: Matrix, n: int) -> Matrix:
    if n < 0:
        a, n = inv(a), -n
    out = IDENTITY
    while n:
        if n & 1:
            out = mul(out, a)
        a, n = mul(a, a), n // 2
    return out


def tr(a: Matrix) -> int:
    return a[0]+a[3]


def add(a: Vector, b: Vector) -> Vector:
    return (a[0]+b[0], a[1]+b[1])


def linear(a: Vector, r: int, b: Vector, s: int) -> Vector:
    return (r*a[0]+s*b[0], r*a[1]+s*b[1])


def det(a: Vector, b: Vector) -> int:
    return a[0]*b[1]-a[1]*b[0]


def act(a: Matrix, v: Vector) -> Vector:
    return (a[0]*v[0]+a[1]*v[1], a[2]*v[0]+a[3]*v[1])


class Seed(NamedTuple):
    p: int
    x: int
    y: int
    A: Matrix
    B: Matrix
    alpha: Vector
    beta: Vector
    path: str
    orientation: int


def iter_seeds(height: int):
    """Both oriented minimum seeds at every ordered Markov occurrence.

    Paths '111', '112', and '125'+binary suffix identify occurrences.
    The stack is pruned by the maximum, whose strict increase is proved
    in the article. Matrix bases and homology columns travel together.
    """
    if height < 1:
        return
    initial = [(E, F, (1, 0), (0, 1), "111"),
               (E, mul(F, E), (1, 0), (1, 1), "112")]

    def emit(u, v, alpha, gamma, path):
        b = mul(u, v)
        p, x, z = tr(b)//3, tr(u)//3, tr(v)//3
        beta = add(alpha, gamma)
        yield Seed(p, x, z, u, inv(b), alpha,
                   (-beta[0], -beta[1]), path, -1)
        yield Seed(p, x, 3*p*x-z, u, b, alpha, beta, path, 1)

    for u, v, alpha, gamma, path in initial:
        if tr(mul(u, v)) <= 3*height:
            yield from emit(u, v, alpha, gamma, path)
    stack = [(E, C, (1, 0), (2, 1), "125")]
    while stack:
        u, v, alpha, gamma, path = stack.pop()
        b = mul(u, v)
        if tr(b) > 3*height:
            continue
        yield from emit(u, v, alpha, gamma, path)
        beta = add(alpha, gamma)
        stack.append((v, b, gamma, beta, path+"1"))
        stack.append((u, b, alpha, beta, path+"0"))


def mechanical_pattern(h: int, t: int) -> tuple[int, ...]:
    if not (h >= 2 and 0 < t < h and gcd(h, t) == 1):
        raise ValueError("require h >= 2, 0 < t < h, and gcd(h,t)=1")
    return tuple((i*t)//h-((i-1)*t)//h for i in range(1, h+1))


def word_matrix(seed: Seed, h: int, t: int, l: int) -> Matrix:
    low = mul(seed.A, pow_mat(seed.B, l))
    high = mul(low, seed.B)
    out = IDENTITY
    for bit in mechanical_pattern(h, t):
        out = mul(out, high if bit else low)
    return out


def ray_trace(seed: Seed, n: int) -> int:
    return tr(mul(seed.A, pow_mat(seed.B, n)))


def isometry_actions() -> tuple[Matrix, ...]:
    g, d = (0, -1, 1, -1), (0, 1, 1, 0)
    out = []
    for i in range(3):
        gi = pow_mat(g, i)
        for a in (gi, mul(gi, d)):
            out.extend((a, tuple(-x for x in a)))
    return tuple(out)


ISOMETRIES = isometry_actions()


def isometry_witness(v: Vector, w: Vector):
    for index, a in enumerate(ISOMETRIES):
        if act(a, v) == w:
            return index
    return None


def bounds(h: int) -> dict:
    if h < 2:
        raise ValueError("occurrence count must be at least two")
    s = 100*h*h
    jp = 2*h+3
    j0 = (h+1)*(s-1).bit_length()+1
    sp, s0 = jp*(jp-1)+4, j0*(j0-1)+4
    phi = sum(gcd(h, t) == 1 for t in range(1, h))
    np = 2*sp*(3*h-1)*phi*(2*h+1)
    n0 = 2*s0*phi*(3*h+1)
    return {"h": h, "positive_height": 4**(h+1)-1,
            "minimum_height": s**(h+1)-1, "positive_word_length": jp,
            "minimum_word_length": j0, "positive_seed_bound": sp,
            "minimum_seed_bound": s0, "primitive_patterns": phi,
            "positive_comparison_bound": np,
            "minimum_comparison_bound": n0, "comparison_bound": np+n0}


def iter_comparisons(h: int, stats: dict | None = None):
    """Stream every comparison; mutate stats only to record actual work.

    Equality rows outside the original k>=1 family are retained, with
    original_k and original_t set to None. The offset is always m-h*l-t;
    it is never the auxiliary offset after a spectral-state interchange.
    """
    b = bounds(h)
    if stats is None:
        stats = {}
    patterns = [t for t in range(1, h) if gcd(h, t) == 1]
    for branch, height, indices, qhi in (
        ("positive", b["positive_height"], range(1, 3*h), (3*h)//2),
        ("minimum", b["minimum_height"], (0,), (5*h)//2),
    ):
        stats[branch+"_seeds"] = 0
        stats[branch+"_comparisons"] = 0
        for seed in iter_seeds(height):
            stats[branch+"_seeds"] += 1
            ray_cache = {}
            for l in indices:
                for t in patterns:
                    target = tr(word_matrix(seed, h, t, l))
                    for q in range(-((h+1)//2), qhi+1):
                        m = h*l+t+q
                        if m < 1:
                            continue
                        for sign in (-1, 1):
                            n = sign*m
                            if n not in ray_cache:
                                ray_cache[n] = ray_trace(seed, n)
                            stats[branch+"_comparisons"] += 1
                            if ray_cache[n] != target:
                                continue
                            v = linear(seed.alpha, 1, seed.beta, n)
                            w = linear(seed.alpha, h, seed.beta, h*l+t)
                            difference = l-n
                            k = (difference if difference >= 1 else
                                 -difference-1 if difference <= -2 else None)
                            original_t = (t if difference >= 1 else h-t
                                          if difference <= -2 else None)
                            yield {"branch": branch, "h": h, "path": seed.path,
                                   "orientation": seed.orientation,
                                   "seed": [seed.p, seed.x, seed.y],
                                   "l": l, "t": t, "q": q, "n": n,
                                   "trace": target, "source": list(v),
                                   "target": list(w), "intersection": abs(det(v, w)),
                                   "original_k": k, "original_t": original_t,
                                   "isometry": isometry_witness(v, w)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h", type=int)
    parser.add_argument("--run", action="store_true", help="execute after the size estimate")
    parser.add_argument("--max-comparisons", type=int, default=2_000_000)
    args = parser.parse_args()
    b = bounds(args.h)
    print(json.dumps({"kind": "estimate", **b}, sort_keys=True), flush=True)
    if not args.run:
        return
    if b["comparison_bound"] > args.max_comparisons:
        parser.error("proved comparison bound exceeds --max-comparisons")
    stats = {}
    matches = scoped_matches = nonisometric = 0
    for row in iter_comparisons(args.h, stats):
        matches += 1
        if row["original_k"] is not None:
            scoped_matches += 1
            nonisometric += row["isometry"] is None
        print(json.dumps({"kind": "equality", **row}, sort_keys=True))
    print(json.dumps({"kind": "completed", "h": args.h, **stats,
                      "equality_rows": matches, "scoped_equality_rows": scoped_matches,
                      "nonisometric_scoped_rows": nonisometric}, sort_keys=True))


if __name__ == "__main__":
    main()
