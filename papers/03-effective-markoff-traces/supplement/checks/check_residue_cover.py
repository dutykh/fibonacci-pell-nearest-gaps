#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Exact marking-orbit coverage through intersection 24, with the boundary at 25.

This checks units and explicit integral marking maps only: there is no label
census. An arbitrary marking is not asserted to be a metric isometry. The
metric conclusion is supplied by the occurrence theorem, the reflection lemma,
or the separate intersection-nineteen proof.

Python 3.10+, standard library. Run from any directory; --write regenerates data.
"""
from math import gcd
import hashlib
import json
from pathlib import Path
import sys

if not __debug__ or sys.argv[1:] not in ([], ["--write"]):
    raise SystemExit("Run without -O; the only optional argument is --write.")


def act(a, v):
    return [a[0]*v[0]+a[1]*v[1], a[2]*v[0]+a[3]*v[1]]


def det(a):
    return a[0]*a[3]-a[1]*a[2]


def marking_map(r, s, n):
    for sign in (1, -1):
        if (sign*s-r) % n == 0:
            a = [1, (sign*s-r)//n, 0, sign]
            assert act(a, [1, 0]) == [1, 0]
            assert act(a, [r, n]) == [sign*s, sign*n]
            assert abs(det(a)) == 1
            return a
        if (sign-r*s) % n == 0:
            a = [s, (sign-r*s)//n, n, -r]
            assert act(a, [1, 0]) == [s, n]
            assert act(a, [r, n]) == [sign, 0]
            assert abs(det(a)) == 1
            return a
    raise AssertionError("Residues do not have a permitted marking map")


rows = []
for n in range(2, 26):
    units = {r for r in range(1, n) if gcd(r, n) == 1}
    remaining = units.copy()
    orbits = []
    while remaining:
        r = min(remaining)
        orbit = {r, (-r) % n, pow(r, -1, n), (-pow(r, -1, n)) % n}
        assert orbit <= remaining
        maps = {str(s): marking_map(r, s, n) for s in sorted(orbit)}
        if (r*r-1) % n == 0:
            reflection = [r, (1-r*r)//n, n, -r]
            assert det(reflection) == -1
            assert act(reflection, [1, 0]) == [r, n]
            assert act(reflection, [r, n]) == [1, 0]
            for basis in ([1, 0], [0, 1]):
                assert act(reflection, act(reflection, basis)) == basis
            cover = {"method": "reflection", "reflection": reflection}
        else:
            hs = [h for h in range(2, 7) if h in orbit and n//h >= 1 and 0 < n % h < h]
            if hs:
                h = min(hs)
                k, t = divmod(n, h)
                bits = [(i*t)//h-((i-1)*t)//h for i in range(1, h+1)]
                gaps = [k+b for b in bits]
                assert sum(gaps) == n and gcd(h, n) == gcd(h, t) == 1
                # Every cyclic factor has the mechanical balanced count.
                for length in range(1, h):
                    counts = [sum(bits[(start+j) % h] for j in range(length)) for start in range(h)]
                    assert max(counts)-min(counts) <= 1
                cover = {"method": "occurrence", "h": h, "k": k, "t": t,
                         "bits": bits, "B_gaps": gaps}
            elif n == 19:
                assert orbit == {7, 8, 11, 12}
                cover = {"method": "intersection_nineteen", "h": 7, "k": 2, "t": 5}
            elif n == 23:
                assert orbit == {7, 10, 13, 16}
                cover = {"method": "intersection_twentythree", "h": 7, "k": 3, "t": 2}
            else:
                assert n == 25 and orbit in ({7, 18}, {9, 11, 14, 16})
                cover = {"method": "uncovered"}
        orbits.append({"residues": sorted(orbit), "marking_maps": maps, **cover})
        remaining -= orbit
    assert set().union(*(set(item["residues"]) for item in orbits)) == units
    assert sum(len(item["residues"]) for item in orbits) == len(units)
    if n <= 24:
        assert all(item["method"] != "uncovered" for item in orbits)
    else:
        assert [sorted(item["residues"]) for item in orbits if item["method"] == "uncovered"] \
            == [[7, 18], [9, 11, 14, 16]]
    rows.append({"intersection": n, "units": sorted(units), "orbits": orbits})

payload = {"status": "PASS", "scope": "Every primitive unit residue at intersections 2 through 24; exact first remaining orbits at 25",
           "intersection_one": "All bases are in the integral reflection case",
           "rows": rows}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))+"\n"
path = Path(__file__).resolve().parents[1] / "data" / "residue-cover.json"
if sys.argv[1:] == ["--write"]:
    path.write_text(encoded)
else:
    assert json.loads(path.read_text()) == payload
print(json.dumps({"status": "PASS", "unit_count": sum(len(row["units"]) for row in rows),
                  "orbit_count": sum(len(row["orbits"]) for row in rows),
                  "nonreflection_representatives": {row["intersection"]:
                      [min(item["residues"]) for item in row["orbits"] if item["method"] != "reflection"] for row in rows}}, sort_keys=True))
print("reference_sha256="+hashlib.sha256(encoded.encode()).hexdigest())
