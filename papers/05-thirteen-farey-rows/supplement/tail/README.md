# Exact unbounded tail checks

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry, France)

This package is the first verification of the tail covers of the article
(parts (ii) and (iii) of its proposition on the finite verifications); the
second, independent one is `../checks/independent_tail_cover.py`. Its finite
trees cover **all** eligible nonnegative integer gap tuples, including
unbounded gaps. It is not a finite-cost sample and does not use a list of
Markov values. Source bands and the six row-tail thresholds are fixed in the
public verifier.

The proof uses exact coordinatewise parity order of the full normalized
ratio. Each parity box is bounded by its two matched corners, where an
infinite gap means the exact limit q^g=0. The entire numerator and
Fibonacci denominator are evaluated at the same corner. A finite tree
closes only through strict exact inequalities or an ineligible-total
check; unresolved singletons and resource exhaustion fail.

## Files and dependencies

- `exact.py`: arithmetic in Z[phi], channel evaluation, source constants and
  seven prescribed normalization controls against physical integer paths.
- `corner_cover.py`: one parity-sector cover, with explicit node and gap
  endpoint caps.
- `sweep.py`: complete parity sweep for one row/defect group; its `--verify`
  mode reconstructs every tree and compares all saved records.
- `verify_all.py`: full fixed-domain replay of all 28 groups, 1976 sectors
  and 1,023,768 boxes, including the separate lower-cost comparison.
- `data/rowP-kK.jsonl`: actual sector records. Each record gives its domain,
  counts, resource diagnostics and a streaming digest of the terminal boxes.
- `data/lower-cost-two-defects.jsonl`: the eight exact sector records
  proving `N/F_(2c+3) < 71/100` for two defects and every `c ⩾ 155`.
- `data/catalog.json`: compact readable inventory of the fixed domains and
  their complete sector manifests, with the group, sector and box counts of
  every cover. This is a catalog, not a proof shortcut.

Only Python 3.10 or newer and its standard library are needed. No package
from outside this directory is imported. The generators in this package are
inverse-Euclidean path maps, not the Cohn trace generators. Gaps are applied
chronologically to (1,2), and the cost is the gap sum plus the number of
gaps. The normalized value is N/F_(2c+3).

## Reproduce the proof

From the manuscript directory, first inspect the fixed domain:

```sh
python3 -B supplement/tail/verify_all.py
```

On Linux, a capped complete replay is:

```sh
timeout 600s prlimit --as=1073741824 --cpu=300 \
  python3 -B supplement/tail/verify_all.py --run
```

The mathematical verification is portable Python; `timeout` and `prlimit`
provide the optional external resource wrapper on Linux. On an idle machine
the replay takes about ten seconds and less than 64 MiB; the generous caps
keep a busy machine from turning a correct replay into a failure. For a single group, for example:

```sh
python3 -B supplement/tail/sweep.py --row 12 --defects 5 \
  --minimum-cost 110 --verify supplement/tail/data/row12-k5.jsonl
```

The verifier recomputes every sign and reconstructs the exhaustive cover.
Checking a saved hash or reading summary counts is never enough for PASS.
An absent, altered, incomplete or extra sector record fails. Ordinary
verification writes no files. Assertions must remain enabled.

The resulting theorem is limited to the stated source bands and cost
thresholds. These checks do not prove Gaussian-fibre completeness,
the bounds on the source labels of a row, an extremizer classification at
smaller costs, or global Markov uniqueness; those are separate mathematical
inputs or open questions.

## The size of the trees

These counts describe the finite trees; they are not part of the argument,
whose completeness comes from the partition of the eligible gap tuples into
root boxes, from the partition of each box at every split, and from the
closing of every leaf. Since every split divides one box into two, each tree
satisfies `visited boxes = 2 × terminal boxes − root boxes`, a convenient
check on the bookkeeping. The middle defect counts are those of the row
thresholds; the final defect count `P − 5` closes already at its root boxes,
all below the source band.

| Row `P` | Middle defects `k` | Visited boxes | Terminal boxes | Final `k` | Its root boxes |
| --- | --- | ---: | ---: | --- | ---: |
| 8 | 2 | 40 | 32 | 3 | 64 |
| 9 | 2 to 3 | 256 | 172 | 4 | 160 |
| 10 | 2 to 4 | 1,188 | 718 | 5 | 384 |
| 11 | 2 to 5 | 12,588 | 6,610 | 6 | 896 |
| 12 | 2 to 6 | 123,236 | 62,382 | 7 | 2,048 |
| 13 | 2 to 7 | 878,276 | 440,926 | 8 | 4,608 |
| Total | | 1,015,584 | 510,840 | | 8,160 |

Together these covers comprise 1,968 parity sectors and 14,256 root boxes;
the largest single stratum is row 13 with six defects. The separate
lower-cost cover for two defects and costs at least 155 adds eight sectors
and 24 root boxes, all closing at the root, for 1,976 sectors and 14,280 root
boxes in all. The order of splitting matters for termination: splitting at
the first non-singleton coordinate instead of the one with the smallest lower
endpoint does not terminate already at row 9, because it keeps refining a
coordinate whose contribution is nearly settled.
