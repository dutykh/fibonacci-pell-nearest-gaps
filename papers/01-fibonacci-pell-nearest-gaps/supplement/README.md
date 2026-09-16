# Reproduction supplement

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

This supplement contains the exact computations behind the
computer-assisted statements of the article. It is self-contained: nothing
below reads a file from outside this directory, and every program uses only
the Python standard library. All decisions that carry part of a proof are made
in exact integer or rational arithmetic, or in outward-rounded fixed-point
arithmetic where an interval is needed.

## Running the checks

From the directory that contains this supplement, run:

```sh
python3 -B supplement/run_all.py
```

The driver stops at the first failure and prints the command and its complete
output. It resolves every executable below `supplement/checks`, so the
manuscript directory can be copied and checked independently of the research
repository. A successful run ends with the line

```text
PASS: every check completed successfully.
```

The complete expected summaries are recorded in `expected-output.txt`.

## Which computation supports which claim

| Statement in the article | Computation | Second, independent computation |
| --- | --- | --- |
| First Matveev comparison and `h < 2.4·10¹⁴ (1 + log(2q))` | Written proof and source constant regression `checks/check_logarithmic_h_bound.py` | Exact rational reconstruction in the logarithmic-bound part of `checks/check_effective_advances.py` |
| Moving coefficient, full norm, second Matveev constants, and `q < 10³²` | `checks/check_adjusted_logarithmic_form.py` | `checks/check_terminal_reduction_second.py` and `checks/check_continued_fraction_closure_second.py` |
| Common rational approximation, `h < 192`, all 190 moving shifts, `q < 90`, and the two-row enumeration | `checks/check_continued_fraction_closure.py` | `checks/check_terminal_reduction_second.py` and `checks/check_continued_fraction_closure_second.py` |
| All-exponent classification and empty odd terminal domain `q = 3, 6, …, 93` | `checks/check_nearest_gap_classification.py --q-max 2000` | `checks/check_nearest_gap_classification_second.py --q-max 2000`; the second checker independently reconstructs the integer arithmetic and does not import the first |
| Exact densities of both, one-sign, and absent candidate windows | Written irrational-rotation proof and `checks/check_window_densities.py --q-max 5000` | The finite counts are regression evidence only; the density theorem is unbounded because of Weyl equidistribution |
| Complete all-anchor two-sign Fibonacci-core orbit list | `checks/check_two_sign_odd_pell_orbit.py` and `checks/check_even_pell_anchor_orbits.py --q-max 5000` | `checks/check_even_pell_anchor_orbits_second.py --q-max 5000`, together with the all-exponent pair; bounded orbit searches are not the completeness proof |
| General-seed semiconjugacy, canonical units, and anchored orbit/gap equivalence | `checks/check_general_seed_orbit.py` | `checks/check_general_seed_orbit_second.py` |
| Uniform simultaneous local clocks, rank-parity criterion, and the obstruction at 241 | Written rank/CRT/equidistribution proof and `checks/check_uniform_local_clocks.py` | The checker verifies exact ranks and witnesses; infinitude rests on the written proof |
| Wider exact nearest-gap regression | `checks/check_even_unit_gap_criterion.py --q-max 2000` | Not used for completeness |
| 7- and 17-clock theorem witnesses and rank data | `checks/check_two_clock_obstruction.py` | `checks/check_two_clock_obstruction_second.py` |
| Two-arms section: Cohn normalisation, arm letter counts, branch recurrence with trace `3m₀`, failure of a half-step on the displayed fixed-34 ray, and the growth-rate trace comparison | `checks/check_two_arms.py` | Self-contained; the two trace-criterion scans are finite verifications, and neither they nor the displayed-ray calculation classify all fixed-34 rays |
| Exact filtering of the published `j = 1` near-square classification | `checks/check_plus_branch_quartic.py` | `checks/check_plus_branch_quartic_second.py` |

Each program can also be run on its own, from this directory; the table above
gives the exact invocation, with its arguments, for every statement. The single
command given under "Running the checks" performs all of them in order.

The independent effective-advances program also reconstructs an earlier
`4 + 4` support exclusion. That extra check is retained, but the present
article uses this program only as an independent check of the first Matveev
comparison.

The wider nearest-gap runs through `q = 2000`, the density run through
`q = 5000`, and the bounded orbit comparisons are regression checks, not the
source of completeness. Completeness comes from the manuscript's unbounded
estimates followed by the exact even domain `q < 90` and odd domain
`q = 3, 6, …, 93`, and from the exact anchored orbit/gap equivalence.

## Optional Magma corroboration

The leading-exponent lemma is proved by the published complete classification
of Alekseyev and Tengely, *Journal of Integer Sequences* 17 (2014), Article
14.6.6, Table 1. Magma is therefore not required for the manuscript proof or
for `run_all.py`.

The `magma` directory preserves two independent Magma V2.29-9 inputs and the
two records containing their complete transcripts,
fixed seeds, version information, integral-point outputs, descent bounds, and
full Mordell--Weil proof flags. See `magma/README.md` for the exact map.

## Arithmetic conventions

- `F(0) = 0`, `F(1) = 1`, and `F(n+2) = F(n+1) + F(n)`.
- `P(0) = 0`, `P(1) = 1`, and `P(n+2) = 2 P(n+1) + P(n)`.
- `C(n) = P(n) + P(n-1)`, with `P(-1) = 1`.
- `λ = 1 + √2`.
- A quadratic integer is represented by its two integer coefficients.
- Every logarithm and square root that decides the effective classification
  is enclosed by rational endpoints or by an outward-rounded fixed-point
  interval. The source first-Matveev regression uses 80-digit `Decimal`
  diagnostics, but its inequalities are separately reconstructed with exact
  rational intervals by the independent effective-advances checker.
- The source local-clock checker uses binary floating point only to locate
  displayed witness candidates and then verifies their nearest windows over
  the integers. The independent local-clock checker reconstructs the displayed
  rows without that locator, while the unbounded clock theorem rests on the
  written recurrence and equidistribution proof rather than on finite search.

## Software

- Python 3.11 or later is recommended; no third-party package is required.
- The optional integral-point reruns use Magma V2.29-9.
- The programs have no nonlocal dependencies, invoke no network service and
  use no nondeterminism. Two of them import named sibling checks in this same
  directory; all the others import only the standard library.

Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon. The package is
distributed under the GNU Lesser General Public License, version 2.1 (see
`LICENSE`); please cite the article when using it (see `CITATION.cff`). The
packages of all the articles of this series are collected at
<https://github.com/dutykh/fibonacci-pell-nearest-gaps/>.
