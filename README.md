<!--
Fibonacci–Pell nearest gaps: an all-exponent classification and
quadratic-unit orbit rigidity

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Fibonacci–Pell nearest gaps

**An all-exponent classification and quadratic-unit orbit rigidity**

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

*MSC 2020*: Primary 11B39, 11D61; Secondary 11J86, 11D09.
*Keywords*: Fibonacci numbers, Pell equations, quadratic units, linear forms in
logarithms, `p`-adic valuations, arithmetic dynamics, Markoff uniqueness
conjecture.

This repository holds the complete LaTeX source of the manuscript, the compiled
PDF, and a reproduction supplement of twenty-one exact certificates that can be
re-run with nothing beyond a Python 3 installation.

---

## What is proved

Let `λ = 1 + √2`, let `(F_q)` be the Fibonacci sequence, and put
`D_{q,σ} = σ F_q + F_{q+1} √2` for `σ ∈ {1, −1}`. Suppose a positive integer `n`
lies in the strict factor-two window `D_{q,σ}² < λⁿ < 2 D_{q,σ}²`, write
`Δ = λⁿ − D_{q,σ}²` for the resulting gap, and let `g` be the gcd of its two
integer coefficients.

**All-exponent nearest-gap classification.** For every `q ⩾ 1`,
every sign, and every positive exponent `n`, the window condition together with
the exact norm–content condition `N(Δ) = −g²` holds in exactly two cases:

| `q` | `σ` | `n` | `g` | `Δ/g` |
| --- | --- | --- | --- | --- |
| 2 | −1 | 2 | 6 | `λ⁻¹` |
| 4 | +1 | 6 | 40 | `λ` |

In particular no solution has an odd exponent. By contrast the underlying
Archimedean windows are abundant: they occur with natural density
`log_λ √2 − η = 0.3300…`, governed by an irrational circle rotation, and that
density survives restriction to any arithmetic progression. Scarcity of windows
is therefore not the exclusion mechanism; the exact norm–content equation is.

**Complete Fibonacci-core Pell-orbit classification.** The map
`S(u,v) = (v, u + 2v)` is integrally conjugate to multiplication by `λ` through
`A_(u,v) = (v − u) + u √2`. For the two Fibonacci cores
`z_{q,σ} = (F_{q+1}, F_{q+1} + σ F_q)`, every hit
`H(Sᵗ z_{q,σ}) = P_{d+2t}` at a positive anchor `d` and time `t ⩾ 0` appears in
exactly four rows: one canonical persistent family `(σ, q, d) = (1, 1, 3)` valid
for all `t`, and three isolated time-zero hits at `(−1, 1, 1)`, `(−1, 2, 3)` and
`(−1, 4, 5)`. There is no hit at a positive even anchor.

**Why `φ` and `λ`, and how far that is proved.** Every branch of the Markoff
tree that keeps one entry `m₀` fixed carries an integral Cohn monodromy whose
dominant eigenvalue is the quadratic Pisot unit solving `X² − 3m₀X + 1 = 0`.
That eigenvalue sets the scale on every branch alike, so on its own it
distinguishes none of them. What the orbit bridge consumes is stronger: an
integral square root of that unit, acting compatibly on two-square data.
Comparing traces shows that any such half-step forces `3m₀ − 2` to be a perfect
square. The Fibonacci and Pell arms, `m₀ = 1` and `m₀ = 2`, carrying `φ` and
`λ`, supply the two explicit half-steps used here, and they are also the two
constant-directive branches. The square condition is necessary, not a
classification. Over the two finite Markoff families scanned in the supplement
it admits exactly one further value, `m₀ = 34`; one displayed ray fixing `34`
is shown to carry no common integral half-step, but the remaining rays fixing
`34`, and any further values meeting the condition, are not classified here.
The constant `γ = log φ / log λ` that governs both the certified reductions and
the density theorem is the exchange rate between the two arms. Transport to
another branch would require both the square condition and a compatible action
on that branch's two-square sequence, and neither is supplied by the general
branch recurrence.

The work is motivated by the Markoff–Frobenius uniqueness conjecture, and the
final sections make the boundary explicit. **No reduction from a hypothetical
Markoff collision to these orbit families is claimed.** The arithmetic theorems
are complete on their stated domains; the upstream reduction remains open.

## Repository layout

```
.
├── DD-LV-Fibonacci-Pell-Gaps.tex        main file: preamble, abstract, \input list
├── DD-LV-Fibonacci-Pell-Gaps.pdf        compiled manuscript (47 pages), tracked
├── references.bib                       bibliography, 23 entries, all cited
├── Makefile                             strict build; `make help` lists all targets
├── CITATION.cff                         citation metadata (Citation File Format 1.2.0)
├── LICENSE                              GNU LGPL v2.1
│
├── sections/                            one file per section, in \input order
│   ├── 01-introduction.tex              problem, both main theorems, logic diagram
│   ├── 02-arithmetic-system.tex         unit dictionary, window densities, rectangle
│   ├── 03-leading-exponents.tex         even branch: which unit exponents can lead
│   ├── 04-height-bounds.tex             even branch: height separation, Matveev form
│   ├── 05-effective-classification.tex  moving form, reductions, enumeration
│   ├── 06-odd-exponents.tex             odd branch: reordered bootstrap, empty list
│   ├── 07-zero-defect-bridge.tex        seed/unit conjugacy, orbit theorem
│   ├── 08-local-clocks.tex              p-adic ranks, Haar recurrence, limits
│   ├── 09-two-arms.tex                  why φ and λ: trace criterion, dictionary
│   ├── 09-scope-outlook.tex             scope, Markoff boundary, next obstruction
│   └── appendix-certificates.tex        certificate architecture and coverage
│
└── supplement/                          reproduction supplement
    ├── run_all.py                       one-command driver for every certificate
    ├── README.md                        claim-to-source map, arithmetic conventions
    ├── expected-output.txt              recorded output of a successful run
    ├── PROVENANCE.md                    research-tree origins, 25 SHA-256 digests
    ├── certificates/                    21 exact standard-library checkers
    └── magma/                           optional Magma inputs and transcripts
```

## Building the manuscript

Requirements: a TeX Live installation providing `amsart`, `latexmk`, `bibtex`,
and the packages listed in the preamble (`mathtools`, `amssymb`, `booktabs`,
`tikz`, `cleveref`, `aliascnt`, `hyperref`, `microtype`).

```sh
make            # build DD-LV-Fibonacci-Pell-Gaps.pdf
make help       # list every target
make rebuild    # force a full rebuild
make check      # alias for rebuild
make clean      # remove LaTeX intermediates, keep the PDF
make distclean  # remove intermediates and the PDF
```

The build is gated on a warning-free log. It fails on any LaTeX, package or
class warning, on any underfull or overfull box, on a missing character, on a
multiply-defined label, on a BibTeX warning, and on any use of the flat
relations `\leq` or `\geq` in place of the slanted `\leqslant` and `\geqslant`.
Five warning texts are whitelisted by exact match, and are documented in the
`Makefile` where the exclusion is written: the pdfTeX font-expansion notice,
and `amsplain`'s missing-pages warning for the four bibliography entries whose
publisher locator is an article number rather than a page range. A successful
`make` leaves the tree clean: `latexmk -c` removes every intermediate and the
PDF is the only build product.

Note that `make distclean` deletes the tracked PDF; run `make` afterwards to
restore it before committing.

## Reproducing the certificates

Requirements: Python 3.11 or later. No third-party package, no network access
and no nondeterminism is involved.

```sh
python3 -B supplement/run_all.py
```

The driver runs all twenty-one checkers in a fixed order, stops at the first
failure, and prints the failing command with its complete output. A successful
run ends with

```text
PASS: every manuscript certificate completed successfully.
```

The complete reference output is recorded in `supplement/expected-output.txt`,
so a run can be compared line by line:

```sh
python3 -B supplement/run_all.py > /tmp/certs.out 2>&1
diff <(tail -n +4 supplement/expected-output.txt | grep -v '^[[:space:]]*$') \
     <(grep -v '^[[:space:]]*$' /tmp/certs.out)
```

`make certificates` runs the same driver, and `make release` rebuilds the
manuscript and then runs every certificate.

Proof-bearing decisions are made with exact integer, rational, or
outward-rounded fixed-point arithmetic. An element of `Z[√2]` is stored as a
pair of arbitrary-precision integers, so norm, conjugation, content and sign
comparisons are exact. Logarithms and square roots that decide the effective
classification are enclosed by rational endpoints, never compared in binary
floating point. `supplement/README.md` gives the full claim-to-source map.

## What the computation does and does not establish

The manuscript is deliberate about this distinction, and so is the supplement.

**Completeness inputs.** The unbounded estimates come first. Two logarithmic
forms and certified rational-interval reductions bound both exponents, after
which the terminal enumeration inspects a complete finite domain: `2 ⩽ q < 90`
with `3 ∤ q` in the even branch, and `q = 3, 6, …, 93` in the odd branch. The
lower boundary `q = 1` is excluded analytically, not by search. Only these
finite domains, together with the exact interval certificates and the published
classifications of Alekseyev and of Alekseyev–Tengely, carry the completeness
claim.

**Regression evidence only.** The wider runs through `q ⩽ 2000`, the density
counts through `q ⩽ 5000`, and the bounded orbit comparisons are consistency
checks. They are not the source of any completeness statement, and no finite
search is extrapolated beyond its proved range.

**Explanatory, and labelled as such.** The two-arms certificate marks every
line it prints. A `PROVED-CHECK` line verifies an instance of a statement
proved in the manuscript; a `VERIFIED` line records a finite scan and nothing
more. The trace criterion `3m₀ − 2 = square` is a theorem, but the assertion
that only `m₀ = 1`, `2` and `34` satisfy it is a scan over the `3502` Markoff
numbers below `10⁶⁰` and over the `4097` produced by the Christoffel tree to
depth `11`. Whether `3m − 2` is a perfect square for infinitely many Markoff
numbers `m` is not known, and nothing in the manuscript depends on the answer.
The failure of a half-step on the displayed fixed-`34` ray is a separate exact
argument about one explicit linear system, and it classifies neither the other
rays fixing `34` nor any further value meeting the criterion.

**Optional.** Magma is never required. The leading-exponent step rests on the
published Alekseyev–Tengely classification. `supplement/magma/` preserves two
independent Magma V2.29-9 inputs with their complete recorded transcripts,
version information, fixed seeds, integral-point output and Mordell–Weil proof
flags, as corroboration that can be inspected without the proprietary
executable.

## A note on the vendored certificates

The files in `supplement/certificates/` and `supplement/magma/` are byte-for-byte
snapshots of the audited research scripts, and their SHA-256 digests are pinned
in `supplement/PROVENANCE.md`. They are intentionally never edited here, not even
to reformat them or to add a header; when one of them changes it is re-copied
from its source and its digest refreshed, which is why `PROVENANCE.md` records
two vendoring dates. Their docstrings still name their original research-tree
paths (`scripts/…`, `attempts/…`); that is provenance, and it is
what makes the digests verifiable. To confirm that nothing has drifted, compare
each file against its recorded digest:

```sh
cd supplement && sha256sum -c <(grep '^| `' PROVENANCE.md \
  | awk -F'|' '{gsub(/[ `]/,"",$2); gsub(/[ `]/,"",$4); print $4"  "$2}')
```

## Verified status

The tree in this repository was built and checked end to end:

| Check | Result |
| --- | --- |
| `make rebuild` | exit 0, 47 pages |
| Strict warning gates | clean |
| Cross-references | 155 distinct labels referenced, no broken target, no duplicate label |
| Citation keys | 23 of 23 entries cited, every cited key present in `references.bib` |
| `python3 -B supplement/run_all.py` | exit 0, output identical to `expected-output.txt` |
| `supplement/PROVENANCE.md` | 25 of 25 SHA-256 digests match |

## Citation

Citation metadata is provided in `CITATION.cff`, from which GitHub renders a
"Cite this repository" button. In BibTeX:

```bibtex
@unpublished{DutykhVuillon2026FibonacciPellGaps,
  author = {Dutykh, Denys and Vuillon, Laurent},
  title  = {Fibonacci--Pell Nearest Gaps: An All-Exponent Classification
            and Quadratic-Unit Orbit Rigidity},
  year   = {2026},
  note   = {Preprint},
  url    = {https://github.com/dutykh/fibonacci-pell-nearest-gaps}
}
```

## License

The entire repository, manuscript sources and certificates alike, is released
under the GNU Lesser General Public License, version 2.1. The full text is in
`LICENSE`.

## Contact

Denys Dutykh, <denys.dutykh@ku.ac.ae>
Laurent Vuillon, <laurent.vuillon@univ-smb.fr>
