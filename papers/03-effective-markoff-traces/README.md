<!--
Effective trace bounds and rigidity for simple curves on the modular torus

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Effective trace bounds and rigidity for simple curves on the modular torus

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

*MSC 2020*: Primary 11D25, 30F60; Secondary 11J06, 68R15.
*Keywords*: Markov equation, simple closed geodesics, Christoffel words, Vieta
involutions, effective Diophantine methods.

This is the third manuscript of the series described in the
[repository README](../../README.md). The directory holds the complete LaTeX
source, the compiled PDF, and the supplementary material carrying the exact
computations of the proofs. Every command below is run from this directory.

## What is proved

The article studies equal traces of primitive simple curves on the modular
once-punctured torus. One curve `A` is taken as the first element of a simple
basis `(A, B)`, and the second curve, of homology `(h, N)`, is written as a
Christoffel word with `h` occurrences of `A` separated by gaps `B^k` and
`B^(k+1)`. The occurrence count `h` is fixed; the gap quotient and the integral
Markov character of the basis are both allowed to vary.

For every `h ⩾ 2` the article gives a complete decision procedure for trace
equality in this family: every primitive mechanical equality lies in an
explicit domain of `O(h⁵)` exact comparisons on `O(h³)`-bit integers, and the
isometry status of every equality is decided by twelve fixed integral actions.
The two parameter bounds behind it are analytic and hold over all positive
integral Markov characters: with `p = tr B / 3` after centring the ray,

```text
p < 4^(h+1)            at positive product index,
p < (100 h²)^(h+1)     at minimum index.
```

The first comes from positive quotient mutation, the second from mechanical
balance and spectral coordinates along the Euclidean algorithm.

Two rigidity statements follow. Every primitive occurrence family with
`2 ⩽ h ⩽ 6` contains only isometric equalities, and every equal-length pair of
simple closed geodesics with intersection number at most 24 is isometric. No
upper bound on the common trace is assumed.

These results leave the occurrence count unbounded and do not resolve the
Markov–Frobenius uniqueness conjecture. A nonisometric equality, if one exists,
has intersection number at least 25.

## Directory layout

```
papers/03-effective-markoff-traces/
├── DD-LV-Effective-Markov-Traces.tex  main file: preamble, abstract, \input list
├── DD-LV-Effective-Markov-Traces.pdf  compiled manuscript (33 pages), tracked
├── references.bib                     bibliography, 14 entries
├── Makefile                           strict build; `make help` lists all targets
├── README.md                          this file
│
├── sections/                          one file per section, in \input order
│   ├── 01-introduction.tex            the question, the main conclusions
│   ├── 00-literature.tex              relation to earlier work, input by the introduction
│   ├── 02-geometry-and-rays.tex       simple curves, mechanical words, Markov rays
│   ├── 03-gap-reduction.tex           positive trace paths and a uniform gap bound
│   ├── 04-positive-vieta.tex          positive product indices, a uniform Vieta estimate
│   ├── 05-spectral-minimum.tex        mechanical minimum patterns and spectral descent
│   ├── 06-finite-decision.tex         the complete finite decision at each count
│   ├── 07-rigidity.tex                rigidity at small occurrence count and intersection
│   ├── 08-outlook.tex                 what remains after effective finiteness
│   └── 09-finite-domains.tex          appendix: finite domains, intersection nineteen
│
├── figures/
│   └── vieta-comparison.tex           the comparison diagram, TikZ
│
├── tools/
│   └── check_build.py                 the strict log, bibliography and style gate
│
└── supplement/                        supplementary material, self-contained
    ├── README.md                      conventions and reproduction
    ├── run_all.py                     every check, under time and memory limits
    ├── requirements.txt               SymPy, pinned
    ├── code/trace_decision.py         the general decision procedure for any h
    ├── checks/                        the six exact programs, with their README
    ├── data/                          equality records, polynomial data, residue cover
    ├── CITATION.cff                   citation metadata for the material alone
    └── LICENSE                        GNU LGPL v2.1
```

## Building the manuscript

Requirements: a TeX Live installation providing `amsart`, `latexmk` and
`bibtex`, together with the packages named in the preamble.

```sh
make            # build DD-LV-Effective-Markov-Traces.pdf
make help       # list every target
make rebuild    # force a full rebuild
make clean      # remove LaTeX intermediates, keep the PDF
```

The build is gated: after `latexmk`, `tools/check_build.py` reads the log, the
bibliography log and the sources, and fails on any LaTeX, package or class
warning, on an undefined or duplicated label, on a BibTeX warning other than
three recorded missing page ranges, and on any use of the flat relations in
place of the slanted `\leqslant` and `\geqslant`.

## Reproducing the computations

Requirements: Python 3.10 or later on a Unix-like system, and SymPy (the
version used is pinned in `supplement/requirements.txt`) for the three symbolic
checks. The occurrence and residue programs use only the standard library. No
network access is needed and nothing outside the supplement directory is read.

```sh
make checks     # runs supplement/run_all.py
```

The runner executes fourteen programs in turn, each under a 30-second limit
and, on Unix, a 128 MiB address-space limit. For each `h = 2, …, 6` a primary
program built on literal integer matrices and an independent scalar program
built on the universal character algebra both reconstruct the whole proved
finite domain and compare every equality record, with its isometry witness,
against `supplement/data/`. Four further programs check the explicit `h = 2`
endpoint controls, the polynomial arguments at intersection numbers 19 and 23,
and the residue cover through intersection 24. The whole suite runs in a few
seconds.

`make release` rebuilds the manuscript, runs every program, and then removes
the build intermediates. From the top of the repository the same work is
reached as `make release-03-effective-markoff-traces`.

`supplement/README.md` and `supplement/checks/README.md` give the data
conventions and the map from each statement of the article to the program that
supports it. The supplement carries its own `LICENSE` and `CITATION.cff` so that
it can be submitted to a journal on its own.

## What the computation does and does not establish

The parameter bounds and the reduction to a finite domain are analytic proofs
in the article. The computation executes that finite domain for `h ⩽ 6`, and
checks the exact polynomial identities and sign conditions used at
intersections 19 and 23. It does not bound the occurrence count, and the two
residue orbits left uncovered at intersection 25 are not counterexamples.

## Citation

Citation metadata for the repository as a whole is in the `CITATION.cff` at the
top of the repository, from which GitHub renders a "Cite this repository"
button. For this manuscript in particular, in BibTeX:

```bibtex
@unpublished{DutykhVuillon2026EffectiveTraces,
  author = {Dutykh, Denys and Vuillon, Laurent},
  title  = {Effective Trace Bounds and Rigidity for Simple Curves on the
            Modular Torus},
  year   = {2026},
  note   = {Preprint},
  url    = {https://github.com/dutykh/fibonacci-pell-nearest-gaps}
}
```

## License

The manuscript sources and the supplementary material alike are released under
the GNU Lesser General Public License, version 2.1. The full text is in
[`LICENSE`](../../LICENSE) at the top of the repository, and a copy travels
with the supplement.

## Contact

Denys Dutykh, &lt;denys.dutykh@ku.ac.ae&gt;
Laurent Vuillon, &lt;laurent.vuillon@univ-smb.fr&gt;
