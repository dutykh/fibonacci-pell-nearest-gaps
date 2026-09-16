<!--
Endpoint costs and balanced near-minimizers under Gaussian reflections

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Endpoint costs and balanced near-minimizers under Gaussian reflections

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

*MSC 2020*: Primary 68R15; Secondary 05C22, 11D25, 11J70, 68Q45.
*Keywords*: Markoff numbers, Gaussian integers, continued fractions, weighted
graphs, balanced words, Sturmian words, return words.

This is the second manuscript of the series described in the
[repository README](../../README.md). The directory holds the complete LaTeX
source, the compiled PDF, and the supplementary material carrying the exact
computations of the proofs. Every command below is run from this directory.

## What is proved

The article proves a quantitative stability theorem for finite graphs carrying
nonnegative integer initial, edge and terminal weights: when each cyclic
component of the zero-weight subgraph has a single letter frequency, the letter
count of every balanced accepted word lies within an explicit multiple of its
completed cost from that finite set of frequencies, the multiple being read off
the condensation graph. Balance is what allows the finite set to replace its
convex hull. The compatible slopes of a finite word are identified with the
projection of a face of the arrangement by which Berstel and Pocchiola
enumerated the factors of Sturmian words, which also measures how thin the
arithmetic input is: of the balanced words of length `n`, only `O(n)` of the
roughly `n³/π²` arise from a Markoff occurrence.

The theorem is then applied to the Gaussian reflection attached to

```text
89  =  8² + 5²
```

acting on primitive sources whose norm has 89-adic valuation exactly one. On an
unrestricted paired-digit language this reflection never decreases subtractive
Euclidean cost, and every word whose cost it preserves is determined. Every
distinct reflected target of a genuine reduced Markoff source therefore has
strictly larger cost, by an amount that grows linearly with the denominator of
the Farey slope outside any neighbourhood of four critical slopes, and an
explicit infinite family of genuine sources realizes that growth.

The proof of the inequality at 89 is in part computer-assisted; the stability
theorem and the infinite family are not. These comparisons concern the complete
norm-89 flip, and they do not establish global Markoff uniqueness.

## Directory layout

```
papers/02-endpoint-costs-balanced-paths/
├── DD-LV-Endpoint-Costs-Balanced-Paths.tex  main file: preamble, abstract, \input list
├── DD-LV-Endpoint-Costs-Balanced-Paths.pdf  compiled manuscript (20 pages), tracked
├── references.bib                           bibliography, 13 entries
├── Makefile                                 strict build; `make help` lists all targets
├── README.md                                this file
│
├── sections/                                one file per section, in \input order
│   ├── introduction.tex                     the question, the theorems, the boundary
│   ├── cost-and-reflections.tex             Euclidean cost and a complete Gaussian flip
│   ├── markoff-sources.tex                  genuine Markoff sources, balanced cores
│   ├── balanced-stability.tex               balanced finite paths near critical frequencies
│   ├── endpoint-inequality.tex              the endpoint inequality at norm 89
│   ├── equality.tex                         the complete equality language
│   ├── markoff-application.tex              quantitative penalties for genuine sources
│   └── conclusion.tex                       further questions
│
├── figures/
│   └── equality-graph.tex                   the zero-slack graph of the equality language
│
├── tools/
│   └── check_build.py                       the strict log, bibliography and style gate
│
└── supplement/                              supplementary material, self-contained
    ├── README.md                            conventions, and which computation supports what
    ├── INDEPENDENT_VERIFICATION.md          how the two computations differ
    ├── verify.py                            the complete computation
    ├── verify_independent.py                the second, independent computation
    ├── verify_dictionary.py                 the source dictionary and the excess
    ├── code/                                the five modules the drivers use
    ├── data/                                the recorded values and the full state table
    ├── CITATION.cff                         citation metadata for the material alone
    └── LICENSE                              GNU LGPL v2.1
```

## Building the manuscript

Requirements: a TeX Live installation providing `amsart`, `latexmk` and
`bibtex`, together with the packages named in the preamble.

```sh
make            # build DD-LV-Endpoint-Costs-Balanced-Paths.pdf
make help       # list every target
make rebuild    # force a full rebuild
make clean      # remove LaTeX intermediates, keep the PDF
```

The build is gated: after `latexmk`, `tools/check_build.py` reads the log, the
bibliography log and the sources, and fails on any LaTeX, package or class
warning, on an undefined or duplicated label, on a BibTeX warning, and on any
use of the flat relations in place of the slanted `\leqslant` and `\geqslant`.

## Reproducing the computations

Requirements: Python 3.10 or later on a Unix-like system. No third-party
package, no network access and no factorization program is involved, and
nothing outside the supplement directory is read.

```sh
make checks     # the three commands below, in order
```

The first command reconstructs the complete finite transition graph, rather
than sampling continued-fraction words, and compares every one of its 11420
states and potentials against the recorded table. The second repeats that
calculation through an independent implementation, which shares no code with
the first beyond its own reducer, and reaches the same constants. The third
checks the two statements the graph takes for granted, rebuilding the
Christoffel word, the Markoff label, the continued-fraction word and the cost
of every reduced slope with denominator at most seventy, and evaluating the
excess directly from its definition on every word of the source language up to
core length eight. Each run bounds itself in memory and in time, and any failed
identity, incomplete graph, exceeded bound or mismatch with the recorded values
fails it.

`make release` rebuilds the manuscript, runs all three, and then removes the
build intermediates. From the top of the repository the same work is reached as
`make release-02-endpoint-costs-balanced-paths`.

`supplement/README.md` gives the arithmetic conventions, the contents of each
module and data file, and the map from each statement of the article to the
computation that supports it. The supplement carries its own `LICENSE` and
`CITATION.cff` so that it can be submitted to a journal on its own.

## What the computation does and does not establish

The finite-balance count argument and the infinite Markoff family are analytic
proofs in the article. What the computation supplies is their finite graph
input and their exact algebraic identities, over a graph that covers inputs of
arbitrary length at this one Gaussian reflection block. It does not substitute
a finite sample for either argument, it does not prove minimality against other
blocks, and it does not settle Markoff–Frobenius uniqueness.

## Citation

Citation metadata for the repository as a whole is in the `CITATION.cff` at the
top of the repository, from which GitHub renders a "Cite this repository"
button. For this manuscript in particular, in BibTeX:

```bibtex
@unpublished{DutykhVuillon2026EndpointCosts,
  author = {Dutykh, Denys and Vuillon, Laurent},
  title  = {Endpoint Costs and Balanced Near-Minimizers under Gaussian
            Reflections},
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
