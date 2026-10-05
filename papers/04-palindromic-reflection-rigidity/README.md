<!--
Rational reflections and palindromic words in the Cohn monoid

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Rational reflections and palindromic words in the Cohn monoid

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

*MSC 2020*: Primary 68R15; Secondary 11D25, 11A55.
*Keywords*: palindromic closure, standard morphisms, Thue–Morse morphism, Cohn
matrices, Markov equation, rational reflections, Frobenius root.

This is the fourth manuscript of the series described in the
[repository README](../../README.md). The directory holds the complete LaTeX
source, the compiled PDF, and the supplementary material carrying the exact
computations that accompany the proofs. Every command below is run from this
directory.

## What is proved

A palindrome `p` defines a rational reflection `R_p` of the plane, which fixes
the terminal column `e = (2, 1)` and is orthogonal for the Gram form of `p`, for
the two Cohn matrices `A = ((2,1),(1,1))` and `B = ((5,2),(2,1))`. The article
classifies every finite word `x` whose terminal column `μ(x)e` the reflection
carries to another word column: this happens exactly when `x` is a concatenation
of the blocks `ab·p` and `ba·p`, and the reflection then exchanges the two
blocks. The projective relation lifts to exact equality of columns, and for a
standard centre a second proof, by projective cylinders, measures how far a
reflected column can move.

For palindromes the word theorem has an arithmetic consequence. Two palindromes
with the same Gram label and a root midpoint at the ratio of the centre are
exactly two palindromes exchanged around it, with no hypothesis on a common
prefix, and every nontrivial pair has nonprimitive completed letter counts, so
that no two distinct central occurrences collide at such a midpoint. An explicit
quadruple of palindromes of label 12 986 074 130 shows that the midpoint
hypothesis cannot be dropped, and exchanges around nested centres compose into
cubes of palindromes with one label.

In the Thue–Morse coding of Reutenauer and Vuillon the Gram label of a
palindrome is the length of an explicit Christoffel word, and the Gram root is
its number of letters `b`, which also fixes where its standard factorization
cuts. The length of a directing word is half the sum of the partial quotients of
the label divided by the root; for Markov words the equality of these lengths
under a coincidence of labels is equivalent to Frobenius–Markov uniqueness.
Away from the midpoint the article derives the exact two-sided reflection
equation and a fourth-power bound on the common-ancestor label when the parent
imbalance is bounded.

An exact enumeration of all 16 743 538 palindromes with label at most 10²⁶
finds no coincidence of labels that the exchanges do not produce, and every
class it finds is a single exchange or a square. This supports the conjecture
on twins stated in the closing section, whose proof would imply the uniqueness
conjecture. These are structural restrictions: no result excludes every
displaced central pair, and the uniqueness conjecture is not resolved.

## Directory layout

```
papers/04-palindromic-reflection-rigidity/
├── DD-LV-Palindromic-Reflection-Rigidity.tex  main file: preamble, abstract, \input list
├── DD-LV-Palindromic-Reflection-Rigidity.pdf  compiled manuscript (38 pages), tracked
├── references.bib                             bibliography
├── Makefile                                   strict build; `make help` lists all targets
├── README.md                                  this file
│
├── sections/                                  one file per section, in \input order
├── figures/                                   the three figure sources, written by supplement/code/make_figures.py
├── tools/
│   └── check_build.py                         the strict log, bibliography and style gate
│
└── supplement/                                supplementary material, self-contained
    ├── README.md                              conventions, and which computation supports what
    ├── run_all.py                             runs the six fixed checks with time limits
    ├── checks/                                the six fixed checks
    ├── code/                                  reference words and Gram matrices, the census, the figures
    ├── data/                                  the worked examples and the recorded output
    ├── requirements.txt                       SymPy 1.14.0
    ├── CITATION.cff                           citation metadata for the material alone
    └── LICENSE                                GNU LGPL v2.1
```

## Building the manuscript

Requirements: a TeX Live installation providing `amsart`, `latexmk` and
`bibtex`, together with the packages named in the preamble, and Python 3.10 or
later for the gate.

```sh
make            # build DD-LV-Palindromic-Reflection-Rigidity.pdf
make help       # list every target
make rebuild    # force a full rebuild
make clean      # remove LaTeX intermediates, keep the PDF
```

The build is gated: after `latexmk`, `tools/check_build.py` reads the log, the
bibliography log and the sources, and fails on any LaTeX, package or class
warning, on a BibTeX warning, on a repeated prose space, and on any use of the
flat relations in place of the slanted `\leqslant` and `\geqslant`.

## Reproducing the computations

Requirements: Python 3.10 or later and SymPy 1.14.0 (`pip install -r
supplement/requirements.txt`). No network access is involved, and nothing
outside the supplement directory is read.

```sh
make checks     # the six fixed checks, each limited in time and memory
```

The six checks verify the symbolic identities and rational constants of the
word decoder and of the decoder at every palindromic centre, reconstruct the
worked examples of the article independently from the literal words, verify the
identities and constants of the moving reflection and of the fourth-power bound
in both parities of the label, verify the Thue–Morse coding and the block
exchange at every palindromic centre of length at most five together with the
four palindromes of label 12 986 074 130, and verify the length of a directing
word, the cubes of nested exchanges up to depth three and the first square. Each
prints a PASS line and a hash that is compared with the recorded output.

`make release` rebuilds the manuscript, runs the six checks, and then removes the
build intermediates. From the top of the repository the same work is reached as
`make release-04-palindromic-reflection-rigidity`.

Two further tools are run on their own, preferably under a memory limit. The
enumeration of twins up to label 10²⁶ visits 16 743 538 palindromes and takes
about 8 seconds and 2 GB; the search for decoder solutions at every palindromic
centre of length at most seventeen takes about 2 seconds:

```sh
python3 -B supplement/code/twin_census.py --bound 1e26
python3 -B supplement/code/centre_decoder.py --max-length 17
python3 -B supplement/code/make_figures.py      # rewrites the three figure sources
```

`supplement/README.md` gives the arithmetic conventions, the contents of each
module and data file, and the commands in full.

## What the computation does and does not establish

Every theorem of the article is proved in the text, by hand. The fixed checks
confirm the displayed identities and rational constants of the proofs and the
entries of the worked examples in exact arithmetic; no enumeration of
directives or of finite words enters any proof. The enumeration of twins is
evidence for the conjecture on twins within its stated range and is used in no
proof. Neither the checks nor the enumeration establish Frobenius–Markov
uniqueness.

## Citation

Citation metadata for the repository as a whole is in the `CITATION.cff` at the
top of the repository, from which GitHub renders a "Cite this repository"
button. For this manuscript in particular, in BibTeX:

```bibtex
@unpublished{DutykhVuillon2026Reflections,
  author = {Dutykh, Denys and Vuillon, Laurent},
  title  = {Rational Reflections and Palindromic Words in the Cohn Monoid},
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
