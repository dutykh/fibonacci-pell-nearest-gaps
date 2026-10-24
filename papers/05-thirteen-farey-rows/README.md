<!--
Euclidean minimality of Markov representations on thirteen Farey rows

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Euclidean minimality of Markov representations on thirteen Farey rows

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

*MSC 2020*: Primary 11D25, 68R15; Secondary 11A55, 11E25, 11J06, 11Y50.
*Keywords*: Markov numbers, Christoffel words, sums of two squares, Euclidean
algorithm, Gaussian integers, palindromic closure.

This is the fifth manuscript of the series described in the
[repository README](../../README.md). The directory holds the complete LaTeX
source, the compiled PDF, and the supplementary material carrying the exact
computations that accompany the proofs. Every command below is run from this
directory.

## What is proved

Markov numbers are indexed by the reduced Farey fractions `P/Q` with
`0 < P < Q`, so that `m_{1/2} = 5`, `m_{1/3} = 13` and `m_{2/3} = 29`. Each such
occurrence induces a primitive representation of its label as a sum of two
squares, `m_{P/Q} = r² + s²`, and the subtractive Euclidean algorithm, which
counts single subtractions until the pair reaches `(1, 1)`, reduces this pair
in exactly `P + Q − 2` steps. The article asks whether another primitive
representation of the same integer, Markov or not, can have a shorter
Euclidean algorithm, and answers it on the first thirteen rows. For coprime
`0 < P < Q` with `P ⩽ 13`,

```text
m_{P/Q} = r² + s²,  0 < r < s,  gcd(r, s) = 1    ⟹    ℓ(r, s) ⩾ P + Q − 2,
```

where `ℓ(r, s)` is the subtractive Euclidean cost, with equality for the pair
induced by the occurrence. There is no bound on the denominator `Q` or on the
label. The minimum is weak, but for the 248 occurrences of the finite windows
the occurrence pair is the only representation of minimum cost. Since equal
labels of equal cost come from the same fraction, the minimality implies that
the Markov numbers of these thirteen rows are pairwise distinct, a
consequence that also follows from the results of Lee, Li, Rabideau and
Schiffler on the ordering of Markov numbers, completed by a few explicit
comparisons of labels; what is new is the minimality itself, a statement
about every primitive representation of the label, and the method that
proves it.

The general tool of the proof is a parity-order theorem. At cost `c` the
largest norm is the Fibonacci number `F_{2c+3}`, every other Euclidean path
departs from the Fibonacci path at a number of defects, and once the parities
of the run lengths between the defects are fixed, the normalized norm
`N/F_{2c+3}` is strictly monotone in each run length, with an explicit
direction, for every positive number of defects. Each infinite box of run
lengths is thereby reduced to two exact corner values in `Q(√5)`. An exact
expansion of the Cohn trace places the labels of each row in a narrow
normalized interval that every shorter path misses, complete Gaussian
factorizations settle the 248 labels of the finite windows, and the
finitely many assertions left to a computer are verified in exact arithmetic,
each twice by independent programs.

Minimality on every row would imply the Frobenius–Markov uniqueness
conjecture, and the article conjectures that it holds. Every one of the 4 145
occurrences with `Q ⩽ 100` or `P + Q ⩽ 160` is strictly minimal, although the
margin can be a single unit, as for `m_{4/5} = 985 = 12² + 29² = 16² + 27²`,
with costs 7 and 8. Among the occurrences with `P + Q ⩽ 160`, a margin at
most `k`, for `k = 1, …, 6`, occurs only when `P + Q ⩽ 24k + 3`, and two Farey
ladders converging to `3/5` show that no smaller threshold is possible: for
every `k ⩾ 1`, conjugating the Gaussian factor above 5
in the occurrence pair of `(9k + 1)/(15k + 2)` or `(9k − 1)/(15k − 2)` gives a
primitive representation whose cost exceeds the occurrence cost by exactly
`k`. A second conjecture, on the directives of the Christoffel words of a
given length, also implies uniqueness: for every Markov number `m ⩾ 5` exactly
one `n` prime to `m` with `1 ⩽ n < m/2` satisfies

```text
m/n = [2; B₁, …, B_j, 2],   B_i ∈ { (2), (1, 1) },
```

or equivalently exactly one word `w` beginning with the letter `a` satisfies
`|Pal(θ(w))| = m − 2`, where `Pal` is the iterated palindromic closure and
`θ` the Thue–Morse morphism. It contains the injectivity of Reutenauer and
Vuillon, and it is verified for all 150 Markov numbers between 5 and 10¹².
The square roots of −1 modulo the label carry its whole two-square fibre,
which recasts the first conjecture in terms of the label and its square roots
of −1 alone. These are
results on thirteen rows and evidence beyond them: the Frobenius–Markov
uniqueness conjecture remains open.

## Directory layout

```
papers/05-thirteen-farey-rows/
├── DD_LV_Thirteen_Farey_Rows.tex  main file: preamble, abstract, back matter, \input list
├── DD_LV_Thirteen_Farey_Rows.pdf  compiled manuscript (40 pages), tracked
├── references.bib                 bibliography
├── Makefile                       strict build; `make help` lists all targets
├── README.md                      this file
│
├── sections/                      one file per section, in \input order
├── figures/                       the TikZ sources of the eight figures and their palette
│
└── supplement/                    supplementary material, self-contained
    ├── README.md                  conventions, file formats, and which computation supports what
    ├── run_all.py                 runs the thirteen fixed checks under time and memory limits
    ├── early/                     the first seven rows
    ├── checks/                    finite windows, row labels, junctions, ladders, directives
    ├── tail/                      the covers of the unbounded tails and their saved records
    ├── code/                      the two censuses, the figures, the build gate, data preparation
    ├── data/                      Gaussian fibres, primality checks, censuses, reference record
    ├── CITATION.cff               citation metadata for the material alone
    └── LICENSE                    GNU LGPL v2.1
```

## Building the manuscript

Requirements: a TeX Live installation providing `amsart`, `latexmk` and
`bibtex`, together with the packages named in the preamble (TikZ among them),
and Python 3.10 or later for the gate.

```sh
make            # build DD_LV_Thirteen_Farey_Rows.pdf
make help       # list every target
make rebuild    # force a full rebuild
make figures    # rewrite the TikZ code of the figures from exact data
make clean      # remove LaTeX intermediates, keep the PDF
```

The build is gated: after `latexmk`, `supplement/code/check_build.py` reads the
log, the bibliography log and the sources, and fails on any LaTeX, package or
class warning, on any underfull or overfull box, on a BibTeX warning other than
the missing page ranges of five articles cited by their article numbers, on a
repeated prose space, and on any use of the flat relations in place of the
slanted `\leqslant` and `\geqslant`.

## Reproducing the computations

Requirements: Python 3.10 or later and SymPy, validated with Python 3.14.4 and
SymPy 1.14.0 (`python3 -m pip install sympy==1.14.0`). The checks need no
network access and read nothing outside the supplement directory.

```sh
make checks     # the thirteen fixed checks, each limited in time and memory
```

The thirteen checks, in order, verify the exact source evaluations and
margins of the first seven rows; verify the algebraic identities of the
occurrence pair and of the path envelopes; verify
the complete Gaussian fibres of the 248 labels of the finite windows, first
from the Gaussian sign cube with complete Lucas primality checks and then
independently by Cornacchia's descent from every square root of −1; verify the
algebra of the tail bounds; replay the covers of the unbounded tails, 28
groups, 1 976 parity sectors and 1 023 768 boxes, in exact arithmetic of
`Z[φ]`, and rebuild every tree independently in `Z[√5]` through the spectral
projectors of the path matrix; verify the six junctions between the finite
windows and the tails, with the signs of the named constants; verify, twice,
the six comparisons of constants and the thirty-eight comparisons of labels
that keep every row in its band up to the second junctions; verify the finite
facts on the two ladders converging to `3/5`, and independently recompute the
whole census of margins; and verify the finite facts on directives for every Markov
number up to 10¹², against the census written by the second, independent
program. Each check ends with a record of its counts and digests, which the
runner compares with `supplement/data/verification-reference.json`; any
difference fails. On an idle machine the thirteen checks take about a minute.

`make release` rebuilds the manuscript, runs the thirteen checks, and then
removes the build intermediates. From the top of the repository the same work
is reached as `make release-05-thirteen-farey-rows`.

Two censuses are run on their own, from the supplement directory, and are not
part of the fixed checks. The census of margins needs PARI/GP 2.15 or later
and takes about ten minutes; the census of directives, the second verification
of the finite facts on directives, needs Julia 1.9 or later and takes about ten
seconds on one thread:

```sh
cd supplement
gp -q -D parisizemax=2000000000 code/margins_census.gp
julia --threads=1 code/directives_census.jl
```

They rewrite `data/margins-census.txt` and `data/directives-census.txt`, which
the fixed checks then recompute line by line with other formulas.
`supplement/README.md` gives the arithmetic conventions, the contents of each
program and data file, the optional preparation of the Gaussian data, and the
commands in full.

## What the computation does and does not establish

The reductions of the proof are established in the text: that the sign cube
lists the whole fibre of a label, that two matched corners bound a whole box,
that the root boxes cover every admissible tuple of run lengths exactly once,
arbitrarily long runs included, and that the labels of a row stay in their
band beyond the junction. The computer contributes only the finitely many
assertions that the article collects in its proposition on the finite
verifications (the factorizations and costs of 1 638 explicit pairs, the
closing of the boxes of finitely many trees, thirty-eight comparisons of
integers and six comparisons of constants), together with the finite facts
on the two ladders and on directives; each is verified twice, by programs that
use different formulas. The census of margins is evidence for the two
conjectures on margins, and the facts on directives support the conjecture on
directives up to 10¹²; neither establishes the Frobenius–Markov uniqueness
conjecture.

## Citation

Citation metadata for the repository as a whole is in the `CITATION.cff` at the
top of the repository, from which GitHub renders a "Cite this repository"
button. For this manuscript in particular, in BibTeX:

```bibtex
@unpublished{DutykhVuillon2026Thirteen,
  author = {Dutykh, Denys and Vuillon, Laurent},
  title  = {Euclidean Minimality of {Markov} Representations on Thirteen {Farey} Rows},
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
