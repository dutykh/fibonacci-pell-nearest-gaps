<!--
Around the Markoff uniqueness conjecture: manuscripts and reproduction
material.

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Around the Markoff uniqueness conjecture

**Manuscripts and reproduction material**

Denys Dutykh<sup>1</sup> and Laurent Vuillon<sup>2</sup>

<sup>1</sup> Mathematics Department, Khalifa University of Science and
Technology, PO Box 127788, Abu Dhabi, United Arab Emirates
&lt;denys.dutykh@ku.ac.ae&gt;
<sup>2</sup> Univ. Savoie Mont Blanc, CNRS, LAMA, 73000 Chambéry, France
&lt;laurent.vuillon@univ-smb.fr&gt;

## The question behind the series

A Markoff triple is a positive integer solution of

```text
x² + y² + z² = 3xyz
```

and every such triple is obtained from `(1, 1, 1)` by repeatedly replacing one
coordinate through the Vieta involution, so that the solutions form a tree. The
largest coordinate of a triple is its Markoff number. Frobenius asked in 1913
whether a Markoff number determines its triple: whether, for each `m`, there is
at most one triple whose maximum is `m`. The question is still open, and it is
the question this series circles.

Each manuscript here attacks one face of it and proves what can be proved
completely, on a domain stated without hidden hypotheses. None of them claims
the conjecture, and each says in its own words where its boundary lies. What
they share is a method: analytic estimates first, cutting an unbounded problem
down to a finite domain, and only then an exact computation over that domain,
carried out in integer or rational arithmetic and reproducible from the
supplement that travels with the manuscript.

This repository is the public home of those manuscripts and of everything
needed to re-run their computations. Every manuscript directory is
self-contained: its own sources, its own bibliography, its own build, its own
supplement. Nothing below `papers/` reaches into anything else.

## The manuscripts

| No | Manuscript | Directory | State |
| --- | --- | --- | --- |
| 1 | Fibonacci–Pell nearest gaps: an all-exponent classification and quadratic-unit orbit rigidity | [`papers/01-fibonacci-pell-nearest-gaps`](papers/01-fibonacci-pell-nearest-gaps) | complete, 47 pages, supplement of 21 checkers |
| 2 | Endpoint costs and balanced near-minimizers under Gaussian reflections | `papers/02-endpoint-costs-balanced-paths` | being prepared for this repository |

The first classifies, for every `q ⩾ 1`, every sign and every positive exponent
`n`, the strict factor-two gaps between a power of the Pell unit `λ = 1 + √2`
and the square of a Fibonacci-built element of `Z[√2]` whose norm is minus the
square of its content. Exactly two such gaps exist, none with an odd exponent,
and a companion theorem lists every Pell-orbit hit of the two signed Fibonacci
cores. Its own README states the results and their boundary in full.

The second proves a stability theorem for Euclidean costs along a finite graph
of continued-fraction states, and applies it to the Gaussian reflection
attached to `89 = 8² + 5²`, showing that every distinct reflected target of a
genuine reduced Markoff source has strictly larger cost.

## Repository layout

```
.
├── README.md      this file: the series, its manuscripts, how to build them
├── LICENSE        GNU LGPL v2.1, covering everything in the repository
├── CITATION.cff   citation metadata for the repository as a whole
├── Makefile       forwards every target to each manuscript below papers/
│
└── papers/
    └── 01-fibonacci-pell-nearest-gaps/
        ├── README.md                     the results, the boundary, the build
        ├── Makefile                      strict build and reproduction targets
        ├── DD-LV-Fibonacci-Pell-Gaps.tex main file
        ├── DD-LV-Fibonacci-Pell-Gaps.pdf compiled manuscript, tracked
        ├── references.bib                bibliography
        ├── sections/                     one file per section
        └── supplement/                   reproduction material
```

## Building

A manuscript is built from its own directory with `make`, exactly as if it were
alone in a repository of its own:

```sh
cd papers/01-fibonacci-pell-nearest-gaps
make            # build the PDF under that manuscript's strict warning gates
make help       # every target that manuscript offers
make checks     # run the reproduction programs of its supplement
```

The `Makefile` at the top of the repository does no work itself. It discovers
the manuscript directories below `papers/` and forwards to each of them, so the
whole series is built and checked with one command:

```sh
make            # build every manuscript PDF
make check      # force a full rebuild of every manuscript
make checks     # run every supplement
make release    # rebuild everything, then run every supplement
make clean      # remove build intermediates everywhere
make list       # the manuscript directories that were discovered
make help       # the above, with the discovered directories listed
```

One manuscript at a time is reached by its directory name, as `make
01-fibonacci-pell-nearest-gaps` to build it, or `check-…`, `checks-…`,
`release-…`, `clean-…` for the other verbs.

Requirements are stated by each manuscript. In general they amount to a TeX
Live installation with `latexmk` and `bibtex`, and Python 3.11 or later with no
third-party package, no network access and no nondeterminism.

## What a manuscript directory contains

A directory under `papers/` is numbered in the order the manuscripts joined the
series, followed by a short name, as `02-endpoint-costs-balanced-paths`. It
holds the LaTeX sources of one manuscript and nothing else: the main file, the
sections, the bibliography, the compiled PDF, a README stating the results and
their boundary, a supplement carrying the reproduction material, and a
`Makefile`.

That `Makefile` is the only thing the series asks of it. Five target names are
common to every manuscript, which is what lets the top-level `Makefile` treat
them alike: `all` builds the PDF, `check` forces a full rebuild under whatever
strict gates that manuscript imposes, `checks` runs the reproduction programs
of its supplement, `release` does both, and `clean` removes the build
intermediates. A manuscript may offer further targets of its own, and those are
reached directly with `make -C papers/<directory> <target>`.

Working notes, plans and internal reviews stay in the research tree and do not
enter this repository. What lands here is the manuscript, the material a reader
needs to check it, and the record of what that material establishes.

## Licence

Everything in this repository, manuscript sources and supplements alike, is
released under the GNU Lesser General Public License, version 2.1. The full
text is in [`LICENSE`](LICENSE).

## Citation

Citation metadata for the repository is in `CITATION.cff`, from which GitHub
renders a "Cite this repository" button. Each manuscript is cited on its own,
and its README gives the BibTeX entry for it.

## Contact

Denys Dutykh, &lt;denys.dutykh@ku.ac.ae&gt;
Laurent Vuillon, &lt;laurent.vuillon@univ-smb.fr&gt;
