# Exact trace comparisons and modular-torus rigidity

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

This package accompanies the article *Effective trace bounds and rigidity for
simple curves on the modular torus* by Denys Dutykh and Laurent Vuillon.
Its purpose is to reproduce the complete finite applications of the analytic
bounds proved in the article. The finite checks do not prove those bounds
or settle Markov uniqueness for unbounded occurrence count.

Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon. The package is
distributed under the GNU Lesser General Public License, version 2.1 (see
`LICENSE`); please cite the article when using it (see `CITATION.cff`). The
packages of all the articles of this series are collected at
<https://github.com/dutykh/fibonacci-pell-nearest-gaps/>.

## Contents and conventions

- `code/trace_decision.py`: the general occurrence-preserving matrix algorithm,
  with the domain bounds, mechanical words, ray traces and all twelve isometries.
- `checks/`: the small-occurrence calculations, an independent scalar
  implementation, the intersection-19 and intersection-23 polynomial proofs,
  and the residue cover.
- `data/`: all equality rows and exact polynomial coefficients used in the paper.
- `run_all.py`: one command reproducing every headline check.
- `LICENSE`, `CITATION.cff`: licence and citation metadata for this package.

Traces in matrix records are ordinary integral matrix traces. A normalized
Markov label is the trace divided by three. Homology uses the root generators
E and F displayed in the paper; vectors are primitive, and the absolute
determinant is the intersection number. Matrix tuples are row-major. All
calculations are deterministic and exact; no random seeds, floating-point
root approximations or numerical tolerances are used.

## Reproduction

Python 3.10 or newer is required. The occurrence and residue programs use
only the standard library; the fixed symbolic controls additionally require
SymPy. The current checks use Python 3.14.4 and SymPy 1.14.0;
the manuscript uses pdfLaTeX and BibTeX through
`latexmk`. From this directory, run:

```sh
python3 -B run_all.py
```

Each subprocess has a 30-second wall-time limit. On Unix the runner also
enforces a 128 MiB address-space limit and a 30-second CPU limit. To check one occurrence count:

```sh
python3 -B checks/classify_small_occurrences.py 6
python3 -B checks/verify_small_occurrences.py 6
```

The first program uses literal integer matrices; the second imports none of
that implementation, generates a scalar Vieta tree and evaluates the universal
four-dimensional character algebra. Both reconstruct the whole proved finite
domain and compare the complete equality records with `data/`. Every equality
retains its tree occurrence, orientation, full character, original gap type,
homology vectors and actual isometry witness. Order-independent fingerprints
summarize all comparisons, including unequal ones. The data are outputs to be
verified, not an oracle controlling enumeration. Add `--write` to the first
program only when deliberately regenerating the reference output.

The general program is separate from the fixed published application:

```sh
python3 -B code/trace_decision.py 6
python3 -B code/trace_decision.py 6 --run --max-comparisons 1000000
```

The first command reports the proved bounds without executing the domain.
Execution requires `--run` and a sufficient explicit comparison budget. The
program streams equality rows, followed by a completion record only after
exhausting the domain. A timeout or interruption has no completion record and
must not be interpreted as a classification. Raising the occurrence count is
a new computation, not part of the published small-occurrence result.

The diagram of the article is drawn from its TikZ source,
`../figures/vieta-comparison.tex`, during the manuscript build. It is a
schematic of the proof, with no numerical source data. All supplement paths
are resolved relative to the scripts, so the check commands can also be called
from a different working directory.
