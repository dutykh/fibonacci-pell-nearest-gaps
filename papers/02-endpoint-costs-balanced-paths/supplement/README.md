# Exact computations for the norm-89 reflection

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

This package contains the exact computations behind the norm-89 endpoint theorem and its quantitative refinement for balanced inputs. It reconstructs a complete finite transition graph rather than sampling continued-fraction words. The graph covers inputs of arbitrary length at this one Gaussian reflection block. It does not prove minimality against other blocks or settle Markoff--Frobenius uniqueness.

Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon. The package is distributed under the GNU Lesser General Public License, version 2.1 (see `LICENSE`); please cite the article when using it (see `CITATION.cff`). The packages of all the articles of this series are collected at <https://github.com/dutykh/fibonacci-pell-nearest-gaps/>.

## Running the checks

Python 3.10 or later and a Unix-like system with Python's standard `resource` and `signal` modules are required. The computations were last run with Python 3.14.4. No third-party Python package, network access or factorization program is needed, and nothing outside this directory is read.

From this directory run:

```sh
python3 -B verify.py
```

The command also works from any other working directory when given the path to `verify.py`. It bounds itself to 128 mebibytes of address space and to 30 seconds of processor and wall-clock time. The fixed graph bound is 24000 quotient states and the sign-resolution prefix bound is 64 nodes. Any failed identity, incomplete graph, exceeded bound, or mismatch with the recorded values fails the run. Optimized Python execution is refused, because the assertions carry part of the argument.

The command prints one record and compares it with [the recorded values](data/reference-output.json). It also compares every entry of [the complete integer table of states and potentials](data/axis89-potential.json), which has 11420 rows. Timings and platform-dependent memory measurements are excluded from the comparison. Both files are regenerated with `python3 -B verify.py --write-reference`, which performs all the same mathematical checks first.

## The source dictionary and the excess

```sh
python3 -B verify_dictionary.py
```

This command checks the two statements that the reduction graph takes for granted, and it uses nothing from that graph. For every reduced slope with denominator at most seventy it rebuilds the Christoffel word, its Markoff label, the continued-fraction word of the occurrence and its cost, and confirms the letter counts and the balance of the core asserted by the proposition on the genuine source dictionary. It then evaluates the excess directly from its definition, as a difference of two subtractive Euclidean costs, on every word of the source language up to core length eight, and confirms that it is nonnegative on each of those whose norm has valuation one at 89. Finally it reproduces the zero-cost light family that closes the article, with its state, its excess and its letter counts.

Both statements are proved analytically in the article. This command confirms them independently; it does not supply a step of either proof.

## Contents and mathematical conventions

- `code/elementary.py` implements integral matrices, primitive reflection numerators, complete sign resolution, and exact batched row subtraction.
- `code/endpoint.py` reconstructs every reachable state and primitive terminal, checks every potential inequality, and settles the full equality language by a complete finite automaton product and two token intertwiners.
- `code/critical.py` checks the zero-slack strongly connected components, all their rational count coboundaries, and every inequality in the weighted condensation longest-path bound.
- `code/witnesses.py` checks the generic Markov-recurrence polynomial identities and the fixed examples used in the analytic arguments. It does not enumerate a modular recurrence period.
- `code/unit_reducer.py` is the second, independent reducer described below; it performs every forced row subtraction separately and imports nothing else of this package.
- `data/reference-output.json` contains the exact reproducible record: graph and potential fingerprints, all seeds, equality transitions, critical components and frequencies, and the fixed identities.
- `data/axis89-potential.json` contains the full sorted table of 11420 pairs `[state, potential]`, one per line. Every state, successor, primitive terminal and potential inequality is regenerated; the fingerprints are an additional control, not a replacement for checking the integer table.
- `data/independent-reference.json` contains the record of the second computation, which compares itself against that file on every run.

Matrices are row-major integer tuples. The source paired letter `j` represents `P(j)²`, where `P(j)` is the matrix with rows `(j, 1)` and `(1, 0)`. Its source cost is `2j`. An edge emits the actual number of unit componentwise row subtractions; the edge weight is emitted cost minus source cost. Row exchange is quotiented by choosing positive determinant. The terminal vector is the column `(2, 1)`, with source terminal cost one. A terminal is admitted only when the raw reflected vector has content exactly 89.

The completed cost difference includes the seed and terminal charges. Zero-cost paths have the complete source language `111 χ((12111 | 21111)*) 2`. The seven cyclic zero-slack components have frequencies `0`, `1/5`, `1/3` and `1`. Their exact weighted condensation bound is `L = 58`, giving the balanced count estimate

```text
min  |m - θ n|  <  59 (Δ + 1) + 2
 θ
```
 The final two accounts for paired letters consumed before graph entry. For genuine occurrence counts, applied at the occurrence's own Farey slope, this becomes

```text
min  |P - θ Q|  ⩽  118 Δ + 121
 θ
```


The finite-balanced count argument and the infinite Markov-family construction are analytic proofs in the article. The code checks their finite graph inputs and exact algebraic identities. It does not substitute a finite sample for either all-parameter argument.

## The second, independent computation

An independent implementation of the complete endpoint and critical-graph calculation is included:

```sh
python3 -B verify_independent.py
```

This command is also independent of the working directory and bounds itself in the same way. It imports only `code/unit_reducer.py`, which performs every forced row subtraction separately. A work-list relaxation, an iterative strong-component search and a sink-first condensation calculation replace the ordered edge sweeps, two-pass component calculation and forward longest-path calculation of the first command. It compares all 11420 supplied state potentials, every critical state, edge and count coboundary, the full equality language and the constants `L = 58`, `K = 59`, and initial prefix allowance two. It does not repeat the separate analytic identity checks of `code/witnesses.py`.

Its record is compared against [data/independent-reference.json](data/independent-reference.json) on every run, and that file is regenerated with `python3 -B verify_independent.py --write-reference`. [INDEPENDENT_VERIFICATION.md](INDEPENDENT_VERIFICATION.md) sets out where the two computations differ mathematically. Running the second is optional; both were run for the recorded values.
