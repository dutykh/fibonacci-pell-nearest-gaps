# Exact computations for thirteen Farey rows

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry, France)

This package accompanies the proof that, on the first thirteen Farey rows, the representation of a Markov label as a sum of two squares induced by its occurrence has minimal subtractive Euclidean cost among all primitive representations of the label. It carries the finitely many assertions that the article collects in its proposition on the finite verifications: the complete Gaussian fibres of the 248 labels of the finite windows, the closed covers of the unbounded tails, and the source comparisons (six constant comparisons and 38 integer comparisons between the junctions and the second junctions, which with the bound on the labels of a row keep every source in its band). It also carries the finite facts on the two Farey ladders below 3/5 used in the last section, the census of margins behind the figures and the conjecture of that section, and the finite facts on the directives of the Christoffel words of a given length: for each of the 150 Markov numbers m with 5 ⩽ m ⩽ 10¹², the root is the only n < m/2 prime to m for which m/n = [2; B₁, …, Bⱼ, 2] with blocks (2) and (1, 1). Each assertion is verified twice, by two independent programs that use different formulas. The tail covers include every eligible gap tuple, arbitrarily large gaps included; they are not finite-cost samples. The mathematical reductions and the covering arguments are proved in the article. Global Markov uniqueness remains outside the result.

## Code and data

- `run_all.py`: fixed runner for the thirteen checks below; it compares the final record of each check with `data/verification-reference.json`.
- `early/check_early_rows.py`: exact source evaluations and arithmetic margins for the first seven rows.
- `checks/check_foundation_identities.py`: SymPy checks of the algebraic identities used for the occurrence pair and the path envelopes.
- `checks/verify_gaussian.py`: first verification of the finite windows, standard library only: complete Lucas primality checks, Gaussian splits, the complete fibre of every label from the Gaussian sign cube, the Euclidean costs, and the occurrence pair recomputed from the central Gram factor; each label is computed twice, by Stern–Brocot mediants and from the floor increments of the Christoffel word.
- `checks/independent_gaussian.py`: second, independent verification of the finite windows: BPSW primality (SymPy) instead of Lucas checks, and every representation by the Cornacchia descent from each square root of −1 instead of the sign cube.
- `data/gaussian/`: finite-fibre records for the reduced pairs with `P = 8, …, 13` and `Q ⩽ 28, 34, 38, 75, 100, 145` respectively; always `P < Q`.
- `data/lucas_primes.json` and `data/gaussian_splits.json`: explicit primality checks and sums of two squares. `data/factor_hints.json` is used only for optional data preparation.
- `data/source-junctions.json`: the six exact source values and positive margins at the junctions between the finite windows and the unbounded tails.
- `data/margins-census.txt`: one line per occurrence, `P Q margin factorization`, for the 4,145 reduced `P/Q` with `Q ⩽ 100` or `P + Q ⩽ 160`; the margin is the cost of the cheapest other primitive representation minus `P + Q − 2`, or `-` when the label has a single representation, and the factorization is written `p^e*q^f`.
- `data/verification-reference.json`: the final record of each of the thirteen checks, with complete counts and result digests.
- `checks/check_tail_bounds.py`: exact symbolic checks of the parity and normalization identities of the tails.
- `tail/`: first verification of the covers, with exact arithmetic in `Z[φ]`, exhaustive parity-box verification and all 28 saved group manifests. See [its README](tail/README.md) for single-group commands, conventions and the counts of the trees.
- `checks/independent_tail_cover.py`: second, independent verification of the covers, written from the text of the article, with arithmetic in `Z[√5]` and the spectral projectors of the path matrix; it rebuilds every tree and compares every count with `tail/data/catalog.json`.
- `checks/check_assembly.py`: exact source-window and cost-threshold comparisons joining the finite and the unbounded parts.
- `checks/check_row_bounds.py`: first verification of the source comparisons beyond the six junctions, in exact arithmetic of Q(√5): the eigenbasis identities behind the proposition on the labels of a row, controls of its trace formula and of its two-sided bound, the six constant comparisons `C_P (1 + φ^(−4 n_P − 6))^P < U_P` (with `n_P − 1` shown not to suffice), `C_7 (1 + φ^(−10))^7 < 1`, the 38 comparisons `m_{P/Q} < U_P F_{2(P+Q−2)+1}` for the reduced `Q` with `Q_P ⩽ Q < P (n_P + 1)`, and the six labels of row 7 below their Fibonacci numbers.
- `checks/independent_row_bounds.py`: second, independent verification of the same comparisons, with labels as continuants of the continued fraction `[2, σ(z), 2]`, Fibonacci numbers by doubling, and the constants bounded by rational enclosures of √5 to sixty digits.
- `checks/check_margin_family.py`: first verification of the finite facts on the two ladders `(1+3j)/(2+5j)` and `(2+3j)/(3+5j)`: the word identity of their central words, the transfer matrices of the two parities, the four six-rung matrices (nonnegative, determinant one, 50 letters in `L_a`, `L_b` against 48 for the cube of the transfer), and the base cases.
- `checks/independent_margins.py`: second, independent verification: the whole census of margins recomputed from the recorded factorizations with other formulas (continuants, BPSW, square roots of −1, Chinese remainders and Cornacchia's descent), the ladders by Gaussian division up to the fortieth rung, and the thresholds of the margin for `P + Q ⩽ 160`.
- `code/margins_census.gp`: the census of margins for every reduced `P/Q` with `Q ⩽ 100` or `P + Q ⩽ 160` (PARI/GP; about ten minutes; complete factorization, whole fibre from the Gaussian sign cube, costs, margins); it writes `data/margins-census.txt`.
- `code/source_words.py`, `code/export_gaussian_data.py`, `code/prepare_gaussian.py`: optional literal-word construction and bounded Gaussian-data preparation.
- `checks/check_directives.py`: first verification of the finite facts on directives: the Markov numbers up to 10¹² from the Markov triples of the Farey tree, their roots from the residues of −m_L/m_R modulo m, and for each label every n < m/2 prime to m of the block form, found by a search over the blocks from the left with the last blocks solved for exactly; the examples 31 and 1130 by a complete scan, and the label of 9/10. It compares everything with `data/directives-census.txt`, and controls three statements proved in the article: the Gaussian recursion along the Farey tree for every reduced P/Q with Q ⩽ 100, the correspondence between the square roots of −1 and the primitive representations for every label up to 20000, and the directives of the Christoffel words of length at most 300.
- `code/directives_census.jl`: second, independent verification of the same facts (Julia): it enumerates all the 528,591,759 continued fractions [2; B₁, …, Bⱼ, 2] with numerator at most 10¹² and takes the labels and roots from the Cohn matrices; it writes `data/directives-census.txt`.
- `data/directives-census.txt`: one line per Markov number m with 5 ⩽ m ⩽ 10¹², `m root n`, where `n` lists the denominators of the block form found for the numerator m; the header records the bound, the number of continued fractions enumerated, the integers 31 and 1130 and the label of 9/10.
- `code/make_figures.py`: writes the TikZ code of the figures of the article from exact data; the two figures on margins read `data/margins-census.txt`, the figure of the Farey tree computes each Gaussian integer by the recursion of the article and checks it against the occurrence pair, and the figure of the label 985 computes the cost of every pair (n, 985).
- `code/check_build.py`: complete LaTeX log, bibliography and source-style gate.

All arithmetic is exact; there are no physical units, random seeds or floating-point tolerances. Gaussian records use positive primitive ordered pairs `(x, y)` with `x < y`, norm `x² + y²` and the subtractive Euclidean cost of the article. In the tail packages, inverse-Euclidean path maps act chronologically on `(1, 2)`; `c` is the gap sum plus the number of gaps, and the normalized norm is `N / F(2c + 3)`. These path maps are distinct from the Cohn word matrices.

## Reproduce the checks

Dependencies: Python 3.10 or newer and SymPy. The validated environment uses Python 3.14.4 and SymPy 1.14.0. Install the symbolic dependency in an existing Python environment with:

```sh
python3 -m pip install sympy==1.14.0
```

From the manuscript directory:

```sh
python3 -B supplement/run_all.py
```

The runner and its children resolve paths from their own locations, so an absolute path to the runner works from any working directory. It executes the thirteen checks in order, prints the final record of each, compares it with the reference record and stops at the first failure or difference. Ordinary verification writes no data. Assertions must remain enabled: optimized Python (`-O`) is rejected.

Each child has a 600-second wall limit and, where the Python `resource` module is available, a 1 GiB address-space limit and a 300-second processor-time limit. These limits are generous on purpose: on an idle machine the thirteen checks take about a minute in total, and a busy machine must not turn a correct verification into a failure. Missing files, nonzero exits, empty output, a result different from the reference record or an exceeded limit cause explicit failure. The domain of the tail replay can be inspected without replaying it:

```sh
python3 -B supplement/tail/verify_all.py
```

To build the manuscript and check the complete build log, use `make rebuild` from the manuscript directory. `make release` also runs the checks and removes auxiliary files while retaining the PDF; `make figures` regenerates the TikZ code of the figures. The build requires LaTeX, BibTeX and `latexmk`, separately from the Python dependencies.

## Optional Gaussian-data preparation

The supplied data suffice for verification. In particular, the first verifier performs no integer factorization and uses no probable-prime test: it checks the saved complete Lucas checks and every Gaussian representation directly.

To regenerate the prescribed Gaussian records on Unix, run:

```sh
python3 -B supplement/code/prepare_gaussian.py --all
python3 -B supplement/checks/verify_gaussian.py
```

Preparation uses SymPy for factor discovery, overwrites the prescribed Gaussian records and updates the shared prime and split records. The optional `--resume` flag skips existing fibre records, so the verification step remains necessary. No enlargement of the finite windows or regeneration of tail manifests is part of the release command.

## The census of margins

The census needs PARI/GP (version 2.15 or later) and is not part of the fixed checks. From the supplement directory:

```sh
gp -q -D parisizemax=2000000000 code/margins_census.gp
```

It writes `data/margins-census.txt` and stops on any inconsistency. The file in this package was produced in four parallel parts (the header of the script shows how) and sorted by `Q` and `P`; `checks/independent_margins.py` recomputes every line of it independently.

## The census of directives

The second verification of the finite facts on directives needs Julia (version 1.9 or later) and is not part of the fixed checks. From the supplement directory:

```sh
julia --threads=1 code/directives_census.jl
```

It takes about ten seconds on one thread, writes `data/directives-census.txt` and exits with an error if some Markov number up to 10¹² had a second admissible integer or a second occurrence; `checks/check_directives.py`, which is part of the fixed checks, recomputes every line of that file by another method.
