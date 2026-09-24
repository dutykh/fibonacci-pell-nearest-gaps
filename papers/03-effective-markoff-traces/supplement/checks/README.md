# Exact classification checks

This package supplies the finite arithmetic used in the article's small
occurrence and small intersection theorems. It is portable: the paper's
supplement directory is sufficient, and needs nothing from outside it.
All calculations use exact integers or symbolic polynomials. No random input,
floating-point root approximation, or numerical-label uniqueness is used.

## Reproduce the checks

Use Python 3.10 or later. The occurrence and residue programs require only its
standard library; the two intersection polynomial checks and the fixed
two-occurrence endpoint controls also require SymPy. From the supplement directory run:

```sh
python3 -B run_all.py
```

The runner is sequential. Each child has a 30-second wall limit; on Unix it also
has a 128-MiB address-space limit and a 30-second CPU limit. A failed assertion,
missing dependency, resource failure, or mismatched data terminates the runner
without a success claim. The commands work from other directories as well.

Individual examples are:

```sh
python3 -B checks/classify_small_occurrences.py 3
python3 -B checks/verify_small_occurrences.py 3
python3 -B checks/check_intersection_19.py
python3 -B checks/check_intersection_23.py
python3 -B checks/check_residue_cover.py
```

The primary occurrence, the two intersection, and the residue programs accept
`--write` to regenerate their reference JSON. The independent occurrence
checker never modifies data. Running with `python -O` is rejected because the
programs' exact assertions form part of their verification.

## What is checked

| Program | Mathematical role | Reference data |
| --- | --- | --- |
| `check_h2_bounds.py` | Verifies the explicit h=2 analytic endpoint identities and rational inequalities. | Fixed formulas in the program |
| `classify_small_occurrences.py` | Imports the public matrix and marked-seed primitives in `../code/trace_decision.py`; evaluates the proved finite domain and records actual isometries. | `../data/occurrences-h.json` |
| `verify_small_occurrences.py` | Independently generates scalar ordered Markov occurrences, uses fundamental ray recurrences and universal character algebra, and reconstructs every comparison and isometry. | The same occurrence files |
| `check_intersection_19.py` | Derives both literal trace formulas, the elimination polynomial, thirteen univariate positivity margins, all four unbounded sign cones, and the actual surviving isometry. | `../data/intersection-19-polynomials.json` |
| `check_intersection_23.py` | Rebuilds both literal traces in the universal character algebra, confirms the elimination polynomial, the fifteen shifted coefficient lists, the two power-trace identities, the four grouped sign estimates and the two exact rational margins. | `../data/intersection-23-polynomials.json` |
| `check_residue_cover.py` | Checks all unit residue orbits through intersection 24, explicit integral marking maps, reflection matrices, and mechanical words; records the two uncovered orbits at 25. | `../data/residue-cover.json` |

The occurrence reference files cover h=2,3,4,5,6. The h=2 domain was executed
after its analytic endpoint refinements had been independently checked. All five counts now have complete primary and independent
checks: 315692 comparisons in total, 57 equality records, and 50 records
in the original positive-quotient scope. Every equality has an actual isometry
witness; duplicated descriptions of a geometric pair are retained.

## Data conventions

An ordinary trace is three times its positive integral Markov label. In a
minimum seed, `seed=[p,x,y]` means `tr(B)=3p`, `tr(A)=3x`, and `tr(AB)=3y`.
The ray is `u[n]=tr(A B^n)/3`, with its minimum at zero. The word has `h` gaps,
of which `t` are long, beginning at product index `l`. Its literal gap bits are
the successive floor differences of `i*t/h`. Each tested single index is
`n=+m` or `n=-m`, with `q=m-h*l-t`.

`path` identifies an occurrence: `111`, `112`, or `125` followed by its binary
child path. Equal numerical maxima at distinct paths are not merged.
`orientation` records the choice `B=(UV)^(-1)` or `B=UV` at that occurrence.
`source` and `target` are primitive homology columns in the root `(E,F)` basis.
`isometry` indexes the twelve actions in the public module, and the independent
checker derives that list from separately verified physical normalizers.
`character` gives the three ordinary traces of the reconstructed original
marking; the original basis direction is reversed when necessary.

The values `original_k` and `original_t` are null for rows outside the paper's
original k>=1 family. Such equalities are retained to make the normalization
transparent. Otherwise the intersection is `h*original_k+original_t`.
When the original type is null, `character` simply records the negative ray
neighbor; no positive-quotient original marking is claimed for that record.
Tied minima can describe the same geometric pair more than once, so the number
of equality records is not the number of distinct isometry orbits.

Every comparison, including inequalities, contributes its complete exact
input and two traces to a SHA-256 sum modulo 2^256. This order-independent
digest permits different seed traversal orders in the independent checker.
It is a compact comparison aid, not the mathematical reason for completeness:
both programs actually regenerate every row of the proved domain.

The polynomial data contain complete integer coefficient arrays. Each cone
coefficient is stored as `[[degree_A,degree_B,degree_C],coefficient]`.
Univariate arrays are in increasing degree, including zero coefficients.
The defining substitutions and exact program regenerate the lists before
checking their signs. They prove positivity on the whole specified orthant;
they are not samples of a polynomial at finitely many points.

## Scope

The analytic finite domains, with their cardinality and bit bounds, are
recorded in the appendix of the article. The programs here need nothing
beyond this directory, Python and SymPy.

The resulting theorems classify the stated occurrence families and all pairs
through intersection 24. They do not prove global Markov--Frobenius uniqueness,
and the residue file's two uncovered orbits at 25 are not counterexamples.

The release checks were run with Python 3.14.4 and SymPy 1.14.0. Each complete
fixed occurrence run took under one second and approximately 21 MiB on the
preparation machine; the two symbolic intersection checks used about 61 MiB and
69 MiB, and the whole suite finished in 6.6 seconds. These timings are
descriptive, not part of the mathematical claims.
