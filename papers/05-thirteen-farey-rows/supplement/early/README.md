# Fixed arithmetic supporting the early-row proofs

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry, France)

The checker in this directory supports the elementary occurrence-cost
dictionary and the fixed source evaluations in Sections 2 and 4.
It evaluates thirteen prescribed reduced slopes, not a denominator range,
and checks their palindromic Gram factors, primitive norms and exact
subtractive Euclidean costs. It also checks the integer margins and
quadratic-field signs printed in the early-row comparisons.

These controls do not prove a symbolic identity or an unbounded envelope.
Those statements are proved in the manuscript. Complete finite Gaussian
fibres and the infinite parity-box covers are supplied by the separate
packages described in the main supplement documentation.

The only dependency is Python 3 and its standard library. Run from this
directory:

    python3 -B check_early_rows.py

The checker accepts no arguments and rejects Python optimization, which
would disable assertions. Its working memory consists of a fixed number
of small integer matrices and thirteen output rows. No binary-word list
or factorization search is generated.

The declared reproduction envelope is 128 MiB and 30 seconds:

    timeout 30s prlimit --as=134217728 python3 -B check_early_rows.py
