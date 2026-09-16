# The two independent computations

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

The two computations of the norm-89 graph and its potential were written
separately and share no code. Running `python3 -B verify.py` performs the
first; `python3 -B verify_independent.py` performs the second, which needs
only the Python standard library, the local `code/unit_reducer.py`, and the
two published data files. It imports nothing from the first computation and
reads no file outside this directory.

They differ at four mathematical steps.

| Step | First computation | Second computation |
| --- | --- | --- |
| Forced reduction | Batched row subtractions | Every unit row subtraction separately |
| Shortest potential | Ordered sweeps over all edges | Relaxation from a work list |
| Zero-slack components | Two-pass postorder | Iterative strong-component search |
| Condensation bound | Forward longest paths | Sink-first remaining budgets |

The second computation also uses the opposite row orientation internally,
converting to the positive-determinant convention only before comparison.
Every reachable state has its two successors regenerated, and the whole
sorted list of states and integer potentials must agree with
`data/axis89-potential.json`; the seed, edge and primitive-terminal slacks
are then checked in full. Fingerprints of the graph and of the potential are
recorded as an extra control, never as the test itself.

The zero-cost language is settled separately by a first-return argument.
Trimming the tight graph leaves six states, one starting seed and one
terminal; removing the terminal leaves an acyclic first-return graph whose
only words are `12111` and `21111`, the initial edge carrying the label
`1`. This covers every finite number of return tokens. Two literal matrix
intertwiners then give the token exchange in the target, and the count
coboundary independently yields `Q - 5P = 0` at every live terminal.

For each cyclic tight component the second computation derives the
rational frequency from fundamental path relations, checks the
count-potential increment on every internal edge, and compares every
normalized rational potential with the first. The resulting `L = 58` gives
`K = 59`, and the full sign-resolution tree confirms that at most two paired
letters are consumed before the graph is entered.

The declared domain is at most 24000 quotient states, each with two
outgoing edges; the relaxation stops inconclusively after one million
removals rather than reporting a partial result. The actual reconstruction
has 11420 states, 22840 edges, 11 seeds, 116 primitive terminals. The exact
mathematical output is in `data/independent-reference.json`.

This settles the finite part of the argument for the fixed reflection
block. It does not establish the balance lemma or the modular-family
argument of the article, both of which are proved there, and it makes no
claim about any other Gaussian block.
