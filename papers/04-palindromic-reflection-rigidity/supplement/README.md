# Exact controls for the palindromic reflection manuscript

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

This package accompanies the proofs of the finite-word decoder, the palindromic midpoint classification and the displaced common-ancestor bound, together with the Thue–Morse coding, the block exchange and the decoder around an arbitrary palindromic centre, the cubes of nested exchanges, the lengths of directing words, the ladder lemma and the three figures. The fixed checks perform prescribed exact algebra and reconstruct the worked examples. Two further tools enumerate twins among palindromes and test the decoding at arbitrary centres, as reported in the closing section of the article; they give evidence within their stated scopes and are not used in any proof.

## Layout and conventions

- `code/palindromic_words.py`: reference word, right standard-morphism and Gram operations using integers and rational fractions.
- `code/make_examples.py`: regeneration of the five prescribed example rows.
- `code/make_figures.py`: writes the TikZ code of the three figures, `cusp-cylinders.tex`, `nested-square.tex` and `twin-ladders.tex`, from exact data, into the `figures/` folder of the manuscript directory that contains this supplement. It is the only tool that reads or writes outside this folder; the checks and the enumerations do not.
- `code/twin_census.py`: exact enumeration of all palindromes up to a label bound, with the twins, their structure, their exchange classes, their letter counts, the squares of nested exchanges formed by the classes of four and, with `--ladder`, the ladder lemma.
- `code/centre_decoder.py`: all solutions of the reflected column equation at a palindromic centre among words up to a length bound, an independent search that the decoder theorem makes a control.
- `checks/check_decoder.py`: exact symbolic identities and strict residual margins for the decoder.
- `checks/check_examples.py`: independent left-morphism reconstruction and Gram/count/reflection checks for the table.
- `checks/check_displaced.py`: exact identities and constants for the moving reflection and the bound in both parities.
- `checks/check_centres.py`: the Thue–Morse coding, the doubled directive, the block exchange and its arithmetic at every palindromic centre of length at most five, the four palindromes of label 12 986 074 130, and the ladder lemma below label one billion.
- `checks/check_centre_decoder.py`: the decoder at every palindromic centre, with the symbolic identities of its proof for a general symmetric centre matrix, every rational constant of the proof, an exhaustive test on the palindromic centres of length at most four and all words of length at most nine, and the classification at every centre on all palindromes of label at most one billion.
- `checks/check_nested_and_lengths.py`: the length of a directing word as a continued-fraction sum, the bridge to the m-value of Lagisquet, Pelantová, Tavenas and Vuillon on all words of length at most ten, the Christoffel and Farey identification on the central words of directive length at most nine, the cubes of nested exchanges up to depth three, and the first square with its two shifted midpoints.
- `data/midpoint_examples.json`: literal words, ancestor and image Gram entries, and completed counts for the table.
- `LICENSE` and `CITATION.cff`: the licence (GNU LGPL v2.1) and the citation metadata of this material.
- `run_all.py`: fixed portable runner with time limits and, on Unix, memory/CPU limits.

Matrices act on columns and products follow the literal word from left to right. The fixed generators are `A=((2,1),(1,1))`, `B=((5,2),(2,1))`; `P=((2,1),(1,0))` and `e=(2,1)`. The Gram normalization is `P mu(z) P`. Completed counts add one to each letter count. Right morphisms compose in directive order; the independent example checker reconstructs the equivalent left action. All arithmetic is exact. There are no units, random seeds or floating-point tolerances.

## Reproduction

Dependencies: Python 3.10+ and SymPy. Validated versions: Python 3.14.4 and SymPy 1.14.0; the reference word module and example programs use only the standard library. To install the validated symbolic dependency in an existing environment:

```sh
python3 -m pip install -r requirements.txt
```

From this directory:

```sh
python3 -B run_all.py
python3 -B code/make_examples.py
python3 -B checks/check_examples.py
```

The scripts resolve paths from their own locations and also run from another working directory. The generator overwrites only the prescribed JSON output. Each check prints a PASS line and a stable check or data hash. Saved reference output is in `data/reference-output.txt`. A nonzero exit indicates failure. Optimized Python (`-O`) is rejected because exact assertions are part of the checks.

Each child process is limited to 30 seconds and, where the Python resource module is available, 128 MiB of address space and 30 CPU seconds. This is ample for the fixed calculations. No broad enumeration or integer factorization is part of `run_all.py`.

The two enumeration tools are run separately, preferably under a memory limit. The scopes quoted in the article are

```sh
python3 -B code/twin_census.py --bound 1e26        # 16,743,538 palindromes; 8 s, 2.0 GB
python3 -B code/twin_census.py --bound 1e21 --ladder
python3 -B code/centre_decoder.py --max-length 17  # 16 centres; about 2 s, 0.1 GB
python3 -B code/make_figures.py                    # rewrites the three figure files in ../figures
```

The census stops with an inconclusive verdict, rather than exhausting memory, when the number of palindromes below the bound exceeds `--max-palindromes`.
