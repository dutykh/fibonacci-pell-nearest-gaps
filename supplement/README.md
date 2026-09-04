# Reproduction supplement

This supplement supports the exact computer-assisted statements in
`DD-LV-Fibonacci-Pell-Gaps.tex`. It is self-contained relative to the
manuscript directory: the driver reads no source, data, or audit file from the
surrounding research repository. Every checker uses only the Python standard
library. The proof-bearing interval and enumeration decisions are made with
exact integer, rational, or outward-rounded fixed-point arithmetic.

## Run the certificates

From `manuscripts/fibonacci-pell-nearest-gaps`, run:

```sh
python3 -B supplement/run_all.py
```

The driver stops at the first failure and prints the command and its complete
output. It resolves every executable below `supplement/certificates`, so the
manuscript directory can be copied and checked independently of the research
repository. A successful run ends with the line

```text
PASS: every manuscript certificate completed successfully.
```

The complete expected summaries are recorded in `expected-output.txt`.

## Claim-to-source map

| Manuscript claim | Complete checker | Independent control |
| --- | --- | --- |
| First Matveev comparison and $h<2.4\cdot10^{14}(1+\log(2q))$ | Written proof and source constant regression `certificates/check_c0082_logarithmic_h_bound_curator.py` | Exact rational reconstruction in the logarithmic-bound part of `certificates/check_c0082_effective_advances_independent_audit_agent_ea.py` |
| Moving coefficient, full norm, second Matveev constants, and $q<10^{32}$ | `certificates/check_c0082_full8_adjusted_form_agent_c.py` | `certificates/check_c0082_full8_complete_exclusion_independent_audit_agent_fra.py` and `certificates/check_full8_adjusted_dp_independent_agent_b.py` |
| Common rational approximation, $h<192$, all $190$ moving shifts, $q<90$, and the two-row enumeration | `certificates/check_c0082_full8_dp_closure_agent_c.py` | `certificates/check_c0082_full8_complete_exclusion_independent_audit_agent_fra.py` and `certificates/check_full8_adjusted_dp_independent_agent_b.py` |
| All-exponent classification and empty odd terminal domain $q=3,6,\ldots,93$ | `certificates/check_all_exponent_nearest_gap_agent_aeng.py --q-max 2000` | `certificates/check_all_exponent_nearest_gap_independent_audit_agent_aei.py --q-max 2000`; the second checker independently reconstructs the integer arithmetic and does not import the first |
| Exact densities of both, one-sign, and absent candidate windows | Written irrational-rotation proof and `certificates/check_candidate_window_density_curator.py --q-max 5000` | The finite counts are regression evidence only; the density theorem is unbounded because of Weyl equidistribution |
| Complete all-anchor two-sign Fibonacci-core orbit list | `certificates/check_two_sign_odd_pell_orbit_independent_audit_agent_tsa.py` and `certificates/check_even_pell_anchor_orbits_agent_epa.py --q-max 5000` | `certificates/check_even_pell_anchor_orbits_independent_audit_agent_eaa.py --q-max 5000`, together with the all-exponent pair; bounded orbit searches are not the completeness proof |
| General-seed semiconjugacy, canonical units, and anchored orbit/gap equivalence | `certificates/check_general_seed_orbit_nearest_gap_curator.py` | `certificates/check_general_seed_orbit_nearest_gap_independent_audit_agent_gsa.py` |
| Uniform simultaneous local clocks, rank-parity criterion, and the obstruction at $241$ | Written rank/CRT/equidistribution proof and `certificates/check_uniform_local_clock_generalisation_agent_lcg.py` | The checker verifies exact ranks and witnesses; infinitude rests on the written proof |
| Wider exact nearest-gap regression | `certificates/check_c0068_nearest_even_unit_gap_curator.py --q-max 2000` | Not used for completeness |
| $7$- and $17$-clock theorem witnesses and rank data | `certificates/check_c0077_imprimitive_two_clock_recurrence_agent_cte.py` | `certificates/check_c0077_imprimitive_two_clock_recurrence_independent_audit_agent_cia.py` |
| Two-arms section: Cohn normalisation, arm letter counts, branch recurrence with trace $3m_0$, failure of a half-step on the displayed fixed-$34$ ray, and the growth-rate trace comparison | `certificates/check_two_arms.py` | Self-contained; the two trace-criterion scans are finite verifications, and neither they nor the displayed-ray calculation classify all fixed-$34$ rays |
| Complete fixed-value enumeration of primitive two-square writings and their candidate Stern--Brocot itineraries | Written Gaussian unique-factorisation and subtractive-Euclid argument in Remark 10.2 | No finite search supplies completeness; Cohn evaluation and the central-word test validate candidates but do not replace the Gaussian enumeration proof |
| Exact filtering of the published $j=1$ near-square classification | `certificates/check_c0068_plus_j1_quartic_curator.py` | `certificates/check_j1_plus_branch_quartic_independent_audit_agent_j1.py` |

The identifier `C0085` in the recorded output of
`check_c0082_full8_adjusted_form_agent_c.py` names the research-tree
derivation of the Matveev prefactor $9.275\cdot10^{12}$. That prefactor is
re-derived here from Matveev's constants with exact rational intervals inside
`certificates/check_c0082_effective_advances_independent_audit_agent_ea.py`,
and the specialisation it rests on is the Matveev lemma of the manuscript, so
the condition recorded in that PASS line is discharged inside this supplement
and does not depend on any file outside it.

The other internal identifiers that survive in shipped file names and PASS
lines are research-tree claim numbers, not manuscript objects. They may be
read as follows: `C0060` and `C0065`, the primitive factorisation and
orientation congruences of the even branch; `C0068`, the nearest even-unit
gap criterion; `C0070` and `C0073`, the half-unit trace conic and its
consequences; `C0077`, the imprimitive two-clock recurrence; `C0080`, the
exclusion of primitive norm one; `C0082`, the effective classification;
`C0088`, the all-exponent classification itself. Nothing in the manuscript
depends on these labels; they are retained only so that the vendored files
remain byte-identical to their sources.

The run instructions and provenance pointers inside individual vendored
files mention `math-sandbox/`, `scripts/` and `attempts/` paths. Those are
historical annotations from the research tree and are not part of the
distributed supplement. The only supported invocation is the one given under
"Run the certificates" above.

The independent effective-advances checker also reconstructs an earlier
$4+4$ support exclusion. That extra check is retained as provenance, but the
present manuscript uses this program only as an independent audit of the
first Matveev comparison.

The wider nearest-gap runs through $q=2000$, the density run through
$q=5000$, and the bounded orbit comparisons are regression checks, not the
source of completeness. Completeness comes from the manuscript's unbounded
estimates followed by the exact even domain $q<90$ and odd domain
$q=3,6,\ldots,93$, and from the exact anchored orbit/gap equivalence.

## Optional Magma corroboration

The leading-exponent lemma is proved by the published complete classification
of Alekseyev and Tengely, *Journal of Integer Sequences* 17 (2014), Article
14.6.6, Table 1. Magma is therefore not required for the manuscript proof or
for `run_all.py`.

The `magma` directory preserves two independent Magma V2.29-9 inputs and the
two audit records containing their complete certificate-relevant transcripts,
fixed seeds, version information, integral-point outputs, descent bounds, and
full Mordell--Weil proof flags. See `magma/README.md` for the exact map.

## Arithmetic conventions

- $F_0=0$, $F_1=1$, and $F_{n+2}=F_{n+1}+F_n$.
- $P_0=0$, $P_1=1$, and $P_{n+2}=2P_{n+1}+P_n$.
- $C_n=P_n+P_{n-1}$, with $P_{-1}=1$.
- $\lambda=1+\sqrt2$.
- A quadratic integer is represented by its two integer coefficients.
- Every logarithm and square root that decides the effective classification
  is enclosed by rational endpoints or by an outward-rounded fixed-point
  interval. The source first-Matveev regression uses 80-digit `Decimal`
  diagnostics, but its inequalities are separately reconstructed with exact
  rational intervals by the independent effective-advances checker.
- The source local-clock checker uses binary floating point only to locate
  displayed witness candidates and then verifies their nearest windows over
  the integers. The independent local-clock checker reconstructs the displayed
  rows without that locator, while the unbounded clock theorem rests on the
  written recurrence and equidistribution proof rather than on finite search.

## Software and provenance

- Python 3.11 or later is recommended; no third-party package is required.
- The optional integral-point reruns use Magma V2.29-9.
- The certificate files are unmodified snapshots of the audited research
  scripts. `PROVENANCE.md` records their source paths and SHA-256 digests.
- The Python certificates have no nonlocal dependencies, invoke no network
  service, and use no nondeterminism. Two new checkers import named sibling
  certificates already vendored in the same directory; all other checkers
  import only the standard library.
- LaTeX compilation is separate and uses `make check` in the manuscript
  directory.
