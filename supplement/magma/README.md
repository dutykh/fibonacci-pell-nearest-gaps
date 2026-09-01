# Optional Magma corroboration

These files preserve two independent Magma V2.29-9 computations for the
quartics induced by the $j=1$ leading-exponent branch. They corroborate the
published Alekseyev--Tengely classification; the manuscript does not depend
on locally rerunning proprietary software.

| Input | Recorded transcript and audit |
| --- | --- |
| `check_c0068_plus_j1_quartic_curator.m` | `2026-08-23-c0068-plus-j1-quartic-curator.md` |
| `check_j1_plus_branch_quartic_independent_audit_agent_j1.m` | `2026-08-23-j1-plus-branch-quartic-independent-audit-agent-j1.md` |

The audit records state the Magma version, fixed random seeds, full input,
certificate-relevant output, Mordell--Weil rank bounds, full-group flags, and
the exact filtering from integral points to Fibonacci--Pell pairs. The two
Python filters in `../certificates` independently verify the returned point
lists and recurrence filtering using exact standard-library arithmetic.

To rerun a local input in Magma, use the ordinary batch form appropriate to
the installation, for example:

```sh
magma -b check_c0068_plus_j1_quartic_curator.m
magma -b check_j1_plus_branch_quartic_independent_audit_agent_j1.m
```

The online runs recorded in the audit files used the University of Sydney's
official Magma calculator. Their completeness claim comes from Magma's
documented `IntegralQuarticPoints` routine together with the recorded full
Mordell--Weil proof information, not from a bounded search.
