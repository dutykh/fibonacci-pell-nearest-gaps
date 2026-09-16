# Optional Magma corroboration

**Authors:** Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)

These files preserve two independent Magma V2.29-9 computations for the
quartics induced by the `j = 1` leading-exponent branch. They corroborate the
published Alekseyev--Tengely classification; the manuscript does not depend
on locally rerunning proprietary software.

| Input | Recorded transcript and check |
| --- | --- |
| `plus_branch_quartic.m` | `plus-branch-quartic.md` |
| `plus_branch_quartic_second.m` | `plus-branch-quartic-second.md` |

The check records state the Magma version, fixed random seeds, full input,
check-relevant output, Mordell--Weil rank bounds, full-group flags, and
the exact filtering from integral points to Fibonacci--Pell pairs. The two
Python filters in `../checks` independently verify the returned point
lists and recurrence filtering using exact standard-library arithmetic.

To rerun a local input in Magma, use the ordinary batch form appropriate to
the installation, for example:

```sh
magma -b plus_branch_quartic.m
magma -b plus_branch_quartic_second.m
```

The online runs recorded in the check files used the University of Sydney's
official Magma calculator. Their completeness claim comes from Magma's
documented `IntegralQuarticPoints` routine together with the recorded full
Mordell--Weil proof information, not from a bounded search.
