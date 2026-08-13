# Audit of the conditioned-row packet-width generalization

The binary64 discovery evaluator is sound for `g=8` under its stated
full-support witness scope.  The packet width affects only the packet moment

```text
R_w(t) = [x^w] (sum_j binom(g,j)t_j x^j)^(64/g) / binom(64,w),
```

the profile dimension and mass, and the optimizer coordinates.  The BCH
split spectra, conditioned-row argument, puncture replacement, and graph
average are independent of `g`.  The shared validator restricts `g` to a
supported divisor of 64.

The bounded regression checks established the following facts.

- The implicit legacy path and explicit `group_bits=4` path are bit-identical.
  Both reproduce the frozen scalar `198830.2490738372144`.
- At zero log variables, the `g=8` evaluator returns
  `1048576.0000000002328` bits.  This agrees with the exact total-mass
  identity `log2(2^K)=K=1048576` within binary64 roundoff.
- The symmetric `g=8` optimizer uses eight free log variables followed by
  `band1` and `theta`.  The asymmetric optimizer uses the same eight variables
  followed by `band1`, the interpolation fraction, and `theta`.  Stubbed
  optimizer returns exercised both unpacking paths.
- The generic CLI derives its default per-profile target from
  `-40-log2(feasible_profile_count(g,N,21))`.  For `g=8`, the default is
  `-168.70099010348517` bits.
- The canonical implementation is
  `scripts/probe_packet_group_conditioned_row_outer.py`.  The historical
  `g4` path is now a compatibility wrapper, and internal imports use the
  canonical module.

The independent outward evaluator now accepts an explicit packet width. Its
`g=8` regression checks five frozen vectors, the exact mass identity, and the
canonical `g=4` leader. A separately written 180-digit implementation also
lies inside every frozen outward interval. Complete domain coverage remains
necessary before any `g=8` result becomes theorem-facing.

The upgrader now handles sparse witnesses explicitly. It rejects negative
fugacities and a zero fugacity on an occupied class. A zero fugacity is valid
on an absent class, and the profile pairing sums only over occupied classes.
This rule avoids the undefined binary64 expression `0*(-infinity)` while
preserving the witness's support boundary.
A bounded `g=8` probe produced a finite combined value when the zero class
was absent and rejected the same fugacity vector when that class was occupied.
