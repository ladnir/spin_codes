# Closing occupancy gaps without changing the construction

The fixed candidate is BCH [256,128] with RM2Sub t64_s20. Setup samples
independent row permutations, independent region permutations, and a fresh
independent nonzero GF(2^20) multiplier at each epoch. The state starts at
zero, output precedes update, and the final state is discarded. All
messages use the same setup. The target remains K=2^20, N=2^21, output
weight cutoff H=209716, and setup-failure probability below2^-40.

Let Z_Q count bad messages with Q nonzero BCH rows. A first-moment proof
needs an upper bound for every Q, not only a set of sampled occupancies.
`close_larger_state_gap.py` selects numerical witnesses and then checks
every integer occupancy in each attempted interval with outward arithmetic.

## The bound evaluated at each occupancy

Use the region matrices R_j(z) and thirteen BCH weight bands from the
existing larger-state calculation. Each R_j bounds a region containing
exactly j input ones, including the separate activation state. The
coefficient calculation samples positions without replacement.

For an ordinary band g, let U_w be the certified BCH shell cap and choose
p_g in(0,1). Define

    Gamma_g = max_(w in g) U_w / (binom(256,w) p_g^w (1-p_g)^(256-w)).

Choose rho_g >= Gamma_g^(1/256). The all-one band instead uses p_g=1
and rho_g=1, since its counting measure is exactly that product measure.
Starting from K_j^(0)=R_j(z), apply the entrywise nonnegative recurrence

    K_j^(r+1) = max_g rho_g ((1-p_g) K_j^(r) + p_g K_(j+1)^(r)).

The established counting comparison gives

    E[Z_Q] <= binom(8192,Q) 13^Q z^(-H) e_Z (K_0^(Q))^256 1.

The thirteen-band assignment count is retained even when the exact
coefficient hull discards redundant lines. The construction and the
probability space do not depend on the selected numerical witness.

## Discovery, certification, and replay

Discovery chooses a tilt and thirteen Bernoulli parameters at the first
uncovered occupancy. The producer attempts a short interval starting at
that occupancy. It recomputes costs and roots in Arb, uses the directed
scaled recurrence, and evaluates the final matrix power in Arb.

Each result is rounded upward to a power of two, retaining at most80
bits of margin. This limits receipt size without weakening the requested
40-bit target materially. Only upper bounds at most2^-70 enter the gap's
coverage set. Nonpassing rows and all discovery choices remain in immutable
shard files. If the anchor fails, the run stops without claiming closure.

The final range receipt is written only after the set of passing integer
occupancies equals the complete requested interval. Replay checks hashes,
reconstructs the same inequalities at512bits, and compares every result
against its saved dyadic upper bound. Replay does not invoke the optimizer.
It also checks the exact interval union and rational aggregation.

`test_gap_dyadics.py` checks rounding at exact powers, between powers, and
below the80-bit cap. Prior exact tests cover the polynomial region
calculation, scaled recurrence, line hull, and all-one endpoint.

## Current completed gaps

The2..511 gap uses40 short witnesses. The1025..2047 gap uses7 witnesses.
Both producers and their512-bit replays passed. Together with the earlier
certificates, every occupancy1..4096 is covered, with union below2^-50.
The exact ledger is `generated/larger_coverage_through4096.json`.

The4097..8192 coarse search certified through8025, then failed at8026.
`finalize_gap_shards.py` replayed all its shards at512bits and extracted
the complete4097..8025 interval. A separate tilt6 witness certifies every
occupancy8026..8192 and passed512-bit replay. The full ledger is now
`generated/larger_coverage_full.json`: all8192 occupancies, union below2^-50.
The complete statement is in `T64_S20_FULL_CLOSURE.md`.

## Prepared fallback for the dense end

`screen_larger_constant_boxes.py` is discovery only. It fixes an interval
lo<=h<=hi for the number h of all-one rows, rather than allowing that band
inside the adaptive maximum. For d ordinary rows, it replaces R_j by

    Rbar_j = max_(lo<=h<=min(hi,8192-j)) R_(h+j).

The ordinary-band recurrence then uses twelve bands. For each d, the count
binom(8192,d) 2^(8192-d) bounds all allowed choices of ordinary positions
and assignments of zero or all-one rows. The moment envelope is applied
to eligible h first; only their count is enlarged. Summing these bounds
over d covers the fixed h interval. A future certificate must check every
d and cover all required h intervals. No such box is currently certified.

`test_constant_box_envelope.py` compares every entry against explicit
maximization, including truncated endpoints where h+j would exceed8192.
