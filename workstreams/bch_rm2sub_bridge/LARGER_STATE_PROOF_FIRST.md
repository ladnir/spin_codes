# Larger state: prove the distance target before optimizing

**Now closed:** `T64_S20_FULL_CLOSURE.md` states the full t64_s20 theorem.
All8192 occupancies are certified and replayed, and the exact union is
below2^-50. The partial milestones below record the development history.

The user authorized larger RM2Sub states on 2026-09-05. The priority is a
full positive distance proof, followed by parameter optimization. The old
t128_s15 second-moment investigation is retained but is no longer the
immediate route. Its obstruction does not transfer to another map.

## Fixed target and changed maps

Keep the fixed BCH [256,128] outer, 8192 rows, independent row and region
permutations, and output length N=2^21. Keep cutoff H=209716 and target
setup-failure probability below 2^-40. Every message uses the same setup.

The inner recurrence remains

    q_0 = 0,  Y_i = X_i + A q_i,  q_(i+1) = alpha_i q_i + B X_i.

Each alpha_i is an independent uniform nonzero element of GF(2^s).
There is no final flush, and state carries across region boundaries.
Changing s changes the field and maps. Changing t also changes epoch count.
No production encoder or other worktree is modified.

Two concrete map pairs were snapshotted from worktree3061 into the
write-once directory `generated/larger_state_inputs_v1/`:

| Configuration | State bits | Epoch bits | Minimum A weight | Minimum kernel weight |
|---|---:|---:|---:|---:|
| t128_s19 | 19 | 128 | 48 | 6 |
| t64_s20 | 20 | 64 | 16 | 8 |

`larger_state_maps.py` independently checks full rank, B=A^T, BA=0,
distinct nonzero B columns, complete A spectra by Gray-code enumeration,
and complete B-kernel spectra by exact MacWilliams arithmetic. Both pass.
These are algebraic checks, not full distance certificates.

## Positive proof route

For each occupancy Q, let Z_Q count nonzero messages with exactly Q
occupied BCH rows and encoded weight at most H. Seek bounds U_Q with

    E_setup[Z_Q] <= U_Q,   sum_(Q=1)^8192 U_Q < 2^-40.

The union bound then establishes the desired distance guarantee. No
second moment is required for this route. The epoch formulas in
`GENERAL_OCCUPANCIES.md` and `tightened_occupancy.py` accept general t,s
and their exact spectra. Numerical certificates must be recomputed for
the new maps; old t128_s15 occupancy bounds are not reused as values.

The first candidate is t64_s20. The t128_s19 map retains the old epoch
size and is a second candidate. Selection is proof-first, not a claim
that either pair optimizes speed, state size, or distance.

## Development record

`screen_larger_state.py` separates one-row discovery from higher-occupancy
discovery. Its outputs are screens, not certificates. The t64_s20 Q1
screen gives approximately 50.439 bits of margin. A first range screen
at tilt -35 covers Q128 numerically but is inadequate at its dense anchors.
Those nonpassing upper bounds do not establish actual failure.

`certify_larger_state_q1.py` provides a 256-bit outward producer and an
independent 512-bit replay using unnormalized supports and binary powering.
For t64_s20, both runs passed all 92 coefficients and exact BCH aggregation.
The Q1 upper bound is below 2^-50, with diagnostic margin50.4390834639.

The first strong-tilt dense screen was limited by the Bernoulli surrogate
for the all-one band. `screen_larger_state_v2.py` keeps this band exact:
set its Bernoulli parameter to one and its counting cost to one. This
degenerate product measure is exactly the single permuted all-one row.
The existing nonnegative adaptive recurrence remains valid at p=1:
the corresponding update is simply K_j -> K_(j+1). All other bands and
the count of band assignments are unchanged. This is a proof comparison,
not a change to the encoder or message family.

`test_larger_state_constant_band.py` checks the endpoint's counting cost,
enumerates every band sequence on a small exact rational example, and
compares the directed scaled recurrence with exact arithmetic.
`certify_larger_state_range.py` recomputes selected witnesses outward;
its range replay uses 512-bit region calculations. Range certificates
retain at most 80 bits per occupancy to avoid oversized experiment files.

The first exact coverage ledger, `generated/larger_coverage_first.json`,
combines Q1 with every Q from512 through1024. Both range producers were
replayed at512bits, and their selected bounds sum with Q1 to less than
2^-50. This is partial coverage for t64_s20, not a full distance theorem.

The next ledger, `generated/larger_coverage_dense_first.json`, additionally
covers every Q from2048 through4096. All three overlapping range receipts
passed512-bit replay. The complete certified set for t64_s20 is therefore
Q1, Q512..1024, and Q2048..4096; its union contribution is below2^-50.
The remaining occupancies are2..511,1025..2047, and4097..8192.

A dense screen also passes at Q7168 with tilt8. At Q8192, all thirteen
pure-band diagnostic bounds pass at
tilt8, but their adaptive combination does not. Keeping the constant band
exact removes its individual defect but does not remove every loss from
the adaptive maximum. A fixed-composition split may still be needed.

A witness optimized for one Q need not work at neighboring Q. In particular,
the anchor1024 certificate fails at some lower Q in its512..1024 range.
The anchor512 certificate covers those missing rows. The coverage ledger
checks the union explicitly and excludes every nonpassing row. Analogous
low-range attempts are retained, including failures; do not infer coverage
from filenames or from the presence of an outward receipt.

The former remaining ranges are now closed by the gap driver and a
finer tilt at the endpoint. Fixed constant-row boxes were not needed.
Next: independent proof review and a benchmark of this fixed baseline,
then parameter optimization. Fresh-multiplier idealization and the
separate SPIN application obligations remain unchanged.
