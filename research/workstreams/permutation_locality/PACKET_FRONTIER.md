# Two-bit implementation and four-bit proof frontier

## Decision (2026-09-28)

Keep the proved two-bit, two-update construction and its padded tiled
implementation as the control. Its confirmed precomputed transpose is
7.422 ms at K=2^20 with 128-bit elements. The whole-code ideal-setup
certificate gives relative distance greater than 9.25% with failure
probability below 2^-49.11; see
[the implementation](TWO_BIT_IMPLEMENTATION.md) and
[the certificate scope](two_bit/FIRST_CLOSURE.md).

Four-bit packets remain a useful candidate: the earlier confirmation
measured 6.442 ms with two updates and 6.727 ms with three. Both use
independent row shuffles, not the older shared-shuffle ensemble.
Neither has a whole-code certificate. A third update improves the
selected proof bounds, but does not resolve the larger occupancies
tested below. Production defaults and the paper are unchanged.
The new same-session confirmation was deferred because another
benchmark held the shared locks. The four-bit timings above therefore
remain the earlier confirmed measurements, not a fresh measurement.

## Two-bit gather experiment

The new diagnostic, mode 6 in `joint.cpp`, retains exactly the same
sampled code as the selected mode 5. Both first write pairs of 128-bit
values into padded buckets. Mode 5 reads each bucket sequentially and
scatters into a local BCH tile. Mode 6 instead uses a precomputed inverse
index to gather four values and write a full cache line into a canonical
four-row tile before applying the same GFNI BCH circuit.

The hypothesis was that sequential full-line writes would outweigh the
local gather cost. The screen rejected it. On Peach's Ryzen 7950X,
GCC 15.2, Release/znver4, CPU 15, K=2^20, two updates, 31 timed calls
after three warmups, the two setup seeds gave:

| Execution schedule | Tile rows | Seed 1 median (ms) | Seed 17 median (ms) |
|---|---:|---:|---:|
| Padded local scatter, selected control | 16 | 7.463 | 7.454 |
| Local gather | 4 | 8.971 | 8.979 |
| Local gather | 8 | 8.882 | 8.951 |
| Local gather | 16 | 9.172 | 9.219 |
| Local gather | 32 | 9.398 | 9.362 |
| Local gather | 64 | 9.988 | 10.008 |
| Local gather | 128 | 13.655 | 13.733 |
| Local gather | 256 | 10.988 | 10.989 |

The phase screen places local restoration plus BCH near 5.44 ms for
eight-row gather, versus 4.08--4.11 ms for the selected scatter. These
phase means include warmups and timer overhead; they are not an exact
decomposition of the process median. No hardware-counter explanation
is claimed. The gather remains opt-in research code, not a replacement.

All 36 new Release cases passed: exponents 14, 18, and 20; seeds 1 and
17; and six tile sizes from 4 to 128 rows. Each checks three input
patterns, dense-inner equivalence, the adjoint identity, complete BCH
output, and suffix preservation. The layout checker also verifies the
inverse gather mapping and rejects malformed routes. The new gather
has not been rerun under sanitizers; the earlier mode-5 sanitizer
results do not constitute such a check.

## Four-bit proof comparison at a common distance

Fix K=2^20, N=2^21, BCH[256,128], and IMT(128,19). A group contains four
adjacent outer rows, each shuffled independently. Here q is the number
of active groups, and u is the common size of the union of their four
shuffled supports. Only these homogeneous support vectors are screened;
they do not exhaust messages at any occupancy q.

The screen regenerates the local operators for each update count. It
uses six-window feedback counts, eight-window output moments,
three-window joint cancellation, and translated-density refinement.
The output tilt is fixed at .056 and the all-one penalty at .9. Only
the auxiliary support tilt is optimized for each point. These are
binary64 proposals for upper bounds on each selected class's expected
bad-message count, including its outer multiplicity and group choices.
There is no outward global replay or complete support cover in this run.

At output cutoff 193986, corresponding to distance 9.25%, the proposed
log2 bounds are:

| Active groups q | Common union size u | Two updates | Three updates |
|---:|---:|---:|---:|
| 80 | 192 | -494.53 | -970.64 |
| 80 | 200 | -288.41 | -774.13 |
| 80 | 208 | -212.92 | -708.81 |
| 96 | 192 | +1321.72 | +795.42 |
| 96 | 200 | +1779.67 | +1240.04 |
| 96 | 208 | +2080.61 | +1527.51 |
| 128 | 192 | +7136.38 | +6501.45 |
| 128 | 200 | +8190.67 | +7537.52 |
| 128 | 208 | +9036.20 | +8364.81 |

The three-update bounds improve by 476--671 bits at these points, but
the tested 96- and 128-group cases remain unresolved. A positive entry
means this upper bound is uninformative; it is not a lower bound on
the number of bad messages, or evidence of an actual low-distance word.

The same run also evaluates cutoff 209715 (10%). With the output tilt
fixed, that change adds exactly .056*(209715-193986)/ln(2), approximately
1270.76 bits, to every entry. This relation is a property of this
Chernoff witness, not a measured distance-versus-margin curve. These
fixed-tilt values must not replace earlier, more heavily optimized
selected-point bounds.

An audit also corrected the root README's stale claim of complete
four-bit coverage through q=128. The documented outward covers reach
q=58, with a combined margin above 43.74 bits for that restricted
message class; see [LOW_OCCUPANCIES.md](independent_rows/LOW_OCCUPANCIES.md).
They concern the two-update four-bit ensemble. The full two-bit
certificate is separate and unchanged.

## Next bounded experiment

The later [four-bit attack](independent_rows/FOUR_BIT_ATTACK.md) retunes
the proof tilt and verifies a three-update selected-event bound below
2^-374.76 at q=96, u=200 and 9.25% distance. This supersedes the positive
fixed-tilt screen for that subcase, not the unrestricted occupancy status.
Larger occupancies and unequal supports remain open.

Update: [SHAPE_POTENTIAL.md](independent_rows/SHAPE_POTENTIAL.md) tests
whole-shape continuation bounds inspired by the companion lifting paper.
The selected 26--29-bit improvements are insufficient for closure. It
also gives an exact two-epoch single-packet return formula and identifies
which packet history the present envelope forgets. A small history-aware
test now precedes the wider retuning proposal below; neither replaces
the proved two-bit fallback.

Retune the output tilt and all-one penalty at q=96 and q=128, at 9.25%
distance, for both two and three updates. This should precede another
full occupancy-cover attempt: increasing updates alone has not closed
the fixed-witness gaps. If the retuned points pass, replay them outward
and extend coverage over heterogeneous supports. If they remain far
from passing, vary the inner batch width or state size and rederive its
operators rather than spending another long run on the same witness.

The four-bit performance advantage is now roughly 9--13% in time
relative to optimized two-bit, not the apparent factor of two from
the obsolete direct-scatter baseline. It is worth investigating, but
does not justify discarding the proved two-bit fallback.

## Reproduction

Run benchmarks serially under the three shared locks. The helper
acquires them before checking or timing; setup and allocation are
outside the measured encoding calls.

```sh
bash research/workstreams/permutation_locality/run_packet_frontier.sh ROOT check
bash research/workstreams/permutation_locality/run_packet_frontier.sh ROOT screen
bash research/workstreams/permutation_locality/run_packet_frontier.sh ROOT confirm 8
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/candidates/round_screen.py --rounds 2 3 --tilts .056 --groups 80 96 128 --supports 192 200 208 --thresholds 193986 209715
```

The candidate suite passes all 49 tests, including update-count
propagation and output-cutoff validation. The cutoff test checks that
changing only the threshold changes only the Chernoff factor.

Raw files remain outside Git. The remote root is `/tmp/spin-joint-9n57TT`;
new Release checks are in `measurements/packet-frontier-check-KZtEV0`
and the screen in `measurements/packet-frontier-screen-NnmGqK`.
The local proof log is `tmp/four-bit-rounds-80-128-d0925-d10.log`.
Screen executable SHA-256:
`dfcbaa8bc6eed5af5c45aedef2239034a58b54dfd1112d4bd5737487f54e47e0`.
Unchanged production library SHA-256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.
