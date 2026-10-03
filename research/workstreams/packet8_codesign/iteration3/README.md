# Iteration goal: a complete byte-packet proof

Active goal, 2026-10-02: obtain a complete 10% distance / 40-bit setup-failure
bound at K=65,536, with competitive precomputed transposed encoding of 128-bit
elements. Preserve the measured 91.120 us byte-packed kernel where possible;
the certified four-bit implementation remains the default.

The preceding [checkpoint](../route_conditioning/README.md) separates rare
routing events from message counting. Its strongest tested q=119 bounds
remain negative: -617.49 bits for the baseline 16-bit state and -258.39 bits
for the actual 24-bit state. These are floating proposal bounds, not observed
code failures. The exact routing witnesses alone are not distance certificates.

## Current experiments

| Track | Question | Decision rule |
|---|---|---|
| Conditioned trajectories | Which potential occupancies, actual nonzero occupancies, and state transitions dominate the remaining H2-bound? | Use those counts to choose the next statistic or map change. |
| H3 routing condition | Does counting the first three packets per step exploit the actual 24-bit feedback rank? | First test q=119; expand occupancies only after a credible positive gate. |
| Fractional moment over routing | Can clipping the conditional failure bound before averaging over routing avoid costly shared events? | Check finite-model inequalities and alpha=1 regression before screening. |

The fractional-moment track is a different probability bound, not an encoder
change. For a fixed outer-group subset, it first conditions on potential
occupancy counts. It must average within-step positions before using the
conditional message bound. Its matrix relaxation may introduce path-counting
loss; a negative screen is a reason to inspect that loss, not evidence of a
low-distance word.

## First results

The H2-conditioned trajectory diagnostic separates potential slots from
nonzero byte labels. In the tilted 24-bit comparison at q=119, all 3,808
potential slots remain present, but only about 1,834 carry nonzero bytes.
There are about 1,683 empty zero-to-zero steps and only 0.059 occupied
zero-to-zero steps. Thus the remaining loss is not mainly cancellation of
occupied feedback inputs. Zero byte labels allow long zero-state runs even
when potential slots are dispersed. These are statistics of the tilted
upper-bound expression, not measured failing codewords.

Counting the first three potential slots per step addresses this loss. For
the actual 24-bit maps at q=119, the condition H3 >= 3,183 has an exact
rational bad-route witness of approximately 61.99 bits. At weight tilt 0.5
and route tilt 4, the floating good-message bound is approximately 168.23
bits. Adding the bad-route term leaves approximately 61.99 bits. An
independent evaluator reproduced the log moment within 3e-12, and an
independent review checked the probability split and strict threshold.
The q=128 point is also positive. Neither point establishes full occupancy
coverage, and the message bounds have not been evaluated outward.

For the unchanged 16-bit kernel, fractional averaging improves the q=119
proposal from -617.49 bits with H2 conditioning to -297.46 bits without a
route threshold. Averaging two steps after forgetting their separate
occupancies is worse (-983.47 bits). Retaining both occupancies, multiplying
the two local matrices, and then taking the fractional power improves the
proposal to -180.32 bits. This tests a proof refinement without changing
the measured encoder. Further grouping and H2-conditioned fractional
bounds are being tested separately.

`local_cache24.py` builds source-pinned local-operator receipts from one
actual 24-bit census for subsequent tilt searches. The cached tilt-0.5
operators are bit-for-bit identical to the preceding authenticated receipt.
The cache uses floating arithmetic and is not a certificate backend.

## Eight-step checkpoint: unchanged fast construction

Retaining the fine potential-occupancy tuple while summing all hidden state
paths inside eight consecutive steps gives a positive selected-occupancy
proposal for the original 16-bit maps. At q=119, weight tilt 0.475 and
fractional power 0.35, the margin is **44.9980768686 bits**. The same saved
operator is positive at q=113--159 and gives at least 50 bits at q=120--151.
These intervals are a single-parameter screen, not an optimized full curve.

The following are selected floating q=119 proposals, with separately chosen
proof parameters. The encoder and its setup distribution are unchanged.

| Bound | Margin (bits) |
|---|---:|
| Ordinary first moment | -4015.97 |
| H2-conditioned first moment | -617.49 |
| Fractional, one-step path blocks | -297.46 |
| Fractional, fine two-step blocks | -180.32 |
| Fractional, fine four-step blocks, refined parameters | -39.49 |
| Fractional, fine eight-step blocks | 45.00 |

The eight-step calculation sums 43,046,721 fine-count tuples. Its alpha=1
control agrees with ordinary placement within 8.5e-14 relative error. An
independent review checked the conditioning, tuple weights, chronological
product order, and all probability prefactors. An independent all-log replay
of the saved operators reproduces the selected q=119 margin. The local and
grouped operators remain floating proposals: this is not an outward certificate.

The declared 64 MiB workspace setting in `long_grouped_gate.py` is an
allocation estimate, not a strict peak limit: during replacement, the previous
product chunk can remain alive alongside the new one. This does not affect
the computed values. Preserve the source-pinned run; a future implementation
should free the old chunk before allocating its replacement.

A cheaper four-step survey already shows why checking only q=119 is
insufficient: selected q=4,16,64 bounds are approximately -12.42, -112.71,
and -215.04 bits. The dense cases q=160 and q=256 are positive. The next
question is whether longer path blocks or the actual 24-bit maps close the
smaller-support gap. A single positive middle-occupancy point does not justify
an encoder promotion or a complete 10% / 40-bit claim.

No new encoder timing was taken in this iteration. The 91.120 us value remains
the earlier matched implementation measurement, and the certified 98.950 us
K16 checkpoint remains intact. All 220 of that checkpoint's source pins and
its archived receipt were rechecked successfully.

## Small-support follow-up

The targeted 16-bit eight-step check at q=64 improves the four-step bound
from about -215.04 to -173.45 bits. Longer fixed physical blocks alone are
therefore not a demonstrated solution to the smaller-support gap.

The potential-birth basis combines the old active-count birth measures
using their weighted masses, after zero-label thinning. A stochastic matrix
Q satisfies `Tnew_j Q = Q Told_j`, `Q 1 = 1`, and `e0 Q = e0`. Thus every
unpowered comparison moment is unchanged. Fractional powers introduce less
branch splitting in the new basis. This improves the baseline16 four-step
q=119 proposal to +1.26 bits, but the q=16 and q=64 proposals remain
approximately -111.73 and -203.93 bits. The construction is unchanged.

Combining fractional grouping with the actual 24-bit maps is stronger:
selected, rebased four-step proposals at q=4 and q=64 exceed 65 and 132 bits.
A finer small-support grid gives q=8,16,32 margins of approximately 30.70,
17.35,34.05 bits. These are still below the complete 40-bit target, which
requires summing all occupancies, not just reaching 40 at each tested point.
See [the numerical record](CAPPED_PROGRESS.md) for scope and receipts.

The activity-aligned proof experiment uses blocks containing two nonempty **potential**
steps, with their intervening empty runs. Fixed eight-step blocks often
contain at most one nonempty step at small q. Activity-aligned blocks can
sum the relevant hidden state paths without enumerating longer fixed
tuples. [The recurrence](EVENT_ALIGNED.md) retains chronological order and
continuous state; its alpha=1 case must recover ordinary regional placement.
Long empty runs are evaluated without underflow before fractional powers.
The recurrence and its alpha=1 controls pass, but the selected bounds are
worse: baseline16 q=16 and q=64 give about -117.28 and -266.47 bits;
actual24 q=16 gives 14.72 bits, below the fixed four-step result of 17.35.
Stop this branch unless a stronger reason for a different partition emerges.
The final bounded test used fixed eight-step grouping at actual24 q=16,
with weight tilt 0.06 and fractional power 0.4. It reached **26.752381 bits**,
a 9.40-bit improvement over four-step grouping, but still below 40. The
saved operator is positive only at q=12--20 and has no 40-bit occupancy at
this one proof setting. Its alpha=1 controls cover q=0--512; maximum regional
and global log discrepancies are 2.5e-12 and 6.5e-11. This receipt remains a
floating proposal, not a partial or complete outward certificate.

The [outward-evaluation design](NUMERICAL_CERTIFICATION.md) records a possible
exact-local/native-macro/interval-global verifier. It is a plan, not a claim
that the current floating receipts have been certified. The
[feedforward analysis](FEEDFORWARD.md) remains an unimplemented construction
backup; proof-only improvements currently have the stronger evidence.

## Handoff and next iteration

The goal remains active. This cycle did not close a complete width-eight
certificate. It did establish a stronger conditional-probability argument,
identify and remove artificial birth-family splitting, obtain independently
checked positive middle-occupancy proposals for the unchanged fast kernel,
and localize the harder cases to smaller outer supports.

Do not expand fixed-block enumeration blindly. The actual24 q=16 check still
needs 13.25 bits merely to reach a per-class margin of 40, before reserving
budget for all occupancies. Test a concrete, performance-conscious change
to the local expansion or a demonstrably tighter local comparison first.
Fixed byte scales in A are one untested option; CA=0 is not required by the
original-order recurrence. Any such candidate needs fresh literal maps and
moments, a scoped gate at the obstructing supports, and a cost assessment.
Keep the measured 16-bit design and certified nibble baselines available.

Verification at the end of the cycle: 97 portable tests pass across the
width-eight study and its follow-ups. Current iteration receipts authenticate
against their source pins. Earlier exploratory parent-directory receipts are
historical snapshots, not substitutes for the current source-pinned witnesses.
No production code, paper, timing claim, Git commit, or published artifact
was changed by this iteration.

## Progress discipline

- Keep all previous authenticated numerical sources and certified checkpoints
  unchanged. New experiments live in this directory.
- Save construction identity, parameters, tested choices, source hashes, and
  numerical scope with every result. Do not commit raw experiment output.
- Record substantive improvements and rejected directions here. Do not replace
  a failed construction with a counterfactual denominator adjustment.
- A positive single-occupancy screen is not completion. A complete result must
  cover all nonzero outer-group occupancies and include the routing-failure union.
- Use outward or exact arithmetic for final probability endpoints and obtain
  an independent review before declaring a certificate.
- Do not implement or benchmark a heavier candidate until its proof gate
  justifies the added work. All encoder benchmarks remain serialized.

No production change, promotion, commit, or publication is part of this iteration.
