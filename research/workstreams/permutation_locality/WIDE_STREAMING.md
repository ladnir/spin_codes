# Wide streaming stores and a late-activation obstruction

2026-09-23. A 16-row candidate improves measured throughput by about 1.22x.
The 2x target and the new distribution's distance certificate remain open.
The faster 32-row candidate fails the required distance. Production code and
the paper are unchanged.

## Implementation result

Retain all 256 regions and share the coordinate permutation within each row
group. Each region independently permutes the groups and independently
permutes the lanes within each group, as defined in ANALYSIS.md. Group sizes
8, 16, and 32 supply contiguous transfers of 128, 256, and 512 bytes.

The new scratch layout interleaves all rows of a group at each BCH coordinate.
The unrolled inner fills an aligned, fixed-size lane buffer. At each group
boundary, full 64-byte streaming stores write this buffer to scratch. A fence
precedes scratch reads. The BCH stage extracts four rows at a time into a
16 KiB tile and calls the existing four-row circuit. All routing metadata and
buffers are prepared before timing; there are no encode-time allocations.

This changes the distribution and the execution layout. The outer generator,
inner maps, recurrence, and parameters remain BCH [256,128] / IMT (128,19).
The kernel's output is checked against the library using the same supplied route.
The production-distribution control is a different sampled linear map.

The initial 51-call screen at K=2^20 gives:

| Group size | Existing tiled kernel (ms) | Wide cached stores (ms) | Wide streaming stores (ms) |
|---:|---:|---:|---:|
| 8 | 8.648 | 14.468 | 8.668 |
| 16 | 8.817 | 13.555 | 7.813 |
| 32 | 8.966 | 12.485 | 7.248 |

The matched production-distribution control takes 9.766 ms. Streaming stores
are essential to the measured result; widening cached stores alone loses.
These timings do not isolate the individual costs of cache allocation,
instruction scheduling, or the final repack.

The 16-row candidate survives the new obstruction below. Confirmation uses
three processes per cell, 101 timed calls after three warmups, and alternating
candidate order. The host is Peach's Ryzen 7950X, pinned to CPU 15. Compiler
flags, timing exclusions, and lock order match COST_AND_CANCELLATION.md.

| Route seed | Production median (ms) | 16-row streaming median (ms) | Throughput ratio |
|---:|---:|---:|---:|
| 1 | 9.592349 | 7.867078 | 1.219x |
| 17 | 9.566902 | 7.860706 | 1.217x |

This is about an 18% latency reduction. Every paired candidate process beats
its control. Two seeds do not establish distribution-wide performance or a
distance guarantee. The matched 2x target remains approximately 4.8 ms.

## Why full zero-state rank was not enough

The 32-row candidate has full rank for the complete zero-state constraints in
the tested setup. That excludes a state that stays zero throughout encoding.
It does not exclude a message that activates the state very late.

The following argument applies to rectangular groups, including c=1. Fix any
realized route and inner randomness. A row group contains g outer rows, each
with 128 message bits. It contributes a contiguous bundle of g*c coordinates
to each of R=256/c macroregions. Assume g*c divides t=128, g*c<=128, and
macroregion boundaries align with epochs. Thus the group occupies one epoch
per macroregion. These conditions hold for the listed experiments.

Restrict messages to this group. For the first ell macroregions, impose
B*x_epoch=0 on its occupied epoch. Each macroregion contributes at most s=19
linear constraints. The resulting message subspace has dimension at least

    h0 = 128*g - 19*ell.

Choose ell<R with h0>0. All messages in this subspace keep the inner state
zero throughout the prefix. Their output there equals their routed outer
word, whose support is contained in ell*g*c positions. In the remaining
macroregions, allow arbitrary output. The union of possible output supports
therefore has size at most

    M = (R-ell)*N/R + ell*g*c.

The encoder is linear and injective. If the actual subspace dimension is h,
each used output coordinate is one for exactly 2^(h-1) of its messages.
Averaging over its 2^h-1 nonzero messages proves that some nonzero codeword
has weight at most

    floor(M * 2^(h0-1) / (2^h0 - 1)).

We use h>=h0 and the fact that the fraction decreases with dimension.
This is an upper bound on minimum distance for every realization in the
specified family, not a probabilistic failure estimate.

At N=2^21, take ell=floor((128*g-1)/19):

| g | c | Prefix macroregions ell | Guaranteed dimension h0 | Relative-distance upper bound |
|---:|---:|---:|---:|---:|
| 4 | 8 | 26 of 32 | 18 | <9.395% |
| 8 | 4 | 53 of 64 | 17 | <8.635% |
| 16 | 2 | 107 of 128 | 15 | <8.286% |
| 32 | 1 | 215 of 256 | 11 | <8.176% |

All four geometries fail the 10% target, including the earlier g=4,c=8
performance candidate. For g=16,c=1 the same bound is about 29.14%, so it
does not decide whether that family meets 10%. It supplies no positive
certificate for the surviving candidate.

The exact local feedback census provides a starting point for that proof.
For one active 16-row group, condition on active-lane count a>0. A uniform
aligned window and uniform lane permutation yield 8*binomial(16,a) equiprobable
supports. None has zero feedback unless a=8. For a=8, exactly two of 102,960
choices have zero feedback, giving probability 1/51,480. This statement concerns
the original uniform lane permutation used in these measurements. It is neither
a bound on cancellation from a nonzero state nor a full-code failure bound.

The obstruction does not depend on a shared coordinate shuffle: the constraint
count also applies when rows have independent coordinate shuffles, provided
the same bundle geometry holds. Enlarging groups cannot be justified by a
full zero-state rank check alone.

## Independently checked witnesses

The experiment now solves the prefix constraints and extracts a nonzero
message. It evaluates this message with the library's packed forward encoder.
It also constructs the BCH outer word independently and checks every prefix
output bit against the routed outer word. For route seed 17:

| g | c | Actual prefix rank | Actual nullity | Verified full output weight |
|---:|---:|---:|---:|---:|
| 4 | 8 | 494 | 18 | 181,460 |
| 8 | 4 | 1,007 | 17 | 174,714 |
| 16 | 2 | 2,033 | 15 | 168,424 |
| 32 | 1 | 3,998 | 98 | 167,816 |

All outputs have length 2,097,152. A second route seed for g=32,c=1 gives
weight 171,056. These concrete witnesses corroborate the obstruction; the
averaging argument establishes it for every route satisfying the geometry.

## Rejected execution schedules

Two other schedules preserve their supplied map but fail the performance gate.

The invalid-distance sequential diagnostic was fused with BCH at input tile
sizes 1024, 4096, and 16384 blocks. It buffers compressed output before copying
it back; otherwise reverse traversal would overwrite unread input. Timings
are 6.438, 6.309, and 6.252 ms, versus 6.166 ms for separate sequential stages.
This simple fusion is not an improvement and is not a valid code candidate.

A state-first schedule saves the exact 19-word reverse state for every epoch,
then gathers contiguous input bundles into BCH tiles. It avoids storing the
full inner output, but takes 27.892, 24.595, and 18.685 ms for g=4 and c=1,4,8.
The corresponding same-route tiled kernels take 8.771, 8.649, and 8.693 ms.
The tested schedule is rejected; setup savings are not counted as encoding wins.

## Replay and next step

Serial scripts are `run_fusion.sh`, `run_state_gather.sh`, `run_wide.sh`, and
`confirm_wide.sh`. Raw receipts remain in the ignored remote measurements
directory. `prefix-rank` uses the same CLI as `rank`, including the optional
final column-bundle argument. For example:

```sh
./build/locality 20 4 1 prefix-rank 17 1 8
./build/locality 20 32 1 prefix-rank 17
```

Every timing process first checks full output and suffix preservation against
the same-route library implementation. `run_sanitize.sh` adds ASan/UBSan checks
at K=2^14 for new kernels and invokes the existing library tests. ASan/UBSan
passes for cached and streaming wide kernels at g=8,16,32 with seeds 1 and 17;
all three fused diagnostics and all three state-gather variants also pass.
The dense transpose and adjoint checks pass for g=16 at both seeds, and all
four existing library tests pass under the sanitizers. The log is
`/tmp/spin-locality-JjY7yR/sanitize.log`. Timing output from these correctness
runs is not performance evidence. Larger-size sanitizer coverage remains
future work; release comparisons cover the measured K=2^20 cases.

Release executable SHA256 after the prefix diagnostic:
`a18fe141d5e29cdf65f6aa73783b64a7367c0a6cb23899a518df0b1d1d6f9a10`.
The production library remains unchanged with SHA256
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

Next, profile the surviving 16-row streaming pipeline and test a SIMD-friendly
within-group shuffle. Changing the lane-permutation distribution requires new
feedback enumeration, even if the one-row marginal stays uniform. Do not spend
more tuning effort on the four geometries excluded above. A full 40-bit margin
certificate for any new family remains missing.
