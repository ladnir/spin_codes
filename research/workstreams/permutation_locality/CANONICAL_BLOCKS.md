# Contiguous column blocks: approximately 1.32x throughput

2026-09-23. A new block-preserving route eliminates the BCH repack and improves
full precomputed transposed encoding to approximately 7.2 ms at K=2^20.
Matched production controls take approximately 9.5 ms. This is progress toward
the 2x target, not its completion. The new distribution has no distance
certificate yet; production code and the paper remain unchanged.

## Confirmed performance

Each entry is the median of three process medians, each with 101 in-place
encodes after three warmups. Candidate order reverses in the middle repetition.
Runs are serial under the shared benchmark locks, pinned to CPU 15 on Peach's
Ryzen 7950X. The workload uses 128-bit elements, K=2^20, BCH [256,128], and
IMT (128,19). Coefficient seed is 2. Setup, allocation, initialization, and
correctness checks are outside timing; inputs are not reset between calls.

| Encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production distribution and kernel | 9.523886 | 9.497617 |
| Previous 16-row streaming candidate | 7.888275 | 7.834886 |
| Four-row/four-column blocks, original BCH circuit | 7.375118 | 7.372763 |
| Same route, nonaliasing BCH emission | 7.175114 | 7.220900 |

The last row gives 1.327x and 1.315x throughput, or 24.66% and 23.97% lower
latency than the respective controls. Its process-median ranges are
7.174073--7.252589 ms and 7.208407--7.226991 ms. The previous candidate's
rank-one partial certificate does not apply to this new route.

## Distribution and layout

Fix groups of four outer rows. Partition the 256 BCH coordinates into 64
canonical consecutive blocks of four coordinates. In each row group, sample
a uniform permutation of these blocks and independent uniform permutations
of the four coordinates inside each block. All four rows share these choices.

There are 64 macroregions. In each macroregion, independently permute the row
groups. Each group contributes its assigned four-coordinate block from all
four rows. Independently permute the row labels inside each of those four
columns. All sampling choices above are independent unless sharing is stated.

Thus each group contributes one 16-input window per macroregion. Crucially,
the corresponding destinations are also contiguous in the canonical four-row
BCH input layout. The inner fills a 256-byte aligned buffer, then writes four
cache lines with streaming stores. The BCH stage reads this layout directly;
there is no repack and no new arithmetic in the inner recurrence.

This is a stricter coordinate-permutation distribution than the previous
rectangular experiment. A fixed word's active macroregions depend on which
canonical four-coordinate blocks it occupies. They are not determined by its
ordinary binary weight alone. Neither the original one-row proof nor the
16-row shared-permutation partial proof can be transferred without a new
argument.

The two-column version was also screened. With eight values per transfer,
it took 7.960 ms, versus 7.340 ms for four columns in the same screen.
Cached stores were substantially slower: 14.372 and 13.814 ms respectively.
The retained candidate uses streaming stores and a fence before BCH reads.

A separate phase diagnostic for the four-column candidate with the original
BCH circuit gave 3.376 ms for inner/routing and 4.007 ms for BCH. These are
instrumented means including warmups, not the confirmation medians. Removing
repacking helps, but the remaining BCH arithmetic and memory traffic are now
the larger measured phase.

## Exact-map BCH scheduling experiments

`schedule_bch.py` retains the existing XOR gates and output expressions. It
verifies all 128 outputs as formal 256-bit linear forms against `BchRows`
before emitting a build-local C++ translation unit. No production source is
rewritten. All variants retain fixed-width SIMD and direct calls in the BCH
loop; the schedule switch occurs outside that loop.

The 101-call, seed-17 screen gives:

| BCH emission within the new block route | Full time (ms) |
|---|---:|
| Original library circuit | 7.410 |
| Greedy peak-live-signal order | 10.379 |
| Greedy release-first order | 10.628 |
| Reverse output order | 7.451 |
| Original order, nonaliasing emission | 7.221 |

The two greedy schedules reduce a symbolic live-signal count but lose badly
in actual execution. Keep them as rejected experiments, not optimizations.
The selected emission uses `__restrict` input/output pointers and emits the
body directly in one non-inlined function. The library version uses a small
wrapper around a non-inlined helper. This screen does not separate the effect
of the alias promise from that emission difference. Scratch input and output
are disjoint in this caller, so the alias promise is valid for this path.

## Proof filters, not a certificate

For one four-row group there are 512 message bits. The sampled full zero-state
constraint matrix has rank 512 at seed 17. Keeping the state zero through the
first 32 macroregions gives rank 505 and a seven-dimensional kernel. A kernel
message independently checked by the forward encoder has output weight
513456 of 2097152, about 24.5%, so this witness does not refute the 10% target.
The tested prefixes of 40,48,52 macroregions all have rank 512.

These tests concern one sampled route and group. They do not bound failure
probability or exclude other kinds of low-weight output.

The exact local cancellation census now conditions on the sorted weights of
the four input columns, each of height four. Uniform column permutations and
independent row permutations are transitive on each such shape. Enumerating
all 65536 masks in each of the eight aligned windows gives:

| Column-weight shape | Zero-feedback choices | Total choices | Probability |
|---|---:|---:|---:|
| (0,2,3,3) | 1 | 9216 | 1/9216 |
| (1,2,2,3) | 1 | 55296 | 1/55296 |

Every other nonempty shape has zero probability of zero feedback. The worst
local probability is larger than the 1/51480 bound for a uniform 16-lane
shuffle. A new proof must count occupied canonical blocks and their column
shapes, not substitute that earlier probability or its BCH support law.

## Reproduction and validation

`run_canonical.sh` screens the layout; `run_schedule.sh` screens BCH emission;
`confirm_canonical.sh` performs the two-seed confirmation. Each script acquires
the same three shared benchmark locks and runs measurements serially. Raw
receipts remain outside version control.

Every timed candidate first compares the complete in-place output, including
the unchanged suffix, with the library using its supplied route. Dense and
adjoint checks pass at K=2^14 for the four-column candidate. The local census
is replayed by `python -B research/workstreams/permutation_locality/canonical_census.py`.
ASan/UBSan passes for both two- and four-column routes at K=2^14, seeds 1
and 17, including the cached and streaming kernels. All four new BCH emissions
also pass at K=2^14, seed 17. All four existing library tests pass under the
sanitizers. These checks do not constitute a full-size sanitizer campaign.
The remote log is `/tmp/spin-locality-JjY7yR/sanitize-canonical.log`.

The confirmed release executable has SHA-256
`fb7fa2001da7a44806a96c37cefc1b3245996798a24b2027c7ca6f754326ee89`.
The production library hash remains
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Canonical-block counting and delayed activation

The next check bounds how long a nonzero message in one four-row group can
leave the state identically zero. It concerns the ideal uniform setup
distribution specified above, not just the two benchmark seeds. The IMT
transvections may be fixed arbitrarily: they preserve a zero state.

Let h be the binary rank of the four BCH words, and u the number of canonical
four-coordinate blocks touched by their union. Write T_h(u) for the number of
ordered four-tuples of rank h touching at most u blocks, including zero rows.
Let k_v upper-bound the dimension of every BCH subcode supported on v fixed
coordinates. The earlier verified shortening bounds give

\[
T_h(u)\leq {64\choose u}R(k_{4u},4,h),\qquad
R(k,g,h)=\prod_{j=0}^{h-1}
\frac{(2^k-2^j)(2^g-2^j)}{2^h-2^j}.
\]

The right side is zero when k_{4u}<h. This counts all tuples in each possible
u-block shortened code; overcounting is harmless. We also take the minimum
with the earlier bit-support CDF bound at 4u bits and the total R(128,4,h).
Monotonicity lets later CDF caps improve earlier ones. We do not subtract CDF
caps to infer shell counts. These bounds require respectively at least
10, 15, 17, and 18 occupied blocks for ranks one through four. They do not
assert that tuples attaining those minima exist in this BCH code.

Rank one has a stronger feedback property. Its rows lie in {0,c} for a
nonzero BCH word c. If a rows equal c, every column has weight either zero
or a. Neither zero-feedback shape in the census has this form. Thus a
nonempty block always activates the state when the incoming state is zero.
For higher ranks, use the uniform upper bound q=1/9216.

Fix a tuple touching u blocks. Its occupied macroregions form a uniform
u-subset of the 64 macroregions. Conditional on that subset, group positions
and the local shuffles in distinct macroregions are independent. If J of
those macroregions occur in a prefix of length ell, remaining identically
zero throughout that prefix requires zero feedback at all J active epochs.
Consequently its probability is at most

\[
f_{\ell,q}(u)=\frac{1}{{64\choose u}}
\sum_j {\ell\choose j}{64-\ell\choose u-j}q^j,
\]

using q=0 for rank one and interpreting 0^0 as one. Each f is decreasing in
u. Summation by parts therefore combines f with upper CDFs using only
nonnegative coefficients. Finally, union-bound over the 2048 possible active
four-row groups and all four ranks. This does not require independence
between events for different messages or groups.

Exact-rational computation gives the following *untruncated* log2 upper
bounds. Positive entries are vacuous probability bounds.

| Prefix in macroregions | All ranks | Rank one |
|---:|---:|---:|
| 48/64 | 68.998 | 28.622 |
| 52/64 | 21.520 | 21.161 |
| 55/64 | -17.975 | impossible |
| 56/64 | -30.427 | impossible |
| 57/64 | -42.886 | impossible |
| 58/64 | -55.354 | impossible |
| 60/64 | -80.307 | impossible |

At 57 macroregions, the code checks the strict inequality against 2^-42
using integers and rational arithmetic; the displayed logarithm is only a
rounded summary. All accepted shortened-code LP witnesses are verified
exactly. The calculation authenticates the same 163 dependency files as the
earlier support calculation; it does not rerun every historical BCH proof.
Small-code exhaustive tests check the block-support count and its weighted
CDF bound. Reproduce with `python -B research/workstreams/permutation_locality/canonical_support.py`.

This is not a 10%-distance certificate. It only bounds a completely zero
state prefix for messages in one active group. Activation after almost 90%
of the encoder, later returns to zero, output weight, and multiple active
groups still require control. The bound has not ruled out this distribution,
but it is substantially weaker than the useful full-distance statement.

## Further BCH sharing experiment

Increasing the minimum reuse threshold from four to six or eight reduces
common intermediate values but adds XOR work. The build generator reuses
the existing Paar synthesizer and verifies all 128 output linear forms
against BchRows. Both variants use the same four-row packing, original
output order, and nonaliasing input/output pointers as the current winner.
There is no runtime dispatch inside the row loop and no change to the code
being evaluated.

The same serial two-seed, three-process, 101-call protocol gives these
medians of process medians:

| Implementation | Seed 1 (ms) | Seed 17 (ms) |
|---|---:|---:|
| Production control | 9.629389 | 9.659025 |
| Canonical blocks, current nonaliasing BCH | 7.238507 | 7.247174 |
| BCH reuse threshold 6 | 7.297558 | 7.347140 |
| BCH reuse threshold 8 | 7.431037 | 7.438801 |

Neither variant improves performance. Retain the current circuit. The
roughly 1.33x ratio confirms the earlier 1.32x result within run variation;
it is not a new optimization gain. Every timed candidate passed the
full-output reference comparison. `run_sharing.sh` reproduces the serial
comparison; measurements remain in the remote scratch directory.

The threshold-six circuit has 4212 XORs and 240 shared gates; threshold
eight has 4539 XORs and 180 shared gates. Their generated stack frames,
excluding alignment and the saved frame pointer, are 10760 and 11592 bytes.
Both are smaller than the nonaliasing threshold-four frame of 12936 bytes,
but the reduction does not offset the added work in these measurements.
ASan/UBSan checks pass for both new circuits at K=2^14, seed 17, together
with the earlier routing checks and all four library tests. The log is
`/tmp/spin-locality-JjY7yR/sanitize-sharing.log`.

The updated release executable hash is
`16759726149e2a2972d5c81f5f4384dada42717d634461419ee85883dfcd0988`;
the production library hash above is unchanged.

The 2x goal still needs roughly another third removed from the candidate's
latency. Small changes to sharing and output order have not delivered that.
Next, screen a different BCH evaluation organization with bounded working
storage before spending more time tuning these circuits. In parallel with
that design work, a full-distance analysis needs a sharper bound on
canonical-block support or a distribution that admits stronger counting
without losing the contiguous transfers. Production defaults and the paper
remain unchanged.
