# One-column groups: stronger bounds at the same encoder speed

2026-09-24. One-column bundles are the leading candidate. They retain the
existing cache-line kernel and take about 5.08 ms at K=2^20. Matched production
controls take 9.58--9.65 ms: a 1.88--1.90x speedup, not 2x.

Outward replays now cover every message with one through thirteen active
four-row groups. Their combined failure contribution is below 3.188222e-15,
giving about 48.1562 bits of margin at 10% relative distance. Fourteen or more
active groups remain uncovered, so this is not a full distance certificate.
The three-/four-group extension is documented in [STATE_MEMORY.md](STATE_MEMORY.md).
The five-/six-group extension is in [ADAPTIVE_OCCUPANCY.md](ADAPTIVE_OCCUPANCY.md);
six groups currently have a 192-bit replay, without the 384-bit repeat.
The seven-through-thirteen extension is in [CDF_COVER.md](CDF_COVER.md).
Nine groups use a 192-bit replay; seven, eight, ten, eleven, twelve, and thirteen also pass
384-bit repeats.

## Change the distribution, not the hot kernel

At K=2^20 there are 8192 BCH[256,128] outer rows. Fix them into 2048 groups
of four. Each group independently samples one uniform 256-coordinate
permutation shared by its rows. Each of the 256 regions independently
permutes the 2048 groups. Each group occupies one four-input window, with an
independent uniform permutation of its four row labels.

Each region contains 64 IMT(128,19) epochs, each containing 32 windows.
The two-column candidate instead places two adjacent columns together in
one eight-input window. Removing that pairing gives each column an
independent group placement across regions.

The physical-column encoder already produces and routes one 64-byte cache
line at a time. It does not need columns to arrive in pairs. The only C++
change enables the existing physical-column schedules for c=1 and rejects
c=1 for the older schedules that require larger bundles. No hot-loop
arithmetic, allocation, or runtime abstraction was added.

## Why this helps the current proof

For the actual feedback map, every pair of distinct four-bit windows has
rank eight: all 496 window pairs are injective on their eight input bits.
Thus two groups sharing an epoch cannot jointly supply a nonzero zero-feedback
input. In contrast, nineteen pairs of eight-bit windows have rank fifteen.

Independent placement also doubles the number of regions in which each
group's support can activate the inner. The verified high-rank bounds improve
substantially without the additional state-memory or row-weight refinements.
These observations do not imply that cancellation from an already nonzero
state is impossible; the transfer bounds still include it.

The calculation uses the same positive seven-coordinate envelope as the
two-column verifier. For two active groups, a common epoch has probability
31/2047. Their local census includes all sixteen ordered pairs of column
weights in {1,2,3,4}. Different-epoch and same-epoch contributions are averaged
with their exact placement probabilities.

The auxiliary Bernoulli coefficient bound in [TWO_GROUPS.md](TWO_GROUPS.md)
now uses degree one per group and 256 regions. A width-two partition of
support sizes covers every support, with additional boundaries at the
minimum supports of the four possible ranks. Each bucket uses the sum of
the four authenticated cumulative rank caps. It does not interpret their
differences as shell counts. The union includes every pair of group locations.

| Covered messages | Failure upper, including all locations | Margin (bits) |
|---|---:|---:|
| Exactly one active group, all ranks | <2.363448e-15 | 48.5880 |
| Exactly two active groups, all rank pairs | <8.179413e-16 | 50.1188 |
| Combined | <3.181389e-15 | 48.1592 |

The combined bound is rounded from the unrounded verified totals. Both
calculations pass independent 192-bit and 384-bit outward replays, with
agreeing bounds. The largest two-group contributions involve supports
around 74--79.

These bounds use ideal uniform setup distributions, not the separate
heuristic permutation bank. The remaining 40-bit failure budget is available
for higher occupancies; it is not a bound on their contribution.

## Matched performance

Peach Ryzen 7950X, GCC 15.2, CPU 15, Release/znver4, 128-bit elements,
K=2^20, IMT(128,19), mask seed 2. Each process performs three warmups and
101 timed in-place encodes without resetting its input. Setup, allocation,
initialization, correctness checks, and checksums are outside timing.

Each entry is the median of three process medians. Mode order reverses in
the middle repetition. All runs are serial under the three shared locks.

| Encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production control | 9.652769 | 9.576125 |
| Two-column candidate | 5.082518 | 5.071166 |
| One-column candidate | 5.078410 | 5.086134 |

One-column process medians range from 5.077298 to 5.100371 ms for seed 1,
and from 5.084061 to 5.120167 ms for seed 17. The measurements do not establish
a speed difference between one and two columns. One column wins on the
current proof coverage. Reaching 2x still needs approximately 5--6% less latency.

Every timed encoder compares its entire output and suffix against a
materialized reference. Raw measurements stay on Peach under
`/tmp/spin-locality-JjY7yR/measurements/one-column`.

ASan/UBSan runs pass at K=2^14 for route seeds 1 and 17. They check the
physical-column encoder, dense transpose, adjoint, direct route, and suffix.
All four existing library tests also pass. The generalized proof helpers
reproduce the preceding two-column single-group bound, 47.3847 bits.

Executable SHA-256:
`147daa285614fc6ddd6ed7661c363271e80472f26c7acdafaf7aaa5d40b8556d`.
Unchanged production-library SHA-256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Negative probes that motivated the change

For two-column bundles, retaining total row weight did not close the
high-rank gap. At support 216 and total row weight 512 per group, the selected
symmetric-point diagnostic still had a log2 contribution near +182.

Remembering the first activation's feedback distribution within a region
improved the coarse full-rank grid from about +203 to +194. Retaining the
number of all-one columns helped further: at support 216 and sixteen such
columns per group, the selected-point diagnostic was about +55. These are
vacuous upper-bound diagnostics, not measured failure probabilities.

The scripts `two_group_weight.py`, `two_group_memory.py`, and
`two_group_flag.py` preserve these experiments. Exact small-code flag counts,
ordered-slot recurrences, and memory-collapse tests pass. Their binary64
output does not certify the missing two-column cases. Smaller bundles give
a simpler and stronger current path.

## Reproduction and next step

From the repository root:

    python -B research/workstreams/permutation_locality/test_two_group.py
    python -B research/workstreams/permutation_locality/two_column_verify.py --precision 192 --columns 1
    python -B research/workstreams/permutation_locality/two_group_verify.py --precision 192 --columns 1 --ranks 0 0 --step 2

Repeat both verifier commands with `--precision 384`. `--ranks 0 0` includes
all ranks, not only rank-zero messages. The scripts authenticate 163 BCH
dependencies and check 67 exact shortening witnesses; they do not rerun
every historical BCH proof.

For the encoder, use `confirm_one_column.sh ROOT` and
`check_one_column.sh ROOT` on the configured Peach build. Production sources
and the paper remain unchanged.

The extensions below and in the linked follow-on notes now cover through
six active groups. Larger occupancies remain open, as does the remaining
encoder-latency gap. Multiple groups in an epoch require joint feedback
bounds; two-window injectivity alone does not close those cases.

The initial exact census for that extension is already available:

| Groups in one epoch | Rank histogram of window subsets | Worst zero-feedback probability over fixed input-weight tuples |
|---|---|---:|
| 3 | rank 11: 27; rank 12: 4933 | 11/1428480 |
| 4 | rank 14: 80; rank 15: 4180; rank 16: 31700 | 7/863040 |

These probabilities condition on all indicated groups sharing one epoch,
with distinct uniform windows and independent uniform lane shuffles at
each fixed column weight. They do not include the placement probability
of that common epoch. The three-group maximum occurs at weights (2,1,1);
the four-group maximum occurs at (4,4,4,2), up to relabeling.

`window_feedback.py` enumerates every window subset and its kernel vectors,
retaining only inputs nonzero in every group. Direct small-map enumeration
checks the kernel construction. The census counts 27 fully active kernel
masks across triples and 3637 across quadruples. These are local activation
counts, not distance bounds for those occupancies.

The follow-on recurrence, direct-coefficient screen, and cancellation
diagnostics are in [MULTI_GROUP.md](MULTI_GROUP.md). That earlier three-group
support cube missed 40 bits after summation. Tracking feedback history
instead closes the full support ranges for three and four groups without
changing the hot kernel; see [STATE_MEMORY.md](STATE_MEMORY.md).
