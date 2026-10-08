# One permutation with a sequential precode

**The striped GRS[128,64] candidate with two scalar weighted scans is
withdrawn.** A deterministic parity-check argument gives minimum distance
at most 832321, below the required 838861 at K=2^20 and N=2K. Thus every
setup fails the 40% target, regardless of its field or coefficients.
The argument is stronger than the earlier probabilistic band screen.

The source is the companion
`gen-BAA/field-universality/expander-lift/CONVOLUTION_WIDTH_ONE_CEILING.md`.
Each inverse scalar scan has at most two entries per column. The two scans,
diagonals, stripe, and intervening permutation therefore leave each final
parity-check column incident to at most four outer check groups. Here each
group contains 64 checks and there are 16384 groups. Selecting 13005 groups
satisfies the lemma's strict counting inequality, giving

    d_min <= 64*13005+1 = 832321 < ceil(0.4*2^21) = 838861.

The retained map and cost analysis below now describe a **rejected control**,
not a candidate certificate. Native encoding also misses the 5ms target:
the completed 101-call baseline median is 12.152ms. Each application record
contains 128 bits. The author of this note ran no builds or benchmarks;
the measurements below were supplied by the parent task's serial run.

## Concrete map

Fix F=GF(2^16), using the existing quadratic extension of AES GF256.
Let R=K/64. Apply the fixed normalized additive-basis GRS[128,64] encoder
to each message row x_i. Write its coordinates in striped order:

    z[j*R+i] = GRS(x_i)[j],  0 <= j <128, 0 <= i <R.

This layout is part of the construction. It is not another sampled
permutation. At K20, each coordinate plane has R=16384 records and 256KiB.
The complete striped vector has N=2^21 records and 32MiB.

Setup samples one uniform permutation pi of {0,...,N-1}. Independently,
it samples d_i,e_i from F minus {0}, separately for every position.
Define the two scans by

    u[-1] = 0;    u[i] = u[i-1] + d[i]*z[i],
    y[-1] = 0;    y[i] = y[i-1] + e[i]*u[pi[i]].

The encoder outputs all N entries of y. The sampled setup remains fixed
for every message. Every setup is injective: the outer has rank K, and
the diagonals, permutation, and triangular scans are invertible.

These are input-weighted scans. Their field multiplications are independent
within a batch; only XOR carries connect adjacent outputs. Do not replace
them by an unproved finite-alphabet feedback-scan ensemble. A gauge change
can relate certain full-field feedback scans, but its distributions and
absorbed diagonals must be derived before importing a certificate.

Eight F coordinates fit in one 128-bit record. All eight use the same
scalar coefficient. The minimum nonzero-record distance equals the scalar
F-code distance: each component is a codeword, and one nonzero component
can realize a minimum-weight scalar codeword.

## Why the layout matters, and what it does not fix

A contiguous-block precode does not cure the earlier window obstruction.
For an outer [b,k] block, forcing an r-state precoder to zero at its end
costs at most r equations. A subcode of dimension at least k-r remains
localized to that block. A subsequent uniform permutation and causal scan
therefore retain vulnerable shortened subspaces.

Striping makes one row occupy b separated positions. However, a band of
a adjacent rows still gives a localized subspace. For the width-one scan
above, the band occupies b intervals, each containing a positions. Force
the scan state to zero at the last position of every interval. The b
constraints leave dimension at least

    d_band = k*a-b.

Those b last output positions are zero by construction. Excluding them
gives a support container of size

    C = b*(a-1).

Let J be a terminal suffix containing W=ceil(delta*N)-1 output positions.
Let T count the container positions mapped by pi into J. Restricting the
band subcode to pi^-1(J) leaves dimension at least

    d_band-(C-T) = T-(b-k)*a.

Thus T >= (b-k)*a+1 forces a nonzero final codeword supported in J.
No final-state equation is needed: a causal postscan preserves suffix
support. This argument applies to every choice of d and e, and to any
causal postscan width. Increasing its field or refresh count cannot remove
this event.

For one band, T has the exact hypergeometric distribution

    Pr[T=t] = choose(W,t)*choose(N-W,C-t)/choose(N,C).

Disjoint row bands give disjoint containers. Their occupancies are not
independent. An elementary lower bound uses only two-term inclusion-exclusion.
For G=floor(R/a), set

    p  = choose(W,t)*choose(N-W,C-t)/choose(N,C),
    p' = choose(W-t,t)*choose(N-W-C+t,C-t)/choose(N-C,C).

Then failure probability is at least G*p-choose(G,2)*p*p'. All these
quantities are rational. The following screen used exact integer arithmetic;
the logarithms are approximate displays. Only one occupancy value t was
included, so a tail sum can strengthen each lower bound.

| Outer | Distance target | Band a | C / t | log2 failure lower bound |
|---|---:|---:|---:|---:|
| [16,8] |40% |5 |64 /41 |-1.36 |
| [32,16] |40% |5 |128 /81 |-10.41 |
| [64,32] |40% |6 |320 /193 |-30.665 |
| [128,64] |40% |6 |640 /385 |-69.616 |
| [128,64] |45% |10 |1152 /641 |-32.519 |
| [256,128] |45% |11 |2560 /1409 |-71.225 |

The first three rows rule out the 40%/40-bit goal for these width-one
striped constructions. The [128,64]/45% row rules out the stretch target
for that candidate. A separate eight-term tail calculation strengthens
the [64,32]/40% failure lower bound to approximately 1.0338e-9.

The [128,64]/40% and [256,128]/45% rows were only inconclusive necessary
screens. The deterministic degree-four ceiling above supersedes them for
this two-scalar-scan topology. A fixed cheap stripe is not equivalent to
the missing first uniform permutation in a two-permutation BAA proof;
even allowing arbitrary permutations does not escape that ceiling.

## Rejected control: schedule and cost

The local GRS outer needs no extension-field products. Its coefficients
are embedded GF256 scalars, so fixed multiplication is bytewise GFNI.
Use four interleaved message rows per ZMM and direct plane stores, as in
`HdFastOuter.h`. For [128,64], two 64-point cosets replace the four cosets
used by the existing [256,64] implementation. No systematic conversion is
required.

The pruned [128,64] circuit has 321 nontrivial GFNI sites per four rows.
At K20 it therefore uses 1,314,816 GFNI instructions. The two weighted scans
each use N/4 ZMM products. `wide16::multiply512` implements each product
with two GFNI instructions, giving 2,097,152 more GFNI instructions.
The total is 3,411,968, excluding prefix shuffles and XORs.

For comparison, the existing [256,64] outer alone uses 2,887,680 GFNI
instructions. The rejected control also halves the internal vector length
from 4K to 2K and removes the four-state Cauchy mixing circuits. These are
operation-count advantages, not measured speedups.

A four-record XOR prefix uses two masked 128-bit lane shifts and XORs.
XOR the incoming carry into all four records, then broadcast the last
record as the next carry. Compute several independent local prefixes
before carrying between vectors. Keep group sizes fixed and explicit.
There is no field-product dependency through the carry.

The first scan reads sequentially. The second scan should issue independent
gathers before consuming their values. Reuse the best validated 128-bit
record gather rather than an unexamined generic gather abstraction.
Prefetch indices must be bounded or explicitly padded with valid entries.

A simple complete schedule is:

1. Evaluate the outer into a 32MiB striped buffer.
2. Apply the first weighted scan in place, sequentially.
3. Gather through pi, apply the second diagonal, and prefix-scan into output.

Ignoring cache-line amplification, that schedule moves 176MiB of payload:
16MiB input reads, 32MiB outer writes, 64MiB for the in-place scan, and 64MiB
for the final read/write pass. Two expanded GF(2^16) coefficient arrays add
16MiB of metadata reads; the 32-bit permutation adds 8MiB. Setup and scratch
allocation stay outside encoding.

Do not extrapolate the existing regional gather time to this global gather.
Its 32MiB source exceeds L2 and competes for the active core's last-level
cache with metadata and output traffic. That stage could consume the
entire 5ms budget even when arithmetic is cheap. Exact source bucketing
would require additional passes and must be charged to the complete call.

Fusing the first scan into outer production requires care. Row-tiled outer
evaluation visits different columns out of chronological order. One can
compute column-local prefixes, then propagate column totals. That still
requires a correction pass, or a per-column correction during later
gathering. It is not a free elimination of the sequential scan.

### Completed native measurements

The full baseline median is 12.152ms over 101 calls. A separate diagnostic
run reports 3.272ms outer, 1.191ms prescan, and 7.746ms global route/postscan.
Diagnostic phase times need not sum to the separately measured full median.

Prefetch leads of 64 and 256 records did not help. Splitting the last route
and scan produced controls around 13.73ms and 14.10ms. Fusing prescan stores
into 256KiB buckets (`bucket_bits=14`), followed by local gathering, gave
12.075ms in a diagnostic run: a negligible change, not a promoted win.
The `bucket_bits=16` diagnostic took 13.926ms.

These results establish neither a certificate nor competitive encoding
cost. Further tuning of this scalar topology is not recommended because
the deterministic distance obstruction already rejects the construction.

## Field choice and alternative directions

GF256 remains an arithmetic cost control, not a way around the obstruction.
Its two scans need 1,048,576 GFNI instructions, for 2,363,392 in total with
the [128,64] outer. Merely storing GF256 operations in 128-bit records does
not enlarge the coefficient field or improve its distance.

For orientation, solving H_q(delta)=1/2 gives random-linear-code benchmark
distances of approximately 38.05%, 43.82%, 46.88%, and 49.22% for field sizes
2^8,2^16,2^32, and2^128. Here
H_q(delta)=delta*log_q(q-1)+H_2(delta)/log_2(q).
These are neither universal upper bounds nor guarantees for this topology.
GF(2^16) is the sensible first 40% arithmetic target; a 45% or near 50% claim
needs separate field and structural analysis.

Adding a reverse scan in the same order does not supply another independent
interleaver. Local interval witnesses can impose both directional boundary
states at constant cost. Such a change should first face a boundary-cut
analysis, not inherit a two-permutation BAA proof.

A cheap implementation control is a weighted two-dimensional prefix.
Arrange the once-permuted input as an H-by-W rectangle. Compute

    h[i,j] = d[i,j]*v[i,j] + h[i,j-1],
    y[i,j] = e[i,j]*h[i,j] + y[i-1,j],

with zero states at row and column boundaries. A row-major implementation
fuses both recurrences using W persistent records, sequential state-array
accesses, and two products per input. It needs no second random permutation
or materialized transpose. The vertical state occupies 16W bytes.

This does not fix the suffix obstruction. Both dependencies precede their
output in row-major order, so the complete inner remains lower triangular.
Used directly after a local outer and one uniform permutation, it is ruled
out by `HIGH_DISTANCE_GLOBAL.md`. Used after the striped precode, it still
inherits the band-to-suffix event above. Reversing one grid axis merely
changes a topological order; a uniform permutation absorbs that relabeling.

A forward/reverse cascade along the same stream has no such simple common
triangular order. It remains a separate, unproved possibility, but its
constant-cost interval boundary constraints are an early rejection test.
More directions are not evidence of another independent interleaver.
Independent diagonals are necessary for any claimed full-field cancellation
law; all-one characteristic-two scans do not manufacture a larger alphabet.

## Reusable sources and next step

| Source | Reusable component | Limitation |
|---|---|---|
| `spin/experiments/char2_field_spin/HdFastOuter.h:32` | Explicit 16-register GF256 FFT leaves and strided stores | Current public wrapper is [256,64], not [128,64] |
| `spin/experiments/char2_field_spin/WideDiagonal16.h:43` | Exact four-record, two-GFNI GF(2^16) product | Expanded setup costs 4 bytes per record |
| `research/workstreams/inner_design/redesign/accumulator_benchmark.cpp:68` | Fixed eight-way XOR scan with precomputed routing | Existing code is transpose scan/scatter and has multiple permutations |
| `research/workstreams/transposed_comparison/BinaryRaa.h:18` | Explicit scan/scatter and independent adjoint pattern | Binary arithmetic and online Feistel timing do not match this proposal |
| `../gen-BAA/field-universality/expander-lift/bench/baa_prime_bench.cpp:129` | Input-weighted scan: independent products before carries | Prime arithmetic and two uniform permutations; transfer scheduling only |
| `../gen-BAA/field-universality/expander-lift/bench/fp128_scan_schedule.h:71` | Bounded prefetch and fixed gather scheduling | Feedback-weighted prime recurrence; not the proposed map |

The companion `gen-BAA/rs-certification` certificates assume two independent
uniform permutations and q>=2^127. They cannot certify this one-permutation,
GF(2^16) construction. Likewise, its integer Newton outer cannot be reduced
modulo two as a substitute for the characteristic-two GRS outer.

Do not continue the scalar candidate's certificate or performance campaign.
The next proposed topology uses a GRS[64,32] outer, a true order-two precode,
one regional routing stage, and a final scalar scan. The precode composes
two independently weighted scans in the same order; it is not two unrelated
parallel state chains. Its inverse has bandwidth two and at most three
entries per column. With the postscan, the check-group degree bound is six.

The degree-six ceiling permits a 40% target but excludes 45% at this scale.
It supplies no distance lower bound or setup-failure guarantee. Before any
optimization campaign, define the actual regional route and coefficient
law, then test their structural cuts and finite-field specialization losses.
Only a complete bound over all nonzero messages can establish the requested
40-bit margin. This note makes no proof or speed claim for that new topology.
