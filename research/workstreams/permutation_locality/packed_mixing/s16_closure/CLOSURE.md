# Packed GL32 with a 16-Bit State

2026-09-30. **Closed: greater than 10% relative distance with more than
55.08 bits of whole-code setup-failure margin**, at K=2^20 and rate
one-half. The construction uses 64-bit inner steps and independent uniform
invertible 16-by-16 state maps. Fresh 384-bit replay covers every occupancy
and passes the exact aggregate and unchanged-source checks.

The checked precomputed transposed encoder takes **5.3445045 ms** on
128-bit XOR elements. Its matched 19-bit-state control takes 5.3241165 ms.
The 0.38% difference does not demonstrate a speedup. The
[19-bit-state certificate](../R4_CLOSURE.md), with more than 52.05 bits of
margin at greater than 10% distance, remains preserved.

## Construction

The forward encoder maps K=2^20 bits to N=2^21 bits. Its fixed outer is
BCH[256,128]. Setup samples the following independent choices once;
every message then uses the same realized encoder.

1. Encode 8,192 BCH rows and collect every four adjacent rows into one
   group. There are 2,048 groups.
2. Divide each group into 32 canonical blocks of eight adjacent columns.
   Apply an independent uniform GL(32,2) map to each four-row/eight-column
   block.
3. Sample a uniform permutation of the 256 columns for each group.
   The four rows share that permutation, so each column forms a four-bit
   packet. Route column j to region j.
4. Independently shuffle the 2,048 packet positions in each of the 256
   regions. Concatenate the regions and retain the order within each packet.
5. Process the resulting bits in 64-bit steps, with a persistent 16-bit
   state. There is no additional GF16 packet randomizer.

For input X_i in F_2^64 and entering state Q_i in F_2^16, the inner is

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i,
    Q_0 = 0.

Each M_i is independently uniform in GL(16,2), independent of all routing
and local block maps. The output uses the entering state; feedback uses
the raw input and is added after M_i. State persists across region
boundaries, and the final state is discarded without a flush.

The fixed expansion A and feedback C=A^T are `Map64S16` in
[SelectedMaps.h](../../../../../spin/src/kernels/generated/SelectedMaps.h).
[kernel_t64.py](kernel_t64.py) reconstructs both maps from the declared
[matrix data](../../../rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json).
It checks their rank, enumerates every expansion word, and verifies all
16 four-column packet ranks. Both maps have rank 16; CA=0. The expansion
has minimum weight 16, and ker(C) has minimum weight 6. Neither fact alone
is the distance claim for the complete encoder.

Every setup realization has rate exactly one-half. The BCH outer is
injective, the local maps and routing are invertible, and the inner is
invertible successively from its output and known initial state.

## Claim and Proof Coverage

The bad-weight cutoff is 209715=floor(N/10). Let Z count nonzero messages
whose output has weight at most this cutoff. A bound E[Z]<=U implies
Pr[d_min<=209715]<=U over the stated setup distribution. Outside that
event, minimum weight is at least 209716, strictly more than 10% of N.
The margin is -log2(U): it bounds setup failure, not an attack work factor.

The proof separates messages by q, the number of nonzero four-row groups.
Sparse bounds cover q=1,...,32 and every possible support size for each q.
The dense comparison covers q=33,...,2048 through an exhaustive
106-interval partition of [2079/7829504,1], with disjoint interiors.
No interpolation between selected points is used.

Fix a message and an output-weight tilt z in (0,1]. Two successive
physical steps form a 128-bit proof block; the midpoint state is retained.
Let L_a(z) be a nonnegative envelope for one physical step's state-mass
transition, weighted by z^wt(Y), with a active four-bit packets.
A valid envelope for j active packets in the proof block is

    K_j(z) = sum_a binom(16,a) binom(16,j-a) / binom(32,j) * L_a(z) L_(j-a)(z).

The product is in chronological order. Conditioning on the outer support
data and packet counts in every proof block leaves uniform local subsets
under the regional shuffle. This conditioning does not fix the sampled
GL32 maps or packet labels. The GL32 maps supply independent nonzero
packet labels conditional on their support. Each physical-step bound
remains valid for the weighted entering-state distribution, so composition
does not assume a new independent state at the midpoint.

Fresh BCH shortening and incidence bounds control expected cumulative
counts by support size, averaging over the GL32 maps. The inner bounds
use the actual selected maps,
including the distribution of states reached from zero. The numerical
replay reconstructs these premises and reevaluates every rational witness
with outward-rounded arithmetic. Saved search scores are not proof inputs.
The coordinator sums all fresh dyadic endpoints exactly, rounds the sum
upward once, and checks the strict whole-code target again.

Both fresh child replays and the whole-code coordinator passed. Independent
exact-arithmetic audits of all 32 sparse and 106 dense endpoints confirm
the aggregate and the 521-file unchanged-source manifest. A separate
review checked the construction, conditioning, and physical-step composition.

| Covered messages | Margin from the outward upper bound (bits) |
|---|---:|
| Sparse, q=1,...,32 | 55.096181354755 |
| Dense, q=33,...,2048 | 62.938762900079 |
| Whole code | 55.089909760705 |

The whole-code upper bound is approximately 2.60786258523e-17. Its exact
upward-rounded dyadic representation is

```text
U = 79502844188172024009604177020280257326396612675787641964254393039785161112775352643955434064789631003832483596073478801761825 * 2^-470.
```

Integer arithmetic verifies U^100 < 2^-5508, giving the conservative
**greater than 55.08-bit** margin without relying on a rounded logarithm.
The sparse contribution now dominates; selected dense points no longer
limit this certificate. The result applies to this K, rate, distance
threshold, and setup distribution, not to unverified parameter ranges.

## Implementation and Matched Performance

[GfniT64.h](../../GfniT64.h) keeps the state packed for GFNI throughout the
inner loop. The 16-bit state fits in 256 bytes for 128-bit payloads.
The matrix-update core uses eight GFNI instructions; the precomputed
update tables occupy 1 MiB across 32,768 physical steps. Fixed-width
expansion and feedback circuits avoid allocations in the hot path.

The [complete encoder probe](../../gfni_t64_probe.cpp) checks explicit
matrix bases, each sampled state map and its transpose, independent inner
references, direct BCH encoding, scalar GL32 mixing, whole-encoder adjoints,
boundary inputs, and preservation of the unused buffer suffix. Checks ran
at K=2^14 with setup seeds 1 and 17 and at K=2^20 with seed 1.
[implementation_t64.py](implementation_t64.py) binds the explicit matrices,
generated source, archived binary, check logs, and timing records.

If E:F_2^K -> F_2^N is the forward encoder, the benchmark applies E^T to
N 128-bit elements and writes K elements. The matched campaign ran
serially on Peach CPU15, alternating candidate
and control processes. Each candidate had four processes of 101 measured
calls. Setup and allocation are excluded.

| Complete in-place transposed encoder | Median of process medians | Range of process medians |
|---|---:|---:|
| t=64, s=16, uniform GL16 | 5.3445045 ms | 5.338589--5.359809 ms |
| Retained t=128, s=19, four transvections | 5.3241165 ms | 5.307380--5.349489 ms |

The benchmark uses deterministic setup generation; the mathematical claim
uses the independent ideal setup distribution specified above. The
implementation checks bind the maps and recurrence, not a formal proof
of the complete software stack. No production default or paper claim
changes in this research result.

## Replay and Retained Evidence

Run from the repository root, choosing a new output directory:

```text
python -B research/workstreams/permutation_locality/packed_mixing/s16_closure/whole_t64.py --dense tmp/s16-closure/t64-dense-merged-p256.json --sparse tmp/packed-hill/sparse-t64-selected-screen-extend-p256.json --output-dir tmp/s16-closure/t64-whole-rerun-p384 --precision 384 --target-bits 40 --dense-workers 4
```

The sparse replay runs first. Four fresh processes then reconstruct their
own models and replay the dense witnesses. Before and after manifests
reject changes to proof sources or inputs during verification. All
numerical receipts remain in ignored `tmp/`; raw experimental data is
not part of the proposed code changes.

Search witnesses:

| Input | SHA256 |
|---|---|
| `tmp/s16-closure/t64-dense-merged-p256.json` | `55ebfe88c86b908240100801fdde8eaf4393b371daf3809f09ce9e7a272e9589` |
| `tmp/packed-hill/sparse-t64-selected-screen-extend-p256.json` | `86dbbfc72505869e76f00bd1df4036e3b5a137a232131d92b6b7c25b54dbc625` |

The final receipt and its fresh children are in
`tmp/s16-closure/t64-whole-p384/`:

| Receipt | SHA256 |
|---|---|
| `sparse.json` | `11930fd644ecb3d309d29af848318847eaa2974d485870ed0b5e0c104ce5b8f2` |
| `dense.json` | `ae61a2ca26424b7cb56b3b718c81d80e4123cfabe6ff356d3a0657b2cc3e6d20` |
| `whole.json` | `698e849a142641f6af3885f72a4ae5b2f285f0100a17051640c6a97bc0409c3f` |

The implementation archive is
`tmp/packed-hill/gfni-t64-performance/records.tar`; its checked source,
map, binary, and timing bindings are in `binding-final.json` alongside it.
That binding receipt has SHA256
`c3e288c88b8ac4915c3870764f074990f74b54ac8e8152adbbb1bf5970c41898`.
All 129 tests in this closure directory pass, including exact small-state
kernel checks, scope/coverage rejection tests, and archived implementation
binding. The tests complement rather than replace the full numerical replay.

The retained S19 whole receipt still has SHA256
`637a1ac049a7aa81441705e5a4fbeebf79e78fa80e63fe54f93cfe6dad05a8f0`.
All 479 source files recorded by that receipt remain byte-identical;
its implementation artifacts are also preserved.

The [progress ledger](PROGRESS.md) retains the unsuccessful t=128, s=16
screens and explains the move to shorter physical steps. A state map
cannot activate zero when CX_i=0. Reducing the step length changes those
activation events; merely strengthening the state map does not. The
shorter-step proof is independent of the earlier t=128 certificates.
