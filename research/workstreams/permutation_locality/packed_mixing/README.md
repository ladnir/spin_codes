# Mixing Inside the Packed BCH Kernel

This research track changes the encoder distribution. It does not change
production defaults or inherit the distance certificates of an earlier route.
The retained constructions and their proofs remain in the
[proof index](../PROOF_INDEX.md).

**Current result (2026-09-30):** fresh 384-bit whole-code replay certifies
**>10% relative distance with >52.05 bits of setup-failure margin** for
canonical GL32 with R4 at K=2^20, rate one-half. The latest persistent-state
GFNI implementation gives **5.334629 ms** precomputed transposed encoding
with 128-bit XOR elements (about 5.32--5.35 ms in matched confirmation).
See the [R4 closure record](R4_CLOSURE.md) for the construction, exact
claim, receipt hashes, and replay command, and the
[matched performance report](../gfni_inner_REPORT.md) for measurements.

The separate [16-bit-state construction](s16_closure/CLOSURE.md) is now
closed at **>10% distance / >55.08 bits of whole-code margin**, using
64-bit physical steps and independent uniform GL16 updates. Fresh
384-bit replay covers all occupancies, with aggregate margin
55.089909760705 bits. Its checked implementation measures 5.344505 ms,
versus 5.324117 ms for the matched R4 control. This is essentially equal
performance; the new result does not replace the retained R4 proof.

The [hill climb](HILL_CLIMB.md) combines stronger canonical BCH block
counts with four actual inner updates. Its sparse q=1..32 and exhaustive
dense q=33..2048 contributions sum to 52.052884355478 bits of margin.
Fusing the four updates preserves the sampled map exactly and saves 8.01%
against its same-build sequential R4 control. The later packed-stage
optimization saves another 4.91% in its own matched comparison and halves
the GL32 coefficient table to 8 MiB.

The [retained R2 certificate](FIRST_CLOSURE.md) remains **>9.5% / >36.68
bits**, with its earlier **6.203224 ms** measurement. R2's stronger sparse
10% bound is still not a whole-code result. No old proof, production
default, or paper claim has been replaced.

## Question and Scope

Can the local bit transpose already used by the GFNI BCH kernel also make a
stronger randomized outer stage affordable? The target workload is the
precomputed transposed encoder at K=2^20, rate 1/2, with 128-bit XOR elements,
four-row routing packets, and IMT(128,19), initially with two updates.
The separate R4 construction uses four updates and preserves the other
parameters; [its record](R4_CLOSURE.md) specifies the complete setup.

The packing is an implementation layout, not field arithmetic on the 128-bit
elements. In each packed byte, eight code coordinates belong to the same
payload bit. Four 128-bit SIMD lanes hold four BCH rows. Every proposed mixer
must preserve all 128 independent payload bit positions.

The new forward construction has the order

    four BCH rows -> local coordinate mixer -> shared column shuffle
                 -> regional packet shuffles -> IMT

Row-local mixers additionally retain the existing independent GF16 packet
randomizers before IMT. A final independent uniform GL(32,2) map on each
four-row/eight-column block supplies the required packet-label distribution
itself, so that candidate omits those randomizers. The transpose implementation
reverses the order and applies the transposed binary maps.

## Candidate Ladder

| Candidate | Random map per local block | Additional packed arithmetic | Proof issue |
|---|---|---|---|
| Dense-half shared GL8 | One map on eight coordinates, reused across four rows; dense half only | One GFNI per prepared vector | Preserves repeated rows and leaves the other half unchanged |
| Full shared GL8 | Same, on all 256 coordinates | One GFNI per vector; also pack the nearly systematic half | Output support depends on local binary rank |
| Full independent-row GL8 | Independent map for each row/octet | One GFNI with lane-specific matrices per vector | Breaks equal rows, but preserves zero row octets |
| Full GL32 | Independent map on four rows × eight coordinates | Four GFNI operations, three lane shuffles, three XORs per vector | Removes local rank dependence, but preserves empty blocks |

These are operation counts, not latency predictions. Full-word variants also
prepare the nearly systematic half of the BCH input and account for its
exceptional parity coordinate. The existing fixed BCH coefficient table is
retained; sampled mixers are factored rather than folded into a different
dense BCH matrix for every group.

For the distribution of one fixed message, a uniform nonzero GF(256) scalar
has the same effect on an active row octet as uniform GL(8,2). Likewise, a
uniform nonzero GF(2^32) scalar has the same effect on a nonempty 32-bit block
as uniform GL(32,2). This first-moment equivalence does not claim equal joint
distributions across different messages. Field multiplication can reduce setup
storage; an online speedup would require a separately tested kernel.

## Proof Boundary

A local mixer does not by itself randomize which of the 32 canonical blocks
are occupied. The shared column shuffle comes **after** mixing in the forward
construction. It randomizes the resulting packet support, not the input
partition seen by the mixer. A kernel for a fresh uniform input partition
therefore describes a different construction and cannot be silently reused.

The full-block map is nevertheless particularly simple. For any fixed nonzero
32-bit input, its output is uniform over the 2^32-1 nonzero vectors. Viewed as
eight four-bit packets, its support is Binomial(8,15/16) conditioned positive.
Conditional on its complete support, the active packet values are independent
uniform nonzero GF16 elements. Independent maps on different blocks preserve
this conditional independence.

Thus the missing outer statistic is the number H of occupied canonical
eight-column blocks, rather than the binary rank inside each block. If d(u)
uniformly bounds the dimension of the BCH code shortened to any u coordinates,
then a direct upper bound on the number of nonzero ordered four-row tuples
with H <= h is

    binom(32,h) * (2^(4*d(8h)) - 1).

Every tuple counted is supported in at least one union of h canonical blocks;
overcounting those unions is safe. This can be intersected with the existing
union-support bounds and refined by tuple rank. It requires no extra random
prepartition or coordinate gather. Local distribution identities and count
bounds are not a complete SPIN distance certificate.

The [local proof note](LOCAL_KERNELS.md) gives the distributions, exact
transport, authenticated counting diagnostics, and independence requirements.
The exact tests include exhaustive small matrix groups, packet-label
factorization, toy-code rank counts, and CDF transport checks:

```text
python -B -m unittest discover -s research/workstreams/permutation_locality/packed_mixing -p "test_*.py" -v
```

The old checker restricts supports to [38,256], whereas the new mixed support
can fall below 38. The new [sparse driver](sparse.py) and
[support-cover adapter](support_cover.py) instead authenticate the full
[5,256] domain and retain rational expected counts. The legacy checker is
unchanged. The new [dense point driver](dense_probe.py) handles all support
values, but selected mean probes do not establish a complete dense cover.
The [resumable dense driver](dense_cover.py) builds exhaustive interval
partitions, using freshly checked pointwise expected-shell bounds and
positive comparison mixtures. Its output still excludes the sparse prefix.
The [progress ledger](PROGRESS.md) records covered scopes, replay status,
remaining work, and the aggregate acceptance criterion.

Search and verification have separate entry points. The dense search can
split pending intervals to the width required by its stronger bound before
optimizing them (`--regional-presplit`). It can also reuse point witnesses
as parameter proposals (`--point-seeds`) and retain one freshly computed
inner transition polynomial in memory (`--region-cache`). None of these
options transfers a saved numerical bound: each accepted interval gets a
new outward evaluation of its counts and geometry. Independent replay
rejects these search options and uses the uncached evaluator.
For independent replay, `--replay-workers 3` optionally assigns the cells
to three fresh processes. Each builds its own ordinary inner model; the
parent checks that exactly the requested cells and witnesses were replayed.
The default remains sequential. A 384-bit, 22-interval regression matched
the sequential dyadic endpoints and aggregate exactly; see the progress
ledger for the retained receipts.

With point seeds, `--seed-width 1/256` first combines adjacent unresolved
cells and then subdivides them to the requested width. It does not combine
accepted cells or increase certified coverage. The regional bound is valid
on these wider intervals; the 1/1024 threshold limits ordinary proposal
search, not the mathematics. If the supplied hints fail on a wider cell,
the driver splits it without invoking the expensive fallback optimizer.
Normal optimization resumes at width 1/1024. A neighboring hint whose inner
operator is already cached goes directly to a fresh exact evaluation,
avoiding a redundant floating screening pass.

The isolated [count-cap prototype](count_caps_fast.py) further reduces
search arithmetic by selecting one independently valid MGF bound per count.
Its floating selector never supplies a numerical endpoint; selected bounds
are recomputed with outward arithmetic. It is tested but not yet connected
to either search or replay. The progress ledger records its limited
actual-witness comparison and checker-component timing.

For a parallel continuation, [the subtree driver](cell_search.py) assigns
disjoint pending roots from one immutable snapshot to fresh search processes.
Each process authenticates the same outer premises and keeps the global model
domain. It need not replay accepted cells outside its assignment. [The
merger](cell_merge.py) checks that every original pending root is assigned
exactly once and that each returned subtree has no gaps or overlaps. Its
output is a search cover, not a certificate; worker numerical endpoints
are discarded, and final replay remains mandatory.

Once the dense partition is complete, [the whole-code replay coordinator](whole_replay.py)
runs the sparse and dense replays in separate fresh processes, checks their
construction, cutoff, provenance, and complementary scopes, and adds all
fresh dyadic upper endpoints. It writes a whole-code certificate only when
that sum is below the requested failure-probability bound. An incomplete
dense partition is rejected before either replay starts. The coordinator
has now completed the [first whole-code closure](FIRST_CLOSURE.md).
Its optional `--dense-workers` setting enables the same parallel evaluator
for the dense stage without overlapping the sparse and dense stages.

### First Inner Diagnostic

At the comparison mean 0.032, for occupancies at least 33, the initial
two-update / 10% point bound failed. Refining the shortened-dimension bounds
reduced its outward log2 upper bound from +12966.65 to **+964.38**, but did not
close that point at 10%. A positive log2 upper bound is uninformative; it
does not exhibit a low-distance codeword.

The refined witness uses output tilt lambda=227094387/2500000000. Its bound
depends on the integer output cutoff T only through exp(lambda*T). Keeping
the same point and witness, the bound at a smaller cutoff T' is therefore
U(T)*exp(lambda*(T'-T)). Evaluating that expression outward at 256 bits gives:

| Requested relative distance | Integer cutoff | Log2 upper bound at this point |
|---|---:|---:|
| 10% | 209715 | +964.38; fails |
| 9.7% | 203423 | +139.81; fails |
| 9.6% | 201326 | -135.01 |
| 9.5% | 199229 | -409.82 |

The lower-distance entries are exact-threshold rescalings of one checked
witness, not newly optimized searches or whole-code margins. They motivated
the original 9.5% closure target, which is now complete. The 9.6% sparse
results and subsequent R2 and R4 investigations are recorded separately
in the progress and hill-climb ledgers.

Reproduce the refined point with:

```text
python -B research/workstreams/permutation_locality/packed_mixing/dense_probe.py --distance .1 --means .032 --last-lp 104 --refined --precision 256 --output tmp/packed-mixing-fresh-point.json
```

The retained local receipt is
`tmp/packed-mixer-results/canonical-r2-d10-m032-refined-p256.json`.

## Implementation and Validation

Research-only sources are [the generator](../packed_mixer_codegen.py),
[the C++ comparison](../packed_mixer_perf.cpp), and
[the serial runner](../packed_mixer_run.sh). The runner uses the existing
reference build and benchmark locks. Coefficient generation and allocation
are outside the precomputed encoding interval. Every timing comparison must
use a matched control from the same build and harness.

## Retained R2 Optimized-Driver Result (2026-09-30)

The [minimal-driver investigation](../packed_driver_REPORT.md) identified
outlined routing callbacks inside the inner loop in the large comparison
build. An explicitly inlined, fixed-width emitter removes that overhead.
Two seeds, both process orderings, and 101 measured calls per process give:

| Full precomputed transpose at K=2^20 | Median | Range of process medians |
|---|---:|---:|
| Retained optimized legacy binary | 5.829 ms | 5.820–5.837 ms |
| Force-inlined legacy, matched new build | 5.805 ms | 5.802–5.814 ms |
| **Force-inlined full GL32** | **6.203 ms** | **6.176–6.229 ms** |

The new GL32 construction costs **0.398 ms (+6.86%)** against the optimized
same-build legacy control. It is not a speedup over the old optimized code;
the purpose is a stronger distance guarantee near the 6 ms target. Independent
full-encoder reference and adjoint checks pass at K=2^14 and K=2^20 for both
seeds. The following earlier measurements remain historical, not the current
baseline. No production kernel or old certificate was replaced.

## Earlier Large-Driver Results

The final matched build includes an exact-code control that moves the old
GF16 packet maps into the packed stage. Reindexing the same maps through the
route preserves the complete encoder; it is not a new distribution. Its
diagonal packed submatrices can use AND masks instead of general GFNI maps.

| Final build, full precomputed transpose at K=2^20 | Median | Range of process medians |
|---|---:|---:|
| Unchanged control | 6.582 ms | 6.554–6.641 ms |
| Same GF16 maps moved into packed GFNI | 6.441 ms | 6.362–6.496 ms |
| Same GF16 maps moved into packed AND masks | 6.438 ms | 6.431–6.451 ms |
| New full GL32 mixing, without GF16 | 6.418 ms | 6.379–6.445 ms |

The exact-code AND variant improves the matched control by **2.19%** and
passes full-buffer equality with the original encoder for both seeds at
K=2^14 and K=2^20. It preserves the existing construction and certificate.
At that earlier checkpoint, the new GL32 candidate was still uncertified.
The later R2 and R4 proof records above supersede that status, not the
historical measurements.

### Earlier Build and Timing Sensitivity

The strongest local randomizer is the best measured candidate because it
replaces, rather than merely supplements, the existing packet randomizers.
The following are medians of four process medians: two seeds, both execution
orders, 101 measured calls per process, with setup excluded.

| Full precomputed transpose, K=2^20 | Median | Range of process medians |
|---|---:|---:|
| Matched unchanged GFNI control | 6.573 ms | 6.544–6.600 ms |
| Full shared GL8, retaining GF16 | 6.856 ms | 6.848–6.963 ms |
| Full independent-row GL8, retaining GF16 | 6.978 ms | 6.966–7.014 ms |
| **Full GL32, removing redundant GF16** | **6.094 ms** | **6.055–6.151 ms** |

In that first build, GL32 saved 0.479 ms (7.29%) against the matched control.
Adding the exact-code control variants and rebuilding changed the GL32
inner/routing time; its packed mixer/BCH phase stayed almost unchanged.
The later minimal-driver investigation above identified routing-callback
outlining as a substantial cause. These earlier **6.09–6.42 ms** measurements
must not be presented as a robust 7% improvement over optimized legacy.
A shorter
31-call screen measured 7.687 ms when it unnecessarily retained GF16. Thus
the useful optimization is the combined distribution change, not the claim
that four-GFNI mixing is free.

Historical controls in a different driver measured about 5.83–5.86 ms.
The earlier [BCH comparison](../bch_compare_REPORT.md) records the original
inner/routing timing discrepancy. The new minimal driver restores that
legacy baseline. No synthetic timing correction is applied. See the
[benchmark report](../packed_mixer_REPORT.md) for scope, checks, and commands.

The full GL32 prototype stores 16 MiB of expanded coefficient data at this
size. All candidate kernels passed comparison against scalar mixing plus
the existing BCH implementation on all 131,072 physical input basis bits
for two setup seeds. Full-encoder checks verify the mixer/BCH result and
preserved suffix, but share the already tested inner/routing implementation.
These checks establish implementation agreement, not a distance certificate.

## Further Directions

The same four-term kernel can also apply a block-diagonal matrix: two
independent GL16 maps on four rows × four columns, four GL8 maps on four
rows × two columns, or eight GL4 maps on individual four-row packets.
This changes the column-group width without changing the general GFNI
instruction sequence. These are different from the **row-pair** GL16
distribution in the local note, which keeps eight columns together.

At width one, the exact GF16 relocation control above preserves the encoder,
not just its first-moment distribution. For the new distributions, a fresh
width comparison using one authenticated dimension table favors width eight:
at support 115 its log2 expected-count cap is 252.31, versus 288.98 for width
four and 344.17 for width two. These are count bounds, not security margins.
See [the width comparison](width_compare.py) and the local proof note.

If empty blocks remain the obstruction, add communication between packed
blocks rather than repeating isolated local maps. Candidates include
row-wise GF256 shear pairs and GF(2^32) MDS pair mixers. They need their own
cancellation model and measurements. Repeated mixing confined to the same
local block cannot activate an empty block.

A single random prepartition shared by all groups could retain one BCH
coefficient table, but introduces common randomness. A proof must condition
on and certify that fixed partition, or analyze the resulting dependencies;
raising an average local enumerator to the number of groups is not valid.
