# Optimized high-distance checkpoint

The current certified frontier is approximately **5.5 ms** for full forward
encoding at K=2^20, rate 1/2, and 40% GF(2^16)-symbol distance. This is about
2.4 times faster than the first 13.249-ms certified native implementation.
The sub-5-ms aspiration was not reached. Neither production nor the paper
has been changed.

The final build retains three fully certified choices:

| Native selection | Emission / feedback / refresh | Whole-code margin | Final reversed medians |
| --- | --- | ---: | ---: |
| `bucket4448-fused-nt` (45) | Global RS / period 4,4,4,8 / rank one | >56.6591 bits | 5.451 / 5.440 ms |
| `bucket8-rank2-nt` (41) | Four local RS maps / eight-check / rank two | >50.9966 bits | 5.457 / 5.460 ms |
| `bucket48-fused-nt` (42) | Global RS / alternating 4,8 / rank one | >51.0426 bits | 5.551 / 5.482 ms |

Select mode45 as the research checkpoint, keeping mode41 as its speed
control. Their ordering varies across runs; the difference is not a robust
speedup. Mode45 has the larger proved margin. Its larger margin is a result
of the particular bounds and witnesses, not a monotonicity claim about
using less feedback.

## Construction and guarantee

Let F=GF(2^16). A 128-bit application record contains eight parallel F lanes;
all lanes use the same scalar code. The message has K=1,048,576 records,
and the output has N=2,097,152 records. Distance counts nonzero F symbols,
not bits. The shared scalar-code guarantee also bounds nonzero records.

The outer consists of 8192 RS[256,128,129] rows. Independently sampled
nonzero F scalars multiply their output coordinates. Collect each outer
column into a region and independently permute its 8192 rows uniformly.
This is one routing stage; it is not a uniform permutation of all N positions.

Process the resulting stream in packets X of 64 symbols, with state Q in
F^16. Mode45 emits Y=X+A Q, where A is RS[64,16] evaluation. Its update is

    Q_next = Q + D2 (H D1 C X + ones (v^T Q)).

The feedback C uses four RS[16,4] transpose maps on three packets out of
four, then two RS[32,8] transpose maps on the fourth. H is the fixed Cauchy
interpolation map from nodes 0,...,15 to nodes 16,...,31. Each packet has
fresh independent nonzero diagonals D1,D2 in setup. Given D2, sample v
uniformly subject to v^T D2 ones != 1. The tuples are independent across
packets and independent of the outer and route. Setup is sampled once and
fixed for all messages. The initial state is zero; the terminal state is
unrestricted. See [the shared-refresh argument](HD_SHARED_DIAGONAL.md),
[local feedback bounds](HD_BLOCK_FEEDBACK.md), and
[period composition](HD_RETURN_DESIGN_SCREEN.md#using-eight-check-feedback-less-often).

For this setup distribution, the probability that **any** nonzero message
has output weight at most 838860 is below 2^-56.6591. Thus, except with that
setup-failure probability, minimum distance is at least 838861, which is
above 0.4N. This is a whole-code statement after setup is revealed, not a
per-message sampling claim.

The 119 pinned rational witnesses and the dense-tail bound cover all 8192
active-row counts. Independent replay at 384-bit directed precision gives

    ideal union log2 upper < -68.92219588903863,
    density loss bits    <  12.26300140707000,
    final log2 upper     < -56.65919448196863.

The displayed final bound comes from the directed computation, not from
subtracting the rounded numbers above. Mode41 has 133 corrected witnesses;
mode42 has 104. Both also pass separate 384-bit replays. The scripts retain
the exact rational witnesses; no numerical optimizer is needed for replay.

## Measurement and validation

Measurements use Peach's Ryzen 9 7950X, CPU2, the established 4.5-GHz
userspace setting, GCC 15.2, and AVX512/GFNI. Each sequence uses three
warmups followed by 501 calls per selection, then reverses the order.
Benchmarks run strictly sequentially under both host locks. No compilation
or proof workload runs on that host during measurement.

Full timing includes the outer, fused outer scaling, routed reads, complete
inner recurrence, output stores, and the final streaming-store fence.
Setup, allocation, input generation, and checksum calculation are excluded.
Input is restored outside every timed call. The deterministic seed is a
reproducible performance fixture, not an information-theoretic certificate
for that seed or its PRNG. The certificate above uses the ideal setup law.

All 45 native variants pass Release and ASAN/UBSAN checks. Checks include
every packet input/state basis direction, independent field arithmetic,
nonzero state, routing, alignment, guards, and exact aliasing. Full K20
checks additionally compare every output record of modes41,42,45 against
the dense scalar inner, including unaligned output and in-place calls.
The final help-text-only rebuild passes the complete Release suite again.
The full high-distance Python regression suite passes all 446 tests.

For seed113 and the benchmark input, the output checksums are:

| Mode | Checksum |
| --- | ---: |
| 45 | 7582624837270167152 |
| 41 | 7812899659201065080 |
| 42 | 16994121132285933612 |

## Why stop this optimization round?

The [chronological ledger](HD_ENCODING_OPTIMIZATION.md) records the complete
search, including proof defects found and corrected. This round exhausted
the concrete candidates with a plausible near-term path to a faster
certified encoder. It does not prove that this family is optimal.

| Direction | Outcome |
| --- | --- |
| Smaller state, fewer refreshes, rank-one feedback | Current bounds fail sparse or middle classes; no replacement certificate |
| Wider packets and states20/24/32 | Validated native implementations are slower; the state20 proof also remains incomplete |
| Smaller RS outer | Current middle-class bounds fail; smaller arithmetic alone does not preserve the certificate |
| Block feedback, shared diagonals, cheaper emission | Yielded the three fully certified alternatives above |
| Carried/common/scattered refresh | About 5.3 ms in the tested cached implementations, but strengthened return bounds still fail middle classes |
| Layout, gather, offsets, prefetch, shared loads | Kept the winning schedule; later controls have no repeatable gain or regress |
| Streaming output and register lifetime | Retained; an observed store-order regression was repaired and checked against the frozen binary |
| Compact setup and non-alias hints | Correct but no repeatable improvement; disabled |
| Alternative dyadic Cauchy factorization | Exact, but the tested schedule needs 32 fixed products versus 20 in the current FFT pair |

Phase controls put the outer near 2.4 ms and the routed inner near 3.1 ms;
pre-routing only lowers the latter to about 2.85 ms. These are diagnostic
medians for the earlier cached rank-two implementation, not additive
components of the final timing. Routing alone is no longer the dominant
obstruction to another large gain.

An F-linear extension-field multiplier remains a possible new feedback
design, not an untested drop-in optimization. Its basis conversion and
product circuit would have to beat the now much cheaper shared refresh.
Arbitrary binary mixing would change the F-linear interface and lose the
existing projective-counting argument. No speed or certificate is claimed
for either option.

Further work should begin with a materially cheaper recurrence and a
credible bound on correlated returns, or a new outer/inner allocation.
More numerical tuning of the failed bounds is not presently a convincing
route below 5 ms. The recommended next step is to review and promote this
checkpoint before starting another structural search.

## Reproduce and preserve

From the repository root, on the stated ISA/compiler:

```sh
cmake -S spin/experiments/char2_field_spin -B out/char2-hd-final \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15 \
  -DSPIN_CHAR2_SCALAR_LOOKUP=OFF -DSPIN_CHAR2_CLMUL_REDUCTION=OFF
cmake --build out/char2-hd-final --target spin_char2_hd_sumrank_forward -j4
out/char2-hd-final/spin_char2_hd_sumrank_forward --check
out/char2-hd-final/spin_char2_hd_sumrank_forward --large-check
flock -n /tmp/spin-char2-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  taskset -c 2 out/char2-hd-final/spin_char2_hd_sumrank_forward \
  --bench 501 full 20 bucket4448-fused-nt
```

Leave optional `HD_SUMRANK_*` macros undefined. Those macros retain negative
controls; they are not enabled by the final measured build. Native source
and headers remain in `spin/experiments/char2_field_spin`.

With Python, python-flint (tested at 0.9.0), NumPy, and SciPy available,
run the independent proof replays from this workstream's `theory` directory:

```sh
python hd_sparse_feedback_witnesses.py --precision 384 --workers 4
python hd_rank_two_block_witnesses.py --precision 384
python hd_bucket_feedback_witnesses.py --precision 384
python -m unittest discover -p 'hd_*_tests.py'
```

Proof workers only parallelize independent arithmetic witnesses. Do not
run proof arithmetic concurrently with a benchmark on the same host.

The final native executable is preserved on Peach as
`/tmp/hd-blockwise-20261006/optimized-checkpoint-69329fef`, SHA-256:

```text
69329fefb8f3e4020488411adfcaa8c0aa694eaf43fdc8b29aaa9e88c57e46c3
```

It contains all three certified choices. Older frozen executables and
witnesses are preserved; nothing has been promoted, committed, or pushed
as part of this optimization round. The source snapshot is local and
ignored by Git; binary hashes identify measurements, not portable builds.

[The source snapshot](../../.proof-archives/gf65536-high-distance-k20-20261007-source.zip)
contains 774 source and documentation files, including the complete native
experimental project, proof scripts, tests, and pinned rational witnesses.
Every archived file was read back and checked against its SHA-256 manifest.
Build products, raw numerical dumps, and measurement logs are excluded.
The archive is 2,554,944 bytes; its SHA-256 is

```text
3690641a7404f5865f1ae745cb508b6698a010c613555008294a9734f69e8101
```

The archive predates this hash paragraph and the final regression-count
update, avoiding a self-referential archive hash. Its native implementation
and proof inputs are the validated versions. Principal native source hashes:

```text
bd6eac925ac6f9b8504ee76708a2be9657bd0821a9ac5b69db2036677a02fc5e  hd_sumrank_forward.cpp
12a91385f98fa6012c53d2fb208dc9b9b413a7b372886a6191fe90072ddb7e7c  HdSumrankInner.h
dc36c46d95e0170f3060eaee11f2e2cbb1acacd90eead651da9141e336496083  HdSumrankOuter.h
```
