# K18: exact certified s22 inner on the optimized small outer

This isolated experiment implements the retained K=262144 certificate's
exact inner maps using the newer outer kernel and routing schedule. It does
not change the library or the frozen K16 checkpoint.

## Measured result

The selected K18 kernel takes **about 0.467 ms** for precomputed transposed
encoding of 128-bit elements, with the existing **10% / 42.09-bit** certificate.
The library control takes about 0.606 ms: the candidate reduces latency by
23%, or increases throughput by about 1.30x.

| Configuration | Median time | K18 certificate status |
|---|---:|---|
| Library control, PacketRsT64S20 | 0.605622 ms | Unchanged by this work |
| Small outer with paired15 | 0.399557 ms | Full certificate remains open |
| Exact certified s22, initial implementation | 0.502955 ms | Retained certificate applies |
| Exact s22, fused feedback | 0.470008 ms | Same code as initial implementation |
| Exact s22, reordered feedback inputs (selected) | 0.466953 ms | Same code as initial implementation |

The certified path is about 17% slower than the lightweight paired15 path.
The larger state is therefore not free, but it closes the proof without
giving up the entire small-outer performance gain.

All runs were serialized on Peach's Ryzen 7950X, CPU15, GCC15.2, with a
configured 4.5 GHz maximum and huge-page preference. Each timed process in
the table used five warmups and 1001 calls, without diagnostic clocks. The
final variant comparison used seeds 41 and 113 in A/B/B/A order; figures
are medians of process medians. The library/paired15 comparison used the
same two seeds and reversed order in the preceding holdout. The initial
implementation used a one-seed A/B/B/A comparison against fused feedback.
One paired15 process measured 0.428 ms; the other three measured
0.398--0.400 ms. The table retains the four-process median rather than
discarding that outlier.

All three s22 implementations have identical output checksums for each seed
and call count. Final diagnostic clocks give about 259 microseconds for
inner/routing and 206 for the outer. These separately measured medians need
not sum to the median total.

Two implementation changes account for the improvement from 0.503 ms:

1. Keep the fixed feedback transform packed and combine it with the random
   update before unpacking. This removes redundant conversion and lane shuffles.
2. Reorder only the feedback transform's input coordinates so existing vectors
   supply them directly. Adjust its constant matrix columns accordingly.

The selected source has SHA-256
`a56e2f4be564b36a316ad01f411398e8d25f63590405fde2ee0289d60d696148`.
All 12 length/seed checks passed in both Release and AddressSanitizer plus
UndefinedBehaviorSanitizer builds. Each checks routed data against the
literal scalar oracle, final output, forward/transpose adjoint, four
alignments, padding and output guards, input preservation, and in-place use.
The exact map and distribution correspondence received an independent review.

Raw timing logs remain ignored under `output/k18-cert22-*.log`.

## Construction and certificate

The message has K=262144 bits and the output has N=524288 bits. The outer
has 2048 independent 128-to-256-bit groups. Each group contains four
systematic RS[16,8] words over GF(16). Each of its sixteen 16-bit symbols
receives an independent uniform nonzero GF(2^16) multiplier and a fixed coordinate
permutation. For each fixed nonzero input symbol this produces a uniform
nonzero 16-bit vector, independently across symbols. This is the joint
fixed-message distribution used by the retained first-moment proof; it does
not assert equality of the entire scalar-multiplier and GL16 ensembles.

Routing independently permutes each group's 64 four-bit packets and the
group slots in each of 64 regions. The inner processes 64 input bits per
step, with 22 state bits, an independent uniform invertible state update
per step, zero initial state, continuous state across regions, and no flush.
The implementation's seeded sampler follows the same ideal-randomness
convention as the retained certificate.

The retained outward-rounded certificate bounds the probability of any
nonzero message having output weight at most 52428 by less than
2^-42.09084006035199. Thus its good setups have minimum distance at least
52429 (at least 10%). The audit authenticates all 103 source pins, four
component receipts, all 2048 occupancy endpoints, and their exact dyadic
union. It is a recheck of an existing certificate, not a new numerical replay.

- Certified map SHA-256: `e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.
- Receipt SHA-256: `20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b`.
- [Certificate audit and basis warning](../../../research/workstreams/k18_state_ladder/README.md#the-certified-s22-fallback).
- [Original replay instructions](../../../research/workstreams/k16_design/README.md#power-of-two-checkpoint-k--218).

## Exact basis change

Let A be the certified 64-by-22 expansion matrix, with feedback C=A^T.
The fast kernel uses a monomial expansion matrix A_m and an invertible
22-by-22 coordinate matrix P satisfying A=A_m P. Writing its state as
b=P a changes the reverse step to

```text
output = input + A_m b
next_b = (P U P^-1) b + (P P^T) A_m^T input.
```

Here U is the literal reverse state-update matrix. The fixed factor P P^T
in the feedback is essential: merely replacing A by another basis for the
same expansion code and resetting C=A^T changes the code. The scalar oracle uses
the original certified coordinates and does not use the fast basis or GFNI
tables. Two extra physical SIMD slots are zero padding, not extra state bits.

## Reproduce

Run from the repository root. The generator's default mode verifies the
checked-in map; only `--write` regenerates it.

```text
python -B spin/experiments/k18_certified22/generate_maps.py
python -B research/workstreams/k18_state_ladder/audit_certified22.py tmp/rs16-s22-k262144-whole-fresh-p256.json
python -B spin/experiments/k16_codesign_100us/checkpoint/verify.py
```

The experimental AVX-512/GFNI build reuses the frozen outer and research
setup. The installed library supplies buffer ownership and CPU detection.

```sh
cmake -S spin/experiments/k18_certified22 -B out/k18-certified22 \
  -Dspin_DIR=/absolute/path/to/spin/package -DCMAKE_BUILD_TYPE=Release
cmake --build out/k18-certified22 -j3
ctest --test-dir out/k18-certified22 --output-on-failure -j1
flock /tmp/prindal-addition-encoder-benchmark.lock \
flock /tmp/bare-spin-benchmark.lock \
flock /tmp/hypercat-benchmark.lock \
flock /tmp/hypercat-global-benchmark.lock \
taskset -c15 out/k18-certified22/spin_k18_certified22 262144 41 1001 0 1
```

Arguments are K, seed, timed calls (zero runs checks only), diagnostic stage
clocks, and huge-page preference. Never run concurrent benchmarks. Timing
excludes setup, allocation, and checks. The page policy permits an ordinary
allocation fallback. Tests at other lengths validate the implementation,
not a distance certificate for those lengths.

The preserved K16 source check still matches all 220 pins. No production
default, K16 source, or certificate source was modified. The next step is
to freeze this K18 implementation/certificate pair as a research checkpoint
before considering a library dispatch change.
