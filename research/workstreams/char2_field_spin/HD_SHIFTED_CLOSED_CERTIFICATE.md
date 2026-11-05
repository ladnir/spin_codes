# Certified paired routing at 35% distance

The independently swapped paired route now has a complete distance
certificate. At K=2^20 and rate 1/2, its setup-failure margin is
55.092738508223 bits for 35% symbol distance. Matched native measurements
give 3.798 ms, versus 4.049 ms for the retained certified reference.
This is a research checkpoint; no production profile or paper is changed.

On 2026-10-09, the user selected this checkpoint as a leading research
design at 35% distance. The faster-design search at this distance was
parked; the next direction considers larger alphabets closer to 50%.
A new 384-bit directed aggregation
of the retained, pinned bounds confirmed complete coverage and the same
55.092738508223-bit ceiling. This was not a fresh coefficient replay or a
new performance measurement. See the [wrap-up](HD_34MS_CLOSURE.md#wrap-up-2026-10-09)
for the preserved alternatives and their remaining gaps.

## Construction and guarantee

Work over GF(256), with K=1048576 and N=2097152. Partition the input into
8192 groups of 128 symbols. Each group uses four RS[64,32] constituents.
Independently in each group and output column, uniformly sample an
invertible 4-by-4 Toeplitz matrix and apply it to the four constituent symbols.

Collect these outputs into 256 regions. Regions alternate the even and odd
adjacent-pair matchings on the cyclic group array. Independently in each
region, uniformly permute its 4096 pairs and independently swap each pair
with probability one half. All swap bits are independent of the permutations,
outer mixers, and refresh diagonals. Setup is sampled once and fixed for
all encodings.

The inner processes 64 symbols per step with a 16-symbol state. Its fixed
expansion map A consists of four RS[16,4] maps; its feedback map is C=A^T,
with CA=0. Starting from Q=0, each step computes

    Y = X + A Q,
    Q_next = D2 H16 D1 (Q + C X).

Here H16 is the fixed GF(256) MDS interpolation matrix from evaluation
nodes 0 through 15 to nodes 16 through 31. Each diagonal entry of D1 and
D2 is uniformly sampled from GF(256) minus zero. These samples are
independent across both maps and every step, and independent of the outer
and routing setup.
The state persists across all regions; it is not reset at region boundaries.
The implementation retains a striped SIMD layout and uses a fixed output
permutation, which does not change Hamming weight. The
[construction analysis](HD_SHIFTED_PAIR_SWAPS.md) specifies the operators.

Let S denote this ideal independent setup and G_S its scalar encoder.
The complete directed calculation gives

    Pr_S[there exists nonzero m in GF(256)^K with wt(G_S m) <= 734003]
        <= U,
    log2 U <= -55092738508223 / 1000000000000.

Thus, except with probability at most U, every nonzero message has output
weight at least 734004, which exceeds 0.35 N. The statement is uniform over
all messages, not a fixed-message guarantee.

The native implementation applies the same scalar map to all sixteen byte
components of each 128-bit record. For every fixed setup, its record-Hamming
distance equals the scalar distance: any nonzero component gives the lower
bound, and embedding a scalar message into one component gives equality.
There is no factor-sixteen failure union. Native seed 113 uses a reproducible
Mersenne Twister setup; the certificate concerns the ideal distribution,
not that particular seed or bit-Hamming distance.

## Complete proof accounting

The [disjoint table](HD_SHIFTED_PAIR_SWAPS.md#disjoint-coverage-checkpoint)
has 38 rows covering every nonzero active-group weight from 1 through 8192.
The manifest combines the two dense rows into one checkpoint, giving
37 entries. Together they cover all 11465830401 feasible matching profiles,
including every even/odd lane split. No feasible profile remains unresolved.

Fresh 256-bit and 384-bit arithmetic replays agree on all 10084 retained
per-weight or box-leaf rational ceilings. The final checker rebuilds
coverage and sums those ceilings with directed arithmetic. An independent
512-bit sum gives -55.09273850822352338... before upward rounding.
Setup-failure events need not be disjoint; the calculation uses a union
bound, not an independence assumption.

The proof combines exact sparse coefficients, fixed common-potential
bounds, and box bounds. The last gap, weights 4068 through 4095, reused a
frozen coefficient bank with finer support boxes. It required no change
to the encoder, setup distribution, or amount of mixing.

Legacy and sharper Toeplitz density bounds both dominate the same outer
counting measure. Using either in different bands is valid. Tail receipts
transfer only through the established unsigned operator domination or the
route-independent terminal bound. Fixed-oriented bounds are not imported.

The frozen manifest is
[hd_shifted_closed_cover_manifest.json](theory/hd_shifted_closed_cover_manifest.json),
with externally recorded SHA-256
`accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.
The final verification receipts are:

| Aggregation precision | Receipt | SHA-256 |
|---|---|---|
| 256 bits | [verified256](theory/hd_shifted_closed_cover_verified256.json) | `bab4c9e1656d86f4c68dfa4b015a70c810f0472a9c19d6d97761666945e3be65` |
| 384 bits | [verified384](theory/hd_shifted_closed_cover_verified384.json) | `446c2ce09488fcaf86c2070d86b8890d510a767a21cd134a2ef9c5c20672ba8e` |

Both report the same rational union ceiling and complete coverage. The
[checker](HD_SHIFTED_VERIFIED_COVER.md) authenticates previously observed
fresh arithmetic, then independently checks geometry and aggregation.
Hashes alone do not establish arithmetic validity. Its selected source
pins identify reviewed helpers, not every transitive numerical dependency.

From the theory directory, the inexpensive composition is:

```text
python -B hd_shifted_verified_cover.py --manifest hd_shifted_closed_cover_manifest.json --expected-manifest-sha256 accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743 --precision 384
```

Add `--fresh` to rebuild all coefficient bounds at both precisions without
using cached numerical bounds. Use `--output` with an unused filename to
retain a new receipt; frozen files are never overwritten.

## Matched native comparison

Each row below uses three serial runs on Peach CPU 2, an AMD Ryzen 9 7950X.
Each run times 501 complete forward encodings of 128-bit records after
warmup, with setup and allocation excluded. The builds use non-temporal
stores and explicit runtime parameters g4/h4. No benchmarks ran concurrently.
The shared host may still introduce noise.

| Construction / layout | Run medians, ms | Median of medians, ms |
|---|---|---:|
| Retained scalar-route reference | 3.996491, 4.118699, 4.048848 | 4.048848 |
| Same distribution, native output layout | 3.870926, 3.881285, 3.870415 | 3.870926 |
| Alternating paired route with independent swaps | 3.798160, 3.800274, 3.781779 | 3.798160 |

The complete candidate takes 6.19% less time than the retained reference,
but only 1.88% less than its native-output-layout refinement. Do not
attribute the whole improvement to the changed route or call it a 2x win.
Earlier shared-host candidate runs measured 3.66--3.78 ms; the matched
comparison above is the preferred current result.

The candidate binary has SHA-256
`103c38e52809ab4a1030150280e55080dc1fdba844b4ee17a7ca95e791ed7312`
and K20 checksum 7701295210980077712. Fresh checks pass 30720 input/state
basis directions, pair-mask and cyclic-boundary cases, small dense-oracle
tests, and the K20 oracle. Earlier sanitizer checks also pass. Build flags
and retained binary identifiers are in the
[iteration log](HD_STRUCTURAL_ITERATION.md#native-reproduction-identifiers).

## Earlier structural follow-up (now parked)

The following proposal records the next step at the time this certificate
was completed. It is historical, not the current work plan.

Preserve this checkpoint and the previous 72.3975-bit reference. The next
candidate is the existing eight-row RS[32,16] outer with normalized
D2(I+11^T)D1 column mixing, using the newly certified route and unchanged
h4 inner. Its static schedule removes about 19.7% of outer GFNI products
and reduces hot keys from 4.5 to 3.75 MiB. These are operation counts,
not an end-to-end latency prediction.

First bound the 17-active-group class with two active constituents per
group under this actual route. Retain constituent-pair correlations across
all 32 columns and charge concentrated pair geometries separately. Scalar-route
receipts and the Toeplitz counting measure do not certify this new outer.
If that bounded proof test passes, run a matched serial native comparison
with identical routing, swaps, inner, and output layout.
