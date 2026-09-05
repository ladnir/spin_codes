# Spectrum-computable outer exploration at B=1024

Status: exploratory. No variant in this directory replaces the current BCH
outer until an end-to-end distance certificate and an encoder benchmark pass.

The downstream construction is the existing bit-transpose design with the
RM2Sub `s=16` inner. This exploration changes only the rate-1/2 outer map from
512 bits to 1024 bits.

## Weight-transition interface

Let `S[h]` be the expected number of outer codewords of weight `h`. Let a
randomized stage `F` have input-output weight enumerator `A_F(h,t)`. A uniform
interleaver before `F` gives the transition

```text
T_F(h,t) = A_F(h,t) / binom(B,h).
```

The expected spectrum after the stage is

```text
S_next[t] = sum_h S[h] T_F(h,t).
```

Independent interleavers make this a weight-only Markov composition. This
interface supports local mixers, accumulators, repeaters, block codes, and
expanders when their input-output enumerators are available.

The computed values are ensemble expectations. They do not certify every
fixed schedule. Markov's inequality converts a cumulative expected mass of
`2^-lambda` into a failure bound of at most `2^-lambda` for a sampled schedule.
A fixed production schedule still requires either this sampling argument or a
separate deterministic certificate.

## Named baselines

### Outer ParityMix8-L

Start from `(message,0)`. Each of `L` layers independently permutes all 1024
coordinates, groups them into sets of eight, and applies

```text
y_i = x_i + sum_j x_j
```

within every group. A group of input weight `r` has output weight `r` for even
`r` and `8-r` for odd `r`. The exact one-layer joint enumerator is the
coefficient table of

```text
(sum_r binom(8,r) x^r y^(r if r is even else 8-r))^(B/8).
```

Each layer costs 1.875 XORs per output bit before permutation overhead.
Four layers therefore cost 7.5 XORs per output bit.

### Outer RA-2 and Outer RAA-3

Repeat every message bit twice. RA-2 applies one independently interleaved
accumulator. RAA-3 applies two. Repetition maps input weight `w` to weight
`2w` exactly.

For a length-`B` accumulator, the exact input-output enumerator is

```text
A_acc(h,t) = binom(B-t,floor(h/2)) binom(t-1,ceil(h/2)-1).
```

One accumulator costs about one XOR per output bit, excluding its interleaver.

### Random BA-3(s) baseline

Split the 512-bit message into `s`-bit blocks. Encode each block with an
independently sampled systematic `[2s,s]` binary matrix. Apply two independently
interleaved accumulators to the resulting 1024-bit word.

The expected local homogeneous spectrum is

```text
L_s(z) = 1 + 2^-s ((1+z)^(2s) - (1+z)^s).
```

The block stage has spectrum `L_s(z)^(512/s)`. A direct dense implementation
costs about `s/4` XORs per output bit. The two accumulators add about two more.
This ensemble remains a comparison baseline. It is not a preferred encoder.

### Outer DBA-3(C)

Fix an injective designed constituent `C:[s]->[2s]` with homogeneous weight
enumerator `W_C(z)`. To encode one 512-bit Riffle outer block, apply `C` to
`512/s` constant-size pieces, then apply two independent length-1024
interleavers and two length-1024 accumulators. Reuse the fixed constituent,
but sample the two BA interleavers independently for each Riffle outer-block
instance.

Before bit transpose, retain the existing independent output-coordinate
permutation for each 1024-bit outer word. Consequently, the downstream proof
uses only the complete weight spectrum of this local BA code. The BA
permutations do not span the full `N=2^21` codeword.

For one Riffle outer block, the constituent stage has exact spectrum
`W_C(z)^(512/s)`. The constituent needs no encoder randomness and does not
grow with the final code length.

### Outer EA-3(d)

Retain the 512 message bits and append 512 parity checks. Each check samples a
size-`d` subset of message positions and XORs those bits. Apply two independently
interleaved accumulators to the resulting 1024-bit word.

This systematic expander costs about `(d-1)/2` XORs per output bit. The two
accumulators raise the nominal total to `(d-1)/2+2`.

## Exact B=1024 results

The table reports the base-2 logarithm of the cumulative expected number of
nonzero codewords through the stated weight. A more negative value is better.
The DBA rows are the complete expected spectra of proposed length-1024 outer
encoders. They are the correct local inputs to the end-to-end Riffle analysis.
They are not end-to-end distance certificates by themselves.

| Outer candidate | First point coefficient at least 1 | First cumulative mass at least 1 | Cumulative through 37 | Cumulative through 95 |
|---|---:|---:|---:|---:|
| RA-2 | 1 | 1 | 54.67 | 140.72 |
| RAA-3 | 43 | 36 | 0.29 | 54.04 |
| ParityMix8-4 | 86 | 86 | -4.27 | 1.83 |
| ParityMix8-5 | 115 | 115 | -11.72 | -5.58 |
| ParityMix8-6 | 115 | 115 | -18.91 | -12.78 |
| ParityMix8-4 + Acc | 115 | 115 | -15.54 | -14.11 |
| ParityMix8-4 + Acc + Acc | 115 | 115 | -23.79 | -22.17 |
| BA-3(16) | 111 | 111 | -13.85 | -11.83 |
| BA-3(32) | 114 | 114 | -27.25 | -23.67 |
| DBA-3(EBCH [8,4,4]) | 95 | 95 | -8.79 | 1.22 |
| DBA-3(EBCH [32,16,8]) | 111 | 111 | -20.47 | -15.66 |
| DBA-3(conjectured EBCH `[64,32,12]`) | 114 | 114 | -33.03 | -25.36 |
| DBA-3(EBCH [128,64,22]) | 114 | 114 | -56.05 | -40.72 |
| EA-3(8) | 115 | 115 | -11.47 | -9.80 |
| EA-3(12) | 115 | 115 | -16.71 | -14.98 |
| EA-3(16) | 115 | 115 | -22.09 | -20.27 |
| EA-3(24) | 115 | 115 | -33.12 | -31.05 |
| EA-3(32) | 115 | 115 | -44.41 | -41.93 |

RA-2 fails because two repeated impulses become an interval under the first
accumulator. Short intervals occur with high multiplicity. RAA-3 mitigates but
does not remove this rate-1/2 failure mode.

ParityMix8 has persistent low-weight orbits. A weight-two word survives a layer
when both bits land in the same group. Later layers are dominated by related
weight-14 spikes. Each extra layer removes about seven bits from the measured
low-weight tail in this range.

BA-3 gives the best measured shape among these baselines. Its low-weight floor
depends on `s`. The `s=16` and `s=32` results agree with the earlier empirical
rule of roughly `3s/4` failure bits.

The local DBA diagnostics show that designed constituents improve the floor
without local encoder randomness. At local dimension 16, the fixed extended
BCH constituent improves the weight-37 diagnostic from `-13.85` for random BA
to `-20.47`. These outer-only values do not replace the end-to-end Riffle
calculation.

## End-to-end one-active gate

The complete B=1024 spectra were composed with the existing per-block output
permutation, bit transpose, region permutations, and RM2Sub `s=16` inner. The
calculation covers messages with exactly one active Riffle outer block.

| BA constituent | One-active margin | Dominant local output weight | Decision |
|---|---:|---:|---|
| EBCH `[8,4,4]` | 2.37 bits | 2 | Reject. |
| EBCH `[32,16,8]` | 19.65 bits | 2 | Reject. |
| Conjectured EBCH `[64,32,12]` | 39.20 bits | 3 | Near miss. |
| EBCH `[128,64,22]` | 80.23 bits | 6 | Continue. |

The two smaller constituents fail because the local BA ensemble retains rare
weight-two outputs after its second accumulator. Across 2,048 possible active
outer-block locations, those events exceed the 40-bit budget. The EBCH
`[128,64,22]` version clears the one-active gate, but it has no end-to-end
certificate yet. Multi-active messages remain open.

The outer is not required to have excellent standalone distance or a negligible
standalone failure probability. Its design target is the composed Riffle bound

```text
sum_h (expected outer multiplicity at weight h)
      * (Riffle-inner bad-event bound conditioned on weight h).
```

Thus minimum distance and outer-only cumulative mass are diagnostics, not the
objective. Weight two is presently expensive because the inner gives it almost
no additional suppression; moderate and large outer weights can be much more
common because the inner rapidly suppresses them. We should choose the cheapest
outer spectrum that makes the complete sum clear 40 bits, with enough reserve
for the still-open multi-active cases. The 80.23-bit DBA result is therefore an
upper anchor, not evidence that Riffle needs an outer that strong.

The length-64 row is an exploratory spectrum proxy, not the enumerator of a
specified fixed code. Start with the exact `[64,36,12]` extended-BCH spectrum
and choose a random dimension-32 subcode conditioned to contain the all-ones
word. Every such subcode retains minimum distance at least 12. Each other
parent codeword is retained with probability

```text
(2^31 - 1) / (2^35 - 1),
```

which is essentially `1/16`. The committed proxy rounds that conditional mean
symmetrically and adjusts the weight-32 coefficient by two so that its total
mass is exactly `2^32`. It should be read as a curve-filling conjecture only.
Its 39.20-bit result, dominated by outer weight three, places the useful BA
threshold very near a 64-bit, distance-12 constituent.

### Multi-active full-spectrum first pass

The BA spectrum can be inserted directly for one active outer block because,
after the block's output-coordinate permutation, its total weight completely
determines the distribution seen by the transpose. This is the 39.20-bit
calculation above.

For multiple active blocks, the first bulk plug-in using the single maximum
spectrum-density ratio is vacuous. The ratio is maximized by the rare BA output
weight three and charges its 435.75-bit density surcharge independently to
every active block. This is an artifact of the compression, not refutation
evidence.

The next first-pass calculation retains every BA spectrum coefficient through
a Holder change of measure from uniform outer words. It is materially stronger
than the maximum-ratio endpoint, but remains vacuous at high occupation: its
dominant row is 2,048 active blocks with a diagnostic margin of about
`-388,035` bits. The Holder moment still permits the rare low-weight tail to be
perfectly correlated with the dense inner bad event. The appropriate next
compression is therefore a two-piece spectrum split: enumerate or sum the rare
low-weight BA tail, and apply the regular bulk bound only to the central body.
This remains entirely an analysis of BA as the local outer inside Riffle.

The endpoint split `min(w,1024-w) <= 97` gives the following current
decomposition:

| Riffle message class | Margin | Qualification |
|---|---:|---|
| Exactly one tail BA block | 39.20 bits | Complete one-active spectrum sum. |
| Body blocks only | 229.10 bits | Complete body-spectrum Holder bound. |
| At least one tail block plus at least one body block | 235.52 bits | Direct fixed-region-weight envelope; permits arbitrary tail supports. |
| Exactly two tail blocks and no body block | 63.34 bits | Fixed-support bound via exact IID conditioning penalty. |
| Three or more tail blocks and no body block | 43.87 bits | Count only; no inner help. |

The direct mixed calculation removes the previous arbitrary-shift hypothesis.
Fix one region and let `F[w]` denote its transfer matrix for a uniformly
located input support of weight `w`. Suppose the region contains `b` body
positions and `j` forced tail ones. The ambient-uniform body values contribute
`Bin(b,1/2)` ones. The region permutation makes the combined support uniform
conditioned on its total weight. Therefore define

```text
G[b,j] = E[F[j + Bin(b,1/2)]].
```

The recurrence

```text
G[b+1,j] = (G[b,j] + G[b,j+1]) / 2
```

computes all such matrices. For each `b`, the analyzer takes the entrywise
maximum over every feasible `j`. This envelope permits the tail count to vary
adversarially between regions. It still retains the independent random
permutation within each region. The resulting mixed-class margin is 235.52
bits after the complete body-spectrum Holder step and the exact ensemble
expectation for all nonempty tail choices.

The recurrence audit compares `G[b,0]` with the existing regular transfer
table for every `b`. Their maximum log-domain difference is `3.83e-12`.
Direct binomial sums for selected `(b,j)` pairs agree with the recurrence to
`1.78e-15`. Thus every message class is covered under the recorded transfer
model. The combined margin remains 39.15 bits because the exactly-one-tail
class contributes 96.2 percent of the bound.

The same one-active composition gives the following comparison. These are
nearest-binary64 diagnostics, not outward-rounded certificates.

| Outer candidate | One-active margin | Nominal XORs/output | Interpretation |
|---|---:|---:|---|
| Random BA-3(16) | 7.15 | about 6 | Reject; weight two dominates. |
| Random BA-3(32) | 21.22 | about 10 | Reject; weight two dominates. |
| EA-3(12) | 9.90 | 7.5 | Reject. |
| EA-3(16) | 15.31 | 9.5 | Reject. |
| EA-3(24) | 26.46 | 13.5 | Reject. |
| EA-3(32) | 37.92 | 17.5 | Just below target. |
| EA-3(34) | 40.82 | 18.5 | Barely clears the one-active gate. |
| EA-3(36) | 43.75 | 19.5 | Clears with modest reserve. |
| ParityMix8-4 + Acc + Acc | 16.92 | 9.5 | Reject; also permutation-heavy. |
| DBA-3(conjectured EBCH `[64,32,12]`) | 39.20 | benchmark open | Near miss; weight three dominates. |
| DBA-3(EBCH `[128,64,22]`) | 80.23 | benchmark open | Strong upper anchor; likely overbuilt. |

EA-3 has the same main distance crossing across the tested degrees. Increasing
`d` suppresses the rare low-weight floor. The end-to-end calculation, rather
than the outer-only weight-95 statistic, places its one-active threshold between
degrees 32 and 34.

## Expander composition

The draft expander enumerator is in `enumerator_paper/Expander.tex`. Suppose
each of `n` output checks samples a size-`d` subset of `k` input positions. For
an input support of weight `w`, define

```text
odd(w)  = sum_j binom(w,2j+1) binom(k-w,d-2j-1)
even(w) = sum_j binom(w,2j)   binom(k-w,d-2j).
```

The normalized expected enumerator is

```text
A_exp(w,h) = binom(k,w) binom(n,h)
             (odd(w)/binom(k,d))^h
             (even(w)/binom(k,d))^(n-h).
```

The local TeX draft omits the denominator `binom(k,d)^n` if its displayed
quantity is intended as an expected enumerator. The numerator is usable after
that normalization audit. An expander can feed one or two accumulator stages
through the same transition interface.

## Artifacts

The exact analyzer is `scripts/analyze_riffle_paritymix8_outer_spectrum.py`.
The fast expander analyzer is
`scripts/analyze_systematic_expander_accumulate_spectrum.py`.
The designed-BA analyzer is
`scripts/analyze_designed_block_accumulate_spectrum.py`.
The one-tail mixed analyzer is
`scripts/analyze_riffle_one_tail_body_shell_conditioned.py`, and the deliberately
coarse all-mixed diagnostic is
`scripts/analyze_riffle_multi_tail_body_shell_holder.py`.
The closing direct analyzer is
`scripts/analyze_riffle_body_arbitrary_tail_direct.py`. Its recurrence audit is
`scripts/audit_riffle_arbitrary_tail_region_recurrence.py`.
The disjoint-class summary is
`scripts/summarize_riffle_tail_body_certificate.py`.
The proof argument is recorded in `direct_arbitrary_tail_lemma.md`.
The `receipts` directory contains the full coefficient tables for every row in
the comparison. The first analyzer uses exact integer polynomials and exact
rational composition. Floating-point arithmetic is used only to serialize its
logarithms. The expander analyzer uses stable binary64 log-domain composition.
At length 64, its output agrees with exact rational composition within
`1.9e-13` bits.

## Explicit shortened-XBCH refinement

The conjectured constituent can be replaced without changing the Riffle
architecture. The replacement is an explicit `[64,32,12]` code generated by
`(I_32 | A)`. The 32 rows of `A` come from the shortened-XBCH diffusion matrix
published in US patent application 20020101986, paragraphs 45--59. The local
files are:

- `scripts/xbch64_32_philips_generator_rows.txt`;
- `scripts/xbch64_32_philips_spectrum.csv`.

The exact enumeration visits all `2^32` codewords. It finds minimum distance
12 and 787 words of weight 12. For comparison, the best self-dual constituent
tested here has 1,312 words of weight 12. The explicit spectrum is not
complement-symmetric, so every calculation uses the complete table rather than
a symmetry shortcut.

The generator now has an independent reconstruction certificate. The script
`scripts/reconstruct_xbch64_32_philips.py` starts from

```text
g(x) = x^27 + x^22 + x^21 + x^19 + x^18 + x^17
       + x^15 + x^8 + x^4 + x + 1.
```

It verifies that `g` divides `x^63+1` and has roots `alpha^1` through
`alpha^10` in `GF(64)`. The BCH bound gives distance at least 11. Shortening
to dimension 32 preserves this distance, and an overall parity coordinate
raises the lower bound to 12. Systematic reduction produces `(I_32 | B)` with
28 parity columns.

The published matrix contains all 28 columns of `B` up to row and column
permutations. Degree invariants uniquely identify zero-based columns
`5,9,28,31` as the four added columns. Two rounds of bipartite incidence
refinement uniquely recover the remaining permutations. Thus puncturing those
four coordinates in the published `[64,32]` code recovers the reconstructed
extended `[60,32,>=12]` BCH code up to a coordinate permutation. The receipt is
`xbch64_philips_reconstruction.json`. The exhaustive spectrum supplies a
weight-12 witness, so the published code has minimum distance exactly 12.

Keep the original two length-1024 accumulators. The following five classes are
disjoint and exhaustive under the endpoint-tail split.

| Message class | Nearest margin | Outward lower margin |
|---|---:|---:|
| Exactly one tail block and no body block | 40.464902278373 | 40.464902278356 |
| Exactly two tail blocks and no body block | 65.828916919490 | 65.828399762063 |
| At least three tail blocks and no body block | 46.912164341081 | 46.912124928575 |
| Body blocks and no tail block | 230.086343012605 | 230.086343012485 |
| At least one body block and one tail block | 237.522306338528 | 237.522280635105 |

The outward calculation uses 256-bit Arb intervals for combinatorial ratios
and transcendental functions. Every nonnegative binary64 array operation is
followed by `nextafter(+infinity)`. Conditioned region matrices carry separate
power-of-two exponents, which prevents underflow without approximate rescaling.

The complete bad-event probability is at most
`2^-40.448462678054`. Thus the explicit constituent clears the 40-bit target
under the recorded inner-transfer model. The retained receipt and its audit
are:

- `receipts/outward_full_endpoint97_dba_xbch64_philips.json`;
- `receipts/outward_full_endpoint97_dba_xbch64_philips_audit.json`.

The outward arithmetic layer and the constituent-provenance layer are closed.
A paper theorem must still import the recorded activation and live-spectrum
lemmas.

## Proof-backed encoder implementation

The implementation uses the explicit constituent above. It does not replace
the constituent with the earlier `[256,128]` benchmark code. Each length-1024
outer block contains 16 copies of the `[64,32,12]` constituent. Two independent
interleaved accumulators follow the constituent layer.

The transposed encoder applies the following data path:

1. The RM2Sub inner emits finalized 128-bit blocks in reverse order.
2. Four sequential route streams partition the blocks into 512-outer-block
   tiles. The route stores each destination in 24 bits.
3. One tile scatter reconstructs an 8 MiB outer tile.
4. The encoder evaluates both accumulators in batches of four outer blocks.
5. A generated scalar circuit evaluates each transposed constituent.

The selected constituent circuit uses 252 block XORs. Its transposed schedule
has peak live count 23 and live area 1779. The generator seeds are 136 and 605
for synthesis and scheduling, respectively. The receipt is
`receipts/xbch64_philips_transpose_circuit_selected.json`.

The benchmark compares the fused encoder with a staged reference. The reference
materializes the inner word, applies the explicit route, evaluates both
accumulators, and uses dense constituent rows. The two implementations produce
the same checksum, `0xae7fddb94ec50842`.

On Peach's Ryzen 9 7950X, GCC 15.2 with `-O2 -march=native` is faster than
`-O3`, `-Ofast`, LTO, and profile-guided optimization. The selected route uses
no software prefetch. Linux huge-page advice for the 32 MiB route buffer and
the 8 MiB tile reduces address-translation overhead. The final uncontended
101-trial run measured a 10.705 ms median and a 10.496 ms minimum. The best
observed 51-trial window measured a 10.606 ms median.

The following alternatives did not improve the complete encoder:

| Alternative | Median time (ms) | Decision |
|---|---:|---|
| Selected circuit, `-O2`, huge-page advice | 10.705 | Retain. |
| Native 32-bit route offsets | 11.396 | Reject; extra metadata traffic dominates. |
| Four-way 24-bit route unpacking | 11.400 | Reject; instruction-level parallelism does not repay the extra unpacking. |
| PGO build | 11.109 | Reject. |
| Explicit reserved huge pages | 10.609 | Reject; no gain over ordinary memory with page collapse. |
| 254-XOR, peak-live-21 circuit | 10.805 | Reject; the lower live count is not an end-to-end win. |
| 4 MiB tile | 10.922 | Reject; eight route streams cost more than the smaller tile saves. |
| 16 MiB tile | 11.158 | Reject; the larger random-write working set costs more than two route streams save. |

## Recommended next step

Register the complete encoder benchmark as an optional libOTe performance
target. Further tuning should target the two route passes; the constituent
circuit and accumulator batch are no longer the dominant unresolved choices.
Retain the three-accumulator 52.09-bit variant only as a fallback with extra
margin.
