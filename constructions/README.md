# Riffle construction registry

This directory records named construction candidates and their status. It is
the authority for candidate identity. A proof or benchmark applies only to the
candidate named in its manifest.

`registry.json` records the active candidate and machine-readable statuses.
Each candidate folder contains its own `manifest.json`.

The registry initially covers the \(g=4\) deterministic-permutation lineage.
It does not attempt to classify every historical experiment in this
repository.

| Candidate | Status | Folder |
|---|---|---|
| Riffle DP g=4 | Frozen reference | `riffle_dp_g4/` |
| Riffle DP-2Lap g=4 | Paused, open | `riffle_dp_2lap_g4/` |
| Riffle PacketMul-2Lap g=4 | Paused, open | `riffle_packetmul_2lap_g4/` |
| Riffle PacketMul-WrapMul-2Lap g=4 | Paused, open | `riffle_packetmul_wrapmul_2lap_g4/` |
| Riffle ParallelAcc g=4 | Diagnostic exploration | `riffle_parallelacc_g4/` |
| Riffle DenseOuter-ParallelAcc g=4 | Proved asymptotic warmup | `riffle_denseouter_parallelacc_g4/` |
| Riffle BCHBlockPerm-ParallelAcc g=4 | Proof-gym exploration | `riffle_bchblockperm_parallelacc_g4/` |
| Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4 | Proof-gym exploration | `riffle_shiftalpha64_bchblockperm_parallelacc_g4/` |
| Riffle RandomStepConv g=4 sigma=20 | Paused: one-lap boundary obstruction | `riffle_randomstepconv_g4_sigma20/` |
| Riffle RandomStepConv-2Lap g=4 sigma=20 | Paused: low compute value | `riffle_randomstepconv_2lap_g4_sigma20/` |
| Riffle RM512-RandomStepConv g=4 sigma=20 | Paused model exploration | `riffle_rm512_randomstepconv_g4_sigma20/` |
| Riffle FrozenRandom512-P2-RandomStepConv g=4 sigma=20 | Active model exploration | `riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/` |
| Riffle FieldPair512-P2-ScalarStepConv g=4 sigma=18 | Superseded entropy reference | `riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/` |
| Riffle FieldPair512-P4-ScalarStepConv g=4 sigma=16 | Active entropy optimization | `riffle_fieldpair512_p4_scalarstepconv_g4_sigma16/` |
| Random-outer parameter landscape | Active diagnostic exploration | `riffle_random_outer_landscape/` |
| Riffle S-Stripe RandomStepConv g=8 | Promising provisional bound | `riffle_striped_randomstepconv_g8/` |
| Riffle S-Stripe RandomStepConv g=4 | Promising provisional bound | `riffle_striped_randomstepconv_g4/` |
| Riffle TransposePacketShuffle-RandomStepConv g | Certified random-ensemble distance | `riffle_transpose_packetshuffle_randomstepconv/` |
| Riffle TransposeBitShuffle-RandomStepConv g | Certified random-ensemble distance | `riffle_transpose_bitshuffle_randomstepconv/` |
| Riffle SpectrumPerm-TransposeBitShuffle-RandomStepConv g | Active outer-replacement experiment | `riffle_spectrumperm_transpose_bitshuffle_randomstepconv/` |
| Riffle BCHPerm-TransposeBitShuffle-StructuredStepConv(t,s) | Active structured design | `riffle_bchperm_transpose_bitshuffle_structuredstepconv/` |
| Riffle FieldCheckpointAccumulate t=32 s=64 K=32 | Active proof exploration | `riffle_fieldcheckpoint_accumulate_t32_s64_k32/` |
| Riffle LDPCSplitState g=4 t=256 s=64 | Paused packet branch | `riffle_ldpcsplitstate_g4_t256_s64/` |
| Riffle LDPCSplitState PacketMul g=4 t=256 s=64 | Paused: low expected value | `riffle_ldpcsplitstate_packetmul_g4_t256_s64/` |

No deployable candidate is currently active. The current direction doubles
the local constituent block size at fixed rate and fixed total packet count.

The active entropy exploration is **Riffle
FieldPair512-P4-ScalarStepConv g=4 sigma=16**. It preserves the random-linear
marginals using entropy-optimal field multipliers. Its central distance
number remains diagnostic pending a complete support cover. Attempts to use
sigma 15 with five or six parity symbols expose a one-active-block
obstruction. The P2, FrozenRandom, RM512, and earlier 256-bit models remain as
reference diagnostics.

The active structural exploration is **Riffle S-Stripe RandomStepConv g=8**.
Exact regional coefficient diagnostics find the first full-stripe contours at
\(B=16384,\sigma=51\) and \(B=32768,\sigma=28\). These cells retain large
margin but still require outward-rounded numerical certification.

The packet-width-four stripe variant is also promising. Its first observed
trade points are \(B=512,\sigma=20\), \(B=1024,\sigma=15\), and
\(B=2048,\sigma=13\). These values use the same random-outer proof model.

The two active structural comparisons are **Riffle
TransposePacketShuffle-RandomStepConv g** and **Riffle
TransposeBitShuffle-RandomStepConv g**. The first retains fixed groups of
\(g\) blocks as packets and shuffles those packets in every row. The second
shuffles individual bits in every row before forming packets. Both now have
complete random-ensemble distance certificates; the bit-shuffle version is
stronger for \(g=8\) but permutes \(g\) times as many items.

The active outer-replacement experiment is **Riffle
SpectrumPerm-TransposeBitShuffle-RandomStepConv g**. It replaces every random
outer matrix by one fixed \([B,B/2]\) code. Each outer block receives an
independent coordinate permutation. The proof then depends on the fixed
code only through its weight spectrum. Goal 01 covers one active outer block
at \(B=256\). An ideal random spectrum gives about 46.62 bits at the higher
states. A modeled distance-38 spectrum gives about 51.95 bits there and 40.55
bits at the adjacent lower states. These are model diagnostics, not distance
certificates.

The active structured design is **Riffle
BCHPerm-TransposeBitShuffle-StructuredStepConv(t,s)**. It replaces the dense
random inner maps with either field multipliers or Toeplitz maps. Both
families preserve the exact one-vector transition used by the random-inner
first moment. At the preliminary point \((t,s)=(32,32)\), the modeled
one-active margin is 80.83 bits. The optimized transposed inner takes 27.14 ms
with FieldMul and 25.83 ms with Toeplitz on the recorded benchmark host. The
complete multi-active spectrum sum remains open.

The current low-cost inner candidate is **Riffle FieldCheckpointAccumulate
t=32 s=64 K=32**. It uses 64 accumulator lanes and one nonzero
\(\mathrm{GF}(2^{64})\) checkpoint multiplier every 1024 bits. The checkpoint
gives an exact zero/live transfer for every epoch occupancy. The one-active
outer-block diagnostic has 81.27 bits of margin at 9% distance and reaches
40 bits at about 14.049% distance. The complete multi-active outer-spectrum
sum remains open.

## Filing policy

Each candidate folder contains a manifest that identifies the exact map,
randomness, inherited mechanisms, and proof status.

Existing artifacts remain in `explorations/` and `scripts/`. Their candidate
manifests link to them. This policy preserves imports, hashes, and historical
references.

New candidate-specific reports and receipts should normally use these paths:

```text
constructions/<candidate>/proof/
constructions/<candidate>/receipts/
```

Shared code may remain in `scripts/`. A shared script must accept or state the
candidate explicitly. Candidate-specific scripts must include the candidate
slug in their filename.

Changing any of the following creates a distinct named candidate:

- the packet width;
- the packet permutation law;
- a packet-local randomizer;
- the convolution map;
- the number of laps;
- the retained-state or termination rule; or
- either outer or inner code.

## Initial transpose-packet result

At \(n=2^{20}\), relative distance 0.09, \(B=1024\), \(g=4\), and
\(\sigma=15\), the packed-group calculation gives provisional margin
\(195.28\) bits. The comparable S-Stripe point gave \(110.03\) bits.

The dominant occupation has one active outer block and is exact under the
random-outer relaxation. It is about 242 bits above the two-active-block
term. Packing domination is proved. A complete outward-rounded calculation
over all 2048 occupations certifies failure exponent
\(195.283900426038\) at distance threshold 188743.

Recommended next goal: audit the ensemble assumptions against the intended
encoder cost model, then determine the largest distance certified by the same
parameters.

## Certified g=4 outer-memory curve

At distance threshold 188743, complete outward-rounded calculations certify
these additional TransposePacketShuffle trade points:

```text
B=1024, sigma=13: lambda 83.087891407589
B=512,  sigma=15: lambda 57.641932041422
B=256,  sigma=18: lambda 45.400358851794
```

The adjacent lower-memory points do not reach 40 bits in the current
one-block bound. Under a dense bit-product proxy, the estimated outer-plus-
inner costs are respectively 1.225 billion, 0.726 billion, and 0.522 billion
products. The smallest outer point is therefore cheapest under this proxy,
despite its larger inner state.

Artifact:

- `constructions/riffle_transpose_packetshuffle_randomstepconv/TRADEOFF_G4.md`.

Recommended next goal: benchmark the three certified cells with common
low-level kernels. For a point that reduces both \(B\) and \(\sigma\), explore
padded intermediate outer lengths between 512 and 1024.
