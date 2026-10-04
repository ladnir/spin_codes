# Current certified research baselines

The selected unified library profile is now
[`PacketRs8T64`](../spin/RS8.md): the same RS8/GF32 outer at K16/K18/K20,
with paired-s15 at K16 and exact s22 above. The promotion preserves the
[three-size research comparison](../spin/experiments/rs8_sizes/README.md)
and all earlier controls below. Its margins are 61.8841/58.2456/51.8225 bits
at 10% relative distance; the public-API timings are 0.100/0.461/2.862 ms.
These are separate certificates, not the RS16 margins in the retained table.

Originally recorded 2026-10-02; K18's unchanged code received a stronger proof
on 2026-10-03. These are the selected closed four-bit research results. Retain
their exact maps, proof sources, and measured kernels as controls. Further K18
small-state and shorter-step construction tuning is parked at the user's request.

All listed codes have rate 1/2, use four-bit packets, and certify at least
10% relative distance. The margin bounds whole-code setup failure under
the specified independent setup ensemble. It is not a per-message bound
or a measured distance of one seeded code.

| Message size K | Outer group, input/output bits | Physical inner t / s | Precomputed transpose | Complete margin | Minimum distance |
|---|---|---|---:|---:|---:|
| 2^16 | 128 / 256, four GF16 RS[16,8] rows | 64 / 15, paired map | **0.098950 ms** | **62.041173 bits** | 13,108 |
| 2^18 | 128 / 256, four GF16 RS[16,8] rows | 64 / 22, exact certified map | **0.466953 ms** | **60.750999 bits** | 52,429 |
| 2^20 | 256 / 512, eight GF16 RS[16,8] rows | 64 / 20, retained extended map | **3.258973 ms** | **68.103757 bits** | 209,716 |
| 2^20, new research winner | 128 / 256, four GF16 RS[16,8] rows | 64 / 22, exact certified map | **2.962029 ms** | **58.750377 bits** | 209,716 |

Here `t` is the number of binary coordinates processed per physical step;
`s` is the persistent state dimension. Encoding carries 128-bit XOR
elements through these binary maps. Timings exclude setup and allocation.
They come from separate serial Ryzen 7950X campaigns, not one same-session
scaling comparison. K16 uses normal pages; K18 and K20 use the recorded
large-page preference. See each report for compiler, repetitions, and
buffer details. No new benchmark was run for this index.

## K16: frozen paired15 winner

- [Source-pinned checkpoint](../spin/experiments/k16_codesign_100us/checkpoint/README.md): `k16-rs16-paired15-v1`, mode 52, 220 pinned files.
- [Construction and complete proof](workstreams/k16_codesign_100us/proof/PAIRED_S15_CONSTRUCTION.md).
- [Implementation and benchmark commands](../spin/experiments/k16_codesign_100us/README.md).
- [Research report](workstreams/k16_codesign_100us/README.md), including the separate 98.950 us holdout.
- Subsequent library integration is documented in [PAIRED15.md](../spin/PAIRED15.md), parameter identifier 6. The table retains the frozen research timing, not a new public-API timing.

The local certificate is
`workstreams/k16_codesign_100us/proof/paired-s15-drop10-whole-p256.json`,
SHA-256 `89bed0c4442349a520d54b1aa6e914eac326beb3530123f8273f78fc764b867d`.
The exact margin recorded by the checkpoint is 62.04117326862128 bits.
The faster approximately 95 us RS8 XOR candidate lacks a complete proof
and is not the selected baseline.

## K18: exact certified22 on the faster small-outer route

- **Proof update, 2026-10-03:** [routing-conditioned clipping](workstreams/k18_margin60/README.md) raises the complete margin to **60.7509985457 bits**, without changing the construction or encoder. The 0.466953 ms measurement therefore remains the matched implementation result; no new timing is claimed.
- [Matched implementation and proof report](../spin/experiments/k18_certified22/README.md).
- Selected kernel: [K18Fast.cpp](../spin/experiments/k18_certified22/K18Fast.cpp), SHA-256 `a56e2f4be564b36a316ad01f411398e8d25f63590405fde2ee0289d60d696148`.
- [Exact-map requirements](../spin/CERTIFIED22.md).
- [Full replay instructions](workstreams/k16_design/README.md#power-of-two-checkpoint-k--218).

The local receipt is `../tmp/rs16-s22-k262144-whole-fresh-p256.json`,
SHA-256 `20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b`.
Its exact dyadic union gives 42.090840060352 bits across all 2048 occupancies.
This original certificate is preserved. The strengthened complete receipt is
`../tmp/k18-margin60-fractional-whole-p256-v2.json`, SHA-256
`44fcdb0389d5cc7be870bdfcbeb0ec5500193706b7e690313054971d052a01d0`.
It replaces 82 intermediate-occupancy endpoints by new outward bounds and
sums all 2048 endpoints exactly; the new bound is 60.75099854572545 bits.
The certified map identity is
`e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.

The same code is now available as library parameter identifier 7; the
[promotion report](../spin/CERTIFIED22.md) records separate public-API timings.
Preserve the actual
expansion and feedback maps: replacing them by another basis with the
same expansion code does not automatically preserve the certificate.
The approximately 0.40 ms paired15/16 K18 kernels and all tested
feedforward/shorter-step candidates remain uncertified at this target.

## K20: larger outer and s20

- [Matched construction, certificate, and 3.258973 ms campaign](workstreams/k16_design/README.md#matched-k20-result-10-distance-in-326-ms).
- [Standalone fresh reproducer](workstreams/k16_design/reproduce_rs16_k20.py).
- Measured driver: `fused-nt` in [rs16x8_border](workstreams/k16_design/implementation/rs16x8_border/), using the fused route and two-plane outer.
- The same construction is integrated as parameter identifier 5; [the library report](../spin/README.md) records separate public-API measurements.

The first full receipt is `../tmp/rs-k20-outer256-s20-whole-fresh-p256.json`,
SHA-256 `97396aee84ff16c05924df80ce6b4c7722bff9ebd16643e4c15bcaf69d508e84`.
An independent complete replay at `../tmp/rs16-k20-reproduction/whole-p256.json`
has SHA-256 `08915a0e2340843999fceff2a94e2ba6927fa1d3375f2408f333c88f4d571642`.
Both give the same 4096 occupancy endpoints and 68.103757217935-bit union.

The matched source archive remains local at
`../tmp/spin-k20-goal-campaign-source.tar.gz`, SHA-256
`3e03aa3308264a5267e0469ca41ed35104301313aa207916f07d7b7da57431d5`.
Do not replace this result with the smaller-outer s16 timing at K20;
that different construction lacks the claimed complete certificate there.

## Preservation and next work

The subsequent [K20 small-outer proof](workstreams/k20_small_outer_screen/README.md)
closes the same 128-to-256-bit outer and exact s22 maps used at K18, now at
K=2^20. Its full 8192-occupancy union gives **58.7503767402 bits** of margin
at 10% distance, meeting the accepted 50-bit target. The whole receipt is
`../tmp/k20-small-outer-whole-v1.json`, SHA256
`90492fc15834d2e28a9ddc8289a68dbf7bd17984146524016f9361a9ec54fe65`.
The [initial matched benchmark](../spin/experiments/k20_small_outer/README.md)
measures 5.213749 ms for the unchanged K18 kernel at K20. The subsequent
[exact-map optimization](../spin/experiments/k20_small_outer/OPTIMIZATION.md)
uses streaming stores, carried state-update outputs, and unrolled outer mixing.
It measures **2.962029 ms**, versus **3.299650 ms** for the larger-outer library
profile in the same holdout: 10.23% less time. The construction and certificate
are unchanged. This is the new K20 research performance winner at the accepted
50-bit target. It is now available as `PacketRsT64S22`, alongside the retained
larger-outer profile. The [promotion report](../spin/CERTIFIED22.md) measures
0.465073 ms at K18 and 2.955466 ms at K20 through the public API. The larger-outer
row retains its original timing and higher-margin certificate.

This index records existing results; it does not rerun or extend them.
Source and receipt identities above were checked when this index was made.
Raw receipts, archives, and timing samples remain ignored, not committed.
The K16 checkpoint has its existing source manifest. K18's selected kernel
and certificate identities are recorded above; this index is not a new
full transitive-dependency manifest. Library validation and source provenance
are recorded separately in the linked promotion report.

Keep the K18 joint-bound and shorter-step investigations as
parked research. Neither should trigger more tuning by default.

## Additional K16 checkpoint: eight-bit packets

The [width-eight study](workstreams/packet8_codesign/iteration5/README.md) now
has a matched full certificate and implementation. It uses the wider
256-to-512-bit outer, a 24-bit state, and eight-bit routing packets.
At K=2^16 it certifies minimum distance 13108 with **68.893748662558 bits**
of whole-code margin. All 256 occupancies are covered by outward bounds;
an independent global replay also passes the exact 40-bit threshold.

Its precomputed transpose takes **98.1625 us** on 128-bit elements. The
four-bit control takes 99.591 us in that same fresh-seed campaign. These
are effectively comparable latencies, not a large demonstrated speedup.
The earlier uncertified byte kernel takes 91.4005 us in the same campaign.
See the [performance record](workstreams/packet8_codesign/iteration5/PERFORMANCE.md).

This additional research checkpoint does not replace the selected baselines
above or promote a production backend. Its receipt hashes, exact dyadic
union, reproduction commands, and implementation audit are in the linked
report. The frozen K16 baseline's 220 source pins and archived certificate
were rechecked successfully after this work.

The subsequent [exact-map optimization](workstreams/packet8_codesign/iteration6/README.md)
measures **93.895 us**, against **97.8175 us** for that frozen byte-packet
control in a new matched campaign. It simplifies the field-randomizer circuit
without changing the code or its 68.893749-bit certificate. A separately
proved seven-diagonal MDS randomizer measures 96.831 us with a transported
68.350804-bit margin; it is not the selected implementation. These remain
research checkpoints, with no production promotion.
