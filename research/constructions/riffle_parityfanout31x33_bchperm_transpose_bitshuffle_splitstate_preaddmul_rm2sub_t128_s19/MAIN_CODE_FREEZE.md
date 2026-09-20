# Main code freeze

Status: frozen **Structured SPIN** implementation candidate as of 2026-08-31.

The frozen construction was developed as Riffle ParityFanout-31x33 and is now
the instance **Structured SPIN (B=256, t=128, s=19)**. It uses a 256-to-128
structured outer, per-block coordinate permutation, bit transpose, independent
region permutations, and the RM2Sub inner with step size 128 and state size 19.
The message and code lengths are `2^20` and `2^21` 128-bit blocks.

The optimized transposed implementation has the following locked results:

- independent dense-inner correctness: pass;
- independent staged end-to-end correctness: pass;
- checksum: `0x95c9d722a9539fef`;
- 21-trial Peach median: 10.823872 ms;
- fine-grid floating proof margins: 66.270087 bits at distance 0.09,
  61.069584 bits at distance 0.10, and 55.864642 bits at distance 0.11.

The benchmark receipt records the exact source hashes and the frozen remote
workspace `/tmp/riffle_parityfanout31x33_s19_v1`. The exact sources have been
copied into `frozen_source/`; all three local hashes match both Peach and the
performance receipt. `frozen_source/SOURCE_MANIFEST.json` records the snapshot.

The recorded correctness and benchmark run used these exact source hashes.
Cleanup must preserve the frozen snapshot and begin from a separate working
copy. The first cleanup checkpoint must reproduce the checksum before any
behavior-preserving refactor begins.

Cleanup is behavior preserving. A change to the construction, schedule,
randomness distribution, constituent, or parameter creates a separately named
variant. It does not modify this frozen record.

The current proof numbers are diagnostics, not the final paper theorem. The
remaining certificate work is:

1. replace nearest-binary64 evaluation with outward-rounded arithmetic;
2. certify a spectrum or sufficient envelope for the actual structured outer
   code, replacing the modeled even-floor spectrum.

Canonical artifacts:

- `README.md`;
- `PROOF_STATUS.md`;
- `PERFORMANCE.md`;
- `benchmark/parityfanout31x33_s19_integrated_peach_7950x.json`;
- `frozen_source/SOURCE_MANIFEST.json`;
- `receipts/floating_distance_curve_d09_d11.json`.
