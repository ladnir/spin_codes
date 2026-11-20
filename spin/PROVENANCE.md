# Kernel generation record

SPIN is an independently licensed MIT library; see [LICENSE](LICENSE). The
package has no runtime or build dependency on its consumer projects or on the
research workspace.

`PacketRs8T64` (identifier 8, Turbo SPIN) promotes the common RS8 outer from
the research experiment `k20_two_ms/rs8`, with the exact s22 inner at every
size. The package-owned `src/rs8` module preserves the field basis, sampling
domain, and unrolled flat-nine-product outer schedule. `src/certified22`
holds the routing and exact s22 inner: its setup and scalar reference preserve
the original research streams, and its packed-state route was ported from the
`k20_small_outer` experiment. The AVX2, AVX2+GFNI, NEON and portable kernels
were written in the package against the scalar oracle. The two packet inner
headers in `src/packet` and the setup sampler in `src/kernels/SetupRandom.h`
are retained from the earlier packet and IMT imports; Turbo's setup streams
depend on them.

`tools/generate_lane_major22.py` regenerates `src/certified22/LaneMajor22.h`
from `productionMap.h`. Ordinary builds need neither Python nor research
sources.

SPIN 0.4 removed the retired families: paper SPIN (BCH outer, IMT inner;
identifiers 1-3), the BCH packet family (4), and the packet profiles s20 (5),
K16 paired-s15 (6), and RS16 s22 (7), together with their importers, overlays,
generic XOR transpose, and research experiments. Commit `584a7fb6` is the last
source tree containing them; their generation records (input hashes, import
commits `3ced6dc0` and `2582a9fb`) are in that tree's `spin/PROVENANCE.md`.
