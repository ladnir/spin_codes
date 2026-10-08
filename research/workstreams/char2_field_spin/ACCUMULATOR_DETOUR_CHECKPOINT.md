# Accumulator detour: preserved research checkpoint

Frozen on 2026-10-06 at the user's request. This branch is no longer the
main high-distance design search. Continue at
[the blockwise SPIN restart](HD_BLOCKWISE_RESTART.md). No source, proof
helper, benchmark receipt, or certified winner was deleted or replaced.

## Construction and measured frontier

The selected research encoder has rate 1/2 and K=2^20. Each 128-bit record
contains eight GF(65536) lanes. GRS[64,32] outer rows feed fixed affine
placement of four-row groups, one weighted prefix, one regional permutation,
and a two-additive-state post-inner. All 2K output records are retained.
The setup is precomputed; the prefix and post states carry globally.

Full one-core forward encoding measured 8.327--8.360 ms on Peach, versus
the retained 16.622 ms FLOCK RS baseline. The instrumented medians were
2.980 ms for outer/placement, 1.711 ms for the prefix, and 3.658 ms for
routing/post. Setup and allocations were excluded. Input restoration was
outside the timed call; required encoding movement was inside it.

These are implementation results, not a 40%-distance certificate. The
selected native source and measurement commands are in
[HD_AFFINE_NATIVE.md](HD_AFFINE_NATIVE.md). The source is
`spin/experiments/char2_field_spin/hd_affine_two_state.cpp`, with
`AF_GROUP=4 AF_NT=1 AF_TILE=16 AF_SPLIT_ROUTE=0`.

The source SHA256 remains

    bc71ce8fe05b6948a69456c6f279f58e58ae1d713a6fe68ddeb05c01ae3338e4

The three-product alternative has a validated local law but was never
implemented or timed. No new timing is implied by its product count.

## Proof results worth retaining

| Result | Scope and reference |
| --- | --- |
| Exact source interface | Prefix-diagonal averaging and component contraction; [source law](HD_AFFINE_PROOF.md). |
| First two classes | Entire permitted-singleton and pair terms below 2^-84.94 and 2^-41.11; [singletons](HD_COMPONENT_RANK.md), [pairs](HD_COMPONENT_PAIRS.md). |
| Dense tail | All 2826 classes a>=29943 below 2^-322.58; [dense analysis](HD_COMPONENT_DENSE.md). |
| First-moment obstruction | Actual-code projective first moment exceeds 2^441457 for group-four placement because a rare route contains a large bad subspace; same reference. |
| Capped event | Contained-component dimension obstruction below 2^-155.23 for all a>=18300 and all low-weight sets with at most 64 runs; [routing gate](HD_ROUTING_CONDITION.md). |
| Exact MDS filtering | Equality-aware component values and degree counters; [MDS interface](HD_MDS_COMPONENT.md). |
| Sparse refinements | Checked interval kernels and failed geometry/offset bounds; [geometry](HD_COMPONENT_GEOMETRY.md), [bounded refinements](HD_SPARSE_REFINEMENTS.md). |
| Cheaper post-inner | Three-product local law, not a complete proof or timing; [inner](HD_THREE_PRODUCT.md). |

The full 40%/40-bit certificate did not close. The first-moment obstruction
does not refute that distance target: its exhibited route event is much
rarer than 2^-40. It does rule out repairing the uncapped first moment
merely by improving source counts or Chernoff marks. The capped event
does not cover all low-weight vectors, especially cancellations and larger
families of output supports.

All 158 high-distance helper tests passed at the final checkpoint.
Independent audits checked the actual-code obstruction and directed bounds.
Native scalar, alias, alignment, guard, and sanitizer checks are recorded
in the measured checkpoint. No production profile adopted this design.

## Why this branch is parked

Its five coefficient products per output record and three accumulator
chains are a substantial cost. The outer and pre-prefix already consume
about 4.69 ms, leaving an implausibly small budget for the remaining phase
if those two costs stay unchanged. Removing one post product alone does
not establish a path below 5 ms.

The restart restores the implementation advantage of SPIN: amortize fixed
maps and state mixing across a packet. It keeps one routing stage and
does not introduce a separate pre-accumulator by default. The original
40%/40-bit and sub-5ms targets are unchanged.

## Verified source snapshot

The local source-only archive contains 578 files (7,116,424 uncompressed
bytes) from the characteristic-two workstream and native experiment tree.
It includes preceding field-code winners, so the detour can be inspected
with its local dependencies. New blockwise restart files are excluded.
Repository HEAD at capture was `3ce07a64936c00e4a9853f1081ab676681bb3382`;
the archive also preserves the relevant uncommitted/untracked source.

Archive, relative to the repository root:

    tmp/hd-accumulator-before-blockwise-20261006.tar.gz
    SHA256 756c949491dd1830e802764d762a42053d480e985598a9ae1b223a135c7826d9

The adjacent `.manifest` records each source path, byte count, and SHA256.
Its SHA256 is

    95fc8eeb47587af0ad7f6613a3a719b566993b473d3bcf5cf3ffe6cc03d32f19

The archive was extracted to a fresh ignored verification directory, and
all 578 file hashes and sizes matched. It contains no generated binaries,
raw timing logs, rendered paper files, or bulk proof dumps. It is a local
checkpoint, not a commit or a published artifact. External dependencies
and the measurement host are described in the original reproduction notes.
