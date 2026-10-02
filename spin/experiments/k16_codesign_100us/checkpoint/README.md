# Frozen K16 research checkpoint

Checkpoint `k16-rs16-paired15-v1` preserves **mode 52**, the fastest fully
certified K16 construction from this search. It is not a library promotion.
Do not replace this baseline with the faster but incompletely certified RS8
XOR candidate, or with the proposed byte-native inner.

The binary code has dimension 65,536 and length 131,072. Each outer group
maps 128 bits to 256 bits using four GF16 RS[16,8] rows and independent
nonzero field16 symbol randomizers. Four-coordinate packets reach 64 regions.
The inner uses 64-coordinate steps and the fixed paired 15-bit map, with
independent GL15 state updates. State starts at zero and is not flushed.

The complete bound gives minimum distance at least 13,108 except with
setup-failure probability below 2^-62.04117326862128. The independent
precomputed-transpose holdout measures **98.950 us** on 128-bit elements.
Setup and allocation are excluded. These are separate mathematical and
performance statements about the construction and implementation.

## Exact implementation

The selected setup is `customizePaired15Shuffle` followed by
`customizeNativeField16`. The allocation-free encoding stages are
`reverseRoutePaired15Fold` and `fieldLoopSharedParity`. Their loop structure,
coordinate permutation, dummy SIMD slot, cached stores, and in-place
contract are part of this checkpoint. Future adaptations should use new
files or a new checkpoint rather than silently changing these sources.

`source_manifest.toml` pins 220 source files: the experiment, independent
oracles and reviews, proof calculation dependencies, and private headers.
Hashes normalize CRLF to LF so Git's Windows checkout convention does not
invalidate an otherwise identical source tree. Raw receipts and timing
samples remain excluded from version control. The original receipt and
measured executable identities are recorded in the manifest.

The experiment links an external `spin::spin` package only for allocation,
capability detection, and the retained comparison driver. The winning code
does not call the production transpose implementation. Its pinned private
headers are unchanged by the pending production work. Compiler and linked
library binaries are not reproduced by the source-identity check.

Run from the repository root:

```text
python -B spin/experiments/k16_codesign_100us/checkpoint/verify.py
python -B spin/experiments/k16_codesign_100us/checkpoint/verify.py --receipt research/workstreams/k16_codesign_100us/proof/paired-s15-drop10-whole-p256.json
python -B -m unittest discover -s research/workstreams/k16_codesign_100us/review -p "test_paired15*.py" -v
```

The second command is optional and requires the original local receipt.
It authenticates that archive, not a new proof. The third command audits
all saved mathematical source hashes, exact union arithmetic, interfaces,
and a fresh direct-physical q1 calculation. See the
[construction record](../../../../research/workstreams/k16_codesign_100us/proof/PAIRED_S15_CONSTRUCTION.md)
for the full replay, and the [experiment README](../README.md) for compilation
and serial measurement commands.

## Larger sizes: probe first, then adapt

A bounded scaling probe ran the unchanged mode 52 executable on Peach,
CPU15, with normal pages, setup excluded, and no diagnostic stage clocks.
It used seeds 1 and 17, increasing and decreasing size order, and 101 calls
per process. Every process first checked scalar equality, adjoints, alignment,
guards, and in-place behavior. The executable still matches the frozen hash.

| K | Median of four process medians | Complete bound for this configuration |
|---|---:|---|
| 2^16 | 99.205 us | 62.04117 bits at 10% |
| 2^18 | 413.813 us | Not yet available |
| 2^20 | 5,103.994 us | Not yet available |

The K16 measurement above is a short scaling control, not a replacement for
the longer independent holdout. No timing or certificate is interpolated.

K18 is the next sensible adaptation: the measured kernel scales smoothly
there. Keep the outer/routing design initially, but screen the actual paired
inner at the new length before choosing a state size. More steps change the
occupancy analysis; the K16 certificate does not extend automatically.

At K20, retain the existing 256-to-512-bit RS/s20 construction, whose complete
bound has 68.10376 bits and whose prior implementation measured about 3.26 ms.
That prior timing is not a same-session A/B comparison. The unchanged K16
kernel's cached-store path is not a competitive replacement in this probe.
Any K20 adaptation needs both a large-working-set implementation and a fresh
complete proof. Earlier small-outer length experiments had strong q1/q2
bounds but failed the middle-occupancy bound; sparse evidence alone is not
sufficient. See the [length-study index](../../../../research/workstreams/k16_design/README.md).
