# Cost of crossed mixing inside an outer superblock

The extra-stage proposal now has a checked isolated implementation.
At K=2^20, rate1/2, on128-bit records, one layer costs approximately
0.47ms beyond the unchanged H* encoder; two layers cost approximately
0.98ms. Neither new ensemble has a whole-code distance certificate.
The [exact local proof](HD_OUTER_PAIR_LAYERS.md) establishes support laws
and improvements for specific quiet-plane inputs, not35% distance.

## Construction and implementation

All variants retain sixteen RS[16,8] rows per superblock, H*, the same
global paired route, and the same64-coordinate inner with16-symbol state.
Independent nonzero D1 and D2 scalars remain. The new layers go between
H* and D2. Each pair uses a fresh independent nonzero alpha sampled once
at setup, then stays fixed for every message. Its map is
`(a,b) -> (a+alpha*b, a+2*alpha*b)` in the AES byte field.

The one-layer and two-layer wiring differ; two layers are not simply a
second pass after the preferred one-layer schedule. One layer uses row
differences8..15 over successive column pairs. Two layers use first row
difference1 and second differences(2,4,8,14) over column quartets.
The lower column is the first endpoint. The exact schedules, gauge
argument, and independence requirements are in the linked proof note.

`cross_pair_native/` is separate from the frozen experiment and production
code. Its blocked kernel reuses consumed RS scratch slots, retains
explicitly unrolled row operations, and applies D2 at the final64-byte
non-temporal stores. There is no new whole-array pass or runtime row
permutation. Alpha and2alpha are prepacked for independent GFNI products.
Each layer adds2MiB of hot keys at K20 and intended32MiB of local scratch
write/read traffic. Actual register spilling was not separately counted.

The original per-column mixing has39 GFNI products and66 XORs including
D1, H*, and D2. Each added layer contributes16 products and16 XORs.
The unchanged RS work is another17 products per column. These are circuit
counts; they do not predict cycle cost by themselves.

## Matched serial measurements

Measurements ran on Peach CPU2, using GCC15, the retained release flags,
and precomputed setup. Each full-encoder entry below is the median of501
calls. Three rounds rotate order: baseline/one/two, two/baseline/one,
one/two/baseline. Both shared benchmark locks were held. No two of our
benchmarks or compiler jobs ran concurrently. Other host activity can
still cause noise, as in earlier measurements.

| Variant | Round1 ms | Round2 ms | Round3 ms | Median of medians ms | Outer-only ms |
|---|---:|---:|---:|---:|---:|
| Unchanged H* |3.373608|3.351397|3.436977|3.373608|1.629893|
| H* plus one crossed layer |3.929185|3.817126|3.847473|3.847473|2.101363|
| H* plus two crossed layers |4.351404|4.339312|4.381160|4.351404|2.530394|

Outer-only measurements are one501-call run per variant, after the full
rounds. The full increases are about14.0% and29.0%. This test does not
establish that these implementations are optimally tuned. It does show
that adding the layers to the existing mixer misses the3.4-ms target.

Raw logs and exact build/benchmark scripts are retained under
`cross_pair_native/`. Remote scratch was
`/tmp/hd-cross-pair-20261009-6HGxEE`; the older scratch was not modified.
WinSCP could not accept the configured host-key fingerprint; the transfer
used the existing binary-safe SSH method with the same pinned SHA256 key.

## Correctness and scope

All three variants pass `--check` and `--large-check`. The new tests cover
8192 binary basis directions,16 random halves, both RS cosets, and the
independent scalar pair oracle. Existing tests retain30720 packet basis
directions, wrapping pair loads, K8192 and K24576 cases, aligned and
unaligned buffers, aliasing, and full K20 scalar-reference agreement.
The baseline full checksum remains6721820652671321300; one layer gives
2389812500930817425 and two layers2597169192294202319.

Independent source review checked pair-key indexing and packing, the
non-clobbering scratch schedule, final D2 exactly once, and the normalized
sampler's ideal distribution. Symbol inspection found no outlined
crossBlock, crossPair, or crossFinal helper. Non-temporal stores remain.
This pass did not run the optional sanitizer build.

The local exact-law suite has12 passing tests, including exhaustiveGF(4)
checks over every input and fresh-key tuple, and all140 affine four-plane
profiles. Those checks validate the local algebra, not the full distance.
The certified fallback manifest is unchanged:
`accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.

Source SHA256 values:

- `HdGroupedOuter.h`: `8df051b4da797597b245711be62a3adec1d279582b6bd58e44fd4439304966d5`.
- `hd_grouped_forward.cpp`: `e8b7eb9257609754b2a7b1195ed9a52bb7756eb45d7e8692b793a3d3b244122c`.
- `build.sh`: `a5f8ddc578e7f6f28bf0a0b2c18aa3f354f606f7dbbba5ae7774b5592d077833`.

Binary SHA256 values are in `cross_pair_native/hashes.log`.

## Next structural test

Prefer one layer for the next proof/cost comparison. To recover the time
budget, test a cheaper first mixer rather than adding another full layer.
One precise candidate is H0=I+11^T:31 XORs and no field products per
sixteen-coordinate column. It still spreads a singleton to15 positions,
but its rank-one parity condition leaves a much larger quiet subspace.
It is not yet implemented or certified. Reject or accept its two-row
quiet contribution under the crossed-layer bound before native work.

An exact intermediate proof gate is also available for H*. For a fixed
labelled input row set A, obtain the full65536-mask output count vector
by zero-set ranks and integer Mobius inversion. At a scalar output marker
t, the one-pair factor matrix is
`[[1,t^2],[t^2,(2t+253t^2)/255]]`. Its16-fold tensor contraction evaluates
the two-column count without enumerating all pairs of masks. Split each
input vector into exact zero and nonzero terms before contraction; retain
the resulting degree0/1/2 original-column activity polynomial. Multiply
the eight actual wiring blocks with the same A and retain degree>=9.
The scalar regional bound and source/RS factors must be supplied once
with their correct scope. This is a proposed gate, not a computed result.

No paper, production profile, commit, or push changed in this pass.
