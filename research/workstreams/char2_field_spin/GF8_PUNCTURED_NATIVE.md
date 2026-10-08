# Three-quarter-rate GF256: validated research encoder

The seven-of-eight construction now has both a complete K20 distance
certificate and a matching native transposed encoder. It improves the
proved rate frontier to 3/4 at roughly the cost of the preceding lower-rate
field encoder. It is not faster than the production binary code.

The [whole-code certificate](GF8_PUNCTURED_WHOLE_CERTIFICATE.md) gives 10%
symbol distance and at least 116 bits of setup-failure margin at K20.
The [construction specification](GF8_PUNCTURED_PACKET.md) fixes the maps,
random variables, state boundaries, and retained output coordinates.

## Implementation and interface

The isolated research driver is `spin_char2_punctured8`, built from
`spin/experiments/char2_field_spin`. No installed library interface or
production profile changes. Setup, allocation, and page preparation are
outside encoding. The timed operation consumes N retained output records
and computes the field transpose into K source records; every record holds
sixteen parallel GF256 symbols. No general GF(2^128) multiplication is used.

The setup plan stores one regional route, its inverse, and independent
input/state diagonal arrays. The benchmark uses reproducible pseudorandom
fixtures; the mathematical probability bound refers to independent uniform
permutations and nonzero scalars. The measured circuit and coordinate
selection match the certified map, but timing a fixture does not certify
that particular sampled code.

The native schedule uses the exact factorization F=242 Aeff C, with
Aeff=2A. Two eight-point XOR-circulant maps replace three dense maps. A
four-packet batch packs the seven retained coordinates directly into seven
ZMM registers; the missing coordinate is filled with zero in the transpose.
The state recurrence runs backward through the four chronological packets.
The kernel uses GFNI byte products, fixed unrolled maps, and no allocation
or runtime callback inside encoding.

Source records use four-row-interleaved additive-polynomial coefficients:
`sourceIndex(row,j)=(row/4)*336+4*j+row%4`.
Retained outputs use `outputIndex(step,j)=28*(step/4)+4*j+step%4` for j<7.
This fixed output-coordinate ordering is part of the interface. No hidden
conversion from packet-major payloads is included in the timings.

Natural requested sizes are shortened from multiples of 84*64. Scratch
holds 128 columns, plus one regional temporary. Exact input/output alias is
supported; precisely K records are replaced and the input suffix is kept.
Column-stride padding is an exact-map experimental control, not a change
to setup, rate, or proof geometry.
The research plan exposes its fields: a custom `planeStride` must be at
least `groups` and divisible by four, with scratch resized accordingly.

## Correctness checks

The independent scalar oracle constructs the full16 Cauchy matrix and
normalizes its state block with generic Gauss-Jordan elimination. It does
not call the factored native maps or generated outer FFT. Its literal
dense transpose and direct additive-polynomial outer evaluate the complete
encoder independently.

Release checks pass for:

- All 36 four-packet input/state coordinates with all 256 byte values,
  plus random cases; both emitted records and incoming dual state agree.
- Full scalar encoders at requested K=1,83,84,85,335,336,337,5376,5377,6000.
- Shortened forward/transpose bilinear identities, including all16 lanes.
- Tight/padded scratch strides, both outer schedules, scatter and regional
  routing, 16/32-bit route tables, aliasing, guards, and unaligned offsets.
- K16/K18/K20 schedule agreement and exact-alias suffix preservation.

The initial and final inlined/stride-control builds pass AddressSanitizer
and UndefinedBehaviorSanitizer, including the final small and large checks.
The default outer schedule is gathered loads, matching the reported winner.
The root's final proof/diagnostic run passes all28 tests, and its full
certificate replay independently reproduces the116-bit aggregate bound.

## Serial measurements

Measured on Peach's Ryzen7950X, GCC15.2, Release, CPU2, boost disabled.
GFNI/AVX512 kernels use 128-bit application records. All calls are serial.
Competing encoders share payload and scratch allocations and equal requested
K. Setup and allocation are excluded. A fresh payload is copied outside
timing before each batch of five calls; candidate order rotates. Checksums
agree between exact-map schedules. Huge pages are preferred for buffers at
least16MiB, as in the previous controls.

The final inlined kernel, with compact columns and gathered outer loads:

| Requested K | Retained N | Requested rate | Full transpose | Matched cubic control |
|---|---:|---:|---:|---:|
| 2^16 | 93,184 | 64/91 | 0.240ms | 0.233ms |
| 2^18 | 351,232 | 256/343 | 0.938ms | 0.930ms |
| 2^20 | 1,404,928 | 256/343 | 5.571--5.816ms | 5.824--5.907ms |

K16/K18 rows are performance and correctness results, not new distance
certificates. The complete certificate in this checkpoint is for K20.
The small K16 request has visibly more padding; do not describe its actual
requested rate as 3/4. The padded construction has rate3/4 at every size.

The K20 range reports three separate 101-call launch medians:
5.816077, 5.571260, 5.815286ms. Their matched lower-rate cubic medians are
5.880287, 5.907137, 5.824202ms. A separate phased launch reached5.366ms;
do not substitute that favorable sample for the repeated full-call range.
In that phased launch, inner plus route was3.277ms and the outer was2.084ms.
The preceding cubic control emits1,622,016 records, so the new output is
13.38% shorter at the same useful K. Against binary rate1/2, it is33.01%
shorter. These are output-length savings, not measured downstream speedups.

### Exact-map experiments

The first non-inlined prototype measured about0.250/0.967/7.210ms at
K16/K18/K20. Its packet hook was emitted as a separate function, forcing
an intermediate output handoff. Forced inlining permits the compiler to
combine that handoff with the diagonal/store sequence. This reduces local
overhead without changing the code. The `SPIN_P8_FORCE_INLINE` CMake option
retains the non-forced control; final measurements above enable it.

An apparent K20 column-stride conflict motivated one-cache-line padding.
It did not help: three final K20 padded-gather medians were6.242,6.247,6.237ms,
versus5.571--5.816ms with compact columns. Padding is therefore not the
default. It gave small K16/K18 wins in one launch (.234/.917ms), but those
are not promoted as a universal layout rule. K20 process/allocation
variation also affected the old control. We do not attribute the entire
first-to-final timing difference solely to inlining.

## Reproduce

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release -DSPIN_P8_FORCE_INLINE=ON
cmake --build build/field --target spin_char2_punctured8
build/field/spin_char2_punctured8 --check
build/field/spin_char2_punctured8 --large-check
taskset -c 2 build/field/spin_char2_punctured8 1048576 101 full
taskset -c 2 build/field/spin_char2_punctured8 1048576 101 phases
```

The next construction experiment should spend some of the sparse margin
on a higher outer dimension while keeping seven retained outputs. That
tests rates above3/4 with the same kernel geometry. First replay the sparse
and middle gates; the116-bit figure alone does not justify such a change.
Alternatively target the measured inner/routing3.3ms without changing the
certified map. Preserve this checkpoint before pursuing either branch.
