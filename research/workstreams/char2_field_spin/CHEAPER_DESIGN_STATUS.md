# Cheaper designs: what passed and what failed

The accepted target is higher rate at 10% symbol distance, or substantially
higher distance at half rate, with competitive total encoding cost.
Half-rate/10% already follows by scalar extension of binary SPIN.
The [GRS128/87 checkpoint](FIELD128_K20_CERTIFICATE.md) supplies a complete
higher-rate proof and a [measured K20 improvement](FIELD128_NATIVE.md).
This note indexes the surrounding experiments without treating every
failed upper bound as a bad-code witness.

## Exact-map and state-coordinate experiments

The residue-coordinate FFT inner computes exactly the certified map.
It reduces GFNI512 products per packet from 70 to 34, but increases shuffles
and dependencies. Release and sanitizer checks pass; serial timing regresses
full encoding by about 8% at K16/K18 and 5% at K20. Keep the dense inner.

A later exact representation change moves its butterflies into 32-bit
lanes and leaves persistent state in that layout. This lane-local FFT also
passes independent native and sanitizer checks, but remains about5% slower
than the dense inner at K16/K18 and loses at K20. Compact route indices are
a separate modest win in direct comparisons; full K20 timings remain
process-sensitive. See [the complete latest measurements](FIELD128_NATIVE.md).
Neither FFT failure changes the certified ensemble.

[Cheaper state coordinates](CHEAPER_INNER_MAP.md) change the diagonal
refresh distribution. Two-input returns can have rank 5 instead of 8;
four selected output zeros can cost one random ratio rather than four.
Those are exact failures of the existing proof properties. They do not
by themselves establish low whole-code distance. A new proof must use
the weaker actual ranks instead of transferring the old certificate.

Restricting the residue map to six state coordinates preserves its
superregularity and zero feedback product. The generic joint-return bound
was independently audited for this restriction. However, the complete
code is not certified. For GRS[128,86], width 16, K20 and a 10% target,
the joint middle-weight screen gives:

| Input density | Sampled margin (bits) |
|---:|---:|
| .03 | 2089.09 |
| .06 | 2351.42 |
| .10 | -18314.37 |
| .14 | -1091.95 |
| .20 | 190997.60 |
| .24 | 394718.97 |

These are floating bounds on individual classes, not certificates or
failure lower bounds. The middle gap outweighs the positive sparse
samples. Six regression tests pass. Replay the sampled screen with:

```sh
python -B research/workstreams/char2_field_spin/theory/symmetric_density_screen.py --outer 128 86 --state 6 --joint --density .03 .06 .1 .14 .2 .24 --iterations 55
```

## Cheaper outers and mixing

[The smaller-outer screen](JOINT_OUTER_SCREEN.md) selected GRS[128,87].
GRS[128,88] has only about 1.27 bits near 240 active rows in the tested bound;
GRS[64,44] has negative sampled multirow bounds. These results alone do
not disprove their distance, but do not justify native implementation.

[The cheap-outer analysis](CHEAPER_OUTER_ENSEMBLE.md) gives stronger
negative results for specific replacements:

- Direct sums of small GRS blocks leave component-supported subcodes that
  fit bad regional windows, independently of the inner coefficients.
- The one-product butterfly has an all-equal-input invariant and retains
  such a subcode.
- For the fixed binary map `I+J` with independent column diagonals shared
  across outer rows, an exact lower bound exceeds 2^-40. This rejects that
  shared-setup ensemble at the desired margin. It does not reject
  independently resampled per-row mixing or all structured mixers.

Eight exact tests pass for this obstruction study. More expensive
two-shear switches have no complete compatible whole-code proof yet.

## Higher-distance track

The new variable-distance screen uses the exact forbidden integer weight
`ceil(delta*N)-1`, not an interpolated distance label. At half rate:

| Outer | Target distance | Active rows | Sampled margin (bits) |
|---|---:|---:|---:|
| GRS[128,64] |30%|1|185.04|
| GRS[128,64] |30%|64|-2846.19|
| GRS[256,128] |30%|64|-5013.75|
| GRS[128,64] |20%|256|-645.27|

The stronger joint middle-weight screen also fails tested 30% classes.
These remain failed upper bounds, not impossibility results for the
regional construction. They identify the multirow/middle-weight proof
problem rather than a need for more singleton margin. Five tests cover
threshold handling, reference agreement, and monotonicity.

```sh
python -B research/workstreams/char2_field_spin/theory/higher_distance_screen.py --outer 128 64 --distance 3/10 --active 1 4 16 64 --iterations 25
```

Do not resume half-rate/10% certification as a substitute for this track.
The useful next step is a construction or rank analysis that improves
the failing multirow classes without making the already costly inner heavier.
