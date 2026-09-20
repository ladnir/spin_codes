# Binary degree-six outer-support certificate

This note records the diagnostics and final certificate for binary
two-sided regular EC with left/right degrees `6/3`.  The early diagnostic
values explain how the parameters were selected.  The final memory-79 values
under "Certified result" use outward rounding.

## Candidate and obstruction

The parameters are

```text
k = 1,048,575
n = 2,097,150
left/right degrees = 6/3
cutoff = 230,729
```

At convolution memory 60, the exact regional transfer gives a weight-64
first-moment term of `2^1.9346367`.  Expander collisions do not cause this
term.  One region retains parity weight 64 with probability
`0.9961618556`, and all six regions do so with probability `0.9771909761`.

The generic positive coefficient saddle is much looser.  At support 64 it
gives about `2^181.97`.  Giving each region its own input marker improves the
bound by only 4.4 bits.  Giving each start/end-state entry its own marker
improves it to about `2^127.32`, which is still far above the exact term.

The inactive-to-inactive regional entry is the largest source of this slack.
Replacing only that entry by its exact coefficient lowers the full bound to
about `2^49.91`.  Replacing the inactive row and column lowers it to
`2^12.86`.  The coefficient distribution for the inactive bridge is
multimodal: paths with different numbers of convolution activations prefer
different input markers.  A single local-limit correction is therefore not
appropriate.

## Memory 80

Memory 80 moves the delicate exact terms below the 20-bit target.  Near the
maximum, fixed near-optimal markers give

```text
support 104   -33.5483 bits
support 112   -31.7177 bits
support 120   -30.7159 bits
support 127   -30.4731 bits
support 128   -30.4845 bits
support 136   -30.9725 bits
support 144   -32.1351 bits
support 160   -36.3299 bits
support 192   -51.3221 bits
support 256  -102.7047 bits
```

The block diagnostic uses one output marker for 16 consecutive supports and
reuses their uniform-slice transfer matrices.  It gives

```text
supports       log2 of summed bound
80..95         -32.6318
96..111        -27.8391
112..127       -25.8801
128..143       -26.1967
144..159       -28.6731
160..175       -33.3347
176..191       -40.1018
combined       -24.7299
```

The next two blocks give `-45.0909` bits for supports `192..223` and
`-68.1582` bits for supports `224..255`.  The original binary64
implementation underflowed after this range.  Carrying one logarithmic scale
per conditional slice removes that implementation limit.  At the fixed
markers used by the scan, the scaled recurrence gives

```text
support 256  -102.7045 bits
support 288  -137.0860 bits
support 320  -176.3113 bits
```

All `81^2` regional entries remain nonzero in the scaled computation.  At
support 320, the smallest is about `2^-1301.54`.  This explains why the
unscaled matrix vanished even though its normalized structure was regular.

A single scaled slice family at output marker `0.992` bounds the sum over
supports `256..383` by `2^-55.1658`.  Its worst support is 256.  Thus the
scaled exact calculation overlaps the positive coefficient blocks.

These diagnostics selected memory 80 for the first certificate.  A later
threshold scan reduced the certified memory to 79.  The scaled exact and
coefficient regimes overlap, so the frozen partition has no support gap.

## Interval audit

A frozen partition covers every nonzero support.  It uses exact
regional blocks through support 383, positive-coefficient blocks on the two
outer ranges, the degree-three local point-mass bound in the central range,
and a separate all-one calculation.  A standard-library checker confirms
that these ranges are disjoint and cover supports `1..k`.

Three preliminary outward-rounded binary64 representations became vacuous
on the last exact block, supports `256..383`: a global relative factor with
one absolute error, a separate scale and residual mass for every row, and an
entrywise upper matrix that preserves structural zeros.  The reason was not
the observed exact term.  Some reachable transfer entries lie outside the
binary64 normal range; an absolute underflow allowance for those entries is
later multiplied by `binom(k,r)`.

The final verifier resolves this problem with exponent layers.  It allows
each reachable entry to use an independent power-of-two scale, performs only
positive binary64 matrix operations, and rounds every operation upward.  It
then combines the resulting matrix bounds with Arb interval arithmetic for
the scalar combinatorial factors.

Increasing convolution memory does not repair the generic coefficient
regime.  On supports `256..281`, its diagnostic stays near `2^141.07` for
memories 80, 100, 120, 160, 240, and 320.  This near invariance confirms that
the gap comes from conditioning six regions on their exact weights, rather
than from insufficient convolution memory.

There is no binary rate-one-half degree-`8/4` fallback.  The all-one message
has even incidence at every right vertex when the right degree is even, so it
is a deterministic kernel word.  The next admissible profile after `6/3` is
the certified `10/5` profile.

## Activation-count refinement

Introduce a marker for transitions that leave the inactive convolution
state.  For a fixed regional input weight, extracting this activation count
separates the modes of the inactive bridge.  In the memory-60, support-64
diagnostic, the exact inactive-to-inactive entry is `2^-279.39`.  A coarse
two-marker grid bounds it by `2^-269.96`, reducing the coefficient gap from
about 42 bits to 9.4 bits.

Representing the structural zeros correctly gives a complete regional bound,
but the resulting support-64 first-moment bound is `2^91.23`.  Most regional
entries retain about nine bits of coefficient slack, and six independent
regions compound that slack.  The activation marker therefore explains the
ordinary saddle's failure but does not repair the memory-60 candidate.

## Fixed-shell saddle audit at memory 80

The difficult memory-80 block begins at message support 256.  At output
marker `0.992`, the scaled exact calculation gives a first-moment term of
`2^-55.4077`.  The older unscaled calculation gives `2^-55.4492`.  Their
largest regional entries agree within `6.1e-10` bits.  The unscaled matrix
loses 3,430 of its 6,561 entries to binary64 underflow; retaining those
entries changes the final exponent by only `0.0415` bits.

Conditioning on the exact parity weight does not by itself yield a useful
analytic bound.  For support 256, an entrywise positive-coefficient bound on
each fixed-weight convolution slice gives `2^176.5213`, which is 231.93 bits
above the scaled exact term.  This loss occurs inside the convolution slice,
after the expander parity weight has already been fixed.

The worst individual slice entries reveal two path families.  For the
transfer from state 75 to the inactive state, the optimizing input marker is
`log(y) = -5.5`.  Under this marker, input weight one has tilted probability
about `2^-0.0235`, whereas input weight 256 has tilted probability about
`2^-60.8185`.  Thus input weight 256 lies in a valley between modes; the
usual one-marker coefficient estimate pays for the much larger low-weight
mode.

Three positive refinements were tested.

- Splitting by the number of transitions that leave the inactive state helps
  entries that start inactive, but not the active-to-inactive entries that
  dominate this obstruction.
- Splitting by whether the path ever waits in the inactive state does not
  improve the well-resolved bound.  Both modes use such a wait.
- Splitting by whether the first inactive wait occurs in the first half of a
  region lowers the complete support-256 bound from `2^176.5213` to
  `2^173.0702`.  Moving the cut to step 350 lowers it only to `2^176.3492`.

These temporal splits are valid, but their gains are too small relative to
the 231.93-bit gap.  The exponent-layer representation used by the final
verifier retains the tiny reachable entries without turning them into a
matrix-wide or row-wide absolute residual.

## Reproduction

For example, the first delicate block is evaluated by

```powershell
python scripts/binary_biregular_outer_scan.py `
  --support-start 80 --support-limit 95 --output-marker 0.99775366
```

Run different blocks serially.

## Certified result

The complete certificate is checked and evaluated with

```powershell
python scripts/check_binary_biregular_ec_d6_certificate.py `
  results/binary_biregular_ec_rate_half_d6_m79_gv.json
python scripts/binary_biregular_ec_d6_certificate.py verify `
  results/binary_biregular_ec_rate_half_d6_m79_gv.json
```

The partition contains 17 exact blocks covering supports `1..383`, 103
positive-coefficient blocks, 28 degree-three central blocks, and a separate
all-one term.  The outward-rounded component bounds are

```text
exact blocks          < 2.829174e-7
coefficient blocks    < 7.292690e-8
central blocks        < 1.397617e-13
all-one term          < 1e-315674
total                 < 3.558444e-7
security exponent     > 21.4222 bits
```

The largest exact block is supports `112..127`, bounded by
`1.229988e-7`.  The difficult final exact block, supports `256..383`, is
bounded by `1.955756e-16`.  The verifier therefore closes the 20-bit target
with more than 1.42 bits of margin.

## Memory threshold scan

Per-support optimization does not close memory 70: supports 80 and 96 have
bounds `2^-14.80` and `2^-6.33`.  At memory 75, the block covering supports
`80..95` is bounded only by `2^-17.66`.

Memory 77 remains inconclusive even after narrower blocks.  Supports
`104..111` contribute about `2^-20.36` before adding any other support.
At memory 78, four-support blocks over `112..127` sum to about `2^-20.74`.
Adding supports `96..111` and the certified outer contribution crosses the
20-bit target.  Thus neither value closes with the current partition and
positive output-tail bound.  This diagnostic does not prove that a different
enumerator decomposition cannot certify memory 78.

## Central and high supports

The first positive coefficient block, supports `384..422`, is bounded by
`2^-23.7090`.  Later low-side blocks decrease rapidly.  Ten-percent blocks
become too wide near `k/2`; this is endpoint slack, since the separate central
Fourier calculation gives `-62.3468` bits for one central support.

For right degree three, conditioning one group on its parity leaves a single
Bernoulli half-count.  The local point-mass argument used for right degree
five therefore applies directly.  One-percent endpoint blocks over supports
`400000..648575`, including the complementary half, sum to `2^-42.7021` in
the diagnostic.  This range overlaps both coefficient tails.  The first
high-side complement block, supports `648576..656477`, is already below
`2^-42901`.

The full outward-rounded sum, including all exact, coefficient, central, and
all-one blocks, gives more than `21.4222` failure bits.  Thus the certified
profile has more than 1.42 bits of margin beyond the 20-bit target.
