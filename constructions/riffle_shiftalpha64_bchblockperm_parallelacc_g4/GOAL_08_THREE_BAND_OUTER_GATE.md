# Goal 08: three-band outer gate

## Target

Goal 07 used one bound for every nonzero BCH weight. That bound mixed common
weight-64 words with rare endpoint words. Goal 08 tests whether three weight
bands preserve enough information without returning to exact profiles.

The bands are:

| band | BCH weights |
|:---|:---|
| light | 22 through 50 |
| normal | 52 through 76 |
| heavy | 78 through 128 |

The normal band contains about 97.94% of the nonzero field values. The light
and heavy bands each contain about 1.03%.

## How the two outer checks enter

Fix an outer support of size (s). Also fix the band of every occupied
coordinate. Choose two coordinates. The two MDS check columns at those
coordinates are invertible.

The remaining (s-2) coordinates are free. Their weighted sums come from the
exact BCH spectrum inside each band. The two chosen coordinates are uniquely
determined. The calculation keeps a term only when each determined value lies
in its declared band. It upper-bounds that term by the largest coefficient
weight in the band.

This gives a rigorous bound for each band pattern. Summing all patterns covers
the complete occupation shell.

## Result

The three-band split improves both tested shells:

| occupied symbols | best scalar bound | three-band bound | improvement |
|---:|---:|---:|---:|
| 3 | (2^{183.223}) | (2^{160.172}) | 23.051 bits |
| 4 | (2^{272.163}) | (2^{231.142}) | 41.021 bits |

The bounds remain trivial. The exact message counts are (2^{103.415}) and
(2^{179.415}). Therefore the trivial statement that every message is bad is
stronger than the new probability calculation.

The dominant patterns are:

| occupied symbols | light | normal | heavy | pattern bound |
|---:|---:|---:|---:|---:|
| 3 | 1 | 0 | 2 | (2^{159.773}) |
| 3 | 2 | 0 | 1 | (2^{158.124}) |
| 4 | 2 | 0 | 2 | (2^{231.142}) |

The all-normal patterns are much smaller but do not close:

| occupied symbols | all-normal pattern bound |
|---:|---:|
| 3 | (2^{18.588}) |
| 4 | (2^{63.043}) |

## Interpretation

The band split recovered a material number of bits. The remaining loss no
longer comes from mixing the central BCH spectrum with the endpoints in one
scalar moment.

The current calculation assumes that free light and heavy values often solve
to the worst allowed dependent values. It does not measure how often the two
field equations permit that outcome. The mixed light-heavy patterns dominate
because of this assumption.

Finer weight bands would not address that assumption. The next useful object
is a compatibility bound. For a fixed pair of outer coefficients, it should
count how often values from two declared bands solve to values in two other
declared bands.

For occupation three, the required object is a multiplier intersection among
three bands. For occupation four, it is a two-input, two-output band table.
The exact weight-22 and weight-24 multiplier scans from Goals 01 and 02 provide
the first authenticated rows of such a table.

## Reproduction

Run the shells sequentially:

```powershell
python scripts/analyze_riffle_shiftalpha64_threeband.py --occupation 3 --optimizer-maxiter 50
python scripts/analyze_riffle_shiftalpha64_threeband.py --occupation 4 --optimizer-maxiter 50
```

The analysis script has SHA-256
`1C9876477FE0A6028D87EE2A77910D175EA0247809627202D69C60929268ED7E`.
The summary receipt is `receipts/goal08_three_band_outer_gate.json`.

## Next goal

Build a small compatibility table for the dominant mixed patterns. Start with
occupation three. Test whether one light value and two heavy values can align
under the actual outer coefficient ratios at the endpoint weights that
dominate the optimized bound. Then extend the test to the two-light,
two-heavy pattern at occupation four.

