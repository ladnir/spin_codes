# Four-window cancellation and output weight

The baseline uses exact joint cancellation/output counts through three
occupied packet windows. At four windows it retains less information.
This candidate extends the same local inequality to four windows without
changing the encoder, setup distribution, or state coordinates.

Fix a packet-weight multiset h of length four. Let X be uniform over
all inputs having that multiset in four distinct windows. Write B for
feedback and E for expansion, using the authenticated Map128S19 maps.
On a lazy return to zero, the entering state equals BX and the emitted
word is E(BX)+X. The census records

    J_h(v,w) = #{X : wt(E(BX))=v, wt(E(BX)+X)=w}.

It also records the fresh-state counts H_(h,a)(w) defined in
[the existing joint-return argument](../JOINT_CANCELLATION.md).
Fresh-state membership is computed using all 32 windows, even in the
small-domain tests. The production map's zero-or-one membership property
is checked explicitly.

For the full domain, every shape has denominator

    D_h = binomial(32,4) 4!/product_b n_b! * product_i binomial(4,h_i),

where n_b counts occurrences of weight b. There are 35 shapes and
binomial(32,4)*15^4 = 1820475000 input patterns in total. These totals
bound every integer histogram accumulator.

The implementation precomputes pairs of windows. It combines pair
feedbacks and emitted words with vectorized XOR and population counts.
In particular, E(BX)+X is linear in X, so the emitted word can be formed
by XORing the precomputed pair words. The calculation includes every
window subset and lane-mask assignment; it does not sample them.

For each shape, the old joint-return formulas give coefficients for
C-to-Z, F-to-Z, and the five U_v-to-Z transitions. They retain the
refresh contribution and the factor 1/4 for two lazy transvections.
Input penalties are applied before maximizing over the complete shape
set. The resulting bounds are intersected with existing scalar bounds
before any coupled mass-based column replacement.

## Checks and integration

`candidates/joint_four.py` checks exact shape totals and output-weight
triangle inequalities. It compares selected window subsets with a
separate scalar enumeration. The small-domain tests compare all joint
and fresh counts for one through four occupied windows. Integration
also checks every expansion-class marginal against the independent
Walsh feedback census.

The old shared-route files are unchanged. The experimental drivers
select this refinement with `--joint-four`; it is not loaded from the
production operator cache. The first full-map replay passes all 35
shape totals, the selected scalar comparisons, and all 35 independent
Walsh expansion-class marginal checks. At each tested tilt it tightens
all seven zero-return coefficients associated with four occupied windows.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 80 --supports 192 200 208 --tilts .068 .072 .076 --exact-feedback 10 --joint-four
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 64 --full-cover --max-splits 300 --tilts .054 .056 .058 .06 --exact-feedback 10 --joint-four
```

The first command addresses selected support vectors only. The second
requires a complete support cover and outward replay before any new
occupancy claim can be made.

The occupancy-64 run has now ended at its 300-split limit, with 8784
leaves and binary64 log2 upper -43.1317304404339. It misses that run's
-54 replay gate, so no outward occupancy certificate was produced.
This numerical improvement does not extend the certified occupancy range.

The first command passes outward replay at 192-bit precision. Its log2
uppers, rounded upward, are:

| Union support per group | Selected tilt | Outward log2 upper |
|---|---:|---:|
| 192 | .068 | +141.953270 |
| 200 | .072 | +164.918511 |
| 208 | .072 | +15.492140 |

These gain about 39.46, 45.12, and 44.88 bits over the matching
exact-feedback results. All remain vacuous probability bounds. They do
not identify low-weight codewords, nor do they cover heterogeneous
support vectors.

The driver can also reuse its operators across an explicit occupancy
batch. Each member still has an independent complete cover and outward
replay. Failed members remain explicitly unverified; the batch aggregate
only sums successful members. A batch from 59 through 63 is running:

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --occupancies 59 60 61 62 63 --full-cover --max-splits 200 --tilts .052 .054 .056 .058 .06 --exact-feedback 10 --joint-four
```
