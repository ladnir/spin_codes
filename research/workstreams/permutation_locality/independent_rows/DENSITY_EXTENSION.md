# Extending feedback-density bounds with upward counts

The baseline density refinement stops at six occupied windows because
its exact convolution requires N*D < 2^63. Here N = 2^19 is the
state count and D is a shape's counting denominator. Exact feedback
counts can now be reconstructed through ten windows, but their larger
denominators can violate that convolution guard.

Fix a shape and let n(s) be its exact feedback count at state s,
with sum_s n(s)=D. Choose an integer r>=0 and define

    n'(s) = ceil(n(s)/2^r),    D' = sum_s n'(s).

The key inequality is pointwise:

    n(s)/D <= (2^r D'/D) * n'(s)/D'.

Consequently, every positive convolution bound for the normalized
counting measure n'/D' remains valid for n/D after multiplication by
2^r D'/D. This factor must not be dropped. The rounded measure is
only a proof upper bound; it is not the encoder's feedback distribution.

The implementation uses r=max(0,bit_length(D)-40). It checks each
integer addition, the normalization, and N*D' < 2^63. Pointwise rounding
also gives

    1 <= 2^r D'/D < 1 + N*2^r/D.

Thus this representation can keep the added mass small even when the
original denominator is large. The existing exact Walsh convolution
routine then computes translated expansion-class counts for n'. It
rounds exponential factors upward and returns positive bounds for the
mature-density and uniform-class contributions to the density coordinate.
The adapter scales those bounds by 2^r D'/D using outward arithmetic.

For each local occupancy, the adapter applies the all-one penalty to
each shape before taking a maximum over all shapes. Only the scalar
C-to-C and U_v-to-C bounds change. No complement involving a rounded
probability is used; zero counts and exact class counts remain separate.

## Implementation and checks

`candidates/density_extend.py` uses the guarded limb inversion from
[SPECTRAL_FEEDBACK.md](SPECTRAL_FEEDBACK.md) to obtain n. It processes
all shapes at occupancies seven through the requested cutoff, at most
ten. It does not change the old density verifier or its cache format.

Tests check pointwise integer domination, including large denominators,
and compare the scaled bounds with direct small-state convolutions.
The complete-map census passes all 791 shapes at occupancies seven
through ten. Its largest mass factor is exactly
15335685049/15335681400, less than 1.00000024. At each tested tilt
(.068, .072, .076), it tightens all 24 scalar entries in the requested
range. This local result alone gives no full-code certificate; the
selected-point transfer replay is a separate calculation.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 80 --supports 192 200 208 --tilts .068 .072 .076 --exact-feedback 10 --joint-four --density-through 10
```

The 192-bit outward replay of this command gives:

| Union support per active group | Selected tilt | Log2 upper, rounded upward |
|---:|---:|---:|
| 192 | .068 | +136.909130 |
| 200 | .072 | +158.598967 |
| 208 | .072 | +9.147664 |

These improve the preceding exact-feedback/joint-four bounds by about
5.04, 6.32, and 6.34 bits, respectively. All three remain positive and
therefore vacuous as probability bounds. They do not increase the
certified occupancy range. This command only bounds the listed
homogeneous support events; even a successful selected event would need
a complete support cover and aggregation with all other occupancies.

The command retains the fixed mass-based replacement at occupancies six
through twelve. In that range it replaces C-to-C completely, so any
improvement to that old scalar coefficient is not used. Improvements
to U_v-to-C remain available. Before judging the extension's full value,
retune the fixed column mixtures with the same density option:

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_optimize.py --groups 80 --support 200 --tilt .072 --maximum 16 --iterations 30 --passes 2 --exact-feedback 10 --joint-four --density-through 10 --epoch-cache tmp/independent-row-local-families
```

This optimizer proposes only rational witnesses and requires their
outward replay. The completed support-200 run gives log2 upper
+158.187720, rounded upward. It selects full replacement of both columns
at local occupancies six through sixteen. Its gain over the fixed mixture
is only 0.412 bits, and the bound remains vacuous.

The optional `--epoch-cache` stores the refined unsplit local matrices
and the three coupled mass coefficients for each local occupancy.
Both experimental drivers use `candidates/local_family.py` to construct
these objects. Column replacement still occurs afterward, through the
same outward `blend` routine. In particular, tuning cannot discard a
coupled coefficient while retaining the benefit of its alternative.

Cache keys bind the generator sources, selected maps, arithmetic precision,
output tilt, all-one penalty, and every refinement option. They also bind
the largest local occupancy available for mass replacement. Occupancy q,
support vectors, and proposed mixture fractions are not dependencies of
these local bounds. Each stored number is an exact dyadic upper endpoint.
Records have checksums and are published atomically; a source change during
construction prevents writing the memo. The memo is not an independent
certificate. Omit the option to regenerate all censuses.

Eight cache tests check exact replay, option and source invalidation,
corruption detection, precision guards, partial cache hits, and a source
change during construction. The combined test suite passes 139 tests
(94 verifier tests and 45 experimental tests). These tests add no
certified occupancies.
