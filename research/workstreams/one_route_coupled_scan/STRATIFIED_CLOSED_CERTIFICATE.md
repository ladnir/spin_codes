# Fixed RS, one stratified route, and uniform recursive feedback

This checkpoint removes the random dense outer and the uniform global route
from the reference construction. It gives finite, all-message distance
certificates over GF(2^128), including a Toeplitz-batched inner. The inner
still uses uniform full-field coefficients. The certificate is independent
of the later [native measurements](NATIVE_BASELINE.md); those results do
not meet the practical speed target.

## Construction and certified parameters

Fix F=GF(2^128), K=2^20, and N=2K. Distance counts nonzero F-symbols.
For an even B dividing N, encode each B/2-symbol message block with the
same fixed RS[B,B/2] code. Its evaluation points are distinct elements of
the GF256 subfield. This code is MDS over F as well as over GF256.

Put M=N/B. Form B consecutive regions of M symbols. Region j contains
coordinate j from each outer block, in a uniform random order. Sample the
B orders independently. This is one global permutation, not B successive
permutations of the whole word.

The scalar inner uses memory L and independent uniform F-coefficients in
each feedback row. The batched inner instead groups b symbols and sets

    Y_t = X_t + T_t S_t.

Here S_t contains the preceding L output symbols, L is a multiple of b,
and T_t is a b-by-L Toeplitz matrix. Each matrix has b+L-1 independent
uniform F-diagonals; different batches use independent matrices. Start
with zero state, append Y_t to the state, and impose no terminal condition.
All routing and inner randomness is sampled once and fixed for every
message. The inner is invertible for every setup, so the rate is exactly 1/2.

The following bounds concern the existence of any bad nonzero message
after setup. They are not fixed-message estimates.

| Outer | Inner memory L | Batch b | Distance greater than | Setup-failure bound | Covered boxes |
|---|---:|---:|---:|---:|---:|
| RS[256,128] | 256 | scalar | 45% | < 2^-66 | 1240 |
| RS[256,128] | 256 | 64 | 45% | < 2^-72 | 2193 |
| RS[128,64] | 128 | 32 | 35% | < 2^-70 | 4440 |
| RS[64,32] | 128 | 32 | 35% | < 2^-70 | 9444 |
| RS[64,32] | 128 | 64 | 35% | < 2^-40 | 19935 |

The tested bad-weight thresholds are 943718 at 45% and 734003 at 35%.
Both equal floor(delta N). Excluding weights through that integer makes
the minimum relative distance strictly greater than delta.

The first four entries use the coarse output-volume bound below. The
last uses the sharper entropy bound. The final
route moments are exact rational numbers. Every box was verified with
outward-rounded logarithms at both 80 and 120 decimal digits. The
reported margins deliberately use the fixed per-term threshold, rather
than extracting every spare bit from the evaluated bounds.

## Why the fixed outer and route suffice

For a scalar coordinate set S, let t_i count the selected coordinates
from outer block i. The routed outer's shortened dimension is exactly

    D(S) = sum_i (t_i-B/2)_+.

At rate one half, the fixed total occupancy also gives

    D(S) = |S|-N/2+D(S^c).

The [fractional-moment proof](STRATIFIED_FRACTIONAL_BOUND.md) bounds the
moment of D(S) under the actual regional permutations. It compares
sampling fixed-cardinality subsets with independent Bernoulli indicators
using a product expansion with nonnegative coefficients. Convex averaging
then bounds the partial-region occupancy by a binomial variable.

For a set that meets at most r partial regions, define

    H_r(mu,theta)
      = max_f E exp(theta*(f+Bin(r,(mu-f)/r)-B/2)_+),

where mu=|S|/M and the integers f satisfy
0<=f<=B-r and 0<=mu-f<=r. Then

    E exp(theta D(S)) <= H_r(mu,theta)^M.

For a union of c intervals, r=min(2c,B) is safe. Increasing r is safe:
pad the original indicators with deterministic zeros or ones before
convex averaging. Increasing mu is also safe: increase indicator
probabilities, then apply the same envelope. These two monotonicities
permit uniform bounds on rectangles of interval sizes and counts.

When |S|>N/2, apply the complement identity before the moment bound. This
preserves the fixed total occupancy. Treating that total as another
binomial fluctuation gives a much weaker bound.

## The finite bound for Toeplitz batches

The [Toeplitz proof](TOEPLITZ_BATCH_BASELINE.md) shows that T_t S is jointly
uniform on F^b for every fixed nonzero S. A triangular selection of b
fresh diagonals proves this directly. Thus the distribution needed by
the inverse-output argument survives batching exactly.

Fix an output interval representation with c intervals and v batches in
their union. Write n=N/b and m=L/b. Let U be its scalar-coordinate union.
Delete the first batch from each interval to obtain A. Then

    |U|=bv,                 a:=|A|=b(v-c).

For a fixed output y, the inverse input is uniform on an affine translate
of F^A. Its intersection with the routed outer has at most q^D(A) points.
If D(U)=0, no nonzero output with this representation is possible. The
first activation input equals its nonzero output batch.

The first c-1 intervals force L(c-1) zero output symbols. Put

    f=bv-L(c-1),            e=floor(delta N).

The number W of eligible outputs is at most the q-ary ball volume on f
positions. The table uses

    log_q W <= w := min(f, e+f/128).

This follows from W<=q^f and W<=2^f q^e. For a fixed route, the probability
of any bad output in this representation is therefore at most

    1[D(U)>0] min(1, W q^(D(A)-a)).

For 0<theta<=ln q, use min(1,z)<=z^(theta/ln q), then average the route.
D(A)=0 does not imply an empty affine intersection, so retain its moment.
Complementing A when needed gives the following valid upper bound after
dropping the indicator:

    exp(theta*(w-min(a,N/2)))
      * H_min(2c,B)(min(a,N-a)/M,theta)^M.                 (1)

The deterministic zero test still uses U, not A. The implementation
discards a class only when bv/M+2c<=B/2, which implies D(U)=0.

The number of possible interval representations at (v,c) is at most

    binom(v-m(c-1)-1,c-1) * binom(n-v+c,c).               (2)

The range is m(c-1)+c<=v<=n and
1<=c<=floor((n+m)/(m+1)). Count (2) permits abutting intervals and includes
some impossible final clusters. It therefore does not exclude any output.
Multiply (1) by (2), then sum over this entire range. This proves a bound
on whole-code setup failure. No source class is omitted.

### Sharper optional output count

For 0<e<f<q, choose t=e/[q(f-e)] in the elementary generating-function
bound V_q(f,e)<=t^-e(1+qt)^f. Since t<=1,

    log_q W <= e + (f ln f-e ln e-(f-e)ln(f-e))/(128 ln 2).

Take the minimum with the coarse bound. At e=0 the logarithm is zero;
at e>=f use f. The bound is nondecreasing in f, so a rectangle may use its
largest f. This option changes only the proof, not the encoder or setup.

## Complete coverage and numerical verification

The proposal pass recursively partitions the entire integer (v,c) range.
It drops only points that violate m(c-1)+c<=v, and bisects all remaining
rectangles until a uniform upper bound passes. A selected grid is not
used as coverage. The scalar control uses the same construction in
scalar interval coordinates and the bound in the fractional-moment note.

Each box fixes a rational marker tau=exp(theta) from

    17/16, 9/8, 5/4, 3/2, 2, 4, 8.

Every feasible integer f in the definition of H is checked. The excess
H-1 is computed as an exact fraction. The directed replay uses enclosing
Decimal logarithms and preserves signs when multiplying negative
coefficients. Binomial counts use the upper entropy bound
n ln n-k ln k-(n-k)ln(n-k). Complement moments retain their zero-dimensional
outcomes. The optional entropy volume uses a downward-rounded denominator.

The first four rows verify that every contained pair contributes at most
2^-100. There are fewer than 2^34 scalar pairs, 2^28 pairs for
(L,b)=(256,64), and 2^30 pairs for (128,32). The final row uses a
2^-69 per-pair threshold and fewer than 2^29 pairs. Union bounding gives the
displayed margins. The two precisions agree on every required inequality;
agreement alone is not the justification for rounding validity.

The largest directed per-pair log2 upper bounds, rounded upward here, are
-101.316543, -132.878987, -200.229442, -100.270196, and -69.262801,
in table order.

Root also passed bounded exact combinatorial checks: 69952 regional
placements, 2766 moment comparisons, 1988 convex equalizations, 2992
partial-count enlargements, 3034 binary outputs, and 3862 small-box
component comparisons. Additional checks cover 2380 GF2 Toeplitz products,
2485 expanded batch-activity patterns, 416 batch-volume comparisons, and
2976 batch-box comparisons. These fixtures check the machinery; the full
directed covers establish the finite certificates. Independent analytic
and code reviews found no gap in the moment, inverse, or coverage steps.

## Reproduction and failed probes

Run serially from the repository root. Python with NumPy is required;
the scripts write no data files or bytecode caches with these commands.

```text
python -B research/workstreams/one_route_coupled_scan/theory/test_stratified_fractional.py --with-float
python -B research/workstreams/one_route_coupled_scan/theory/stratified_certificate.py
python -B research/workstreams/one_route_coupled_scan/theory/stratified_batch_certificate.py
python -B research/workstreams/one_route_coupled_scan/theory/stratified_batch_certificate.py --B 128 --L 128 --batch 32 --distance 35/100
python -B research/workstreams/one_route_coupled_scan/theory/stratified_batch_certificate.py --B 64 --L 128 --batch 32 --distance 35/100
python -B research/workstreams/one_route_coupled_scan/theory/stratified_batch_certificate.py --B 64 --L 128 --batch 64 --distance 35/100 --entropy-volume --term-bits 69
```

The proposal pass is capped at 120 seconds and 50000 rectangle checks.
Each directed replay is capped at 180 seconds. Caches are bounded and
remain in memory. Root ran all jobs serially with more than 60 GiB free.

Two smaller-parameter probes did not complete a 40-bit cover with this
envelope, even after sharpening the output count:

| B,L,b | Distance | Per-pair target | First rejected (v,c) | Floating log2 bound |
|---|---:|---:|---:|---:|
| 128,128,32 | 45% | -70 | (34261,1193) | -69.6463 |
| 64,64,16 | 35% | -72 | (46087,54) | -62.9247 |

These are first rejected points, not global maxima or low-distance
witnesses. A point exceeding the uniform per-pair target does not rule
out a tighter summation, more markers, or a better route envelope.
The failed probes do not weaken any completed certificate above.

## Remaining implementation question

At (L,b)=(256,64), Toeplitz setup stores 319/64 field coefficients per
output symbol, instead of 256 for scalar feedback. At (128,32), it stores
159/32 instead of 128. The final (128,64) profile stores 191/64,
about three coefficients per output symbol. Setup is fixed and may be
preprocessed. Its batch/window ratio is 1/2: the finite certificate proves
this instance directly, without invoking the stricter asymptotic corollary
in the Toeplitz note.

A direct Toeplitz product still performs L full-field multiplications
per output symbol. Reducing coefficient storage is not a runtime result.
A polynomial middle-product kernel is required to exploit the arithmetic
structure. Its multiplication, scratch, transformed setup, and route
traffic must be measured together. The small RS coefficients do not make
the uniform GF(2^128) inner coefficients cheap.

There is an exact kernel candidate for L=128,b=64. In the Toeplitz
polynomial representation, D has degree at most 190 and S has degree
at most 127. Only product coefficients 127 through 190 are needed.
Evaluate both polynomials at all 256 GF256 subfield elements, multiply
pointwise in F, and interpolate the degree-at-most-255 polynomial.
This computes D*S modulo Z^256-Z. Product terms of degrees 256 through
317 reduce to degrees 1 through 62, so they do not affect the requested
middle coefficients. No approximation or change of setup is involved.

Precomputing D's evaluations leaves one evaluation transform, 256 full-field
products, and one interpolation transform per batch: four full-field
products per output symbol, plus transform work, including any monomial-to-
additive-FFT-basis conversions and their inverses. The transformed coefficient
tape stores four F-elements per output rather than 191/64. The transforms
use GF256 coefficients, but an efficient representation for both those
operations and full-field multiplication was unresolved at this proof
checkpoint. In particular,
conversion between a byte-oriented tower basis and a CLMUL-friendly basis
must not be treated as free. This is an algebraic kernel plan, not a
measured implementation or a claim of four total operations per output.

The next implementation gate is the certified RS[64,32], L=128, b=64
profile at 35%, with b=32 as a control and the 45% profile retained as a
high-distance control.
Compare against the preserved 3.798 ms GF256 code under a matched harness.
Do not transfer that timing to these constructions. The older coupled-scan
proofs and all production and paper files remain unchanged.

The native gate has now been run. The existing sparse-byte F128 arithmetic
avoids a field-coordinate conversion, but the monomial transform interface
is expensive. The [spectral replacement](SPECTRAL_BATCH.md) changes the
inner distribution while preserving its exact fixed-output inverse law.
It retains this row's finite bound without another numerical cover.
See the [separate measurement record](NATIVE_BASELINE.md) for the remaining
performance gap; none of these timings replace the preserved GF256 result.
