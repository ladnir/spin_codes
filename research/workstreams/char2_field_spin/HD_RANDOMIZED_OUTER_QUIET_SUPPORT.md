# A global input-support mark for quiet words

A bounded experiment with the outer support floor gives no new certified
quiet occupancy range when the local moments are averaged over all packet
subsets before the fractional power. At the three selected witnesses, an
exact-coefficient derivative proves that every input mark v>1 makes this
particular upper bound worse. This is a limitation of the proposed local
surrogate, not of the randomized-outer construction or of support marks in
other bounds.

The qualification matters. A separate calculation that retains tile
profiles before the fractional power improves the r=4096 marker bound by
about 92.62 bits using v=257/256. That bound is still positive by more than
62000 bits. No whole-code distance certificate follows from either result.

## The support floor and its positive mark

Use the construction and probability space of
[the randomized-outer proof](HD_RANDOMIZED_OUTER_PROOF.md), and the
coefficient notation of
[the exact regional calculation](HD_RANDOMIZED_OUTER_QUIET_EXACT.md).
There are R=8192 outer message rows, 256 coordinate regions, 128 packets
per region, packet width 64, and scalar field size q=65536. The total
packet count is P=32768 and the output-weight threshold is D=838860.

Fix an active-row set of size r. Each nonzero outer message row encodes to
an MDS[256,128,129] word. The independent nonzero output scalars and the
subsequent route preserve its support. The routed input X to the inner
encoder therefore satisfies

    wt(X) >= 129r.

The artificial uniform measure in the outer comparison does not itself
satisfy this constraint. Retain the constraint in that comparison and,
for any v>=1, use the pointwise inequality

    1[wt(X)>=129r] <= v^(wt(X)-129r).

The probability still ranges over the construction's sampled scalars and
routes. The artificial uniform amplitudes are only a device for bounding
that probability. The support mark does not assert that these amplitudes
are actual outer codewords.

For a fixed packet subset S of size j, let K_S(v,z) be the uniform-input
moment of v^wt(X) z^wt(FX) restricted to CX=0. Let T_S(v,z) omit the
feedback condition. Here 0<z<=1. The zero input has probability q^-j and
contributes exactly that amount to either moment. Every other input has
weight at most j. Consequently,

    K_S(v,z) <= q^-j + v^j (K_S(1,z)-q^-j),
    T_S(v,z) <= q^-j + v^j (T_S(1,z)-q^-j).

Average over a uniform j-subset S and substitute the existing quiet upper
bound Khat_j(z) and exact free-final moment T_j(z). Define

    Ktilde_j(v,z) = q^-j + v^j (Khat_j(z)-q^-j),
    Ttilde_j(v,z) = q^-j + v^j (T_j(z)-q^-j).

The subtraction removes only the known exact zero-input contribution.
Khat_j remains an upper bound, rather than an exact quiet moment. For
j<=4, the feedback map is injective on the chosen positions and the
quiet moment is exactly q^-j, so the helper sets its nonzero part to zero
directly.

For 0<alpha<=1, define positive polynomials

    U_v(x) = sum_j binom(64,j) Ktilde_j(v,z)^alpha x^j,
    V_v(x) = sum_j binom(64,j) Ttilde_j(v,z)^alpha x^j,
    A_r(v) = [x^r] U_v(x)^128 / binom(R,r),
    H_r(v) = [x^r] U_v(x)^127 V_v(x) / binom(R,r).

The free final packet remains in H_r. The route-fixed clip precedes the
fractional power; concavity is applied to the within-packet subset
average afterward. With beta=q^256/(q-1)^128, the resulting class bound is

    B_r(v) = binom(R,r)
             [beta^r z^(-D) v^(-129r)/(q-1)]^alpha
             A_r(v)^255 H_r(v).

The outer-row binomial remains outside alpha. The support-floor penalty
and projective division remain inside alpha. This is the same regional
coefficient calculation as before, with marked local moments.

## A bounded numerical comparison

The helper screened r in {1024,2048,4096}, using

    z = 2^-h,       h in {0,2,4,8,12,14,16},
    log2(v) in {0,1/2,1,2,4,8,12,16},
    alpha in {1,3,8,12,16,32}/32.

For each tuple, a floating saddle-point calculation proposed the regional
coefficient marker. The best marker-bound proposal for each r used v=1.
The selected proposals were replayed with 384-bit Arb arithmetic and
exact regional coefficient extraction. Displayed decimals below are
rounded for comparison, not outward-rounded certificate endpoints.

| Active rows r | z | alpha | Selected v | Exact-coefficient log2 upper bound |
| ---: | ---: | ---: | ---: | ---: |
| 1024 | 1 | 1/32 | 1 | 5397.046715 |
| 2048 | 1/4 | 1/32 | 1 | 52528.659853 |
| 4096 | 1/4096 | 1/32 | 1 | 75014.791267 |

The nonzero-input term can be much smaller than q^-j. For example, at
z=1/65536 and j=5, converting Khat_j and q^-j separately to ordinary
floating point loses their difference. The helper subtracts them in Arb
first and only then converts the logarithm of the positive difference
for witness selection. Directed evaluations never use the floating
difference. An unresolved positive difference causes an error rather
than silently being replaced by zero.

This finite screen alone does not exclude an improvement at another
input mark or another z and alpha. The next argument gives a stronger,
but explicitly fixed-witness, stopping rule.

## Convexity gives a fixed-witness stopping rule

Put t=log2(v). Each local surrogate has the form

    a_j + b_j 2^(jt),       a_j=q^-j, b_j>=0.

It is log-convex in t. Raising it to alpha, multiplying positive factors,
and adding them preserve log-convexity. A polynomial coefficient here is
a finite sum of such products. Thus log2 B_r(2^t) is convex: its support
penalty contributes only the affine term -129r alpha t.

The derivative at t=0 can be evaluated with the same exact polynomial
coefficients. For the quiet factors define

    s_j = j (Khat_j-q^-j)/Khat_j,
    dotU(x) = sum_j binom(64,j) Khat_j^alpha s_j x^j.

Define dotV analogously using T_j. All polynomials in the following
expressions are evaluated at v=1. Set

    s_A = 128 [x^r] U^127 dotU / [x^r] U^128,
    s_H = [x^r](127 U^126 dotU V + U^127 dotV)
          / [x^r] U^127 V.

The factors alpha and ln(2) cancel as appropriate when differentiating
the base-two logarithm, giving

    d/dt log2 B_r(2^t) at t=0
        = alpha (255 s_A + s_H - 129r).

Positive Bernoulli scaling and normalization of U,V,dotU,dotV leave these
coefficient ratios unchanged. Directed 384-bit evaluation gives the
following strictly positive derivatives, shown approximately:

| Active rows r | z | alpha | Derivative at t=0 |
| ---: | ---: | ---: | ---: |
| 1024 | 1 | 1/32 | 3746.631236 |
| 2048 | 1/4 | 1/32 | 8125.656313 |
| 4096 | 1/4096 | 1/32 | 8596.095312 |

Convexity and each positive derivative prove that no v>=1 improves the
corresponding fixed z,alpha exact-coefficient witness. They do not exclude
jointly changing z or alpha, retaining additional local information, or
using a more faithful input-weight moment. This stopping rule applies to
the upper-polynomial calculation, not to the actual failure probability.

## Retaining tile profiles changes the comparison

The independent
[tile-conditioned calculation](HD_RANDOMIZED_OUTER_QUIET_TILE_CONDITIONING.md)
averages the fractional power only after conditioning on the eight tile
occupancies. It applies the same exact zero-input subtraction separately
to each profile. Its directed optimized-marker comparisons are:

| Active rows r | z | alpha | v=1 log2 bound | v=257/256 log2 bound |
| ---: | ---: | ---: | ---: | ---: |
| 1536 | 1/4 | 1/32 | 8584.890630 | 8616.685303 |
| 4096 | 1/4096 | 1/32 | 62367.062196 | 62274.439422 |

These are marker bounds, not exact-coefficient bounds. The second row
improves by about 92.62 bits. The derivative certificates in the preceding
section do not apply to these different polynomials. No additional tile
derivative or broader input-mark sweep was run for this note.

## Tests, reproduction, and next step

The implementation is
[theory/hd_randomized_outer_quiet_support.py](theory/hd_randomized_outer_quiet_support.py).
Its five tests all pass. They include exhaustive GF(4) input enumeration
for the support-mark inequality, exact zero-input packet behavior, an Arb
subtraction that ordinary floating point loses, a small exact-coefficient
log-convexity check, and agreement of the directed derivative with a small
forward difference. From the theory directory:

    python -B -m unittest hd_randomized_outer_quiet_support_tests -v
    python -B hd_randomized_outer_quiet_support.py

The second command repeats the bounded screen and replays only its three
selected witnesses and derivatives. These are proof computations; no
native code, performance measurement, or remote execution is part of this
experiment.

The next useful proof iteration is to retain joint input/output weights
before applying the input mark and fractional power. For a selected
a-column Cauchy tile, the systematic map [I_a;F_U] generates an
MDS[a+8,a,9] code. Its bivariate enumerator can retain the actual input
weight instead of charging every nonzero input the maximal weight a.
Composing that information with the cross-tile quiet-feedback condition
requires a separate argument. No such implementation or certificate is
claimed here; the present experiment stops at the bounded comparison and
the fixed-witness convexity certificates above.
