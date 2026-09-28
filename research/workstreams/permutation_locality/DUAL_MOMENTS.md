# Average shortening counts from the dual code

The complementary-dimension bound in [EVEN_COLUMNS.md](EVEN_COLUMNS.md)
uses the worst dual shortening dimension. This refinement instead bounds
averages over all coordinate sets. It improves outer counts without
changing the encoder or its setup distribution.

## Positive identity

Let C be a binary [n,k] linear code and D=C^perp. For a coordinate set S,
write C_S for the subcode supported on S. Fix t>=n-k and set a=t-(n-k).
Complementary shortening gives

    dim C_S = a + dim D_(S^c),       |S|=t.

Write G(d,h) for the number of h-dimensional subspaces of F_2^d. Set
G(d,h)=0 when h>d. For nonnegative a,b, counting subspaces by projection
onto the b-dimensional summand gives

    G(a+b,h) = sum_{r=0}^h 2^(r(a-h+r)) G(a,h-r) G(b,r).

Terms with G(a,h-r)=0 are zero. Otherwise the exponent is nonnegative.
To see the count, choose the intersection with the first summand and the
projection onto the second. A linear map from that projection to the
quotient of the first summand then specifies the subspace.

Define M_C,h(t)=sum_{|S|=t} G(dim C_S,h), and define M_D,r similarly.
Taking sums in the identity yields

    M_C,h(t) = sum_{r=0}^h 2^(r(a-h+r)) G(a,h-r) M_D,r(n-t).

Every coefficient is nonnegative. Upper bounds on the dual moments
therefore give an upper bound on the primal moment. The r=0 moment is
exactly binom(n,n-t).

Let B_D,r(v) count r-dimensional dual subcodes with support size at most v.
Then

    M_D,r(s) = sum_{v=0}^s (B_D,r(v)-B_D,r(v-1)) binom(n-v,s-v).

An upper CDF may replace B in this expression: the binomial weights are
nonincreasing, so summation by parts has nonnegative coefficients.
Every primal h-subcode supported on at most u coordinates belongs to at
least binom(n-u,t-u) sets of size t. Dividing the moment upper by this
number bounds its support CDF. Multiplication by the number of spanning
ordered four-tuples converts this into the outer count used by SPIN.

All these operations use integers or rationals. The method can alternate
between C and D: each finite iteration uses only previously valid bounds.
The implementation performs three iterations, not an assumed fixed point.

A later finite-iteration check used the repaired sandwich dual shells and
compared three iterations with twenty. It tightened 210 integer CDF entries,
but the log2 totals at supports 96, 112, 128, 144, 160, 176, and 192 were
unchanged to ten decimal places. Thus more iterations of these same
inequalities did not materially improve the tested intermediate supports.
The default remains three; no fixed-point assumption is introduced.

## Dual inputs for the implemented BCH code

The program replays dual distance at least 30 and checks Q<=C<=P against
the production generator. The same check verifies that C is even and
contains the all-ones vector. Thus D is also even and contains all ones.
The primal distance lower bound 38 makes D an orthogonal array of strength
at least 37.

Let K_j be the binary Krawtchouk polynomial of degree j. For each weight w,
the degree-18 kernel gives the exact shell upper

    A_D(w) <= floor(2^128 / sum_{j=0}^{18} K_j(w)^2/binom(256,j)).

Indeed, square the kernel centered at w and normalize its value at w to
one. Its expectation under uniform D equals its binomial expectation,
because its degree is at most 36. Orthogonality makes that expectation
the reciprocal of the displayed denominator. The square is nonnegative
at every other weight, so the shell inequality follows.

These bounds are intersected with the distance, parity, all-ones, total-size,
and binomial bounds. The existing basis and containment-moment arguments
then bound dual subcodes. The 375 exact tests include Gaussian identities,
averaged complement identities, exhaustive small-code tuple counts, and
kernel shell bounds. They supplement the algebra above, not replace it.

## Effect on outer counts

The table gives log2 upper counts of nonzero ordered four-tuples with union
support at most u. Both columns already use positive-polynomial shortening,
dual-complement dimensions, and primal containment moments.

| u | Before averaged dual identity | After three iterations |
|---|---:|---:|
| 96 | 254.2645 | 254.2645 |
| 112 | 303.6766 | 303.6766 |
| 128 | 348.9427 | 341.7826 |
| 144 | 374.9980 | 358.5181 |
| 160 | 396.1708 | 379.4003 |
| 176 | 419.8266 | 406.8611 |
| 192 | 450.4205 | 444.8091 |

The underlying bounds are exact; the table logarithms are rounded.
These are outer counts, not setup-failure probabilities.

    python -B research/workstreams/permutation_locality/dual_moments.py

## Occupancy integration

The optional `--dual-shortening` and `--dual-moments` flags enable the
refinements in `occupancy_cdf_cover.py`. Existing defaults are unchanged.
At q=64, tilt .032, all-one penalty .75, and two updates, the screen with
dual dimensions but without averaged dual moments gives:

| Equal support in each active group | Without averaged dual moments | With averaged dual moments |
|---|---:|---:|
| 96 | +1840.4742 | not rerun |
| 128 | +3298.6653 | +3046.1426 |
| 144 | +3373.8855 | +2593.1769 |
| 160 | +3224.6792 | +2396.4237 |
| 176 | +3274.8276 | +2335.6182 |

The screen includes exact feedback distributions through six active
windows and the existing state refinements. Values are log2 contribution
uppers from binary64 diagnostics. These selected points remain
open; positive values provide no nontrivial failure bound.

Reproduce the screen, adding averaged dual moments, with:

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-moments --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties .75 --tilts .032 --probe-supports 128 144 160 176

Omit `--dual-moments` to reproduce the earlier column. No complete occupancy
cover, full-code certificate, production default, or paper claim is changed.

The support-128 improvement is smaller than 64 times its unweighted count
gain because the inner bound also charges for all-one columns. That charge
uses separate joint-count upper bounds. Thus an improved total count does
not translate directly into the same weighted-count gain.

Next investigate tighter low-weight dual spectra or joint all-one-column
counts. The current dual shell caps use only the distance-derived moments;
they do not exploit the full existing BCH sandwich LP. Any reused LP must
exclude empirical orbit-search lower bounds and produce a checked rational
witness. The present diagnostic still needs substantial improvement, so
another broad witness grid is not the immediate next step.
