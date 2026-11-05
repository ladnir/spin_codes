# Retaining fixed-input locations and active outer rows

Follow-up: [count-region certificates](LAYERED_COUNT_REGIONS.md) extend this
single-profile result to gapless regions across neighboring suffix lengths.
This note preserves the original proof and its integer-power replay.

This continues the [all-pair-type count](LAYERED_ALL_PAIR_TYPES.md) with the
outer, regional permutations, and layered inner unchanged. The test profile is

    R=192, z=49136, r=24568, f=589632, e=368520, B=49136,
    W=838860, d_bar=153246.

It is the retained large-empty-pair obstruction at the exploratory 40%
threshold. A result for this profile is not full count-space coverage or a
near-half distance theorem. The GF32 sandwich remains the whole-code control.

**Result:** the profile now closes over every active-block count, all support
and exceptional placements, and all amplitudes. The exact good-route
contribution is at most 2^-5153. This replaces the previous failed bound
+491591 for the same aggregate profile. The construction is unchanged.

## What is lost by counting only zeros

Fix a routed history, including supports, exceptional branches, and internal
zero-output positions. Let Z be its linear state space. The fixed-input map
L maps Z to the fixed input coordinates; its kernel has dimension kappa.

Choose a basis matrix for Z and form the matrix U of L in that basis. A basis
H of the left kernel of U gives every linear relation on the fixed inputs.
If G_fix is the outer generator restricted to those coordinates, then

    K_compatible = K_R - rank(H G_fix),
    dim{(state,message): L(state)=G_fix(message)} = kappa + K_compatible.

This identity retains locations, temporal components, and all outer relations
together. Zero rows of U are forced-zero inputs. Other rows of H can express
relations between nonzero fixed inputs. Thus the zero rows alone need not
capture the full rank saving.

The toy verifier checks the identity on 128 routed GF4 fixtures. In 63 cases,
nonzero relations add rank beyond the forced zeros; the largest extra rank is
four. These fixtures validate the algebra, not the K20 class bound. A simple
negative control also shows why zero locations matter: four zeros in one
RS[4,2] block remove two dimensions, while one zero in each of four blocks
removes four.

## A tractable location-aware bound

Write q=2^32, M=8192, and k=R-128. Each RS block shortened to the R-region
suffix has dimension k. Partition outer words by their exact set S of active
blocks, and put s=|S|. Let J and U_+ be the fixed input positions known to be
zero and nonzero, respectively.

Every position in U_+ must route into S. A position in J that routes into S
imposes a zero in an active RS block. If j_i zeros hit active block i, it has
at most q^(k-j_i) nonzero words when j_i<k, and none when j_i>=k. Consequently,

    #compatible outer words with active set S
        <= q^(ks-|J inside S|) * 1[U_+ inside S].

The expression remains an upper bound when a row has too many zeros: its
actual contribution is zero. Positions outside S impose no dimension cost.

The first inverse bound uses d_S, the shortening dimension on the free
coordinates within S. On the existing good-route event,

    d_S <= min(d_bar,ks).

For one history, interpolate the state-count and compatibility bounds with
theta=m/32, where m is an integer from 0 through 32. Before summing histories,
the contribution for S is at most

    eta * q^[theta*(a+d_S)+(1-theta)*(kappa+ks-|J inside S|)-D-e]
        * 1[U_+ inside S],

where eta<2. The two bounds apply to the same history and the same active
set; no independence between their constraints is asserted.

## Averaging the actual regional permutation

For fixed S, the positions routed into S form a uniform s-subset in each
region. Different regions are independent. To bound this distribution, assign
independent Bernoulli colours of probability p to all positions, then condition
each region to have exactly s coloured positions. For any 0<p<1, this gives
exactly the required uniform-subset distribution.

The probability of the conditioning event in one region is

    pi_s(p) = binom(M,s) p^s (1-p)^(M-s).

For a nonnegative weight, dropping the conditioning event in the numerator
costs at most pi_s(p)^(-R). Thus every forced-nonzero position contributes p,
and every forced-zero position contributes

    beta = 1-p+p*2^(-(32-m)).

Include binom(M,s) for the choice of S. The value p is a proof parameter,
not a change to the permutation distribution. It need not equal s/M.

The route envelope is used only on the good-route event. After replacing
d_S by its simultaneous envelope, the positive right-hand side can be
averaged over all routes. The resulting bound controls the expected bad-word
count intersected with that event; it is not a conditional expectation without
the appropriate normalization. The existing bad-route probability is added
once in a future whole-code certificate.

## Keeping input-zero choices in the component transfer

The existing local fixed inputs have form X=Q_next+gamma Q. For gamma=0,
the next support determines whether X is zero. For gamma=1, unequal adjacent
supports force X nonzero. A zero output also forces X nonzero when both states
are active, because Q_next=delta Q and delta!=1.

The remaining case has gamma=1, both adjacent states active, and a nonzero
output. Split it into X=0 and X!=0. The former adds Q_next=Q to the actual
state constraints. These new equations occupy distinct temporal edges and
are independent of the prescribed zero-output edges, so each reduces a by
one. They do not reduce kappa again: the auxiliary kernel already imposes
the same gamma=1 equation.

After interpolation and route averaging, the combined local factor is

    p + beta*2^(-m).

The two summands retain the nonzero and zero alternatives. This is not a
claim that X is uniform. Initial fixed inputs receive their zero/nonzero
factors explicitly; terminal input classifications are omitted conservatively.

These factors multiply the component-completion, exceptional-update, birth,
and output-weight factors in the same positive transfer. Both the one-update
operator and the existing two-successive-update operator are supported.
The latter retains pair correlation on four coordinate paths. No encoder
operation or setup distribution changes.

## Interval verification

For fixed rational markers, p, m, and local potentials, the operator norms
do not depend on s. Above s=5900, d_S's upper bound is the constant d_bar.
All remaining dependence on s has the form

    constant + slope*s + (1-R)*log2 binom(M,s).

The binomial logarithm is discretely concave; R>1 makes this expression
convex. Its maximum on an integer interval is attained at an endpoint.
The verifier uses exact rational row maxima and integer-power logarithm
bounds at both endpoints. Floating proposals only choose intervals and
weights; verified endpoints and a gap check establish the coverage.

For 1<=s<=5900, the previous two-update bound suffices without location
weights. Replace K_R by ks in its consistency branch and union over at most
2^M active sets. Its log2 upper bound is

    556423 - 4*64*(8192-5900) + 8192 = -22137.

This includes all choices of S in that range, conditional on any good route.
The location-aware intervals handle the remaining range. At integer-power
precision 2048, their verified per-s log2 upper bounds are:

| Active blocks s | Exact upper bound |
| --- | ---: |
| 5901–6903 | −220882 |
| 6904–7636 | −10905 |
| 7637–7990 | −5166 |
| 7991–8061 | −5229 |
| 8062–8125 | −6687 |
| 8126–8192 | −8549 |

The union over all s, including the already-unioned small-s range, costs
at most 13 additional bits. The resulting good-route upper bound is −5153.
Precision 512 was insufficient to certify two of the narrow intervals;
the improvement to 2048 only tightens the exact logarithm rounding error.
Some bounded proposals stopped at their time cap. Every retained proposal
is validated by recomputing its exact local row maxima, independently of
optimizer termination or claimed feasibility.

The executable verifier reports uncovered intervals explicitly and never
substitutes passing sample points for coverage. It conditions no new random
encoder maps and introduces no additional bad-route event.

The checks also cover 972 exact regional-colour inequalities, 528 actual
GF4 fixed-input classifications, 1,188 exact/float local rows, and 11 added
input-zero equations in the routed matrix fixtures. The latter lower the
actual state dimension without double-counting a reduction of kappa.

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_outer_locations.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/verify_layered_location_cover.py --verify

Only compact rational proof inputs are retained. No benchmark, raw timing
dump, production change, or paper claim is part of this investigation.

Next apply this bound to nearby empty-pair count profiles and then bound
regions of the remaining count domain. The proof still omits some temporal
relations among nonzero fixed inputs; the exact quotient test records those
as a possible further improvement. No full layered-code certificate follows
from closing this single count profile.
