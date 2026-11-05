# Outer compatibility for one-diagonal feedback

Status: an additional counting interface, not a distance certificate.
It applies to the existing unshuffled packet route. No construction
change or additional randomness is required.

The inverse input is uniform only on selected coordinates. Its remaining
coordinates must already agree with some outer codeword. Earlier bounds
discarded this deterministic compatibility condition. Retaining it gives
a complementary bound indexed by the union of output and inverse-randomness
supports. Masks with the same union can be grouped before taking a
fractional moment.

## Exact compatibility and joint dimension

Fix a route and let P be a full-row-rank parity-check matrix for its routed
outer code C=ker P. Fix the deterministic strictly causal map L taking an
output Y to the state observations V=LY in SINGLE_DIAGONAL.md. The inverse is

    X=Y+diag(r)LY,

where the coordinates of r are independent uniform elements of F. Put
q=|F|. This setup is sampled once and fixed for all messages.

Fix exact supports W=supp Y and A=supp(LY), and put S=A union W and w=|W|.
For a coordinate set T, let P_T contain the columns of P indexed by T;
write r_T=rank(P_T). For each fixed eligible Y, diag(r)LY is uniform on
F^A. Therefore

    Pr_r[X in C]
      = q^(-r_A) if P Y belongs to image(P_A),
        0 otherwise.                                      (1)

Let Q_A be any matrix whose kernel is image(P_A). The compatibility
equations on the free output coordinates are Q_A P_W Y_W=0. Their rank is

    r_out=rank(Q_A P_W)=r_S-r_A.                            (2)

For the direct sum of MDS [B,B/2] outer rows, this rank is explicit. If
a_i=|A intersect row i| and u_i=|(W minus A) intersect row i|, then

    r_out=sum_i min(u_i,(B/2-a_i)_+).                       (3)

Indeed each parity-check restriction has rank min(B/2,number of columns).
The formula applies to actual routed scalar coordinates, including partial
packets. It needs no randomized constituent generator.

The state-zero and compatibility equations act on the same output variables.
Their joint solution space has dimension

    d(W,A)=w-rank([ L_(A^c),W ; Q_A P_W ]).                 (4)

The brackets denote vertical concatenation. This counts a linearized set;
required nonzero coordinates can only reduce it. If d(W,A)=0, the class
contains no nonzero output. Otherwise its expected number of eligible
projective output lines is at most

    ((q^d(W,A)-1)/(q-1)) q^(-r_A).                         (5)

All conditions are invariant under common nonzero scaling, which justifies
counting lines. Equation (5), capped at one, bounds the class's failure
probability over setup. It does not concern a message fixed before setup.

If E is a previously established dimension bound using only state-zero
equations, (4) gives d<=min(E,w-r_out). Consequently the two simpler
expected-word bounds are

    q^(E-r_A),             q^(w-r_S).                     (6)

Take their minimum; do not add r_out to the earlier state-zero rank.
Those equation families can overlap. A larger improvement requires a
bound on their joint rank in (4).

## A valid additive rank from overlap columns

Let O=A intersect W. The compatibility matrix Q_A P_W is zero on every
column indexed by O, since Q_A annihilates P_A. Define

    R_overlap=rank(L_(A^c),O).

First select R_overlap independent state-equation pivots on those columns.
Compatibility rows have zero entries there and retain their rank on the
remaining columns. A block-triangular minor therefore gives

    rank([ L_(A^c),W ; Q_A P_W ]) >= R_overlap+r_out.        (7)

Equivalently, restrict the joint row space to its O columns. The image
contains R_overlap dimensions, and the kernel contains the r_out-dimensional
compatibility row space. This argument does not require generic minors
across several observation times.

The selected-observation bound in SINGLE_DIAGONAL_RANK.md applies to
L_(A^c),O. Use the same observation deficits b-h_t, but replace source
counts w_t by o_t=|O intersect batch t|. Any valid selected-time bound
R_o<=R_overlap then strengthens the second bound in (6) to

    q^(w-R_o-r_S).                                        (8)

For example, maximize sums of

    min(b-h_(t_i), sum_(t=t_(i-1))^(t_i-1) o_t)

over ordered observation times, with t_0=0. A prescribed closed terminal
can also contribute a full-syndrome observation on the overlap columns.
No such observation is available at a free terminal. Bounds from several
observations must use the proved partition rule, not a generic-minor claim.

## Grouping by the union support

Fix per-batch statistics h_t=|A_t|, w_t=|W_t|, and o_t=|A_t intersect W_t|.
Let s_t=h_t+w_t-o_t, s=sum_t s_t, and w=sum_t w_t. All batches contain
complete routed packets. For each fixed union support S with these s_t,
the number of mask pairs in the complete class is

    L_S=product_t binom(s_t,h_t) binom(h_t,o_t).            (9)

The bound (8), using an R_o determined only by h_t,o_t, is common to
every such pair. Sum the expected
projective counts over those pairs before taking a fractional power.
For 0<nu<=1 and tau=q^nu, the fixed-S failure probability is at most

    (q-1)^(-nu) L_S^nu tau^(w-R_o-r_S).                   (10)

This is an upper bound even if the mask-pair counts overlap. Exact support
classes in fact partition the eligible outputs.

There are M_S=product_t binom(b,s_t) union masks in the complete family.
They are closed under independent packet-lane permutations. Thus the
mask-orbit identity and moment M_p from SHUFFLED_PACKET_ROUTE.md apply
to their sum under the original unshuffled route. Since r_S=s-D_out(S),

    Pr[an eligible output in the class]
      <= min(1, (q-1)^(-nu) M_S L_S^nu
                    tau^(w-R_o-s) M_p(s,tau)).            (11)

Only the union-mask multiplicity remains outside the fractional power.
This is different from grouping W at fixed A: neither grouping dominates
universally. Compare complete bounds for the same class, with their actual
multiplicities. Neither bound permits discarding an uncounted support class.

Taking R_o=0 retains the earlier bound that drops all state-zero ranks.
The overlap refinement charges only ranks proved independent of outer
compatibility. It complements, rather than replaces, the full delayed-rank
bound. Nearly full S has few union masks, which is the first class where
this grouping may help. It proves no
finite margin until the full class sum is evaluated soundly.

## Retain fully selected and empty regions

The arbitrary-set moment averages away some information available from a
support profile. In particular, a dense suffix can select complete regions.
Every outer row then contributes a coordinate in each of those regions,
regardless of the sampled packet route.

Fix a coordinate set T. Let f be the number of regions wholly contained
in T, and let t be the number of regions intersecting T. Put M=N/B,
the number of outer rows, and k=B/2. For every route, each row's selected
occupancy i lies in [f,t], and the sum of occupancies is |T|. Therefore

    D_out(T) <= C(T),

where

    C(T)=0,                                      if t<=k;
    C(T)=|T|-Mk,                                 if f>=k;
    C(T)=floor((|T|-Mf)(t-k)/(t-f)),               if f<k<t.       (12)

The second case is exact. In the third, convexity bounds (i-k)_+ by the
chord (i-f)(t-k)/(t-f). Sum over rows and round down, since the shortened
dimension is an integer. At t=f one of the first two cases applies.
This argument uses no independence assumption and no coefficient
randomness. It is valid for the physical unshuffled packet route.

For H- or S-grouped profile bounds, f and t are determined by the
per-batch support counts. The complete mask family has the same full
and empty regions, even after packet-lane symmetrization. Thus its route
moment can be capped by tau^C(T), as well as by tau^(|T|/2) and the
previous activity-union cap. Use the cap separately for T=A and T=S;
the coarser whole-batch activity set U can have neither full-region
distinction. No rank is added to the inner in this refinement.

## Zero-dimensional classes and checks

Discard a class for zero dimension only after proving that its actual
joint linearized space in (4) is zero, or a valid upper bound is zero.
For example, if A is empty and D_out(W)=0, compatibility alone forces Y=0.
In contrast, D_out(A)=0 does not generally exclude an affine intersection.
Nor does a small probability in (1) imply a zero-dimensional output class.

`theory/test_affine_compatibility.py` provides bounded GF16 checks of (1),
the MDS rank difference, the joint matrix, exact linearized output counts,
projective counting, the overlap-column rank addition, and the
zero-dimension guard. A negative control rejects adding the separate ranks
blindly. Root's bounded serial replay passes 10240 joint-rank fixtures,
including the overlap-column and selected-observation inequalities. It
finds 2294 cases where adding the unrestricted ranks would overcount.
It also checks 13225 randomizer assignments, 40 exact joint spaces,
6280 output assignments, and 32 zero-dimensional classes. These local
checks are not a K20 certificate.

The following diagnostic compares (11) with the A-grouped delayed-rank
bound on the same structured support classes. It retains the free
terminal and actual packet moment, without assuming independence between
the outer and state equations.

## Bounded profile diagnostic

`theory/screen_single_profiles.py` initially evaluates 18 named absolute support
profiles at K20, B64, packet width four, batch size 128, and q=2^32.
It retains activation, observed resets, the free terminal, per-batch mask
multiplicities, and the actual partial-packet moment. It compares complete
H-grouped and S-grouped bounds, rather than adding incompatible ranks.
Late-only activity touches one outer region, so D_out(U)=0 excludes it.

The run first passes 64 exact comparisons of the selected-time rank
algorithm with its direct quadratic recurrence. All profiles then finish
within the 30-second/256-MiB caps. Scores use ordinary floating point and
ten fixed fractional markers. They are not directed certificates, and
the script does not sum other profiles, translations, or placements.

Examples below give uncapped log2 upper-bound expressions. Each live
period follows an initial (h,w,o)=(0,128,0) activation batch. The entire
sequence has 16384 batches; the terminal is free.

| Repeated profile (h,w,o) | Earlier ell-only H bound | Delayed H bound | Overlap S bound |
|---|---:|---:|---:|
| (86,44,44) | +1716192 | -18286748 | -18287452 |
| (96,44,44) | +1053565 | -15721603 | -15722627 |
| (128,128,128),(128,0,0),(1,0,0) | -22238605 | -22238605 | -22240653 |

The ell-only relaxation loses the fresh-column rank in the first two
profiles. Delayed ranks recover that loss. The overlap bound offers a
complementary grouping, but it is not responsible for every improvement.
For example, the resetting profile `(0,64,0),(65,65,65),(0,0,0)` has an
uninformative positive S-overlap expression; the H-grouped expression is
about -4146486. Neither grouping may be discarded globally.

All eligible initial profiles have a negative expression under at least
one valid bound. One deliberately heavier profile lies above the 35%
threshold and is labelled ineligible. These results remove the named
profiles as diagnostic obstructions; they do not establish a minimum
distance or a margin for the full code.

### A phased profile exposes a regional loss

A targeted nineteenth profile is more informative. Start with the same
activation batch, followed by p=10734 copies of (65,1,1), then m=5649
copies of (128,128,128). Its terminal state is free. This gives

    w=733934 <= floor(0.35 N)=734003,
    H=1420782,     R_partition=p+62=10796,
    E<=723138.

The simpler cap D(A)<=H/2 leaves E-H/2=12747>0. Even the full
dimension-only relaxation permits E=723076, so retaining its known
column core does not resolve this global inequality. This is an
admissible relaxed support profile, not a realized low-weight codeword.

The previous H/S route-moment expressions also fail: their best log2
value on the fixed grid is about +1174431. However, the dense suffix
contains 22 whole regions. Every row contributes to A in those regions.
Equation (12) gives, for every route and every mask in this class,

    D(A)<=floor((1420782-22*32768)*32/42)=533246.

Hence E-H+D(A)<=-164398. This removes the failed bound without changing
the construction or asserting stronger multitime inner ranks. The new
floating H-regional expression is about -3863508. The bounded replay
also passes 48105 exact pointwise secant checks and 144315 integer
aggregation checks.

This one absolute profile admits a simple exact class bound. The sparse
phase has at most binom(128,65)*65 mask pairs per batch; all other masks
are fixed. Since binom(128,65)<=2^128 and 65<=2^7, the number of pairs is
at most 2^(135p). Summing the conditional inverse-word bounds therefore
gives

    Pr[an eligible nonzero output in this absolute profile]
      <= 2^(135*10734) (2^32)^(-164398)
       = 2^(-3811646).                                  (13)

This bound holds uniformly over the packet route and averages only the
independent inner coefficients. It omits the optional projective saving.
It does not union over other profiles, phase lengths, or placements and
is not a full-code margin. Actual realization of the sparse-phase
syndrome alignment is unnecessary for this upper bound.

The next refinement should retain regional occupancy information in a
complete count. For near-full regions, one valid candidate is to fill
their missing positions virtually, apply (12) to that larger support,
and use monotonicity of shortening. Another is to retain regional
Bernoulli marginals before the equalization in the route-moment proof.
Neither refinement is implemented here. Do not infer complete coverage
from the now-handled finite list of support patterns.
