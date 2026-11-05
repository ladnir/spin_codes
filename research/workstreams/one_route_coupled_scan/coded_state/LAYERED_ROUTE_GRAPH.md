# A routing bound for the layered inner's actual free coordinates

The affine proof in [LAYERED_INNER_PROOF.md](LAYERED_INNER_PROOF.md) leaves a
specific gap: its free coordinates could concentrate in outer blocks.
For dense states, those coordinates are not arbitrary. They occur in pairs.
This note uses that restriction without changing the encoder or its keys.
The GF32 sandwich remains the whole-code certified control.

The new result is restricted. It covers all low-weight J-valued outputs
whose state becomes dense within a bounded budget of pair locations.
Long or frequently recurring sparse states remain open. Here J=GF(2^32).
All probabilities concern setup, with the bad message chosen after setup.

## The route supplies random matching graphs

Keep N=2^21, M=8192 outer blocks, B=256 regions, and RS dimension k=128.
Each region independently permutes the block labels uniformly. The inner's
pair wiring partitions its 8192 physical positions into pairs, even though
the chosen index bit changes between batches.

Consequently, pulling the pairs back through the regional route gives an
independent uniform perfect matching of the M outer-block labels in each
region. This is a property of the existing route, not a new permutation pass.

Consider the last R regions. If both incoming state coordinates of a pair
are nonzero, the regular branch has two free inverse-input coordinates and
the exceptional branch has none. Thus a free set F_0 consists of endpoints
of selected matching edges. Its size is D=2m when m edges are selected.
The degree f_i of outer block i in the selected multigraph is its number
of free coordinates. Therefore

    d_pi(F_0)=sum_i(f_i-128)_+.

Let S_0 be the set of vertices with f_i>128, and put s=|S_0|. If E_R(S_0)
counts all matching edges internal to S_0 before any edges are selected,
then the selected degree sum on S_0 is at most each of

    2m,       m+E_R(S_0),       R*s.

The middle bound counts a selected internal edge twice and any other selected
edge at most once. A uniform bound E_R(S_0)<=e_R(s) therefore gives

    d_pi(F_0) <= d_R(m)
      := min(m, max(0, max_s[min(2m,m+e_R(s),R*s)-128s])).

The cap m also uses the route-independent RS inequality d_pi(F_0)<=D/2.
Crucially, this inequality holds simultaneously for every selected-edge mask.
There is no separate routing union over those masks.

## An elementary uniform internal-edge bound

For a fixed s-vertex set in one uniform perfect matching, let Z be its number
of internal edges. Counting specified collections of l disjoint edges gives

    E binom(Z,l)
      = [s!/(s-2l)!] / [2^l l! product_(i=0)^(l-1)(M-2i-1)].

Set n0=floor(s/2). For even s, put p=(s-1)/(M-1); for odd s, put p=s/(M-1).
Comparison of successive factors proves

    E binom(Z,l) <= binom(n0,l)*p^l.

For example, in the even case the residual ratio is
(s-2i-1)/(M-2i-1)<=(s-1)/(M-1). The odd case instead uses
(s-2i)/(M-2i-1)<=s/(M-1). These inequalities hold for every allowed i.
The empty set is deterministic and needs no tail bound.

Expanding exp(theta*Z) in powers of exp(theta)-1, for theta>=0, shows that
its moment is bounded by that of Bin(n0,p). Independent regions give the
same bound with n=R*n0. Two elementary Chernoff consequences are

    Pr[E_R(S_0)>=np+t] <= exp(-2t^2/n),
    Pr[E_R(S_0)>=np+t] <= exp(-t^2/[2(np+t/3)]).

For the first, the second derivative of the centered binomial log moment
is at most n/4. Integrate twice, then optimize theta. For the second, use
log(1+p(exp(theta)-1))<=p(exp(theta)-1) and
exp(theta)-1-theta<=theta^2/[2(1-theta/3)] for 0<=theta<3.
The latter follows term by term from j!>=2*3^(j-2), j>=2.

Let c_s be the bit length of binom(M,s), and set L=(7/10)(c_s+77).
Since ln(2)<7/10, either of the following integer-rounded slacks makes
the tail smaller than 2^(-c_s-77):

    ceil(sqrt(nL/2)),
    ceil(sqrt(2npL)) + ceil(2L/3).

The implementation takes the smaller, adds ceil(np), and caps by n.
All square-root ceilings and mean ceilings are computed with integers.
For s>M/2, use the identity

    E_R(S_0)=R*(s-M/2)+E_R(S_0^c).

Union over at most 2^c_s sets of each size, fewer than 2^13 sizes, and at
most 2^8 suffixes. Except with probability at most 2^-56, every required
e_R(s) bound holds simultaneously. The route event is independent of keys.
This proof uses no floating-point logarithms or statistical sampling.

## Combine the inverse and forward bounds on a dense suffix

First suppose the output is zero before the start of the last R regions.
Ignore the first batch's keys, since its incoming state is zero. Suppose
every subsequent incoming state coordinate is nonzero. Put

    u=R*M,       P=(u-512)/2.

If e of the P pair updates are exceptional, there are m=P-e regular pairs.
For a good route, the inverse count gives, for outputs of weight at most W,

    Pr[bad output in this class]
       < 2^(u+P+1) * q^(W-P-m+d_R(m)).                 (I)

Here 2^u bounds output-support choices, q^W bounds amplitudes, and 2^P
bounds exceptional-pair masks. The extra factor two bounds the actual key
density inflation, as in LAYERED_INNER_PROOF.md. No independence between
outputs and exceptional masks is assumed; the union includes both.

There is a complementary forward bound. A suffix U has shortened outer
dimension max(0,u-N/2), for every route. Charge e fresh exceptional roots
for each fixed outer message and union over their locations:

    Pr[bad output in this class]
       < 2^(P+1) * q^(max(0,u-N/2)-e).                (F)

Take the smaller complete bound for each e. A suffix of at most 128 regions
cannot support a nonzero outer word at all, so only R=129,...,256 matter.

## Complete coverage of exception counts, without a huge enumeration

The verifier requires each conditional class bound to be at most 2^-80.
Bound (F) identifies an initial interval of m values that already pass.
For the remaining interval [lo,P], maximize the W allowed by (I).

For fixed s, the function

    min(2m,m+e_R(s),R*s)-128s

has slopes 2,1,0, with breaks at e_R(s) and R*s-e_R(s). In particular,
minimizing m-d_R(m) reduces exactly to minimizing

    max(0, max(-m,-e_R(s),m-R*s)+128s)

jointly over s and m in [lo,P]. For each s, choose m by clamping e_R(s)
to that interval. This attains its minimum, including the possible flat
segment. Thus all integer e values are covered analytically, not sampled.
An independent heap-based sweep and small exhaustive graph checks test this
piecewise reduction.

Complete exact-integer results over J:

| Suffix regions R | Exception-count classes covered | Certified weight threshold W |
|---:|---:|---:|
|160|655105|868139|
|192|786177|856971|
|224|917249|886155|
|256|1048321|950021|
|All 129..256|100892800|854892|

The minimum occurs at R=182. Its threshold is about 40.7644% of N.
The total conditional union is less than 2^27*2^-80=2^-53. Adding the
route failure gives less than 2^-52 for this entire restricted family.
The percentages are descriptive; the integer thresholds define the results.

At R=192,e=314572, the new dimension bound is 208514, rather than the
coarse free-set bound 471604. That selected class now closes even at 45%,
with conditional log2 failure below -1018495. Covering all exception counts
is harder, as the complete thresholds above show.

## Allow warm-up, arbitrary starting batches, and occasional sparse states

For a fixed Y, call an incoming pair nonfull if at least one coordinate is
zero. Apart from the first batch, suppose there are at most h such pairs.
Condition on their keys when bounding this particular Y. Those inverse-input
coordinates then become fixed; the other pair keys remain independent.
This is legitimate because Q(Y), and hence the locations, are key-independent
in the inverse direction. The bound is uniform over the conditioned values.

If e full pairs are exceptional, the actual number of free edges is
m=P-h-e. Put m'=P-e. The bound d_R is nondecreasing, so d_R(m)<=d_R(m').
The inverse exponent is therefore at most the old exponent at m' with W
replaced by W+2h. The forward exceptional-root bound is unchanged. Since
the route bound holds for all free sets, no union over nonfull-pair masks
is needed in addition to the existing union over outputs.

Allow h<=6144 after the first batch. The complete threshold becomes

    W=854892-2*6144=842604, approximately 40.1785% of N,

with the same failure bound below 2^-52. The first nonzero output may occur
anywhere in its first region: choose the suffix starting at that region,
and include its initial zero batches in the nonfull-pair budget. Up to 24
full batches' worth of additional pair locations fit in this budget. This
allows a warm-up period and occasional holes, not arbitrary sparse histories.

The restricted theorem is over J. A general F-valued state can be dense
while each J-component is sparse. Full scalar-extension transfer requires
covering all J-valued bad outputs, not just this family.

## Checks, remaining gap, and next step

layered_route_gate.py runs 203 exact factorial-moment comparisons and 6144
selected-edge graph checks. The analytic minimization agrees with exhaustive
small sweeps; its four displayed large thresholds also reproduce the earlier
full integer sweeps. The complete suffix pass covers 100892800 (R,e) classes
through the proved piecewise reduction. No raw masks, graph dumps, or timing
claims are retained.

This is progress on the actual route/inner interaction, not a full 40% code
certificate. Histories with more than 6144 nonfull pairs remain uncovered,
unless another restricted result applies. The fixed-small-spatial-support
result in LAYERED_INNER_PROOF.md handles one such extreme, not every pattern.

Next: retain the one-sided pair's free coordinate instead of conditioning
away both coordinates. Combine that smaller loss with coordinate-excursion
rank, then extend to long sparse periods and multiple live episodes. The
graph bound and key-root accounting can stay; the construction need not change.
