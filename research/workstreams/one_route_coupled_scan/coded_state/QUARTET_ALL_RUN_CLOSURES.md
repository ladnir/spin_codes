# Global availability closures and joint profile boxes

This checkpoint keeps the construction, setup distribution, outer, route,
field and free terminal fixed. It continues the global low-rank proof;
local seven-type refinements remain parked. The full 45%-distance,
40-bit certificate is still open. Notation is inherited from
[the forced-run checkpoint](QUARTET_GLOBAL_RUN_COST.md).

## What the joint profile calculation resolves

For a canonical footprint with rank R, available-position count L and f
forced births, use the grouped conditional field bound

    Phi(R,L,f) = min(1, (gamma/q) q^(-R)
      inf_(0<x<=1) x^(-W) c^f lambda^(L-f) e^max(0,f-512)).

Here c=(q+1)x, lambda=1+c, e=min(1,lambda/(q+1)), and gamma<2.
For a box of integer profiles, Phi decreases with R and f and increases
with L, so its upper endpoint is (R_min,L_max,f_min). Each forced birth
starts a distinct run; all but at most 512 close. Hence 2f-512>W is
infeasible. The algebraic x optimizer is implemented in
`quartet_rank_boxes.py` and checked against direct evaluations.

The geometric quantity G_B is the expected NUMBER of footprints in the
box, not the probability that a geometry has a footprint. Its contribution
is at most G_B sup_B Phi. Capping G_B at one before this multiplication
would be invalid. H is counted once over the complete code, not once per
region. The positive formal annotations of the existing transfer do not
have to represent physical availability profiles; apply the bound first
to the canonical coefficients, then enlarge positively.

Selected h=2048 boxes, each of half-width 1024 in (R,L,f,F), give the
following FLOAT diagnostics with one shared geometric moment:

| Box center (R,L,f,F) | log2 G_B | log2 Phi | log2 product |
|---|---:|---:|---:|
| (903775,891035,71426,1132763) | 155388 | -797226 | -641838 |
| (900000,1000000,40000,1200000) | -36169 | 0 | -36169 |
| (900000,1100000,20000,1200000) | -201169 | 0 | -201169 |
| (900000,900000,20000,1200000) | 111331 | 0 | 111331 |

These are selected class bounds, not a cover or directed certificate.
Correction from the subsequent canonical-feasibility audit: the last
three boxes are EMPTY at h=2048. Their F_min-R_max=297952 exceeds
128h=262144. Their loose upper bounds therefore do not identify a
remaining class or a low-birth obstruction. The first box survives the
aggregate screen and its displayed conditional bound is negative.
The replacement test uses F=1150000, with the other coordinates of the
last row unchanged; it is not excluded by the aggregate tests, but is
not asserted to be realizable. Finer boxes alone need not fix a genuine
remaining gap: for a fixed positive
generating function, unconstrained multivariate Chernoff bounds can merely
recover the same common-moment variational bound. Nonconvex class
exclusions, switching field branches, a tighter geometry count, or the
quota-tilt upper envelope can give genuine improvements.

## Charge every closed availability run

For one coordinate, a run is a maximal consecutive interval in A_t.
Recall A_T=V_T; no terminal reset is imposed. Write D_A for the number of
closed runs. The relaxed state/output matrix is

    H = [[1,qx],[x,1+qx]].

A forced-start closed run of length ell has count at most
c lambda^(ell-1) e. An unforced closed run has count

    U_ell = e_0^T H^ell (1,x)^T.

Its normalized count U_ell/lambda^ell decreases with ell. A direct proof
uses u=e_0^T H^ell=(a,b). The ratio b/a stays in [0,q+1]: its update
is (qx+(1+qx)r)/(1+xr), an increasing map preserving that interval.
For w=(1,x)^T,

    u(Hw-lambda w)
      = ax[-(q+1-qx)+(b/a)(1-x)] <= -a x^2 <= 0.

Therefore U_ell <= (1+qx^2) lambda^(ell-1). Define

    d = max(e, (1+qx^2)/lambda),  0<d<=1.

Multiplying the forced and unforced run bounds gives the single joint
inequality

    state sum <= c^f lambda^(L-f) d^D_A.              (1)

It covers all run lengths and every terminal endpoint. Do not multiply
(1) by the earlier forced-birth closure bonus: that would spend some
closures twice. At x=1, d=1. Closures do not supply an unconditional extra
q^(-D_A) rank factor: a closed run of ell states already has only ell
free state values. Their benefit here comes from the output-weight bound.

Two independent algebraic audits agree. Exact rational checks cover both
closed and free forced/unforced runs, including q=2^32.

## A new complete-history class certificate

The all-run lemma yields an additional exclusion without evaluating the
tagged transfer. Set x=2^-16, so

    lambda=2^16(1+2^-16+2^-32),    d<2^-15.

Also (1+2^-16+2^-32)^N<2^47. To verify the constant rationally, use
ln(1+u)<u and exp(32+1/2048)=[exp(1)exp(1/65536)]^32. The bounds
exp(1)<87/32 and exp(1/65536)<=65536/65535 give a product less than 11/4;
finally 11^32<2^111. The first exponential bound follows by summing its
series through degree five and bounding the tail by 7/4320.

Drop the birth saving c^f/lambda^f<=1 from (1). For each fixed footprint,
the conditional probability is strictly less than

    2^[16+16W-32R+16L-15D_A].

Condition on arbitrary route/partition geometry. There are at most 2^N
footprints; no additional state, support, row-set, or tag union remains.
Consequently the whole-history class

    32R-16L+15D_A >= 16W+N+92 = 17196732             (2)

has failure probability below 2^-76. Average over the fixed sampled
geometry afterward. Every terminal endpoint remains included. This is
a J-field class theorem: its particular class label need not be preserved
by extension-field basis projection. A completed proof covering every
J-word would transfer the distance guarantee in the usual way.

The integer remainder satisfies 32R-16L+15D_A<=C=17196731. For eta>=1,
multiply its canonical count by eta^(C-32R+16L-15D_A)>=1. The transfer
substitution is

    m,b -> eta^16(m,b),    alpha -> eta^-15 alpha,
    rank weight a -> eta^32 a,
    global factor eta^(C-4096h).

Apply this before relaxing footprints or rounding tags. The helper uses
`closuregap=log2(eta)`. Two independent audits confirm (2), its single
footprint union, free terminal, and this substitution. Exact rational
constant checks pass in `quartet_rank_boxes.py`. It is an additional
class certificate; the full remaining-event certificate is still open.

## Retain a coarse random subset of closure obligations

After taking a fractional moment theta, let alpha=d^theta. Independently
tag each position in A_(t-1) with probability 1-alpha, and require every
tag to remain in A_t. Conditional on the complete canonical history,
the probability of satisfying these requirements is exactly alpha^D_A.
Tags are proof variables only; the encoder and its randomness do not change.

To keep the transfer small, choose fixed levels J containing zero. If s
tags were sampled, retain a uniformly selected j(s)-subset, where
j(s)=max{j in J:j<=s}. Requiring only retained tags enlarges the probability
pointwise. Given ell available positions, the retained count h has law

    w_h(ell) = Pr[j(Bin(ell,1-alpha))=h].

Conditional on its count, the retained subset is uniform. Current birth
labels F and carried tag labels G are disjoint. Both must be available
next, but G receives bulk weight b=lambda^theta rather than birth weight
m=c^theta. Tags have no required survival beyond the free terminal.

The reduced states are (i,g), with 0<=i<=128, g in J and g<=4i.
J={0,32} gives 250 states; J={0} recovers the 129-state transfer exactly.
For c_(k,ell)=[s^ell]((1+s)^4-1)^k, a=mp chi and o=4i-g,

    C_((i,g),(k,h)) = binom(128,i)(r chi)^(4i)
      / [binom(512,4i) binom(4i,g)]
      * sum_(ell>=g) c_(k,ell) w_h(ell) binom(ell,g) b^g
        [z^o](a+bz)^(ell-g)(1+z)^(512-4k).

All terms are nonnegative. This factors as C=VU with U independent of
the routing-quota tilt p. To justify entrywise quota minima, insert an
eligible-label counter y into V. Apart from its binomial factors, V then
contains

    (p chi y)^f [chi(py+(1-p)/a)]^(4i)
      [1-p+py]^(512-4i-f),

where a here denotes the rank weight q^theta, not the coefficient
substitution mp chi used in the preceding matrix formula. Extract y^h
across all 16 V factors in a region. Dividing by
Pr[Bin(8192,p)=h] cancels p^h(1-p)^(8192-h), leaving the same conditioned
reduced operator for every p. Positivity bounds its nonzero-history part
by each quota-adjusted unconditioned operator. Their entrywise minimum
is therefore valid in the common basis. The last region has 16 V factors
and 15 U factors, ending in the free terminal vector. Independence of U
alone would not justify moving an entrywise minimum through a factor.

The normalization follows from

    multinom(512;f,g,4i-g)
      = binom(512,4i) binom(4i,g) binom(512-4i,f).

The terminal vector is independent of g:

    t_(i,g)=binom(128,i)(r chi b)^(4i)(1+mp chi)^(512-4i).

The full transfer is e_(0,0) C^4095 t. The zero operator remains
e_(0,0)e_(0,0)^T, so the existing positive regional and occupied-region
recurrences apply. Establish these statements for canonical histories
before adding the old positive fully-forced-quartet relaxation.

`quartet_available_tags.py` implements the factorization through positive
Taylor convolutions. Independent rational partition enumeration checks
308 matrix entries, 24 composed paths and 34 terminal entries. Additional
checks cover 135 thinning inequalities, 739 uniform-subset probabilities,
13 exact old-transfer entries and 384 run cases. All pass on Peach.

## Sharper tag rounding with the same latent counts

Flooring the sampled tag count can discard much of the closure penalty.
An audited alternative uses barycentric rounding to adjacent permitted
counts. For ell available positions and D closures, the probability that
a uniform g-subset survives is

    P(g)=binom(ell-D,g)/binom(ell,g)
        =binom(ell-g,D)/binom(ell,D).

This sequence decreases and is discretely convex. For D>=2 its second
difference is binom(ell-g-2,D-2)/binom(ell,D)>=0; D=0,1 give a constant
or linear sequence. Given a sampled count s, choose adjacent a<=s<=b
in J with b<=ell, then round to a,b with mean s. Convexity gives

    P(s) <= E[P(rounded s)] <= P(a).

If there is no upper level at most ell, use the floor. After choosing the
rounded count, select a uniform subset from ALL ell positions. It is not
a subset of the original Bernoulli tags when rounding upward. Averaging
over the sampled count proves that this gate is at least alpha^D and no
greater than the floor gate. The same positive transfer and quota proof
apply. Matrix entries need not decrease: mass moves between tag levels.

The implementation retains both modes. Both pass independent rational
partition comparisons: 616 matrix entries, 48 composed paths, 68 terminal
entries, and 26 old-transfer identities. There are also 1794 rounded
probability comparisons and 2184 exact gate-sandwich checks.

At h=2048, theta=1/16 and log2 z=-7/16, nine quota tilts and levels
{0,8,16,32,64,128} give FLOAT log2 upper bounds +111544.715 (floor) and
+105850.123 (barycentric), compared with +117834 for the prior forced-run
control. These are still failed upper bounds, not certified margins.
Changing log2 z to -5/8,-7/8,-1 at five quota tilts gives respectively
+190358,+397421,+507512. A stronger output penalty alone does not repair
the tested count. A finer level grid is evaluated separately.

This checkpoint is not a full-code certificate. Jobs are serial, bounded,
single-threaded and file-free on Peach. No encoder, construction, paper,
benchmark, commit, or push changes accompany these proof diagnostics.
