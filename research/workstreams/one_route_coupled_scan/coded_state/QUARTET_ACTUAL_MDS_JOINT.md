# Keep outer-zero constraints inside the global actual-state count

The construction and 45%/40-bit target are unchanged. This is a global
support enumeration, using the actual-state statistics in
[the actual-state checkpoint](QUARTET_ACTUAL_STATE_PROFILES.md).
It does not resume the local seven-type analysis. The full certificate
remains open.

The separate routing and support bounds permit combinations that the
accepted quartet maps cannot realize. We can exclude some of those
combinations before counting. This note specifies that joint count;
it does not claim that the additional restriction suffices for closure.

## A necessary condition on each live quartet

Fix a complete candidate actual-state support history S and a candidate
eligible outer-row set H. In one quartet let O be the incoming state
support, S_new the outgoing support, and H_local the eligible positions
after routing. Write r=|O| and v=|S_new|. The input is zero outside H.

If r=0, the new state equals the input, so S_new is contained in H_local.
If r>0, the new state outside H_local equals A Q_old. Superregularity
implies that a nonzero r-supported vector has at most r-1 zero coordinates
in its image. Therefore

    |H_local^c \ S_new| <= r-1.                          (1)

This applies also at the final update. It needs no availability fixed
point. It removes only transitions with r+v<=4: when r+v>=5 there are
already at most r-1 new zero positions in the whole quartet. A probe
must determine whether such transitions matter in the remaining bound.

## A positive quartet polynomial

For an auxiliary Bernoulli model, label every position eligible with
probability p, independently. Exact regional quotas will be restored
below; this does not change the routing distribution in the theorem.

Use positive markers s for new state positions, t for actual coordinate
deaths, z for the old-hit union G, nu for forced births f, and eta for
covered positions outside H. Put

    A = p+(1-p)eta,
    D_(r,v)(t) = [y^v](t+y)^r(1+y)^(4-r),
    T_(r,v) = sum_(j=0)^min(4-v,r-1)
                binom(4-v,j) p^(4-v-j) ((1-p)eta)^j.

For r>0 the entire quartet is covered. New nonzero positions permit
either label, contributing A^v. Among its new zero positions, at most
r-1 may be outside H by (1); T_(r,v) enforces this condition. For r=0,
every new nonzero position is a forced birth in H. Thus the polynomial

    P(U,V) = (1+p nu s V)^4
       + z^4 sum_(r=1)^4 binom(4,r) U^r
                sum_(v=0)^4 D_(r,v)(t) (sV)^v A^v T_(r,v)    (2)

has nonnegative coefficients. The variables U,V count old and new
support size. No accepted-map probability is multiplied into (2): it
retains a necessary condition of the already accepted setup.

For n=512, define the physical count-basis matrix

    K_(u,v) = [U^u V^v] P(U,V)^128 / binom(512,u).        (3)

The denominator averages the uniform quartet partition of a fixed
u-position old support. The outgoing masks and Bernoulli labels are
summed. Starting at the empty state and retaining every free endpoint
gives M=e_0 K^4096 1. There is no expansion after the final state.

## Restore quotas and the global field bound

For a fixed profile (h,R,Ls,f,F), put G=F-f and C_out=R-128h. In every
8192-position region the true H-label set is a uniform h-subset. It is
the Bernoulli-p distribution conditioned on its exact count, whose
probability is

    pi_h(p) = binom(8192,h) p^h (1-p)^(8192-h).

The 256 conditioning events are independent. Dropping their indicators
from a positive weighted sum costs at most pi_h(p)^(-256). Union over
the candidate H sets contributes binom(8192,h) once.

Choose 0<x<=1 and 0<=theta<=1. Put lambda=1+qx and
t=(qx^2/lambda)^theta, where q=2^32. The actual-state field lemma and
its cap at one give the profile bound

    binom(8192,h) pi_h(p)^(-256)
      s^(-Ls) z^(-G) nu^(-f) eta^(-C_out) M
      * [(gamma/q) q^(-R) x^(-W) lambda^Ls]^theta.      (4)

The terminal factor (qx/lambda)^|S_T| is at most one and has been dropped.
Apply the field bound and (1) on the actual eligible-row class first;
then drop additional row-class restrictions in the positive sum. Fixing
S and H before averaging avoids conditioning partitions on a support
chosen from their realized values.

The joint label/support count in (4) REPLACES the separate routing ratio.
Multiplying that ratio into (4) would count the same routing randomness
twice without independence. Likewise, eta is only a counter here: q^(-R)
occurs once in the external field factor. A separate probability budget
is required when excluding previously certified classes before (4).

At any fixed p in (0,1), (4) is valid. Different p or marker choices may
be minimized within a fixed profile. A finite numerical marker search
is a diagnostic until its chosen bound is verified with directed or
exact arithmetic; it never proves coverage of unsampled profiles.

## Bounded implementation check

The intended implementation first compares (2) with explicit labelled
(O,S_new,H_local) enumeration and with width-four/eight partition sums.
The full-width probe then tests the recorded low-birth profile
(h,R,Ls,f,F)=(4096,786432,786432,4,1048576).

Two meaningful label-marker choices are p=1/2, eta=1/4, and p=1/4,
eta=1/9. The latter matches the profile's uncovered eligible fraction
and its covered outside-H fraction in the unconstrained Bernoulli model.
The MDS restriction may change the optimizer. No performance benchmark
or construction change is part of this test.

The local condition, polynomial, exact-quota prefactor and probability
ordering have passed independent mathematical audit. Source checks and
the bounded numerical result are to be recorded after execution.

## Later high-birth diagnostic with independent geometric markers

The low-birth controls below do not settle the later high-birth remainder.
[quartet_actual_profile_potential.py](quartet_actual_profile_potential.py)
therefore reuses the same physical kernel with independent markers and
theta=0. For the positive potential alpha^u, one quartic coefficient pass
computes all 513 ratios
[X^u]P(X,alpha)^128 / (binom(512,u) alpha^u).
It retains actual live hits, forced births, deaths, and the accepted-MDS
support condition. It includes the free endpoint and exactly one H union
and regional-quota conditioning cost. It uses no separate route factor.

The independent source audit passed. Exact checks cover 18 physical
partition rows, 28 coefficient rows, 24 derivative coordinates, 28 free
endpoint inequalities, and four finite-box dual comparisons.

At (h,R,L,f,F)=(4096,819200,838861,65536,1081344), D<=104857,
and p=1/2, a bounded 256-evaluation control returned a FLOAT exponent
of +869385, down from its initial +1074107. No barrier stage was centered;
the dual lower diagnostic was -2252280. Thus this run establishes neither
closure nor optimizer optimality. It does not test the stronger count of
run-saturated supports in [the saturation reduction](QUARTET_RUN_SATURATION.md).
The latter removes within-run pattern multiplicity while retaining hits.
No additional parameter grid or code-construction change accompanies this
diagnostic.

## Whole eligible-row classes without a profile union

There is a useful specialization of (4). Fix h and a positive birth cap
f0. Choose x=2^-16, theta=1/16 and the rational upper markers

    s=2+2^-19,  t=1/2,  eta=1/4,  z=1,  0<nu<=1.

Concavity gives 65537^(1/16)<=2+2^-19, and
65537^(-1/16)<1/2. Thus the state marker dominates lambda^theta, and
the death marker dominates the fractional field death factor. On the
actual eligible-row class, R=128h+C_out, so the outside marker cancels
the variable part of the field rank. For all histories with 1<=f<=f0,
the entire class has failure bounded by

    binom(8192,h) 2^(W-256h-1) nu^(-f0) M_h.          (5)

Here M_h already restores the exact regional quotas, or bounds the
resulting quota-conditioned operator entrywise. The factor 2^-1 bounds
(gamma/q)^(1/16). Neither R,Ls,F,D nor f is separately union-bounded
in (5); the moment has already summed their histories. Only the chosen
h and birth-cap classes require accounting when several bounds are
combined. The all-zero path can be removed positively, and every
nonzero outer word occupies at least 129 regions.

This cancellation is independently audited. It reduces the next
coverage problem to whole h/birth-cap classes when the specialized
moment is strong enough. It is not itself a numerical certificate.

## A quota-approximation obstruction, before further tuning

The unrestricted Bernoulli moment includes the path S_1=...=S_T=[512].
Its exact contribution is

    (p nu s)^512 [s z (p+(1-p)eta)]^(512(T-1)).        (6)

There are no deaths and the MDS restriction is vacuous on this path.
At p=1/2, eta=1/4, s=2, z=1, the bulk factor is 1.25^(N-512).
Even nu=2^-18 leaves the primary profile's FLOAT lower bound on this
particular upper-count method at about +570897. The path is not a
codeword or a witness to the target profile. It proves that this scalar
quota relaxation cannot close that profile at these markers, so no
full-width run of this doomed baseline is needed.

Restore quotas after each 16-update region instead. Let Z=e_0 e_0^tr
be its exact empty-footprint operator. For each p,

    Q_true <= (K_p^16-Z)/pi_h(p)

entrywise on the same physical population basis. Consequently the
entrywise minimum over several p values is also a valid upper bound
Q_min, and the global region operator is Z+Q_min. The external
pi_h(p)^(-256) factor must then be removed. Q_min(0,0) is generally
positive: a nonempty region can return to the empty state. It must not
be erased as if it were the wholly empty region.

For the primary h=4096, eta=1/4 comparison, p in {1/5,1/3,1/2}
matches fully covered, half-covered and empty regions respectively in
the independent-label approximation. This finite family is a justified
control for the obstruction in (6), not an assertion of optimality.

The base implementation passes 150 exact local coefficient checks,
106 exact matrix entries, 23 endpoint checks, 225 common-basis quota
conditioning entries, 50 positive nonzero-power entries and 81 floating
comparisons against exact matrices. Independent static audit passes.
The regional-minimum implementation also passes 100 exact comparisons,
25 floating regional comparisons against exact rationals, and 25 checks
of the explicit positive-part power. Independent source audit confirms
that composition precedes quota adjustment and the common-basis minimum.

At the primary profile, with x=2^-16, theta=1/16, s=2, z=1,
nu=2^-18 and eta=1/4, the three-tilt regional probe gives a FLOAT
log2 moment bound of 771645.686 and a combined profile exponent of
+675044.004. This does not certify that profile. In particular, the
regional upper kernel has log2 Q_min(0,0)=2029.557, so its repeated
zero-boundary excursions alone give a moment lower bound exceeding
519566. This is a floor on the relaxed counting method, not on code
failure. Each such excursion needs a forced birth; repeating it in
256 regions is incompatible with the target f=4. The finite nu marker
does not adequately suppress that incompatible family.

The stronger theta=1/4 and theta=1/2 controls initially stopped at
strict positive-arithmetic range guards. They yielded no bound and
must not be treated as failures or successes of those mathematical
weights. Any relaxed numerical range handling must add positive mass,
never silently discard underflowed terms.

For the whole-h specialization, M_h is nondecreasing in h when
0<eta<=1. Couple uniform regional h-subsets by adding eligible slots.
This relaxes both the birth constraint and (1), and increases eta^C_out.
The coupling also preserves restrictions depending only on the state
history and partitions. Therefore an h interval [a,b] can use M_b.
The remaining factor c_h=binom(8192,h)2^(-256h) satisfies
c_(h+1)/c_h<2^-243, so summing that interval costs less than one bit
above c_a M_b. This is an audited coverage reduction; the numerical
moment is not yet small enough to apply it as a new certificate.

## What the stronger-weight controls actually show

The helper now has optional comparison-operator diagnostics and an
explicit FLOAT range control. The latter rounds tiny positive entries
up to a relative floor of 2^-4000 and uses logarithms when combining
regional scales. It adds mass in the mathematical comparison, rather
than dropping underflowed paths. Ordinary floating arithmetic is still
not a directed certificate. Exact toy-reference checks and 131 floating
majorant comparisons pass; the default strict range guards remain.

The following controls use the same primary profile, x=2^-16, nu=2^-18,
and no occupied-region marker. Each minimizes over three regional label
tilts. All values are FLOAT upper-bound diagnostics, not failure rates.

| theta | s | z | eta | label tilts | profile exponent |
|---|---|---|---|---|---:|
| 1/16 | 2 | 1 | 1/4 | 1/5, 1/3, 1/2 | +675044 |
| 1/8 | 4 | 1 | 1/16 | 1/17, 1/5, 1/2 | +695800 |
| 1/4 | 16 | 1 | 1/256 | 1/257, 1/17, 1/2 | +1883688 |
| 1/4 | 4/3 | 32/3 | 1/256 | 1/257, 1/17, 1/2 | +1108803 |

The last two use the positive floor. Their large bounds do not establish
that the corresponding unfloored or exact-quota calculations fail.
The last control is motivated by the population-3/4 stationary vector
of [[1,4/3],[1/16,4/3]], whose eigenvalue is 3/2, together with the
full-coverage quota factor. It is a targeted comparison, not a claimed
optimizer for the quartet process. The actual t is slightly below 1/16.

In the natural theta=1/4 comparison, the approximate left/right
eigenvector diagnostic concentrates near population 510. In the last
comparison it concentrates near populations zero and one. Both use
the p=1/17 upper matrix at every entry. These are properties of the
comparison matrices, not distributions of realized code histories.
Thus changing weights shifts the overcount between dense paths and
short regional excursions; it does not close the remainder.

An [exact four-birth persistent family](QUARTET_LOW_BIRTH_MOMENT_OBSTRUCTION.md)
proves that an exact birth count alone cannot repair the first row's
envelope. It has nearly full coverage, outside the actual remainder.
Also, [two tempting output-MDS refinements](QUARTET_MDS_REFINEMENT_LIMIT.md)
are ineffective at that row's weights. These results justify retaining
the omitted global restrictions before expanding the numerical search.

## Insert the remainder restrictions before the positive relaxation

Bounds for previously certified classes cannot simply be subtracted
from this moment. Their histories must be excluded or penalized inside
the sum. For example, the remaining coverage satisfies F<=F_max=1736703.
For 0<chi<=1, multiply its indicator by

    chi^(F-F_max) >= 1.

This replaces z by z chi and nu by nu chi inside the kernel, with the
external factor chi^(-F_max). Other external marker corrections keep
their original definitions. At a fixed exact F this is already an
ordinary coverage tilt; the new point is to retain it when covering
whole h/birth-cap classes rather than dropping F.

The actual-state census complements also have positive multipliers.
For mu>=1, the complement

    32(R-Ls)-2jD <= c_j+44

permits s -> s mu^32, eta -> eta mu^-32, t -> t mu^(2j),
with external factor mu^(c_j+44-4096h). The closure complement

    32R-16Ls+(16-2j)D <= 16W+c_j+91

instead permits s -> s mu^16, eta -> eta mu^-32,
t -> t mu^(2j-16), with external factor
mu^(16W+c_j+91-4096h). These follow by substituting R=128h+C_out;
they use actual-state statistics and the earlier exact constants.
The multipliers must be identical across label tilts before taking a
regional minimum. Prior-class failure allowances remain separate.

These are valid next interfaces, not evaluated class certificates.
Retaining D also permits choosing the output marker for its actual
remaining range; exact binomial extraction alone is not an exponential
improvement over a correctly chosen death-conditioned Chernoff bound.
The full 45%/40-bit goal remains open.
