# Cover actual deaths within the live-selector quotient

The live-selector quotient groups amplitudes before applying a fractional
moment. Its volume-only bounds can still count support histories whose
deaths are incompatible with the target output weight. This note retains
actual deaths through two global interfaces: a field bound for each death
interval, and a relation between actual deaths and closed allowed runs.
Neither interface enumerates actual support masks outside the quotient.

The construction and the 45%/40-bit target are unchanged. The results below
are counting interfaces and a proposed finite-cover strategy. They do not
close the remaining multirow class. Use the field-rank proof order from
[the live-selector quotient](QUARTET_LIVE_SELECTOR_QUOTIENT.md) and the
integrated birth count from
[the selector-hit note](QUARTET_SELECTOR_HIT_MOMENT.md).

## Actual deaths and holes in the allowed support

Write q=2^32, n=512, T=4096, N=nT, and W=943718. For an actual state history,
let S_t=supp(Q_t), with S_0 empty. Its actual volume, deaths, and terminal
population are

    L=sum_t |S_t|,
    D=sum_(t=1)^T |S_(t-1)\S_t|,
    a_T=|S_T|.

There is no death after update T. The output is
Y_t=Q_t+delta Q_(t-1), with delta nonzero.

Let U_t be the union of the current quartets meeting S_(t-1), and let
B_t=S_t\U_t. Set U_1=empty and use the formal endpoint U_(T+1)=[n]. Then

    B_t subset U_(t+1)\U_t,
    A_t=(U_t intersect U_(t+1)) union B_t,
    S_t subset A_t.

The allowed support A is a relaxation of the actual support S. Write

    f=sum_t |B_t|,             G=sum_t |U_t|,       F=G+f,
    A0=sum_(t=1)^T |U_t intersect U_(t+1)|,
    m=L-f,                    H_A=A0+f-L=A0-m.

Here H_A counts allowed positions where the actual state is zero. Let C_A
be the number of maximal allowed coordinate runs ending before T. Thus

    C_A=sum_(t=1)^(T-1) |A_t\A_(t+1)|
       =sum_(t=1)^(T-1) |A_t\U_(t+2)|.                 (1)

For the second equality, A_t is contained in U_(t+1), whereas B_(t+1) is
disjoint from U_(t+1). The formal full endpoint makes the final summand
zero. It does not introduce a closing output.

The two-sided relation is

    |D-C_A| <= H_A.                                    (2)

To prove D<=C_A+H_A, consider an actual death. If it occurs when the allowed
run closes, charge it to that run. Otherwise its following allowed
position is a hole. Distinct remaining deaths use distinct holes.
Conversely, the last position of each closed allowed run is either
actually nonzero or a hole. A nonzero last position produces an actual
death at the next update. Distinct closed runs use distinct endpoints.
This proves C_A<=D+H_A.

Consequently a death interval d0<=D<=d1 imposes both

    A0-C_A >= m-d1,       A0+C_A >= m+d0.              (3)

The first inequality retains a small-death restriction. The second can
retain a large-death restriction. Both concern allowed runs, so they can
be inserted before enlarging the amplitude count.

The first inequality can be strengthened if allowed-run lengths are
retained. For an integer k>=1, define

    R_k=sum_closed_allowed_runs min(k,run_length).

Every closed allowed run containing an actual nonzero state contributes
at least one actual death. An entirely zero allowed run consumes its full
length in holes. Summing these alternatives gives

    kD+H_A >= R_k,
    D<=d1  implies  A0-R_k >= m-kd1.

The choice k=1 recovers C_A. Larger k records that removing a long allowed
run costs more than one hole. The current nested population helpers do
not retain these capped run lengths; this is a possible source of new
information if improving their potential does not suffice.

## Output weights at fixed actual support

Fix routing, quartet partitions, a candidate pair (U,B), and an actual
support S contained in A. Retain the exact-live conditions when applying
the field-rank argument. Let R be the outer parity rank of the footprint
U union B. The accepted-field density loss is at most gamma<=1024/989.

The actual support has D+a_T coordinate runs. Every B position starts a
distinct such run and has a prescribed nonzero value. For 0<x<=1, the
weighted number of amplitude assignments at this support, for any fixed
nonzero vector on B, is exactly

    (q-1)^(D+a_T-f) x^(2D+a_T)
       * [1+(q-2)x]^(L-D-a_T).                         (4)

Each unprescribed run start has q-1 nonzero choices. A closed run has two
mandatory nonzero boundary outputs; an open run has one. At a transition
between two nonzero states, exactly one next value gives zero output and
q-2 values give nonzero output. Scaling each forced run transports any
prescribed nonzero start to another, without changing output weight.
Thus (4) is independent of the nonzero values on B.

The residual RS equations admit at most q^(f-R_B) compatible birth
vectors, and the raw-field probability on live coordinates is q^(-R_U),
where R_U+R_B=R. Apply these factors while the live conditions are exact.
Then enlarge the amplitude family. After summing compatible birth
vectors, projective normalization divides by q-1. The accepted setup
introduces gamma only once.

Define

    lambda_x=1+qx,
    t_x=qx^2/lambda_x,
    B_x=(gamma/(q-1)) q^(-R) x^(-W) lambda_x^L.

Replacing q-1 by q and 1+(q-2)x by lambda_x in (4) gives the fixed-support
contribution bound

    B_x t_x^D (qx/lambda_x)^a_T <= B_x t_x^D.          (5)

Both omitted terminal factors and t_x lie in (0,1]. The factor x^(-W)
converts the output enumerator into an upper bound for weight at most W.

For fixed (U,B), there are at most C(A0,m) supports with actual volume L.
First sum the supports in the interval [d0,d1]; only then apply a
fractional moment. Their field contribution is bounded by

    min(1, B_x t_x^d0 C(A0,m)).                        (6)

Equation (6) remains uniform in the compatible nonzero birth vector.
It does not add a support union outside the quotient. In particular,
the binomial factor is raised to theta when 0<theta<=1.

## A finite death-cover bound

Let h be the number of eligible footprint rows and put

    C_H=F-R+128h,
    Q_route=min(1, C(G,C_H-f) C(N-F,256h-C_H)
                    / C(8192,h)^255).

This route bound applies uniformly to each fixed selector and birth
configuration at the specified profile. It already includes the union
over eligible row sets. Apply it before dropping the fixed G, f, and C_H
restrictions from the positive moment.

Define Xi(beta,d,w,z) as follows. Sum over candidate selector histories U
and birth sets B, and average over the random quartet partitions. Require
B_t subset U_(t+1)\U_t. Each selected new quartet must meet U_t union B_t;
the formal full endpoint has no such requirement. Give a history weight

    z^G beta^A0 d^C_A w^f.

The positive moment Xi includes all G, f, and allowed volumes. Its formal
endpoint includes arbitrary births outside U_T.

Choose b,k_minus,k_plus>=1, y,w,z>0, and 0<theta<=1. Set

    beta=b(1+y)^theta k_minus k_plus,
    d=k_plus/k_minus.

The binomial estimate and the necessary conditions are

    C(A0,m)^theta <= y^(-theta m)(1+y)^(theta A0),
    1_(A0>=m) <= b^(A0-m),
    1_(A0-C_A>=m-d1) <= k_minus^(A0-C_A-m+d1),
    1_(A0+C_A>=m+d0) <= k_plus^(A0+C_A-m-d0).

Here each indicator inequality is used on the corresponding restricted
family. Combining them with (6) gives the whole death-interval bound

    P_profile,[d0,d1]
      <= Q_route z^(-G) w^(-f)
         * B_x^theta t_x^(theta d0)
         * b^(-m) y^(-theta m)
         * k_minus^(-(m-d1)) k_plus^(-(m+d0))
         * Xi(beta,d,w,z).                             (7)

For theta=0, omit the field and y factors and set
beta=b k_minus k_plus. This is a geometry bound on the same interval.
There is no separate birth binomial or separate empty-quartet hit factor:
Xi counts the birth choices directly.

Low-death intervals can use k_minus>1. High-death intervals gain the
factor t_x^(theta d0), and may also use k_plus>1. Taking a minimum over
several marker choices for one interval costs no union factor. Disjoint
intervals covering the allowed death range must be summed explicitly.

For example, if at most 2^r intervals each have a certified bound
2^(-172-r), their profile sum is at most 2^-172. The existing strict
2^96 union over (h,R,L,f,F) then gives a qualifying-class bound below
2^-76. A direct certified sum can avoid this uniform allocation.
Previously closed classes retain their separate failure allocations.
Neither sampled profiles nor a finite marker family alone proves a full
profile cover.

## A scalar interface for closed runs

Xi has a nested state (U_t,A_(t-1)), with A_(t-1) subset U_t. A scalar
supersolution avoids storing all pairs of population counts. Choose
positive potential parameters alpha and rho, and define

    V(U,Aprev)=alpha^|U| rho^|Aprev|.

In a fresh quartet let r=|U intersect quartet| and
a=|Aprev intersect quartet|, so 0<=a<=r<=4. If the new quartet is
unselected, a allowed coordinates close, giving d^a. If it is selected,
all r old-live coordinates belong to the new allowed support. The other
coordinates may be births; at least one birth is required when r=0.
The local action, including the future potential, is

    h(a,r)=d^a
       + z^4 alpha^4 (beta rho)^r
           * [(1+w rho)^(4-r)-1_(r=0)].                (8)

If a positive Lambda satisfies

    h(a,r) <= Lambda alpha^r rho^a
        for every 0<=a<=r<=4,                          (9)

then each partition sends V to at most Lambda^(n/4)V. No average over
partitions is needed for this implication. The initial potential is one.
The terminal weight is beta^u(1+w)^(n-u), where u=|U_T|. Its ratio to V
is at most

    E_end=max(1+w, (beta/alpha)max(1,rho^(-1)))^n.

Consequently

    Xi(beta,d,w,z) <= Lambda^((n/4)(T-1)) E_end.       (10)

The endpoint includes B_T and has no final closure charge. For fixed r,
the normalized expression in (9) is a sum of two exponentials in a.
It is log-convex, so a=0 and a=r suffice: nine distinct cases in total.
Checking all fifteen cases is equally valid.

If d<=1 and rho>=1, a=0 always maximizes the normalized expression.
The closure marker then disappears from the bound. A useful small-death
control within this scalar family must use rho<1 and account for the
resulting endpoint cost. Equations (8)--(10) do not assert that such a
control has a favorable numerical exponent.

The integrated birth transition preserves a smaller population domain.
Every selected new quartet meets the new allowed support: it contains
either an old-live intersection or a required birth. Therefore every
reachable nested state satisfies

    |U|/4 <= |Aprev| <= |U|,

including the initial state (empty,empty). This condition concerns the
allowed support, so it survives amplitude enlargement. The fresh partition
is independent of the preceding partition; its averaged action still
depends only on these two populations.

A scalar multiplier can exploit this condition before averaging. For
0<t<=1, replace (9) by

    h(a,r) <= Lambda alpha^r rho^a t^(a-r/4).

After multiplication over all quartets, the extra factor is
t^(|Aprev|-|U|/4)<=1. The normalized local log action gains
(r/4-a)log2(t). Convexity and the nine-case reduction are preserved; the
terminal bound is unchanged. Choosing t=v^4 makes the new exponents
integers for rational verification. This is an additional valid interface,
not a recorded numerical result.

The regional capacity can be retained in the same potential. For c>=1,
let Lambda_0 and Lambda_1 satisfy (9) at z and zc, respectively. A region
of ell new-selector layers sends V to at most

    min(Lambda_1^((n/4)ell), c^h Lambda_0^((n/4)ell)) V.

Use ell=15 for the first region and ell=16 for each of the next 255.
Multiply the terminal factor and the external correction c^(-(C_H-f)).
Both comparisons must use the same alpha and rho.

The pointwise bound (9) replaces the actual quartet mixture by its worst
local case. If that replacement is too loose, average the fresh partition
before testing the same potential. Fix |U|=u and |Aprev|=a. Let X mark
Aprev, Y mark U\Aprev, and let the constant term mark positions outside U.
Define the positive polynomial

    P(X,Y)=(1+dX+Y)^4
       + z^4 alpha^4 [(1+w rho+beta rho(X+Y))^4-1].

The exact expected one-step action divided by the old potential is

    [X^a Y^(u-a)] P(X,Y)^(n/4)
       / [multinomial(n;a,u-a,n-u) alpha^u rho^a].     (11)

The first term of P describes an unselected quartet. The second describes
a selected quartet; subtracting one removes exactly the old-empty,
birth-empty selection. All coefficients remain nonnegative. The
multinomial denominator averages the locations of the two nested old
sets relative to the fresh random partition.

It suffices to bound (11) uniformly over u=0,4,...,512 and u/4<=a<=u.
There are 24,897 reachable population rows. The larger domain
0<=a<=u has 33,153 rows and remains a conservative comparison. Each
coefficient can use its own positive Chernoff variables X,Y, followed by
a maximum over the rows. The endpoint and regional comparison still use
the same alpha and rho. This verifies a scalar potential for the
population process; it does not require a transition matrix or a
refinement by local state types.

The local endpoint reduction in a does not extend to the averaged rows.
For a counterexample, take n=8, u=4, and a small positive rational epsilon.
Choose

    d=epsilon^8, beta=epsilon^(-8), alpha=1, rho=epsilon,
    w=(beta-1)/epsilon, z=epsilon/beta.

These are valid low-death markers with b=1 and k_minus=1/d. A selected
quartet with r>0 old-live coordinates has weight epsilon^(r+4), including
the future potential. At r=0 its weight is
epsilon^4(1-epsilon^32). Let R(4,a) be the exact normalized row in (11).
Then

    lim_(epsilon->0) R(4,a)/epsilon^4 = C(4,a)/35,
        a=1,2,3,4.

The leading terms select one quartet containing all a old allowed
coordinates and no other old-live coordinates. Their total coefficient is
2 C(4,4-a)/C(8,4): either of the two quartets may be selected. All other
terms have higher epsilon degree after normalization. Hence a=2 strictly exceeds
both reachable endpoints a=1 and a=4 for sufficiently small epsilon.
The parent verified this strict comparison at epsilon=1/8 using exact
Fraction arithmetic and the physical coefficient formula; no floating
decision was used.
All reachable a values require comparison unless another envelope is
proved.

Perfect partition averaging does not remove every source of large mass.
For any positive potential satisfying K V<=mu V on the reachable states,
three exact transition families imply

    mu >= 1,
    mu >= (z beta)^n,
    mu^2 >= [1+z^4((1+wd)^4-1)]^(n/4).

The first two are the empty and full-state self-loops. The third follows
all two-step excursions from the empty state back to itself. Each selected
first-step quartet requires a nonempty birth set; all those allowed
coordinates close on the second step. In fact the last right-hand side
is exactly the corresponding entry of K^2. These are lower bounds on the
relaxed operator's growth. They do not establish a route-compatible
support geometry or a codeword.

These floors give a parameter-free lower estimate for the geometry-only
uniform-step comparison when d1>=f. Put

    C=(n(T-1)-m)/8,       k=(G-m)/4,
    Z=z^4,               V=wd.

A convex combination of the full-state and two-step floors cancels the
overlap marker. The resulting profile exponent is at least

    log2 Q_route + C log2[1+Z((1+V)^4-1)]
       - k log2 Z - f log2 V.

The omitted death-marker term is (d1-f)log2(k_minus)>=0. The omitted
initial-potential times endpoint factor is also at least one. If
k<=f<=2k and C>=k, set n0=C-k, n1=2k-f, and n2=f-k. The variational
inequality for the positive polynomial 1+Z(4V+6V^2+4V^3+V^4) gives the
further lower estimate

    log2 Q_route + C log2 C - sum_(i=0)^2 n_i log2 n_i
       + 2n1 + n2 log2 6,

with 0 log2 0 interpreted as zero. This is a test for a possible
limitation of the uniform-step bound, not a positive lower bound on the
code's failure probability. No sign is asserted here. Regional-capacity
or exact-quota operators need their own analysis; this argument concerns
K followed by the external uniform route factor.

A short collection of positive potentials can address the loss from one
exponential. Fix common geometry markers and let
V_j(u,a)=alpha_j^u rho_j^a. Suppose physical row bounds F_j>=K V_j are
available. For c_j>=0 with at least one positive coefficient, the sum
V=sum_j c_j V_j satisfies

    (K V)(s)/V(s) <= sum_j c_j F_j(s)/sum_j c_j V_j(s).

For two strictly positive coefficients, the physical minimum
V=min(c_1 V_1,c_2 V_2) instead satisfies

    (K V)(s)/V(s)
       <= min(c_1 F_1(s),c_2 F_2(s))
            / min(c_1 V_1(s),c_2 V_2(s)).

In either case, maximize the row ratio over reachable populations and
multiply mu^(T-1) by V(0,0) and max_s terminal(s)/V(s). The initial factor
is sum_j c_j for the sum and min(c_1,c_2) for the minimum. These interfaces
reuse the averaged actions and require no transition matrix. A minimum
of physical functions is valid by positivity of K; an entrywise minimum
of coefficient matrices in an overlapping positive basis is not generally
valid.

## Which regimes require geometry

The field penalty cannot replace a necessary geometric exclusion. Fix
any support mask history S, not necessarily realizable as an exact state
support. Its supported state space has dimension L. Its outer parity
constraints have rank at most R, so the compatible subspace has dimension
at least L-R. The triangular map from states to outputs is injective and
has output support length at most L+D.

If L>R, the Singleton bound therefore supplies a nonzero compatible output
of weight at most

    (L+D)-(L-R)+1=R+D+1.

Thus a genuine mask geometry with

    L>R and D<=W-R-1                                   (12)

would force a bad codeword for every field map. The codeword may have a
smaller actual support than the mask used in this argument. Proving that
no such geometry exists is a structural task. A positive relaxed upper
bound is not evidence that a mask satisfying (12) exists.

The remaining diagnostic profile

    (h,R,L,f,F)=(4096,819200,838861,65536,1081344)

has L-R=19661 and D>=f-512=65024. For D<=W-L=104857,
every supported output already has weight at most W. The larger interval
D<=W-R-1=124517 is also structural by (12). Deaths beyond that threshold
can benefit from the grouped output penalty in (7).

A finite cover should split at these two thresholds, then refine only
intervals left unresolved by the combined bound. Rank surplus R>L is a
different regime, already addressed by grouped-rank classes. The earlier
large birth penalties were partly counting losses; integrated births
remove those losses but do not remove condition (12).

## Optional joint regional quota transfer

One can retain H labels together with selectors and births instead of
multiplying by Q_route. This is an alternative to the external route
bound in (7), not an additional factor.

Fix an eligible row set H of size h. Each region has a uniform h-subset
of its 8192 coordinate slots labelled by H. Introduce a positive formal
quota variable x, a covered-H marker a_H>0, and positive birth and volume
markers w,z. At a current update, set

    nu=1+x,       mu=1+a_H x,       wprime=w a_H x/nu.

A slot outside the current live set contributes nu if it is not a birth.
A current live slot contributes mu. A birth outside the live set must
have an H label, so its contribution relative to nu is wprime.

For old population 4i and new population 4j, define K_x by the row factor

    nu^n (z mu/nu)^(4i)

times the integrated-birth selector kernel with overlap beta and birth
marker wprime, with no new-volume marker. Its selected-quartet polynomial
is (1+wprime+beta X)^4-1. This transfer charges the CURRENT live volume.
The endpoint is

    t_x(i)=nu^n (z mu beta/nu)^(4i)
                  (1+wprime)^(n-4i).

Although this notation uses rational functions, every transfer entry is
a polynomial in x with nonnegative coefficients. Construct that formal
generating operator first. Extract the quota coefficient only after a
whole region; current labels are not independent after conditioning on
the quota.

For fixed marker values beta>=1 and w,z,a_H>0, the exact unnormalized
regional operators are [x^h]K_x^16. They admit entrywise upper bounds

    B=min_(x in finite positive family) x^(-h) K_x^16.

The final vector admits the componentwise upper bound

    t=min_(x in finite positive family) x^(-h) K_x^15 t_x.

The common population basis permits these minima after composition.
There are 255 ordinary regions, followed by 15 transitions and the last
endpoint. The first 255 blocks charge updates 1 through 4080; the final
block charges updates 4081 through 4096.

The geometry-only profile bound is

    C(8192,h)^(-255) a_H^(-C_H) w^(-f) z^(-G)
       * beta^(-(L-f)) * e_0 B^255 t.                 (13)

The denominator accounts for 256 conditional quotas and one H-set union.
There is no further Q_route. Equation (13) still relaxes individual row
degrees, but retains regional quotas jointly with birth feasibility and
selector geometry.

The parent's bounded controls were positive: approximately +38809 bits
at h=2048, +157151 at h=4096, and +125804 at h=3584. They did not improve
the corresponding external-route controls. The h=4 proposal hit a zero
probability parameter and a division guard, so it produced no bound.
These are diagnostic outcomes, not exact certificates or codeword
witnesses. Further quota tuning is parked until closed runs are retained.

## Rank and row information that the refinement must preserve

Fix the field maps and an allowed state mask S. Let E_S embed its L state
coordinates, let T_field map states to inputs, and let H_RS be the outer
parity-check matrix. The columns of T_field E_S are supported on V(S), so
T_field E_S=E_V C for a suitable matrix C. Therefore

    rank(H_RS T_field E_S)
       <= rank(H_RS E_V)
        = sum_rows min(128,|V intersect row|)=R.

The rank is also at most L. Elimination or a change of basis cannot create
more than R independent outer parity equations on this state space.
Additional support-feasibility or probability information may still help;
it is not an increase in the rank of this linear system.

Exact row classification is also nearly automatic at the primary profile
`(4096,819200,838861,65536,1081344)`. Here C_H=786432 and C_outside=294912.
Fix any H-slot masks and regional covered-slot counts c_r,o_r with these
totals. Conditional on these masks, the independent within-group row
permutations give each fixed H row independent Bernoulli incidences across
regions, with mean degree 192. An outside row has mean degree 72.

Concavity of log(1-p+pz) shows that equal regional densities maximize
each such exponential moment. Chernoff therefore gives, uniformly in the
fixed regional geometry,

    Pr[H-row degree<=128] <= (3/4)^128 < 2^(-48),
    Pr[outside-row degree>=129]
       <= Pr[outside-row degree>=128]
       <= (207/256)^128 < 2^(-38).

The markers are z=1/3 and z=23/9, respectively. For the rational bounds,
(3/4)^8<1/8, 207^4<7*2^28, and 7^16<2^45 suffice. A union over 4096 rows
in each group gives error below 2^(-36)+2^(-26)<2^(-25). Thus imposing
the correct row classification retains at least 1-2^(-25) of any
nonnegative geometric count whose weights are independent of these row
identities. Its loss is less than 2^(-24) bits. This ratio is not asserted
after imposing amplitude equations that depend on the row identities.

Regional footprint capacity is different: it constrains which slot
geometries can realize the profile. Let v_r be the full footprint volume
in region r, including births. Exactly h slots in each region belong to H,
so every matching route satisfies

    C_H <= sum_r min(h,v_r),
    sum_r (v_r-h)_+ <= R-128h.                         (14)

This is stronger than the live-volume condition that treats all births
as a separate allowance. It also implies the top-128-regions constraint:
for any set I of 128 regions,

    sum_(r in I) v_r
       <= 128h+sum_r(v_r-h)_+ <= R.

Unlike a general threshold representation of the top-128 sum, (14) uses
the fixed profile threshold h. No union over thresholds is needed.
The external route bound remains compatible:

    Pr[required route profile | fixed U,B]
       <= 1_(sum_r(v_r-h)_+<=R-128h) Q_route.

This uses a necessary condition and a uniform conditional bound, not an
independence assumption. There is no additional H-set count.

For eta>=0, the positive regional factor is

    2^(-eta(v_r-h)_+)=min(1,2^(eta(h-v_r))).

An upper regional operator takes the physical entrywise minimum of an
ordinary block and its volume-tilted block times 2^(eta h), before
composing regions. Both live and birth markers must be tilted by
2^(-eta), since v_r includes both. The external correction is
2^(eta(R-128h)). With current-update charging, use 255 blocks of 16
transitions followed by 15 transitions and the free endpoint. The exact
regional quota operator already contains (14); a scalar comparison can
lose its dependence on the regional history. The new 129-population
regional filter in
[quartet_regional_rank_cap.py](quartet_regional_rank_cap.py) retains this
necessary condition without changing the construction. Static source
review passed current-volume charging, positive birth coefficients,
terminal births, 16/15 powers, both tilted markers, physical minima, and
fixed-h versus arbitrary-threshold scope. Parent checks passed 13 current
kernel entries, five terminal entries, 72 regional comparisons,
34 two-region moments, and 4401 capacity identities.

The parent tested four bounded centers. At beta=3, w=1/2, and
z=2^(-1.01)/(3/2), eta=0 gave +148157.034 bits. This improves the earlier
+148463.722 by about 307 bits through the exact endpoint; every tested
positive eta from 1/32 through 2 was worse. The other controls were:

| log2 z | Best tested eta | FLOAT exponent |
| ---: | ---: | ---: |
| approximately -1.5649625 | 1/32 | +151836 |
| approximately -1.5449625 | 1/32 | +163721 |
| approximately -1.5249625 | 1/16 | +162661 |

At the last center, the unfiltered bound was approximately +199198 bits,
so the capacity filter reduced that bound. None of these four centers
improved on +148157.034. These are failed probability upper bounds,
without an optimizer or method-impossibility claim.

## Verification and next step

Independent derivation audits passed (2), (4)--(7), the scalar transition
and endpoint, and the formal-quota proof order. The positive averaged-row
formula (11) also passed a derivation audit. This note reports no new
numerical probability certificate.

[quartet_nested_selector_scalar.py](quartet_nested_selector_scalar.py)
implements the pointwise scalar bound for a low-death interval, with
integrated births and a bounded six-variable FLOAT search. Static source
review passed the affine local cases, endpoint, profile corrections,
smoothing error, and finite-box lower diagnostic. The parent-run checks
passed 45 exact local cases, 1593 physical nested-partition cases,
63 endpoint cases, 12 temporal moments, and six gradient/Hessian
coordinates.

The parent ran two bounded controls at
`(h,R,L,f,F)=(4096,819200,838861,65536,1081344)`:

| Death upper bound | FLOAT exponent | Evaluations | Iterations | Finite-box lower diagnostic |
| ---: | ---: | ---: | ---: | ---: |
| 104857 | +177477.46 | 979 | 117 | approximately -1901157 |
| 124517 | +191951.90 | 744 | 100 | approximately -2973278 |

Both controls failed to certify their death intervals. The large gaps to
the lower diagnostics establish neither optimizer optimality nor an
impossibility result for the scalar family.

The subsequent averaged implementation is
[quartet_partition_averaged_potential.py](quartet_partition_averaged_potential.py).
Its static audit passed positive coefficient assembly, face extraction,
normalization, reachable rows, terminal potential, and the pointwise
fallback. Parent-run checks passed 42 physical rows, 29632 birth
transitions, 42 face extractions, 42 rational Chernoff comparisons,
42 FLOAT row comparisons, and two gradient coordinates.

At the preceding six-variable marker choice, a complete 24,897-row scan
took 0.106 seconds. Its raw one-step Chernoff growth upper was about
40.391 bits, above the pointwise comparison of about 35.577 bits. The
fallback therefore retained the profile bound of about +177477 bits.

[quartet_partition_averaged_search.py](quartet_partition_averaged_search.py)
uses 49 positive monomials and approximate envelope derivatives to choose
outer markers. Boundary faces use their reduced tilt coordinates.
Gradients and covariance Schur complements require optimized inner tilts;
with finite inner iterations they guide candidate selection only. Every
saved global bound comes from a complete reachable-row scan.

Static review passed this separation. Parent checks passed 49 monomials,
30 exact collapsed coefficients, six face/row comparisons, and six
derivative coordinates. A bounded search used 247 outer evaluations and
19 complete scans in 1.9 seconds. It returned the preceding +177477 bound;
its current candidate was about +184782. Competing large rows occurred
near the small-population boundary and near (u,a)=(508,503). The search
established no optimum or limitation of the exponential family.

The subsequent sum/minimum control is
[quartet_partition_potential_mixture.py](quartet_partition_potential_mixture.py).
Parent checks passed 17 sum comparisons, 17 minimum comparisons,
34 endpoint comparisons, 34 FLOAT comparisons, four normalization checks,
and eight temporal checks. The bounded control considered 56 candidates
from 19 components, with complete 24,897-row comparisons. Its best sum
bound was +170880.799 bits, compared with the old +177477.464. This is a
modest improvement within the uniform-step family. It remains weaker
than the earlier all-death integrated-birth bound of approximately
+148464 bits on this profile; that whole-profile bound also applies to
the low-death subset.

A joint-quota control retaining C_A used a common potential and 21 quota
marker probes. Its best exponent was +177497.738 at a_H approximately
7.0277778, so it gave no improvement. The low-death external correction
in this control includes +d1 log2(k_minus): it is
-m log2(b)-(m-d1)log2(k_minus), not merely -m log2(beta).
The control is recorded in
[quartet_quota_potential.py](quartet_quota_potential.py), with separate
[checks](quartet_quota_potential_tests.py). These are parent-reported FLOAT
controls, not exact probability certificates.

The return-floor test also failed to establish a limitation. The parent
checked the rational candidate z=21/20 and V=2/25 with directed logarithms;
the reported upper value for the test was -52568 bits. Its rational
full-floor dominance check passed. Thus these growth floors do not prove
that a uniform-step potential cannot close the class. The negative test
value is not itself a successful potential or a failure-probability bound.

This checkpoint makes no optimizer, scalar-family impossibility, or bad
codeword claim. Further tuning of the tested controls is parked. The
next review should retain regional evolution across several updates or
use a comparison that avoids a single uniform one-step maximum. Retain
the field factor in (7) when assembling high-death slices. The low-death
structural regime requires a successful geometry bound. The full
fixed-construction distance certificate remains open.
