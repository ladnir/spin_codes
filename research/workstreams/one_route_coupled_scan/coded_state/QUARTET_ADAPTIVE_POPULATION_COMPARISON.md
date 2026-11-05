# A positive bin comparison for the adaptive Hall population operator

2026-11-01. The adaptive Hall normalization has an exact positive
operator on 33,153 safe population states. This note gives a comparison
on disjoint population bins. It requires exponential-potential actions,
not a dense physical transition matrix. It preserves a common lift
through regional quota conditioning. No numerical bound or certificate
is claimed.

The physical operator and its local probabilities are defined in
[QUARTET_ADAPTIVE_HALL_NORMALIZATION.md](QUARTET_ADAPTIVE_HALL_NORMALIZATION.md).
Use `n=512`, `m=n/4`, 16 updates per region, quota `h=4096`,
`C0=binom(8192,4096)`, and 256 regions. The marker `y>1` is fixed
throughout a complete comparison.

## Physical states and an exact exponential action

Immediately before processing a reverse-time update, let

```
A=V_(t+1), W=V_(t+1) intersect V_(t+2),
s=(a,w)=(|A|,|W|).
```

The safe state set is `0<=w<=a<=n`, with `a` divisible by four.
The next state is `(a',w')=(|V_t|,|V_t intersect A|)`. Keeping
unreachable states in this set weakens a comparison but does not
invalidate it. The initial state is `(n,n)`.

Fix a positive label evaluation `x`. For local old counts
`(w,e,z)` with `w+e+z=4`, put `u=w+e` and use the adaptive probabilities
`p(j,r,u),g(j,r,u)` from the normalization note. Define

```
h_(u,w)(x;alpha,beta)
 = sum_(j=0)^4 sum_r binom(w,r)binom(4-w,j-r)
   *x^j*y^(w-r)*[p(j,r,u)+g(j,r,u)*alpha^4*beta^u],

F_(x,alpha,beta)(X,Y,Z)
 = sum_(w+e+z=4) multinomial(4;w,e,z)
   *h_(w+e,w)(x;alpha,beta)*X^w Y^e Z^z.
```

Here `X,Y,Z` mark the old types `W,A\W,A^c`. The selected current
quartet has next potential one; a hole quartet has next potential
`alpha^4 beta^u`. Consequently, for
`phi_(alpha,beta)(a,w)=alpha^a beta^w`, the exact action is

```
(K_x phi_(alpha,beta))(a,w)
 = [X^w Y^(a-w) Z^(n-a)] F_(x,alpha,beta)^m
   /multinomial(n;w,a-w,n-a).                         (1)
```

The denominator averages the fixed old population over the ordered
quartet slots. It must not be replaced by a new-population denominator.
All coefficients are nonnegative. No maximum over the 15 local types
is taken in (1).

Dehomogenizing one variable leaves at most `binom(n+2,2)=131841`
coefficients. Repeated multiplication by the 15-term quartic gives all
safe row actions in roughly 85 million positive scalar updates per
potential at this width, without a 33,153-square matrix. This is an
operation count, not a measured runtime. A fixed row ordering with
arrays for `a,w` and either the absolute action logs or the action-ratio
logs suffices for the comparison below. If a helper returns
`log2(K_x phi/phi)`, recover the absolute log action by adding
`a*log2(alpha)+w*log2(beta)` before changing the lift.

## One fixed disjoint lift

Partition the safe physical states into nonempty, disjoint bins
`B_1,...,B_d`. One possible partition retains every value of `a` and
uses a few intervals for `w`; coarser intervals for both coordinates
reduce the comparison dimension further. Fix positive `alpha0,beta0`
and put

```
phi0(a,w)=alpha0^a beta0^w,
H[s,I]=phi0(s)*1[s in B_I].                           (2)
```

The bins and `phi0` remain identical for every quota evaluation `x`
used in one comparison. Auxiliary tilts used to bound individual
destination bins do not change (2).

For a destination bin `B_J` and a real tilt `xi=(xi_a,xi_w)`, define

```
d_J(xi)=min_((a',w') in B_J) (xi_a*a'+xi_w*w'),
phi_xi(a',w')=phi0(a',w')*2^(xi_a*a'+xi_w*w').
```

Pointwise on all destination states,

```
H[.,J] <= 2^(-d_J(xi))*phi_xi.
```

Apply positivity of `K_x`, then divide by the fixed source lift:

```
(K_x H)[s,J]/phi0(s)
 <= 2^(-d_J(xi))*(K_x phi_xi)(s)/phi0(s).             (3)
```

The action on the right is (1) with
`alpha=alpha0*2^xi_a`, `beta=beta0*2^xi_w`. A finite tilt family gives
the valid entry bound

```
M_x[I,J]
 = max_(s in B_I) min_(xi in Xi_J)
   2^(-d_J(xi))*(K_x phi_xi)(s)/phi0(s).              (4)
```

Thus `K_x H<=H M_x` entrywise. Each tilt supplies a simultaneous
physical-row upper bound; taking its minimum incurs no union factor.
Take that minimum before the source-bin maximum. Different destination
bins may use different finite tilt families. If the destination bin
has no state satisfying
`max(0,a+a'-n)<=w'<=min(a,a')`, its entry for that physical source row
is exactly zero and may be removed before the maximum.

The order in (4) matters for sharpness: reversing the maximum and
minimum is valid but weaker. Taking minima between coefficients of
different non-disjoint potential bases is not justified by this proof.

## Regional quotas come after composition

Let the exact conditional regional operator be

```
T_h=[x^h]K_x^16/C0.
```

Although (1) comes from a polynomial, `M_x` in (4) is only a positive
pointwise comparison at a real evaluation `x`. Do not extract a formal
coefficient from `M_x`. Instead, positivity and repeated application of
`K_x H<=H M_x` give, separately for every `x>0`,

```
T_h H <= x^(-h) K_x^16 H/C0
      <= H [x^(-h) M_x^16/C0].                      (5)
```

For a finite set of quota evaluations, define the entrywise minimum

```
Mbar[I,J]=min_x x^(-h)*(M_x^16)[I,J]/C0.             (6)
```

Then `T_h H<=H Mbar`. To see why the minimum is valid, fix a physical
source `s in B_I`. Its row of `H` has exactly one nonzero entry,
`phi0(s)`, so every right side in (5) bounds the same quantity
`(T_h H)[s,J]/phi0(s)`. Their minimum still bounds it. This argument
would fail for a general overlapping positive mixture basis.

The complete order is therefore:

1. Bound a destination-bin indicator with exponential tilts.
2. Minimize those bounds separately for each physical source row.
3. Maximize over each source bin to form `M_x`.
4. Compose all 16 updates using the same regional quota evaluation.
5. Multiply by `x^-h/C0`, then take the entrywise quota minimum.

Taking the quota minimum at each individual update would impose a
different condition and is not allowed. Changing the bins or the lift
between quota evaluations also invalidates (6); such choices require
separate complete bounds.

## Endpoint and full-length comparison

The exact input-one terminal vector is

```
q_x(a,w)=(y+x)^w(1+x)^(n-w).
```

For each quota evaluation, put

```
e_x[I]=max_(s in B_I) q_x(s)/phi0(s),
bbar[I]=min_x x^(-h)*(M_x^15 e_x)[I]/C0.              (7)
```

The same positivity and disjoint-lift argument proves that the exact
first-region terminal vector is at most `H bbar`. The terminal maximum
is over the actual states of each bin; for rectangular bins it can be
computed from the feasible boundary points of the affine log ratio.
The evaluation `x` is common to its fifteen kernel actions and its
input-one labels before taking the minimum in (7).

Let `I_full` contain `(n,n)`. The adaptive majorant and one row-set
union give the whole half-row-class bound

```
Pr[some H of size4096 has integer deficit>=d_star]
 <= C0*y^(-d_star)*phi0(n,n)
    *e_(I_full) Mbar^255 bbar.                       (8)
```

The matrices and vector in (8) already include all 256 quota
normalizers. There is no additional factor `y^512`: the adaptive local
cost charges the last source layer, and (7) charges input one.
The threshold may be the dimension-obstruction threshold or one of the
stronger low-death thresholds in
[QUARTET_FLOW_FIELD_BRIDGE.md](QUARTET_FLOW_FIELD_BRIDGE.md).

An arbitrary positive bin potential `v` may replace direct powering:
if `Mbar v<=rho v` and `bbar<=c v`, then the full probability bound is
`C0*y^(-d_star)*phi0(n,n)*v[I_full]*rho^255*c`. This permits a
non-exponential shape across bins while preserving the same physical
lift and quota proof.

## Limits and the first bounded implementation

The known `y=2` coefficientwise 15-type scalar comparison fails by an
exact positive-term lower bound. That result does not prove failure of
(1), a one-bin partition-average comparison, or (8). A one-bin control
can still retain exact partition averaging; it merely discards the
population state between updates. It is a conservative diagnostic, not
an established no-go for every scalar potential.

Fine bins reduce source variation but increase matrix-composition cost.
Auxiliary tilts also overcount states outside each destination bin.
Both losses must be reported separately from the adaptive majorant
itself. A small bin partition and a fixed bounded tilt family are the
first appropriate control; this note makes no claim that such a control
will close the tail.

Before a production action, small-width checks should compare (1) with
physical partitions and labels, verify (3)--(5) on every toy row, and
compare regional coefficients and the input-one endpoint with exact
enumeration. Positive arithmetic still needs explicit range and
rounding control for a certificate. Ordinary floating evaluation,
clipping negative outputs, or silently dropping positive underflows
does not provide one. No helper or numerical job is added by this note.

## Verified single-potential action checkpoint

Independent mathematical review passed (1)--(8), including the old
multinomial denominator, destination-indicator bound, physical-row
minimum before the source-bin maximum, transition-feasibility zeros,
common disjoint lift, regional quota minimum, and both endpoints.

The separately implemented
[quartet_adaptive_population.py](quartet_adaptive_population.py) exposes
`adaptive_population_action(y,x,alpha,beta,width)`. It retains all
33,153 safe physical rows in order `a=0,4,...,512`, then `w=0,...,a`,
and returns absolute action, potential, and action-ratio logarithms.
It performs no source-bin maximum or regional composition. Its source
received an independent static audit.

The parent ran the exact small checks serially. They passed 180 local
comparisons, 360 power coefficients, 126 physical partition actions,
252 floating action/ratio comparisons, and 42 unit-potential identities,
including all 35 width-eight quartet partitions.

One production-width action with `alpha=beta=1`, `y=2`, and `x=1`
then completed in 1.590864 seconds. All 33,153 rows agreed with the
closed-form identity

```
K_x 1(a,w)=(y+x)^w(1+x)^(n-w).
```

The maximum absolute logarithmic discrepancy was `3.21965e-15`.
The action logarithms ranged from 512 to `811.500800369232`.
The implementation uses three full square coefficient buffers, roughly
13 MB, rather than a packed triangular layout. This is a validated
single-potential operator evaluation in floating arithmetic. It does
not implement the bin comparison, exact regional quotas, or a global
tail certificate.

The next implementation target is the actual disjoint-bin comparison
with a common lift and the reviewed low-death integer deficit threshold
`60518`. The current helper and checks are frozen; no numerical job is
running at this checkpoint.

## A nonlinear comparison retaining each physical row's total mass

The entrywise source-bin maximum in (4) may choose a different physical
source for each destination column. It also discards the exact total
mass of the row. Both losses can be reduced without another physical
potential action.

Retain the destination-tilt bounds before the source maximum:

```
m_(s,J)=min_xi 2^(-d_J(xi))*(K_x phi_xi)(s)/phi0(s),
t_s=(K_x phi0)(s)/phi0(s).
```

The second quantity is the action at zero auxiliary tilt. It is the
exact physical row mass, not the source-bin maximum. If
`r_(s,J)=(K_x H)[s,J]/phi0(s)`, the disjoint bins give

```
0<=r_(s,J)<=m_(s,J),
sum_J r_(s,J)=t_s.
```

For `t_s>0`, define the individual caps

```
c_(s,J)=min(1,m_(s,J)/t_s).
```

Then the true vector `r_s/t_s` lies in the capped probability simplex.
For a nonnegative bin vector `v`, let

```
F_s(v)=max {sum_J p_J v_J:
           p_J>=0, sum_J p_J=1, p_J<=c_(s,J)},
N_x(v)[I]=max_(s in B_I) t_s*F_s(v).                 (9)
```

The maximum in `F_s` is a fractional-knapsack problem: sort the common
values `v_J` in decreasing order, then fill a unit mass budget using
the caps in that order. The same sorting order serves every physical
row. The caps sum to at least one in exact arithmetic. They must not
be renormalized to sum to one; doing so could remove the true row.

Equation (9) is the sharp upper bound obtainable from these individual
caps and the row mass. It gives

```
K_x H v <= H N_x(v).                                (10)
```

The same physical source is used for its entire weighted destination
sum before the source-bin maximum. The map `N_x` is monotone and
positively homogeneous, though generally nonlinear. Positivity of
`K_x` and induction therefore imply

```
K_x^ell H v <= H N_x^ell(v).                        (11)
```

Here `N_x^ell` means repeated application of the map, not entrywise
exponentiation. The quota evaluation `x` remains fixed in all sixteen
applications. Define

```
T(v)[I]=min_x x^(-h)*N_x^16(v)[I]/C0.                (12)
```

Positive coefficient extraction and the common disjoint lift prove
`T_h H v<=H T(v)`. This is a vector inequality for each fixed `v`;
it does not assert that `T` is a linear regional matrix. If a positive
`v` satisfies `T(v)<=rho*v`, then positivity of the physical regional
operator alone gives

```
T_h^255 H v <= rho^255 H v.
```

For the endpoint, retain `e_x` from (7) and set

```
b_NL[I]=min_x x^(-h)*N_x^15(e_x)[I]/C0,
mu=max_I b_NL[I]/v[I].                              (13)
```

The complete probability bound is consequently

```
C0*y^(-d_star)*phi0(n,n)*v[I_full]*rho^255*mu.        (14)
```

No linearity, independent source choice, extra row-set union, or extra
endpoint factor is used in (11)--(14).

For a fixed `x`, `N_x(v)<=M_x v`. Across regional quota evaluations,
however, the old matrix `Mbar` can choose a different evaluation for
each destination column. Thus (12) need not dominate `Mbar v`.
The componentwise minimum of the two valid bounds,
`min(T(v),Mbar v)`, remains valid under the same lift and prevents a
regression. Likewise `min(b_NL,bbar)` is a valid endpoint bound.

A certified implementation may know only an upper bound `tbar_s` on
the physical row mass. It can instead maximize over
`sum_J p_J<=1`, with caps `min(1,m_(s,J)/tbar_s)`, and multiply the
answer by `tbar_s`. If those caps sum to less than one, use all available
caps; do not manufacture the missing mass. For nonnegative `v` this
is a valid upper comparison. The equality-budget proof above requires
the actual mass `t_s`, and floating approximations alone are not a
certificate of either input bound.

This nonlinear interface received an independent mathematical audit.
It is a proposed refinement; no implementation or numerical result for
it is recorded at this checkpoint.

## First disjoint-bin control

The separately implemented population-bin comparison received an
independent static source audit. The parent ran its checks serially:
378 destination-indicator comparisons; 216 each for physical-row,
source-maximum, and composition inequalities; 108 quota comparisons;
21 endpoint checks; two full comparisons; 92 floating comparisons;
one workspace guard; and all 35 width-eight partitions passed.

The first production control used 46 bins, `y=2`, integer deficit
threshold `60518`, `beta0=2`, bin width 64, auxiliary tilts `+-0.5`,
and quota evaluations `x=1,2,4`. It completed in 42.44 seconds and
returned logarithmic tail upper bound `+275797.785214`. It therefore
did not close the tail.

Every regional comparison entry chose `x=2`. The dominant comparison
bin was `448<=a<=508`, `384<=w<=447`, with reported mass `0.9989`.
At `x=2`, destination-tilt excess was at most `2.8074` bits, whereas
source-bin variation reached `49.7586` bits, with mean `20.134`.
The maximum row-mass inflation was `2.8074` bits. These are diagnostics
of the comparison, not of the exact adaptive moment or code failure.

The parent then repeated the bounded comparison with `beta0=4/3`,
matching the endpoint shape at `y=2,x=2`. Its logarithmic tail upper
bound was `+355171.026705`, worse than the first control. At `x=2`,
source-bin variation fell to maximum `21.11` and mean `5.70`; the
dominant comparison bin moved to `448<=a<=508`, `448<=w<=508`.
Changing the lift also shifts the absolute auxiliary-action tilts, so
this comparison does not isolate a pure gauge effect. Neither result
is an optimization claim. The nonlinear row-budget proposal above
remains separate from this first implementation.

## Row-budget control

The nonlinear helper `quartet_adaptive_budget.py` received independent
source review, including the optional retained-row outputs in the bin
helper. Serial checks passed 135 exact capped-simplex vertex comparisons,
42 physical-row and iterated-composition comparisons, 21 quota and
endpoint comparisons, two Collatz checks, 36 floating comparisons, and
the workspace guard. These validate the implementation on small inputs;
they do not certify production floating arithmetic.

The fixed 46-bin control at `y=2`, threshold 60518, and the same nine
auxiliary modes and three quota evaluations completed in 50.58 seconds.
The original direct-power bound was +275797.785214 bits. Two nonlinear
positive-witness evaluations gave +275791.055784 and +275790.249049.
The best gain was only 7.54 bits; the largest regional gain was 0.0264
bits and the endpoint gain 1.10694 bits. No caps were renormalized and
no cap shortage was reported. Every regional output selected `x=2`.

Thus row-mass preservation does not resolve the dominant loss in this
tested comparison. Further tuning of this refinement is parked. All
numbers in this paragraph are FLOAT upper-comparison diagnostics; no
new geometry-tail class or full-code certificate is established.
