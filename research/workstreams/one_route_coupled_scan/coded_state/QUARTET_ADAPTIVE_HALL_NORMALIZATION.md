# Adaptive reverse-time normalization of the Hall maximum

2026-11-01. A uniform flip penalty ignores predecessor sources already
covered by future quartets. This note uses that known coverage in each
reverse-time decision. It proves a positive majorant for the exact Hall
maximum, gives its finite population kernel, and states a 15-type
coefficient comparison. No numerical or full-code certificate is claimed.

Use the exact short matching and saturated dual witnesses from
[QUARTET_SHORT_MATCHING_HALL_COUNT.md](QUARTET_SHORT_MATCHING_HALL_COUNT.md).
Fix every quartet partition and every routed row label. Let `J_t` be
the non-`H` inputs and assume `2*sum_t|J_t|=nT`. Inputs have capacity
two and sources have capacity one. A saturated witness has `C_1=J_1`
and, for `t>=2`,

```
C_t=U_t intersect J_t,
U_t is a union of current quartets, each meeting J_t.
```

Its score is `2|C|-|N^-(C)|`, and the maximum score is the integer flow
deficit `delta`. Put `V_t=U_t^c`. Both formal future sets
`V_(T+1),V_(T+2)` are full.

## One decision with the future fixed

Process updates from `T` down to two, using any fixed order of quartets
within an update. At update `t`, the future hole sets give

```
A=V_(t+1),       W=V_(t+1) intersect V_(t+2),
B=W^c=U_(t+1) union U_(t+2).
```

For a current quartet `Q`, define

```
j=|J_t intersect Q|,
r=|J_t intersect W intersect Q|,
u=|A intersect Q|,
w=|W intersect Q|,
b=2j-r,
a=max(0,4+u-b).                                     (1)
```

The local current-time score is `-|B intersect Q|` if the quartet is
unselected, and `b-|B intersect Q|` if it is selected. Its immediate
gain is therefore exactly `b`.

Hold every remaining earlier choice fixed. Selecting `Q` can add
predecessors only at times `t-1` and `t-2`. At `t-1`, the future set
`U_(t+1)` already covers all positions outside `A`, so at most `u`
new predecessors are possible. At `t-2`, at most four are possible.
Missing early time layers only reduce these costs. Hence the difference
between the two complete scores lies in

```
[b-u-4,b] subset [-a,b].                             (2)
```

The same interval bounds the difference between the two maxima over
remaining earlier choices. Both branches optimize over the same finite
choice space, and (2) holds for every point of that space. Earlier
selected quartets are required only to meet their fixed `J` labels;
their feasible choices do not depend on the current branch.

If `j=0`, selection is forbidden. If `a=0` and `j>0`, selecting cannot
decrease the optimized score, so the unselected branch can be discarded.
The values in (1) do not depend on the decisions for other quartets at
the current update.

## Weighted maxima and the complete decision tree

Fix `y>1`. For `a,b>0`, define

```
w0=(1-y^-a)/(1-y^(-a-b)),
w1=(1-y^-b)/(1-y^(-a-b)).                            (3)
```

They satisfy

```
w0+y^-a*w1=1,          y^-b*w0+w1=1.                 (4)
```

If two branch maxima `F0,F1` obey `-a<=F1-F0<=b`, then

```
y^max(F0,F1) <= w0*y^F0+w1*y^F1.                    (5)
```

When `a=0,b>0`, use `w0=0,w1=1`. When `j=0`, use `w0=1,w1=0`.
These cases agree with the valid forced choices above.

Apply (5) at the root of the reverse-time decision tree. Apply it again
to each remaining branch maximum, with its own fixed future. Induction
on the number of remaining decisions gives a weighted sum over every
complete saturated witness that upper-bounds `y^delta`. The weights
may depend on prior decisions in this reverse order; the induction does
not require independent or fixed weights across the entire tree.

Define the conditional hole and selection probabilities by

```
g=y^-b*w0,             p=w1.                         (6)
```

Equation (4) gives `g+p=1`. At `j=0`, set `g=1,p=0`; at `a=0`, set
`g=0,p=1`. Within a given update, these probabilities depend only on
the fixed future sets, current labels, and current partition. Its
quartet decisions can therefore be sampled independently conditional
on that information.

## Exact transformation to an adaptive cover moment

For the current quartet, the weighted raw-score branches are

```
unselected: w0*y^(-|B intersect Q|),
selected:   w1*y^(b-|B intersect Q|).
```

Factoring out `y^(b-|B intersect Q|)` leaves the probabilities `g,p`.
Summing that common exponent across the update gives

```
2|J_t|-n+|H_t intersect W_t|,
W_t=V_(t+1) intersect V_(t+2).
```

At update one, the forced choice `C_1=J_1` has exactly the same score
expression. There is no `U_1` choice or cover coin. Exact total balance
cancels `sum_t(2|J_t|-n)`. The decision-tree inequality becomes

```
y^delta <= E_adaptive coins[y^(sum_t|H_t intersect W_t|)].  (7)
```

This statement holds for each fixed balanced routing and partition
realization. The expectation in (7) is only over the conditional coins
in (6). An individual sampled cost need not upper-bound the deficit.
The result is a positive moment majorant, not a sampled feasible flow.

The reverse-time start has `A=W=[n]`, so the last source layer is already
charged by the local costs. There is no additional factor `y^n`.
At the other end, input-one labels contribute their own `H_1 intersect
W_1` cost. Exact regional quota conditioning is still required when
averaging the labels.

## A positive polynomial and the full population interface

Let `x` mark the number of current non-`H` inputs. The old coordinate
labels are `W`, `E=A\W`, and `Z=A^c`. If the current quartet is a hole,
the new hole set contains all four positions and its intersection with
old `A` has size `u`. If the quartet is selected, the new hole set is
empty. Thus the next population state is

```
(|V_t|,|V_t intersect A|).
```

An exact local polynomial uses old variables `X_W,X_E,X_Z` and new
variables `X'_W,X'_E,X'_Z`. For each `w+e+z=4`, put `u=w+e`. Its term is

```
multinomial(4;w,e,z) X_W^w X_E^e X_Z^z
 * sum_(j=0)^4 sum_r binom(w,r)binom(4-w,j-r) x^j y^(w-r)
     * [p(j,r,u)(X'_Z)^4
        +g(j,r,u)(X'_W)^u (X'_E)^z].                 (8)
```

Use zero for invalid binomial coefficients. The probabilities in (8)
come from (1), (3), and (6), including the forced cases. All coefficients
are nonnegative. Raising the sum of (8) to `n/4`, extracting old and
new population coefficients, and dividing by the old multinomial
averages the fresh uniform partition exactly.

The physical state is `(a,w)`, where `a=|A|` is divisible by four and
`0<=w<=a<=n`. At width 512 there are 33,153 safe states. The endpoint
vector is

```
q_x(a,w)=(y+x)^w(1+x)^(n-w).                         (9)
```

It is independent of `a`. At an old `W` position, a first-layer `H`
label costs `y` and a non-`H` label contributes `x`; other positions
contribute `1+x`.

## A 15-type coefficient comparison

For positive `alpha,beta`, use

```
phi(a,w)=alpha^a beta^w.
```

For one old quartet type `0<=w<=u<=4`, its exact action on this potential
is the positive polynomial

```
h_(u,w)(x)
 = sum_(j=0)^4 sum_r binom(w,r)binom(4-w,j-r)
       *x^j*y^(w-r)
       *[p(j,r,u)+g(j,r,u)*alpha^4*beta^u].           (10)
```

Define its coefficient comparison by

```
lambda_j=max_(0<=w<=u<=4) [x^j]h_(u,w)/(alpha^u beta^w),
lambda(x)=sum_(j=0)^4 lambda_j*x^j.                  (11)
```

For any current partition, multiplying the local inequalities gives
`K_x phi<=lambda(x)^(n/4) phi`, coefficientwise in `x`. Thus this is
also valid after averaging the fresh partition. It takes a maximum
over local types; it does not retain the full partition-average gain.

For the endpoint, define

```
E0=max(1,alpha^-1,y/(alpha*beta)),
E1=max(1,alpha^-1,1/(alpha*beta)),
E(x)=E0+E1*x.
```

On the three old coordinate types, `q_x/phi` contributes respectively
`(1+x)`, `(1+x)/alpha`, and `(y+x)/(alpha*beta)`. Therefore
`q_x<=E(x)^n phi` coefficientwise.

For the production parameters `n=512`, regional length `R=16`, quota
`h=4096`, `m=n/4=128`, and `C0=binom(8192,4096)`, put

```
b=[x^4096]lambda(x)^2048,
e=[x^4096]lambda(x)^1920 E(x)^512.                    (12)
```

There are 255 ordinary reverse-time regions and one region containing
fifteen cover updates plus input one. Since the initial state is full,
`phi(n,n)=(alpha*beta)^n`. For any integer threshold `d_star`, (7)--(12)
give the complete fixed-row-class bound

```
Pr[some H of size4096 has integer deficit>=d_star]
 <= y^-d_star (alpha*beta)^512 b^255 e C0^-255.       (13)
```

The factor `C0^-256` conditions the 256 regional quotas, and the single
row-set union contributes `C0`. No separate terminal factor `y^512`
appears. Formula (13) may target either the dimension-obstruction
threshold or the stronger low-death thresholds in
[QUARTET_FLOW_FIELD_BRIDGE.md](QUARTET_FLOW_FIELD_BRIDGE.md).

This comparison preserves regional coefficients before applying the
global bound. It still relaxes the old quartet populations independently
in (11), so a positive result need not diagnose the adaptive majorant
itself. The exact polynomial (8) supplies the stronger reference
interface. Any bounded control must distinguish those two scopes and
make no optimizer or full-code claim from isolated marker values.

## A certified limitation of the coefficientwise scalar comparison

At `y=2`, comparison (13) cannot close the obstruction threshold for
any positive `alpha,beta`. For old type `(u,w)=(4,4)` and `j=r=1`,
the hole probability is `127/255`. Its unchanged full-hole branch
cancels the potential and gives

```
lambda_1 >= 4*2^3*(127/255)=4064/255 > 125/8.
```

For old type `(0,0)` and `j=2,3`, selection is forced, giving
`lambda_2>=6` and `lambda_3>=4`. In the coefficient of `x^4096` in
`lambda(x)^2048`, retain counts `(683,682,683)` of those three terms.
Its multinomial coefficient exceeds
`3^2048/binom(2050,2)>2^3220`. Its remaining weight exceeds
`375^682*(125/2)>2^5802`, since `375^2>2^17`. Therefore

```
b>2^9022.
```

Also `E(x)>=(2+x)/(alpha*beta)` coefficientwise. This cancels the
initial potential in (13), and the remaining endpoint coefficient is
at least one. Since `C0<2^8192`, the scalar expression at
`d_star=104860` is strictly larger than

```
2^[255*(9022-8192)-104860] = 2^106790.                (14)
```

The proof is independent of the potential parameters. It diagnoses the
coefficientwise maximum in (11), which can mix returns from incompatible
population phases without paying for transitions between them. It does
not refute the exact adaptive polynomial (8), another marker, or the
distance target. No scalar search at `y=2` is warranted for (11).

Independent source audit of
[quartet_adaptive_cover_checks.py](quartet_adaptive_cover_checks.py)
passed the future-dependent branch probabilities, physical predecessor
expansion, raw-score and unbalanced-offset identities, forced cases,
15 local types, and first-layer endpoint. Its separate exact-integer
scalar obstruction retains stronger multinomial and endpoint terms than
(14). The parent ran the checks serially in 0.518 seconds. They passed
1,007 canonical maxima, 259,201 input subsets, 3,850 physical branches,
3,907 flip bounds, 2,014 decision-tree majorants, 2,006 balanced cover
identities, 450 local coefficients, 450 endpoint coefficients, and
162 coefficientwise partition checks. Its exact scalar floors are
`9091` bulk bits and `9285` initial-times-endpoint bits, with quota
upper exponent `8186`. Hence the stronger certified bound is

```
255*9091+9285-255*8186-104860 = 135200.
```

Thus the coefficientwise scalar tail expression is at least
`2^135200` at `y=2`, uniformly over positive `alpha,beta`.

## A one-mode population action and the remaining regional comparison

For a fixed positive `x`, the action on one exponential potential can
retain every physical old population without a dense matrix. Form the
homogeneous quartic

```
P_(alpha,beta)(W,E,Z)
 = sum_(0<=w<=u<=4) multinomial(4;w,u-w,4-u)
                       h_(u,w)(x) W^w E^(u-w) Z^(4-u).
```

Then

```
(K_x phi)(a,w)
 = [W^w E^(a-w) Z^(n-a)]P_(alpha,beta)^(n/4)
                     /multinomial(n;w,a-w,n-a).      (15)
```

Repeated positive multiplication by this quartic requires only fifteen
shifts of a triangular two-dimensional coefficient array per power.
At width 512, the largest array has 131,841 entries. Two extended-precision
arrays use about 4 MB, and the 128 powers require roughly 85 million
coefficient additions and multiplications. These are operation and
storage counts, not measured runtime claims. Positive scaling and
underflow bounds remain necessary before treating a numerical action
as a certificate.

Equation (15) provides all 33,153 physical rows for one exponential
potential. It is not an action on an arbitrary population vector, so
it does not by itself compute sixteen regional steps. A finite positive
mixture of exponential potentials would need a proved domination after
each action, in the physical-state basis. To retain quota extraction,
that domination must be coefficientwise in the label variable, or
otherwise be applied only after a valid complete-block comparison.
An entrywise minimum in overlapping coefficient coordinates is not
valid merely because the represented functions are positive.

The bounded helper
[quartet_adaptive_population.py](quartet_adaptive_population.py) now
implements this single-mode action. Read-only source audit passed its
positive probability formulas, fifteen shifts, common normalization,
structural triangle zeros, one old-multinomial division, separate
potential-ratio calculation, and fixed physical row order. It uses
three contiguous full-square buffers, about 13 MB at width 512, rather
than the packed two-triangle storage estimated above. This permits
contiguous vectorized shifts while performing additional zero work.

The parent ran its exact checks serially: 180 local coefficients,
360 positive power coefficients, 126 physical partition actions,
252 floating action/ratio comparisons, and 42 exact unit-potential
identities passed. The physical checks include all 35 width-eight
partitions. The identity

```
(K_x 1)(a,w)=(y+x)^w(1+x)^(n-w)
```

follows from `g+p=1` and supplies an independent all-row check at
production width. This helper is neither a bin engine nor a quota or
global-probability computation.

The parent's single width-512 control at `y=2,x=1,alpha=beta=1`
passed all 33,153 rows in 1.590864 seconds. Its maximum absolute
logarithmic error against that closed form was `3.21965e-15`; the
minimum and maximum log actions were `512` and `811.500800369232`.
This was one proof-operator evaluation, not an encoder benchmark or
a probability certificate. No repeated control or marker search was
needed for this checkpoint.

The next unimplemented task is a population-aware positive comparison
that uses (15) while preserving the regional quota. The scalar family
excluded by (14) and the earlier pair-cover comparison are not proposed
as further tuning targets.

The mathematical interface for disjoint bins is given in
[QUARTET_ADAPTIVE_POPULATION_COMPARISON.md](QUARTET_ADAPTIVE_POPULATION_COMPARISON.md).
Independent audit passed its common lift, physical-row minimum before
source-bin maximum, sixteen-step composition before quota minimum,
and fifteen-step terminal comparison. Implementing and evaluating that
interface remains a separate checkpoint. Its first useful global target
is integer deficit `60518`, which would support the reviewed
`D<=65536` field bridge at eligible count 4096.

## Appendix: mandatory predecessors from the fixed earlier setup

The following deterministic refinement uses earlier input labels and
partitions that are part of the fixed setup. It is not implemented and
is not a drop-in replacement for the population kernel above.

Fix the future holes `A=V_(t+1)`, `W=V_(t+1) intersect V_(t+2)` and a
current quartet `Q` at time `t>=2`. Keep `u,j,r` as in (1), and write
`B=2j-r` for the raw immediate selected-minus-unselected score. Define

```
m1=|Q intersect A intersect H_(t-1)|,
m2=|Q intersect H_(t-2) intersect Gamma_(t-1)(J_(t-1))^c|,
m=m1+m2.                                            (16)
```

Here `Gamma_(t-1)` is closure in the partition of update `t-1`.
Set `m2=0` when `t=2`, since time-zero sources do not exist.
Every feasible earlier selected union satisfies
`U_(t-1) subset Gamma_(t-1)(J_(t-1))`: a selected quartet must meet
the fixed non-`H` labels.

At time `t-1`, positions in `Q intersect A` are outside the fixed
future union `U_(t+1)`. Other current quartets cannot cover them.
The positions counted by `m1` also cannot be covered by own inputs
`C_(t-1) subset J_(t-1)`. Selecting `Q` therefore introduces all
`m1` of these predecessor sources for every earlier choice.

At time `t-2`, the positions counted by `m2` are outside both possible
own-input coverage `J_(t-2)` and the largest feasible union
`Gamma_(t-1)(J_(t-1))`. Other current quartets again cannot cover
positions in `Q`. These `m2` predecessors are also necessarily new.
The two counts concern different source layers and add without
overlap. The former upper bounds on newly covered predecessors remain
valid. Consequently every fixed earlier witness has complete score
difference in

```
[B-c_max, B-m],
c_max=u+4  for t>=3,       c_max=u for t=2.           (17)
```

Both branches optimize over the same earlier witness space: its
selected-quartet eligibility depends only on the fixed input labels.
Thus (17) also bounds the difference of their optimized scores. This
argument does not assume any independence between the earlier setup
and the fixed future holes. At `t=2`, the forced choice `C_1=J_1`
in fact makes the predecessor increment exactly `m1`, so the complete
branch difference is exactly `B-m1`.

If `j=0`, selection is forbidden. Otherwise, if `B-m<=0`, one may
force the hole branch; if `B-c_max>=0`, one may force selection.
In the remaining case put

```
a_loss=c_max-B>0,       b_gain=B-m>0
```

and use weights (3) with `(a,b)=(a_loss,b_gain)`. The weighted maximum
induction is unchanged, but its conversion to a probability kernel is
different. The raw immediate score still differs by `B`, not by
`b_gain`. Factoring out the same raw selected score therefore gives

```
g_sub=y^-B*w0,       p_sub=w1,
g_sub+p_sub
 =1-(1-y^-m)*y^-b_gain*w0 <=1.                      (18)
```

These are substochastic weights. A forced hole has
`g_sub=y^-B,p_sub=0`; a forced selection has `g_sub=0,p_sub=1`.
At `j=0`, `B=0` and the hole weight is one. Ties may use either valid
forced branch. One must retain the mass lost in (18). Equivalently,
normalize the surviving two probabilities and multiply their common
factor by `g_sub+p_sub`; silently renormalizing without that factor
discards the improvement.

The raw common exponent is still the one preceding (7). For a fixed
balanced setup the complete decision tree thus gives a positive
substochastic sum of `y^(sum_t |H_t intersect W_t|)` that majorizes
`y^delta`. It can be represented by conditional coins with an absorbing
cemetery outcome of weight zero. This establishes a deterministic
branch/majorant lemma, not an averaged finite-state operator.

The distinction matters. A future decision at time `t+1` now consults
`H_t` and `Gamma_t(J_t)`, so the future hole sets can reveal information
about the current labels and partition. Conditional on those holes,
the current setup no longer has the independent distribution used to
derive the existing 33,153-state kernel. A valid averaged construction
must retain the relevant setup information or prove another exact
conditioning argument. Reusing that kernel with only changed local
weights would be unjustified. This refinement is reserved for later
work if the current row-budget comparison remains insufficient; no
engine or numerical search is proposed for it here.
