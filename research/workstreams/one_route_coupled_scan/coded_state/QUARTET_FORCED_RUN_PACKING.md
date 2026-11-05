# Charge forced runs before packing the remaining state volume

2026-10-31. A declared birth makes its allowed run nonempty. Every closed
such run therefore consumes a death, even when that run has length one.
The capped-run inequality in [the death cover](QUARTET_GLOBAL_DEATH_COVER.md)
does not fully retain this fact when the cap exceeds one. This note gives
an exact support-packing criterion and a stronger positive-moment filter.
The construction and the 45%/40-bit target remain unchanged. No production
bound is asserted here.

Fix a selector and birth pair `(U,B)` from
[the live-selector quotient](QUARTET_LIVE_SELECTOR_QUOTIENT.md). Its allowed
sets are

```
A_t = (U_t intersect U_(t+1)) union B_t,
U_1 = empty,                 U_(T+1) = [512].
```

An allowed run is a maximal consecutive run of one coordinate in `A`.
It is forced when its first position belongs to `B`. Every birth starts
a distinct allowed run, and no run contains another birth. A run is open
when it reaches update `T`; otherwise it is closed. Write `F_c,F_o,U_c,U_o`
for the forced closed, forced open, unforced closed, and unforced open
runs, respectively. Let `ell_r` be a run's positive integer length, and set

```
f     = |F_c|+|F_o|,
c_f   = |F_c|,
C_A   = |F_c|+|U_c|,
A_tot = sum_r ell_r = A0+f.
```

In the quotient, `f>0`. The combinatorial statements below also hold for
`f=0`, but that case does not justify projective normalization by a forced
nonzero coordinate.

The enlarged support family consists of all `S_t subset A_t` containing
every position of `B`. Its volume and deaths are

```
L(S) = sum_t |S_t|,
D(S) = sum_(t=1)^T |S_(t-1) \ S_t|,       S_0=empty.
```

There is no death after `T`. These supports need not realize every selected
live quartet. All sufficiency statements below concern this enlarged
family only.

## Exact volume interval under a death budget

For a finite multiset `E` of nonnegative integers and an integer `b>=0`,
write `Top_b(E)` for the sum of its `min(b,|E|)` largest elements. In
particular, `Top_0(E)=0`.

Fix an integer budget `d>=0`. If `d<c_f`, no support in the enlarged family
has at most `d` deaths. Otherwise put `b=d-c_f` and define

```
L_min(d) = c_f + sum_(r in F_o) ell_r
                 - Top_b({ell_r-1 : r in F_o}),

L_max(d) = sum_(r in F_c union F_o union U_o) ell_r
                 + Top_b({ell_r : r in U_c}).             (1)
```

Then a support of exactly `L` positions and at most `d` deaths exists if
and only if

```
d>=c_f  and  L_min(d)<=L<=L_max(d).                       (2)
```

To prove necessity, every forced closed run contains a nonzero state and
contributes at least one death. These deaths are distinct across runs.
Only `b` deaths remain. An occupied unforced closed run consumes another
death, so at most `b` such runs can contribute volume. Filling all forced
and open runs gives the upper endpoint in (1).

A forced closed run contributes at least one position. An open forced run
that is not completely filled contributes a death: its first position is
present, and some later allowed position is absent. At most `b` open
forced runs can therefore be shortened. Shortening a run saves at most
`ell_r-1` positions. This proves the lower endpoint in (1).

For sufficiency, use prefixes in closed runs. A nonempty prefix of a closed
run has exactly one death, including when it fills the run. Fill each open
forced run either completely, with no death, or with a nonempty proper
prefix, with one death. An open unforced run can have any volume from zero
to its length by choosing a suffix; this costs no death.

Start with one position in each forced closed run, every open forced run
filled, and the unforced runs empty. Its volume is

```
L_base = c_f + sum_(r in F_o) ell_r,
```

and it has `c_f` deaths. Every integer volume between `L_min(d)` and
`L_base` is obtained by shortening at most `b` open forced runs. Indeed,
the `b` largest savings give integer intervals from zero to each saving,
whose sum is the whole interval from zero to their total. Every volume
between `L_base` and `L_max(d)` is obtained by extending forced closed
prefixes, adding open unforced suffixes, and occupying at most `b`
unforced closed prefixes. The same integer-interval argument applies.
This proves (2), including exact volume and the free endpoint.

In particular, `L_min(d)>=f`. The weaker conditions `L>=f` and
`L<=L_max(d)` do not replace the lower endpoint. For one forced open run
of length five and `d=0`, only volume five is possible, although `f=1`.

## A stronger capped-run inequality

For an integer `k>=1`, define

```
R_k^u = sum_(r in U_c) min(k,ell_r),
E_k   = k*c_f + R_k^u,
H_A   = A_tot-L = A0-(L-f).
```

Every actual support in the enlarged family satisfies

```
H_A + kD >= E_k.                                        (3)
```

At most `D-c_f` unforced closed runs are occupied. Every other unforced
closed run consists entirely of holes and contributes its full length to
`H_A`. Charging at most `k` to each occupied unforced closed run proves
`R_k^u<=H_A+k(D-c_f)`, which is (3).

Let `R_k=sum_(r in F_c union U_c) min(k,ell_r)` be the older statistic.
The improvement is exactly

```
E_k-R_k = sum_(r in F_c) (k-ell_r)_+.                   (4)
```

Thus a short forced run receives its full mandatory-death charge. At
`k=1`, (3) is the old inequality `H_A+D>=C_A`. At `k=2`, if `C_(u,1)`
counts unforced closed runs of length one, then

```
E_2 = 2*C_A-C_(u,1),
H_A+2D >= 2*C_A-C_(u,1).                               (5)
```

Forced closed runs of length one receive weight two in (5). The older
`R_2` assigns them weight one. For example, when `F_c` consists of `f`
length-one runs, the gain in (4) is `(k-1)f`.

The family of upper-capacity tests is exact. For `b>=0`,

```
Top_b({ell_r : r in U_c})
  = min_(integer k>=0) [k*b + sum_(r in U_c)(ell_r-k)_+].
```

Consequently `L<=L_max(d)` is equivalent to `L<=A_tot` and all inequalities
`H_A+kd>=E_k` for `k>=1`. The separate conditions `d>=c_f` and
`L>=L_min(d)` remain necessary. The capped tests alone do not assert exact
live-set realization or field compatibility.

## Positive marking without a support union

Let `m=L-f`. On the class `D<=d1`, (3) gives

```
A0-E_k >= m-k*d1.
```

For any `kappa>=1`, the indicator of this necessary condition is at most

```
kappa^(A0-E_k-m+k*d1).                                 (6)
```

The marked geometry weight factors over whole allowed runs. Charge each
closed forced run `kappa^(-k)`, each closed unforced run of length `ell`
the factor `kappa^(-min(k,ell))`, and each open run the factor one. The
product is `kappa^(-E_k)`. All factors are positive. At `k=2`, closed
unforced singletons receive `kappa^-1`; every other closed run receives
`kappa^-2`. No terminal run is closed artificially.

Define `Omega_k(beta,kappa,w,z)` by summing over the same admissible
selector and birth histories as the death-cover moment `Xi`, and averaging
over the random quartet partitions. Give each history weight

```
z^G w^f beta^A0 kappa^(-E_k).                          (7)
```

There is one summand per `(U,B)`, not per support. The usual selected-quartet
hit requirement remains in this geometry sum. Its enlargement from exact
live histories is the same enlargement used in `Xi`.

For comparison with equation (7) of the death-cover note, retain its
fixed-profile route factor `Q_route`, field term `B_x`, death factor `t_x`,
and interval `[d0,d1]`. Choose `b>=1`, `y,w,z>0`, `kappa>=1`, and
`0<theta<=1`, and put

```
beta = b*(1+y)^theta*kappa.
```

The resulting fixed-profile bound is

```
P_profile,[d0,d1]
 <= Q_route*z^(-G)*w^(-f)
      * B_x^theta*t_x^(theta*d0)
      * b^(-m)*y^(-theta*m)
      * kappa^(-m+k*d1)
      * Omega_k(beta,kappa,w,z).                       (8)
```

Apply the route factor while the profile is fixed. Sum the support
amplitudes before applying `theta`, as in the death-cover argument. Then
insert the necessary geometry filter (6) and enlarge the positive moment.
For `theta=0`, omit the field and `y` factors and use `beta=b*kappa`.
At `k=1`, (8) is the previous low-death interface with
`Xi(beta,1/kappa,w,z)`.

Equation (8) uses `d0` in the amplitude penalty and `d1` in the support
capacity filter. These are two necessary restrictions on the same support.
There is no second probability penalty per mandatory death. In particular,
`c_f` must not also be added to `d0` without proving that stronger lower
bound on the actual total `D`.

An exact regional-quota operator can replace the external route factor.
A full-footprint-cap operator instead supplies an additional necessary
geometry filter and retains `Q_route`. Either operator must retain the
same whole-run marking across region boundaries. A run crossing a boundary
is not closed there. Equation (7) specifies the statistic without
prescribing a large state representation or claiming that a scalar
comparison will suffice.

## An exact support polynomial, kept inside the quotient

There is also a small whole-run interface for the exact support count.
Let `s,v,e` mark state volume, deaths, and terminal population. Define

```
T(s,v) = [[1,s], [v,s]],
c_closed = column(1,v),       c_open = column(1,e).
```

Rows and columns of `T` are the preceding and next occupancy bits. A
transition from one to zero contributes one death. For a run of length
`ell`, its support polynomial is

```
unforced: row(1,0) * T(s,v)^ell * c_endpoint,
forced:   s*row(0,1) * T(s,v)^(ell-1) * c_endpoint.      (9)
```

Use the appropriate closed or open endpoint. Let `P_(U,B)(s,v,e)` be the
product of (9) over all allowed runs. Its coefficients are nonnegative,
and `[s^L v^D e^a]P` counts exactly the enlarged supports with those three
statistics. Distinct maximal allowed runs have disjoint state positions,
and their death contributions are disjoint, which proves the product.
The closed endpoint contributes its final death; the open endpoint does
not. A forced first position is counted exactly once.

This polynomial can sharpen the amplitude count without an external
support enumeration. For `0<x<=1`, set

```
c_x = 1+(q-2)x,
v_x = (q-1)*x^2/c_x,
e_x = (q-1)*x/c_x.
```

For any fixed nonzero vector on `B`, the exact weighted amplitude count
over the enlarged supports is

```
(q-1)^(-f) * P_(U,B)(c_x,v_x,e_x).                    (10)
```

Indeed, a support of volume `L`, deaths `D`, and terminal population `a`
has amplitude weight

```
(q-1)^(D+a-f) * x^(2D+a) * c_x^(L-D-a).
```

Every closed state run has two nonzero boundary outputs, every open state
run has one, and each remaining live-to-live transition has weighted
count `1+(q-2)x`. This is uniform in the prescribed nonzero birth values.

For fixed `L` and `D<=d1`, any `s>0` and `0<v<=1` give the upper count

```
(q-1)^(-f) * s^(-L)*v^(-d1)
                 * P_(U,B)(c_x*s,v_x*v,e_x).           (11)
```

The field-rank step still precedes this enlargement. Its residual birth
count is at most `q^(f-R_B)`, its raw live-map probability is `q^(-R_U)`,
and projective normalization divides by `q-1` because `f>0`. Thus (10) or
(11), multiplied by `gamma*q^(f-R)/(q-1)` and `x^(-W)`, gives a fixed-pair
field bound before capping at one. Apply any fractional power only after
this complete amplitude count. The packing filter can instead accompany
the older binomial count, as in (8); the two support counts are alternatives,
not factors to multiply together.

## A scalar restart budget and the next comparison

For an actual support, let `S_t^u` contain its positions lying in unforced
allowed runs, and set `s_t^u=|S_t^u|`. Deaths decompose across allowed runs,
so

```
D >= c_f + sum_(t=1)^(T-1) (s_t^u-s_(t+1)^u)_+.        (12)
```

The sum counts downward crossings of the population levels `1,...,512`.
This is a scalar consequence that could be coupled to regional boundary
populations. It is weaker than retaining coordinate runs. Replacing
`s_t^u` with total population would double-count some deaths already
charged to forced runs, so that replacement is not justified.

The next comparison is the full regional footprint constraint together
with (5), or with the whole-run polynomial. Short forced excursions then
spend their mandatory deaths before other closed runs can supply state
volume. This removes one counting relaxation; it does not yet prove a
negative exponent at the remaining profile or a whole-code certificate.

## Verification

The parent ran the exact checker
[quartet_forced_run_checks.py](quartet_forced_run_checks.py) serially. The
first replay passed 2,376 support families, 19,720 volume intervals,
183,354 capped inequalities, 2,376 support polynomials, and 167 fixed-birth
amplitude counts over GF(4), including the free endpoint. A separate
multi-run replay passed 968 families and 7,808 intervals. It used up to
three runs of lengths one through four, with every forced/open
combination, and included multiple open forced runs.

Two independent derivation audits passed the exact interval criterion,
capped inequalities, positive marking, support polynomial, and amplitude
substitution. The review correction distinguishing a regional quota
operator from a footprint-cap filter is incorporated above. These checks
support the stated counting interfaces; they do not evaluate a production
geometry bound or close the fixed-construction certificate.

## Appendix: a capped-age comparison with few local cases

The statistic `E_k` has a finite label representation. At a selector
transition, write `U` for the old live set and `Aprev` for the preceding
allowed set. Label each coordinate as follows:

* `Z`: outside `U`;
* `O`: in `U` but outside `Aprev`;
* `J_j`, for `1<=j<k`: in an unforced allowed run whose current length is
  exactly `j`;
* `C`: in a forced allowed run, or an unforced allowed run of length at
  least `k`.

There are `k+2` labels. For `k=1`, the `J_j` labels are absent. Selecting a
new quartet retains every old `U` coordinate in the next allowed set.
Its old `O` coordinates become `J_1`, each `J_j` advances one age, and `C`
remains `C`. Interpret `J_k` as `C`. Every birth is forced and enters `C`.
An unselected quartet closes each `J_j` run with weight `kappa^-j` and each
`C` run with weight `kappa^-k`. Thus these transitions retain `E_k` exactly.

Choose positive `alpha,rho,sigma_1,...,sigma_(k-1)`, and set `sigma_k=1`.
For population `u=|U|`, `a=|Aprev|`, and `n_j` coordinates labelled `J_j`,
use the positive potential

```
V = alpha^u rho^a product_(j=1)^(k-1) sigma_j^n_j.
```

In one quartet let `o,c,n_1,...,n_(k-1)` be its inside-label counts, and put
`r=o+c+sum_j n_j` and `a=c+sum_j n_j`. The local action including the future
potential is

```
h = kappa^(-k*c-sum_j j*n_j)
      + z^4 alpha^4 (beta*rho)^r
          * sigma_1^o * product_j sigma_(j+1)^n_j
          * [(1+w*rho)^(4-r)-1_(r=0)].                (A1)
```

For `k=1`, the convention `sigma_1=1` applies. The old local potential is
`alpha^r rho^a product_j sigma_j^n_j`. If their ratio is at most `Lambda`
for every quartet, positivity gives

```
Omega_k <= Lambda^(128*(T-1)) * E_end,

E_end = max(1+w, beta/alpha, beta/(alpha*rho),
                   max_(1<=j<k) beta/(alpha*rho*sigma_j))^512.   (A2)
```

Omit the last maximum when `k=1`. The actual terminal weight is
`beta^u(1+w)^(512-u)`; (A2) bounds its ratio to the potential. The endpoint
has no closure charge and no artificial age transition.

For `r>0`, the logarithm of the normalized expression (A1) is a
log-sum-exp of two affine functions of the inside-label counts. The birth
factor is then `(1+w*rho)^(4-r)`, without a subtraction. The polytope of
nonnegative counts with `1<=r<=4` has vertices consisting of either one or
four coordinates of a single inside label. Convexity therefore reduces the
comparison to these vertices. Handle `r=0` separately, retaining its
subtraction. Only

```
1 + 2*(k+1)
```

local cases are needed. In particular, `k=2` requires seven cases. This
pointwise convexity argument does not extend to a population row after
averaging over partitions.

The age factors also admit a backward feasibility test. Fix the other
markers and a candidate `Lambda>0`. For `r` equal to one or four, put

```
C_r = z^4 alpha^4 (1+w*rho)^(4-r).
```

The monochromatic `J_j` comparisons are equivalent to

```
sigma_j >= max_(r in {1,4}) Lambda^(-1/r)
            * [(kappa^(-j)/(alpha*rho))^r
                       + C_r*(beta*sigma_(j+1)/alpha)^r]^(1/r). (A3)
```

Starting from `sigma_k=1`, choose each preceding factor at its lower
bound. The right side increases with the following factor, so this choice
minimizes `sigma_1` and hence the `O`-case action. Check the `O`, `C`, and
empty cases afterward. Thus bulk feasibility does not require a search
over the `k-1` age factors. To impose a chosen endpoint bound `E_end<=E`,
also require `sigma_j>=beta/(alpha*rho*E^(1/512))` at each backward step,
and check the other terms in (A2).

The substitution `tau_j=kappa^j sigma_j` makes the recurrence independent
of `j`. Its boundary is `tau_k=kappa^k`, and its minimal backward values
are repeated applications of

```
Phi(t) = max_(r in {1,4}) Lambda^(-1/r)
           * [(1/(alpha*rho))^r
                         + C_r*(beta*t/(alpha*kappa))^r]^(1/r). (A4)
```

This is a finite iteration of one positive function, rather than a larger
set of independent tuning parameters.

The added information is not removable by changing the potential. For a
physical width-four example, concatenate two separated full-live phases
of lengths two and three, returning to the empty selector between phases
and at the end. Seed the phases with birth subsets of sizes one and three.
A phase of live length `a` seeded by `b` births has `b` forced allowed runs
of length `a` and `4-b` unforced allowed runs of length `a-1`.

Assign the one-birth seed first to the shorter phase, then swap the two
seed sizes. Both histories have

```
(G,A0,f,C_A) = (20,12,4,8)
```

and empty endpoints. Their unforced singleton counts are three and one,
so their `E_2` values are thirteen and fifteen. Their `k=1` moment weights
are identical, whereas their `k=2` weights differ by `kappa^2`. Potential
ratios telescope on a closed excursion and cannot remove this difference.

This distinction does not guarantee a better optimized scalar bound at a
particular profile. Another local case can dominate both comparisons.
Combining the two valid inequalities retains the earlier information.
The bounded scalar outcomes are recorded below.

## Appendix: actual live hits and compulsory positions

Packing the allowed runs does not ensure that the packed support realizes
every selected quartet. Exact realization additionally requires

```
S_t intersect e != empty
```

for each selected quartet `e` of `U_(t+1)`, for `t<T`. Its candidate
positions are `A_t intersect e`. The selector moment requires this candidate
set to be nonempty, but the packed actual support may omit every candidate.
There is no such requirement at the formal endpoint `U_(T+1)`.

For a small example, take width eight and three updates. Set `B_1={a}`,
`U_2=Q_a`, and `U_3=[8]`, where `Q_a` is one quartet containing `a`.
Set the other birth sets empty. The partition at update three splits
`Q_a` two-and-two between its two quartets. Then

```
A_1={a},       A_2=Q_a,       A_3=[8].
```

The constant support `{a}` has volume three and no deaths, so the packing
criterion passes. Exact realization of `U_3` requires at least two states
at update two. With no deaths, both persist to update three, so every exact
support has volume at least five. No selected quartet has a sole candidate
in this example. Thus singleton forcing alone cannot recover all live hits.

A second example also passes an aggregate hit-volume bound. Use width
eight, six updates, `B_1={1}`, and all other birth sets empty. Set

```
U_1=empty,
U_2={1,2,3,4},
U_3=[8],
U_4=U_5=U_6={1,5,6,7}.
```

At update two use quartets `{1,2,3,4}` and `{5,6,7,8}`. At update three use
`{1,2,3,5}` and `{4,6,7,8}`. At updates four through six use `{1,5,6,7}`
and `{2,3,4,8}`. These choices give

```
A_1={1},       A_2={1,2,3,4},       A_3=...=A_6={1,5,6,7}.
```

Keep coordinate one at all six updates and coordinate five at the final
two updates. This support has `L=8,D=0,a_T=2`, contains every birth, and
satisfies `L=G/4+a_T=24/4+2`. It nevertheless misses the selected quartet
`{4,6,7,8}` at update three. Its sole allowed incoming position is
coordinate four at update two. That position disappears from `A_3`, so
every exact-live support has at least one death.

These examples do not invalidate the amplitude quotient. Its field-rank
step retains exact live hits before enlarging the support family. They
identify a separate loss in a geometry bound for that enlarged family.

Call a position compulsory when it belongs to `B`, or is the sole member
of `A_t intersect e` for a selected quartet `e` at update `t+1<=T`. Call
an allowed run mandatory when it contains a compulsory position. Let
`c_m` count mandatory closed runs, and let `N_c` be the other closed runs.
Every exact-live support satisfies

```
D>=c_m,
H_A+kD >= k*c_m + sum_(r in N_c) min(k,ell_r).          (B1)
```

Each mandatory closed run consumes a distinct death. The proof of (3)
then applies to the remaining closed runs. Several compulsory positions
in one run contribute only one mandatory-run death. An interior compulsory
position need not be the first position of its run, so the exact lower
endpoint in (1) cannot simply treat every mandatory run as birth-forced.

The new compulsory positions have unknown amplitudes. They do not change
`f`, the residual birth equations, or projective normalization. Only the
original `B` positions have the prescribed nonzero amplitudes used in the
field argument.

For `k=2`, a mandatory unforced singleton receives charge two instead of
one. Longer unforced runs already receive charge two. Thus the four-label
automaton needs only one changed transition. When a fresh selected quartet
has exactly one old `U` coordinate, that coordinate has label `O`, and no
birth is chosen, its new allowed run is a compulsory singleton. Put it in
`C` instead of `J_1`. A birth singleton is already in `C`, and a continuing
`J_1` already matures into `C`.

Writing `sigma=sigma_1`, the selected term for `r=o=1` becomes

```
z^4 alpha^4 beta*rho * [1+sigma*((1+w*rho)^3-1)].       (B2)
```

The zero-birth term has lost its factor `sigma`; terms containing births
retain it. All coefficients in the displayed form remain nonnegative.
There is no compulsory-hit promotion at the formal endpoint.

The seven-case convex reduction needs care after this isolated change.
When `sigma>1`, (B2) lowers the `r=1,O` vertex, so the old convex argument
cannot use that lowered value to bound the unchanged interior cases.
One safe reduced list contains eight cases: the modified `r=1,O` case,
the unchanged `r=2,O` case, and the other six old endpoint cases. To see
this, remove the exceptional integer composition and write `a=r-o`.
The remaining integer compositions have convex hull cut out by
`1<=r<=4` and `o+2a>=2`; its new vertex is `r=o=2`. The original
log-sum-exp convexity applies throughout that remaining domain. Keeping
all 35 quartet compositions also remains valid.

In the six-update example, the allowed runs consist of one forced open
run of length six, three unforced closed singletons, and three unforced
open runs of length four. Thus `A_tot=21`, `H_A=13`, and the original
`E_2=3`. Promoting coordinate four's compulsory singleton raises the
statistic to four. The finite-cap inequality still passes at `D=0`,
because `13>=4`; the separate condition `D>=c_m=1` rejects it. Therefore
(B2) is a valid strict refinement, but does not by itself enforce every
mandatory-run constraint.

More generally, assign each selected quartet the set of allowed runs
containing its candidate positions. Every exact-live support must occupy
at least one run in each such set. For a lower bound on deaths, declare all
mandatory and open runs available, and discard sets already meeting them.
The remaining sets contain only nonmandatory closed runs. If nonnegative
weights on these sets sum to at most one over every remaining run, their
total is a lower bound on the additional occupied closed runs. Consequently

```
D >= c_m + sum_(remaining candidate sets e) weight_e.  (B3)
```

Pairwise disjoint sets of candidate runs give an integer version of this
certificate. This is a global covering constraint, not a change to the
amplitude enumeration. It can detect compulsory deaths without singleton
candidate sets. It does not by itself control the minimum volume forced
by open runs, as the three-update example illustrates.

A whole-run Lagrangian can retain both volume and deaths. Assign a
nonnegative price `lambda_e` to each selected quartet, and let
`P_hit=sum_e lambda_e`. Give an allowed position `p` at time `t<T` the
price `c_p=lambda_e` of its selected quartet at time `t+1`. Set `c_p=0`
at time `T`. Each such earlier position belongs to exactly one selected
quartet, because `A_t subset U_(t+1)`. For `eta>=0`, define

```
M_r = max_[binary support pattern s in allowed run r]
          [sum_(p in r) (c_p-1)*s_p - eta*D_r(s)].
```

Retain the required first position when the run starts in `B`. Use the
actual closed or free endpoint in `D_r`; no other live-hit constraints are
imposed inside this maximum. This is a finite optimization within one
run. The actual pattern is among its candidates. Exact live hits give
`sum_p c_p*s_p>=P_hit`, and summing the run maxima therefore proves

```
L+eta*D >= P_hit-sum_r M_r.                            (B4)
```

If the right side exceeds `L+eta*d1`, no exact support has that volume and
at most `d1` deaths. In the three-update example, price both update-three
quartets by two and every other quartet by zero, and take `eta=1`. The
forced run has maximum `-1`; each other run has maximum zero. Equation
(B4) gives `L+D>=5`, recovering the missing volume restriction.

Equations (B3) and (B4) are geometry filters for a fixed `(U,B)`. They do
not count supports outside the quotient or prescribe additional field
amplitudes. No numerical bound for their global geometry sums is asserted.

## Scalar implementation audit and bounded outcomes

[quartet_forced_run_scalar.py](quartet_forced_run_scalar.py) implements the
second cap, a combined first-and-second-cap comparison, and optional
compulsory singleton promotion. Source review passed the modified
zero-birth monomial, its independent physical label transition, the
combined marker corrections, and the free endpoint. The implementation
retains all 35 quartet compositions, so it does not rely on the modified
vertex reduction. It preserves the original birth count and amplitude
normalization.

The parent's compulsory-hit checks passed 70 local comparisons, 266
independent singleton classifications, two modified monomials, and eight
gradient/Hessian coordinates. The checks also passed the explicit
six-update example, its charges three and four, its free endpoint, and
the assertion that thirteen holes still satisfy the refined second-cap
inequality at zero deaths. They therefore do not misidentify that finite
cap as enforcing the separate bound `D>=c_m`.

At the diagnostic profile
`(h,R,L,f,F)=(4096,819200,838861,65536,1081344)`, the bounded scalar searches
reported the following approximate logarithmic upper bounds:

| Comparison | Returned exponent in bits |
| --- | ---: |
| Second cap alone | +182765 |
| First and second caps together | +177477 |
| Combined caps with compulsory singleton promotion | +177477 |

The combined interfaces retain the earlier first-cap control. These
positive results establish no qualifying probability bound for the
profile, no optimizer optimum, and no impossibility of the scalar family.

The parent subsequently ran an epigraph barrier diagnostic for the
combined caps with compulsory hits and death upper bound `104857`. After
259 evaluations and 14 centered stages, it reported

```
FLOAT primal exponent:             +177152.546568 bits
FLOAT dual lower bound in the box: +177140.016189 bits
Reported gap:                           12.53038 bits
```

The narrow reported gap strongly indicates that the scalar relaxation is
too loose in the tested parameter box. Incomplete optimization alone does
not explain the large positive exponent there. The calculation is not a
rational certificate and makes no claim about parameters outside that
box or an impossibility of the fixed construction.

The next proof step is to couple actual live-hit constraints with the
regional count. Further tuning of this scalar comparison is not the
recommended next step.
