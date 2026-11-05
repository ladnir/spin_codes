# Remaining-quota envelopes for the exact population operator

2026-11-01. This note gives two positive comparisons that retain the
remaining regional label quota jointly with the physical population.
They use the fixed construction and the original `(1,1,2)` adaptive
Hall operator. Neither comparison is a completed tail certificate.
The next discriminating check is a small restricted-state comparison,
before any full-size implementation.

Let `n=512`, `R=16`, `h=4096`, and `C0=binom(nR,h)`. Write the positive
one-update operator from
[QUARTET_EXACT_POPULATION_HORNER.md](QUARTET_EXACT_POPULATION_HORNER.md)
as

```
K_x=sum_(j=0)^n x^j K_j,
B=[x^h]K_x^R/C0.
```

The physical state is `s=(a,w)`, with `a` a multiple of four and
`0<=w<=a<=n`. Every comparison below retains these physical rows.
For a fixed positive vector `phi`, define the raw coefficient vectors

```
f_t(r)=[x^r]K_x^t phi,
I_t=[max(0,h-n(R-t)), min(h,nt)] intersect integers.       (1)
```

Only degrees in `I_t` can contribute to degree `h` after all `R`
updates. The recurrence is

```
f_t(r)=sum_(j=0)^n K_j f_(t-1)(r-j).                     (2)
```

For computing the final coefficient, vectors outside the feasible
interval may be set to zero after every update. This discards no path
that reaches degree `h`. The regional quota remains one condition;
the recurrence imposes no separate quota on an individual update.

## A fixed set of quota slopes

Choose positive, distinct numbers `x_1,...,x_L`. Keep positive physical
vectors `d_(t,l)` such that, for every integer `r in I_t`,

```
f_t(r,s)<=U_t(r,s):=min_l d_(t,l)(s)*x_l^(-r).            (3)
```

Initially `I_0={0}` and every `d_(0,l)=phi`. Given (3), set

```
c_l=K_(x_l) d_(t-1,l),
U(r,s)=min_l c_l(s)*x_l^(-r),   r in I_t.                (4)
```

Each term of (2) is bounded by the same slope before summation:

```
sum_j K_j f_(t-1)(r-j)
 <=x_l^(-r) sum_j x_l^j K_j d_(t-1,l)
 =x_l^(-r) c_l.
```

The zero extension outside `I_(t-1)` also satisfies this inequality.
Taking the pointwise minimum proves (4). Thus each step requires one
ordinary Horner action per quota slope, applied to a physical vector.
It does not require a polynomial-valued population vector.

Canonicalize the intercepts over the new feasible interval:

```
d_(t,l)(s)=max_(r in I_t) x_l^r U(r,s).                  (5)
```

This preserves the represented envelope exactly on `I_t`. Indeed,
`d_(t,l)<=c_l`, while the definition of (5) gives
`U(r,s)<=d_(t,l)(s)x_l^(-r)` for every slope. Taking minima gives both
directions of equality. The canonicalization can reduce intercepts
whose slopes are redundant on the feasible interval. It introduces
no additional error on that interval.

In logarithms, (5) maximizes a concave piecewise-linear function of
one variable. Its continuous maximum occurs at an interval endpoint
or a pairwise line intersection. This gives a conservative maximum
over integer degrees and still preserves the envelope on the whole
continuous interval. For the exact integer maximum, inspect the
integer endpoints and the floor and ceiling of each in-range
intersection. Parallel slopes can first be merged by taking the
smaller intercept.

At the last update, `I_R={h}`. Therefore

```
B phi <= U_R(h)/C0.                                     (6)
```

The losses in this comparison come from passing the positive
operator through a quota-dependent minimum and extending each
exponential bound outside the preceding feasible interval. They do
not come from a source-population maximum. A single envelope also
has a structural limitation: for fixed `s`, `log U_t(r,s)` is concave
in `r`. It cannot retain a low valley between two separated quota
modes. Feasibility truncation may help, but does not by itself prove
that this limitation is harmless.

## Binomial-reference bands

The following stronger interface can retain different quota modes in
separate bands. Partition `I_t` into fixed disjoint integer intervals.
Let

```
b_t(r)=binom(nt,r),
f_t(r,s)<=b_t(r)*v_(t,I)(s),   r in band I.               (7)
```

At `t=0`, use the singleton band `{0}`, `b_0(0)=1`, and `v_(0,{0})=phi`.
For an old band `J`, a target degree `r`, and a positive evaluation
`x`, define

```
J(r)=J intersect [r-n,r],
m_(J,r)(x)=max_(u in J(r)) binom(n(t-1),u)*x^u.           (8)
```

Omit an empty intersection. Positivity gives the sufficient
one-update inequality

```
f_t(r,s)/binom(nt,r)
 <= [sum_J min_(x in X_J) x^(-r)*m_(J,r)(x)
                           *(K_x v_(t-1,J))(s)]
       /binom(nt,r)
 =:W_t(r,s).                                            (9)
```

To prove it, put `u=r-j` in the contribution from one old band. For
every `u in J(r)`, (8) gives
`binom(n(t-1),u)<=m_(J,r)(x)x^(-u)`. Substitute this into (2), extend
the restricted positive sum over `j` to all `0<=j<=n`, and then take
the minimum over the chosen evaluations. Sum the separate band
bounds only afterward.

A valid next band vector is

```
v_(t,I)(s)=max_(r in I) W_t(r,s).                        (10)
```

The minimum over evaluations in (9) precedes the band sum, and the
maximum in (10) follows that whole sum. Other orders require their
own inequality. Each expensive `K_x` action depends only on the old
band and evaluation, not on the individual target degree. The
remaining operations combine physical vectors with scalar factors.
The maximum in (8) is a clipped binomial mode: successive terms have
ratio `x*(n(t-1)-u)/(u+1)`, so it requires no population computation.
For rational `x=A/B` and `M=n(t-1)`, one mode is
`floor((M+1)A/(A+B))`, clipped to `J(r)`. If the fraction is an
integer, the preceding integer is also a mode; clipping the displayed
mode still attains the maximum on an interval.

For the label-only scalar kernel `(1+x)^n`, the unbanded factor in
(9) is exactly

```
max_(u in J(r)) Pr[Bin(M,p)=u] / Pr[Bin(M+n,p)=r],
p=x/(1+x).
```

At balanced degrees and `p=1/2`, these factors telescope between
successive updates. Thus the reference does not inherently charge a
new square-root coefficient loss at every update. Banding and the
physical transition can introduce additional losses.

For the final singleton band `{h}`, (7) already divides by the exact
regional normalizer, giving `B phi<=v_(R,{h})`. The binomial reference
also separates the ordinary label-count profile from population
growth. It does not assert a particular number of bands, runtime,
or successful numerical bound.

One audited compression option permits geometric variation within a
band. If the old bound is

```
f_(t-1)(u,s)<=binom(N,u)*min_l d_(J,l)(s)*z_l^u,
```

then that band's term in (9) can be replaced by

```
min_(l,x) x^(-r)*m_(J,r)(x*z_l)*(K_x d_(J,l))(s).
```

After summing bands and dividing by `binom(nt,r)`, a valid new
intercept is `d'_(I,l)(s)=max_(r in I)W_t(r,s)/z_l^r`. Unlike (5),
this last projection need not preserve `W_t` exactly. Its loss must
be measured before using wider bands. The slopes `z=1` and `z=1/2`
represent the flat and `2^(-r)` binomial-normalized unit-input rows
exactly; this motivates a small compression test but proves no
successful full-population bound.

## Retaining the current-label interval

The extension from a restricted `j` sum to all labels in (9) can be
reduced by a uniform point-mass estimate. This is an entrywise
inequality for the operator evaluated at a fixed positive `x`, not a
coefficientwise inequality between polynomials.

Fix one nonzero physical entry of `K_x`. Condition on the fresh
partition and on which quartets take the hole or selected branch.
These choices fix the destination population independently of the
current label counts. Their contribution is a product of 128 local
polynomials `T_uw(x)` or `S_uw(x)`, with a nonnegative coefficient
independent of the labels. Thus the normalized label-degree law of
each physical entry is a positive mixture of sums of 128 independent
quartet label counts.

For a nonzero local factor `F(x)=sum_(j=0)^4 F_j x^j`, define

```
p_j=F_j x^j/F(x),
alpha_x=min_F sum_(j=0)^3 p_j p_(j+1).                   (11)
```

The minimum ranges over all local factors that may appear. For the
full original operator, these are the thirty `T/S` factors from the
fifteen old quartet types. Every `T` has positive degree-zero and
degree-one coefficients; every `S` has positive degree-one and
degree-two coefficients. Hence `alpha_x>0` for every positive `x`.
A restricted operator may use its smaller valid factor set. If a
different model has `alpha_x=0`, the bound below remains valid but
gives no gain. Also `alpha_x<=1/4`, by separating the even and odd
label probabilities.

For the characteristic function `chi(theta)=sum_j p_j exp(ij theta)`,

```
|chi(theta)|^2
 =1-4 sum_(d=1)^4 (sum_j p_j p_(j+d))*sin^2(d theta/2)
 <=1-4 alpha_x sin^2(theta/2).
```

Products, Fourier inversion, and the triangle inequality therefore
bound every point mass of each product law, and hence of every
mixture, by

```
eta_x=(1/(2*pi))*integral_(-pi)^pi
             (1-4 alpha_x sin^2(theta/2))^64 dtheta

     =sum_(k=0)^32 64!/[k! k! (64-2k)!]
           *alpha_x^(2k)*(1-2alpha_x)^(64-2k).           (12)
```

The last expression is the constant coefficient of
`(1-2alpha_x+alpha_x z+alpha_x z^-1)^64`. It has only nonnegative
terms and is an exact rational number when the local coefficients
and `x` are rational. No real-rootedness, unimodality, or independence
between different physical entries is assumed.

For any set of `d` label degrees, the corresponding truncated
evaluated kernel consequently satisfies

```
sum_(j in J) x^j K_j <= min(1,d*eta_x)*K_x              (13)
```

entrywise. In (9), the current label interval is
`J={r-u:u in J(r)}` and has `d=|J(r)|` integer points. The stronger
band inequality replaces each evaluation's term by

```
x^(-r)*m_(J,r)(x)*min(1,|J(r)|*eta_x)
                    *(K_x v_(t-1,J))(s).               (14)
```

Apply this factor before the minimum over evaluations and before
summing old bands. It retains some of the restriction that was
dropped in (9). For the three-state restriction, every nonzero entry
is already a product of 128 copies of one local factor, so the same
proof applies directly. A floating implementation is still a
diagnostic unless all rounding needed for an upper certificate is
controlled; the exact rational formula alone does not certify later
floating actions.

A stronger use of the same point-mass bound retains the variation of
the binomial weights inside an old band. For a fixed positive vector
`d`, the normalized probabilities

```
p_j(s)=x^j*(K_j d)(s)/(K_x d)(s)
```

sum to one and satisfy `p_j(s)<=eta_x`: they are positive mixtures of
the normalized physical-entry laws. Let `J(r)` contain `d_J` quota
values, sort the weights `b_u=binom(M,u)*(xz)^u` on that interval as
`b_1>=...>=b_(d_J)`, and put `q=min(d_J,floor(1/eta_x))`. The exact
capped-mass optimum is

```
Psi(J,r,xz;eta_x)
 =eta_x*sum_(i=1)^q b_i
   +1_(q<d_J)*(1-q*eta_x)*b_(q+1).
```

It is the maximum of `sum_u p_u b_u` subject to `0<=p_u<=eta_x` and
`sum_u p_u<=1`. For the within-band geometric representation, the
sufficient update becomes

```
f_t(r,s)/binom(nt,r)
 <= [sum_J min_(l,x) x^(-r)*Psi(J,r,x*z_l;eta_x)
                         *(K_x d_(J,l))(s)]/binom(nt,r).
```

As before, form the whole band sum before projecting to the next
band's geometric intercepts. This bound is never worse than replacing
all weights by their maximum and using `min(1,d_J*eta_x)`. The
binomial weights are unimodal, so the largest required weights can
be taken consecutively around the clipped mode, resolving ties in
either direction. At the current seven evaluations, at most fourteen
weights are needed. The refinement permits a small scalar calculation
even for wide quota bands, while preserving arbitrary physical
vectors. Its algebra has passed independent review; a useful
compression and a full-population implementation remain to be shown.

## Endpoint and whole-region composition

For the full adaptive population operator, the branch probabilities
sum to one and hence

```
K_x 1=q_x,
q_x(a,w)=(y+x)^w(1+x)^(n-w).
```

Consequently the exact free endpoint is

```
b=[x^h]K_x^(R-1)q_x/C0=B1.                              (15)
```

It can use a separate run of either comparison with initial vector
`1`. Alternatively start at label layer `t=1` from the explicit
endpoint coefficients

```
q_r(a,w)=sum_j binom(w,j)binom(n-w,r-j)y^(w-j),           (16)
```

and apply the remaining fifteen updates. For (3), initialize
`d_(1,l)(s)>=max_(r in I_1) q_r(s)x_l^r`; the simpler choice
`d_(1,l)=q_(x_l)` is valid but weaker. For (7), initialize each band
by `max_(r in I)q_r(s)/binom(n,r)`. Formula (16) depends only on `w`,
so the endpoint initialization needs no physical transition action.

Given valid physical bounds `Bphi<=F` and `b<=E`, let
`rho=max_s F(s)/phi(s)` and `nu=max_s E(s)/phi(s)`. The original
half-row geometry comparison is

```
C0*y^(-d_star)*phi(n,n)*rho^255*nu.                      (17)
```

There is one row-set union and no additional terminal `y^n` factor.
This is a sufficient fixed-witness inequality; it requires no
iteration or convergence assumption. Bounds for other row counts
and the field-count interfaces remain separate obligations in the
full proof.

## Restricted checks and scope

Replacing `K_x` by a positive three-state restriction `L_x` is useful
for a cheap diagnostic of the quota comparison. An envelope upper
bound for `[x^h]L_x^R` is an upper bound for that restricted operator
only. It is neither a full-operator upper bound nor a lower bound for
the full operator. The full identity `K_x1=q_x` must not be assigned
to `L_x`: discarded destinations remove terms from its row sums.

Both interfaces have passed independent mathematical review. A small
quota-envelope source,
[quartet_three_state_quota_envelope.py](quartet_three_state_quota_envelope.py),
passed independent static review. The parent ran its checks serially:
16 each for final coefficients, continuous bounds, scale invariance,
and damping; and 364 each for pruned coefficients, exact integer
envelopes, and projection identities passed.

The restricted three-state control completed in 0.03 seconds. Its
initial maximum comparison-ratio log was `454.8944`, with minimum
`277.1229`. After two half-updates the maximum was `425.67291` and
the minimum `408.0235`; after twelve, they were `425.67291` and
`425.65567`. These floating values show that this single-envelope
comparison remains too loose even in the restricted diagnostic.
They are not certified lower bounds for the full physical operator.
No full-size implementation or further iteration of this comparison
is justified by the test.

The next small comparison uses the binomial-reference bands. A
constant band vector needs fine initial resolution: for unit input
at the full state, (16) divided by `binom(n,r)` is exactly `2^(n-r)`.
A band of width `Delta` can therefore lose `Delta-1` bits at one
endpoint before any population approximation. Width eight is the
initial diagnostic choice; widths 128 or 256 would conceal a material
part of the regional growth budget. This is a restricted-state test,
not authorization or a feasibility claim for a full-size engine.

The baseline source
[quartet_three_state_quota_bands.py](quartet_three_state_quota_bands.py)
passed independent static review of mode clipping, physical action
orientation, feasible slices, the order of minima and sums,
binomial normalization, and the final singleton. The parent ran its
small checks serially: 990 mode comparisons, 112 coefficient-band
inequalities, 12 each for normalization, floating rows, and scale
invariance, and 52 binomial-log comparisons passed.

The width-eight initial control took 1.067 seconds and returned
maximum/minimum comparison-ratio logs `401.9365` and `118.2494`.
Three explicitly bounded half-updates took 4.275 seconds and ended
at `288.8389796` and `263.0044313`. This baseline remains too loose
in the restricted diagnostic. The interval factor (14) is a new
structural correction to that baseline. Its exact helper and
optional caller integration passed separate independent static review.
The helper is
[quartet_quota_anticoncentration.py](quartet_quota_anticoncentration.py).
Its exact checks passed 90 local laws, six central coefficients, 180
homogeneous products, 90 each for mixed pairs, mixed four-factor
products, and positive mixtures, 36,270 interval bounds, three
production constants, and seven argument guards.

The optional caller delta also passed: 990 mode checks, 196
coefficient-band inequalities, 18 normalization rows, 12 each for
floating rows, scale invariance, and the neutral `eta=1` comparison,
52 binomial logs, 702 actual-quartet coefficient atoms, and six
refined floating rows. The default comparison remains unchanged.

The parent then completed the width-one control with the interval
factor and three half-updates in 41.848 seconds. The seven exact
rational `eta` values were approximately
`0.0728337, 0.0734126, 0.0738484, 0.0741146, 0.0745928, 0.0760282,
0.0827526` at the fixed evaluations. Initially the maximum/minimum
comparison-ratio logs were `202.6749196` and `-31.06075`, with
restricted comparison-tail log `+238.95`. After one half-update the
maximum was `161.245197` and the tail log was `-10408.915`. After
three, the maximum/minimum were `161.0313548587` and `142.278699`,
and the tail log was `-10463.875533`. The final witness logs relative
to `H` were `(-644.923415,-503.265671,0)`.

This is the first passing restricted-operator upper diagnostic for
the remaining-quota comparisons. The transition and interval-bound
algebra is proved, but the displayed control uses floating arithmetic
and only three physical states. It is not a full physical-operator
bound or a code-distance certificate.

The direct full-size lift is not yet a practical implementation
interface. Width-one bands peak at 4097 quota values; paired with all
33,153 physical states, one layer has 135,827,841 long-double entries,
about 2.024 GiB before another layer or a cache. The literal
seven-evaluation recurrence would require 229,488 Horner actions per
witness: 32,784 incoming feasible quota values across the sixteen
steps, times seven evaluations. The next scaling objective is a
verified compression of quota profiles, such as the within-band
geometric representation above, while preserving the interval factor
and physical-state inequalities. No full-size engine is claimed.

The first fixed compression test used width-256 bands, slopes
`z={1/2,1,2}`, the same seven evaluations, capped binomial mass, and
the nearby dyadic witness `(log2(phi_O),log2(phi_P),log2(phi_H))=
(0,142,645)`. Its frozen source,
[quartet_three_state_sloped_bands.py](quartet_three_state_sloped_bands.py),
passed independent static review of the capped-mass normalization,
top-weight search and fractional last weight, the distinction between
tilt `xz` and atom cap `eta_x`, the minimum before band summation,
binomial normalization, and reconstructed-intercept bounds. The
optional dictionary diagnostic uses streaming distances, with no
quadratic all-pairs allocation, and does not replace any action.
The parent ran the exact checks serially: 960 each for capped-mass
vertices and top-weight logarithms, 98 each for coefficient and
projection bounds, and eight each for normalization, floating rows,
and scale invariance passed.

The control completed in 3.501 seconds and failed: maximum/minimum
comparison-ratio logs were `3455.8609` and `3172.9278`, with restricted
tail log `+833012.56`. Its largest first-layer projection loss was
`416.997` bits at degree 255 in state `H`; the sum of the per-stage
maximum projection losses was `4480.566`. These losses are measured
against the computed degreewise upper, not against unknown exact
coefficients.

This rules out that particular coarse basis as a useful comparison.
It does not rule out sloped bands or capped-mass bounds. The large
`H` input coefficient favors the all-hole factor, whose leading
binomial profile varies roughly as `2^(-2r)`, outside the chosen
slopes. In addition, the fixed grid `x>=1` overbounds the zero-degree
endpoint; projecting that loose value across a wide band spreads the
error. These are distinct losses and should be diagnosed before
another compression control.

The full physical operator has an exact boundary that needs no
coefficient engine:

```
(K_0 phi)(a,w)=2^w*phi(n,a).
```

With no non-H labels, every quartet is a forced hole; the destination
is `(n,a)` and the cover cost is `w`. This can remove the initial
zero-degree Cauchy overbound. It does not, by itself, repair the
missing within-band slopes. The next credible compression objective
is to retain the observed first-layer quota profiles and structural
endpoints with a verified reconstruction loss and an explicit full
action budget; no general population-curvature shortcut is available.
See [FULL_POPULATION_SCALING.md](FULL_POPULATION_SCALING.md).

The full proof objective remains to bound the exact quota-conditioned
operator while retaining every physical population state. Neither
interface nor the present diagnostics close a codeword class.
