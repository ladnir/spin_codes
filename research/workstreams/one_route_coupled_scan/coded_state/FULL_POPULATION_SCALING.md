# Conditions for scaling the remaining-quota comparison

2026-11-01. The width-one, anti-concentrated quota-band diagnostic
passed its restricted three-state target: after three half-updates,
the parent reported growth log `161.0313548587`, comparison floor
`142.278699`, and restricted-moment tail log `-10463.875533` in
41.848 seconds. These are noncertifying floating results for three
states. They authorize a compression test, not a full physical proof
engine or a code-distance claim.

## Population curvature cannot be assumed

The exact adaptive action preserves neither log concavity nor log
convexity on physical populations, even for a log-affine input.
Fix `y=2`, `x>0`, and old hole population `a=n`. Use the positive
input `f(a',w')=c^(a'/4)`. For one quartet with `w` old intersection
positions, its action is

```
F_w=S_4w+c*T_4w=A*2^w+B*r^w,
A=256(c-1)(1+x/4)^4/255,
B=(256-c)(1+x)^4/255,
r=(2+x)/(1+x)<2.
```

For `c=2`, both coefficients are positive and `log F_w` is strictly
convex. For `c=512`, `A>0>B`, the action remains positive, and its
logarithm is strictly concave. This follows either by differentiation
or from `F_0 F_2-F_1^2=A*B*(2-r)^2`.

The effect survives uniform partition averaging on the actual width
512. Let `J_w=(K_x f)(512,w)`. Every partition has local intersection
counts summing to `w`. Multiplying the local chord inequalities, then
averaging, gives

```
J_256 <sqrt(J_0*J_512)  for c=2,
J_256 >sqrt(J_0*J_512)  for c=512.
```

Strictness follows because a uniform partition has positive probability
of a quartet with one, two, or three intersection positions. All three
population states used here are reachable nested-set types. Thus a
general curvature argument cannot reduce all physical rows to a few
extremes. This does not exclude a separately proved invariant for a
more restricted potential family.

## The next quota compression and its cost

The tested tiny compression uses width-256 quota bands, three within-band
slopes `z={1/2,1,2}`, the same seven quota evaluations, and one initial
witness. It extends the independently audited band inequalities in
[QUARTET_REMAINING_QUOTA.md](QUARTET_REMAINING_QUOTA.md).
For an old band `J`, represent its coefficient bound by
`binom(M,u)*min_l d_(J,l)(s)*z_l^u`. Replace the maximum binomial
weight by its support function under the proven atom cap `eta_x`:

```
Psi_eta(w)=max{sum_u p_u*w_u : 0<=p_u<=eta, sum_u p_u<=1},
w_u=binom(M,u)*(x*z_l)^u,  u in J intersect[r-512,r].
```

Fill the largest weights up to the cap, with one possible fractional
last weight. Binomial unimodality locates them around the clipped
mode. The update takes the minimum over `(l,x)` before summing bands,
then forms new intercepts by maximizing over the target degrees.
Report the resulting intercept-envelope inflation relative to the
computed degreewise upper bound, not relative to unknown exact
coefficients. No endpoint-only reduction in the quota degree is
justified by the physical-population calculation above.

Sharpness alone is insufficient. At the measured cached full-action
time of about 1.95 seconds, the exact action counts are:

| Quota representation | Full actions per witness | Action time alone |
| --- | ---: | ---: |
| Width 1, one slope, seven evaluations | 229488 | about 5.18 days |
| Width 128, three slopes, seven evaluations | 5712 | about 3.09 hours |
| Width 256, three slopes, seven evaluations | 3024 | about 1.64 hours |

These omit preparation and scalar band work. Width-one full physical
storage is also unsuitable: one layer of 4097 vectors on 33153 states
already exceeds 2 GiB in 16-byte arithmetic. Streaming degrees helps
memory, but not the action count or the degree-by-population scalar work.

## A concrete action-reuse gate

Before any full run, test whether the required physical input vectors
admit a small dictionary with a verified multiplicative distortion.
For each input `v`, select one representative `R` and scalar `c` such
that `v<=cR<=delta*v`, checking all physical rows. Positivity then gives
`K_x v<=c K_x R<=delta K_x v`; one cached representative action serves
every assigned input. Logarithmic distortion adds across update steps.
The required dictionary size and distortion must both be measured;
neither is assumed small. If the same distortion bound is verified for
a pointwise minimum of representatives, applying `K_x` to that minimum
preserves it. Replacing that action by the minimum of separate
representative actions need not preserve the distortion bound.

The continuation gate is therefore: a three-state compression
must retain useful slack, and its projected full action count must fit
an explicit budget. Only then assess representative reuse or a faster
native arithmetic kernel. Changing language alone does not resolve
229488 actions. No implementation or numerical job was performed for
this note; the full-code certificate remains open.

Representatives can be shared globally across all sixteen updates,
because `K_x` is time-homogeneous. Once a representative and its seven
actions are stored, later stages reuse those same arrays. A total
dictionary of eight or twelve representatives therefore needs 56 or
84 actions: about 109 or 164 seconds at the measured action time,
before cache preparation and scalar band work. These are conditional
costs, not evidence that either dictionary size is sufficient.

A verified positive combination can improve the input fit without
new actions. For `r=sum_j b_j R_j`, `b_j>=0`, verify
`d<=c*r<=delta*d` on every physical row. Then
`K_x d<=c*sum_j b_j K_x R_j<=delta*K_x d`. A uniform per-stage
log-distortion `epsilon_t` adds at most `sum_t epsilon_t` to a regional
bound; separate bulk and endpoint errors contribute
`255*E_bulk+E_endpoint` to the complete logarithmic comparison.
Representatives may be added adaptively but remain fixed once cached.
Reject an inadequate dictionary before launching its costly actions.
The eventual full-state verifier must stream row comparisons; a tiny
three-state all-pairs broadcast is not a scalable implementation.

## First-layer failure and the revised gate

The parent executed the sloped-band helper after its independent source
audit and exact toy checks passed (960 support-function checks, 960
top-weight log checks, 98 coefficient and 98 projection checks, and
eight normalization/FLOAT/scaling checks). The 3.501-second restricted
control failed: growth log `3455.8609`, comparison minimum `3172.9278`,
and tail log `+833012.56`. Its first-layer projection alone lost
`416.997` bits at degree 255 in state `H`; the sum of the sixteen
stagewise maximum projection losses was `4480.566` bits. These are
FLOAT comparison diagnostics, not a statement about the code.

The restricted witness `(1,2^142,2^645)` exposes two distinct losses.
Its `H` coefficient profile is within a fixed factor of
`2^(1157-2r)`, whereas the tested slopes have minimum `1/2`. Its `O`
profile has a support activation at degree 128. In fact its adjacent
ratios satisfy `q_126/q_127=965/4>2^7` and
`q_128/q_127=(4/1925)(1+2^114)>2^105`. Any log-concave majorant spanning
these three degrees inflates the middle value by more than 56 bits.
Adding quota slopes cannot remove that particular loss without changing
the band or representation. These facts depend on the restricted
kernel and the stated witness.

The full kernel includes mixed destinations, so those restricted
support gaps are not a sound basis for its compression. Consider
instead the full positive witness `phi(a,w)=2^(a/4+w)` and source states
`O=(0,0)`, `P=(512,0)`, and `H=(512,512)`. Their exact first-layer
normalized coefficient profiles are

```
q_s(r)=[x^r] F_s(x)^128 / binom(512,r),
F_O=(1+x)^4+(1+4x/5),
F_P=A*(1+x)^4+B*(1+x/4)^4,
F_H=16*[A*(1+x/2)^4+B*(1+x/4)^4],
A=224/255, B=7936/255.
```

These are cheap positive one-dimensional coefficient calculations.
For a quartet polynomial `F=sum_j binom(4,j)*a_j*x^j`, double counting
an `r`-subset and one added position shows that `q(r+1)/q(r)` is a
positive weighted average of the local ratios `a_(j+1)/a_j`.
Consequently the full-source ratios lie in `[3/5,1]` for `O`,
`[23/85,85/116]` for `P`, and `[131/510,45/152]` for `H`.
There is no large one-degree support jump. Nevertheless `O` satisfies
the exact strict convexity test
`q(0)*q(2)/q(1)^2=1+7/3066>1`. The parent's
[`full_source_profile_geometry_main`](quartet_quota_profile_geometry.py)
replay checked curvature by exact integer comparisons: `O` has 510
strictly convex points and one flat point; `P` and `H` each have 511
strictly convex points. Exact integer chord tests gave full-interval
loss lower bounds of 46, 88, and 12 bits for `O`, `P`, and `H`;
the corresponding width-256 bounds were `(15,18,1)` on degrees
0 through 255 and `(5,24,3)` on degrees 256 through 511. A general
log-concave quota envelope is therefore still the wrong default.

The next bounded test should stop at this first-layer representation.
Use the exact full-source profiles to construct common adaptive quota
intervals. On each interval, form a secant for each source and raise
its intercept by the maximum residual over every integer degree;
this gives an upper bound without assuming curvature. Split at the
largest reconstruction error until the verified error is at most one
bit. Report the interval count and the implied action-reuse budget
before considering another sixteen-step control. A positive sum of
exponential representatives is also admissible if the same pointwise
distortion test passes. The local two-phase formulas above do not,
by themselves, imply a small exponential mixture after coefficient
extraction. Three source profiles cannot establish a small dictionary
on all 33153 physical states; that remains a separate gate.

This compression preflight is parked pending the parent's full
population control at marker `y=3/2`. No further compression caller or
full-population job is authorized by these first-layer results.
The complete first-layer section and the linked integer geometry
checks passed independent read-only review.

## Why the smaller marker is a distinct test

Changing `y` changes the proof's exponential moment, not the encoder.
For the original `(1,1,2)` flow capacities, a dense quartet admits an
exact two-channel expression. Let `u` be its number of old hole
positions, `w` its number of old intersection positions, with
`0<=w<=u<=4`. Use the positive physical potential
`phi(a,w)=alpha^a*beta^w` and put `c_u=alpha^4*beta^u`. For a label
assignment with `j` nonH positions and `r` of those in the intersection,
the adaptive hole probability is

```
b=2j-r,  a_local=4+u-b,
g=(y^a_local-1)/(y^(4+u)-1)  when a_local>=0.
```

This includes `g=0` at `a_local=0` and `g=1` at `b=0`. If
`u+w>=4`, all label assignments have `a_local>=0`, so no forced-branch
truncation changes the formula. Define the label polynomials

```
Q_uw=(y+x)^w*(1+x)^(4-w),
R_uw=y^(u+w-4)*(y^2+x)^4.
```

Summing `x^j*y^(w-r)*g` over the labels gives the exact hole polynomial
`T_uw=(R_uw-Q_uw)/(y^(4+u)-1)`. The quartet action on the potential is
therefore

```
h_uw=Q_uw+(c_u-1)*T_uw
    =[(y^(4+u)-c_u)*Q_uw+(c_u-1)*R_uw]/(y^(4+u)-1).
```

At `y=3/2`, `alpha=2^(236/512)` and `beta=2^(232/512)`, both channel
weights are positive. Indeed `c_u>1`, and raising the other inequality
to power 512 gives
`c_u^512=2^(944+232u)<2^(256*(4+u))<(3/2)^(512*(4+u))`.
The last strict inequality follows from `2<9/4`. Thus `h_uw` is a
positive convex combination of `Q_uw` and `R_uw` for these dense types.

At half label quota, a factor `(c+x)^4` prefers evaluation `x=c`.
The displayed channels therefore range from evaluations 1 to `y^2`:
1 to 2.25 at the new marker, compared with 1 to 4 at `y=2`.
This reduces one identified source of Cauchy comparison loss. It does
not bound the regional growth below the target, and the formula must
not be applied without the qualification `u+w>=4`.

There is a countervailing change. For fixed positive integers
`a_local,b`, the ratio `(y^a_local-1)/(y^(a_local+b)-1)` decreases as
`y` increases: it is the probability of the initial `a_local` indices
under geometric weights on `a_local+b` indices. Hence lowering `y`
increases the base hole probability. It simultaneously reduces the
reward `y^(H intersect W)` and its bias toward H labels. Under the
reference label weight `Q_uw`, an intersection position is nonH with
probability `x/(x+y)`, which increases as `y` decreases. The full test
measures this balance; neither the comparison moment nor its typical
hole population is the actual maximum flow deficit.
