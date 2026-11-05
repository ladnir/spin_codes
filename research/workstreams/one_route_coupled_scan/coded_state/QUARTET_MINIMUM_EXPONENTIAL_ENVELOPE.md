# A minimum-exponential comparison on physical populations

2026-11-01. This comparison retains every safe physical population
while using finitely many exact exponential actions. It avoids source
bins and an overlapping linear lift. Its two relaxations are moving a
minimum outside a positive operator and projecting to a fixed family
of concave logarithmic envelopes. Neither relaxation is assumed small.

Let `K_x` be the positive adaptive population operator from
[QUARTET_ADAPTIVE_HALL_NORMALIZATION.md](QUARTET_ADAPTIVE_HALL_NORMALIZATION.md).
Its physical states are `s=(ell,omega)`, where `ell` is divisible by
four and `0<=omega<=ell<=512`; there are 33,153 safe states. Fix
positive exponential modes

```
phi_l(s)=alpha_l^ell beta_l^omega,
```

using the same modes at every quota evaluation `x`. For positive
coefficients `c_l`, define the physical function

```
V_c(s)=min_l c_l phi_l(s).
```

## One update and the least fixed-slope majorant

The exact exponential actions `K_x phi_l` are supplied by the positive
quartic-power engine. Positivity gives

```
K_x V_c <= F_x(c),
F_x(c)(s)=min_l c_l (K_x phi_l)(s).                   (1)
```

Define a coefficient map

```
P_x(c)_l=max_s F_x(c)(s)/phi_l(s).                   (2)
```

Each term `P_x(c)_l phi_l` majorizes `F_x(c)`, so

```
K_x V_c <= F_x(c) <= V_(P_x(c)).                    (3)
```

Equation (2) is the smallest coefficient for each fixed slope. Indeed
`min_l d_l phi_l>=F` requires every `d_l phi_l>=F`, hence
`d_l>=max_s F/phi_l`. No concavity of `F` is needed. The logarithm
of the resulting envelope is a concave piecewise-affine function on
the physical populations; this shape restriction may be costly.

The map `P_x` is monotone and positively homogeneous. More directly,
positivity of `K_x` and induction in (3) give

```
K_x^r V_c <= V_(P_x^r(c)).                           (4)
```

All updates in (4) use one fixed evaluation `x`.

## Regional quota and endpoints

For the half-row problem put `h=4096`, `C0=binom(8192,4096)`, and
`B=[x^h]K_x^16/C0`. Positive coefficients and (4) imply

```
B V_c <= x^-h V_(P_x^16(c))/C0.
```

For a fixed finite set of positive quota evaluations, define

```
R(c)_l=min_x x^-h P_x^16(c)_l/C0.                   (5)
```

The minima over evaluations and modes commute as physical function
minima. Therefore `B V_c<=V_(R(c))`. This is not an entrywise minimum
in an overlapping linear basis, and does not require disjoint bins.
The starting coefficients and modes must be common to all evaluations.

The first-layer endpoint is
`q_x(s)=(y+x)^omega(1+x)^(512-omega)`. First majorize it using

```
e_x,l=max_s q_x(s)/phi_l(s).
```

Then propagate fifteen updates at the same `x`, and set

```
e_l=min_x x^-h P_x^15(e_x)_l/C0.                    (6)
```

For the exact endpoint `b=[x^h]K_x^15 q_x/C0`, this gives
`b<=V_e`. The physical start is the full-hole state `(512,512)`;
there is no additional factor `y^512`.

For any positive coefficient vector `c`, compute physical ratios

```
rho=max_s V_(R(c))(s)/V_c(s),
mu=max_s V_e(s)/V_c(s).                             (7)
```

Then positivity of the exact regional operator gives the whole bound

```
Pr[exists half-row H with integer deficit>=d_star]
 <= C0*y^-d_star *V_c(full)*rho^255*mu.              (8)
```

The maximum coefficient ratio is a sufficient but possibly weaker
comparison; (7) deliberately uses actual physical functions. A common
rescaling of `c` cancels between `V_c(full)` and `mu`, and leaves
`rho` unchanged. One may therefore replace a witness by a normalized
`R(c)` and evaluate another bounded candidate without assuming
convergence or optimality.

## Implementation and verified scope

[quartet_adaptive_envelope.py](quartet_adaptive_envelope.py) implements
(1)-(8) in logarithmic floating arithmetic. The frozen default uses
nine common modes around `(alpha,beta)=(1,2)`, exponent offsets
`-1/2,0,1/2`, evaluations `x=1,2,4`, `y=2`, and target 60518.
It records projection loss and mode activity, and also reports the
single-mode controls from the same exact population actions. It does
not measure the loss in the first inequality of (1).

Independent mathematical and static source audits passed the common
mode requirement, scale handling, physical minima, least majorant,
sixteen/fifteen-step composition before quota minima, physical ratios,
start and endpoint, and the single row-set union. Tiny Fraction tests
use crossing modes with strict `K(min)<min(K)` and compare against
the exact physical kernel. The parent ran these checks serially:
42 physical-minimum and 42 projection comparisons, 21 quota and 21
endpoint comparisons, two Collatz checks, and 42 floating comparisons
passed, including all 35 width-eight partitions.

An independent full-note review passed (1)-(8), including the common
mode requirement and the order of the quota and physical minima.

The parent then ran one bounded production control serially in
42.91 seconds. Its best logarithmic tail upper bound was
`+218575.870`, supplied by single mode 6; the second envelope witness
gave `+218580.985`. Its regional growth ratio had logarithm
`1056.323`, and all quota choices were `x=2`. Projection excess
reached `82.74` bits per update at physical row `(400,211)`; the
dominant regional-ratio row was `(312,309)`. The bound did not close.

These are noncertifying floating results for the tested mode family.
The projection statistic does not include the loss from moving the
minimum outside `K_x`. The control neither establishes a construction
defect nor excludes every minimum-exponential family. Broad slope
search is parked.

The next authorized control minimally generalizes the local action
and endpoint to `(a,b,c)=(2,1,4)`, preserving the old defaults. Its
specified rational marker is `y=3/2`, evaluations are
`x=1,9/4,81/16`, target is `184116`, and the nine modes are centered
at the endpoint gauge `(alpha,beta)=(15/13,6/5)`. The field target
and one-step scalar obstruction are documented in
[QUARTET_PARTIAL_TEMPORAL_FLOW.md](QUARTET_PARTIAL_TEMPORAL_FLOW.md).
That one bounded tradeoff control subsequently remained positive at
`+258810.31584`; its full diagnostics are recorded in the partial-flow
note. Further capacity and slope tuning in this envelope family is
parked. No numerical job was run by the note's author.
