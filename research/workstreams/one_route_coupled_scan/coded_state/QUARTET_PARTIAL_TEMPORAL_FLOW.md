# Partial temporal capacity and an all-death field target

2026-11-01. Partial temporal capacity trades a stronger field census
against a weaker flow network. This note gives two explicit half-row
targets, a direct reuse of the population action, and checks that
should precede an implementation. No required deficit tail is proved.
These networks are proof comparisons; the encoder remains fixed.

Use the notation and exact matching theorem from
[QUARTET_GENERAL_FLOW_COVER.md](QUARTET_GENERAL_FLOW_COVER.md), with
`h=4096`, `c=2a`, `0<b<a`, `N=2097152`, and `T=4096`.
The physical deficit is `tau_H=delta_H/c`. A uniform geometry bound
`tau_H<=52429-s` gives a field rank margin `s` for this eligible-row
class.

## Exact all-death census

Set

```
A=1+2^-16,
eta=(a-b)/(2a),       z=q^-eta,
t=sqrt(z)=2^(-8(1-b/a)).
```

The fixed-support field bound is

```
gamma/(q-1) * q^-s * A^(L-D-a_T) * z^D.
```

Summing every support, including its free endpoint, gives the exact
per-coordinate factor

```
Z=e_0 [[1,A],[z/A,A]]^T (1,1/A)^transpose
 =e_0 [[1,t],[t,A]]^T (1,1/t)^transpose.              (1)
```

Thus the entire conditional half-row field class, with no death or
volume restriction, has probability at most
`gamma/(q-1)*q^-s*Z^512`. The maximum-row-sum estimate gives

```
Z^512 <= t^-512*(A+t)^N.                            (2)
```

For `(a,b,c)=(2,1,4)`, `t=1/16` and `A+t<16/15`. Using
`ln2>56/81` and `N*81<202226*56*15` gives
`Z^512<2^204274`. For `(4,1,8)`, `t=1/64` and
`A+t<64/63`. Here

```
N*81=169869312<169869672=48149*56*63,
```

so `Z^512<2^(3072+48149)=2^51221`. Together with
`gamma/(q-1)<2^-31`, these yield the following sufficient margins:

| Network `(a,b,c)` | Census exponent upper bound | Margin `s` | Conditional field bound | Physical cap | Integer cap |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(2,1,4)` | 204274 | 6400 | `<2^-557` | 46029 | 184116 |
| `(4,1,8)` | 51221 | 1602 | `<2^-74` | 50827 | 406616 |

The physical cap is `52429-s`; the integer cap multiplies it by `c`.
Using that integer cap itself as a tail threshold is conservative for
strict violation of the cap. The first row's scalar inequalities were
checked by the parent with exact integers. The second row follows
from the displayed integer comparison and
`51221-31-32*1602=-74`.

There is no row-set union in these conditional field bounds. The
geometry statement must already hold for every half-row set; its
separate probability estimate pays that union once. Other eligible
row counts remain separate classes.

## The same 15-shift population action

Only the local coefficient table changes in the existing homogeneous
quartic action. For an old quartet type `W=w`, `E=u-w`, `Z=4-u`,
sum non-`H` counts `r,s,v` in those three parts. Put

```
j=r+s+v,
B=cj-(a-b)(r+s)-br,
A_loss=max(0,4a+bu-B),
cost=(a-b)(u-r-s)+b(w-r).
```

The coefficient contributing to an exponential population potential
`phi(ell,omega)=alpha^ell beta^omega` is

```
binom(w,r)binom(u-w,s)binom(4-u,v)
 *x^j*y^cost *[p_select+g*alpha^4*beta^u].            (3)
```

The probabilities `g,p_select` use the adaptive loss/gain rule of the
general-flow note, including its forced cases. Raising the resulting
15-term homogeneous quartic to power 128 and dividing by the old
population multinomial gives the exact action on all 33,153 safe
physical states. The positive-shift engine, state ordering, and
regional comparison framework need no mathematical change.

For a comparable marker per unit of source capacity, choose `y^a`
near two. Convenient rational initial choices are `y=3/2` for
`(2,1,4)` and `y=6/5` for `(4,1,8)`. These are bounded calibration
choices, not an optimized family. Write `r0=y^a` and use quota
evaluations `1,r0,r0^2` to represent the empty, central, and full
label phases.

At the central evaluation `x=r0`, the exact first-layer endpoint is
matched by the exponential gauge

```
alpha=(y^(a-b)+r0)/(1+r0),
beta=2r0/(y^(a-b)+r0),
q_x(ell,omega)=(1+x)^512*alpha^ell*beta^omega.        (4)
```

This gives `(alpha,beta)=(15/13,6/5)` for the first network and
`(2376/1921,12/11)` for the second. An endpoint-matching gauge is
only a principled starting point; it is not necessarily a bulk
eigenfunction. Changing it also changes a relative auxiliary-mode
family unless the absolute modes are held fixed.

For an exact regional operator `B`, endpoint `b`, and positive
physical potential `psi`, a bound `B psi<=rho psi`, `b<=mu psi`
gives tail upper bound

```
C0*y^-d_star *psi(full)*rho^255*mu,
C0=binom(8192,4096).                                (5)
```

Consequently the useful bulk range is determined by

```
255 log2(rho)+log2(psi(full)*mu)
 <d_star log2(y)-log2(C0)-desired_failure_bits.       (6)
```

For the displayed markers and targets this permits only about 390
bits of growth per bulk region, before endpoint costs. A large
per-step scalar bound can therefore erase the field-census gain.

## A pure-state obstruction to a one-step scalar maximum

For `(2,1,4)` and `y=3/2`, put `r0=9/4` and `s0=3/2`.
The exact one-step self-loop polynomials at the empty and full-hole
states are respectively

```
S(x)=4*r0^2/(1+r0^2)*x+6x^2+4x^3+x^4,
T(x)=r0^4*[(r0*s0)^4*(1+x/r0^2)^4-(1+x/r0)^4]
             /[(r0*s0)^4-1].                        (7)
```

The corresponding diagonal matrix entries are `S(x)^128` and
`T(x)^128`; no partition averaging changes these pure states.
For every positive potential, its global one-step action ratio is
at least the larger of these diagonal entries.

The ratio `S(x)/(16x^2)` increases for `x>=r0`, while
`T(x)/(16x^2)` decreases for `x<=r0`. For the latter claim, at
`x=r0` the coefficient distribution is a fair binomial law multiplied
by a strictly decreasing hole probability, so its mean degree is
less than two. The derivative of `T/x^2` is therefore negative there;
its second derivative is positive everywhere. The analogous argument
for `S` uses its nondecreasing selection probabilities and `r0>1`.
Direct rational evaluation gives

```
S(r0)/(16r0^2)=33441/24832>4/3,
T(r0)/(16r0^2)=182080305/135000320>4/3.               (8)
```

It follows that any one-step global scalar maximum, followed by a
single common regional Cauchy evaluation, has bulk factor strictly
larger than

```
(2^8192/C0)*(4/3)^2048 >2^(2048/3).                 (9)
```

The last inequality uses `(4/3)^3>2`. This already exceeds the
growth range in (6), independent of the positive potential. It rules
out this organization of the estimate, not the exact population
operator or the fixed-x sixteen-step comparison followed by physical
quota minima. The latter can retain distinct label phases. Pure
constant-state exact-quota paths themselves have conditional weight
at most one; (9) is a loss from the scalar maximum and quota order.

## Explicit Hall witnesses and the choice of temporal fraction

The short matching has `a-b` short copies per source, using only own
and next inputs, and `b` copies that can also use the two-ahead
quartet. Restrict to short copies whose own input is in `H`. For each
next quartet, with `r` current `H` labels and `j` next non-`H` labels,
these give the disjoint Hall deficit

```
((a-b)r-2aj)_+.                                    (10)
```

At `j=1`, this deficit vanishes for every `r<=4` exactly when
`b>=a/2`. Thus `(2,1,4)` is the smallest integral intermediate
network that removes the one-target overload term responsible for
the no-temporal obstruction.

Under independent fair labels, the per-quartet means in (10) are
`1/8` for `(2,1,4)` and `1/2` for `(4,1,8)`. At time `T-1`, all
`a` copies have only own and next inputs, so the last such layer
instead has mean `a*7/32` per quartet. Terminal own-`H` sources
contribute mean `256a`. The resulting disjoint-witness expectations
are 66072 and 263152, respectively. Under the exact regional quotas,
the eight-label pattern argument from
[QUARTET_NO_TEMPORAL_FLOW.md](QUARTET_NO_TEMPORAL_FLOW.md) gives
rigorous lower expectations

```
512 +(505/512)*65560    for (2,1,4),
1024+(505/512)*262128   for (4,1,8).                  (11)
```

Neither exceeds its proposed threshold. These witnesses therefore do
not reproduce the no-temporal contradiction. They do show that the
second network spends a larger fraction of its allowed deficit on a
known local family. A lower witness below the target does not prove
that the remaining matching deficit is small.

Next work should first complete the bounded minimum-exponential
comparison for the existing operator. If a partial-capacity control
is then warranted, `(2,1,4)` is the principled first candidate: it
removes the identified one-target overload while retaining a small
unrestricted field census. Such a control must preserve population
information through the regional quota calculation; (9) prevents a
useful one-step scalar shortcut. No partial-capacity engine or
numerical job was added for this note.

## Exact arithmetic verification

The parent checked the scalar comparisons serially with exact
arithmetic: both transfer row-sum inequalities, field exponents
`-557,-74`, integer caps `184116,406616`, both rational values in
(8), and the Hall means and boundary totals `66072,263152` passed.
These checks verify the stated census targets and lower-witness
calculations. They compute no partial-network deficit tail.

An independent full-note review also passed the census similarities,
endpoint gauges, pure-state polynomials and monotonicity argument,
regional scalar floor, and short-star boundary counts. Its scalar
obstruction has only the scope stated after (9).

## Implemented bounded control

[quartet_partial_temporal.py](quartet_partial_temporal.py) now supplies
the local coefficients and endpoint for the single specified
`(2,1,4)` control. The shared population helper accepts a validated
positive quartic; its 15-shift power loop and old-multinomial division
are unchanged. The envelope helper accepts explicit action and
endpoint callbacks, including both its envelope and single-mode
endpoint consumers. Original defaults remain available.

Independent static audit passed those hooks, the general gain/loss
and forced cases, unit-potential identity, scale handling, requested
marker/threshold, and baseline preservation. The parent ran the exact
checks serially: 900 local coefficients, 126 physical partition
actions, 63 baseline equivalences, 42 unit identities, 997 balanced or
unbalanced decision-tree instances with 3728 branches, and five
callback/validation checks passed. The original population and
envelope checks also passed again after the hooks were added.

The parent then completed the single full control in 42.99 seconds.
Its best logarithmic tail upper bound was `+258810.31584`, from the
second envelope witness. The regional growth logarithm was
`1397.99436`, compared with a useful range of roughly 390 bits.
The dominant physical ratio row was `(172,171)`, and every mode chose
the central quota evaluation `x=9/4`. Maximum per-update projection
excess was `80.1266` bits at `(396,123)`. The central endpoint
projection error was about `2.22e-16`, consistent with the exact
gauge identity (4).

No class closed. These are floating results for the specified
minimum-exponential comparison; they do not establish a deficit-tail
obstruction for the partial network or a construction defect.
Further capacity and slope tuning inside this envelope family is
parked. The next issue is to remove or certify the population
projection loss and retain actual support eligibility, rather than
repeat controls with nearby parameters.
