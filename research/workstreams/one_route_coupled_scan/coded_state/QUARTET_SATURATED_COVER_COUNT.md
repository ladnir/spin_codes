# Count saturated covers with actual quartet hits

2026-10-31. The [run-saturation theorem](QUARTET_RUN_SATURATION.md) reduces
the structural obstruction to supports that occupy whole allowed runs
and fill every open allowed run. This note counts those supports directly.
Every selected quartet must have an actual incoming state. The count is
for geometry; it does not assign field amplitudes or prove a whole-code
distance certificate.

Use `n=512`, `T=4096`, `N=n*T`, and the unchanged outer rows and regional
route. For a selector and birth pair `(U,B)`, the allowed sets are

```
A_t=(U_t intersect U_(t+1)) union B_t,
U_1=empty,                  U_(T+1)=[n].
```

A saturated cover chooses which allowed runs to occupy. Every forced run
and every open run is occupied. In addition, each selected quartet of
`U_(t+1)` must meet the chosen `S_t`, for `t<T`. Thus the resulting support
realizes `U` exactly. Its state volume is the total length of the occupied
runs; its deaths equal the number of occupied closed runs.

Write `G=sum_t|U_t|`, `f=sum_t|B_t|`, and `F=G+f`. For a fixed `(U,B)`,
let `c_f` and `c_u` count forced and unforced closed allowed runs. A death
budget `D0` leaves at most

```
sum_(j=0)^min(c_u,D0-c_f) binom(c_u,j)                  (1)
```

cover choices, with zero choices when `D0<c_f`. This counts choices of
whole runs. It replaces the `2^A0` unrestricted binary support patterns,
where `A0=sum_t|U_t intersect U_(t+1)|`. It does not assume that a feasible
cover is unique. Actual hit and volume restrictions can only reduce (1).

## Four labels and one positive transition

At a transition from `U_t` to `U_(t+1)`, use the old sets
`U=U_t`, `Aprev=A_(t-1)`, and `Sprev=S_(t-1)`. Label coordinates by

* `Z`: outside `U`;
* `O`: in `U` but outside `Aprev`;
* `E`: in `Aprev` but outside `Sprev`;
* `S`: in `Sprev`.

The old `E` coordinates lie in unoccupied allowed runs. An old `S`
coordinate lies in an occupied run. If a fresh quartet is selected, old
`S` must continue occupied and old `E` must continue unoccupied. Each old
`O` coordinate starts an unforced allowed run and chooses its occupancy.
Each old `Z` coordinate is either absent or a forced birth. Require at
least one resulting occupied coordinate in the selected quartet. If the
quartet is unselected, all old occupied coordinates close their runs.

These rules enforce both run saturation and exact next-quartet hitting.
They do not impose an artificial condition at the formal endpoint.

Use positive markers `z,w,v,d` for `G,f,L,D`. For positive potential
parameters `alpha,beta,gamma`, define

```
V(U,Aprev,Sprev)=alpha^|U| beta^|Aprev| gamma^|Sprev|.
```

In a fresh quartet let `o,e,s` be the three inside-label counts and put
`r=o+e+s<=4`. Its action, including the future potential, is

```
h(o,e,s) = d^s
  + z^4 alpha^4 beta^r
      * [(gamma*v)^s (1+gamma*v)^o
           (1+w*beta*gamma*v)^(4-r) - 1_(s=0)].       (2)
```

The first term closes the old occupied runs. In the second term, the
subtraction removes exactly the selection with no occupied coordinate.
Equivalently, expand the selected term over `i` occupied old `O` positions
and `j` births, retaining only `s+i+j>=1`. Its coefficients are
`binom(o,i)binom(4-r,j)`, so the action is a positive sum of monomials.

Normalize by the old local potential:

```
ell(o,e,s)=h(o,e,s)/(alpha^r beta^(e+s) gamma^s).        (3)
```

There are 35 compositions. Let `Lambda_pt` be their maximum. Multiplying
the local comparison over a partition bounds one step by
`Lambda_pt^128`. This baseline needs no partition averaging.

For the all-open-filled normal form, the exact terminal weight is

```
t(U,Aprev,Sprev)=1_(|Aprev\Sprev|=0)
                  * v^|U| (1+w*v)^(n-|U|).           (4)
```

An old empty allowed run cannot become occupied only at its last position.
Hence the indicator in (4). All old `O` positions start occupied open
runs, old `S` positions continue, and births outside `U` remain optional.
There is no final death, new quartet, or added live volume. A safe terminal
comparison is

```
E_end=max(1+w*v, v/alpha, v/(alpha*beta*gamma))^n.       (5)
```

## The fixed-profile geometry bound

Fix `(h,R,f,F)`, with `h` eligible outer rows, parity rank `R`, and
`G=F-f`. Put

```
C_H=F-R+128h,
Q_route=min(1,
  binom(G,C_H-f)*binom(N-F,256h-C_H)/binom(8192,h)^255).
```

This is the existing route bound, including the eligible-row union once.
It is applied only to birth-compatible masks. The saturation note explains
the pruning step for arbitrary dimension obstructions; that step can
change the profile. Saturation itself preserves this profile and regional
footprint volumes.

Let `M_sat(z,w,v,d)` sum `z^G w^f v^L d^D` over saturated covers obeying
the transition rules, and average over the fresh quartet partitions.
Include the terminal rule (4). Each cover is counted once; there is no
additional support union or amplitude count.

For `v>=1`, `0<d<=1`, and positive `z,w`, the fixed-profile event with
`L>=L0` and `D<=D0` has probability at most

```
Q_route*z^(-G)*w^(-f)*v^(-L0)*d^(-D0)*M_sat.           (6)
```

Apply `Q_route` while the profile is fixed, then enlarge the remaining
positive sum by the markers. The pointwise comparison gives

```
M_sat <= Lambda_pt^(128*(T-1))*E_end.                 (7)
```

For the dimension/Singleton obstruction, take `L0=R+1` and
`D0=W-R-1`. This lower-volume class is appropriate for saturation;
an exact-volume class is not preserved. Equations (6)--(7) are geometry
bounds, without a projective field factor. Constraints involving the
accepted field map at update `t` cannot simply be attached to the fresh
partition at update `t+1`.

## Average the partition through a quartic

Introduce four formal variables `x_Z,x_O,x_E,x_S`. Normalize the local
action as in (3) before forming the homogeneous quartic

```
Q(x)=sum_(o,e,s>=0,r<=4)
       multinomial(4;4-r,o,e,s)*ell(o,e,s)
                 * x_Z^(4-r)*x_O^o*x_E^e*x_S^s.       (8)
```

Fix an old population vector `c=(c_Z,c_O,c_E,c_S)`, whose entries sum to
`n`. Fresh partition independence gives the exact averaged potential ratio

```
Row(c) = [x^c]Q(x)^128 / multinomial(n;c).             (9)
```

Every reachable state has `u=c_O+c_E+c_S` divisible by four and
`c_S>=u/4`: the old occupied support hits each old selected quartet.
There are 3,219,969 such population rows. The following comparison avoids
enumerating them.

Set `p=c/n`, with `0^0=1`. Positivity of (8) implies

```
Row(c) <= Q(p)^128 / pi(c),
pi(c) = multinomial(n;c)*product_i p_i^c_i.            (10)
```

The denominator is the probability of the population `c` under its own
multinomial proportions. For `n=512`, its exact uniform reciprocal bound is

```
C_mode = 4^512*(128!)^4/512!.                         (11)
```

To prove (11), write
`log pi(c)=log(n!)-n log n+sum_i g(c_i)`, where
`g(a)=a log a-log(a!)` and `g(0)=0`. Its successive differences are
`g(a+1)-g(a)=a log(1+1/a)` for `a>=1`, with initial difference zero.
They increase with `a`, so `g` is discretely convex. Balancing the four
counts minimizes their sum and hence `pi(c)`. The balanced counts are
128 each and are reachable. Thus (11) is sharp for this population domain.
Its logarithm is approximately 13.5 bits, charged once per update,
not once per quartet.

## Thirty-five positive Bernstein rows on the reachable cone

The reachable proportions satisfy `3p_S>=p_O+p_E`. Their simplex is the
tetrahedron with vertices

```
Z,       S,       P=(3O+S)/4,       Qv=(3E+S)/4.
```

Here `Z,O,E,S` also denote the four coordinate-unit vectors. For a
composition `(k_Z,k_S,k_P,k_Q)` of four, define

```
B_k = sum_(o=0)^k_P sum_(e=0)^k_Q
        binom(k_P,o)*binom(k_Q,e)*3^(o+e)/4^(k_P+k_Q)
        * ell(o,e,k_S+k_P+k_Q-o-e).                    (12)
```

These are the degree-four Bernstein coefficients of (8) on that
tetrahedron. To derive them, assign the four samples to their vertices.
A sample at `P` is `O` with probability three quarters and otherwise `S`;
a sample at `Qv` similarly chooses `E` or `S`. The two binomial counts give
(12). Its weights sum to one and are nonnegative. Each coefficient is
therefore a positive mixture of the original local actions.

Let `Lambda_B=max_k B_k` over these 35 compositions. Bernstein weights
on the tetrahedron are nonnegative and sum to one, so

```
sup_reachable_p Q(p) <= Lambda_B,
Row(c) <= C_mode*Lambda_B^128.                        (13)
```

The global moment bound becomes

```
M_sat <= (C_mode*Lambda_B^128)^(T-1)*E_end.            (14)
```

For the same potential, both (7) and (14) hold. Their one-step growth
bounds may therefore be replaced by their minimum, without a union cost.
The 35 positive mixtures also permit the existing log-sum-exp/epigraph
comparison; the mode cost is an additive constant independent of its
markers.

Rational subdivision of the tetrahedron can tighten the Bernstein bound.
Each child coefficient remains a positive rational mixture of the
original coefficients. Take the maximum over all child coefficients;
there is no sum over events and no probability union factor. This changes
the quartic supremum bound, not the mode argument.

No smaller uniform mode factor follows from positivity alone. For
`Q(x)=24*x_Z*x_O*x_E*x_S`, equality holds in (10) at balanced counts, and
the quartic maximum occurs there. Strictly positive perturbations approach
the same example. A sharper factor for the particular saturated kernel
would require a quantitative property of its coefficient distribution,
such as nondegeneracy; no such additional bound is asserted here.

## Current scope and next comparison

At `(h,R,f,F)=(4096,819200,65536,1081344)`, the full structural thresholds
are `L0=819201` and `D0=124517`. The parent's pointwise baseline reported a
floating exponent near `+247479` bits. This does not close that class or
establish an optimized lower limit for the new count.

The parent's averaged-cone checks passed. Its bounded search for this
full structural profile returned

```
FLOAT primal exponent:             +299226.248959 bits
FLOAT dual lower bound in the box: +299225.406721 bits
Reported gap:                            0.842238 bits
Total mode-factor cost:                55203.738372 bits
```

The endpoint and once-per-update mode accounting were checked. The narrow
reported gap indicates that this tested scalar comparison remains too
loose; incomplete optimization does not explain its large positive value.
The numbers are floating diagnostics, not rational certificates or bounds
over an unrestricted parameter domain. The pointwise comparison remains
the better bound at this profile.

The next proof direction is a global rank or expansion argument retaining
the shape of a violating support. Further local scalar refinement is not
the recommended next step. The other profiles and the separate
field-sensitive regimes also remain to be covered.
