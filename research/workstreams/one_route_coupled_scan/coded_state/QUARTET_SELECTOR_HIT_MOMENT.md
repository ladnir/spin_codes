# Birth-supported live selectors

The live-selector quotient counts selected quartet unions without enumerating
actual state masks.  A selected quartet disjoint from the preceding live union
must contain a forced birth.  This condition improves the quotient geometry
bound while retaining its 129 population states.  It does not change the
construction.  The full distance certificate remains open.

The objects and field-rank grouping are those of
[the live-selector quotient](QUARTET_LIVE_SELECTOR_QUOTIENT.md).  The production
parameters are `n=512`, `T=4096`, and `N=nT=2097152`.  Write
`U_1=empty`, and use the formal endpoint `U_(T+1)=[512]`.  A birth set satisfies
`B_t subset U_(t+1) \ U_t`; its total size is `f`.  The selected volume is
`G=sum |U_t|=F-f`, and the actual state volume is `L`.  Let
`A0=sum |U_t intersection U_(t+1)|`.  Necessarily `A0 >= L-f`.

## A separate hit marker

Let `J0` count selected quartets in `U_(t+1)` disjoint from `U_t`, over the
`T-1` random selector transitions.  Each such quartet needs a distinct point
of `B_t`, so `J0 <= f`.  For `0<u<=1`, inserting `u^J0` and multiplying the
result by `u^-f` gives an upper bound.

For width `n`, `m=n/4`, old population `4i`, and new population `4j`, the exact
positive kernel with overlap marker `beta` and volume marker `z` is

```
K[i,j] = C(m,j) z^(4j) / C(n,4i)
         * [X^(4i)] (u+4 beta X+6 beta^2 X^2+4 beta^3 X^3+beta^4 X^4)^j
                      (1+X)^(n-4j).
```

The old-population denominator averages over the positions of the old set
relative to the new partition.  The new-selector factor `C(m,j)` counts which
quartets are selected.  The formal endpoint has weight `beta^(4i)`, with no
hit or volume charge.  The implementation uses its upper bound `beta^512`.

The regional row capacity gives a separate necessary inequality

```
C_H-f <= sum_regions min(h,G_region),   C_H=F-R+128h.
```

For `c>=1`, an upper operator for a region with `s` random selector layers is
the entrywise minimum of `K(z*c)^s` and `c^h K(z)^s`.  Both matrices use the
same population basis.  The first region has 15 layers, and each of the
remaining 255 has 16.  Their positive composition, multiplied by
`c^(-(C_H-f))`, bounds the required regional-capacity marker.  The existing
fixed-profile route factor already includes the union over `H`.

[quartet_selector_hit_moment.py](quartet_selector_hit_moment.py) implements
these positive convolutions and regional minima as a floating diagnostic.
The parent-run checks passed: 13 physical kernel comparisons, 26,912 birth
feasibility checks, 19 floating moment/neutral-marker comparisons, and two
exact-plus-floating regional-cap comparisons, in 0.044 seconds.  Independent
static review passed the kernel normalization, 15/16 layer accounting,
regional minima, and free endpoint.  The source is frozen.  These checks do
not certify production floating-point rounding.

## Exact dyadic replay

[quartet_selector_hit_certificate.py](quartet_selector_hit_certificate.py)
replays the low-birth control using integer arithmetic and directed
logarithms. It uses

```
beta=35/16, z=459/1024, c=33/32, u=1/256,
v_i=16^i.
```

The fixed positive similarity basis balances the matrix entries. Exact
quartic coefficients are rounded upward to a dyadic grid; every matrix
product and scalar multiplication is also rounded upward. Both regional
operators stay in this basis, and their exponents are aligned exactly
before taking the entrywise minimum. The free endpoint is bounded by
`max_i beta^(4i)/16^i=(beta/2)^512`. The capped route factor includes
the H-set union once.

The parent-run exact checks passed 13 physical kernel entries, 26 dyadic
power entries, 52 regional comparisons, and four global moments against
independent Fraction references. Two independent static source audits
passed the coefficient normalization, similarity direction, outward
rounding, regional minima, 15/16 composition, and final logarithm signs.

The two replayed actual profiles were

```
(h,R,Ls,f,F)=(2048,819200,819200,4,950272),
(h,R,Ls,f,F)=(2560,819200,829030,4,983040).
```

| Global mantissa bits | First profile exponent | Second profile exponent |
| ---: | ---: | ---: |
| 60 | +23394 | +11056 |
| 128 | +6014 | -6325 |
| 256 | -311 | -17026 |
| 512 | -4453 | -18217 |

Every row is a valid upper bound. At low precision, rounding very small
transition weights upward added too much mass; the positive results do
not indicate a mathematical defect or a bad codeword. The 256-bit first
profile already met the required threshold, but the 512-bit replay gives
a larger margin. The production default is now 512 bits. Its two-profile
serial replay took 10.427 seconds.

Both 512-bit exponents are below the existing `-172` threshold, so both
whole-profile classes qualify for the `2^-76` allocation after the
conservative 96-bit global profile union. Each class includes all death
counts and free terminal states. These selected classes do not complete
the full-code certificate. Load the exact and FLOAT hit helpers in
separate namespaces, since some public helper names overlap.

## Integrating the birth subsets

The separate hit marker still pays for birth masks independently.  Instead,
introduce a positive birth marker `w` and sum the admissible birth subsets
inside each transition.  If a selected new quartet contains `c` old
coordinates, its `4-c` new-only coordinates contribute `(1+w)^(4-c)`.
When `c=0`, exclude the empty birth subset.  The resulting positive
bivariate polynomial is

```
P(X,w) = (1+w+beta X)^4 - 1.
```

The subtraction removes only the constant coefficient.  Thus the integrated
kernel is obtained by replacing the selected-quartet polynomial above with
`P(X,w)`.  The endpoint allows arbitrary `B_T` outside `U_T` and has weight
`beta^(4i) (1+w)^(512-4i)`.  It has no selected-quartet condition.

For fixed `w>0`, set

```
b = beta/(1+w),  zprime = z(1+w),  u = 1-(1+w)^(-4).
```

Each integrated transition equals the raw hit transition at
`(b,zprime,u)`.  The endpoint contributes the additional global factor
`(1+w)^512`.  Consequently the exact integrated moment is
`(1+w)^512` times the raw hit moment with exact endpoint `b^(4i)`.
For `b>=1`, the existing raw hit moment with endpoint upper `b^512` is a
valid upper bound.  The production helper currently requires `b>1`.

At fixed profile, define the capped route factor

```
Q = min(1, C(G,C_H-f) C(N-F,256h-C_H) / C(8192,h)^255).
```

It applies to each fixed selector and birth configuration.  In the
geometry-only branch, use `beta^(A0-(L-f))>=1` and the positive coefficient
bound `[w^f] P <= w^-f P(w)`.  The resulting logarithmic bound is

```
log2 M_hit(b,zprime,u) + log2 Q
 - G log2 zprime - (L-f) log2 b
 + (F-L+512) log2(1+w) - f log2 w.
```

If the raw moment includes the regional-capacity marker, add
`-(C_H-f) log2 c`.  There is no separate `C(N-G,f)` birth factor and no
separate `u^-f` correction.  In particular, use `selector_hit_moment` or
`selector_hit_regional_moment`, not `selector_hit_bound`.

[quartet_selector_birth_checks.py](quartet_selector_birth_checks.py) is a
separate exact checker.  It enumerates width-four and width-eight physical
partitions, selected quartets and birth subsets; compares every birth-degree
coefficient; and checks the reparameterization, free endpoint and positive
coefficient bound. Independent source review passed. The parent-run checks
passed 101 physical birth coefficients, 10,032 valid birth choices,
39 reparameterizations, 15 endpoints, 18 moments, and 234 coefficient
bounds in 0.034 seconds. No production numerical job was run by the
helper's author.

The parent's bounded floating controls with coarse markers remained
positive at all five tested profiles: approximately 36,037 bits at `h=2048`,
148,464 at `h=4096`, 120,236 at `h=3584`, 200,203 at `h=5120`, and 6,429 at
`h=4`.  These are upper-bound failures for those selected markers, not
codeword witnesses or impossibility results for the integrated count.

The next structural step is to retain output constraints from closed runs
or small actual death counts in this quotient.  The joint sparse-row
growth-score interface remains a separate possible extension of
[the one-active-row theorem](QUARTET_ONE_ACTIVE_ROW.md); it is not implemented
here.  The full multirow remainder is still unresolved.
