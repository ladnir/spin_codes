# Retain actual populations through the regional footprint cap

2026-10-31. A regional comparison retains population transitions that a
uniform scalar potential can lose. The bounded control below improves
its uncapped counterpart by `75681.858270` bits, but its resulting bound
is still positive. This is a diagnostic checkpoint, not a distance
certificate or a claim about an optimized family.

The construction and once-sampled setup remain fixed. Use `n=512`,
`T=4096`, `M=8192`, `N=nT`, and `W=943718`. Each region contains 16 updates.
The probability space consists of the independent uniform quartet
partitions and regional routing permutations in the existing setup.
The calculation averages this setup; it does not resample setup for
different support histories.

## The class and its route constraint

For a support history `S=(S_1,...,S_T)`, set `S_0=empty` and define

```
U_t=Gamma_t(S_(t-1)),        B_t=S_t\U_t,
V_t=U_t union B_t,
G=sum_t|U_t|,               f=sum_t|B_t|,
F=G+f,                     L=sum_t|S_t|,
D=sum_t|S_(t-1)\S_t|.
```

There is no death after `T`. Let `H` be the outer rows with covered
degree at least 129, let `h=|H|`, and let `R` be the footprint rank cap
`sum_rows min(128,covered degree)`. For fixed `(h,R,f,F)`, put

```
Cout=R-128h,       C_H=F-Cout,       G=F-f.
```

We count histories with `B subset H`, `L>=L0`, and `D<=D0`. The defaults
`L0=R+1` and `D0=W-R-1` give the sufficient dimension-obstruction class.
The [pruning and saturation theorem](QUARTET_RUN_SATURATION.md) explains
why a global search for that obstruction may require `B subset H`.
Pruning can change the profile; it is not a fixed-profile transformation.

Fix a support history and the quartet partitions before averaging the
routing. The existing forced-birth route bound for its profile is

```
Q=min(1, binom(G,C_H-f) binom(N-F,256h-C_H)
         / binom(M,h)^255).                          (1)
```

This bound already includes the union over `H`. Let `F_j` be its full
footprint volume in region `j`. Each region has exactly `h` slots routed
to `H`, so every compatible routing satisfies

```
C_H <= sum_j min(h,F_j),
sum_j(F_j-h)_+ <= Cout.                              (2)
```

Consequently the route probability is at most `Q` times the indicator
of (2). The cap uses the full footprint, including births outside `U`.
Its threshold is the fixed value `h`; no union over thresholds is needed.

## The actual-support kernel

Take positive markers `z,w,v,d` for `G,f,L,D`. For old/new populations
`a,b`, define the positive quartic and its population kernel by

```
P(X,Y)=(1+w*v*Y)^4
       +z^4 sum_(r=1)^4 binom(4,r) X^r
                    (d+v*Y)^r (1+v*Y)^(4-r),

K_ab=[X^a Y^b] P(X,Y)^128 / binom(n,a).               (3)
```

When a quartet has no old state, every new state position is a forced
birth. Otherwise the quartet is live, costs `z^4`, and has no forced
birth. Each disappearing old state costs `d`; each new state costs `v`.
Thus (3) retains exact live hits, forced births, deaths, and state volume.

For a fixed old support, a fresh partition is uniform. Summing over all
new supports of size `b` gives (3), independently of the old support's
labels. This permutation symmetry justifies the common population basis
`0,...,512`, including for a composed regional operator.

Every transition charges the current update. The unrestricted moment is
`e_0 K^4096 1`, with no additional endpoint update or terminal death.
Unlike the [saturated-cover count](QUARTET_SATURATED_COVER_COUNT.md),
this kernel also counts unsaturated supports. Neither comparison is
claimed to dominate the other.

## Compose the cap inside each region

Fix `eta>=0`. Let `C_eta` be the exact 16-update population operator,
with each regional history multiplied by `2^(-eta*(F_j-h)_+)`.
Since `F_j` is the sum of its live volume and forced births,

```
C_eta <= B_eta
      := min(K(z,w)^16,
             2^(eta*h) K(z*2^-eta,w*2^-eta)^16),      (4)
```

entrywise in the same population basis. All other markers in (4) are
unchanged. For each history, its cap factor is at most both factors
used on the right. Summing proves (4); positivity permits composition.

On the target class, with `v>=1` and `0<d<=1`,

```
v^(L-L0) d^(D-D0) >= 1,
2^(eta*(Cout-sum_j(F_j-h)_+)) >= 1.
```

Apply these inequalities after (1)--(2), then sum histories and average
the quartet partitions. The resulting whole-class upper bound is

```
Q z^(-G) w^(-f) v^(-L0) d^(-D0) 2^(eta*Cout)
  * e_0 B_eta^256 1.                                (5)
```

There is no extra eligible-set union, quota divisor, birth binomial,
amplitude factor, or field estimate in (5). Different choices of `eta`
bound the same event, so their final scalar bounds may be minimized.
Their regional matrices cannot be minimized without also retaining the
corresponding global correction.

The implementation forms `K^16-Z` by positive arithmetic and restores
`Z=e_0 e_0^T` once. Regional `(0,0)` entries also contain nonempty
excursions and remain present in all global powers. The terminal vector
is all ones.

## Verification and paired controls

The frozen implementation is
[quartet_actual_regional_rank.py](quartet_actual_regional_rank.py).
It reuses positive power and validation primitives from
[quartet_actual_mds_moment.py](quartet_actual_mds_moment.py); its new
fixed-buffer convolution assembles the polynomial in (3).

Parent-run checks passed: 25 exact local coefficients, 106 physical
matrix entries at widths 4 and 8, 10 free-endpoint checks, 80 regional
population checks, 281 capacity checks, and 132 floating toy comparisons.
The physical references enumerate masks and quartet partitions. Exact
regional tests compare capped histories with (4), then compare physical
and population compositions. An independent static source audit passed.

The parent ran two controls serially on Peach, with the same profile
and markers:

```
(h,R,f,F)=(4096,819200,65536,1081344),
L0=819201,        D0=124517,
(log2 z,log2 w,log2 v,log2 d)=(-1,0,1,-1).
```

| eta | log2 moment | external correction | raw log2 bound | seconds |
| --- | ---: | ---: | ---: | ---: |
| 0 | 890248.134690 | -46404.341678 | +843843.793012 | 7.702 |
| 1/8 | 777702.276420 | -9540.341678 | +768161.934742 | 12.215 |

Both controls used the 175-second, 2-GiB memory, and 20-GiB free-space
guards. All arithmetic in these production controls is floating, not
directed. The paired comparison isolates the capacity penalty at fixed
markers. It establishes neither marker optimality nor a lower bound on
the true failure probability. Both reported upper bounds are vacuous
after capping at one.

The measured gain confirms that regional capacity removes substantial
mass from this comparison. It does not close the profile, and no further
marker probe is planned for this checkpoint. The next global count should
extend the [synchronous obstruction bound](QUARTET_SYNCHRONOUS_OBSTRUCTION.md),
which treats an entire small family of histories through routing rank.
A useful extension would count deviations or contours of partial states,
retaining that family's small description count and the regional routing
concentration. The full construction-level proof remains open.
