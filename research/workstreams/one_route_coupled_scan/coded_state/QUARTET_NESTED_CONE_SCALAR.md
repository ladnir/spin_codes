# A reachable cone for the nested scalar comparison

The integrated-birth process has a population invariant that is absent
from the unrestricted local maximum.  Every selected quartet contains
at least one allowed coordinate: an intersection with the old live set,
or a mandatory birth when that intersection is empty.  Therefore each
intermediate nested state `(U,Aprev)` satisfies

```
Aprev subset U,       4|Aprev| >= |U|.                 (1)
```

The initial state is empty and satisfies (1).  Every transition preserves
it.  The formal full endpoint is not an intermediate selected state and
does not need to satisfy a new hit condition.

Use the definitions and unchanged profile correction from
[the death-cover interface](QUARTET_GLOBAL_DEATH_COVER.md).  Its potential
is `V(U,Aprev)=A^|U| B^|Aprev|`.  Let `h(a,r)` be its local weighted action,
where `a<=r` counts old allowed and old live positions in a fresh quartet.
For a new parameter `0<t<=1`, suppose

```
h(a,r)/(A^r B^a) * t^(r/4-a) <= lambda
    for all 0<=a<=r<=4.                              (2)
```

For any fixed fresh partition, multiplying (2) over its quartets yields

```
K V(U,Aprev) <= lambda^128 V(U,Aprev) t^(|Aprev|-|U|/4)
              <= lambda^128 V(U,Aprev).              (3)
```

The second inequality uses (1).  The initial potential and free endpoint
are unchanged, so the full moment is bounded by the same endpoint factor
times `lambda^[128(T-1)]`.  There is no external cone-marker correction.
Taking `t=1` recovers the earlier scalar comparison.  For fixed `r`, the
logarithm of the left side of (2) is convex in `a`; consequently the nine
endpoint cases `a=0,r` still suffice.

[quartet_nested_cone_scalar.py](quartet_nested_cone_scalar.py) wraps the
frozen six-variable actions and profile accounting.  It adds `log2 t` as
a seventh affine coordinate and uses a bounded convex Newton diagnostic.
Exact small tests choose rational `v` and `t=v^4`, so (2) uses only integer
powers `v^(r-4a)`.  They check reachable transitions, partition actions,
the unchanged endpoint, temporal composition, neutral-marker agreement,
and the seven gradient/Hessian coordinates.

Independent source review passed.  The parent-run checks passed 60 exact
local comparisons, 54,480 reachable-transition checks, 1,700 partition
actions, 68 endpoint comparisons, 16 temporal moments, and all seven
gradient/Hessian coordinates.  The neutral cone marker agrees with the
frozen six-variable control.  No jobs were run by the helper's author.

At `(h,R,L,f,F)=(4096,819200,838861,65536,1081344)` and death upper bound
`104857`, the parent's bounded 1,024-evaluation search returned a FLOAT
exponent of approximately `+178772.94`, with `log2 t=0`.  It did not improve
the previously available scalar bound `+177477.46`.  The finite-box lower
diagnostic was approximately `-482844`, so this result establishes neither
optimizer optimality nor a limitation of the entire cone family.  The
family contains the old control at `t=1`; taking their minimum retains the
earlier bound.  Further cone searches are parked.

Production arithmetic is floating point; neither a successful optimizer
return nor its finite-box lower diagnostic is a probability certificate.
The fixed construction and full multirow distance goal are unchanged.

The next comparison, if this inexpensive cone control remains positive,
is to average the actual quartet partition before maximizing over the
reachable populations.  The valid rows then satisfy `u=0,4,...,512` and
`u/4<=a<=u`.  This keeps 24,897 population rows and does not require a
two-dimensional transition matrix.  It retains the same scalar potential,
profile correction and free endpoint.

This is the full one-step reachable count domain for the integrated-birth
relaxation.  Starting from empty `U`, select `u/4` quartets and choose `a`
birth coordinates hitting each quartet.  Every listed `(u,a)` is obtained.
Any further universal restriction on the current two counts alone would
exclude valid relaxed transitions.  Further gains must use cumulative
profile information, transition dependence, or a different potential.
