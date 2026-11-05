# Fill occupied allowed runs before counting geometric witnesses

2026-10-31. In the structural low-death class, state patterns within one
allowed coordinate run need not be counted separately. Filling an
occupied run preserves its exact live selectors and births, increases
volume, and does not increase deaths. This gives a normal form for the
geometric witnesses that would force a bad codeword by dimension.

The construction, once-sampled setup, and free terminal remain unchanged.
This note is deterministic after fixing the routing and quartet
partitions. It transforms support masks, not field-compatible amplitudes.

## The normal form

Fix all routing permutations and quartet partitions. For a state support
history `S=(S_1,...,S_T)`, set `S_0=empty` and define

```
U_t = union of the update-t quartets meeting S_(t-1),
B_t = S_t \ U_t,
U_(T+1) = [512],
A_t = (U_t intersect U_(t+1)) union B_t.
```

The final set is formal: it adds no partition, update, or live-volume
charge. These definitions imply

```
S_t subset A_t,       A_(t-1) subset U_t,
A_t \ U_t = B_t.
```

An allowed run is a maximal consecutive interval of one coordinate in
`A`. It is closed if it ends before `T`, and open otherwise. It is forced
if its first position belongs to `B`. Every birth starts a distinct
allowed run, and every forced run meets `S`.

Let `Sat(S)` fill every allowed run that meets `S`. Let `Sat_open(S)`
also fill every open allowed run, including those disjoint from `S`.
Both transformations use the allowed sets of the original exact pair
`(U,B)`.

Write

```
L(S) = sum_t |S_t|,
D(S) = sum_(t=1)^T |S_(t-1) \ S_t|.
```

There is no death after `T`.

**Saturation theorem.** For either transformed support `S+`,

```
S subset S+ subset A,
U(S+) = U(S),       B(S+) = B(S),
L(S+) >= L(S),      D(S+) <= D(S).                    (1)
```

Every allowed run in `S+` is either completely occupied or completely
empty. For `Sat_open(S)`, every forced or open run is occupied. The
transformation is idempotent because `(U,B)`, and hence `A`, is preserved.

To prove selector preservation, consider any update `t`. Every selected
quartet still meets `S+_(t-1)`, since it already met `S_(t-1)`. Conversely,
`S+_(t-1) subset A_(t-1) subset U_t`, so an unselected quartet cannot
acquire a nonzero incoming position. This proves exact equality of the
live selectors, not merely containment.

Every birth remains present because `S subset S+`. The equality
`A_t\U_t=B_t` prevents any new birth, proving preservation of `B`.

Each closed allowed run that meets the original support contains at
least one death. After filling, it contains exactly one death, at its
closing edge. These deaths are distinct across allowed runs. Each open
run contributes zero deaths after filling, even if it was originally
empty. Thus

```
D(S+) = number of closed allowed runs meeting S <= D(S).
```

Volume can only increase. This proves (1), including the free endpoint.

## The sufficient dimension obstruction is preserved

Let

```
V_t = U_t union B_t,
F = sum_t |V_t|,          f = sum_t |B_t|,
R(V) = sum_rows min(128, number of covered positions in that row).
```

The eligible row set `H` consists of rows with at least 129 covered
positions; write `h=|H|`. Since saturation preserves `U` and `B`, it
preserves `V,F,f,R,H,h` and every regional footprint volume exactly.

Recall the [dimension interface](QUARTET_DIMENSION_CUT_INTERFACE.md).
For any support mask, the supported state space has dimension `L`,
while the outer parity equations have rank at most `R(V)`. The triangular
state-to-output map is injective and has output support length at most
`L+D`. Therefore

```
L>R(V)  and  R(V)+D+1<=W                              (2)
```

force a nonzero output of weight at most `W`, for every choice of the
field maps. Its actual support may be smaller than the mask.

If `S` satisfies (2), then `S+` satisfies (2) with the same `R(V)`.
Consequently, a mask satisfying this sufficient obstruction exists if
and only if a mask in the saturated normal form satisfies it. The same
equivalence holds after requiring every open allowed run to be filled.

More generally, for fixed `(U,B)`, any thresholds `L0,D0` give the exact
existence equivalence

```
exists exact-hit S with L(S)>=L0 and D(S)<=D0
  iff
exists saturated exact-hit S with all open runs filled,
       L(S)>=L0 and D(S)<=D0.                         (3)
```

Only the forward direction needs proof; the reverse direction is set
inclusion. Equation (3) also holds when `(U,B)` ranges over any class
defined using the preserved geometry.

This does not preserve exact state volume, an upper-volume constraint,
or the output weight of a particular amplitude assignment. No claim is
made that filling a codeword's support gives another codeword. The
dimension argument establishes the sufficient obstruction separately
for the transformed mask.

At the current diagnostic profile

```
(h,R,L,f,F) = (4096,819200,838861,65536,1081344),
```

all original masks with `D<=124517` satisfy (2). A geometric exclusion
may therefore count saturated masks with the same `(h,R,f,F)`, volume
at least `838861`, and at most `124517` deaths. For the smaller original
slice `D<=104857`, use the smaller death threshold. The transformed
volume need not remain `838861`.

## Remove ineligible births before applying a global route bound

The route factor used in the quotient requires every birth to lie in an
eligible row. Actual bad-codeword supports already satisfy `B subset H`.
An arbitrary dimension mask need not satisfy that condition. The following
pruning step supplies it when proving absence of all obstructions (2).

Suppose a birth position `(t,i)` lies in a row whose current footprint
degree is at most 128. Delete that position from `S_t`, and recompute the
exact selectors, births, and footprint. Denote the new statistics by
`L',R',D'`.

Since `i` is outside `U_t`, its preceding state position is zero. Removing
this `0`-to-`1` start cannot add a death; it removes the following death
if the next state position is zero. The same conclusion holds at the
free terminal. Consequently `D'<=D`.

The new support is a subset of the old support, so the footprint can only
shrink. At update `t`, the deleted birth slot leaves the footprint because
`U_t` is unchanged and did not contain it. Its row had degree at most 128,
so that row's contribution to `R` drops by at least one. Thus

```
L' = L-1,       R' <= R-1,       D' <= D,
L'-R' >= L-R,   R'+D' <= R+D.                         (4)
```

Both inequalities in (2) persist. Repeat this step whenever an ineligible
birth remains. Each step strictly decreases `L`, so the process terminates.
The preserved positive surplus `L-R` prevents termination at the empty
mask. Recomputing is essential: removing one position can create later
births or make another row ineligible.

At termination every birth lies in the current eligible set. Saturation
then preserves that inclusion and all current footprint statistics.
Therefore an arbitrary obstruction (2) exists if and only if a
birth-compatible saturated obstruction exists, with every open allowed
run filled.

Pruning may change `h,R,L,f,F`; it does not preserve a specified profile.
A global union over the final profiles can use the birth-compatible route
factor. A fixed-profile argument starting from an actual bad-codeword
support needs no pruning and retains its original `h,R,f,F` throughout
saturation.

## What remains to be selected

For a fixed `(U,B)`, let `F_c,F_o,U_c,U_o` be its forced closed, forced
open, unforced closed, and unforced open allowed runs. Let `ell_r` be
the length of run `r`. The normal form fills the baseline

```
J = union of all runs in F_c union F_o union U_o,
c_f = |F_c|,                L_base = |J|.
```

Its only remaining choices are a subset `C subset U_c`. The resulting
support is determined by this subset, and its statistics are

```
S(C) = J union union_(r in C) r,
L(S(C)) = L_base + sum_(r in C) ell_r,
D(S(C)) = c_f + |C|.                                 (5)
```

Exact hitting must still be imposed. For every selected quartet at
update `t`, its incoming positions at time `t-1` must meet `J` or a run
in `C`. Quartets already met by `J` need no further condition. Every
other quartet gives a clause involving at most four closed unforced
runs: at least one of those runs must belong to `C`.

Thus (3) becomes a selection problem on whole runs: satisfy these hit
clauses, choose at most `D0-c_f` closed unforced runs, and obtain volume
at least `L0-L_base` from them. There is one support per chosen subset;
patterns within an occupied run create no further multiplicity. Dropping
the hit clauses would return to the packing relaxation and is not part
of this reduction.

At most 512 forced runs can be open. Hence `c_f>=f-512`, and the number
of chosen closed unforced runs is at most `D0-f+512`. At the primary
profile this is 39,833 for `D0=104857`, and 59,493 for the full structural
threshold `D0=124517`. These bounds concern whole-run choices, not the
number of positions those runs may occupy.

## A positive transfer can enforce the normal form

The chronology agrees with the four coordinate labels proposed for the
next counting operator. Before choosing `U_(t+1)`, the labels are

```
Z: outside U_t,
O: in U_t outside A_(t-1),
E: in A_(t-1) outside S_(t-1),
S: in S_(t-1).
```

If the new quartet is selected, old `S` positions remain occupied and
old `E` positions remain empty. Each old `O` position starts a new
unforced allowed run and chooses its occupancy once. A birth outside
`U_t` starts an occupied forced run. Requiring at least one resulting
occupied position enforces the exact new live selector.

If the new quartet is unselected, its allowed runs close. Charge one
death for each old `S` position and none for an empty run. These updates
preserve the chosen occupancy throughout each allowed run.

For the version that fills every open run, the terminal rejects any old
`E` position and fills every old `O` position. Old `S` positions remain
occupied; births outside `U_T` remain optional. For a state-volume marker
`s>0` and birth marker `w>0`, populations `(u,o,e,a)` therefore have
terminal weight

```
1_(e=0) * s^(a+o) * (1+w*s)^(512-u).                  (6)
```

Here `a` is the old `S` population. The terminal has no hit condition,
new live-volume charge, or death. Allowing empty open runs instead is a
valid looser normal-form family, with a different terminal.

The geometric moment uses a lower-volume threshold, so its volume
Chernoff marker must satisfy `s>=1`; a death upper bound uses `0<d<=1`.
It does not reuse a fixed-exact-volume field exponent for the original
support. Routing conditions based only on the preserved footprint can
be applied before enlarging this positive count.

## Status

Two independent derivation reviews passed saturation, exact selector and
birth preservation, death monotonicity, and the free terminal. The
positive run-cover count and its averaged operator are separate work.
No numerical job or probability certificate is reported here. The next
step is to bound this smaller family while retaining its exact hits and
regional routing constraints.
