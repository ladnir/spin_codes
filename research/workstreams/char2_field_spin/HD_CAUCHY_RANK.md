# State-rank retention for the cheaper Cauchy inner

The quiet/singleton refinement is composable, but it does not close the
failed 2K control. At the two saved witnesses, it gains about 80 bits for
A=16 and 7 bits at raw density .03. Both bounds remain negative. The
refinement and its tests are useful proof components, not a distance
certificate or a performance claim.

The control has GRS[256,128] outer rows, one independent uniform permutation
per coordinate region, and a width-four/state-four full-Cauchy inner.
All four outputs are retained. Input and state diagonals are independent
uniform nonzero GF65536 scalars. Initial state is zero and terminal state
is free. The requested source dimension is K=2^20, with rate 1/2 and a
target of 40% symbol distance. A 128-bit record carries eight scalar lanes.

## Local objects and the continuation contract

Fix an outer support, routing, and all earlier setup randomness. Impose the
selected zero equations accumulated so far. State coordinates and current
inputs are linear forms on the remaining source space V. They need not be
independent. Each input prescribed live by the support must remain nonzero
on a surviving witness. Earlier actual-nonzero state requirements are also
retained; a later restriction cannot restore a killed requirement.

Write the packet as

    X = D_x x;
    R = D_s S;
    Y = F X + A R;
    S' = R + C X.

The fixed maps are obtained from a full 8-by-8 Cauchy matrix by normalizing
all four bottom rows. This normalization makes the state-to-state block
identity. The emission block A, feedthrough F, injection C, and return
matrix F-AC have the required nonzero minors. Selected top rows remain
Cauchy rows up to invertible coordinate scalings. No claim about arbitrary
subsets of normalized bottom rows is needed.

Fix marks 0<y<=1 and 1<=tau<=q-1, with q=65536. An imposed equation space
of rank r receives factor tau^(-r). Output weight is expanded positively:

    y^wt(Y) = sum over Z subset [4]
                 y^(4-|Z|) (1-y)^|Z| [Y_Z=0].

An unselected output has no nonzero requirement. In particular, the
y branch of one coordinate does not assert that this coordinate is live.

State 0 means the witness's state is zero and those equations have been
imposed. State c in {1,2,3,4} means the witness's state is nonzero and its
state-form rank is at most c. These are caps, not disjoint rank events.
The operator acts on continuations satisfying

    v_j >= 0 for j=0,...,4;   v_1 <= v_2 <= v_3 <= v_4.

There is no ordering requirement between v_0 and v_1.

## One simultaneous state-rank envelope

Suppose the incoming state-form rank is H. Let z selected equations have
a superregular coefficient matrix on the refreshed state D_s S. Let h be
the state-form rank remaining after those equations. Then

    0 <= h <= min(H,4-z),
    E[q^h-1] <= (q^H-1)/(q-1)^z.

For each fixed nonzero state vector, superregularity gives probability at
most (q-1)^(-z) of satisfying these equations. Summing over that image
space gives the moment inequality. Cases with insufficient nonzero
coordinates have probability zero and obey the same upper bound.

A branch with h=0 cannot retain the incoming actual-nonzero state. For
each threshold j, maximize

    E[tau^(h-H) [h>=j]]

over H<=c and all laws satisfying the displayed constraints. Call the
result T_j(c,z), and set T_5(c,z)=0. These maxima decrease with j and
increase with c. Therefore

    mu_h(c,z) = T_h(c,z)-T_(h+1)(c,z)

is one positive measure whose upper tails dominate every admissible
tilted law. Summation by parts proves domination for every nonnegative
nondecreasing continuation of h. The implementation reuses the exact
two-point moment calculation in `state_rank_cap_transfer.py`.

This construction matters: unrelated upper bounds on individual successor
rank masses are not substituted for a joint law. The same measure works
simultaneously for every continuation in the stated cone.

## Quiet packets

With no live input, S'=D_s S. If z selected output equations impose rank
H-h, the successor state-form rank is exactly h. The actual-state-zero
branch remains zero with coefficient one. A live state cannot return to
zero without violating its earlier liveness requirement.

For each z=0,...,4, the live continuation receives the positive measure
mu_h(c,z), multiplied by C(4,z)y^(4-z)(1-y)^z. This retains rank depletion
across repeated quiet restrictions. It does not refresh the image to a
new unconstrained four-dimensional state.

## Singleton continuation and birth

Let the only live input be in column j. From incoming zero state, both the
emission column F_j and injection column C_j have no zero entries. Thus
all four outputs are live and the new state-form rank is one. The birth
entry is exactly y^4 into cap1.

For an incoming live state and no selected output equation, the successor
rank is at most min(4,H+1). The continuation coefficient is y^4.

Now select z>=1 zero outputs. Use one selected row i to eliminate the
diagonally scaled input X_j. Its coefficient F_ij is nonzero. The other
z-1 equations become

    (A_l - (F_lj/F_ij) A_i) D_s S = 0.

Every square minor of this residual state matrix is nonzero. Its
determinant is a full-Cauchy minor containing the pivot row and column,
divided by F_ij. Expose D_s and these state equations first. Their tilted
residual-rank law is dominated by mu_h(c,z-1).

Next expose the one fresh input diagonal in the pivot equation. If the
input form and the remaining state expression are independent, its
equation imposes rank one. If they are dependent, at most one nonzero
diagonal value preserves input liveness. The expected surviving factor
is therefore at most tau^(-1), since tau<=q-1.

The pivot expresses X_j using residual state forms. Thus the successor
cap is h, not h+1. Further history restrictions can only discard branches.
The resulting term is

    C(4,z)y^(4-z)(1-y)^z * mu_h(c,z-1)/tau

into successor cap h. This argument does not assume the input form is
independent of the old state.

## A singleton return uses the incoming cap

Analyze return constraints before output constraints for this branch.
Condition on the fresh input diagonal and call the resulting required-live
form X. Let c_j be the relevant column of C. Return requires

    D_s S = -c_j X.

Every entry of c_j is nonzero. Set

    h = dim(span(S_1,...,S_4,X)/span(X)).

On the affine slice X=1, the image of S is an affine space of dimension h.
Projection onto some h coordinate positions is injective. Consequently it
contains at most (q-1)^h points with every coordinate nonzero. Meanwhile,
-D_s^(-1)c_j is uniform on (F*)^4. Every feasible diagonal choice imposes
exactly h equations and leaves X live. Hence the return rank moment is
at most

    (q-1)^(h-4) tau^(-h).

Because h<=H<=c and (q-1)/tau>=1, a uniform cap bound is

    (q-1)^(c-4) tau^(-c).

On return, Y=(F-AC)_j X. All four coordinates are nonzero, so multiply
by y^4. This return coefficient is nondecreasing in c. At cap4 it is the
old y^4/tau^4 bound; at smaller caps it can be much tighter.

The return and live-continuation bounds use different exposure orders.
They bound separate actual-zero and actual-nonzero successor covers.
Adding them is valid; no independent randomness between these covers is
asserted or required.

## Multi-input packets and composition

For two or more live inputs, the implementation uses the existing safe
two-state full-Cauchy envelope. Every nonzero successor is placed in cap4.
Its live-input rows are identical for all incoming positive caps. This
loses information but needs no new multi-input rank argument.

Quiet and singleton live continuations increase with the incoming cap,
by the simultaneous tail construction. Their return coefficients also
increase with cap. The multi-input fallback has constant positive-cap
rows. Therefore every local matrix maps the stated continuation cone
into itself and has only nonnegative entries.

Products, nonnegative sums, and positive matrix-polynomial coefficients
preserve this property. The all-ones terminal continuation belongs to the
cone. Thus the same local domination can be composed through regional
occupancy sums and the existing source/routing bounds. Free terminal
state is included; no terminal-zero constraint is inserted.

## Exact checks and bounded replay

Four tests passed. The local test uses a full-Cauchy 2/2 analogue over
GF8, exhausts its current nonzero input/state diagonals, and enumerates
every selected-output subset. It checks 1,224 joint cone comparisons.
Fixtures include rank-one and rank-two state images, zero state
coordinates, inputs inside and outside the state span, and an extra
earlier live-history requirement. Marks include tau=1 and tau=q-1.

Other tests check the required small Cauchy and Schur minors, positivity
and cone preservation of the actual 4/4 operator at GF256/GF65536 marks,
and matrix-polynomial multiplication against direct positive products.
The analytic proof supplies the general-field claim; the finite tests
are regressions, not exhaustive GF65536 evidence.

The following comparison reuses the saved natural-log marks in
`HD_FAST_PROOF.md` without optimization. These are floating upper-bound
margins, not directed certificates or observed distances.

| K20 class | Old margin | Refined margin | Gain |
|---|---:|---:|---:|
| A=16 | -2515.661491 | -2435.316855 | 80.344636 |
| Raw weight 62915, approximately .03 density | -220206.823065 | -220200.205239 | 6.617826 |

The sparse marks are (6.217618013797049, 2.0114849606107517).
The middle marks are (-2.513585167041394, -0.07889607306903107,
2.805601752915663, 0.27462567642365154). The source factors and forbidden
output cutoff are unchanged.

```text
python -B research/workstreams/char2_field_spin/theory/hd_cauchy_rank_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_cauchy_rank.py --replay
```

No cone obstruction appeared for this bounded refinement. The remaining
gap is numerical, and large: quiet/singleton rank retention alone does
not repair the low-middle bound. Multi-input packets still erase the rank
history by sending every continuation to cap4. A further proof step would
need a joint multi-input successor law or a source/state cut argument,
not another independent per-entry maximum. The present receipts do not
justify promoting this 2K construction as a 40%-distance code.

Only these new proof files were changed. No native build, benchmark,
remote command, parameter optimization, or frozen-certificate edit was run.
