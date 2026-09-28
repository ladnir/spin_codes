# Excluding impossible mixture points from the moment bound

The affine moment argument bounds a convex function at a cell's vertices.
Some rectangular vertices cannot represent any permitted mixture. Including
them is valid but can force much finer subdivision than the actual domain
needs. This refinement evaluates vertices of a tighter rational enclosure.
It changes only the proof, not the SPIN construction.

## A polytope containing every permitted mean

Use the component features f_i and active labels b_i from
[LIFTED_MIXTURE.md](LIFTED_MIXTURE.md). The only inactive component has
feature zero. Let rho=q_min/G. For any composition with active fraction
a>=rho, its mean can be written as m=a*v, where v is a convex combination
of active component features. Therefore m belongs to

    P = convex hull { rho*f_i, f_i : b_i=1 }.

Indeed, for rho<1 write a as a convex combination of rho and 1, then apply
the same combination to each active feature in v. For rho=1, m=v directly.
This argument covers all allowed integer compositions and also includes
fractional compositions, which only weakens the bound.

The implementation obtains candidate supporting planes from a floating
convex-hull calculation. It constructs each proposed plane from three
source points with exact rational arithmetic. It retains the inequality
only after checking it against every source point. Thus their intersection
H contains P; facet completeness and numerical hull accuracy are not
assumptions of the proof. Missing or rejected facets make H larger.

For each rectangular cell C, compute the vertices of C intersect H.
The implementation enumerates every triple of constraint planes with
independent normals, solves their intersection exactly, and checks every
inequality. The box makes the intersection bounded. Every vertex has
three independent active constraints, including for a lower-dimensional
intersection. Degenerate cells and empty intersections are covered by tests.

## The same affine moment argument on the smaller domain

Fix one output tilt lambda, the positive anchors, and the affine upper
weights W(m) from [AFFINE_MOMENT.md](AFFINE_MOMENT.md). The function
phi(m)=log H_lambda(W(m)) is convex. For any rational slope u and outward
vertex bounds B_v, let A=max_v(B_v-u dot v). Then

    phi(m) <= A + u dot m

throughout C intersect H. Every actual composition in C lies in this
intersection. Thus the same outer generating-function sum remains valid;
only its inner affine upper bound has changed.

The floating least-squares fit merely proposes u. Replay reconstructs
all supporting inequalities and vertices, and recomputes B_v with Arb.
New witnesses mark this option with `clipped`; old witnesses keep their
original rectangular vertices. Positive anchors are also chosen using
the clipped vertices, avoiding large values caused by impossible corners.

This refinement is opt-in via `--clip-domain` and remains under evaluation.
On the larger q_min=512 diagnostic rectangle from `AFFINE_MOMENT.md`,
the first clipped proposal was slightly worse than the existing proposal
(log2 upper +5101 versus approximately +5037). Thus it has not replaced
the default search. The smaller diagnostic rectangle still closes below
-293.04967. Neither diagnostic certifies the full dense range.
The completed 8% certificate still replays with its original witnesses.
