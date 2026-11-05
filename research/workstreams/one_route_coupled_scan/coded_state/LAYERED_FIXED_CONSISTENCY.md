# Layered inner: fixed-coordinate consistency

Status: local counting lemma and exact toy checks; **no whole-code certificate**.
The layered construction, its setup distribution, and the certified GF32
sandwich control are unchanged. This continues [mixed histories](LAYERED_MIXED.md).

## The joint route bound does not fix the retained class

For a fixed subset S of s outer blocks, condition on the regional head sets.
If region j contains h_j heads from S, its number of internal matching edges
has distribution Hypergeom(M/2, s-h_j, h_j). Regional matchings remain
independent. Writing h=sum_j h_j and n=M/2, their total conditional mean is
at most (R*s*h-h^2)/(R*n). Thus head bias and internal edges cannot independently
attain their worst values.

To bound the conditional tail, reveal either all selected heads or all selected
tails, whichever gives fewer draws in total. Write d=min(h, R*s-h). The sum of
squared martingale range lengths is at most both d and
((n+1)*R*d-d^2)/(R*(n-1)). Indeed, for d_j draws in one region the exact range
sum is sum_{i=1}^{d_j} ((n-d_j)/(n-i))^2 when d_j<n, and zero when d_j=n.
It is at most d_j*(n+1-d_j)/(n-1); sum and use convexity. The conditional tail
bound is uniform over the conditioned head sets, so no additional union over h
is needed. The integer slack in `layered_joint_route.py` allocates the same
2^-56 union budget over block subsets and suffixes as the matching envelope.

The retained profile has R=192, f=589632, r=0, e=368520, and B=49136.
Combining this conditional tail with the head bound leaves its dimension
upper bound unchanged at 163923 and its best log2 failure upper bound at
+734997. The maximizing subset size is 2562. This is a failed bound, not a
low-weight codeword. We therefore do not run broader coverage with this bound.

## What the inverse count was discarding

Fix the route, incoming-state support masks, exceptional-pair mask, and
internal output-zero mask. Let U be a regional suffix of length u=RM,
containing L=u/b batches. Work over J=GF(2^32). The shortened outer code C_U
has dimension K_R=max(0,u-K). The first incoming state is zero; the terminal
state is free.

The inverse local law partitions input coordinates into a free set F and
fixed coordinates F^c. Put D=|F|, let e count exceptional full pairs, and put
d=dim(C_U intersect J^F). The existing probability bound for each assignment
of states is at most eta*q^(d-D-e), where eta<2. It allowed the fixed values
to lie in any affine outer-code coset, even when that coset was empty.

The fixed values are especially simple. For an internal pair with incoming
state (x,y), each fixed coordinate has form Q_next+gamma*Q, gamma in {0,1}:

| Pair type / branch | First coordinate | Second coordinate |
| --- | --- | --- |
| 00 | fixed, gamma=0 | fixed, gamma=0 |
| 10 | fixed, gamma=1 | free |
| 01 | free | free |
| 11, regular | free | free |
| 11, exceptional | fixed, gamma=0 | fixed, gamma=1 |

All first-batch coordinates are fixed with X=Q_1. These formulas follow from
X=Q_next+A Q; they do not replace the actual encoder by a different map.

Let Z be the linear space of state assignments after imposing inactive-state
zeros and the prescribed internal equations Q_next=delta*Q at zero outputs.
We may omit nonzero inequalities to obtain an upper bound. Write a=dim Z.
Let L_fix:Z -> J^(F^c) extract the fixed input values, and define
kappa=dim ker L_fix. This kernel is an auxiliary counting object; actual
fixed input coordinates are not assumed to be zero.

**Consistency lemma.** The number of assignments in Z whose fixed values
can belong to the outer code is at most

    q^min(a, kappa+K_R-d).

Proof: the projection of C_U onto F^c has dimension K_R-d. Every nonempty
fiber of L_fix has dimension kappa. Intersecting its image with that outer
projection has dimension at most K_R-d, and the preimage is also contained
in Z. This proves the bound.

Multiplying by the per-assignment inverse probability gives the improved
amplitude bound

    eta * q^min(a+d-D-e, kappa+K_R-D-e).

Support, branch, and zero-mask counts still have to be included. This is a
fixed-mask lemma, not a minimum-distance theorem. The route event applies
simultaneously to the relevant masks; the remaining mask sum cannot be skipped.

If only d<=d_bar is known, use

    eta * q^min(a+d_bar-D-e, kappa+K_R-D-e).

Do not insert d_bar into the apparent rank saving
max(0,a-kappa-K_R+d): that would turn an upper bound on d into an unjustified
lower bound on the saving. The second exponent above avoids that mistake.

## Computing kappa without large matrices

All equations defining ker L_fix involve a single spatial coordinate at
adjacent times. Inactive vertices and fixed gamma=0 equations pin vertices
to zero. A zero-output edge identifies adjacent vertices with ratio delta;
a fixed gamma=1 equation identifies them with ratio 1. If both occur on the
same edge, their component is pinned because delta!=1. Otherwise the graph
is a collection of time paths, with no cycles other than these parallel edges.
Therefore kappa is exactly the number of unpinned connected components.

One can scan each path with a single live-component bit. On an unlinked edge,
finalize the previous component and start the next; on a linked edge, preserve
the live bit unless a pin kills it. Count the final component at the end.
This yields a two-state counting interface for summing histories. It does not
yet provide an upper bound across all admissible support/exception masks.

The elementary cap kappa<=D follows by extending the fixed-coordinate map
to all coordinates: on free coordinates set the auxiliary value to Q_next.
The full map from states to auxiliary inputs is triangular with identity
diagonal. Setting its fixed coordinates to zero leaves D dimensions before
the additional support and output-zero constraints.

Fixed equations cannot simply be counted as independent constraints. Along
one zero-output run, two fixed coordinate-zero equations can impose only one
additional constraint after the temporal equations. The test includes this
negative control.

## Numerical target and exact checks

For the retained hard class, a<=790492, D=638768, e=368520, K_R=524288,
and d_bar=163923. The existing support/branch/zero-mask count contributes
at most 2426932 bits. Thus a **uniformly proved** kappa<=407155 for that
class would give

    2426932 + 1 + 32*(407155+524288-638768-368520) = -107.

This is a required cap, not an established one. It addresses one retained
class at the exploratory 40% threshold, not all classes or the near-half goal.

Run serially with Python -B:

    python -B research/workstreams/one_route_coupled_scan/coded_state/layered_joint_route.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/layered_fixed_rank.py

The first checks 42 exact range inequalities and 783 conditional
hypergeometric probabilities, then reproduces the unchanged hard-class bound.
The second uses GF4 and routed RS[4,2] blocks: 256 fixtures check the fixed
input formulas, actual outer projections, and consistency dimension bound.
Two fixtures have a positive additional rank saving. Another 864 exhaustive
path patterns compare component counts and the two-state scan with Gaussian
rank. These checks test the lemma and its implementation, not K20 distance.
Both commands are bounded and retain no raw output files.

Next: count or bound unpinned temporal components jointly with exceptional
pair locations and support changes. Test the retained hard class first.
If the resulting count cannot meet the cap, inspect which long components
survive instead of treating every fixed equation as a new independent root.
