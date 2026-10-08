# Return bounds without enumerating zero state coordinates

The new [partial-live moment helper](theory/partial_live_rank.py) bounds a
packet's return to zero by counting vectors whose input coordinates are
all nonzero. It allows any pattern of zero state coordinates. This removes
the zero-flat enumeration from the return bound; it does not change the
encoder or provide a whole-code certificate.

## The local interface

Fix a source-form image space \(U\leq\mathbb F_q^{w+s}\) of dimension
\(d\), before exposing the current packet's independent nonzero input
and state diagonals. Its first \(w\) coordinates are the active input
positions. The other \(s\) coordinates are the incoming state. The
dimension of this joint image need not equal \(w+s\).

Let \(K\leq U\) be the kernel of the equations under consideration after
sampling those diagonals, and write \(j=\dim K\). Let \(E\) mean that
\(K\) contains a vector with all \(w\) input coordinates nonzero. For a
rank tilt \(1\leq\tau\leq q-1\), the required quantity is
\[
  \mathbb E[\tau^{j-d}\mathbf 1_E].
\]
An actual live-state return also requires a nonzero incoming state.
Dropping that requirement only enlarges the event we bound.

The helper takes two facts from the packet analysis: each fixed vector
with live inputs satisfies the equations with probability at most
\(p_0\), and every kernel has dimension at most \(J\). Numerical matrix
rank alone does not prove the first fact when state coordinates may vanish.

## Counting and the rank moment

Let \(M\) count vectors in \(K\) with all input coordinates nonzero.
Linearity of expectation gives
\[
  \mathbb E[M]\leq(q^d-1)p_0.
\]
On \(E\), each of the \(w\) designated coordinate functionals is
nonzero on \(K\). Its zero set has \(q^{j-1}\) vectors. All these zero
sets contain zero, so a union bound gives
\[
  M\geq L_j:=q^j-1-w(q^{j-1}-1)>0\qquad(1\leq w<q).
\]

Use the points \((0,0)\) and
\((L_j,\tau^{j-d})\), for \(1\leq j\leq\min(d,J)\). The largest
expected second coordinate subject to an expected first coordinate at
most \((q^d-1)p_0\) is an upper bound on the required moment. There is
one mean constraint and one total-mass constraint, so an extremizer uses
at most two points. The helper checks all such pairs and all feasible
vertices. This also handles the possible nonconcavity of the first
segment when \(\tau\) is close to \(q-1\).

A bound on \(d\), rather than its exact value, is handled by maximizing
over every permitted dimension. No monotonicity assertion is needed.
All equations using the current diagonals must be included in this one
moment. Multiplying by another charge derived from those same diagonals
would reuse their randomness.

## Application to packet returns

For a pure return, the equations are \(R+CX=0\), with superregular
\(C\in\mathbb F_q^{s\times w}\). Fix an input-live vector whose state
support has size \(h\). The restricted matrix \([C\mid I]\) has rank
\(\min(s,w+h)\). If \(w+h\leq s\), the equations cannot hold for this
vector. Otherwise, conditioning on the nonpivot diagonal entries fixes
each of \(s\) pivot entries to at most one value. Hence
\(p_0=(q-1)^{-s}\). The unrestricted matrix has kernel dimension \(w\),
so take \(J=w\). This argument needs no cover of the zero state positions.

For the [full-Cauchy packet](GF8_MDS_PACKET.md), adding \(z\) selected
output-zero equations gives the joint rank
\(\min(s+z,w+h)\). Thus take
\(p_0=(q-1)^{-(s+z)}\) and \(J=w-z\); when \(z\geq w\), no live
solution exists. A joint image-rank cap \(c\) on the incoming state gives
\(d\leq c+w\). This joint bound can be used directly in the positive
output-weight expansion.

## Validation and scope

The [exact tests](theory/partial_live_rank_tests.py) enumerate every
subspace of \(\mathbb F_4^3\), every nonzero diagonal, both splits into
input and state coordinates, and three rank tilts. All 258 resulting
moment comparisons pass. Additional tests cover kernel caps, the
near-endpoint interpolation, and a low-rank return case. A second exhaustive
test uses a two-input, one-state full-Cauchy packet over \(\mathbb F_8\).
It checks 1,503 moments including joint output-zero and return equations,
all source-image subspaces, and every nonzero diagonal; all pass.

These are local algebraic checks. Sparse and middle class calculations,
outward rounding, complete coverage, and an independent audit of the
assembled certificate are still required for the three-quarter-rate goal.
