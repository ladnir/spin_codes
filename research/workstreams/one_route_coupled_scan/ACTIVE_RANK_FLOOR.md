# Retain a lower bound on the outgoing state rank

The active-step continuation currently caps the outgoing rank from above.
Allowing every smaller rank loses information even when no zero equation
was imposed. A fresh pair, for example, cannot produce rank one without
eliminating at least one of its independent input images.

The bounds below apply before maximizing a terminal function over outgoing
ranks. They change only the proof continuation, not the construction.

## 1. Count zero equations, not only coefficient charges

Fix one step, its current coefficients, and one retained support branch.
Write I,H,J,H' for the input, incoming-state, output, and next-state supports.
Let k=|I|. The recurrence is

    Y=X+DQ,             Q'=RQ+Y=(R+D)Q+X.

Assume R+D is invertible, as in the construction. Let W be the common span
of zero forms before this step, U=span(Q), and

    r=dim((U+W)/W),
    c=|I minus (J intersect H')|,
    d=|(H minus I) minus H'|.

The letter d here counts possible old-state zero equations; it is not the
global remaining message dimension from GLOBAL_RANK_GATE.md. All earlier
nonzero requirements remain enforced. In particular, the declared nonzero
coordinates of X and Q cannot be killed by a retained branch.

Let Z be the span of the zero equations on Y and Q' imposed in this step,
W'=W+Z, and g=dim(W'/W). Then

    0 <= g <= c+d.                                                (1)

At an active coordinate j, at least one of Y_j,Q'_j is required zero
exactly when that coordinate contributes to c. If both are zero, their
difference is (RQ)_j. A retained support branch requires (RQ)_j=0 there,
so the two equations are identical modulo W. They contribute at most one
rank, not two.

At an inactive coordinate, Y_j=D_jj Q_j. Its zero/nonzero status is already
determined because D_jj is nonzero. Only Q'_j can impose a new equation,
and this requires j in (H minus I) minus H'. Coordinates outside I union H
give no new equation: Q'_j is then a rotated old coordinate. This proves (1).

## 2. Arbitrary input forms

Put b=dim((U+span(X)+W)/W), and let s be the outgoing state rank modulo W'.
Every equation in Z lies in U+span(X). Therefore its residual dimension is
b-g. Invertibility gives the exact identity

    span(Q',X)=span(Q,X),

both before and after quotienting. Since the residual span of X has rank
at most k, b-g<=s+k. Also b>=r. Consequently every retained branch satisfies

    s >= max(0, r-k-c-d),
    s <= min(|H'|, r+k).                                         (2)

These bounds hold for fixed input forms and for inputs exposed through fresh
outer coefficients. Fresh coefficients alone do not imply fresh source
dimensions. In particular, c must not be subtracted from the upper cap in
that latter setting.

## 3. Inputs independent of the old source span

Now assume the k active input forms are linearly independent modulo U+W.
This holds in the proved complete-column fresh prefix. Then b=r+k.

Project the current zero equations onto (U+span(X)+W)/(U+W). The equation
at each charged active coordinate has its unique input image as pivot.
Different such coordinates have independent images; inactive equations
project to zero. Thus their projected rank is exactly c, regardless of
dependencies among the old-state equations. In particular g>=c, and the
surviving input-image rank modulo U+W' is exactly k-c.

The identity Q'=(R+D)Q+X shows that Q' contains these surviving images.
Together with (1), this gives

    s >= max(0, k-c, r-c-d),
    s <= min(|H'|, r+k-c).                                      (3)

One may also record the exact old-state rank after this step as r-(g-c).
The quantities c,d are not substitutes for the unknown actual gain g.
Equation (3) needs only their proved bounds on g and the projected rank c.

For a live outgoing support, maximize the terminal function only over

    max(1, lower bound) <= s <= upper bound.

An empty interval of ranks contributes zero. For H'=empty, permit the
quiet terminal only when the lower bound is zero. From a quiet incoming
state, a nonempty fresh input gives Y=Q'=X and rank exactly k. An arbitrary
input instead has some rank from one to k; it never remains quiet.

## 4. Apply the floor before aggregating branches

The current fresh and lazy-outer term lists retain c for each output-mask
branch, even when several terms share H'. Compute d from H,I,H', apply
(3) or (2), and only then sum their terminal bounds.

The legacy ordinary matrix has already summed over J. Its stored zero in
the gain field is a placeholder, not evidence that c=0. Either retain c
when forming those terms or use the safe aggregate upper bound

    c_max=|I intersect (complement(H') union (H intersect R(H)))|.  (4)

Indeed, an active coordinate with nonzero Q'_j can have Y_j=0 only when
both Q_j and (RQ)_j are nonzero. Replacing c by c_max in (2) gives a valid,
possibly weaker floor for the aggregated ordinary entry.

[test_active_rank_floor.py](theory/test_active_rank_floor.py) checks these
relations on bounded exact GF16 quotient fixtures, using the independent
Gaussian oracle from the paired-active tests. Such checks are not a
whole-code certificate. The next diagnostic is the change in the actual
five-rank regional continuation, with all other parameters fixed.
