# Expose the regional route online

The fixed-route envelopes collapse the live-state information before
averaging insertion positions. That collapse can be avoided within a region.
The route can instead be exposed batch by batch, and the suffix bound can
retain the current residual rank and support.

This note proves the online recurrence for ordinary chronological equations
and the complete-column fresh prefix, then checks its lazy-outer exposure
conditions. It does not assert a full-code certificate. The combined numerical
screen remains floating-point diagnostic evidence only.

## 1. Fix the source class, not a message chosen after setup

First fix the outer realization and an exact outer-support class, independently
of the route. Let V be its shortened message space, of dimension D over F.
The outer source forms and their declared zero/nonzero coordinates are now
fixed. Initial nonzero requirements may include every declared active source
form: they do not depend on the route.

Retain the common zero-form space W and all nonzero requirements from
[PROOF.md](PROOF.md). Reject a branch if W contains a required nonzero form.
Use the same global potential

    Phi = 1[branch retained] tau^(-dim W),  tau>=1.

A support history records the zero/nonzero choices for the successive output
and next-state coordinates. The proof sums over all such histories. It does
not condition on the history of a message selected after seeing the setup.

For each fixed prefix of support labels, its W, exposed state forms, and
retained-prefix indicator depend only on the revealed route prefix, exposed
inner coefficients, and the fixed outer data. Ordinary zero equations use
no unconsumed source position or future inner coefficient. Future source
nonzero requirements, when included initially, are fixed independently of
the route and do not change this fact.

Consequently, conditional on this full revealed prefix, the remaining pair
permutation is uniform over the remaining labeled source pairs. Conditioning
only on its coarse rank/support summary is unnecessary. The recurrence below
first establishes a uniform bound for every full prefix with that summary.

This handles messages chosen after setup. For each complete setup, every
low-weight nonzero message supplies one of the enumerated support histories.
That history has rank at most D-1. Therefore, for a free terminal and output
threshold e, the same pointwise implication gives

    1[there exists a qualifying message of weight <=e]
      <= tau^(D-1) z^(-e)
           sum_histories Phi_final z^(output weight),   0<z<=1.       (1)

Expectation and the finite nonnegative sums can be taken in either order.
No argument assumes that the witness message was fixed before the setup.

## 2. A suffix bound uniform over the remaining forms

Within a width-four region, let m be the number of complete batches still
unconsumed. There are 2m remaining physical pair slots. Write

    n=(n_0,n_1,n_2),  n_0+n_1+n_2=2m,

for the counts of remaining source pairs with zero, one, and two declared
active coordinates. Keep the current residual state rank r and support H.

For a nonnegative terminal function f(r,H), define V_m(n;r,H) to upper-bound
the conditional expected sum over all suffix support histories, weighted by

    tau^(-(dim W_final-dim W_now))
       z^(suffix output weight) f(r_final,H_final).                  (2)

The bound must hold for every retained full prefix having the indicated
counts, rank, and support. In particular, it is uniform over:

- the actual remaining labeled source forms and all their correlations;
- the current W and every earlier nonzero requirement;
- the revealed route labels, swaps, and inner coefficients.

The source forms are not resampled to become independent or uniform. Their
values can depend on the past. Only the unrevealed route labels and fresh
inner coefficients retain their original conditional distributions.

At m=0, n=(0,0,0), set V_0(n;r,H)=f(r,H). Invalid count profiles contribute
zero. The initial regional profile in the present two-group gate has either
one full pair or two singleton pairs, so every batch has at most two active
coordinates. The existing active-step operator covers those masks.

Let T_I be a sound one-step continuation operator for fixed input support I.
For I=empty, use the empty-step operator. Otherwise use either the ordinary
active bound or the [paired-rank operator](PAIRED_RANK_GATE.md). These operators
are uniform over the actual forms and common quotients in their stated scope.

## 3. Exact next-batch probabilities

Expose the two ordered source-pair labels assigned to the next batch, then
their independent swap bits. For categories a,b in {0,1,2}, the probability
of that ordered category pair is

    n_a (n_b-1[a=b]) / ((2m)(2m-1)).                               (3)

There is no replacement. A category-zero pair has local mask 0, a singleton
has local mask 1 or 2 with probability 1/2 each, and a full pair has mask 3.
Let S_0={0}, S_1={1,2}, and S_2={3}. For u in S_a and v in S_b, the
width-four input mask and remaining profile are

    I = u + 4v,
    n' = n - e_a - e_b,

and their joint probability is (3) divided by |S_a||S_b|. Terms with
insufficient category counts are omitted. Individual pair identities have
the same number of permutation completions; averaging over their values
requires no additional coefficient after applying a uniform local bound.

For the full-pair profile n_2=1,n_1=0:

    Pr[I=empty] = (m-1)/m,
    Pr[I=3] = Pr[I=12] = 1/(2m).                                  (4)

The active event leaves n_2=0. Formula (4) implies a uniform active batch,
but the online recurrence does not first choose that whole future position.

For the two-singleton profile n_1=2,n_2=0, the next batch consumes both
singleton pairs with probability

    1/(m(2m-1)),                                                  (5)

one with probability

    4(m-1)/(m(2m-1)),                                            (6)

and neither with the remaining probability. Conditional on (5), masks
5,6,9,10 are equiprobable. Conditional on (6), masks 1,2,4,8 are equiprobable.
After one singleton is consumed, its remaining partner is still routed by
the exact remaining-slot rule (3).

Across the entire region, (3) reproduces the same-batch collision probability
1/(2m_initial-1) from [TWO_GROUP_ROUTE.md](TWO_GROUP_ROUTE.md). The online
formulation neither removes collisions nor replaces the slots by independent
samples with replacement.

## 4. Conditional induction gives the nonlinear recurrence

For each joint mask/count outcome above, let p_m(I,n'|n) denote its exact
probability. Define

    V_m(n) = sum_(I,n') p_m(I,n'|n)
                  T_I(V_(m-1)(n')).                             (7)

The vectors on the right are indexed by the next residual rank and support.
Their operators preserve the same global rank potential; no rank is reset
when a pair or batch is consumed.

To justify (7), fix an arbitrary retained full prefix and expose the next
two labels and swaps. Their specific input forms are now fixed. Conditional
on every resulting inner-step branch, the induction hypothesis bounds the
remaining expected history sum by V_(m-1)(n'), uniformly over its new forms
and W. Apply the fixed-input one-step bound T_I to this terminal function.
Then average over the revealed labels and swaps. Its uniformity removes
the dependence on the selected form values, leaving exactly p_m(I,n'|n).

The next-mask law is unchanged by the current W or support labels because
these use only the revealed prefix. This assertion is about each fixed
enumerated prefix; it is not conditioning on the existence of a future
message satisfying that prefix. Additional branch rejection can only lower
the history sum.

This is not an interchange of T with a fixed-route expectation. In the
earlier argument, E_route[T(f_route)] cannot generally be replaced by
T(E_route[f_route]). Here the induction hypothesis supplies one suffix
bound valid for every possible actual prefix, after conditional expectation
over its still-unrevealed route. Equation (7) follows directly from that
uniform statement and the tower property. The two arguments have different
intermediate objects.

After the region, one may collapse this rank/support operator to a common
two-state linear envelope by evaluating the quiet and live basis terminals
separately. Monotonicity, positive homogeneity, and subadditivity justify
their nonnegative linear combination. This delays the scalar collapse until
the whole region. It does not justify an arbitrary nonlinear matrix product
between regions or an unproved exchange with the outer-support sum.

## 5. Fresh-prefix scope and forbidden future dependence

For the complete-column source class, the first h-32 occupied columns give
independent source coordinates. These are the first occupied columns, not
the first physical columns. Their unused input forms stay independent of
every earlier source form and ordinary zero equation, including equations
formed during preceding empty batches. Thus the fresh-input version of T_I
can be used throughout those columns in (7).

The input counts always describe the declared exact source support. If an
exposed form is forced zero in W despite its required nonzeroness, reject
that branch. Do not silently remove it from the active-pair counts and
continue with a different support class.

The present filtration does not permit lookahead equations involving an
unconsumed input or a future route choice. Such equations can make W depend
on the unrevealed permutation. They would require an enlarged exposure state
and new conditional probabilities. Likewise, one cannot select a source
support after inspecting the route without separately counting that choice.
The current proof fixes each support class before route exposure and sums
its alternatives outside the class bound.

## 6. Lazy outer coefficients require a separate exposure argument

The ordinary proof above conditions on the full outer setup. That convention
cannot simultaneously be used to call an outer coefficient fresh. The
[lazy-outer proposal](LAZY_OUTER_COEFFICIENTS.md) instead leaves later raw
Toeplitz diagonals unrevealed and introduces future mixed-input nonzero
requirements only when those forms become exposed.

An online version of that proposal would need the suffix bound to be uniform
also over the pending outer data: retained pivot forms, revealed diagonals,
and the unexposed parts of each column. It must establish that:

- the shortened source space is fixed independently of the unexposed route
  and raw diagonals;
- each column's source partition is performed before its scalar regions,
  without conditioning on future route positions or future source support;
- the stated fresh diagonal remains independent of the full exposed prefix
  until its particular group/component input is consumed;
- every lazy active-step bound is uniform over this enlarged prefix history;
- the final invertibility-conditioning factor is applied once, as specified
  by the raw experiment, rather than absorbed into a local route probability.

The four component regions share one Toeplitz matrix. They are not fresh
independent matrix samples. Deferred exposure must use the proved diagonal
order, not infer independence from a component's new physical position.

For the shared-full-column class, the lazy experiment satisfies these
conditions. The absent columns define V independently of the route and the
pending raw diagonals. Each pivot partition uses pre-mixer RS forms and the
current W before that column begins. On each fixed enumerated branch, it
therefore consults no future route label. Earlier components use only lower
Toeplitz diagonal indices, while different groups have independent matrices.
Thus the current component's fresh diagonal has not been used by an earlier
batch, even when the other group's input occurred first in the same region.

Enlarge the full prefix to include these retained pivots, all revealed
diagonals, and the pending labeled forms. Conditional on that prefix, both
the remaining route labels and the required fresh diagonal have the stated
distributions. The local lazy bound is uniform over this information, so the
induction in Section 4 applies. Retain the first b genuinely invertible
mixers and apply the final conditioning factor exactly once. Do not add
future mixed-input nonzero requirements before they are exposed.

An independent proof/code review found no discrepancy in this restricted
combination. The first screen in `theory/online_route.py` remains a failed
floating-point bound, not a directed numerical certificate. Exact tests now
cover 12 small-route laws/ordered linear recurrences, 30 local terminal and
boundary comparisons, and 396 directed entries against exact three-batch
nonlinear tables. They do not assume that nonlinear averaging commutes.
The next implementation gate is a tighter regional boundary continuation.
