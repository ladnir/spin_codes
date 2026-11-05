# The retained paired-route gate transfers to full GF32 diagonals

The three-/four-group obstruction in the
[analytic sparse gate](SMALL_STATE_PROOF_GATE.md) is a limitation of its
deterministic interval bound. The earlier paired-route calculation already
retains enough placement and source-rank information. Its positive operator
also bounds the proposed full-GF32 refresh, without a new spacing estimate.

This note transfers that operator, not the entire restricted-coefficient
proof. In particular, the later scalar-shadow and minimal-polynomial
arguments are not covered by this transfer.

## Construction and event

Use exactly the construction in [the sparse gate](SMALL_STATE_PROOF_GATE.md):
K=2^20, N=2^21, 8192 groups of four RS[64,32] rows, independent invertible
GF256-Toeplitz column mixers, and the alternating paired regional route.
The inner emits 64 symbols from a 16-symbol state:

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H16 D1_t (Q_t + C X_t).

Here A, C, and the superregular H16 are the retained GF256-valued maps.
Both diagonals have independent entries uniform in F minus zero, where
F=GF(2^32). The initial state is zero and the terminal state is free.
Setup is sampled once for every message.

For an integer g, let Bad_g mean that some message supported on exactly
g actual outer groups has output weight at most 734003. This is an
existence event over all F-valued messages chosen after setup. It is not a
fixed-message average or a count restricted to GF256-valued messages.

## What the direct gate needs

The relevant retained argument is
[the direct two-through-four-group gate](../char2_field_spin/CHEAP_LANE_RANK_INTERFACE.md#extending-the-direct-rank-count-to-two-through-four-groups).
It fixes an actual source set and uses one common F-linear quotient of its
message space. Its rank price is tau=65280, independent of the field size.
The potential at residual dimension d is

    P_tau(d) = (tau^d-1)/(tau-1).

A surviving nonzero witness has potential at least one. An independent
pivot lowers the dimension and contracts this potential by at least tau.
A dependent pivot survives for at most one value of its fresh D2 entry.
Under the new distribution its probability is at most 1/(2^32-1), which
is smaller than 1/tau. Thus the same pivot contraction remains valid.
There is no multiplication by the number of F-valued messages.

Fix any D1 realization. Multiplication by nonzero diagonal entries only
rescales columns of H16. Consequently the MDS support constraints,
extinction constraints, output-weight floors, and pivot-rank lower bounds
used by the direct gate are unchanged. The sole averaged local factor
that needs replacement is the activation test.

## Full-field activation is no worse

Let q be a field size, let c be the dimension of a fixed possible input
image, and let m be the number of requested zero coordinates after H16 D1.
Union over projective nonzero inputs in that image. For each input, the
superregular diagonal-pivot argument costs at most (q-1)^(-m). Hence the
activation probability is bounded by

    f_q(c,m) = min(1, P_q(c) (q-1)^(-m)),
    P_q(c) = 1+q+...+q^(c-1).                                  (1)

The old bound used q=256 because its D1 and input matrix were GF256-valued.
The new bound uses q=2^32 directly; it does not project an F-valued input
to a subfield witness.

For every c>=1 and m>=0, f_q(c,m) is nonincreasing in q>1. If m<=c-1,
the expression inside the minimum is at least one. Otherwise every term
q^i/(q-1)^m, for 0<=i<c, decreases with q, since

    d/dq log(q^i/(q-1)^m) = i/q-m/(q-1) < 0.

It follows that f_(2^32)(c,m)<=f_256(c,m). Combining the activation event
with the pointwise rank-potential cap before averaging D1 is essential.
This is not multiplication of two unrelated marginal probability bounds.

## The source ranks and exact route remain the same

For g<=4, every allowed input submatrix C_I has full column rank.
A quiet-to-quiet batch therefore forces its allowed inputs to zero.
If q_j of the four regions in outer column j are wholly quiet, their
constraints have rank g q_j on that column's evaluations. Any 32 RS
columns can independently realize all their evaluations. The invertible
column mixers preserve this fact after extending scalars to F.
Thus source rank is at least

    g * sum_largest32(q_0,...,q_63).                             (2)

No scalar shadow is needed for (2). The argument remains valid for
F-valued equations: constraints on disjoint evaluations within a column
have additive rank, and any 32 complete columns are jointly surjective.

The retained operator averages the actual paired route. In each region
it extracts the coefficient for the specified singleton/full-pair counts,
with the exact without-replacement normalizer. It retains one fixed pair
profile in each parity throughout the 64 columns, then maximizes over
those profiles and includes all source-set choices.

The best-32-column identity partitions paths into four threshold classes.
Each class retains the strict source-rank cutoff, which excludes the zero
message. All transfer matrices, polynomial coefficients, cutoff weights,
and final factors are nonnegative. Replacing (1) by a smaller entry and
retaining all other entries therefore cannot increase the complete bound.
The first batch has zero incoming state; no charge uses the final refresh.

## Inherited finite bounds

The retained [directed receipts](../char2_field_spin/CHEAP_LANE_MIXING_GOAL.md#two-through-four-active-groups-also-close-at-35)
use tau=65280 and distance 7/20. They already include all source supports
and matching profiles. By the preceding domination, the same upper bounds
apply to the new full-GF32 diagonals:

| Actual active groups g | Inherited bound on Pr[Bad_g] |
|---|---:|
| 2 | <2^-450 |
| 3 | <2^-250 |
| 4 | <2^-100 |

The fixed verifier is
[cheap_lane_few_group_gate.py](../char2_field_spin/theory/cheap_lane_few_group_gate.py).
The existing ledger records independent 192-/256-bit directed replays and
their rational markers. No new numerical run or performance measurement
was performed for this transfer. The analytic domination, rather than a
new calculation, connects those receipts to the changed distribution.

Together with the new one-/two-group analytic bound, groups 1 through 4
have combined failure probability below 2^-99. This remains a partial
minimum-distance result, not a whole-code certificate.

## The same comparison covers the direct five-/six-group gates

The [five-group gate](../char2_field_spin/CHEAP_LANE_RANK_INTERFACE.md#five-groups-retain-rank-even-when-quiet-inputs-cancel)
does not assume quiet inputs vanish. It retains deterministic input rank
max(0,c-h-j), where c is the actual feedback-image rank and h,j are state
support sizes. This rank lower bound holds for every nonzero D1 diagonal.
The only averaged activation factor is again (1). Within-column rank
additivity and the best-32-column cutoff still hold over F. Its retained
bound therefore gives Pr[Bad_5]<2^-50 for the new diagonals.

The [six-group gate](../char2_field_spin/CHEAP_LANE_RANK_INTERFACE.md#six-groups-preserve-activation-ranks-and-the-source-cutoff)
tags an activation's input rank d jointly with its probability. Its cap,
for 0<=d<c, has the form

    g_q(c,m,d) = min(1, P_q(c)(q-1)^(-m)/P_q(c-d)).              (3)

This cap is also nonincreasing in q. If m<=d it is one. Otherwise put
z=1/q and rewrite its unclipped part as

    z^(m-d) (1-z^c)/(1-z^(c-d)) (1-z)^(-m).

Each factor is nondecreasing in z. For the middle ratio, write a=c-d;
its derivative has the sign of a-(a+d)z^d+d z^(a+d), which is nonnegative
on 0<z<1. The projected live-transition cap has the same form.
Eliminating incoming columns uses their span, which is unchanged by D1's
nonzero column scalings. The selected projected subspace is fixed before
the remaining D1 entries are exposed. Thus the retained joint rank and
probability accounting also survives the change of field.

The existing five-/six-group receipts are implemented by
[the rank-debit verifier](../char2_field_spin/theory/cheap_lane_rank_debit_gate.py)
and [the joint-rank verifier](../char2_field_spin/theory/cheap_lane_joint_activation_gate.py).
Their inherited bounds are 50 and 45 bits, respectively. Combined with
the classes above, all messages on one through six actual groups have
failure probability below 2^-44. Independent analytic review covered
the direct transfer and both monotonicity arguments.

## Remaining work

This closes the proposed bounded three-/four-group gate by reusing the
actual route and consumed source ranks. It does not establish a generic
stochastic domination between whole code ensembles. In particular, later
proofs that expand an E-valued kernel into auxiliary E messages can change
when D1 becomes F-valued. Filters using theta's minimal polynomial also
need a different argument.

The next useful task is a large-source interface for the full-field DHD
refresh, retaining deterministic source ranks without that E-kernel
expansion. Matching more of the old small-support margins is not the
priority. No new kernel, runtime claim, or production promotion follows
from this proof transfer.
