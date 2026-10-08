# Preserved rank-based proof work for one regional route

The complete one-route result now uses GF(16) coefficients and four
state symbols; its proof is in [GF16_PROOF.md](GF16_PROOF.md). This
note preserves the distinct GF(256), two-state rank-based attempt.
No change to its construction or proof status is implied by the new
GF(16) certificate.

## Objective and fixed model

The GF(256) construction has a repeated GRS[16,8] outer, sixteen
independently permuted column regions, independent nonzero input
diagonals, and pre-output uniform SL(2,256) refresh. Its maps have
A_j=(1,j), C=A^T. The full fixed-word first moment is obstructed by
the rare large bad subspaces in [SUBFIELD_PROOF.md](SUBFIELD_PROOF.md).
The rank-based attempt instead counts prescribed zero equations.

Fix a source support and an output support path before exposing setup.
The source is shortened to a vector space V of dimension D. A branch
remains alive only if its current quotient admits a nonzero witness
satisfying all prescribed nonzero forms. Exposing a new independent
zero equation multiplies a rank potential by tau^-1. An alive branch
has at most D-1 independent equations.

The fresh-diagonal pivot lemma costs

    g=1/(q-1)+(1-1/(q-1))/tau

per independent coefficient pivot. It remains valid for dependent
input forms because every input form is required nonzero on an alive
branch. The state rank must remain explicit: collapsing ranks one
and two repeatedly charges a rank loss that can occur only once.

## A sharper active-state return bound

Partition lanes into canceled active I, uncanceled active P, inactive
zero U, and inactive nonzero lanes. Write R=M Q and X=diag(d)f. On
a return to zero state, canceled inputs satisfy X_I=A_I R, and

    F_I R=C_P X_P,       A_U R=0,
    F_I=I+sum_(j in I) C_j A_j.

Expose the fresh M first. Let K stack A_U and the product of F_I
with a left annihilator of the span of C_P. The required nonzero R
must satisfy K R=0. Rank(K)=2 kills the branch. Rank(K)=1 costs
1/(q+1) for an incoming rank-one state, or one independent equation
for an incoming rank-two state. Rank(K)=0 imposes no zero equation;
the rank-one state must still avoid every required-nonzero lane.

Only after these matrix-only conditions are imposed do we expose the
fresh diagonal. Its remaining affine equations have coefficient rank

    |I|+min(2,|P|).

The diagonal lemma charges g to this power. This exposure order does
not charge the same randomness or rank loss twice. It improves a
singleton return from a coarse whole-state bound to a projective-line
check followed by one diagonal pivot.

The two-state maps give particularly small geometry. If p=|I| modulo
two and x is the XOR of its lane labels, then

    F_I=[[1+p,x],[x,1+x^2]].

This matrix is singular exactly when p=1 and x=0, in which case it
is diag(0,1). Its image then excludes every nonzero multiple of
C_j=(1,j), so a return with exactly one uncanceled active lane is
impossible. Otherwise that lane fixes the projective direction of R.
Directions that kill a required-nonzero canceled input or the required
nonzero uncanceled output are removed. A subset-count dynamic program
over cardinality and XOR aggregates these conditions exactly.

## Composition and numerical status

`theory/one_route_rank.py` composes the local bounds with regional
binomial normalization. Its linear `sum` mode upper-bounds an unknown
nonzero next-rank continuation by the sum of rank-one and rank-two
continuations. This is valid but can be loose. The source envelope
uses tau^(w-8) per outer row; a separate four-row calculation retains
every outer support size and removes all impossible sizes below nine.

At K=65536 and 10% distance, floating searches with exact source
support sizes give the following rank-moment screens:

| Active outer rows | Screened margin |
|---:|---:|
| 1 | 42.04 bits |
| 2 | 8.17 bits |
| 3 | -37.31 bits |
| 4 | -77.86 bits |

The sharper return lemma improves the former coarse four-row bound
by over one hundred bits. However, retaining exact source support
sizes changes the four-row result by less than one bit. These screens
do not establish a complete distance certificate.

The optional `max` mode is only a diagnostic after source aggregation.
A maximum is sound for a fixed future source type. Moving it outside
a sum over source types is not justified: those types can change the
incoming quotient and its actual rank. For example, the continuation
vectors (1,0) and (0,1) can select different ranks and sum to two,
although the maximum coordinate of their sum is only one. No finite
certificate uses that mode.

A useful next refinement would retain the shortened-source dimension
and stop a path when its zero equations reach full rank. The present
scalar tau moment permits such already-dead paths in its envelope.
For a small number of active rows, exact rank truncation may be more
effective than another scalar tilt or outer-support count.

## Reproduction

```sh
python -B research/workstreams/char2_field_spin/theory/one_route_rank.py --rows 1 2 3 4 --steps 12 --exact-source
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_one_route_rank.py -v
```

Five tests check the GF4 return-mask aggregation, singleton and empty-
input cases, regional normalization, numerical continuation comparisons,
and invalid geometry. Independent GF4 checks in
`lane_mixing/test_rank_lemmas_audit.py` verify the exposure-order lemmas,
including dependent input forms and their nonzero requirements.
