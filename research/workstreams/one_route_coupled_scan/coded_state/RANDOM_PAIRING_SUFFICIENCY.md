# Is one random pair layer enough?

2026-10-26. Design analysis following the
[conditional audit](RANDOM_PAIRING_CONDITIONAL_AUDIT.md).
The implemented construction, existing certificates, and paper are unchanged.

## Recommendation

Reconsider one pair layer per update as the final inner. This is a
proof-and-design recommendation, not a demonstrated distance obstruction.
Random pairing could still have the desired distance. However, a stronger
joint count remains substantially inconclusive, whereas a four-coordinate
comparison passes the same difficult family and its selected neighbors.

Keep the RS outer, regional route, 512-coordinate state, and recurrence.
Use a randomized quartet mixer as the next stronger control. Establish
broader coverage and a cheap realization before adopting it. Two pair
layers do not automatically implement the quartet distribution analyzed
here. There are no native timings for either this comparison or the
current random-pairing construction.

## What seems weak, and what does not

The accepted pair map already has almost the best possible point bound
on the image of a fixed nonzero pair: approximately q^-2. More random
coefficient bits inside that pair do not supply substantially more than
two random coordinates. The question is how many coordinates a state
reaches before the next external input can cancel or steer it.

With zero local input, a superregular pair map takes one active coordinate
to two. Two active coordinates can produce only one active coordinate.
Re-pairing prevents a permanent fixed partition, but it does not make one
update dense. Inputs arrive between pair layers, so an empty-input
spreading argument cannot simply be applied across many updates.

For a superregular quartet map, an incoming support of size a>0 gives
at least 5-a nonzero next-state coordinates on zero local input. The same
statement holds for the output if A+delta I is also superregular. In
particular, two active coordinates now force at least three outputs,
rather than permitting one. This addresses the partially active groups
directly, before another input arrives.

These are local statements. The outer prevents arbitrary independent
choices of cancelling inputs, and that restriction must remain in the
global count. Neither local example is a low-weight codeword.

## Controlled comparison

Retain q=2^32, K=2^20, N=2^21, 8192 RS[256,128] rows, 256 uniform regional
permutations, and the recurrence

    Q_t=A_t Q_(t-1)+X_t,   Y_t=Q_t+delta Q_(t-1),   Q_0=0.

The current model samples a uniform perfect matching at each update and
uses the accepted SL2 pair matrices. For the comparison, partition the
512 coordinates uniformly into groups of four at each update. On each
quartet, sample a uniform 4-by-4 matrix over GF(q), rejecting unless both
A and A+delta I are superregular. The new partitions and matrices are
sampled once at setup and fixed for all messages. They are independent
of the route. This is a local state operation, not another global route.

The comparison has a justified point bound. Before rejection, A v is
uniform over GF(q)^4 for each fixed v!=0. A j-by-j minor is a nonzero
polynomial of degree j in the raw matrix entries. The union bound for
all minors of A and A+delta I costs at most

    (2 sum_(j=1)^4 j binom(4,j)^2)/q = 280/q.

There are N/4 matrices. Their total density inflation is at most

    1/(1-(N/4)*280/q) = 1024/989 < 2.

Here the denominator also follows directly from the union bound on raw
matrix rejection anywhere in setup. Consequently the same global factor
two used for the current pair model is valid for this quartet comparison.
This does not justify substituting an arbitrary fast structured matrix.

## Additional constraints used in the new count

Fix a candidate incoming support U, next support V, and a set Z of
selected zero outputs inside U intersect V, within one mixing group.
Let T be the group coordinates whose outer rows are inactive, so their
input coordinates are zero. Write a=|U|.

If a=0, a necessary condition is T intersect V empty: an empty state
can only activate where its input is nonzero.

If a>0, superregularity gives two further necessary conditions:

    |T minus V| < a,
    |T intersect ((group minus (U union V)) union Z)| < a.   (1)

The first counts forced zeros of A Q. The second counts forced zeros of
(A+delta I)Q. Any a such zeros would give an invertible a-by-a submatrix
annihilating the nonzero vector on U, which is impossible.

For pairs, these conditions include the useful special case: a one-active
pair with zero input on a coordinate cannot have either next state or
output zero there. The previous coarse count omitted some of these
conditions. Their effect is measurable, but too small to close its gap.

## Joint count and probability space

The tested family has 256 active state coordinates at every update.
At each of the 4095 transitions after activation, exactly b coordinates
leave and b enter. Fix r inactive outer rows and union over their choices.
For W=floor(0.46 N)=964689, put

    A=4096*256,
    D=4095*b,
    z=A+D-W.

The displayed tests have 0<=z<=4095*256-D. Each candidate history needs
at least z zero outputs in consecutive-support intersections.

We use the conditional pin-and-shortening argument from the audit, with
live groups in place of live pairs. It gives a field exponent

    E=A-K+128r-C_T-z,

where C_T counts covered coordinates belonging to inactive rows. This
identity does not require every outside crossing to belong to an active
row. That extra pair-specific condition cannot be imposed on quartets.
Conditions (1) supply the justified replacement.

For a fixed geometry and history, let C count zero selections satisfying
(1). The conditional bad-message probability is bounded by

    min(1,2 C q^E).

The count includes all messages on the exact active rows. If an active
row has fewer than 129 covered positions, its contribution is zero.
Dropping that test only increases the bound. For 0<theta<=1, use
min(1,x)<=x^theta after counting messages, before averaging geometry.

The following positive generating functions retain local compatibility
jointly with matching, support changes, and selected zeros. They replace
the audit's separate history factor and worst-case geometry cutoff.

Let m=32 theta, and let ell,u,v be positive markers for selected zeros,
next-state weight, and deaths. For fixed U,V,T within a group, define

    P_(U,V,T)(ell)=sum_(Z allowed by (1)) ell^|Z|.

It is zero for impossible local support choices. Zero selections across
groups satisfy

    C <= ell^-z product_groups P_(U,V,T)(ell).

For a fixed set of r inactive rows, represent their regional positions by
independent Bernoulli labels of any probability p, conditioned to contain
exactly r labels in each region. This is the actual uniform regional
distribution after conditioning. The conditioning probability per region is

    pi=binom(8192,r) p^r (1-p)^(8192-r).

Before conditioning, define L_a as the sum over V,T for a fixed incoming
support U of size a in a group of size g, with summand

    p^|T| (1-p)^(g-|T|) 2^(-m*c)
      u^|V| v^|U minus V| P_(U,V,T)(ell)^theta.

Here c=|T| for a live group, and c=|T intersect V| for an empty group.
The latter is zero on allowed choices. A uniform partition gives the
following upper bound on one transition with next weight 256 and b deaths:

    K_g = [h^256](sum_(a=0)^g binom(g,a) L_a h^a)^(512/g)
          / (binom(512,256) u^256 v^b).                    (2)

The coefficient in (2) averages the actual incoming group occupancies.
It is not an independent-partner approximation. Positive u and v bound
the required next-weight and death coefficients. Symmetry makes this
transition bound independent of the incoming support's labels, so it
composes over all histories in the specified family.

Initial activation has binom(512,256) choices and requires its 256 input
positions to be on active rows. The complete family bound is

    binom(8192,r) binom(512,256) (1-p)^256 pi^-256
      * 2^theta * 2^[m(A-K+128r-z)] * ell^(-theta*z)
      * K_g^4095.                                        (3)

This includes all row sets, supports, permitted selected zeros, and
messages in this particular family. It is not just a fixed-message test.
Setup is not refreshed per message. The Bernoulli experiment is only a
proof device, with its full conditioning cost present in (3).

## Results at the original difficult point

Use b=32 and r=1974. Logarithms below are base two. A positive upper bound
is inconclusive; it is not a lower bound on the failure probability.

| Local mixer and count | Log2 upper bound | Status |
|---|---:|---|
| Pairs, coarse local compatibility | +286860 | FLOAT proposal |
| Pairs, both constraints (1) | +282899 | FLOAT proposal |
| Quartets, both constraints (1) | -99576 | FLOAT proposal |
| Quartets, rational positive-weight replay | -44413.390625 | Directed restricted bound |

The pair improvement is about 3961 bits, against a remaining positive
bound of about 283000 bits. The quartet change is much larger. We did not
prove that the optimized pair bound is globally minimal; two numerical
starts agreed closely. More importantly, even failure of the best such
upper bound would not prove that the pair code lacks distance.

The directed replay uses m=6, theta=3/16, p=8/11, ell=5/16, u=5/18,
v=1/24. It upper-bounds the incoming coefficient in (2) with h=12.
Each P^(3/16) is rounded upward to a 40-bit dyadic using exact integer
sixteenth-power comparisons. Rational logarithm intervals come from
sixteenth powers and integer bit lengths. The resulting log2 upper bound
is exactly -2842457/64. No floating-point optimizer result is trusted in
this replay. The global accepted-map density factor is already included.

This proves a restricted failure bound for the quartet comparison at
46% distance: all full-length half-full histories with exactly 32 swaps
per transition and exactly 1974 inactive outer rows. It is not a
whole-code certificate, a 49% result, or a statement about the current
pair construction.

## Neighboring profiles and occasional stronger mixing

A bounded FLOAT screen used b in {16,32,64} and r in {1024,1974,2560,4096}.
All twelve quartet bounds are negative; five pair bounds remain positive.
The worst selected quartet result is the original -99576 case. Other
activity levels, turnover patterns, and target-row counts were not screened.

For the original b=32,r=1974 case, we also replaced a fixed fraction of
updates by the quartet comparison, retaining pairs at the others. The
same proof uses the product of the respective transition bounds.

| Fraction of post-activation quartet updates | FLOAT log2 upper bound |
|---|---:|
| 25%, rounded to 1024 updates | +192759 |
| 50%, rounded to 2048 updates | +98523 |
| 75%, rounded to 3071 updates | +960 |
| 100% | -99576 |

The markers are shared between update types, so this need not be the best
mixed-schedule bound. Nor does this family cover adversarially concentrated
turnover around weaker updates. These results do not support a claim that
rare quartet updates alone fix the proof.

## Decision and next step

No actual distance obstruction was found for random pairs. The earlier
bounded rank search also found none, but neither search proves absence.
Thus the mathematical answer to whether pairs suffice remains open.

The research decision is clearer: do not invest in completing the old
flow enumeration now. The joint count tested here still struggles with
pairs, and a modestly larger local mixing group changes the same test
from failure to a directed restricted success. This is a better-founded
next design step than adding coefficient randomness inside a pair.

Next examine quartet coverage across varying activity and turnover,
alongside a cheap quartet realization. Preserve the uniform accepted
quartet law as the proof control. A factorized circuit, two pair layers,
or subfield multiplier needs its own image law and compatibility proof.
Do not inherit (3) or its point probability merely because the circuit
touches four coordinates. Dense quartet arithmetic may cost more; no
end-to-end performance claim is made before implementing and measuring it.

## Validation and storage

`diagnose_pairing_sufficiency.py` retains the FLOAT screens and the directed
quartet replay. Its 10740 checks include local constraints on 1302 accepted
GF8 pair maps with selected state/input vectors, coefficient-power identities,
empty-state constraints, and exact small random-partition averages. These
checks support the implementation; the general arguments above establish
the stated conditional inequalities.

One initial FLOAT run rejected an underflowing coefficient evaluation.
Positive polynomial recentering fixed the evaluation; it preserves the
coefficient exactly as an algebraic identity. No failed run was counted as
a small bound. The selected screens subsequently terminated in under five
seconds each. All runs were serial, file-free, and guarded by a 60-second
watchdog, fixed small problem sizes, and more than 20 GiB free space.
No raw output, native benchmark, production edit, commit, or promotion.
