# Finite small-group gates for the small-state GF32 proposal

The full-GF32 refresh proposed in [the cost review](STRUCTURAL_COST_REVIEW.md)
admits a small analytic distance gate. Conditional on any realization of
the retained outer mixers and paired route, the probability that a nonzero
message in one outer group has output weight at most 734003 is less than
2^-411, including all 8192 group choices. A separate two-group extension
below bounds their combined union by 2^-405. Probability is over the inner
diagonals only. These analytic bounds do not cover larger source sets.
The subsequent [paired-route transfer](SMALL_STATE_ROUTE_GATE.md) extends
coverage to one through six actual groups, below 2^-44 combined.

The earlier restricted-coefficient construction already has a
[one-group partial certificate](../char2_field_spin/CHEAP_LANE_RANK_INTERFACE.md#a-closed-bound-for-one-actual-source-group).
The new point is a short full-field argument with no state-mask operator,
outer-density substitution, or numerical optimization. It provides a
bounded proof gate for a cheaper implementation direction, not evidence
that the remaining whole-code proof is automatic.

## Construction held fixed

Let E=GF256, F=GF(2^32), q=2^32, K=2^20, and N=2^21. Use exactly the
[historical outer and route](../char2_field_spin/HD_SHIFTED_CLOSED_CERTIFICATE.md):
8192 groups of four RS[64,32] constituents, independent invertible
E-Toeplitz mixers on their four-symbol columns, and 256 regions with
alternating adjacent-pair matchings, uniform pair permutations, and fair
independent swaps. Fix any realization of these invertible mixers and route.

The inner has T=N/64=32768 batches and a 16-symbol state:

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H16 D1_t (Q_t + C X_t),
    Q_0 = 0.

A consists of four RS[16,4] evaluation maps; C=A^T and CA=0. H16 is the
fixed superregular E-valued interpolation map from nodes 0..15 to 16..31.
Every entry of every D1_t,D2_t is independently uniform in F minus zero.
The final state is discarded. In particular, the final refresh cannot
contribute a probability charge to the output-distance argument.

The retained route uses a striped local layout. Each physical pair feeds
two distinct RS[16,4] blocks, not adjacent coordinates of one block.
For 16 pairs, with probabilities p0, p1, p2 of zero, one, or two active
members, the two-block occupancy polynomial is
`(p0+p1*(x+y)/2+p2*x*y)^16`. Two independent halves supply the four
blocks. This is the layout in `hd_shifted_alias_certificate.odd_coefficients`
and the historical encoder's `Plan<Group>` route assignment. Sorting
the four occupancies requires summing all ordered configurations; it does
not make the two members of a pair independent.

The proof uses two facts about these actual maps. Every column of C is
nonzero. Also, for a set Z of requested output zeros,

    rank_F(A_Z) = sum_(j=1)^4 min(4, |Z_j|) >= |Z|/4,             (1)

where Z_j is the restriction to one 16-coordinate block. This is a rank
statement with all four state columns present. It does not assert that
arbitrary row/column restrictions of the coefficient-basis A are superregular.

## The live-state bound that survives

For every fixed nonzero v in F^16 and every nonzero target w,

    Pr[D2 H16 D1 v = w] <= (q-1)^(-16).                           (2)

To see this, let i be the input support size and j the target support size.
The 16-j requested zeros impose independent equations on the active input
coordinates until their rank reaches i. If 16-j>=i the event is impossible.
Otherwise, condition on i-(16-j) nonpivot nonzero inputs. The pivot values
are then determined. The zero constraints cost at most (q-1)^(-(16-j));
D2 reaches the prescribed nonzero target amplitudes with probability
(q-1)^(-j). Multiplying gives (2). Invertibility excludes a zero image.

Define gamma=(q/(q-1))^16. Given a live pre-refresh syndrome and any fixed
current input X, (2) and (1) imply

    Pr[(X+A Q)_Z=0] <= gamma q^(-rank A_Z).                       (3)

There are at most q^(16-rank A_Z) solutions for Q; counting a zero solution
too is harmless. This is a pointwise bound conditional on the past, not an
assertion that Q is exactly uniform.

Expose the refresh of batch t-1 immediately before emitting batch t.
The syndrome Q_(t-1)+CX_(t-1) is known before that refresh. Its being nonzero
is exactly the event Q_t!=0. Thus the live indicator is predictable before
the fresh randomness used in (2)--(3). For fixed proposed live and zero
masks, intersect their events sequentially. Do not condition and renormalize
on the entire future mask. The inequalities remain valid after discarding
other path restrictions.

## Outer dimension and activity geometry

Fix one actual active group. Let H be the exact set of its nonzero RS
columns, and write h=|H|. Invertible column mixers preserve whether all four
evaluations at a column vanish. The RS distance gives 33<=h<=64. Requiring
the other 64-h columns to vanish leaves a shortened source space of dimension

    d_h = 4(h-32).                                                (4)

Only the exact-H subclass is under consideration. Its number of F-projective
message lines is at most P_h=(q^(d_h)-1)/(q-1). The shortened superspace can
also contain smaller supports; those words are not assigned the geometry of H.

One group contributes at most one input coordinate in each region, and
region boundaries align with batches. Hence every batch has at most one
nonzero input. Since C has no zero column, a nonzero input activates a
quiet state. A live state can return to zero only in a nonempty input batch.

Let r count completed live excursions before the final batch. Let U be
the union of their batch intervals and the possible final excursion,
including activation batches. The final excursion is unrestricted. Write
c for the number of intervals, u for their total length, and a for the
number of batches with nonzero incoming state. Then

    c <= r+1,             a = u-c.                               (5)

Adjacent excursion intervals are allowed; a final activation may have
length one. Every nonzero input lies in U. Four consecutive regions form
one RS-column superregion, containing 512 batches. An interval of length
ell touches at most ell/512+2 such superregions. Consequently

    h <= u/512+2c,
    a >= 512(h-2c)-c.                                            (6)

No route probability or typical-spacing assumption is used here.

## Split at a bounded number of returns

For h=33..64 set

    r_h = floor((h-32)/4)+2,       a_h = 512(h-2r_h)-r_h.          (7)

These integers obey r_h<=10, 16r_h-d_h>=20, and a_h>=14846.
For example, writing h-32=4k+j with 0<=j<=3 gives
a_h=14334+1023k+512j. The minimum is h=33, k=0, j=1.

### Many returns

A live return at time t requires Q_t=-CX_t, a prescribed nonzero vector.
It costs at most (q-1)^(-16) by (2). Chosen return times use distinct fresh
incoming refreshes. For a fixed message, the probability of at least r_h
returns is at most binom(T-1,r_h)(q-1)^(-16r_h). The binomial overcounts
the eligible nonterminal times, which is safe. Projective union gives

    U_hi <= 8192 sum_(h=33)^64 binom(64,h) P_h
                         binom(T-1,r_h) (q-1)^(-16r_h).           (8)

This branch discards the output-weight restriction.

For a compact exact bound, Bernoulli's inequality gives

    (q/(q-1))^m <= q/(q-m) < 2          when 0<=2m<q.             (9)

Here m<=16T=524288. Using P_h<=q^d_h, binom(T-1,r_h)<=T^10=2^150,
q^(d_h-16r_h)<=q^-20=2^-640, and sum_h binom(64,h)<=2^64 yields

    U_hi < 2^(13+64+150-640+1) = 2^-412.                         (10)

### Few returns

If r<r_h, (5)--(7) give a>=a_h. Fix an incoming-live index set of size a.
Output weight at most e=734003 forces at least z_a=64a-e zeros among its
64a output coordinates. For any selected z_a coordinates, sequential use
of (3) bounds the joint probability by

    gamma^a q^(-ceil(z_a/4)).

The fixed live-index event remains an intersection inside the expectation,
not a conditioning on future randomness. The exact conservative union is

    U_lo <= 8192 sum_(h=33)^64 binom(64,h) P_h
               sum_(a=a_h)^(T-1) binom(T-1,a) binom(64a,64a-e)
                                  gamma^a q^(-ceil((64a-e)/4)).  (11)

All z_a are positive. The zero-mask binomial is at most 2^(64a), and
q^(-ceil(z_a/4))<=2^(-8z_a). Therefore their product is at most
2^(8e-448a), decreasing in a. At a=14846 its exponent is -778984.
Equation (9) gives gamma^a<2. Bounding projective lines by q^128=2^4096
and all live sets by 2^T now yields

    U_lo < 2^(-778984+4096+32768+64+13+1) = 2^-742042.             (12)

Combining (10) and (12) gives U_hi+U_lo<2^-411. This includes every
nonzero one-group message and all group choices, with a free terminal.
The bound is uniform over the fixed outer/route setup, so averaging that
setup preserves it. Four SIMD field lanes do not require four independent
certificates: any nonzero component of a record-valued message gives a
scalar witness with contained output support.

## A separate two-group theorem

Fix two actual active groups, and let H be the union of their nonzero
RS-column sets. Both groups vanish outside H, so their combined shortened
dimension is at most d_h=8(h-32). Count only the exact union-H subclass;
its projective line count is at most (q^d_h-1)/(q-1). There are at most
binom(8192,2)<2^26 choices of groups and at most 2^64 column sets.

Each batch contains at most two possible input coordinates. Any two columns
of C are independent: inside one block this is the RS[16,4] evaluation
rank, and different blocks occupy disjoint state coordinates. Thus quiet
nonzero input still activates the state. All exposure, return, emission,
interval, and free-terminal arguments above remain valid.

Set r_h=floor((h-32)/2)+2. The exact expressions (8) and (11) apply with
the new d_h and group count. Now r_h<=18, 16r_h-d_h>=24, and

    a_h = 512(h-2r_h)-r_h >= 14318.

For h-32=2k+j, 0<=j<=1, this is a_h=14334+512j-k; its minimum occurs
at h=64. The same conservative integer bounds give

    U_hi,2 < 2^(26+64+15*18+1-32*24) = 2^-407,
    U_lo,2 < 2^(8*734003-448*14318+8192+32768+64+26+1)
           = 2^-501389.

The one-group and two-group unions together are therefore less than
2^-405. Neither statement uses typical placement or averages the route.

The [integer checker](theory/check_small_state_one_group.py) checks all
64 combinations of group count in {1,2} and h=33..64, their finite union
inequalities, and bounded toy interval geometry. It does not replace the
field/rank and exposure-order proofs. The root's serial replay passed all
64 column classes and 8144 tiny interval fixtures under a 60-second timeout
and a 256-MiB memory cap. It reproduced the 411-bit one-group, 406-bit
two-group, and 405-bit combined bounds using integer arithmetic only.
Run `python -B theory/check_small_state_one_group.py` from this workstream.
The source writes no files. Independent analytic review covered the
one-group proof and the two-group extension's injectivity and geometry;
the root additionally reviewed the full formulas and checker.

## The next obstruction is not more singleton margin

This proof does not use random routing, which is a benefit for this sparse
gate but cannot hold uniformly for all source sizes. It also avoids the
old ambient scalar-shadow operator's artificial activations by counting
actual F messages on an actual source support.

The same coarse strategy already has an explicit limitation at three
groups. If their common occupied-column union has h=64, the shortened
source dimension bound is 384. The analogous return cutoff is 26, and
(6) then permits only

    a >= 512(64-2*26)-26 = 6118.

Since 64*6118<734003, the low-return branch need not force any output zeros
on live batches. Its present inequality therefore cannot close that class.
This is an obstruction to the bound, not an exhibited low-weight word.
For at most four groups, C is injective on every possible batch input
support: each emission block contains at most four selected RS evaluation
columns. For five or more groups, a second issue appears: five positions
in one emission block can contain a nonzero input in ker C, so quiet
nonzero input need not activate the state. The retained all-quiet receipt
does not automatically bound trajectories that are quiet only part of the time.

A useful next lemma must retain the paired-route probability of long gaps
and concentrated reset coverage, or retain source-rank constraints consumed
by resets. Later gates must also charge constraints that keep inputs quiet.
These arguments must use the old paired route and grouped E-valued outer rather
than silently import the new packet-four route or a uniform-F density.
Closing that bounded multigroup gate is more informative than increasing
the already excessive one-group margin or launching a large mask operator.
