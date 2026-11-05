# Every message on one active outer row

For the unchanged construction, the probability that any nonzero message
supported on exactly one outer RS row produces output weight at most
W=943718 is less than 2^-82. This is a whole input-support class, with
every input weight, route, state history, death count and terminal state
included. It does not close the full distance goal.

An **active input row** contains a nonzero outer codeword. It is different
from a footprint-eligible row, which has at least129 covered positions.
The theorem allows the resulting footprint to make arbitrarily many rows
eligible. In particular, it covers the entire one-eligible-row remainder,
because any nonzero message there has exactly one active input row.

## Fixed-message counts and three exceptional events

Let the unique nonzero RS row have weight 128+d, where 1<=d<=128.
For any prescribed support of that size, MDS shortening leaves dimension
d. There are at most8192 choices of row and at most2^256 choices of its
regional support. Thus the number of fixed messages in this d-class is
at most 2^269 q^d. The accepted-setup density gamma<2 makes the common
union prefactor less than 2^(270+32d).

Fix one of these messages and a route before revealing the random inner
maps and partitions. Let E count state zeros at live coordinates, let r
count transitions from a nonempty state to the empty state, and let B
count the bad partial-state partitions defined below. Set

    e_d = ceil((32d+362)/14),
    r_d = ceil((32d+362)/116),
    b_d = max(48,ceil((32d+362)/12)).

Each of the three events E>=e_d, r>=r_d and B>=b_d has probability less
than 2^-92 after the fixed-message union for this d-class.

For E, the [adaptive zero argument](QUARTET_ADAPTIVE_ZERO_BUDGET.md)
gives binom(N,e)q^(-e). Since e_d>=24,
binom(N,e_d)<=(3N/e_d)^e_d<=2^(18e_d). This leaves a cost of14 bits per
state zero after the subset count.

For an extinction, a nonempty incoming state makes at least one quartet
live. Its four next-state coordinates must all be zero, with conditional
probability at most q^-4. Any nonzero input outside the live set makes
extinction impossible. For fixed prescribed extinction updates, tower
conditioning therefore gives q^(-4r). The update subset count is at most
binom(4096,r)<=2^(12r), leaving116 bits per extinction. All later input
restarts remain in this count.

The growth lemma below gives probability at most2^-20 for a bad partition
at any prescribed partial-state update, conditional on the entire past.
Since b_d>=48,
binom(4096,b_d)<=(3*4096/b_d)^b_d<=2^(8b_d). This leaves12 bits per bad
partition. The definitions of e_d,r_d,b_d make each total exponent at
most-92. There is no state-support union.

## A finite uniform growth lemma

For state population u define the integer level

    ell(0)=0,
    ell(u)=1+floor(log2 u),                  1<=u<=256,
    ell(u)=17-ceil(log2(512-u)),             257<=u<=511,
    ell(512)=18.

It is nondecreasing and1-Lipschitz on integer populations. For a partial
incoming support, call its new uniform quartet partition bad when the
covered population g=|Gamma| does not increase the level. Since g>=u,
this means ell(g)=ell(u).

Only thirteen exact integer comparisons are needed to prove the uniform
bound Pr[bad partition | incoming support]<=2^-20. For
u in {4,8,16,32,64,128}, let k=u/2-1. A failure is contained in the event
that the old support lies in some k quartets, giving the union bound

    binom(128,k) binom(4k,u) / binom(512,u).

For a hole population v in {4,8,16,32,64,128,256}, let k=floor(v/8)+1.
A failure needs k quartets consisting entirely of holes, giving

    binom(128,k) binom(512-4k,v-4k) / binom(512,v).

All thirteen ratios are at most2^-20 by exact integer comparison.
Within a level interval, a nested-prefix coupling of uniform supports
or holes makes the listed endpoint the worst case. The remaining
small-population cases have deterministic level growth.

Let P count partial incoming states over all updates. If a step has e_t
state zeros at live coordinates, its new population is at least g-e_t.
The level can lose at most e_t by the Lipschitz property. Every good
partial-state partition gains at least one level before these deletions;
all other partitions have nonnegative gain before deletions. Telescoping
from level0 to a final level at most18 proves

    P <= E+B+18.                                      (1)

This includes every extinction and subsequent restart.

## The output and regional occupation contradiction

The already certified fixed-one-row joint-zero budget bounds the event
K_bad>=254 by less than2^-92. Its proof counts messages on one prescribed
input row; it does not require that row to be the only footprint-eligible
row. Work on its complement, K_bad<=253.

Each fully active incoming state has512 live coordinates. At most K_bad
of their outputs can be zero over the whole computation. An output of
weight at most W therefore permits at most

    floor((W+253)/512)=1843

fully active incoming updates. If A counts nonempty actual states, (1)
and the one free terminal state give

    A <=1843+P+1 <=1862+E+B.                           (2)

The binary stream recording nonempty states has r closed runs and at
most one terminal run. Its covered stream is the union of each state
indicator and its predecessor, so it has exactly A+r occupied updates.
At most2r+1 regions can be partially occupied; every other occupied
region consumes16 occupied updates. Hence the number Q of occupied
regions satisfies

    Q <=floor((A+r)/16)+2r+1.

The nonzero RS input row occupies128+d distinct regions. Each such input
requires either the incoming or the new state to be nonempty, so
Q>=128+d. Combining this fact with (2) gives the necessary inequality

    E+B+33r >=170+16d.                                (3)

This uses a conservative endpoint allowance; no terminal reset is added.
The complements of the three exceptional events instead give

    E+B+33r <=(e_d-1)+(b_d-1)+33(r_d-1)
             <170+16d.                               (4)

All128 integer comparisons in (4) pass. For d>=7 they also follow from
the upper bound (32d+362)*535/1218: subtracting this from170+16d leaves
(2368d+13390)/1218>0. The six smaller d values are direct comparisons.
Thus a bad one-row message must belong to an exceptional family.

## Probability budget and verification

There are3*128 fixed-weight exceptional families and one joint-zero
family. Each contributes less than2^-92. Their total is therefore less
than385*2^-92<2^-82. Neither a further route union nor a state-history
union is required.

[quartet_one_active_row_certificate.py](quartet_one_active_row_certificate.py)
rechecks the thirteen growth ratios, the128 budget contradictions, and
the previously verified joint-zero threshold. It also checks the level
and region-occupation identities on independent small binary examples.
The complete serial replay passed512 level checks,33023 endpoint-tail
coverage checks,1616 binary occupation checks, all thirteen growth
comparisons, and all128 budget comparisons. It confirmed the joint-zero
threshold254 and a minimum integer contradiction margin of12. The
smallest certified inverse growth-ratio exponent is24, leaving slack
over the chosen20-bit per-step bound. The source and derivation also
passed independent static review.

The remaining global proof still needs messages with two or more active
outer rows; the one-row theorem must not be described as a full-code
distance certificate. The next step is to retain both the total shortened
dimension and the union of input-occupied regions when extending this
growth argument beyond one active row.

## Sparse-row extension and its frontier

The same argument gives a uniform class criterion beyond one row. Let a
be the number of true active input rows, let their row weights be
128+d_i, and put

    d=sum_i d_i,       s=max_i d_i.

Here a<=d<=128a and1<=s<=128. The count of all such messages at fixed
(a,d) is at most

    binom(M,a) binom(256a,128a+d) q^d.                 (5)

Indeed, choose the active rows and the input support positions, then use
the exact shortened dimension sum_i d_i=d. Dropping the individual row
weight requirements only enlarges the support count. There is no union
over row-weight tuples in (5).

Let k_a be a certified joint-zero threshold for messages supported on a
input rows. The fixed-row/message proof of the adaptive-zero bound applies
with this meaning of a, independently of the footprint's eligible count.
Outside that exceptional class put

    n_a=min(4096,floor((W+k_a-1)/512)).

This bounds fully active incoming updates. Every input row of maximum
weight128+s forces that many occupied regions. Repeating (1)--(3) gives

    E+B+33r >=2013+16s-n_a.                           (6)

For a concrete all-class interface, take tau=128 and define

    A_(a,d)=32d+1+ceilLog2 binom(M,a)
                    +ceilLog2 binom(256a,128a+d),
    e_(a,d)=max(24,ceil((A_(a,d)+tau)/14)),
    r_(a,d)=ceil((A_(a,d)+tau)/116),
    b_(a,d)=max(48,ceil((A_(a,d)+tau)/12)).

Each corresponding exceptional family has probability less than2^-tau
after the message union. Thresholds above a statistic's maximum describe
empty exceptional events. Thus all bad messages satisfying

    (e_(a,d)-1)+(b_(a,d)-1)+33(r_(a,d)-1)
                         <2013+16s-n_a               (7)

are excluded outside those families. For fixed(a,d), this covers every
qualifying s simultaneously; there is no s union. There are fewer than
2^32 admissible(a,d) pairs. Their three exception families cost less than
2^-94 in total. Including the all-a joint-zero budget below2^-79 keeps
this entire extension below2^-78.

Criterion (7) favors dimension concentrated in a few rows. It does not
uniformly close two active rows. For example, a=2,d=256,s=128 gives
A_(a,d)=8218 and thresholds(e,r,b)=(597,72,696). Its complementary score
upper bound is3634, whereas the right side of (7) is at most2218. This
is a limitation of the aggregate growth argument, not a bad codeword.
Replacing s by the total d in (6) would be incorrect.

A sharper interface can retain j, the number of regions containing any
nonzero input. Then Q>=j, and the message support count in (5) can be
replaced by

    binom(M,a) binom(256,j) binom(a*j,128a+d) q^d.

Use16j-n_a-35 on the right side of (7). The extra at-most256 choices of
j cost eight union bits, still within the stated conservative extension
budget. This retains overlap of the active rows' input supports without
enumerating state histories. It is a theorem interface awaiting a
coverage argument; no all-sparse-row conclusion is asserted.
