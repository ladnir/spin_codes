# XOR2 column potentials with two active constituents

An eight-region potential closes the classes with one through four active
groups, each containing exactly two active RS constituents. The directed
log2 bounds are approximately -140.9034, -108.4888, -76.4422, and -44.7680.
Their union has directed log2 bound below -44.7679755806.

Each class includes every group choice and all 28^r constituent-pair choices.
These are class certificates for the unchanged scalar-route XOR2 construction.
Groups with more active constituents, five or more two-constituent groups,
and mixtures of constituent counts still require coverage. A bounded
five-group FLOAT screen below does not reach 40 bits; it does not establish
that this potential family or the construction must fail.

## Column operators retain the cancellation mixture

Use the D2 H D1 construction and scalar inner from
[HD_XOR_MULTIGROUP.md](HD_XOR_MULTIGROUP.md), with H=I+11^T over F=GF256.
There are 8192 groups, eight regions per outer column, and 32 columns.
Each active group in this note has two nonzero RS[32,16] constituents.

Fix the message and its active constituent pairs. At one column, let S_g
be the subset of those two constituents whose input coordinates are nonzero.
The nonzero D1 diagonals give the following exact output-support law:

    S_g empty:     T_g=empty;
    S_g={i}:       T_g=[8]\{i};
    S_g={i,j}:     T_g={i,j} with probability 1/255,
                   T_g=[8]   with probability 254/255.

The two alternatives in the last line describe equal-value cancellation
and its complement. The laws are independent across groups and columns.
Conditioned on the output supports, D2 makes every occupied amplitude an
independent uniform element of F*. These are setup probabilities for each
fixed message; no message is sampled from a surrogate distribution.

Fix an active-group count r. For 0<=m<=r, let R_m be the two-state operator
averaged over all subsets
of m positions in one region. It uses the nonzero-amplitude packet caps
and feedback cancellations proved in HD_XOR_MULTIGROUP.md. For output masks
T_1,...,T_r, define the chronological column operator

    Q(T_1,...,T_r) = R_(m_0) R_(m_1) ... R_(m_7),
    m_i = sum_(g=1)^r 1_(i in T_g).

For source supports S_1,...,S_r, define

    M(S_1,...,S_r) = E[Q(T_1,...,T_r) | S_1,...,S_r].

This expectation is evaluated before the potential inequality below.
For example, one double input averages its two complete eight-region
operators with weights 1/255 and 254/255. Two double inputs retain all
four cancellation outcomes with their product probabilities.

## A common potential for every source-pair choice

Choose 0<t<=1 and positive column marks f_1,f_2. Put v=(1,t)^T.
For one column, let n_1 count groups whose source support has size one,
and let n_2 count groups whose source support has size two. Verify

    M(S_1,...,S_r) v <= f_1^(n_1) f_2^(n_2) v.             (1)

The inequality is componentwise and must hold for every source-pair choice
and every source-support tuple. Group labels can be exchanged in this
check: adding their output masks produces the same regional counts.
Region order is preserved throughout.

The original two-group checker uses the following cases:

| (n_1,n_2) | Source-support cases |
| --- | ---: |
| (0,0) | 1 |
| (1,0) | 8 |
| (0,1) | 28 |
| (2,0) | 64 |
| (1,1) | 224 |
| (0,2) | 784 |

Thus 1109 checks cover every one of the 28^2 active-pair assignments and
all 16 source-support tuples within each assignment. The last category
explicitly includes every ordered pair of constituent pairs. No multiplicity
from this table is added to the eventual message union; these are universal
inequality checks.

For three and four groups, the checker uses the exact group-permutation
reduction in the next section. This reduction gives the same two-group
support laws as the original ordered checks.

For a fixed message, multiply (1) over its 32 columns. The initial state is
zero. Comparing the terminal all-ones vector with v costs 1/t once. If
N_1 and N_2 count its singleton and double source columns across all groups,
then

    E[z^W] <= (1/t) f_1^(N_1) f_2^(N_2).                  (2)

The expectation includes outer diagonals, independent regional permutations,
and inner coins. The regional counts within a column may be correlated.
The definition of M preserves that correlation and the eight-region order.

## Exact RS overlap count and projective union

Let B_w be the number of RS[32,16] words on one prescribed support of size w.
It vanishes for 1<=w<17 and has the MDS support formula recorded in
HD_XOR_MULTIGROUP.md. For two active constituents of weights w_1,w_2 and
support overlap h, the number of ordered support pairs is

    binom(32,w_1) binom(w_1,h) binom(32-w_1,w_2-h).

The pair has w_1+w_2-2h singleton columns and h double columns. Consequently
the exact marked message count for one prescribed constituent pair is

    T(f_1,f_2) = sum_(w_1,w_2=17)^32 B_(w_1) B_(w_2)
        sum_(h=max(0,w_1+w_2-32))^min(w_1,w_2)
          binom(32,w_1) binom(w_1,h) binom(32-w_1,w_2-h)
          f_1^(w_1+w_2-2h) f_2^h.

Both support floors remain inside this sum. In particular, every valid
constituent pair overlaps in at least two columns.

For r active groups, applying Markov to projective message classes
and using (2) gives

    failure probability <=
       binom(8192,r) 28^r T(f_1,f_2)^r z^(-734003) / (255 t).    (3)

All nonzero global scalar multiples have identical output weight. Therefore
the denominator contains one factor of 255 for the complete message.
Each active group independently chooses one of 28 unordered constituent
pairs. Its two constituent messages remain ordered by their row labels.

## Exact multiset reduction and labelled-region factorization

Within a column, exchange group labels but keep all eight region labels.
A singleton source chooses one missing output region. Record these choices
by a multiset with multiplicities s_i. A double source chooses one of the
28 region pairs. Record these choices by a multigraph G on the eight
labelled region vertices, allowing repeated edges.

Let C be the submultiset of double groups whose two nonzero D1 values
cancel. Each edge instance belongs to C independently with probability
1/255. Its output occupies its two endpoints; an uncancelled double
occupies all regions. Therefore the complete occupancy vector is

    m_i = n_1+n_2-|C|-s_i+deg_C(i),   0<=i<8.              (4)

This expression is unchanged by a permutation of group labels. Thus one
check for each singleton-row multiset and double-edge multiset covers
every ordered source-support assignment. For fixed (n_1,n_2), the number
of checks is

    binom(7+n_1,n_1) binom(27+n_2,n_2).

Summing over n_1+n_2<=r gives binom(36+r,r) checks. There are 37, 703,
9,139, and 91,390 checks at r=1,2,3,4. The r=2 count is smaller than
1,109 because this reduction also exchanges singleton or double groups.
It does not identify any region permutations.

Repeated edge instances retain independent cancellation coins. For n_2
double groups, the implementation sums all 2^n_2 branches, giving a branch
with k cancellations probability 254^(n_2-k)/255^n_2. If several branches
produce the same occupancy vector, their probabilities add. For example,
k repeated copies of one edge have exactly c cancellations with probability
binom(k,c)254^(k-c)/255^k. These binomial multiplicities belong inside M.
The multiset-case count is not an additional factor in the message union (3).

Directed replay stores four-region products P=R_(m_0)...R_(m_3) and
four-region vectors h=R_(m_4)...R_(m_7)v. It evaluates Qv as P h,
then averages the cancellation mixture and checks (1). Both tables have
(r+1)^4 entries: only 625 matrices and 625 vectors at r=4. The first four
labelled regions are always on the left. This factorization changes neither
chronology nor the inequality being checked.

## Saved witnesses and validation

The state potential and column marks are defined by rational base-two
exponents. All z values below are rational.

| Active groups | z | log2(t) | log2(f_1) | log2(f_2) | Directed log2 bound |
| ---: | --- | --- | --- | --- | ---: |
| 1 | 499451537/500000000 | -3991/64 | -343125/4096 | -342223/4096 | -140.90342915 |
| 2 | 199748791/200000000 | -499/8 | -12225/256 | -194985/4096 | -108.48876086 |
| 3 | 62411847/62500000 | -1997/32 | -145985/4096 | -72745/2048 | -76.44222390 |
| 4 | 998438053/1000000000 | -3995/64 | -120989/4096 | -120555/4096 | -44.76797558 |

The saved replay checks every case of (1) and evaluates (3) with Arb.
The receipts agree at 256 and 384 bits. Their FLOAT proposals were about
-140.9223, -108.5923, -76.4674, and -44.7894. No floating optimizer runs
during saved replay. Summing the four directed receipts at 384 bits gives

    [-44.76797558061644571538022952297023344339665534064019 +/- 3e-51]

for the log2 bound on their union. This interval is displayed at reduced
precision; the script's individual receipts retain their full precision.

The implementation is
[theory/hd_xor_column_potential.py](theory/hd_xor_column_potential.py).
Its [sixteen tests](theory/hd_xor_column_potential_tests.py) pass. They check
all 12,544 source-pair/subset combinations against the independent support
kernel, preserve noncommuting region order, exhaust a small-field RS overlap
enumerator, verify mixture placement, and replay the directed witnesses.
An invalid column mark is rejected. Further tests compare the exact reduced
and ordered support-law sets, check repeated-edge binomial multiplicities,
compare all 3^8 noncommuting test products with the four-plus-four factorization,
and compare streamed count codes and FLOAT costs with materialized cases.

Run from the theory directory:

```text
python -B -m unittest hd_xor_column_potential_tests -v
python -B hd_xor_column_potential.py --groups 1 2 --replay --precision 384
python -B hd_xor_column_potential.py --groups 1 2 3 4 --multiset --replay --precision 384
python -B hd_xor_column_potential.py --groups 5 --stream-mark .99829 --stream-potential -62.375
```

## Five-group boundary and next proof step

The streaming FLOAT screen checks all 749,398 five-group multiset cases.
It stores 1,296 four-region prefixes and processes at most 1,024 source
cases at a time, with at most 32 branches each. It does not materialize
the 6^8 eight-region word table or all five-group source cases.

The following bounded set of five points was evaluated; it is not a
certified optimization over z or t.

| z | log2(t) | FLOAT log2 bound |
| --- | ---: | ---: |
| .99824 | -62.4375 | 26.51213015 |
| .99829 | -62.4375 | -13.45161484 |
| .99834 | -62.4375 | -2.79563670 |
| .99829 | -62.375 | -13.61840264 |
| .99829 | -62.5 | -13.28068107 |

At the central point, the pure five-singleton and five-double categories
set log2(f_1)=-25.83179428 and log2(f_2)=-25.73756119. The latter worst case
includes all five double groups using the same constituent pair. The
separable marks charge the resulting per-group bound even to columns with
fewer simultaneously active groups. Although T retains both constituents'
support floors within each group, T^5 does not retain support coincidences
between groups. This identifies a concrete loss in the bound; it is not a
demonstrated obstruction to the unchanged construction.

The next proof step should keep the category contractions c_(n_1,n_2)
instead of replacing them by f_1^n_1 f_2^n_2. For example, introduce one
support variable for each of the 2r constituent words. A one-column
polynomial sums the 4^r source-support patterns, weighted by the appropriate
c_(n_1,n_2) and by the variables for nonzero source coordinates. Its 32nd
power counts complete support assignments, including cross-group overlap.
Multiply the coefficient at a weight vector by the product of its exact
prescribed-support counts B_w, retaining each w>=17. This is an exact,
finite support-count formulation after the category potential bound.
A direct 2r-variable coefficient table is too large; a useful implementation
needs symmetric-state compression or a rigorous envelope for the count of
fully occupied columns. A two-column potential could also retain more state
chronology, but it should be tested against a small worst-case family before
expanding the full support enumeration. No constructor alteration follows
from these failed FLOAT points.

For arbitrary r, the same multiset proof remains valid, but its check count
binom(36+r,r) and literal 2^n_2 cancellation expansion do not scale to the
whole family. Repeated-edge branches can instead be aggregated by their
cancellation-count vectors, with the binomial weights above, or convolved
directly into occupancy-vector probabilities using (4). This preserves all
eight labelled regions but can still require up to (r+1)^8 occupancy states.
The present code deliberately caps directed checks at four groups and
streaming FLOAT checks at five. No encoder change or benchmark was made.
