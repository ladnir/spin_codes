# Independent outer scalars: rigorous distance classes

Independent nonzero GF65536 scalars restore the outer measure comparison
used by the binary SPIN proof. They also give concrete global bounds for
the current shared beta64 construction. The following ranges are certified
by directed arithmetic, including their stated unions.

| Class | Number of active outer rows | Failure probability |
| --- | --- | --- |
| Any word of weight at most 838860, [state-density refinement](HD_RANDOMIZED_OUTER_STATE_PROFILE.md) | 1 | below 2^-391.9739 |
| Any word of weight at most 838860, [exact sparse ladder](HD_RANDOMIZED_OUTER_SPARSE_LADDER.md) | 1 through 5 | below 2^-56.5182 |
| Any word of weight at most 838860 | 7548 through 8192 | below 2^-326.6848 |
| Quiet word of weight at most 838860 | 6566 through 8192 | below 2^-1174.2244 |
| Any nonzero quiet word, regardless of weight | 1 through 953 | below 2^-142.1972 |
| Quiet word of weight at most 838860, [joint-moment refinement](HD_RANDOMIZED_OUTER_JOINT_QUIET.md) | 5940 through 8192 | below 2^-3832.4504 |
| Any nonzero quiet word, [exact-coefficient extension](HD_RANDOMIZED_OUTER_QUIET_EXACT.md) | 954 through 962 | below 2^-226.2046 |
| Quiet word of weight at most 838860, [exact-coefficient extension](HD_RANDOMIZED_OUTER_QUIET_EXACT.md) | 5865 through 8192 | below 2^-355.6458 |
| Quiet word of weight at most 838860, [tile and coefficient refinements](HD_RANDOMIZED_OUTER_JOINT_INPUT.md) | 1508 through 7574 | below 2^-72.6143 |

For a quiet word, the state is zero before every emitted packet. These
results leave quiet words with 963 through 1507 active rows and nonquiet
words with 6 through 7547 active rows. They are not a complete 40%/40-bit
distance certificate.

## Setup and probability space

Put q=65536, R=8192, n=256, k=128, N=nR=2^21, and P=N/64=32768.
One scalar lane contains R independent outer message rows in F_q^128.
Encode each row with the same fixed MDS[256,128] generator over the GF256
subfield. Extension to F_q preserves the nonzero MDS minors.

For every row i and outer coordinate j, independently sample

    d_(i,j) <- F_q \ {0}.

Multiply that outer output symbol by d_(i,j). Independently sample one
uniform permutation of the R rows in each of the n coordinate regions.
Route the scaled symbols through those permutations, retaining each
region's chronological order of 128 packets.

The inner recurrence is

    Y_t = F X_t + A Q_t,
    Q_(t+1) = M_t Q_t + C X_t,       Q_0=0.

Here F is the shared beta64 feedthrough. The required identities are
F^2=I, CF=TC, CA=0, and T invertible. The map C checks MDS[64,60,5].
There is no terminal-state constraint. Every M_t is chosen independently
of the outer scalars and route. The current construction also samples its
mixers independently between packets. The original dense and quiet bounds
permit arbitrary fixed mixers. The sparse state-density refinements use
the actual independent rank-one law. No uniform-GL assumption is used.

Setup is sampled once and reused across all messages. Failure means that
some nonzero message gives an output of weight at most D=838860.
The probability is over that setup. All whole-class events below quantify
over every message in the stated class, not one presampled message.

The native representation packs eight such F_q lanes into a 128-bit record.
When all lanes use the same scalar code, the code's symbol distance is the
minimum scalar-lane distance, without an eight-fold union. If independent
lane setups are selected instead, union their failure probabilities.
An aggregate length of eight lanes has a different distance normalization;
the N used here is the length of one scalar lane.

## Expected outer measure

Fix a nonzero scalar row output x with support S of size h. Shortening the
outer code to S has dimension h-128 when h>=129. Any h-128 coordinates
form an information set for that shortened code. An exact-support word
has nonzero values in those coordinates, so at most (q-1)^(h-128)
outer words have exactly support S.

For every such fixed outer word, its independent scalars send its active
symbols to independent uniform elements of F_q\{0}. Its probability of
becoming x is (q-1)^(-h). Consequently the expected number mu(x) of
nonzero row messages producing x satisfies

    mu(x) <= (q-1)^(-128).

For x=0 or 0<h<129, mu(x)=0. Thus, with U uniform in F_q^256,

    mu <= beta Law(U),       beta = q^256/(q-1)^128.

This is domination of expected counting measures. It neither asserts
independence between messages nor bounds every sampled enumerator.
The comparison deliberately includes artificial zero words.

Fix an active-row set S with |S|=r. Independent scalar setups between rows
give product domination by beta^r times uniform values in those rows'
256r potential positions. It holds for every fixed regional route.
All other positions are zero. The uniform comparison values include zero;
an active row does not imply 256 genuinely nonzero symbols.

The exact pointwise measure is also available. For a prescribed support of
size h, put a_0=1, a_h=0 for 1<=h<=128, and

    a_h = (q-1) sum_(j=0)^(h-129) (-1)^j binom(h-1,j) q^(h-129-j).

The expected multiplicity at a particular nonzero word of weight h is
a_h/(q-1)^h. This permits sharper sparse support folds later.

## Dense classes by counting inverse outputs

Fix the route and every mixer, independently of the outer scalars. The
entire inner map is invertible: each packet recovers X_t from Y_t and
the already recovered state using invertible F. Hence each output Y
determines one input X and one set of active outer rows after inverse
routing. Call the size of that set r(Y).

For this fixed inverse input, outer membership under the random scalars
has probability at most (q-1)^(-128r(Y)). There is no binomial union over
active-row subsets: Y already determines the only possible subset.

Let B be a scalar-invariant set of output words. Every sampled code is
F_q-linear, so each nonzero witness comes with all q-1 scalar multiples.
A union bound over these projective classes gives

    Pr[there is a nonzero Y in B with r(Y)>=r_0]
        <= |B| / (q-1)^(128r_0+1).

Counting zero in |B| only enlarges the right side. The bound holds for
each fixed route and mixer, and therefore after averaging them. This
step needs independence of the inverse inner map from the outer scalars.

For the ordinary weight ball and any rational 0<z<1,

    |B_D| <= z^(-D) (1+(q-1)z)^N.

At z=2/(3(q-1)), its base-two logarithm is at most

    15457972.0464719948269747264767736661610447706111726505271.

The inequality above certifies the entire r>=7548 class with logarithmic
failure upper bound -326.6848629258387367777791982780827979651.
For exactly r=8192, the bound is below 2^-1319236.8702.

## A smaller ball for quiet outputs

A word is quiet when Q_t=0 for t=0,...,P-1. This is equivalent to
CX_t=0 for t<P-1. The final packet has no feedback constraint because
Q_P is not emitted. Since CF=TC, quiet output packets satisfy CY_t=0
for t<P-1, while their last packet is arbitrary.

Let A(z) be the exact weight enumerator of MDS[64,60,5] over F_q:

    A(z) = sum_(w=0)^64 binom(64,w) b_w z^w,
    b_0=1,       b_w=0 for 1<=w<=4,
    b_w=(q-1) sum_(j=0)^(w-5) (-1)^j binom(w-1,j) q^(w-5-j).

The number of quiet outputs in the weight ball is at most

    z^(-D) A(z)^(P-1) (1+(q-1)z)^64.

The final factor includes the unconstrained final packet. At the rational
mark z=133/(128(q-1)), the logarithm of this count is at most

    13445991.2739180345763294430889312649831395793686508130473.

The dense inverse-output argument then certifies all quiet words with
r>=6566, with logarithmic failure upper bound -1174.22448200137953.
The GF256 quiet descent used for the earlier fixed outer does not apply
to this setup: the full-field diagonal scalars change those equations.
This calculation uses F_q throughout.

## Independent amplitudes and quiet syndromes

There are two useful local laws, under different probability spaces.
First fix one outer message and a realized route. Its genuinely nonzero
input positions receive independent uniform nonzero amplitudes from the
fresh coordinate scalars. If a packet has h such positions, then

    Pr[CX=0] = 1                              if h=0,
             = 0                              if 1<=h<=4,
             = b_h/(q-1)^h                    if h>=5,

where b_h is the prescribed-support count for MDS[h,h-4,5]. In particular,
the value at h=5 is (q-1)^(-4), close to 2^-64. Different packets use
disjoint scalar coins. This statement is made before conditioning on
their output weights or syndrome values.

For the second law, fix a potential-slot subset U of one packet and let
j=|U|. Under the uniform outer comparison, X_U is uniform in F_q^j.
Every restriction C_U has rank min(j,4), so

    Pr_uniform[CX=0 | U] = q^(-min(j,4)).

This law includes zero amplitudes. Its probability depends only on j,
even though U may have any positions. Therefore no averaging over the
within-packet subsets is required for this unweighted quiet bound.

## Quiet words at sparse row occupancies

Fix an active-row set S of size r. Let J_t count its potential slots in
packet t after routing, including positions whose comparison values are
zero. In each region, the r potential positions are a uniform r-subset
of R positions. The n regions are independent.

For a fixed route, the uniform comparison has quiet probability

    q^(-sum_(t=0)^(P-2) min(J_t,4)).

Let Z_S count projective nonzero quiet witnesses with active set S. The
outer envelope and scalar-orbit division give

    E_scalars[Z_S | route]
       <= beta^r/(q-1) * q^(-sum_(t=0)^(P-2) min(J_t,4)).

For 0<alpha<=1, conditional Markov and min(1,x)<=x^alpha bound failure
by the alpha power of the right side. Crucially, this cap is applied
before averaging rare route outcomes. No inverse expectation or
conditioning on quiet outputs is used.

Put

    H_alpha(x) = sum_(j=0)^64 binom(64,j) q^(-alpha min(j,4)) x^j.

One ordinary region contributes

    [x^r] H_alpha(x)^128 / binom(R,r).

The final region contributes

    [x^r] H_alpha(x)^127 (1+x)^64 / binom(R,r),

because its final packet is free. Their product, multiplied by
binom(R,r) beta^(alpha r)/(q-1)^alpha, bounds the whole r-row quiet class.
The union over active-row sets remains outside alpha.

For a rational marker p in (0,1), define

    S_p = sum_(j=0)^64 binom(64,j) p^j(1-p)^(64-j)
                               q^(-alpha min(j,4)),
    c_r = binom(R,r) p^r(1-p)^(R-r).

Nonnegative coefficient domination gives the computable bound

    Pr[there is a nonzero quiet word with r active rows]
       <= binom(R,r) beta^(alpha r)/(q-1)^alpha
          * S_p^(P-1) / c_r^256.

The marker introduces no change to the routing distribution. In particular,
the binomial normalization must be paid once per region, not once globally.

The helper selects rational witnesses from alpha=a/32, 1<=a<=32,
p=b/128, 1<=b<=127, and p=2^-j for 8<=j<=17. Floating arithmetic only
proposes each witness. Directed Arb evaluation verifies every integer
occupancy from 1 through 953 and their final sum. The total logarithmic
failure upper bound is -142.1972142213165099582265531398421676211.
The weakest occupancy is r=953, with alpha=3/32 and p=1/8.
This result excludes nonzero quiet words at any output weight in that range.

## Exact interface for the remaining composition

For a weight mark 0<z<1 and fixed potential subset U, define

    K_U(z) = q^(-|U|) sum_(x in F_q^U) 1[CX=0] z^wt(FX).

This is a joint quiet-and-output moment. Its input amplitudes are averaged
before applying alpha. Averaging K_U(z) over the uniform j-subsets of a
packet gives K_j(z). The final packet uses the same definition with the
syndrome indicator removed.

Condition on all potential-slot counts J_t. Within-packet subsets are
independent uniform subsets of their declared sizes. After uniform outer
comparison, the artificial values in different packets are independent.
Thus the conditional quiet first moment is bounded by

    beta^r z^(-D)/(q-1) * product_t K_(J_t)(z),

with the final-packet replacement. Conditional Markov now gives a regional
coefficient construction with K_j(z)^alpha in place of
q^(-alpha min(j,4)). It can tighten the unweighted sparse result.
Factoring K_j(z) into a quiet probability and an unconditional weight
moment would require a separate proof; those events share amplitudes.

Nonquiet words require a state envelope. An exact, possibly impractical,
interface uses states a,b in F_q^4 and defines the positive operator

    P_j(z)[a,b] = E_(U,X_U,M)[ z^wt(FX+A a) 1[M a+CX=b] ],

where U is a uniform j-subset, X_U is uniform in F_q^j, and M has the
actual rank-one mixer law. The state begins at zero; sum all final states.
This operator averages uniform amplitudes, subset positions, and the fresh
mixer before any fractional power. A compressed operator must dominate
these weighted measures, including correlations between output and state.

Fix the active-row set and condition on all packet counts J_t. Its
conditional expected projective witness count is bounded by the product
of these operators, times beta^r z^(-D)/(q-1). Apply conditional Markov
and the alpha power at this point. Average the count sequence afterward.

For any valid nonnegative matrix envelope P_j(z), group consecutive packet
operators into blocks. Expand scalar matrix paths, then use subadditivity
of x^alpha to replace each block product by its entrywise alpha power.
Average potential-slot placements only afterward. In a two-packet block,

    T_l = sum_j [binom(64,j) binom(64,l-j)/binom(128,l)]
                    * (P_j(z) P_(l-j)(z))^(entrywise alpha).

The sum is over feasible j. The hypergeometric factors stay outside the
fractional power. The regional operator is

    R_r = [x^r](sum_(l=0)^128 binom(128,l) T_l x^l)^64 / binom(R,r).

The whole r-row contribution is then at most

    binom(R,r) beta^(alpha r) z^(-alpha D)/(q-1)^alpha
        * e_0 R_r^256 1.

This preserves state across region boundaries. It identifies precisely
which part of the binary proof transfers: outer domination and the
route-conditioned composition do transfer. Its old binary local operator
does not transfer to GF65536 and the rank-one mixer without a new proof.
Independent outer scalars do not make the state uniform after conditioning
on earlier weighted outputs.

## Reproduction and next step

The new helper and seven tests are in
[theory/hd_randomized_outer_proof.py](theory/hd_randomized_outer_proof.py) and
[theory/hd_randomized_outer_proof_tests.py](theory/hd_randomized_outer_proof_tests.py).
Tests exhaust a small scalar ensemble, check exact MDS counts and both
syndrome laws, and verify regional normalization with the free final packet.
All displayed endpoints use 192-bit Arb enclosures; regression tests use
160 bits. Floating witness selection cannot turn an invalid endpoint into
an accepted one.

From the theory directory:

    python -B -m unittest hd_randomized_outer_proof_tests -v
    python -B hd_randomized_outer_proof.py

The next proof step is a useful compressed nonquiet operator for the
actual mixer, together with joint quiet-and-output moments. Evaluate every
remaining occupancy before claiming a complete certificate. No benchmark,
native implementation change, or commit was made for this note.

An independent mathematical audit checked the outer atom proof, inverse-output
union, sparse quiet cap, regional normalization, free final packet, and
weighted matrix composition. It found no mathematical blocker. The audit
did not regenerate numerical endpoints; the executable checks do that.
