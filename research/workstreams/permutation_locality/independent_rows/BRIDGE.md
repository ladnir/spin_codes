# From independent row shuffles to an occupancy certificate

This note connects the averaged outer counts to the conditional inner bound.
The connection is specific to the independent-row construction. It does not
transfer a certificate from the shared-row-shuffle construction.

The parameters are BCH[256,128], message length K = 2^20, output length
N = 2^21, and the two-update IMT(128,19) inner. There are G = 2048 groups,
each containing four adjacent BCH rows. A group is active when its four
input messages are not all zero. The bad event is an output of weight at
most D = 209715 from some nonzero message.

Setup independently samples a uniform coordinate permutation for each
BCH row, a uniform group permutation in each region, and a uniform lane
permutation for each four-bit packet. It also samples the two independent
transvections at every inner step. Setup is sampled once and then used for
every message. Probabilities below are over this ideal setup distribution;
they are not assertions about a particular sampled seed.

The formulas below use input-weight tilt one. An additional input-weight
tilt requires its reciprocal factor in the outer measure and valid caps
for that modified measure. The low-occupancy replays need no such factor.

## Conditional input distribution

Fix a message before sampling setup. For each active group g, let
H_g = (n_0,...,n_4) be its histogram of packet weights across the B = 256
regions, before lane shuffling. Thus

    U_g = B - n_0,       J_g = n_4

are its union-support size and its number of all-one packets.

For fixed row weights, the four independently shuffled row supports are
independent uniform subsets of the corresponding sizes. Their joint
distribution is invariant under a common permutation of region labels.
Consequently, conditional on any positive-probability histogram H_g, its
packet-weight sequence is uniform among the arrangements of that histogram.
The active-region set is therefore uniform among the U_g-subsets.

An independent lane permutation maps each packet to a uniform four-bit
mask of its weight. Conditional on the packet-weight sequence, these masks
are independent across packets. Histograms of different groups depend on
disjoint row-permutation randomness. Conditioning on all group histograms
preserves independence between groups.

This is the input distribution obtained by taking a fixed four-row array
with histogram H_g, applying a shared uniform column permutation, and then
applying independent lane permutations. This equality concerns the
conditional input distribution, not the outer tuple-counting measure.

The universal inner envelope is valid for every fixed sequence of packet
weights. Within a region, condition on the assignment of active labeled
groups to epochs, but leave window positions unexposed. Each epoch then
has an independent uniform injection into its 32 windows. These window
choices are independent of the entering state. The envelope averages those
choices, lane permutations, and inner updates before maximizing over
packet-weight shapes. This ordering is the reason its reuse is valid.

## The conditional inner inequality

Fix a tilt lambda > 0 and a penalty 0 < rho <= 1. The local operator
T_j bounds an epoch containing j active packets, after multiplying its
tilted transition by rho for each all-one packet. The multiplier is
applied before maximizing over shapes. The operators are nonnegative
11-by-11 matrices on the existing state-envelope coordinates.

Let e_0 select the zero initial state. Let tau select the mass coordinates
at termination, excluding the density coordinate and the two mature-tail
coordinates. Those three coordinates are auxiliary bounds, not additional
probability mass. The inner has no flush.

For r active groups in one region, define

    R_r = [x^r] (sum_{j=0}^{32} binomial(32,j) x^j T_j)^64
                    / binomial(2048,r).

The product follows epoch order. Coefficient extraction averages the
uniform placement into distinct windows; it does not replace placement
by independent slots. Labels can be removed because T_j bounds every
assignment of packet weights to those labels.

For q active groups, introduce formal variables x_1,...,x_q and define

    P(x_1,...,x_q) = sum_{S subseteq {1,...,q}} R_|S| product_{g in S} x_g.

Fix group histograms H_1,...,H_q, and write u_g and j_g for their union
sizes and all-one counts. Let W be the final output weight. The conditional
input argument and the shape-uniform local envelope give

    rho^(sum_g j_g) E[exp(-lambda W) | H_1,...,H_q]
      <= [x_1^u_1 ... x_q^u_q] e_0 P(x_1,...,x_q)^B tau
                    / product_g binomial(B,u_g).                 (1)

The coefficient enumerates every possible active-region arrangement.
The packet weights within an arrangement need not be independent: the
operator bounds all compatible weight shapes. Thus (1) does not pay an
additional probability denominator for the finer histograms.

For arbitrary witnesses 0 < p_g < 1, let independent auxiliary variables
B_g have Bernoulli distributions with parameters p_g. Define

    R(p) = sum_{r=0}^q Pr[sum_g B_g = r] R_r,
    beta_u(p) = binomial(B,u) p^u (1-p)^(B-u),
    M(p) = e_0 R(p)^B tau.

All coefficients in (1) are nonnegative. Equivalently, average its
coefficient numerator under independent Bernoulli support indicators and
discard the nonnegative terms of other degrees. This yields

    Pr[W <= D | H_1,...,H_q]
      <= exp(lambda D) M(p)
           product_g (rho^(-j_g) / beta_u_g(p_g)).                (2)

The variables B_g are coefficient witnesses, not replacement randomness
in the construction. For a coordinate restricted to u_g = B, p_g = 1 is
also valid: evaluate that deterministic coordinate directly and use
beta_B(1) = 1. No zero binomial mass may be divided into a bound.

## Average over row setup and sum messages

For one group, let m range over its 2^512 - 1 nonzero ordered message
tuples. The corresponding four BCH words are fixed once m is fixed.
Define the weighted averaged shell count

    a_rho(u) = sum_{m != 0} E_row[rho^(-J(m)) 1_{U(m)=u}].

The expectation here uses only this group's four row permutations.
These counts may be rational; they need not count words in any one
realized setup. The identities and caps in [README.md](README.md) bound
this measure. In particular, the reciprocal penalty is retained inside
the expectation; it is not replaced by a typical value of J.

Fix a set of q active groups and a labeled support box
I_1 x ... x I_q. Let Z_I count messages on these groups whose output is
bad and whose shuffled support vector lies in that box. Sum (2) over
messages and average over the histograms. Independence of row setup
between groups and the Cartesian product of group message spaces give

    E[Z_I] <= exp(lambda D) M(p)
                product_g sum_{u in I_g} a_rho(u)/beta_u(p_g).    (3)

This factorization uses independence between groups for a fixed message.
It does not assume that two messages have independent setup. Linearity
of expectation permits the sum over all messages despite their shared
permutations and inner maps.

All sets of q active groups have the same bound, because the region
permutations treat their labels symmetrically. Multiply (3) by
binomial(G,q) to include their locations.

## Fold CDF caps and cover all support vectors

Suppose the nondecreasing function C_rho(u) bounds
sum_{v<=u} a_rho(v). For an interval I = [lo,hi] and its chosen p, put

    w(u) = 1/beta_u(p),
    wbar(u) = max_{u<=v<=hi} w(v).

The restricted prefix sum from lo through u is at most C_rho(u).
Summation by parts, with the decreasing majorant wbar, gives

    sum_{u=lo}^{hi} a_rho(u) w(u)
      <= C_rho(lo) wbar(lo)
         + sum_{u=lo+1}^{hi} (C_rho(u)-C_rho(u-1)) wbar(u).       (4)

The first term does not subtract C_rho(lo-1). Differences in (4) belong
to an extremal upper calculation; they are not asserted shell bounds.
Rounding a rational nondecreasing CDF upward to integers preserves its
validity in (4).

Individual shell caps S_rho(u) >= a_rho(u) give another valid bound,
sum_{u in I} S_rho(u) w(u). The smaller of this bound and (4) remains
valid. The shell caps must come from the positive weighted subset measure,
not from differences of C_rho.

The optional prefix-rank fold combines these constraints more directly.
For every index set A contained in I and every cut k in I,

    sum_{u in A} a_rho(u)
      <= C_rho(k) + sum_{u in A, u>k} S_rho(u).

The shell-only bound sum_{u in A} S_rho(u) is also valid. Let r(A) be the
minimum of these bounds. Order the interval's indices as i_1,...,i_m so
that w(i_1) >= ... >= w(i_m), and put w(i_{m+1}) = 0. Then

    sum_{u in I} a_rho(u) w(u)
      <= sum_{j=1}^m (w(i_j)-w(i_{j+1})) r({i_1,...,i_j}).

This follows by expanding the weighted sum as the corresponding positive
combination of nested-set masses. It needs no interpretation of CDF
increments as shell counts. In outward evaluation, first replace each
w(u) by an exact upper endpoint, then order these endpoints. The set caps
remain exact rationals; outward arithmetic encloses the final weighted sum.

Each nonzero group contains a BCH word of weight at least 38. Therefore
its union size belongs to {38,...,256}. The cover partitions the full
ordered support domain {38,...,256}^q into labeled boxes. It stores a
sorted box with multiplicity L when L distinct assignments of its
intervals to group labels have the same bound. Multiplying (3), with
(4) substituted, by L and binomial(G,q) accounts for all these assignments
and all group locations.

Each box may choose its own lambda, rho, and p. Equation (3) is valid
separately for every such choice, so adding their bounds is valid.
Selected parents replace their entire subtrees: a parent and its children
are never both included. The split construction, not volume equality
alone, establishes disjointness and completeness of the partition.

The sum U_q over a complete cover bounds the expected number of bad
messages with exactly q active groups. Markov's inequality gives

    Pr[there is a bad message with exactly q active groups] <= U_q.

Finally, a full-code certificate requires bounds for every q = 1,...,2048
and an outward aggregate sum_q U_q below the claimed failure budget.
For example, a 40-bit claim requires that aggregate to be below 2^-40.
A certificate at one occupancy does not cover intervening occupancies.
Selected homogeneous support vectors also do not form a support cover.

## What the verifier checks

`verify.py` uses the new averaged outer counts, the universal inner
operators, and the existing support-cover helper. Its search arithmetic
only proposes witnesses. A successful run reconstructs the complete
selected cover using outward Arb arithmetic, including binomial masses,
folded counts, location factors, label multiplicities, and the final sum.

The local census tests and cover regression tests check implementation
identities. They do not substitute for the conditioning, positivity, and
independence arguments above. This lemma justifies the reduction; it does
not assert that the complete occupancy sum has already met its budget.
