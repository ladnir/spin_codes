# Output-weight analysis for a single grouped outer

2026-09-23. A new outward-rounded calculation bounds the 10%-distance failure
contribution from all rank-one messages confined to one 16-row group by
approximately 3.03570e-13, or 41.583 bits. The union includes all 512 possible
groups at K=2^20 and every nonempty subset of rows within each group.

This improves the previous delayed-activation result: the calculation includes
emitted weights, later cancellations, and repeated returns to zero. It still
does not certify the candidate's full distance. Higher-rank tuples within one
group and messages spanning multiple groups remain unresolved. The best
confirmed implementation time remains approximately 7.86 ms against 9.57 ms;
no new benchmark or production change was made in this iteration.

## What failed, and what changed

The first output-weight calculation replaced every active column by the worst
allowed column shape. Combined with the shortened-support counts, its full
one-group union bound remained vacuous (binary64 log2 bound about +941.9).
This is a limitation of that relaxation, not evidence of a bad codeword.
It ignores constraints linking column shapes in a tuple of BCH words.

Rank-one tuples have a simpler description: a nonzero BCH word c is repeated
in some nonempty subset of the 16 rows. For a subset of size a, every active
column contains exactly a ones. Keeping this size fixed improves the bound.
However, summing all row subsets immediately gives only about 38.134 bits.

All these subsets share the same support placement for c. A late first active
region is therefore one common event, not a separate event for every subset.
Condition on that first region, bound the union over row subsets by at most
one, and only then average over the first region. This yields the 41.583-bit
result without changing the construction or its randomness.

## State and output envelopes

Each group shares one uniform coordinate permutation of the 256 BCH positions.
In each region, the group labels and the lanes within each group are permuted
independently and uniformly. These region draws are independent of each other
and of the coordinate permutations and transvections.

Use the existing expansion matrix E and feedback matrix B. For epoch input x
and entering state q, the output is x+Eq and the next state is Fq+Bx. The
transvection F is independent in each epoch. For each fixed nonzero q, Fq has
the exact distribution: keep q with probability 1/2, and otherwise sample a
uniform nonzero state. All probabilities here concern the uniform grouped
route and those transvections, with the BCH code and message fixed.

The current expansion map has nonzero output weights 48,56,64,72,80. Let S_v
count the nonzero states with expansion weight v, and put m=2^19-1. The
weighted-measure envelope has seven coordinates:

- Z bounds mass at state zero.
- D bounds arbitrary mass supported on nonzero states.
- L_v bounds a measure pointwise by L_v/S_v on states with expansion weight v.

The L_v coordinates are upper density bounds, not assertions that the true
conditioned state is uniform. The sum of all coordinates bounds total mass.
All transfers below weight an epoch by z to its emitted Hamming weight,
where 0<z<1. Their rows need not sum to one because they are upper envelopes.

For a fixed active-row count a, let X_a contain all weight-a patterns in the
eight aligned 16-position windows. Its size is n_a=8 binomial(16,a). The
active input is uniform in X_a. Define

    f_v(z) = (1/(S_v n_a)) sum_{q:wt(Eq)=v} sum_{x in X_a} z^wt(Eq+x),
    c_v(z) = (1/(S_v n_a)) sum_{x in X_a:Bx!=0, wt(E Bx)=v}
                                          z^wt(E Bx+x).

The script evaluates f_v through exact counts of ones in each 16-position
window. It evaluates c_v by enumerating the local input patterns. For arbitrary
nonzero state mass, use f_D(z)=z^(48-a). A bound on c_D is the minimum of f_D
and M_a z^b_a/n_a, where M_a is the largest feedback-syndrome multiplicity in
X_a and b_a is the least cancellation-output weight over its nonzero syndromes.

For an entering D or L_v coordinate, write f and c for the corresponding
functions. An active epoch has the following valid outgoing bounds:

    Z: c/2 + f/(2m),    D: f/2,    L_w: f S_w/(2m).

The first term in Z bounds lazy cancellation. The second bounds cancellation
after uniform refresh. The surviving lazy mass is sent to D. After adding Bx,
the refreshed part has pointwise density at most f/(2m) on each nonzero state.
The envelope can count some mass twice; this only increases the bound.

Starting from zero, let p_a=Pr[Bx=0]. The outgoing bounds are p_a z^a in Z
and (1-p_a)z^a in D. The exact census gives p_8=1/51480 and p_a=0 otherwise.

For an empty epoch, zero remains zero. An entering L_v contributes z^v/2
back to L_v and z^v S_w/(2m) to each L_w. An entering D uses the same refresh
rule with z^48, and keeps its lazy part in D. Thus long empty gaps retain
the actual mixing effect instead of repeatedly assuming a worst-case state.

The script extracts the production E and B columns and checks their identity.
It reconstructs the complete expansion spectrum and all local cancellation
counts. No older feedback-map certificate is substituted for the current map.

## Conditioning before the union bound

Each region contains 64 epochs. Let W_0 and W_a be the preceding empty and
weight-a active epoch envelopes. The region envelopes are

    R_0 = W_0^64,
    R_a = (1/64) sum_{j=0}^{63} W_0^j W_a W_0^(63-j).

The first occupied region leaves some l regions after it. Before that region,
the state and input are zero. Conditional on an outer support of size w and
this first occupied region, the remaining w-1 occupied regions form a uniform
subset of those l regions. Consequently the conditional moment is bounded by

    M_{a,l,w}(z) = e_0 R_a [x^(w-1)](R_0+x R_a)^l 1 / binomial(l,w-1).

Here coefficient extraction selects the indicated power of the formal marker
x; e_0 starts at zero state, and the final vector of ones sums envelope mass.
The scalar z weights output bits, while x counts occupied regions.

For D=209715, Chernoff's inequality bounds the probability of output weight
at most D by z^(-D) M_{a,l,w}(z). The implementation minimizes over seven fixed
positive tilts lambda=-ln(z), separately for each a,l,w. It also caps each
probability by one. Let p_{a,l,w} denote the resulting upper bound.

The probability that the first occupied region leaves l following regions is
binomial(l,w-1)/binomial(256,w). A union over every rank-one message therefore
has the upper bound

    512 sum_{w>=1} A_w sum_{l=w-1}^{255}
        [binomial(l,w-1)/binomial(256,w)]
        min(1, sum_{a=1}^{16} binomial(16,a) p_{a,l,w}).

The factor A_w uses the authenticated BCH shell cap. No independence between
row subsets, BCH words, or groups is required by these union bounds. Shared
first-support placement is used only before summing the row subsets for the
same fixed BCH word and group.

## Verification and limitations

Independent runs at 192-bit and 384-bit precision both give

    failure upper < 3.035696e-13 < 2^-41.

The numerical verifier constructs transfer entries and all polynomial
coefficients with Arb interval arithmetic. It rounds intermediate results
upward and checks the final inequality directly, not via a rounded logarithm.
The principal shell in the binary64 diagnostic is weight 38.

Run from the repository root:

    python -B research/workstreams/permutation_locality/test_group_moment.py
    python -B research/workstreams/permutation_locality/group_rank_one_verify.py --precision 192
    python -B research/workstreams/permutation_locality/group_rank_one_verify.py --precision 384

The first command exhaustively verifies the transvection law at state sizes
2 through 5. It checks the conditional coefficient recurrence against both
the unconditional recurrence and exhaustive placements at lengths 1 through 8.
It also compares independently constructed binary64 and outward transfers
and checks the empty-input invariants. Each verifier command authenticates the
BCH caps and rebuilds the production-map census and rank-one bound. Precision
is set after importing the historical authentication modules, which otherwise
set a default precision at import time.

This is a partial certificate for a specific new distribution, not reuse of
the paper's independent-row certificate. The failed worst-column calculation
suggests retaining more joint column information for higher ranks. Separately,
the implementation still needs a substantial improvement to reach 2x; closing
this partial proof would not by itself satisfy the performance goal.
