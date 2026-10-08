# An exact source law for the affine one-prefix candidate

The fixed affine placement has a useful proof interface: integrate out the
independent prefix weights exactly, then bound a rank generating function of
the remaining map. This avoids treating the fixed Cayley neighbors as
independent coordinates. It also retains the actual regional routing law.

The bounded calculations below validate this interface. They do not certify
40% distance, find a counterexample to the frozen K20 stride family, or make
a native performance prediction.

Update: the [component-class analysis](HD_COMPONENT_RANK.md) proves that
one common Chernoff mark cannot close the relaxed bound below. The
identities remain valid. Splitting its positive expansion by permitted-row
count permits separate marks and now gives rigorous bounds for the entire
singleton term and [pair term](HD_COMPONENT_PAIRS.md).

A subsequent [actual-code obstruction](HD_COMPONENT_DENSE.md) shows that
even per-class or per-row-set marks cannot close the full uncapped first
moment. Rare routes contain large low-weight subspaces satisfying the actual
outer checks. The identities below remain valid, but full certification now
requires a capped or routing-conditioned rank argument. The first
[conditioning gate](HD_ROUTING_CONDITION.md) bounds one such route family;
it is not a whole-code certificate.

## The map and the remaining randomness

There are R=32768 source rows. Each row is a GRS[64,32] word over GF65536.
Thus K=32R and N=64R. The fixed placement S sends coordinate c of row r to

    c*R + ((a_c^(-1)*(r-b_c)) mod R).

The odd strides are the frozen family from `HD_TWO_STATE_REGIONAL.md`.
There is one input-weighted prefix, with independently sampled nonzero
weights at all N positions. Its state carries across region boundaries.

Every identity in this note also applies to any other fixed bijective placement
S, including the native four-row-group placement. That placement is a changed
construction. The individual-row Cayley spectral receipt does not transfer to
it, and the small examples below certify neither placement's distance.

Next, Pi applies one independent uniform permutation within each coordinate
region. The post-inner B has two globally carried additive states:

    y_i = s1 + f_i*s2 + c_i*x_i;
    s1' = s1 + d1_i*x_i;
    s2' = s2 + d2_i*x_i.

Each tuple consists of nonzero field elements and is conditioned on
c_i+d1_i+f_i*d2_i being nonzero. Tuples at different positions are independent.
All prefix weights, post tuples, and regional permutations are independent.
Both states start at zero; terminal states are free.

Let Delta be the global first-difference operator, with the predecessor of
position zero fixed to zero. For a fixed realization of B and Pi, define

    A = S^(-1) Delta Pi^(-1) B^(-1).

This is an invertible N-by-N matrix. It maps an output word to unweighted
prefix differences in source-row order. The actual inverse source word is
diag(alpha) A y, where the alpha entries are independent uniform nonzero
field elements. Conjugating the prefix diagonal by S only reorders them.

The argument can be applied to scalar GF65536 codewords. Using the same map
on all field lanes in a 128-bit record gives the same minimum record distance:
a nonzero lane gives a scalar word with no larger support, and one live lane
attains the scalar minimum.

## Exact membership after the random diagonal

Fix y before sampling alpha, and write v=Ay. Let w_r be the number of nonzero
coordinates in source row r of v. Every nonzero coordinate of diag(alpha)v is
independent and uniform over F*. Coordinates that are zero remain zero.

For an MDS[b,k] row over a field of size q, let rho=b-k. The number of codewords
with one prescribed exact support of size w is

    M_0 = 1;
    M_w = 0                                      for 1 <= w <= rho;
    M_w = (q-1) sum_{j=0}^{w-rho-1}
                  (-1)^j C(w-1,j) q^(w-rho-1-j)  for w > rho.

The usual MDS shortening dimension is w-rho. Inclusion-exclusion over zero
coordinates in that shortened code gives the formula. Therefore

    Pr_alpha[diag(alpha)v belongs to the outer code]
      = product_r M_(w_r)/(q-1)^(w_r).

This is an exact source law, including infeasible row supports. It does not
depend on an independent-neighbor approximation or a global permutation.

There is a simpler upper bound. Projection onto any w-rho information
coordinates is injective on the shortened code. A word with full prescribed
support projects into (F*)^(w-rho), so M_w <= (q-1)^(w-rho). Consequently every
active source row contributes at most

    z = (q-1)^(-rho).

Rows with 1 through rho active coordinates still have exact probability zero.
Dropping that feasibility condition gives the block-rank envelope below.

## A positive rank contraction

Fix an output support T and let A_T contain its columns. For a source-row set U,
let A[U,T] contain all b coordinates from each row in U and the columns in T.
The rows in U will be constrained to zero; U is not an active-row set.

For a vector y supported on T, write a(y) for the number of active source rows
in A_T y. Expanding each factor

    z^[source row is active] = z + (1-z)[source row is zero]

and then counting solutions to the selected zero equations gives

    sum_{projective nonzero y on T} z^a(y)
      = 1/(q-1) sum_{U subset [R]}
          z^(R-|U|) (1-z)^|U|
          (q^(|T|-rank A[U,T]) - 1).

“On T” permits zeros inside T. One representative of each nonzero scalar line
is counted, because scalar multiples describe the same codeword event. The
subtracted one removes the zero vector exactly. The right side is an upper
bound on the expected number of projective outer-code witnesses for T.

A second expansion retains exact output weights without summing overlapping
support events. Fix 0<v<=1. Sample auxiliary sets U and T independently:
include each source row in U with probability 1-z, and each output position in
T with probability v. These sets are algebraic devices, not routing choices.
Then

    sum_{projective nonzero y} v^wt(y) z^a(y)
      = E_(U,T)[q^(|T|-rank A[U,T]) - 1] / (q-1).

For a fixed word, the probability that T contains its support is v^wt(y).
This establishes the second identity from the first. All summands are
nonnegative, including the rank-defect term after subtracting one.

Let e be a forbidden output-weight cutoff. Averaging over the actual
conditioned post setup and the independent regional permutations gives

    Pr[distance <= e]
      <= v^(-e)/(q-1)
           E_(B,Pi,U,T)[q^(|T|-rank A[U,T]) - 1].

This is a valid certificate target, not a successful bound yet. The exact
source law can sharpen it by retaining infeasible row supports. The displayed
bound already uses the conditioned post ensemble; it has no extra conditioning
factor. An analysis performed in the independent unconditioned post ensemble
must instead pay the previously established factor Z^(-N), once.

## Contracting equalities before analyzing the post-inner

The selected source-zero constraints have a concrete implementation. For each
r in U and every coordinate c, mark the physical prefix-difference position
corresponding to (r,c). A marked position imposes u_i=u_(i-1), with u_-1=0.

These are edges of one global path, including its anchored initial vertex.
Contract the marked edges. The remaining free components number exactly

    m = N-b|U|.

Every component carries one field value. A component can cross a coordinate
boundary, so independent regional resets would give the wrong space. Apply
the actual regional permutations to these component labels, then feed their
values into B. Let W_U map the m component values to the final output.
For every T,

    |T| - rank A[U,T]
      = m - rank W_U[outside T, all components].

Both sides are the dimension of the same output subspace. The right side
describes repeated component values entering the two-state inner. It is a
more concrete next proof object than an arbitrary dense inverse matrix.

The affine geometry determines the joint component lengths along all 64
paths. Regional permutations randomize their positions but do not resample
component values. A proof must retain these repeated-value constraints and
the globally carried post states. Replacing labels by independent live flags
would discard the property this contraction is intended to capture.

## The generic-rank test is necessary but insufficient

For a fixed T, let H be the block-diagonal MDS parity-check matrix. The actual
check matrix on that support is H diag(alpha) A_T. Its generic rank is

    min_{J subset [N]}
       (sum_r min(rho, |J intersect row_r|) + rank A[J^c,T]).

Indeed, Cauchy-Binet expresses each minor as a sum of distinct diagonal
monomials. A monomial has a nonzero coefficient exactly when its coordinate
set is independent in both the MDS column matroid and the row matroid of A_T.
The common-independent-set rank formula yields the displayed minimum.
For q>2, a nonzero multiaffine polynomial cannot vanish on all nonzero diagonal
assignments. Thus the generic rank is attainable over the specified field.

A value below |T| proves a deterministic failure on T. A value equal to |T|
does not provide a 40-bit setup guarantee. In particular, a one-minor
Schwartz-Zippel estimate cannot support a union over all low-weight supports.
The rank moment, not just generic full rank, remains the certificate target.

## A dense-row shortcut that cannot work

For comparison, suppose every inverse word of weight at most e had at least
a active source rows. The projective count would give the upper bound

    e C(N,e) (q-1)^(e-1-rho*a).

At K20 and e=floor(0.4N)=838860, a directed-Arb calculation finds that a=30192
rows suffice for a bound below 2^-40. This is 92.1387% of the source rows;
the displayed bound is below 2^-317 at that integer threshold.

That uniform premise is impossible, regardless of affine expansion. Fix any
19661 source rows and any e output positions. Vectors whose source differences
vanish outside those rows form a 64*19661-dimensional space. Its output image
intersects the e-dimensional supported-output space nontrivially because

    64*19661 + 838860 > N.

Thus some nonzero output on those positions has at most 19661 active source
rows. This observation is not a codeword counterexample: it precedes the
random diagonal MDS membership test. It explains why a proof must retain
dimensions or weighted counts rather than require every inverse word dense.

## Bounded receipts

The new helper and ten tests use exact rational arithmetic and GF4 operations.
The only K20 numerical calculation here is the directed dense-row threshold.

For a complete small map with GRS[4,2], R=2, and output support {0,3,6}, the
test exhausts all 3^8 prefix diagonals. Exactly 459 have rank two; the remaining
6102 have rank three. Consequently the failure probability is 17/243. Every
failing kernel is one-dimensional, so the exact expected number of projective
kernel lines is also 17/243. The block-rank envelope gives 29/81.

The tests separately compare regional value profiles. At R=4, each region
has two live positions in both examples, but its prescribed value multiset is
either (0,0,1,1) or (0,0,1,2). Their exact average source-membership probabilities
under common placement are respectively

    5645/86093442 and 21797/1377495072.

Thus live occupancy alone does not determine the inverse source law.

At R=8 with two repeated nonzero values per region, the toy strides (1,3,5,7)
give fewer feasible arrangements than common order: 76 versus 193 out of
614656. Nevertheless, the weighted source probability is larger:

    affine: 113/1296243648;
    common: 5627/108884466432.

A rare one-active-row arrangement accounts for the direction of this change.
This example does not evaluate the frozen K20 family. It rules out replacing
the weighted proof object by a count of feasible routed masks.

Other tests exhaust every GF4 row vector and nonzero diagonal, check the MDS
support enumerator at both field sizes, verify global prefix carry, and compare
both rank identities with projective enumeration. All 1024 pairs of source-row
and output-position subsets in the small map satisfy the contraction identity.

```text
python -B research/workstreams/char2_field_spin/theory/hd_affine_proof_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_affine_proof.py
```

The original sufficient all-support condition was to find a mark 0<v<=1
and prove, for the chosen fixed placement S and actual conditioned setup,

    E_(B,Pi,U,T)[q^(|T|-rank A[U,T]) - 1]
      <= (q-1) * 2^(-40) * v^838860.

The random sets U and T have the distributions defined above. This one
inequality would cover every output support up to the target cutoff.
The subsequent component-class analysis rules it out for a common mark:
sparse and dense terms demand different marks. It replaces this condition
with a positive class-specific sum, preserving all terms. Checking a few
supports, mean ranks, or generic full rank cannot establish that sum either.

The next proof step is a bound on this rank-defect moment after component
contraction, with the chosen fixed paths and actual regional label routing.
A useful first checkpoint would retain a bounded family of component-length
profiles and identify which profiles dominate the moment. A full certificate
must then bound the omitted profiles and all relevant row-set sizes; the
quarter-row spectral receipt alone does neither.

No native implementation, build, benchmark, remote call, or existing-file edit
was performed. The changes are this note and the two new `hd_affine_proof`
helpers.
