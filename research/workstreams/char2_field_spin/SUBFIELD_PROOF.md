# Proof status for the batched GF(256) candidate

The exact one-message transfer does not suffice for this construction.
It certifies the classes with one or two active outer rows, but a rare
routing event contributes a huge subspace of bad messages. Counting
every bad projective line makes the complete first moment unusable.
The construction itself is not disproved. The next proof must count
constraints or minimal witnesses rather than all messages.

## Fixed candidate and partial certificates

The current proof model has K=65536, N=131072, GRS[16,8] outer blocks,
L=8192 outer rows, and t=16-symbol inner steps. Each outer column is
permuted independently within its own L-symbol region. Independent
uniform GF(256)* coefficients multiply the input symbols. The inner
uses a two-symbol state, pre-output refresh by a uniform SL(2,256)
matrix, and the A_j=(1,j), C=A^T maps of
[BATCHED_PROOF.md](BATCHED_PROOF.md). The target excludes output
weights at most e=13107, giving distance at least 13108.

The whole generator has GF(256) coefficients. Its minimum symbol
distance is preserved under the declared scalar extension; the proof
below concerns base-field words, not arbitrary extension-field states
treated as uniform GF(256) pairs.

`theory/subfield_certificate.py --certify` evaluates rational witnesses
with 80-digit directed interval arithmetic. Integer powers use binary
exponentiation of outward intervals. The proof calculation uses no
transcendental functions. It verifies the following upper bounds:

| Active outer rows | Tilt z | Failure contribution, rounded upward here |
|---:|---:|---:|
| 1 | 8183/8192 | 1.031016e-15 |
| 2 | 511/512 | 7.830859e-14 |
| Both classes | Separate tilts above | 7.933961e-14 < 2^-43 |

The first class uses exact MDS support multiplicities and exact regional
transfers. The second uses the source envelope defined next. These are
certified contributions from two classes, **not a complete distance
certificate**. The record explicitly leaves active-row counts 3 through
8192 uncovered. Its verifier rejects a claim of complete coverage.

## A coupled bound and why improving its source count is insufficient

Fix z in (0,1]. Let R_h(z) be the normalized two-state transfer of one
region with h active input positions. Its definition includes the exact
binomial denominator for that region. For support masks of size at
least two, the implementation uses an entrywise upper envelope for
the joint return moment:

    p_ret z^t <= J_S(z) <= p_ret z^max(0,t-1-|S|),
    p_ret = (1-p0(|S|))/(256^2-1).

The lower bound follows because an output has at most t nonzero symbols.
The upper bound follows because A times a nonzero state has at most
one zero coordinate. The exact singleton transfer is retained.
In the nonreturn entry, subtract the lower bound on J_S from its exact
output moment. Thus every entry remains an upper bound before region
matrices are multiplied.

For a prescribed support of w>=9 coordinates in one GRS[16,8] row,
the number of codewords with exactly that support is at most q^(w-8),
where q=256. Assigning the same positive envelope to impossible
nonzero-row supports of size below nine only increases the sum.
Fix a set of a outer rows and define

    S_a(z) = sum_(h=0)^a binom(a,h) q^h R_h(z).

Summing the relaxed binary a-by-16 support matrices column by column
gives the complete a-class upper bound

    binom(L,a) / (q-1) * q^(-8a) * z^(-e)
      * [1,0] S_a(z)^16 [1,1]^T.

The state is not reset at region boundaries. The single factor q-1
removes scalar multiples of a complete source word; it must not be
applied separately to the a rows.

A separate small-a calculation retains each row's exact support
weight and MDS multiplicity. It uses a multivariate polynomial with
one support variable per selected row, so its cost is exponential in a.
The floating-point screens at 10% are:

| Active rows a | Coupled source envelope | Exact source-support weights |
|---:|---:|---:|
| 1 | 37.33 bits | 49.79 bits |
| 2 | 43.53 bits | about 43.76 bits |
| 3 | 20.89 bits | about 21.19 bits |
| 4 | -4.90 bits | about -4.55 bits |

The searches use different finite grids, so these rounded values are
diagnostic rather than certified optima. Beyond a=1, removing the
source envelope recovers less than half a bit in these tests. The
missing margin is not primarily a count of invalid outer supports.

## An exact obstruction to a complete word first moment

The following argument applies more generally to a fixed outer [b,k]
code and an s-symbol linear state. Suppose each of the b regions has
length L and ends with an aligned window of W symbols. Fix a set of
a outer rows. Require the regional permutations to place all a labels
in the final W positions of every region. This event has probability

    P_window = (binom(W,a)/binom(L,a))^b.

On this event, all outer words supported on those rows have inputs
confined to the b windows. Their message space has dimension ka.
Impose zero state at the first b-1 region ends. These impose at most
(b-1)s linear equations. Every resulting word is zero before each
region's window and has output weight at most bW. No equation is
needed for the discarded final state. Hence, if

    D = ka-(b-1)s > 0 and bW <= e,

the expected number of bad projective source lines is at least

    P_window * (q^D-1)/(q-1).

This fixes only one set of a rows. There is no multiplication by
binom(L,a), and therefore no overlap assumption about different row
sets. The argument holds for every diagonal and every inner refresh
on the specified routing event.

For the present parameters choose a=32 and W=816. Then D=226 and
bW=13056, below the forbidden threshold. Exact integer comparison gives

    (binom(816,32)/binom(8192,32))^16 * (256^226-1)/255 > 2^83.

The routing event itself has probability about 2^-1716.53. Thus the
large expectation does not imply that bad codes occur frequently.
It shows that a complete union bound over field-valued projective
messages cannot certify this ensemble, even with exact transfer
matrices and exact outer multiplicities. Increasing a fixed state size
only moves the onset of this phenomenon; it does not repair the
underlying count.

## A rank-potential route

The companion's weighted-scan proof instead fixes source and output
support paths and counts ranks of the required zero equations. Let V
be the source space shortened to a prescribed outer support, with
dimension D_V. As coefficients are exposed, let W be the span of
the zero equations in V*. A branch remains alive only if ker(W)
contains a message satisfying all required nonzero conditions already
defined on that branch. For tau>=1, assign potential tau^(-rank W)
to an alive branch, and zero to a dead branch. A nonzero witness
requires rank W<=D_V-1. This replaces a q^(D_V-1) message count by
a tunable tau^(D_V-1) rank factor.

The fresh input diagonal supplies a useful local lemma. Fix nonzero
required input forms f_1,...,f_w in the current quotient, and a
coefficient matrix B of rank r. These objects, and any prescribed
output support, must be fixed before the fresh diagonal is exposed.
Selecting the output support after seeing the diagonal requires a
separate union over support paths. Expose independent d_j from F_q*.
For the equations

    B diag(d_1,...,d_w) (f_1,...,f_w)^T = 0,

the expected alive potential is multiplied by at most g^r, where

    g = 1/(q-1) + (1-1/(q-1))/tau.

Choose r independent columns of B and row-reduce them to identity.
Expose the nonpivot coefficients first. Each pivot equation then has
the form d_j f_j+ell_j=0, where ell_j is fixed before d_j is exposed.
At most one d_j makes this equation redundant. Every other choice
increases the rank or kills the required-nonzero branch. If an earlier
pivot makes f_j zero in the quotient, that branch is already dead.
Multiplying these conditional bounds proves the lemma without assuming
the f_j are independent.

For an active-state return to zero, write X=diag(d)f. The recurrence
forces R=C X and therefore Y=(I+A C)X. If Z is the zero-output set,
first impose

    (I+A C)_(Z,S) diag(d) f = 0.

These equations do not depend on the current fresh M. Apply the
diagonal lemma first, then impose M H=C X, where H contains the two
incoming state forms in the resulting quotient. This exposure order
avoids reusing the same randomness for separate rank charges.
The required nonzero incoming state rules out rank(H)=0 on an alive
return branch.

For uniform SL(2,q), a fixed rank-one H reaches a specified complete
target with probability at most 1/(q^2-1). Its return-potential factor
is therefore at most

    tau^-1 + (1-tau^-1)/(q^2-1).

For rank-two H, the probability that rank(M H+T)<=1 is at most
1/(q-1). Union-bound over the q+1 left projective lines: for each
line, its product with M is uniform over the nonzero two-vectors.
The probability of rank zero is at most 1/(q(q^2-1)), since H is
injective on the matrix M. Thus a valid rank-two factor is

    tau^-2 + (tau^-1-tau^-2)/(q-1)
      + (1-tau^-1)/(q(q^2-1)).

This is a local return bound, not a complete transfer operator.
In particular, output constraints can change the incoming state rank.

The state rank must be retained. For an empty input step, rank zero
stays zero and emits zero. A rank-one state has the fixed-vector
output law G_0(z) and stays rank one on surviving branches. For a
rank-two state, the all-nonzero output pattern contributes z^t and
retains rank two. A prescribed single zero-output coordinate imposes
one independent equation and drops the state rank to one; summing
its possible positions contributes t*tau^-1*z^(t-1). Two zero
coordinates kill the required-nonzero state because two rows of A
are independent. Thus the empty-input upper transfer, with states
ordered by ranks zero, one, two, is

    [[1, 0,                         0],
     [0, G_0(z),                    0],
     [0, t*tau^-1*z^(t-1),         z^t]].

Collapsing ranks one and two would let a one-time rank loss recur at
every empty step and destroy this bound. The next subsection gives
conservative active-input bounds. Their composition with the complete
regional support count remains unfinished.

### Conservative active-input transitions

The following bounds give a starting operator, not an evaluated
whole-code certificate. Fix an input support S and a zero-output
set Z before the current independent coefficients are exposed. Put

    w=|S|, u=|Z\S|, r=|Z intersect S|, p=w-r.

Thus p active input lanes also have nonzero outputs. Let V_0,V_1,V_2
be nonnegative upper bounds for the remaining computation, conditional
on the next state-form rank. Each transition below also receives the
output marker z^(t-|Z|).

First consider a nonzero incoming state and a required nonzero next
state. Expose M, impose the inactive zero-output equations, and only
then expose the diagonal. If u>=2, the branch is dead. For u=0 or 1,
the inactive-lane factors and intermediate ranks are:

| Incoming rank | u | Factor alpha | Rank h' after these equations |
|---:|---:|---:|---:|
| 1 | 0 | 1-(t-w)/(q+1) | 1 |
| 1 | 1, prescribed lane | 1/(q+1) | 1 |
| 2 | 0 | 1 | 2 |
| 2 | 1, prescribed lane | tau^-1 | 1 |

For rank one, the inactive outputs are scalar multiples of the one
required nonzero state form, so these are exact probabilities for
their prescribed zero pattern. For rank two, one zero equation is
independent and reduces the state rank. Ignoring additional nonzero
conditions gives the stated upper factors.

Each of the r active zero-output equations has its own fresh diagonal
pivot and costs at most g. On a surviving branch, a canceled active
input is a linear function of the refreshed state. Only the p
uncanceled inputs can add new state forms. Consequently a safe
nonreturn bound is

    alpha * g^r * max { V_j : 1 <= j <= min(2,h'+p) }.

The maximum is over the unknown but actual next rank; summing the
continuations would introduce an avoidable overcount. This bound
does not assert that every listed rank is achievable. Further equations
may lower h' or kill the branch, which only decreases the potential.

For a required return to rank zero, use the different exposure order
proved above: diagonal-only output equations first, then the fresh M.
Let r_B=rank((I+A C)_(Z,S)), and call the two return factors above b_1
and b_2. A safe return factor for old rank one is g^(r_B)*b_1.
For old rank two it is

    min { g^(r_B)*b_1, b_2 }.

The second term needs a separate justification. If the diagonal
constraints leave the state rank two, the M equations cost b_2.
If they reduce it to one, they have already increased the zero-equation
rank by at least one; the combined cost is at most tau^-1*b_1<=b_2.
If they reduce it to zero, the required nonzero incoming state kills
the branch. The first term follows from b_2<=b_1 and the diagonal
lemma. Taking the smaller of these two independently valid upper
bounds avoids charging the same rank decrease twice. A product
g^(r_B)*b_2 based solely on the old rank is not justified.

A zero incoming state has Y=X, so its output support must equal S.
For w=1 the next rank is exactly one. For w>=2, a nonzero next state
can conservatively use max(V_1,V_2). Returns with w=1 or 2 are
impossible because any two columns of C are independent and active
inputs are required nonzero. For w>=3, a return imposes C_S diag(d)f=0
and the diagonal lemma gives the factor g^2 V_0. Empty input simply
stays at rank zero.

These bounds retain source-form dependencies through the alive-branch
condition. They do not treat input values as independent after a
quotient. The remaining task is to aggregate the rank transitions
over lane masks and regional occupancies, while preserving the
maximum-continuation semantics and complete source-support coverage.

## Reproduction

```sh
python -B research/workstreams/char2_field_spin/theory/subfield_certificate.py --certify
python -B research/workstreams/char2_field_spin/theory/subfield_certificate.py --coupled --screen 1 2 3 4
python -B research/workstreams/char2_field_spin/theory/subfield_certificate.py --exact-sparse --screen 1 2 3 4
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_subfield_certificate.py -v
```

Eleven tests replay the partial certificate, reject overstated coverage
and bounds, check interval operations against exact rationals, verify
regional normalization, and check the exact window obstruction.
GF4 tests also stress dependent input forms in the diagonal lemma,
enumerate all fixed two-by-two targets for the SL2 return bounds,
and check the empty-input rank operator. Independent rank tests are in
`lane_mixing/test_rank_lemmas_audit.py`; the whole-code counting audit is
in [SUBFIELD_AUDIT.md](SUBFIELD_AUDIT.md).
No encoding benchmark accompanies these proof calculations.
