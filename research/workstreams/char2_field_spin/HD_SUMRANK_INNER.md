# Sum-rank-separated inner maps: proof-side candidate

This is a distinct candidate, not an implementation change
for the existing encoder. It keeps the randomized RS[256,128] outer,
one regional routing stage, 64-symbol packets, and four GF65536 state
coordinates. It changes the fixed emission and feedback maps. Its complete
quiet-word class is now certified below 2^-93.9900 at 40% symbol distance.
The full distance certificate remains open because nonquiet histories remain.
The separate [sixteen-state feedback refinement](HD_SUMRANK_FEEDBACK_CERTIFICATE.md)
now has a complete 40%-distance certificate with more than 59.87 bits of
margin. That theorem does not close the four-state construction here.

## Construction

Write B=GF16 and F=GF65536, with F of degree four over B. Choose an additive
subspace H of B of size eight and a shift delta outside H. Let E0 and E1
evaluate degree-below-eight polynomials on H and delta+H. Define

    F0=E1 E0^-1,
    A0=the first four columns of E1,
    C0=the last four rows of E0^-1.

Thus F0 is a Cauchy map up to nonzero coordinate scalings. Translation by
delta preserves the degree-below-four subspace, so C0 A0=0 and
C0 F0^-1 A0=0.

Choose a B-basis beta_0,...,beta_3 of F and eight distinct nonzero a_i in B.
Set R_i[j,l]=a_i^j beta_l^(16^j), for j,l=0,...,3. The concatenation
[R_0|...|R_7] is a dimension-four linearized Reed--Solomon generator.
The nonzero a_i have distinct norms a_i^4 and hence distinct Frobenius
conjugacy classes. Theorem 4 and Section 4.2 of
[Martinez-Penas, Skew and linearized Reed--Solomon codes](https://arxiv.org/abs/1710.03109)
give maximum sum-rank distance for this choice.

In particular, for arbitrary B-subspaces W_i of B^4,

    rank_F [R_i basis(W_i)]_i = min(4, sum_i dim_B W_i).

To derive this consequence, extend each local subspace basis to an
invertible B-matrix. The transformed generator is MDS in Hamming metric:
a nonzero word has sum-rank weight at least 29, so it also has Hamming
weight at least 29. Any four transformed columns are independent.

Use the packet maps

    F=diag(F0,...,F0),
    C=[R_0 C0 | ... | R_7 C0],
    A=vertical_stack(A0 R_0^T,...,A0 R_7^T).

For the principal candidate the recursion remains Y=FX+AQ,
Q_next=M Q+CX with the existing independent rank-one mixers.
Both CA=0 and C F^-1 A=0 hold term by term. Therefore F+lambda AC is
invertible for every scalar lambda: F^-1 AC has square zero. The existing
weighted-state column argument requires this identity, not a shared
four-dimensional translation matrix.

## Exact local rank formulas

Let a_i count potential input positions in tile i, and let b_i count
output positions selected to be zero. Put r0=sum_i min(a_i,b_i).

For a quiet packet the remaining polynomial space in tile i has dimension
d_i=max(a_i-b_i,0). Its highest-four-coefficient image is a B-defined
subspace of dimension min(4,d_i). Sum-rank separation consequently gives

    rank [F_(zeros,U); C_U] = r0 + min(4,sum_i d_i).

For an arbitrary state, eliminate the a_i local input variables from
the selected output constraints. The remaining constraints on the local
four state coefficients have rank min(4,max(b_i-a_i,0)) over B.
This is the extended Cauchy--Vandermonde rank: multiply the rational
functions with poles at the potential input nodes by their common
denominator to obtain all polynomials of degree below a_i+4.
The transposed R_i blocks again give sum-rank separation. Hence

    rank [F_(zeros,U), A_zeros]
        = r0 + min(4,sum_i max(b_i-a_i,0)).

These are universal formulas, not rank estimates inferred from samples.
Expanding z^wt(Y) over selected output zeros gives exact rational quiet,
free, and uniform-nonzero-state moments. All depend only on the eight
occupancies, not the particular potential slots inside a tile.

For a fixed nonzero state, an output-zero mask is impossible when
sum_i min(4,max(b_i-a_i,0))>=4. Otherwise its probability over uniform X
is either zero or q^(-r0). Summing these positive mask bounds gives a
pointwise live-state moment, tighter than multiplying V by q^4-1.

The same formula gives the exact minimum weight over X and nonzero state:
maximize sum_i b_i subject to that rank cost being at most three.
With one potential input in every tile, each nonzero state forces at least
53 nonzero outputs, even though no tile is untouched.

## Complete quiet-class certificate

The [directed helper](theory/hd_sumrank_quiet_certificate.py) covers all
8192 possible positive active-row counts. The old sparse quiet bound
transfers because C is MDS[64,60,5]. The inverse-output dense bound transfers
because the encoder remains invertible on its full input space. Eleven
new intervals cover the middle counts using the exact moments above.
Convexity and directed endpoint evaluation cover every integer in each
interval. No gap is inferred from sampled row counts.

At 192-bit Arb precision, the combined logarithmic failure bound is

    [-93.9900198408693736724657884346656242590532030707230 +/- 4.35e-50].

This bounds the existence of any nonzero quiet message of output weight
at most 838860, including all outer-row choices. The final packet's
successor state is unrestricted. These quiet bounds do not depend on M_t.
They also apply to the independent invertible
[feedback refresh](HD_FEEDBACK_REFRESH.md), because R_t CX=0 if and only
if CX=0 when the previous state is zero.

The same property permits an exact input-weight mark. For any such moment
H_a, marking by v^wt(X) replaces it by

    sum_(b_i<=a_i) product_i[binom(a_i,b_i) v^b_i (1-v)^(a_i-b_i)
                              q^(-(a_i-b_i))] H_b.

This is an algebraic inclusion/exclusion identity. The helper evaluates
its signed integer numerators exactly; it never treats negative terms as
probabilities. Final moments remain nonnegative.

## Status and cost considerations

Exact finite-field regressions check the field embedding, nilpotence,
local Cauchy minors, base-field subspace transforms, and concrete mask ranks.
The mathematical guarantee comes from the rank argument above, not from
the sampled checks. Initial global screens improve some dense classes but
do not close intermediate occupancies. Local rank additivity alone is
therefore insufficient with the current packetwise fractional bound.

The factorization R_i=diag(1,a_i,a_i^2,a_i^3) M uses one shared Moore matrix.
It permits feedback and emission to share their extension-field products:
16 fixed F products for feedback and 16 for emission per packet, with
the remaining coefficient multiplications in B. This is an operation
count, not a native performance result. No code has been promoted.

A separately labeled `uniform_pre` diagnostic replaces M by a uniform
transitive state action before emission: H=M Q, Y=FX+AH, Q_next=H+CX.
It has a two-state transfer because H is uniform nonzero whenever Q is
nonzero. Such an action can be multiplication by a uniform nonzero element
of GF(65536^4), viewed F-linearly. This is a different ensemble and is only
a proof control; its bounds do not apply to rank-one refresh.

Sources: [exact maps](theory/hd_sumrank_maps.py),
[map regressions](theory/hd_sumrank_maps_tests.py), and
[moment/composition screens](theory/hd_sumrank_transfer.py).
