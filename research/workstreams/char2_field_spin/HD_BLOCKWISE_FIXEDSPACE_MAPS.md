# Sharing the beta64 residual with state feedback

The current shared circuit needs six independent local forms to produce its
rank-four feedthrough residual and four feedback values. This proposal changes
the state maps so that both operations use the same four forms. The beta64
feedthrough itself is unchanged. This is an exact construction proposal, not
a whole-code distance certificate or a performance result.

## Construction

Work over GF65536, with the AES GF256 subfield used for fixed coefficients.
The packet has 64 input coordinates, arranged in eight consecutive tiles.
The state Q has four coordinates and starts at zero. Set

    z(x)=x^2+64*x,             A[x,j]=z(x)^j,  j=0,1,2,3,

for input labels x=0,...,63. Multiplication and addition in these formulas
are field operations, not integer arithmetic. Let

    K = [17 114   0 240]
        [114  0 240   0]
        [ 0 240   0   0]
        [240  0   0   0].

For tile i, write A_i for its eight rows and set C_i=K A_i^T.
Let C=[C_0 ... C_7]. For each packet, compute

    c_i=C_i X_i,
    Y_i=X_i+A_i(Q+c_i),
    Q_next=M Q+sum_i c_i.

Here M is the independently sampled invertible rank-one mixer of
[the transfer note](HD_BLOCKWISE_SHARED_TRANSFER.md). The outputs use the
incoming Q, not the refreshed state. In particular, each output tile uses
its own c_i; replacing all c_i by their sum would change the construction.

## Exact identities

The map z is additive. Its kernel in GF256 is {0,64}, so it is injective on
the 64 packet labels. Thus A is a Vandermonde generator for an MDS
[64,4,61] code. K is invertible, so C has the same column dependencies as
A^T: every four columns are independent, and ker C has distance five.
Both facts extend to the GF65536 coefficient field.

The beta64 tile feedthrough F_i satisfies

    F_i-I=A_i K A_i^T=A_i C_i.

One exact derivation writes its entries as g(z(x)+z(y)), where

    g(Z)=17+114 Z+240 Z^3.

Expanding this cubic gives K. The translation-circulant tile matrix is
independent of the tile offset; adding the same z-offset to both arguments
does not change their sum.

Each tile's eight z-nodes form an additive three-dimensional affine space.
Their power sums through degree six vanish. Therefore

    C_i A_i=0,  F_i A_i=A_i,  C_i F_i=C_i.

Globally, CA=0, FA=A, and CF=C. The packet remains invertible for every M:
given Y and old Q, recover X=F Y+A Q and use CY=CX for the state update.
Starting from zero, the full inner encoder is again an involution for each
fixed mixer sequence. No independent-amplitude assumption is used.

## The local support tradeoff

The old local joint [16,12,5] description does not survive this change.
For a singleton x in tile i, take Q=C_i x. The tile emits y=x. This gives
a joint input/output word of weight two.

More generally, the tile's pairs (x,y) are exactly

    C_i(x+y)=0.

Indeed, ker C_i has dimension four and equals image A_i. With free state,
any difference in that kernel can be obtained. For input support S and
required output zeros Z, the exact constraint rank is

    |S intersect Z| + min(4, |Z minus S|).

To prove the formula, replace the free local state by Q+C_i x. The
constraints become x_Z+(A_i)_Z Q=0. Coordinates in S intersect Z supply
distinct input pivots. The remaining state rows have Vandermonde rank.

This state substitution is local. It does not create eight independently
adjustable states in a packet. A whole-packet proof must retain the single
shared Q and all eight feedback dependencies.

## A surviving whole-packet bound

For any nonzero incoming Q, the packet satisfies

    wt(X)+wt(Y) >= 37.

First, each tile has combined weight at least two. If a tile has combined
weight at most four, C_i(x+y)=0 and the four-column independence imply
x=y. Such a tile has input weight at most two, and Q=C_i x.

Two distinct tiles cannot both have combined weight at most four. Their
two representations of Q would yield a dependence among at most four
distinct columns of the global C. Consequently, at most one tile has
combined weight below five. The minimum possible total is 2+7*5=37.

For Q=0, every nonzero tile retains the unchanged feedthrough branch
bound wt(x)+wt(F_i x)>=9. These are deterministic packet statements, not
global relative-distance bounds.

## Quiet words and the remaining proof question

A quiet packet only imposes sum_i c_i=0. It does not impose c_i=0 in every
tile, and it does not force Y=X. The helper gives an exact five-symbol
input with zero global feedback and nonzero local feedback; its output
differs from its input.

The stronger condition (F-I)X=0 does imply CX=0. The stacked map [C;F-I]
has rank 32, versus 34 for the previous state maps. Thus the sparse-row
fixed-subcode dimension argument no longer gives a useful upper bound
below rate-one support. This does not prove that all other quiet traps
are absent. The existing sparse quiet routing gate still applies because
ker C retains distance five.

The rank-one mixer's deterministic joint-image return mechanism also
remains: a diagonal subspace {(h,h)} of dimension at least two always
contains a nonzero returning image vector. The new fixed maps do not
repair that issue automatically.

The next proof step is to combine the exact local rank formula with the
shared-state/feedback image, outer equations, and regional routing.
The saved rank-transfer operator applies to that joint image after output
constraints. It cannot be replaced by independent tile amplitudes or by
the old local-GRS ranks.

The helper and five tests check exact matrices, local constraint ranks,
singleton cancellations, and a quiet/nonfixed example. The Vandermonde
argument and the weight-37 argument are proofs, not conclusions from
sampled tests. No native kernel or timing is included in this checkpoint.
