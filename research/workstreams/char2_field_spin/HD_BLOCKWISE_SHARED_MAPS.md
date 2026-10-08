# Shared transforms for nonidentity feedthrough

The identity feedthrough has a deterministic quiet-subcode obstruction at
width 64/state 4. This proposal replaces it with eight-point MDS maps and shares
their transforms with state emission and feedback. It is a new construction,
not an optimization of the same code and not a distance certificate.

The exact matrix checks pass. The quiet-subcode obstruction does not transfer
unchanged, but partial source supports and cancellation histories remain open.
The later [native comparison](HD_BLOCKWISE_SHARED_NATIVE.md) records timings
for the beta64 variant; no complete distance certificate is claimed.

The construction below preserves the original translation-by8 checkpoint.
The next working choice uses translation by64: its output evaluation labels
lie outside the packet's input labels. This removes a specific neighbouring-
tile cancellation without changing the tile size or state dimension. The
exact comparison is recorded in [the joint-support analysis](HD_BLOCKWISE_SHARED_SUPPORT.md).

## Construction and coordinate order

Work over F=GF(65536), with the AES GF256 subfield used for fixed coefficients.
Each 128-bit record holds eight F-valued lanes. The same linear encoder acts
independently on those lanes. Distance refers to F-symbol positions, not bits.

Retain the GRS[256,128] outer and the single regional route. A packet has 64
input positions, partitioned into eight consecutive tiles of eight positions.
The state has four F coordinates. There is no preliminary accumulator, input
diagonal, puncture, or state flush.

For j=0,1,2, let

    L_j(X) = product_(v=0,...,2^j-1) (X+v),
    b_i(X) = product_(j: bit_j(i)=1) L_j(X)/L_j(2^j),   0<=i<8.

The integers denote the polynomial-basis labels in AES GF256. These normalized
novel-basis polynomials have degree i. Define the evaluation matrix

    E_alpha[u,i] = b_i(alpha XOR u),    0<=u,i<8.

The packet state Q contains the coefficients of a polynomial of degree below4
in this **global** basis. It is not a vector of evaluations on the old state
nodes. For a tile at alpha in {0,8,...,56}, compute

    z_alpha = E_alpha^(-1) X_alpha;
    h = sum_alpha (z_alpha[4],...,z_alpha[7]);
    Y_alpha = E_(alpha XOR8) (z_alpha + (Q,0,0,0,0));
    Q_next = M Q + h.

All eight tiles use the old Q. Update the state once after the packet. The
output from tile alpha stays in that tile's eight storage positions; the
shifted evaluation labels define a value-mixing map, not an extra route.
Initialize Q=0 and leave the final state free.

For the current proposal, each packet independently samples M=I+uv^T using
the exact law in [HD_BLOCKWISE_INNER.md](HD_BLOCKWISE_INNER.md): u is uniform
over F^4 excluding zero, and v is uniform conditioned on v^T u !=1.
Sampling and routing are independent. The matrix is invertible.

## Exact fixed maps

The tile computation has the packet form

    Y = F_64 X + A Q,        Q_next = M Q + C X.

Here F_64 is block diagonal. Each eight-point block is
E_(alpha XOR8) E_alpha^(-1), and the blocks are identical. Its first row is

    [52,172,158,17,209,32,250,29].

The block is an XOR circulant. The corresponding rows of A are the first four
columns of E_(alpha XOR8). The corresponding columns of C are the last four
rows of E_alpha^(-1). Thus feedback is the sum of the high coefficients already
computed by the inverse transforms. It needs no additional per-tile linear map.

Exact arithmetic gives C=H A^T, with the invertible matrix

    H = [245 29  9 1]
        [ 29  0  1 0]
        [  9  1  0 0]
        [  1  0  0 0].

In particular, C is not literally the old A^T. The output shift is included
in this displayed identity.

The four state polynomials have degree at most3. Products of any two have
degree at most6. Their sum over the 64-point additive space is zero, so

    C A = 0.

The matrix A generates an MDS[64,4,61] code. Since H is invertible, C checks
an MDS[64,60,5] code. These statements hold over F as well as the coefficient
subfield: the required minors remain nonzero under field extension.

One can also see why the feedback map is independent of a tile's coset.
The Lagrange polynomial for a node u in alpha+[0,8) is

    (L_3(X)+L_3(alpha))/((X+u)L_3').

Its coefficients of degrees4 through7 depend on u but not on L_3(alpha).
They are polynomials of degree at most3 in u. Conversion to the normalized
novel basis applies the same invertible triangular change of coordinates to
every tile. This proves the fixed relationship between C and low-degree
evaluation moments; the helper also verifies H entry by entry.

## Feedthrough and the state identities

The eight-point feedthrough evaluates a degree-below8 polynomial on a disjoint
coset. Its matrix is a row- and column-scaled Cauchy matrix. Every square minor
is nonzero. Consequently

    wt(x) + wt(F_8 x) >=9       for x!=0.

Applying the same translation twice is the identity in characteristic two.
Thus F_8^2=I. The exact state translation is

    T = [1 8 28 224]
        [0 1  0  28]
        [0 0  1   8]
        [0 0  0   1],

and the packet maps satisfy

    F_64 A = A T,     C F_64 = T C,     T^2=I.

Hence C Y=T C X. Given a fixed output stream, its inverse state recurrence
uses Q_next=M Q+T C Y, followed by X=F_64 Y+A T Q. This inversion does not
depend on the route. The full encoder is invertible before the outer-code
restriction.

Equivalently, changing the inverse state to TQ gives the same packet form with
mixer TMT. Its distribution equals the original mixer distribution, because
u maps to Tu and v maps to T^(-T)v. A fixed sampled encoder need not be its
own inverse; this is an equality of ensembles.

## What happens to the quiet-subcode obstruction

Imposing CX=0 still keeps the state zero. It now emits F_64 X, rather than X.
An input coordinate can spread to eight output coordinates within its tile.

One possible replacement obstruction is to impose both CX=0 and
(F_64-I)X=0. Exact ranks are

    rank(F_64-I)=32,
    rank([C; F_64-I])=34.

For P packets, fixing the output on this subspace costs at most 34(P-1)+32
equations: the last packet needs no feedback constraint. The earlier row-set
dimension argument would asymptotically require a fraction 2*34/64=1.0625 of
all source rows at rate 1/2. Thus it guarantees no admissible row set and no
low-distance word. This does not exclude a more selective bad subspace.

The involution alone is not a sufficient obstruction. Its fixed space has
dimension32, but imposing it uses substantial source dimension. Within one
eight-point tile, the fixed polynomials are polynomials in X^2+8X of degree
below4. Their evaluations form an MDS[8,4,5] code, because X^2+8X is injective
on [0,8). The fixed space is not a collection of sparse coordinate vectors.

A correction of the form F=I+B C would not help: it is the identity on ker C.
The present F does not have that form.

## Outside-packet evaluation labels

Translation by8 pairs adjacent input cosets: one tile's output evaluation
labels are another tile's input labels. The joint-support analysis finds
duplicated state constraints and an exact zero-feedback swap on these pairs.
This is a local cancellation mechanism, not a complete distance attack.

Instead, replace alpha XOR8 by alpha XOR64 in the shared circuit. The eight
output cosets now cover64,...,127, disjoint from all input labels0,...,63.
Every tile still has branch number9. The first row of its feedthrough is

    [16,105,170,110,165,217,144,80].

The same identities F^2=I, FA=AT, CF=TC, C=HA^T, and CA=0 hold, with the
new translation matrix T and corresponding H. All12,869 tile minors and
the displayed matrix identities have been checked exactly. This changes
the fixed map, not the regional route or memory positions.

There is a wider cost-compatible choice. For each tile, choose any aligned
output coset delta distinct from its input coset alpha. Require the eight
input cosets to be distinct, and likewise the eight output cosets. Set

    F_j = E_delta E_alpha^(-1),
    A_j = E_delta[:,0:4],
    C_j = E_alpha^(-1)[4:8,:].

Each C_j A_j is zero individually: translating a degree-below4 polynomial
still gives degree below4. Also F_j^2=I. Swapping the input and output cosets
gives inverse-family maps

    A'_j = F_j A_j = E_alpha[:,0:4],
    C'_j = C_j F_j = E_delta^(-1)[4:8,:].

Consequently C F A=0 and the inverse state recurrence remains
Q_next=M Q+C F Y. The stacked A and A' are MDS evaluation maps; the stacked
C and C' are their fixed degree-below4 moment checks on distinct nodes.
Independent tile shifts need not admit one common T or satisfy C=HA^T.
The inverse belongs to the swapped-coset family, not necessarily to the
same sampled construction. This generalization is an algebraic interface,
not a new random-label distribution or a distance claim.

Translation by64 is the current simple choice. Arbitrary distinct cosets
remain available if the joint-rank analysis identifies a useful change.

## The support statistic that measures the benefit

For a tile with w actual nonzero input symbols, define

    beta(0)=0,       beta(w)=9-w for 1<=w<=8.

The MDS branch bound gives wt(F_8 x)>=beta(w). Singletons therefore emit eight
symbols without fresh multiplications from the full field. Cancellations in
multiply occupied tiles are allowed by this bound.

As a diagnostic, independent Bernoulli occupancy with probability p gives

    E[beta(W)]/g = ((g+1)/g)*(1-(1-p)^g)-p,

for a g-point MDS feedthrough. At the minimum complete-row quiet threshold
p=2s/t=1/8, g=8 gives approximately0.61344. For comparison, g=4 gives0.39227.
These are means under a diagnostic occupancy model, not distance bounds.

There is also an exact without-replacement interface. Fix a region of R
positions, with a fixed set of a genuinely nonzero positions before its
uniform permutation. Let n_w count g-point tiles containing w selected
positions. Then

    Pr[(n_0,...,n_g)]
      = ((R/g)! / product_w n_w!)
        * product_w C(g,w)^(n_w) / C(R,a),

provided sum n_w=R/g and sum w*n_w=a. Otherwise the probability is zero.
Equivalently, the branch-weight generating polynomial is

    E[z^(sum beta(W_tile))]
      = [x^a](sum_(w=0)^g C(g,w)x^w z^beta(w))^(R/g) / C(R,a).

For a fixed source word, apply this formula separately to each region's actual
support count. The regional permutations remain independent. A row allowed to
carry a nonzero outer word need not be nonzero in every coordinate; substituting
the number of selected rows for every regional support count is invalid.

The exact local value count is also available. Since (x,F_8x) forms an
MDS[16,8,9] code, prescribed input and output supports of sizes w and u have
the usual MDS prescribed-support count at total weight w+u. This count is
implemented in the helper. It is not a conditional independence statement
about routed outer symbols.

## Remaining proof work and operation scope

For each fixed packet input X, set u=wt(F_64 X) and c=CX. The previous
three-component envelope applies with shift F_64 X, because im A remains
MDS and C F_64 X=0 exactly when c=0. A global proof must retain enough joint
information about raw support, mixed support, outer constraints, and the
syndrome to compose this local statement. Replacing the mixed support by an
independent random support would change the construction.

The next useful enumeration is therefore the eight-tile input/output support
profile together with the rank of the summed four feedback coefficients.
The outer MDS equations must be charged jointly with those constraints. The
rare-route first-moment issue also remains; improving local distance alone
does not justify an uncapped whole-code first moment.

The proposed arithmetic uses one inverse FFT8 and one forward FFT8 per tile.
State injection uses four coefficient XORs per tile; feedback accumulates
four coefficients per tile. The existing rank-one mixer adds eight full-field
multiplications per 64-position packet. Fixed transform coefficients are in
GF256. These are algebraic operation counts, not measured instruction counts
or evidence that the end-to-end encoder meets a timing target.

### A possible later mixer change

One separate candidate uses M(a)=Theta+aI, where Theta represents
multiplication by a root of an irreducible degree-s polynomial over F and
a is sampled from F. Each factor is invertible for s>1. A companion basis
uses a shift, fixed feedback, and s full-field scalar products per packet.

For a **zero-input** window of s fresh factors and a fixed nonzero entering
state, every point probability is at most s!/q^s, where q=|F|. The product
polynomial product_i(X+a_i) is monic of degree s. Its remainder modulo the
fixed monic irreducible polynomial determines that polynomial uniquely, and
an ordered root tuple has at most s! permutations. At s=4 this gives 24/q^4.

This statement does not apply across arbitrary input syndromes or after
conditioning on output cancellations. Those dependencies need a different
proof interface. No specific irreducible polynomial, implementation, or
whole-code claim is selected by this observation.

## Reproduction

[theory/hd_blockwise_shared_maps.py](theory/hd_blockwise_shared_maps.py) defines
the exact matrices and structural receipts. Its seven tests check all12,869
square minors of F8 at both shifts8 and64, all68 input/state basis cases of
the saved circuit, the displayed identities and ranks, exact small regional
histograms, and the MDS prescribed-support count. They also check several
outside-packet translations and the swapped-coset inverse family. All pass.

    python -B -m unittest hd_blockwise_shared_maps_tests -v
    python -B hd_blockwise_shared_maps.py

The prior construction and its proofs remain unchanged.
