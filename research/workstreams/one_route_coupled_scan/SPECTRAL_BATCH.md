# Evaluation-based batches with the same inverse distribution

This inner replaces the Toeplitz product by two fixed transforms and a
random pointwise product. It retains the exact distribution used by the
[closed finite certificate](STRATIFIED_CLOSED_CERTIFICATE.md). The argument
below certifies the replacement at memory 128 and batch size 64; it makes
no timing claim or claim about the older small-coefficient constructions.

The change is structural: state entries are interpreted directly as
coefficients in the transform's basis. They are not converted from
monomial coefficients merely to reproduce the old Toeplitz map.

## Construction

Fix E=GF256 as a subfield of F=GF(2^128). Fix an ordering of all 256
elements of E. Also fix polynomials P_0,...,P_255 in E[X], with
deg(P_i)=i and nonzero leading coefficients. A degree-graded novel basis
used by an additive transform is permitted; normalization need not make
the polynomials monic. All choices in this paragraph are deterministic.

Define two F-linear maps:

    G : F^128 -> F^256,
    (G S)_a = sum_(i=0)^127 S_i P_i(a),       a in E;

    H : F^256 -> F^64.

To compute H(v), interpolate the unique polynomial Q_v of degree below
256 with Q_v(a)=v_a for every a in E. Expand Q_v in the basis P_i and
return its coefficients with indices 192 through 255, in that order.
Thus G and H have coefficients in E, although their inputs lie in F.

Fix a set J subset E of active evaluation points, with |J|>=191.
The first version uses J=E. For each batch t, independently sample r_(t,a)
from all of F for every a in J, including zero. Set r_(t,a)=0 for a outside J. These coins
are independent of the outer and route. Define

    T_t := H diag(r_t) G.                                  (1)

Setup samples every r_t once and fixes it for all messages. It does not
sample a new diagonal during each encoding.

Retain the certified outer and route: K=2^20, N=2K, fixed RS[64,32] outer
blocks over E extended to F, and 64 independent uniform regional orders.
Each region receives one coordinate from every outer block. This is the
same single structured global permutation as in the finite certificate.

Divide its routed input into batches X_t in F^64. Start with two zero
output batches. Let S_t in F^128 concatenate Y_(t-2) and Y_(t-1), oldest
first, and compute

    Y_t := X_t + T_t S_t.                                 (2)

Append Y_t and discard the oldest batch. There is no feedback within a
batch and no terminal-state constraint. The transforms operate only on
this local state; they add no global routing stage.

For every setup, (2) is an invertible F-linear map: reconstruct S_t from
the preceding outputs and recover X_t=Y_t-T_t S_t. Hence the full encoder
has rate exactly one half.

## Every 64 columns of H are independent

First extract the highest 64 monomial coefficients instead of the highest
64 P_i coefficients. The two extraction maps have the same row space.
Indeed, monomial coefficients of degrees 192 through 255 depend only on
P_i coefficients with i>=192. Their change-of-basis matrix is triangular
with nonzero diagonal, and is therefore invertible.

It remains to compute the monomial extraction matrix. The polynomial

    Z(X) := product_(a in E) (X-a) = X^256-X

has derivative Z'(a)=1 in characteristic two. Its Lagrange polynomial at
a is Z(X)/(X-a). For 0<=j<=63, the coefficient of X^(255-j) in that
polynomial is a^j. One direct identity is

    Z(X)/(X-a) = sum_(j=0)^255 a^j X^(255-j) + 1,

where a^0=1, also when a=0. The extra constant does not affect the
coefficients under consideration.

Consequently the monomial extraction matrix is a Vandermonde matrix
with its rows reversed. Any 64 columns have nonzero determinant because
their evaluation points are distinct. The same holds over F, and the
invertible row transformation proves the claim for H.

## Exact joint uniformity

Fix any nonzero S in F^128 before sampling the current diagonal. Its
polynomial

    p_S(X) := sum_(i=0)^127 S_i P_i(X)

is nonzero and has degree at most 127. It vanishes at at most 127 points
of E. Therefore at least |J|-127 entries of G S inside J are nonzero.

If |J|>=191, at least 64 active coordinates remain. Select any 64 of
them. Their columns of H are independent, and scaling those columns by
their nonzero entries of G S preserves independence. Thus the map

    (r_a)_(a in J) -> H diag(G S) r

has rank 64 over F. Its input is uniform on F^J, so its output is exactly
uniform on F^64. In particular, for every z in F^64,

    Pr_(r_t)[T_t S=z] = |F|^-64.                          (3)

This is joint uniformity, not just uniformity of each coordinate. For
S=0, the output is deterministically zero. Fresh independent diagonals
give independent feedback vectors at different batches when their
states have been fixed in advance.

## Transfer of the finite all-message certificate

Fix a candidate output y and a realization of the route. Its preceding
output states S_t(y) are now fixed, independently of the inner setup.
Its inverse satisfies

    X_t = Y_t - T_t S_t(y).

Equation (3) makes X_t jointly uniform on F^64 whenever S_t(y) is
nonzero. When S_t(y)=0, X_t=Y_t. Independence across batches therefore
gives exactly the affine-coordinate distribution proved for the
[Toeplitz inner](TOEPLITZ_BATCH_BASELINE.md).

Let A contain all scalar coordinates of batches with nonzero preceding
state. Let U also include every coordinate of each nonzero output batch.
Write F^A for the vectors supported in A. If C is the fixed routed outer and

    D(A) := dim_F(C intersect F^A),

then the inverse belongs to C with probability at most |F|^(D(A)-|A|).
The same activation intervals, forced zero tails, and D(U)=0 test apply.
The conditional counting bound and its route-averaged fractional moment
are consequently unchanged for every candidate output class.

This transfers the completed finite bound, including its full coverage:

    Pr_setup[d_min <= floor(0.35 N)] < 2^-40,             (4)

at K=2^20, rate 1/2, outer RS[64,32], memory 128, and batch size 64.
Distance in (4) counts nonzero F-symbols. The probability is over the
ideal regional permutations and the independently uniform diagonals.
It concerns the existence of any bad nonzero message after setup.

The Toeplitz and evaluation-based maps need not have the same matrix
distribution or the same joint action on several candidate outputs.
The certificate uses the fixed-output inverse distribution followed by
counting, so equality of those stronger distributions is unnecessary.
No new numerical cover is needed for this exact distributional transfer.

## Using 192 or 191 active points

The construction and proof permit any fixed J of size 192 or 191, with
all other diagonal entries zero. They retain the full 256-point maps G
and H. In particular, H still interpolates a vector on all 256 points;
inactive pointwise products are zero, not omitted interpolation data.

For 192 active points, a nonzero state leaves at least 65 nonzero active
evaluations. For 191, it leaves at least 64. Both satisfy the same rank
argument and inherit (4) without a probability loss.

The value 191 is minimal for this fixed evaluation-and-diagonal template.
If |J|<=190, a nonzero polynomial of degree at most 127 can vanish on
min(127,|J|) chosen active points. At most 63 active evaluations then
remain. The diagonal-to-output map has rank at most 63 and cannot be
uniform on F^64. This is not a lower bound on every possible inner design.

Pointwise full-field products per output symbol are 4, 3, and 191/64 for
|J|=256,192,191. Fixed transforms and state updates remain additional work.
Masked points may permit transform pruning, but skipping their pointwise
products alone does not establish a cheaper transform.

### Exact three-coset kernel

Choose a Cantor basis beta_0,...,beta_7 for E, with beta_0=1 and
beta_i^2+beta_i=beta_(i-1). Order points by their binary coordinates
in this basis. Put V_6=span(beta_0,...,beta_5) and
s_6(X)=product_(v in V_6)(X-v). Use normalized novel polynomials P_i.
Their normalizations are one in this basis.

Write the state polynomial as A(X)+s_6(X)B(X), with A and B expressed
in P_0,...,P_63. Choose J to be the first three cosets of V_6.
On those cosets, s_6 is respectively 0, 1, and beta_1. Thus the forward
transform consists of three offset-specific 64-point transforms on the
coefficient vectors A, A+B, and A+beta_1 B.

After the 192 independent pointwise products, apply the corresponding
three inverse transforms. XOR their coefficient vectors to obtain H's
64 outputs. To see this, the high half of each inverse additive-FFT
butterfly is the XOR of its two inputs, independently of its twiddle.
The highest quarter of the full 256-point inverse is therefore the XOR
of all four coset inverse transforms; the inactive fourth coset is zero.
This computes (1) exactly, rather than approximating a truncated inverse.

### A limited arithmetic lower bound

Consider an F-bilinear implementation whose random map has the form

    T_r(S) = T_0 S + sum_(j=1)^R v_j ell_j(S) m_j(r),

where T_0 in F^(b by L), vectors v_j in F^b, and linear forms ell_j on F^L
are fixed.
The implementation uses R variable scalar F-products; fixed linear
operations are not charged in this model. If T_r(S) is uniform on F^b
for every nonzero S, then R>=L+b-1.

Indeed, the ell_j must span the dual of F^L. Otherwise a nonzero common
kernel vector gives deterministic feedback. Choose L-1 independent
forms and a nonzero S annihilating them. The random part of T_r(S) now
lies in a span of at most R-L+1 output vectors. Uniformity requires this
span to have dimension b, proving the bound. For L=128,b=64 this gives
R>=191.

This is not a CPU-instruction lower bound. It does not rule out a better
field representation, a bit-level algorithm, or a weaker proof interface.

## Implementation and scope

The intended kernel evaluates the state directly from P_i coefficients,
multiplies active evaluations by the precomputed r_t, and extracts the
top 64 P_i coefficients after interpolation. No monomial-basis conversion
is part of the construction. Any pruned transform must compute these
same maps exactly.

This proof still requires independently uniform full-F diagonal entries.
Subfield entries, shared diagonals across batches, and a compressed scalar
hash are different distributions and do not inherit (3). The earlier
cheap-lane families and the memory-256/45% profile are not certified by
this particular 128-to-64 construction.

The native implementation now passes independent polynomial oracles for
G, H, and (2), including the pruned three-coset kernel. See the
[native record](NATIVE_BASELINE.md) for serial measurements and limitations.
Neither (3) nor (4) supplies an encoding time.
