# Exploratory nonbinary syndrome mixer

This proposal removes all quiet vectors of weights four and five using eight
syndrome coordinates. It is an algebraic fallback, not an implemented mixer or
a distance certificate. The existing H* construction and all its receipts remain
unchanged. In particular, the current arity-five proof failure does not establish
that quiet vectors are its dominant obstruction.

## Map and scope

Work over the AES field F of order 256. The constants 1, 2, and 4 denote the
corresponding byte-basis elements, which are linearly independent over GF(2).
Index the sixteen coordinates by b=(b0,b1,b2,b3) in GF(2)^4. Let V have rows
`(1,b0,b1,b2,b3)` and retain the existing invertible coefficient matrix B in

`H* = I + V B V^T`.

Define three additional column vectors by evaluation on these labels:

\[
g_1(b)=b_0b_1,\qquad g_2(b)=b_0b_2,\qquad
q(b)=b_1b_2+2b_1b_3+4b_2b_3,
\]
\[
f(b)=b_0b_3+q(b),\qquad W=[V,g_1,g_2,f].
\]

The proposed column mixer is

\[
H'=H_*+g_1g_1^T+g_2g_2^T+ff^T
   =I+W\operatorname{diag}(B,I_3)W^T.
\]

The intended outer column map is `D2 H' D1`, with the same independent uniform
nonzero diagonal entries as the existing construction. An input permutation
schedule, if used, remains deterministic and must be specified separately.
The RS rows, group size, regional route, and inner encoder need not change.
This is a new ensemble; existing numerical receipts do not transfer by relabelling.

Call a vector x quiet when H'x=x. The eight columns of W are independent:
their multilinear monomials have distinct nonzero coefficients, and f contains
the monomial b0*b3 absent from the other seven columns. Therefore

\[
H'x=x\quad\Longleftrightarrow\quad W^Tx=0.
\]

The quiet code has dimension eight. The following results concern this code,
not the support of every nonquiet output.

## Invertibility

Summing the products of the displayed evaluation functions over all sixteen
binary labels gives

\[
V^TV=0,\qquad V^Tg_1=V^Tg_2=V^Tf=0,
\]
\[
G:=\begin{pmatrix}g_1^T\\g_2^T\\f^T\end{pmatrix}
       [g_1,g_2,f]
 =\begin{pmatrix}0&0&4\\0&0&2\\4&2&0\end{pmatrix}.
\]

For example, the only term of g1*q whose sum is nonzero is
`4*b0*b1*b2*b3`, so g1^T f=4. The analogous term gives g2^T f=2.
The diagonal entries vanish because the sums of g1, g2, and f vanish;
in characteristic two, the sum of squares is the square of the sum.

The determinant identity `det(I+UV)=det(I+VU)` now gives

\[
\det H'=\det(I_3+G)=1+4^2+2^2=\mathtt{0x15}\ne0.
\]

The last equality uses AES byte arithmetic: 4^2=0x10 and 2^2=0x04.
Unlike the binary star variants, H' is not an involution. Indeed, its
correction has square `W diag(B,I3) (W^T W) diag(B,I3) W^T`, which is nonzero
because W is injective, W^T is surjective, and G is nonzero.

## No nonzero quiet vector has weight at most five

Fix x in ker(W^T). Split its coordinates into the two b0 halves. The first
seven checks, namely V, g1, and g2, force each half separately to satisfy
the three checks `(1,b1,b2)`. They also impose one global b3 check.

Within one half, a nonzero support of size one is impossible. A support of
size two must use the same `(b1,b2)` value twice, hence the two b3 values,
with equal amplitudes. A support of size three is impossible with all three
amplitudes nonzero. Three distinct projections give independent columns
`(1,b1,b2)`; a repeated projection forces the third amplitude to zero.

Suppose first that only one half is active and its weight is at most five.
The global b3 check then supplies the fourth affine check in that half.
Four dependent affine columns are a binary affine plane, and their unique
dependence has four equal amplitudes. Every five-element subset of GF(2)^3
contains such a plane. It spans affine rank four, so its one-dimensional
kernel is supported on that plane and has no word of exact weight five.

For independent plane directions u,v in GF(2)^3, summing q on their affine
plane gives

\[
(u_1v_2+u_2v_1)
+2(u_1v_3+u_3v_1)
+4(u_2v_3+u_3v_2).
\]

At least one of these three binary coefficients is nonzero. Independence of
1,2,4 over GF(2) makes the displayed sum nonzero. The b0*b3 part of f is
affine within a half and has zero plane sum. Thus the f check excludes
every remaining weight-four word in one half.

If both halves are active and total weight is at most five, their weights
must be two and two; the three-coordinate case was excluded above. Let the
pairs have fixed projections `(a1,a2)` and `(c1,c2)`. The global b3 check
forces their nonzero amplitudes to be equal, say alpha. Their combined
f check is

\[
\alpha\bigl(1+2(a_1+c_1)+4(a_2+c_2)\bigr),
\]

which is nonzero. Consequently every nonzero quiet vector has weight at
least six. This argument applies over F, not only to binary amplitudes.

## Circuit budget and untested obligations

A direct circuit evaluates H* on x, then adds the three rank-one corrections.
Each g1 or g2 correction costs three XORs to form its four-input sum and
four XORs to apply it. For f, write `z_j=x_(0,j)+x_(1,j)`, where
`j=b1+2*b2+4*b3`. Only z3,z5,z6,z7 are needed, and

\[
q^Tx=(z_3+z_7)+2(z_5+z_7)+4(z_6+z_7).
\]

Forming f^Tx costs thirteen XORs and two fixed products, including the
four-input b0*b3 sum. Applying f times this scalar costs thirteen XORs and
two fixed products: f has nine nonzero entries, and its coefficients lie
in `{1,2,3,4,5,6,7}`. Therefore this explicit, unshared correction circuit
adds at most forty XORs and four fixed field products to H*.

The resulting budget is at most 106 XORs and twelve fixed field products
per sixteen-symbol column, versus H*'s 66 and eight. Diagonal keys and
variable-product counts do not change. Shared intermediate sums may reduce
the XOR count, but no such native circuit has been implemented or timed.

Invertibility and quiet minimum weight do not establish nonquiet spreading.
The nonbinary correction can cancel entries of H*. The next local checks,
if this proposal is pursued, are all matrix entries, labelled pair laws,
and actual cross-fiber ranks. The inverse must be computed anew; arguments
that identify it with H' would be invalid. No exact enumeration, regional
provider, benchmark, or certificate was run for this proposal.

## Why a schedule alone has a limit

With a active RS[16,8] rows, the message space has dimension 8a. Under any
deterministic input permutation schedule and any fixed nonzero diagonals,
all sixteen rank-five H* quiet conditions impose at most eighty linear
constraints. Thus the all-quiet RS subspace has dimension at least 8a-80.
It is nonzero for a>=11, although a witness may use fewer than all a rows.
A sixth syndrome merely moves this dimension obstruction to a>=13.

This does not refute the current code's distance target: quiet is not the
same as low output weight. It does show why a schedule cannot universally
remove the quiet branch. The present proposal changes the quiet code
itself, at an explicit arithmetic cost. Continue the existing inverse-schedule
coverage first; revisit this candidate only if quiet low-weight words,
rather than looseness in the nonquiet bounds, remain decisive.
