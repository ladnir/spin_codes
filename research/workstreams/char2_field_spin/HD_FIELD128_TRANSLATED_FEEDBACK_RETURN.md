# Return counts for translated feedback at state 32

Changing the feedback interpolation nodes can strengthen the return bound
without changing the emission code. This note gives a uniform rank bound
for one such change. It does not establish a distance certificate or a
performance claim. The existing gamma30 certificates concern a different
construction and remain unchanged.

## Maps and the missing full-rank premise

Let \(E=\operatorname{GF}(256)\subset F=\operatorname{GF}(2^{128})\), and
write \(q=|F|\). Byte labels below use the AES polynomial basis of \(E\)
(modulus 0x11b); addition is XOR. Set \(U=\{0,\ldots,31\}\),

\[
L_U(X)=\prod_{u\in U}(X+u),\qquad \ell=L_U'(X)\in E^*.
\]

For each \(\delta\in E\), let \(T_\delta\) translate the evaluations on
\(U\) of a polynomial of degree below 32. Thus
\(T_\delta T_\epsilon=T_{\delta+\epsilon}\). These matrices are symmetric.
For \(\delta\notin U\), their entries are

\[
(T_\delta)_{x,y}=\frac{L_U(\delta)}{\ell(\delta+x+y)},\qquad x,y\in U.
\]

Let \(A\) interpolate from nodes 128..159 to nodes 0..63. Let \(B\)
interpolate from nodes 64..95 to those same output nodes, and set \(C=B^T\).
In output halves \(H_0=U\) and \(H_1=32+U\),

\[
A=\begin{pmatrix}T_{128}\\T_{160}\end{pmatrix},\quad
B=\begin{pmatrix}T_{64}\\T_{96}\end{pmatrix},\quad
P:=AC=\begin{pmatrix}T_{192}&T_{224}\\T_{224}&T_{192}\end{pmatrix}.
\]

Consequently \(CA=0\), \(P^2=0\), and \(\operatorname{rank}P=32\).
Both \(A\) and \(B\) generate the same MDS evaluation code. Every set of
at most 32 rows of either matrix is independent. Also \(B=A T_{192}\).

The entries of \(P\) are

\[
P_{x,y}=\frac{L_U(192+x+y)}{\ell(192+x+y)},\qquad x,y\in H_0\cup H_1.
\]

All entries are nonzero, but this is not one globally scaled Cauchy matrix.
The numerator differs between equal-half and unequal-half blocks. Mixed
minors can be singular: the exact GF256 probe reports rows \(0,32\) and
columns \(6,44\), with matrix
\(\left(\begin{smallmatrix}145&85\\136&32\end{smallmatrix}\right)\)
and determinant zero. The probe source is
[hd_field128_feedback_rank_probe_v41.py](theory/hd_field128_feedback_rank_probe_v41.py).
Density alone therefore does not justify full rank for every support submatrix.

The encoder uses the recurrence

\[
Y=X+AQ,\qquad Q_{\rm next}=D_2H_{32}D_1(Q+\gamma CX).
\]

Use the \(H_{32}\) refresh from the [random-gamma construction](HD_FIELD128_RANDOM_GAMMA32.md),
with fresh independent nonzero diagonal entries in \(D_1,D_2\).
The outer ensemble, paired route, initial zero state, and unrestricted final
state are unchanged. Sample one global scalar
\(\gamma\) independently of those random choices from

\[
\Gamma_{16}:=F\setminus\operatorname{GF}(2^{64}).
\]

These are exactly the elements of degree 16 over \(E\). The local bounds
below are pointwise in every such scalar satisfying the stated setup
conditions; no scalar is resampled at a packet or region boundary.

## Principal invertibility

Define \(\mathcal G\) by requiring
\(M_D(\gamma):=I_{32}+\gamma C_DA_D\) invertible for every coordinate set
\(D\subseteq H_0\cup H_1\). Here \(A_D\) selects rows and \(C_D\) selects
columns. The diagonal of the symmetric matrix \(P\) is the nonzero constant

\[
a=\frac{L_U(192)}{192\ell}.
\]

For \(c=|D|\), Sylvester's identity and cancellation of reversed cycles in
characteristic two give

\[
p_D(x):=\det(I_c+xP_{D,D})=\det M_D(x)
       =(1+ax)^{c\bmod2}R_D(x)^2,\qquad R_D(0)=1.
\]

The matching expansion is the same as in the
[fixed-gamma proof](HD_FIELD128_GAMMA_RETURN.md), with \(1+x\) replaced by
\(1+ax\). Since \((I+xP)^{-1}=I+xP\) and \(\det(I+xP)=1\), Jacobi's
identity gives \(p_D=p_{D^c}\). Hence

\[
\deg p_D\le\min(c,64-c,32).
\]

If \(c\ne32\), then \(\deg R_D\le15\), so no element of \(\Gamma_{16}\)
is a root. The possible root \(a^{-1}\in E\) is also excluded. Only sets
of size 32 remain. Pairing them with their complements yields

\[
\Pr_{\gamma\gets\Gamma_{16}}[\neg\mathcal G]
\le \varepsilon_{\rm princ}:=
\frac{8\binom{64}{32}}{q-2^{64}}<2^{-64}.
\]

The looser setup bound in the [random-gamma note](HD_FIELD128_RANDOM_GAMMA32.md)
would also suffice. Neither argument certifies an arbitrary fixed tower
generator at state 32 without checking its setup conditions.

## A return constraint and its baseline rank

Fix a potential input support \(S\) of size \(j\). Let \(X\gets F^S\)
mean that its coordinates on \(S\) are independent uniform field elements
and its other coordinates are zero. Fix active output-zero coordinates
\(D\subseteq S\) and inactive output-zero coordinates \(Z\subseteq S^c\).
Write

\[
c=|D|,\quad e=|Z|,\quad T=S\setminus D,\quad t=|T|=j-c,
\quad u=\min(32,t).
\]

For the return contribution put \(R=CX\ne0\); the required incoming state
is \(Q=\gamma R\). Under \(\mathcal G\), the output-zero equations become

\[
X_D=\gamma A_DR,\qquad
R=M_D(\gamma)^{-1}C_TX_T,\qquad
A_ZR=0.
\]

Thus the remaining constraint on \(X_T\) is the matrix

\[
\mathcal B_{Z,D,T}(\gamma):=A_ZM_D(\gamma)^{-1}C_T.
\]

At \(\gamma=0\), this matrix equals \(P_{Z,T}\). Put
\(e_i=|Z\cap H_i|\), \(t_i=|T\cap H_i|\). If all selected rows lie in
one half, the resulting matrix is a column-scaled Cauchy matrix. If all
selected columns lie in one half, it is a row-scaled Cauchy matrix.
The scaling factors are nonzero. The MDS factorization \(P=AC\) also gives
the last term in the following lower bound:

\[
h_{\rm half}:=\max\{0,\min(e_0,t),\min(e_1,t),
                    \min(e,t_0),\min(e,t_1),
                    \min(e,32)+\min(t,32)-32\}
\le\operatorname{rank}P_{Z,T}.
\]

For an occupancy-only formula, the pigeonhole principle gives the weaker
bound

\[
h:=\max\{0,\min(\lceil e/2\rceil,t),
             \min(e,\lceil t/2\rceil),
             \min(e,32)+u-32\}.
\tag{1}
\]

Either bound holds for every placement of the indicated coordinates. A
computed exact rank of \(P_{Z,T}\) can replace them, but cannot be inferred
from a few nonsingular examples.

## Degree-16 feedback loses at most one baseline dimension

Fix any positive integer \(h\) for which \(P_{Z,T}\) has a nonsingular
\(h\)-square submatrix. Fix such row and column subsets \(Z',T'\) before
choosing \(\gamma\). Consider the polynomial

\[
F_{D,Z',T'}(x):=
\det\begin{pmatrix}I_{32}+xC_DA_D&C_{T'}\\A_{Z'}&0\end{pmatrix}.
\]

Its value at zero is \(\det P_{Z',T'}\ne0\). Its degree is at most

\[
d_h:=\min(c,64-c,32-h).
\tag{2}
\]

The variable part has rank at most \(c\); determinant multilinearity
prevents a higher degree. Also \(C_DA_D=C_{D^c}A_{D^c}\) because \(CA=0\),
giving the bound \(64-c\). Finally, the last \(h\)
columns must use \(h\) of the first 32 rows. At most \(32-h\) remaining
rows can contribute variable entries.

Where \(M_D(x)\) is invertible, the Schur complement identity gives

\[
F_{D,Z',T'}(x)=\det M_D(x)\,
\det\bigl(A_{Z'}M_D(x)^{-1}C_{T'}\bigr).
\]

Suppose the full matrix \(\mathcal B_{Z,D,T}(\gamma)\) has rank \(b<h\).
Constant invertible row operations make at least \(h-b\) rows of the
selected Schur complement vanish at \(\gamma\). Its determinant has a
zero of multiplicity at least \(h-b\). Because \(\det M_D(\gamma)\ne0\),
the same holds for \(F_{D,Z',T'}\). Its coefficients lie in \(E\), so the
degree-16 minimal polynomial of \(\gamma\) divides it to that power.
Therefore

\[
\operatorname{rank}\mathcal B_{Z,D,T}(\gamma)
\ge h-\left\lfloor\frac{\min(c,64-c,32-h)}{16}\right\rfloor.
\tag{3}
\]

For \(h>0\), the loss is at most one. It is zero if \(c\le15\) or
\(64-c\le15\), or \(h\ge17\). This bound requires no setup event beyond \(\mathcal G\).
For \(h=0\), use the trivial zero lower bound instead of (3).

There is a stronger degree bound for the return problem. Set
\(V=(D\cup T\cup Z)^c\) and \(v=64-j-e\). The original constraints say
that \(X\) is supported on \(D\cup T\), while
\(Y=(I+\gamma P)X\) is supported on \(T\cup V\). Since

\[
(I+\gamma P)^2=I,\qquad C(I+\gamma P)=C,
\]

exchanging \(X\) and \(Y\) bijects this solution space with the reversed
problem. That problem has active zero set \(V\), inactive zero set \(Z\),
and the same overlap \(T\). The bijection also preserves the syndrome.
Under \(\mathcal G\), elimination in the two directions gives dimensions
\(t-\operatorname{rank}\mathcal B_{Z,D,T}\) and
\(t-\operatorname{rank}\mathcal B_{Z,V,T}\), respectively. Thus these ranks
are equal. Their baseline matrix is the same \(P_{Z,T}\).
Applying (3) in both directions gives

\[
\operatorname{rank}\mathcal B_{Z,D,T}(\gamma)
\ge h-\left\lfloor\frac{\min(c,v,32-h)}{16}\right\rfloor.
\tag{3a}
\]

This is a bijection of finite solution spaces, not a replacement of the
input law. Probabilities below still use the original denominator \(q^j\).

## An implementable joint-return cap

If \(e\ge32\), then \(A_ZR=0\) forces \(R=0\), and the return event is
empty. If \(t=0\), it is also empty. Otherwise assume \(e<32\) and define

\[
b_{\rm old}:=\max(0,e+u-32),\qquad
b:=\max\left\{b_{\rm old},\,
h-\left\lfloor\frac{\min(c,v,32-h)}{16}\right\rfloor\right\},
\tag{4}
\]

where the second argument is defined as zero when \(h=0\). Here \(h\)
may be (1), \(h_{\rm half}\), or an established exact baseline rank.
The old lower bound follows from the dimensions of \(M_D^{-1}\operatorname{im}C_T\)
and \(\ker A_Z\), independent of the Cauchy structure.

For the actual rank \(b_*\), the number of allowed \(X_T\) is
\(q^{t-b_*}\). Of these, exactly \(q^{t-u}\) have \(C_TX_T=0\).
Those vectors have \(R=0\) and \(X_D=0\), so they must be removed.
Dividing by \(q^j\), and using \(b_*\ge b\), gives

\[
\Pr_X[R\ne0,\ Y_{D\cup Z}=0]
\le p_{j,c,e}:=q^{-c-b}\bigl(1-q^{-(u-b)}\bigr).
\tag{5}
\]

The right side is exactly zero when \(b=u\). This is a zero-syndrome
subtraction, not an independence assumption about active and inactive zeros.
It improves or equals the old dimension-only count term by term.

For fixed \(S,\gamma\), define the true return moment
\(J_{S,\gamma}(z):=\mathbb E_X[z^{\operatorname{wt}(Y)}1_{CX\ne0}]\).
Using the occupancy-only (1), its positive expansion for \(0<z<1\) gives
the uniform upper cap

\[
J_{S,\gamma}(z)\le\overline J_j(z):=
\sum_{c=0}^{j-1}\sum_{e=0}^{\min(31,64-j)}
\binom jc\binom{64-j}e\,
z^{64-c-e}(1-z)^{c+e}p_{j,c,e}.
\tag{6}
\]

Set \(\overline J_0=0\), and take the minimum with the existing valid basic return
caps. For a positive nonzero term, its probability logarithm is
\(-(c+b)\log_2q+\log_2(1-q^{-(u-b)})\); skip the exact zero case.

Only the local entry \(\overline J_j/(q^{32}-1)\) changes. The quiet, birth, and
live moments and the refresh density retain their formulas: both feedback
kernels are MDS, while emission is unchanged. The incoming-state atom bound
still applies because this estimate holds for every input support and every
fixed good scalar. Multiply chronological packet matrices before applying
the fractional power. Keep route multiplicities and support counts outside
that power. A half-specific version must retain the existing 16 pair slots
in each half; it cannot replace routing by independent half draws.

## What reversal removes, and when another setup event is needed

For \(j=32,c=e=31\), one has \(t=u=h=v=1\). The forward-only degree
bound permits losing the sole rank dimension. Reversal instead gives
degree at most one, so \(b=u=1\). The former probability
\(q^{-31}(1-q^{-1})\) is replaced by exact zero.

More generally, suppose \(c+e\ge60\) and \(t\ge1\). Then
\(t+v\le4\), so \(v\le3\) and there is no degree-16 rank loss.
Whenever the baseline bound satisfies \(h=t\), the entire return term
is zero. In particular, all such shapes that would incur the forward-only
loss have \(c\le48\), hence \(e\ge12\) and \(t\le4\). Formula (1)
then gives \(h=t\). These high-zero terms need no targeted gamma exclusion.

For other shapes, an additional setup condition can preserve a chosen
baseline rank when \(d:=\min(c,v,32-h)\ge16\). Fix one nonsingular
baseline minor for every coordinate placement, choosing the smaller
elimination direction before sampling \(\gamma\). Each resulting nonzero
polynomial has degree at most 31, hence at most 16 roots in \(\Gamma_{16}\).
For one shape \((j,c,e)\), the additional bad probability is at most

\[
\frac{16}{q-2^{64}}\frac{64!}{c!(j-c)!e!(64-j-e)!}.
\tag{7}
\]

This preserves rank \(h\), not necessarily \(\min(e,t)\). Sum (7) only
over an explicitly selected list of shapes. A union over all four-way
partitions is too large for this argument to give a useful security bound.
No additional setup event is used by formulas (1)--(6).

If a later distance proof gives conditional failure at most \(\varepsilon\)
for every scalar satisfying \(\mathcal G\), its unconditional bound is
\(\varepsilon+\varepsilon_{\rm princ}\). Any separately chosen shape
costs are added once, not per source support, packet, or region. None of
these setup conclusions transfers a previous numerical receipt.

## Smallest discriminating screen

At the saved state32, 45% literals and original profile \((2026,2026)\),
first replace only the return entry by (6) with (1)--(5). Retain all saved
markers, original normalization, projective divisor, and outer prefactor.
Compare against the untouched full-return control. This is a conditional
FLOAT profile diagnostic, not a directed or all-profile certificate.

Compare the forward-only loss from (3) with the reversal loss from (3a).
The former leaves the \((32,31,31)\) term intact; the latter deletes it
without another setup cost. This separates the benefit of dense blocks
from the benefit of preserving the input/output symmetry. If the gain is
material, inspect the remaining return summands before expanding a
protected shape list or implementing full half-aware coefficients.

There is an unavoidable floor on the true return moment. For \(j>0\),
\(\operatorname{rank}C_S=\min(32,j)\) and \(\operatorname{wt}(Y)\le64\), so

\[
J_{S,\gamma}(z)\ge z^{64}\bigl(1-q^{-\min(32,j)}\bigr).
\]

Divide by \(q^{32}-1\) for the corresponding return-entry floor. This is
a lower bound, not an upper cap. Replacing the upper cap by this floor is
an optimistic diagnostic only. Once an upper cap is close to the floor at
fixed \(z\), sharper ranks offer little further gain for that occupancy;
delaying the fractional power is a separate proof refinement.
