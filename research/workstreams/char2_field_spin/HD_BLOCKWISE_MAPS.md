# Fixed maps and a local bound for blockwise field-SPIN

## Result and exact maps

The first native candidate uses a 64-symbol packet and a 16-symbol state.
The fixed maps below have fast additive-FFT circuits and satisfy \(CA=0\).
That relation preserves the feedback syndrome through emission. It also
supports a positive local bound that retains arbitrary state dependence.
This note proves the local bound, not a whole-code distance certificate.

Let \(\mathbb F=\mathrm{GF}(2^{16})\), with \(q=65536\). Embed the AES
\(\mathrm{GF}(256)\) coefficient field in \(\mathbb F\). Byte labels refer
to its polynomial basis modulo \(X^8+X^4+X^3+X+1\). Define

\[
 U=\{0,\ldots,63\},\quad V=\{128,\ldots,143\},\quad
 Q_V(X)=\prod_{v\in V}(X+v),\quad D_V=Q'_V.
\]

Both node sets are additive affine spaces, so \(D_V\) is a nonzero constant.
For \(u\in U\) and \(v_j=128+j\), set

\[
 A_{u,j}=\frac{Q_V(u)}{D_V(u+v_j)},\qquad C=A^T.
\]

Thus \(A\) interpolates a polynomial of degree below 16 from its values on
\(V\), then evaluates it on \(U\). The 16 state coordinates remain these
original values. No hidden change to coefficient-state coordinates is made.
In the additive-FFT notation, \(A\) is a stack of four maps

\[
 \operatorname{FFT}_{16}(b+128)\operatorname{IFFT}_{16}(0),
 \qquad b\in\{0,16,32,48\}.
\]

Injection uses the literal transpose circuit. This is the same map-building
pattern as the existing `Cauchy256InnerFft.h`, with larger rectangular maps.
The row-major bytes of \(A\) have SHA256
`4a8399e103c2fb2ca5824883d7da1f6d3e235c4e5873c8467701a6cde071586a`.

The maps are row-scaled Cauchy matrices, so every square minor of admissible
size is nonzero. Consequently, \(\operatorname{im}A\) is an MDS
\([64,16,49]\) code and \(\ker C\) is an MDS \([64,48,17]\) code.
These claims remain valid over the extension field \(\mathbb F\).

For \(CA=0\), write the columns of \(A\) as evaluations of Lagrange
polynomials \(\ell_j\), each of degree at most 15. Every product
\(\ell_j\ell_k\) has degree at most 30. Its sum over \(U\) is zero:
Lagrange interpolation on \(U\) identifies this sum, divided by the constant
derivative of its node polynomial, with the coefficient of \(X^{63}\).
That coefficient vanishes. Therefore

\[
 CA=0,\qquad (AC)^2=0.
\]

The blockwise recurrence remains

\[
 Y=X+AZ,\qquad Z'=MZ+CX.
\]

The initial state is zero; the final state is free. In particular,
\(CY=CX\). The emitted weight and the return event cannot be analyzed as
independent events. The following bound uses this equality explicitly.

## Deterministic emission constraints

For this section allow general packet width \(t\), state dimension \(s\),
and the same MDS and \(CA=0\) properties. Fix an input packet \(x\) of weight
\(w\), and let \(c=Cx\). Put \(d_A=t-s+1\).

For any live state \(z\ne0\), the triangle inequality gives
\(\operatorname{wt}(x+Az)\ge\max(0,d_A-w)\). If \(c\ne0\), then
\(C(x+Az)=c\), so the emitted packet cannot be zero. Define

\[
 m_1(w)=\max(1,d_A-w).
\]

If \(c=0\), the emitted packet belongs to \(\ker C\). For \(w<d_A\), it
cannot vanish: \(x=Az\) for a live state would contradict the distance
of \(\operatorname{im}A\). Hence define

\[
 m_0(w)=
 \begin{cases}
 \max(s+1,d_A-w),&w<d_A,\\
 0,&w\ge d_A.
 \end{cases}
\]

These minima cover all fixed packet amplitudes. Also, \(c=0\) is impossible
when \(1\le w\le s\), by the distance of \(\ker C\). At \(w=0\), only
\(c=0\) is possible. Both cases remain possible in the relaxed bound for
\(w\ge s+1\).

## Emission from a uniform live state

Fix \(0<\theta\le1\). We need an upper bound for

\[
 E_x=\frac1{q^s-1}\sum_{z\ne0}\theta^{\operatorname{wt}(x+Az)}
\]

without assigning random amplitudes to \(x\). Expand each output coordinate
using \(\theta^{\mathbf1[y\ne0]}=\theta+(1-\theta)\mathbf1[y=0]\).
For a selected set of zero coordinates, let \(a\) lie inside
\(\operatorname{supp}x\) and \(b\) outside it.

If \(b\ge s\), the homogeneous equations outside the input support force
\(z=0\); there are no live solutions. For \(b<s\), the following count is
valid uniformly over the selected coordinates and the fixed amplitudes:

\[
 n(a,b)=
 \begin{cases}
 q^{s-b}-1,&a=0,\\
 q^{\max(0,s-a-b)},&a>0.
 \end{cases}
\]

For \(a=0\), this is the exact kernel size with zero removed. For \(a>0\),
the system cannot contain the zero state. Its rank is \(\min(s,a+b)\).
The bound drops only the consistency condition when there are at least
\(s\) equations.

It follows that \(E_x\le E_w\), where

\[
 E_w:=\min\left\{
 \theta^{\max(0,d_A-w)},
 \frac1{q^s-1}
 \sum_{a=0}^w\sum_{b=0}^{\min(t-w,s-1)}
 {w\choose a}{t-w\choose b}n(a,b)
 \theta^{t-a-b}(1-\theta)^{a+b}
 \right\}.
\]

All summands are nonnegative. For syndrome flag \(j\in\{0,1\}\), write
\(L_j=\theta^{m_j(w)}\) and \(E_j=\min(E_w,L_j)\). Then \(E_x\le E_j\)
whenever \(\mathbf1[c\ne0]=j\). For \(w=0\), the positive expansion is
exact for the nonzero MDS emission code.

## A composable bound on weighted state measures

Use the rank-one setup distribution from
[HD_BLOCKWISE_INNER.md](HD_BLOCKWISE_INNER.md), with \(r\) independent
refresh factors before the single feedback addition. Its nonzero-state
kernel is

\[
 P^r=\alpha I+\beta L+\gamma J,\qquad
 \alpha=a_*^r,\quad
 \beta=q^{-r}-a_*^r,\quad
 \gamma=1-q^{-r},\quad
 a_*=\frac{q^{s-1}-1}{q^s-1}.
\]

Here \(L\) samples within the same nonzero scalar line and \(J\) samples
all nonzero states. Put \(\varepsilon=\alpha+\beta=q^{-r}\) and
\(D=q^s-1\). The maximum return probability from any fixed live state to a
specified nonzero syndrome is

\[
 \rho=\alpha+\frac\beta{q-1}+\frac\gamma D.
\]

For one round, \(\rho=q^{s-1}/D\). This worst-case value is not the return
bound for a uniform incoming state.

Let \(\mu\) be the nonnegative state measure after weighting previous
packets by \(\theta\) to their emitted weights. Represent an upper bound by
three nonnegative numbers \((Z,B,U)\):

\[
 \mu(0)\le Z,\qquad
 \mu(z)\le\nu(z)+U/D\quad(z\ne0),\qquad
 \sum_{z\ne0}\nu(z)\le B
\]

for some nonnegative measure \(\nu\). Thus \(B\) allows arbitrary live
state dependence, while \(U\) is a pointwise uniform-live density bound.
These are overlapping upper-bound components, not probabilities of a
partition of the setup space.

For each fixed input \(x\), the output-weighted update is a positive linear
operator on \(\mu\). Fresh \(M\) is independent of this incoming measure.
The following matrices update the row vector \((Z,B,U)\); rows and columns
are ordered as zero, arbitrary-live, uniform-live.

If \(c=0\), a valid matrix is

\[
 T_0(w)=
 \begin{pmatrix}
 \theta^w&0&0\\
 0&\varepsilon L_0&\gamma L_0\\
 0&0&\varepsilon L_0+\gamma E_0
 \end{pmatrix}.
\]

If \(c\ne0\), put \(H_1=\varepsilon L_1+\gamma E_1\). A valid matrix is

\[
 T_1(w)=
 \begin{pmatrix}
 0&\theta^w&0\\
 \rho L_1&\varepsilon L_1&\gamma L_1\\
 H_1/D&0&H_1
 \end{pmatrix}.
\]

For the arbitrary-live row, the point and line branches retain arbitrary
dependence and have total weighted mass at most \(\varepsilon L_jB\).
The full-refresh branch, followed by translation by \(c\), has density at
most \(\gamma L_jB/D\) at every live successor. Returns have probability
at most \(\rho\), and \(CY=c\ne0\) supplies the factor \(L_1\).
If \(c=0\), invertibility of \(M\) forbids returns.

For the uniform-live row, the input density after emission is at most
\(L_jU/D\). Both \(I\) and \(L\) preserve uniform-live density; translation
preserves a pointwise upper bound. Their live-successor density is therefore
at most \(\varepsilon L_jU/D\). The full-refresh contribution instead uses
the average emitted weight: its density is at most \(\gamma E_jU/D\).
The same bound at the zero successor gives \(H_1U/D\). Zero is impossible
when \(c=0\).

This proof retains the joint emission/return event. It does not multiply an
independent output probability by an independent return probability.

### A second choice for the uniform-live row

There is an alternative that can be tighter when \(L_j\) greatly exceeds
\(E_j\). Place the weighted point/line branches into the arbitrary-live
component by total mass, rather than pointwise density. Replace only the
last row by

\[
 (0,\varepsilon E_0,\gamma E_0)\quad\text{if }c=0,
 \qquad
 (H_1/D,\varepsilon E_1,\gamma E_1)\quad\text{if }c\ne0.
\]

The nonzero-state mass before these branches is at most \(UE_j\), so their
weighted mass is at most \(\varepsilon UE_j\). The return entry remains a
separate upper bound. Each version is composable; a prescribed packet may
use either version. Taking entrywise minima between the two is not justified,
because their different component allocations certify different measures.

If only \(w\) is retained, take the entrywise maximum over feasible syndrome
flags within one chosen version. This loses information but remains valid for
every fixed input of that weight. No input diagonal, independent input
amplitude, or fresh dense-GL assumption is used. Consequently, one may first
condition on a routed outer word, apply these operators, and only then average
over the route. The source/route averaging still requires its own proof.

## Checkpoint and next bound

The helper [theory/hd_blockwise_maps.py](theory/hd_blockwise_maps.py) uses exact
integer and rational arithmetic. Five
[tests](theory/hd_blockwise_maps_tests.py) pass, covering the native matrix
digest and \(CA=0\), every small-map minor, all 4096 small-field inputs for
the emission bound, 2994 local measure checks, and ten consecutive packets.
The last test alternates the two uniform-live allocations and verifies the
pointwise residual-mass invariant, not only the final scalar moment.

As a numerical illustration, for \(t=64,s=16,q=65536,\theta=1/2\), the
uniform emission bound is approximately \(5.4263\cdot10^{-20}\) at the
tested representative weights. This is about \(2^{-64}\). The return bound
from the uniform component also pays the \(q^{16}-1\) denominator. These
are local generating-function bounds, not distance or failure guarantees.

The decisive next step is to insert the three-component operator into a
regional source bound for the fixed outer. Retain both zero-syndrome and
live-syndrome information where affordable. The exact matrices, packet
chronology, and mixer distribution are now fixed for that calculation.
Native timing remains a separate validation step owned by the implementation
task; this checkpoint changes no production code and asserts no full margin.
