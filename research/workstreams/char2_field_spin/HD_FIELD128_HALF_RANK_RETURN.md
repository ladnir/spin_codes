# A return-cap refinement from the two coordinate halves

The current gamma return cap counts active and inactive output-zero constraints jointly, but uses only their separate MDS ranks. The fixed Gram matrix provides an additional cross-rank bound. This note derives that bound and a uniform, occupancy-only return cap. It changes neither the construction nor the route.

Work with the construction in [the gamma return note](HD_FIELD128_GAMMA_RETURN.md). Let \(K=\operatorname{GF}(256)\), \(F=\operatorname{GF}(2^{128})\), and \(q=|F|\). Fix \(s\in\{24,28,30\}\). The matrix \(A\in K^{64\times s}\) evaluates the degree-less-than-\(s\) Lagrange basis on nodes \(0,\ldots,63\), with interpolation nodes \(128,\ldots,127+s\). Put \(C=A^{\mathsf T}\). Fix \(\gamma\in F\) of degree 16 over \(K\).

The existing principal-nonsingularity lemma applies: \(CA=0\), every row of \(A\) sums to one, and

\[
L_D(\gamma):=I_s+\gamma C_DA_D
\]

is invertible for every coordinate set \(D\subseteq\{0,\ldots,63\}\). Here \(A_D\) contains the rows in \(D\), and \(C_D=A_D^{\mathsf T}\). No randomness is added to \(A\), \(C\), or \(\gamma\).

## The dimension that controls a returning state

Fix a potential input support \(S\) of size \(j\). The comparison input \(X_S\gets F^S\) is uniform and may have zero coordinates; \(X_{S^c}=0\). Define

\[
R=C_SX_S,\qquad Y=X+\gamma AR.
\]

The return moment to be bounded is

\[
J_S(z):=\mathbb E_{X_S}\left[\mathbf1_{R\ne0}\,z^{\operatorname{wt}(Y)}\right],
\qquad 0<z<1.
\]

This is the same returning-state moment used by the current local cap. The packet transfer entry is the return cap divided by \(q^s-1\).

Fix output-zero sets \(D\subseteq S\) and \(Z\subseteq S^c\). Write

\[
c=|D|,\quad e=|Z|,\quad T=S\setminus D,\quad t=|T|=j-c,
\qquad u=\min(s,t).
\]

The active-zero equations determine \(X_D=\gamma A_DR\). Substituting in \(R=C_DX_D+C_TX_T\), and using the inactive-zero equations, gives

\[
L_D(\gamma)R=C_TX_T,\qquad A_ZR=0.
\tag{1}
\]

Define the matrix and its rank

\[
B_{Z,D,T}(\gamma)=A_ZL_D(\gamma)^{-1}C_T,
\qquad r_\gamma=\operatorname{rank}_F B_{Z,D,T}(\gamma).
\]

The map \(C_T\) has rank \(u\). Since \(L_D(\gamma)\) is invertible, the allowed syndrome space in (1) has dimension exactly \(u-r_\gamma\). Every allowed syndrome has \(q^{t-u}\) preimages \(X_T\), and then \(X_D\) is unique. Consequently,

\[
\Pr_{X_S}[R\ne0,\ Y_D=0,\ Y_Z=0]
=\frac{(q^{u-r_\gamma}-1)q^{t-u}}{q^j}.
\tag{2}
\]

For \(e\ge s\), \(A_Z\) has full column rank, so this probability is zero. For \(e<s\), the existing cap uses only

\[
r_\gamma\ge\max(0,e+u-s),
\qquad u-r_\gamma\le\min(u,s-e).
\tag{3}
\]

The refinement below supplies another lower bound on \(r_\gamma\).

## Rank stability at a degree-16 field element

Let \(x\) be an indeterminate over \(K\), and put

\[
B(x)=A_Z(I_s+xC_DA_D)^{-1}C_T.
\]

Suppose a \(k\times k\) minor of \(B(0)=A_ZC_T\) is nonzero, with \(k\ge1\). Let \(Z'\subseteq Z\) and \(T'\subseteq T\) index this minor. Its bordered numerator is

\[
N(x)=\det\begin{pmatrix}
I_s+xC_DA_D&C_{T'}\\
A_{Z'}&0
\end{pmatrix}
=\det L_D(x)\,\det B(x)_{Z',T'}.
\tag{4}
\]

Signs disappear in characteristic two. In particular, \(N(0)\ne0\).

There are three useful degree bounds. In the determinant on the left of (4), the \(k\) bottom rows must occupy \(k\) state columns. At most \(s-k\) remaining state columns can contribute an \(x\), so \(\deg N\le s-k\). The Woodbury form gives the same numerator as

\[
N(x)=\det\begin{pmatrix}
I_c+xA_DC_D&A_DC_{T'}\\
xA_{Z'}C_D&A_{Z'}C_{T'}
\end{pmatrix}.
\tag{5}
\]

Only the first \(c\) columns depend on \(x\), so \(\deg N\le c\). Finally, \(CA=0\) implies \(C_DA_D=C_{D^c}A_{D^c}\). Applying (5) with \(D^c\) gives

\[
\deg N\le\nu(k,c):=\min(c,64-c,s-k).
\tag{6}
\]

Let \(m_\gamma(x)\in K[x]\) be the degree-16 minimal polynomial of \(\gamma\). The denominator \(\det L_D(\gamma)\) is nonzero by the existing principal lemma. If \(r_\gamma\le k-\delta\), the selected \(k\times k\) minor matrix has rank deficiency at least \(\delta\) at \(\gamma\). Its determinant, and hence \(N\), vanishes there to multiplicity at least \(\delta\). One way to see this is to use invertible constant row and column operations at \(\gamma\): each determinant term then uses at least \(\delta\) entries vanishing at \(\gamma\).

The same multiplicity holds at each conjugate of \(\gamma\). Finite fields are perfect, so these 16 conjugates are distinct. Thus \(m_\gamma(x)^\delta\mid N(x)\). Since \(N\ne0\), (6) yields

\[
\boxed{\quad
r_\gamma\ge k-\left\lfloor\frac{\min(c,64-c,s-k)}{16}\right\rfloor.
\quad}
\tag{7}
\]

This statement applies to any certified nonzero minor size \(k\), not only the exact rank of \(A_ZC_T\). If \(c\le15\), the minor loses no rank at \(\gamma\). It also loses no rank when \(s-k\le15\) or \(64-c\le15\). For \(s\le30\) and \(k\ge1\), the loss is at most one.

The square-factor lemma is used here only to establish denominator nonvanishing. It does not halve the degree of the arbitrary bordered numerator \(N\).

## A fixed Gram cross-rank from the coordinate halves

Let \(H=\{0,\ldots,31\}\subset K\), an additive subgroup. The evaluation coordinates split into \(H_0=H\) and \(H_1=32+H\). The full interpolation coset is \(128+H=\{128,\ldots,159\}\). Define

\[
\begin{aligned}
d&=32-s,\\
p(X)&=\prod_{b=128}^{127+s}(X+b),\\
m(X)&=\prod_{b=128+s}^{159}(X+b),\\
P(X)&=p(X)m(X)=\prod_{b=128}^{159}(X+b).
\end{aligned}
\]

The affine linearized polynomial \(P\) has constant derivative \(a=P'\ne0\). It is constant on each evaluation half:

\[
P(x)=\kappa_i\quad(x\in H_i),
\qquad \kappa_0,\kappa_1\ne0,\quad\kappa_0\ne\kappa_1.
\]

For distinct evaluation coordinates \(x,y\), the Gram kernel \(G(x,y)=(AA^{\mathsf T})_{x,y}\) satisfies

\[
G(x,y)=\frac{p(x)m(y)+m(x)p(y)}{a(x+y)}.
\tag{8}
\]

To verify (8), the Lagrange coefficient at node \(b\) is \(p(x)/((x+b)p'(b))\). At such a node, \(p'(b)m(b)=a\). Since \(d<s\), partial fractions give

\[
\sum_b\frac{1}{p'(b)^2(X+b)}=\frac{m(X)}{a p(X)}.
\]

Take the sum of this identity at \(x\) and \(y\), divide by \(x+y\), and multiply by \(p(x)p(y)\) to obtain (8).

The sets \(Z,T\) are disjoint, so all denominators \(x+y\) below are nonzero. Diagonal scaling by \(a m(x)\) on rows and \(m(y)\) on columns preserves the rank of \(G[Z,T]=A_ZC_T\). The scaled entries are

\[
\begin{aligned}
\widetilde G(x,y)
&=\frac{\kappa_x m(y)^2+\kappa_y m(x)^2}{x+y}\\
&=\frac{(\kappa_x+\kappa_y)m(x)^2}{x+y}
 +\kappa_x\frac{(m(x)+m(y))^2}{x+y},
\end{aligned}
\tag{9}
\]

where \(\kappa_x=\kappa_i\) for \(x\in H_i\). The first term is zero within a half. Across halves, it is a nonzero row scaling of a Cauchy matrix. With

\[
e_i=|Z\cap H_i|,\qquad t_i=|T\cap H_i|,
\]

its rank is

\[
k_{\rm cross}=\min(e_0,t_1)+\min(e_1,t_0).
\]

Both rectangular Cauchy blocks have full rank: their row and column nodes are distinct and disjoint. The second term of (9) has rank at most \(2d\). Indeed, \((m(X)+m(Y))^2/(X+Y)\) is a polynomial of total degree at most \(2d-1\). Its evaluation matrix is a sum of at most \(2d\) rank-one matrices, one for each power of \(X\). Therefore

\[
\operatorname{rank}_K A_ZC_T\ge k_{\rm cross}-2d.
\tag{10}
\]

This is not a generic-transversality claim. If \(Z,T\) lie in the same half, the Cauchy term vanishes and the cross-rank is at most \(2d\). At \(s=30\), that upper bound is four, even for disjoint eight-element row and column sets. A blanket lower bound \(\min(e,t)\) is therefore false for the current maps.

## An occupancy-only return cap

For \(e<s\), define

\[
\begin{aligned}
k_0&=\max\{0,\ e+u-s,\ \min(e_0,t_1)+\min(e_1,t_0)-2d\},\\
\eta&=\left\lfloor\frac{\min(c,64-c,s-k_0)}{16}\right\rfloor,\\
\underline r&=\max\{0,\ e+u-s,\ k_0-\eta\},\\
\overline\delta&=u-\underline r.
\end{aligned}
\tag{11}
\]

When \(k_0=0\), the last expression involving \(k_0-\eta\) cannot improve the trivial bound; (7) need not be invoked. Otherwise, (10) and the ordinary Sylvester rank bound certify a nonzero minor of size \(k_0\le\min(e,u)\). Equations (3) and (7) then justify (11). In particular,

\[
0\le\overline\delta\le\min(u,s-e),\qquad
\Pr[R\ne0,Y_D=Y_Z=0]
\le\frac{(q^{\overline\delta}-1)q^{t-u}}{q^j}.
\tag{12}
\]

To obtain a cap depending only on \(j\), first fix a support split \(j_i=|S\cap H_i|\), with \(j_0+j_1=j\) and \(0\le j_i\le32\). For counts

\[
0\le c_i\le j_i,\qquad 0\le e_i\le32-j_i,
\]

put \(c=c_0+c_1\), \(e=e_0+e_1\), \(t_i=j_i-c_i\), and compute (11). There are exactly

\[
M(j_0,j_1;c_0,c_1,e_0,e_1)
=\prod_{i=0}^1\binom{j_i}{c_i}\binom{32-j_i}{e_i}
\]

choices of \(D,Z\) with these counts. Expanding

\[
z^{\operatorname{wt}(Y)}
=\prod_{x=0}^{63}\bigl(z+(1-z)\mathbf1_{Y_x=0}\bigr)
\]

and using (12) yields the positive cap

\[
U_{j_0,j_1}(z)=
\sum_{\substack{0\le c_i\le j_i\\0\le e_i\le32-j_i\\e_0+e_1<s}}
M(j_0,j_1;c_0,c_1,e_0,e_1)
z^{64-c-e}(1-z)^{c+e}
\frac{(q^{\overline\delta}-1)q^{t-u}}{q^j}.
\tag{13}
\]

Terms with \(e\ge s\) vanish. Terms with \(c=j\) also vanish because \(u=0\). Every term is nonnegative, and no cancellation or heuristic rank assumption is used.

For every support \(S\) with the fixed split, \(J_S(z)\le U_{j_0,j_1}(z)\). Thus a uniform occupancy-only replacement is

\[
J_j^{\rm new}(z)=\min\left\{J_j^{\rm old}(z),\quad
\max_{\max(0,j-32)\le j_0\le\min(32,j)}U_{j_0,j-j_0}(z)\right\}.
\tag{14}
\]

Here \(J_j^{\rm old}\) is the current returning-state cap. The maximum in (14) is essential: the proof has fixed the support, not sampled its half split. No support split is averaged, and no route randomness enters this step. A later refinement could average the split through a larger route coefficient table, but (14) does not require that change.

## Scope and next step

Only the returning-state entry changes, from \(J_j^{\rm old}/(q^s-1)\) to \(J_j^{\rm new}/(q^s-1)\). The quiet, birth, and live entries, fractional exponent, all-pair composition, support counts, and route remain unchanged. The penalty \(2d\) is 16, eight, and four for states 24, 28, and 30, respectively.

This note supplies an algebraic candidate cap, not numerical evidence that it closes any distance gap. No provider, receipt, gate, benchmark, or fixed map is changed. Before numerical use, an independent review should check the minor multiplicity argument and a tiny exact oracle should verify (11)--(13). Any implementation must retain exact multiplicities and directed arithmetic where a certificate is claimed.
