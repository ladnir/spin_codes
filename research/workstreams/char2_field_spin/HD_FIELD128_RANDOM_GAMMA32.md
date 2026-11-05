# Optional random-gamma variant at state dimension 32

This note describes a future construction variant, not the fixed-gamma construction in current receipts. It gives a setup-failure bound, not a distance certificate. No implementation, selected-weight gate, or whole-code cover is supplied here.

## Construction and objective

Let \(E=\operatorname{GF}(256)\), \(F=\operatorname{GF}(2^{128})\), and \(q=|F|\). Let \(A\in E^{64\times32}\) be the Lagrange evaluation map from nodes \(128,\ldots,159\) to nodes \(0,\ldots,63\). The evaluation nodes form a 64-element additive subspace, so \(A^\mathsf TA=0\). Each row of \(A\) sums to one.

Sample one global scalar \(\gamma\), independently of the outer randomizers, route, and refresh diagonals. Fix this scalar for the complete encoder. Replace the feedback map by \(\gamma A^\mathsf T\):

\[
Y=X+AQ,\qquad
Q_{\rm next}=D_2H_{32}D_1(Q+\gamma A^\mathsf TX).
\]

The outer encoder and paired route remain unchanged. Each refresh uses fresh independent nonzero full-field diagonal entries. The initial state is zero; the final state is unrestricted.

The desired setup property is

\[
\mathcal G(\gamma):\quad
\det(I_{|D|}+\gamma A_DA_D^\mathsf T)\ne0
\quad\text{for every }D\subseteq\{0,\ldots,63\},
\]

where \(A_D\) contains the rows indexed by \(D\). This property supports the joint return count from [the fixed-gamma note](HD_FIELD128_GAMMA_RETURN.md). Unlike states at most 30, state 32 does not follow deterministically from the degree-16 condition alone.

## One polynomial for each complementary pair

Fix a coordinate set \(D\), and write \(c=|D|\) and \(M_D=A_DA_D^\mathsf T\). The matrix \(M_D\) is symmetric with diagonal one. In characteristic two, inverse-cycle cancellation in the determinant leaves only matchings. Consequently,

\[
p_D(t):=\det(I_c+tM_D)
=\begin{cases}
R_D(t)^2,&c\text{ even},\\
(1+t)R_D(t)^2,&c\text{ odd},
\end{cases}
\qquad R_D(0)=1.
\]

The matching expansion establishing this identity is given in the fixed-gamma note. Sylvester's determinant identity bounds \(\deg p_D\le32\). Thus \(\deg R_D\le16\) for even \(c\), and \(\deg R_D\le15\) for odd \(c\).

These polynomials also agree on complementary coordinate sets. Over the rational-function field \(E(t)\), define

\[
P(t)=I_{64}+tAA^\mathsf T.
\]

The identity \(A^\mathsf TA=0\) gives \(P(t)^2=I_{64}\). Sylvester's identity gives \(\det P(t)=1\). Jacobi's complementary-minor identity therefore yields

\[
p_D(t)=\det P(t)\det(P(t)^{-1})_{D^c,D^c}=p_{D^c}(t).
\]

There are \(2^{63}\) complementary pairs of subsets. The pair consisting of the empty and complete sets has polynomial one. Hence at most \(2^{63}-1\) nontrivial polynomials need a root bound.

## Two admissible sampling rules

First, sample \(\gamma\) uniformly from \(F\setminus\{0,1\}\). After excluding \(1\), each nontrivial polynomial has at most 16 distinct roots in this sampling set. The union bound gives

\[
\Pr[\neg\mathcal G(\gamma)]
\le\frac{16(2^{63}-1)}{q-2}
<2^{-61}.
\]

Alternatively, sample \(\gamma\) uniformly from the elements of degree exactly 16 over \(E\). Every proper divisor of 16 divides 8. Thus this sampling set is

\[
\Gamma_{16}=F\setminus\operatorname{GF}(2^{64}),
\qquad |\Gamma_{16}|=q-2^{64}.
\]

For odd \(|D|\), the nonzero polynomial \(R_D\) has degree at most 15 and cannot vanish on \(\Gamma_{16}\). Only even sets remain. Complementation preserves their parity, leaving \(2^{62}-1\) nontrivial complementary pairs. Therefore

\[
\Pr_{\gamma\gets\Gamma_{16}}[\neg\mathcal G(\gamma)]
\le\frac{16(2^{62}-1)}{q-2^{64}}
<2^{-62}.
\]

These are bounds for one global setup event. They must not be charged separately for each source support, packet, or region.

## Use in a later distance proof

For every fixed \(\gamma\) satisfying \(\mathcal G(\gamma)\), the joint return count holds uniformly over active supports and output-zero sets. If a later proof bounds the distance-failure probability by \(\varepsilon\) for every such fixed scalar, then

\[
\Pr[\text{distance failure}]
\le \Pr[\neg\mathcal G(\gamma)]+\varepsilon.
\]

The setup term is added once to the final union bound. Existing fixed-gamma receipts do not establish the conditional state-32 bound \(\varepsilon\); fresh local coefficients and coverage are required.

A generally sampled scalar also loses the sparse multiplication rule available for the fixed tower generator. The implementation would need to account for applying that scalar to the feedback vector. No performance transfer from the fixed-generator implementation is claimed.
