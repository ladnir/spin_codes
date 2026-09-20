# Pairwise-systematic constituent: finite certificate

## Result

This certificate replaces the uniform 256-by-512 outer matrix with a
512-bit sampler and an algebraic encoder. Setup samples two elements of
\(\mathbb F_{2^{256}}\). The same sampled constituent is used in all outer
rows.

Set

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\lceil0.109N\rceil=228{,}590.
\]

For the construction below,

\[
 \Pr[d_{\min}<D]
 <2^{-42.2593129905}<2^{-40}. \tag{1}
\]

The probability in (1) includes the outer constituent, both permutation
layers, and every inner step map. Setup performs no rejection test.

## Outer constituent

Represent a 256-bit row as an element of

\[
 \mathbb F:=\mathbb F_2[X]/
 (X^{256}+X^{10}+X^5+X^2+1).
\]

The verifier applies Rabin's criterion to the displayed polynomial. Setup
samples independent

\[
 a,b\gets\mathbb F.
\]

For every row, the outer encoder computes

\[
 \mathsf{Enc}_{a,b}(u)
 :=\bigl(u,au+bu^2\bigr)
 =\bigl(u,u(a+bu)\bigr)\in\mathbb F_2^{512}. \tag{2}
\]

Frobenius squaring and multiplication by a fixed field element are
\(\mathbb F_2\)-linear. Hence (2) is a binary-linear map. Its first half is
\(u\), so every sampled map has dimension 256.

The sampler consumes 512 unbiased bits. Evaluation uses two
\(\mathbb F_{2^{256}}\) multiplications and one field addition per row. This
is an operation count, not a timing claim.

## Pairwise-independence lemma

Fix distinct nonzero \(u,v\in\mathbb F\). The two parity halves satisfy

\[
 \begin{pmatrix}
  au+bu^2\\
  av+bv^2
 \end{pmatrix}
 =
 \begin{pmatrix}
  u&u^2\\
  v&v^2
 \end{pmatrix}
 \begin{pmatrix}a\\b\end{pmatrix}. \tag{3}
\]

The determinant of the matrix in (3) is

\[
 uv^2+vu^2=uv(u+v)\ne0. \tag{4}
\]

Thus (3) is a bijection on \(\mathbb F^2\). Under the setup distribution,
the parity halves for \(u\) and \(v\) are independent and uniform in
\(\mathbb F\). This statement holds for every fixed pair of distinct
nonzero messages; it is not an asymptotic approximation.

## Spectrum event

For \(w\in\{1,\ldots,512\}\), define

\[
 A_w(a,b):=
 \left|\left\{u\in\mathbb F\setminus\{0\}:
       \operatorname{wt}(u)+\operatorname{wt}(au+bu^2)=w
 \right\}\right|.
\]

For a fixed message \(u\) of weight \(i\), its parity half is uniform.
Therefore

\[
 \mu_w:=\mathbb E_{a,b}[A_w(a,b)]
 =2^{-256}
 \begin{cases}
  \binom{512}{w}-\binom{256}{w},&w\le256,\\
  \binom{512}{w},&w>256.
 \end{cases} \tag{5}
\]

Equation (5) follows from Vandermonde's identity after removing the zero
message. The verifier checks this identity for every shell.

Let \(I_{u,w}\) indicate that message \(u\) produces weight \(w\). Equations
(3)--(4) make \(I_{u,w}\) and \(I_{v,w}\) independent for distinct nonzero
\(u,v\). Consequently,

\[
 \operatorname{Var}_{a,b}(A_w)
 =\sum_{u\ne0}\operatorname{Var}(I_{u,w})
 \le\mu_w. \tag{6}
\]

Let \((T_w)_{w=1}^{512}\) be the integer cap vector from the uniform-matrix
certificate. Its nonzero support is \(42\le w\le470\). Define

\[
 \mathcal E:=\{A_w(a,b)\le T_w\text{ for every }w\}. \tag{7}
\]

For a zero cap, the verifier applies Markov's inequality. For a positive
cap, it applies Cantelli's inequality using (5)--(6). Exact rational
arithmetic gives

\[
 \Pr_{a,b}[\neg\mathcal E]
 \le1.896832534308249\cdot10^{-13}
 <2^{-42.2614729204}. \tag{8}
\]

Event \(\mathcal E\) is only a proof event. The sampler does not evaluate
the shell counts.

## Routing and inner encoder

View the complete message as a matrix
\(U\in\mathbb F_2^{4096\times256}\). Apply (2) to every row using the same
sampled pair \((a,b)\).

For each row, setup independently samples a uniform permutation of its 512
coordinates. The encoder then transposes the 4096-by-512 array. For each of
the 512 resulting regions, setup independently samples a uniform
permutation of its 4096 positions.

The inner encoder is RandomStepConv with memory \(M=22\). It starts from
\(\sigma_0=0\in\mathbb F_2^{22}\). For each
\(t\in\{0,\ldots,N-1\}\), setup independently samples

\[
 H_t\gets\mathbb F_2^{23\times23}.
\]

For routed input bit \(x_t\), the encoder computes

\[
 (y_t,\sigma_{t+1}):=H_t(x_t,\sigma_t)
\]

and discards \(\sigma_N\). Every setup object is sampled once and shared by
all messages. Thus each setup defines one binary-linear code.

## Conditional transfer

Fix any dimension-256 constituent whose shell counts satisfy (7). The
sparse and dense outward verifiers use only the caps \((T_w)\); they do not
use the distribution that produced the constituent.

The sparse verifier covers occupations \(1\le Q\le159\). The dense convex
cover handles every integer composition for \(160\le Q\le4096\). Together
they prove

\[
 \Pr[d_{\min}<D\mid\mathcal E]
 <2^{-51.6439589890}. \tag{9}
\]

The probability in (9) is over the routing permutations and inner step maps
after fixing \((a,b)\in\mathcal E\).

Combining (8) and (9) gives

\[
 \Pr[d_{\min}<D]
 \le \Pr[\neg\mathcal E]
    +\Pr[d_{\min}<D\mid\mathcal E]
 <2^{-42.2593129905},
\]

which proves (1).

## Scope

This is a finite theorem for \(k=2^{20}\). It uses the information-theoretic
RandomStepConv-M22 inner ensemble, not RM2Sub. The outer sampler and outer
map are compact, but this certificate makes no compact-sampler claim for the
inner step maps. An asymptotic family requires a separate parameter
schedule.
