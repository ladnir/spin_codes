# EBCH-direct-sum outer with a linearized mixer: finite certificate

## Result

Let \(C_0\subseteq\mathbb F_2^{512}\) be the direct sum of four fixed
extended BCH \([128,64,22]\) codes. Apply one sampled linearized field map
to every word of \(C_0\), and repeat the resulting constituent in all 4,096
outer rows.

Set

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\lceil0.109N\rceil=228{,}590.
\]

With the structured SPIN permutation and RandomStepConv-M22 inner,

\[
 \Pr[d_{\min}<D]<2^{-42.2593129905}<2^{-40}. \tag{1}
\]

The construction uses no accumulator stage inside the outer constituent.
The linearized mixer is an additional randomized layer; therefore (1) is not
a certificate for the BA-only construction.

## Linearized mixer

Identify \(\mathbb F_2^{512}\) with

\[
 \mathbb F:=\mathbb F_2[X]/
 (X^{512}+X^8+X^5+X^2+1).
\]

The verifier applies Rabin's criterion to the displayed polynomial. Setup
samples independent \(a,b\gets\mathbb F\) and defines

\[
 L_{a,b}(x):=ax+bx^2=x(a+bx). \tag{2}
\]

The map in (2) is \(\mathbb F_2\)-linear. The complete local encoder maps a
256-bit row \(u\) to

\[
 L_{a,b}(\mathsf{EBCH}^{\oplus4}(u)). \tag{3}
\]

Setup samples 1,024 bits for \((a,b)\). After the fixed EBCH encoding,
evaluation of (3) uses two \(\mathbb F_{2^{512}}\) multiplications and one
field addition. This is an operation count, not a timing claim.

## Pairwise independence and rank

Fix distinct nonzero \(x,y\in\mathbb F\). The coefficient matrix that maps
\((a,b)\) to \((L_{a,b}(x),L_{a,b}(y))\) has determinant

\[
 xy^2+yx^2=xy(x+y)\ne0. \tag{4}
\]

Thus the two images are independent and uniform in \(\mathbb F\). This holds
for every fixed pair of distinct nonzero words in \(C_0\).

The map \(L_{a,b}\) can have a one-dimensional kernel. Let
\(M=|C_0|=2^{256}\) and \(q=|\mathbb F|=2^{512}\). The restriction of
\(L_{a,b}\) to \(C_0\) fails to be injective with exact probability

\[
 p_{\mathrm{rank}}
 =\frac{1+(q-1)(M-1)}{q^2}<2^{-256}. \tag{5}
\]

Indeed, \((a,b)=(0,0)\) gives the zero map. If \(a,b\ne0\), the nonzero
kernel element is \(a/b\), which is uniform in
\(\mathbb F\setminus\{0\}\). The cases with exactly one of \(a,b\) equal
to zero are invertible.

## Spectrum event

For \(w\in\{1,\ldots,512\}\), let \(A_w(a,b)\) count nonzero words in
\(L_{a,b}(C_0)\) with weight \(w\), with messages counted before a possible
collision. Equation (4) gives

\[
 \mu_w:=\mathbb E[A_w]
 =(2^{256}-1)\binom{512}{w}2^{-512},
 \qquad
 \operatorname{Var}(A_w)\le\mu_w. \tag{6}
\]

Let \((T_w)\) be the integer caps from the uniform-matrix certificate. Define
\(\mathcal E\) to require injectivity on \(C_0\) and
\(A_w(a,b)\le T_w\) for every \(w\). Markov's inequality handles zero caps.
Cantelli's inequality handles positive caps using (6). Exact rational
arithmetic and (5) give

\[
 \Pr_{a,b}[\neg\mathcal E]
 <2^{-42.2614729204}. \tag{7}
\]

The rank term in (7) is less than \(2^{-256}\); the shell caps dominate the
event failure. Setup does not test \(\mathcal E\).

## Conditional SPIN transfer

For each outer row, setup independently samples a uniform permutation of its
512 coordinates. The encoder transposes the 4096-by-512 array. Setup then
samples one independent uniform permutation of 4,096 positions in each of
the 512 regions.

The inner encoder is RandomStepConv-M22. At each of the \(N\) positions,
setup samples an independent binary 23-by-23 linear map. Every sampled object
is fixed and shared by all messages.

Condition on any \((a,b)\in\mathcal E\). The existing outward sparse and
dense verifiers depend only on the caps \((T_w)\). They cover every
occupation \(1\le Q\le4096\) and prove

\[
 \Pr[d_{\min}<D\mid\mathcal E]
 <2^{-51.6439589890}. \tag{8}
\]

The probability in (8) is over the routing permutations and inner maps.
Combining (7) and (8) proves (1).

## Scope

The proof of (4)--(8) applies to any fixed binary \([512,256]\) base code;
the EBCH direct sum is the concrete instance. The same proof permits any
number of invertible BA stages before the mixer, but those stages do not
improve the bound. Consequently, zero BA stages minimize this modified
outer.

This theorem is finite at \(k=2^{20}\). It concerns RandomStepConv-M22, not
RM2Sub. The outer encoder is linear-time in the number of fixed-size rows,
but the stated field-operation count has not been benchmarked.
