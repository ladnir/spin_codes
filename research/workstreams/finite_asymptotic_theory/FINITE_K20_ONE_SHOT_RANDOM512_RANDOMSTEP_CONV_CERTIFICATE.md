# One-shot random-constituent finite certificate

## Result

This certificate removes the infeasible spectrum acceptance test from the
earlier repeated-random-outer construction. Setup samples one uniform
binary (256)-by-(512) matrix. The encoder repeats that same matrix in all
outer rows.

Set

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\lceil0.109N\rceil=228{,}590.
\]

For the construction below,

\[
 \Pr[d_{\min}<D]<2^{-42.2593129905}<2^{-40}.
 \tag{1}
\]

The probability in (1) includes the outer matrix, both permutation layers,
and every inner step map. No setup test, rejection sampler, or codeword
enumeration occurs.

## Construction

Let (B=512), (K=256), and (L=4096). View a message as a matrix
(U\in\mathbb F_2^{L\times K}). Setup samples

\[
 G\gets\mathbb F_2^{K\times B}.
\]

The outer encoder maps row (U_i) to (U_iG). Every row uses the same
sampled matrix (G).

For each row \(i\), setup independently samples a uniform permutation
\(\rho_i\in S_B\) and permutes that row's coordinates. The encoder then
transposes the \(L\)-by-\(B\) array. For each of the \(B\) resulting regions,
setup independently samples a uniform permutation
\(\pi_j\in S_L\). Concatenating the permuted regions produces
\(x\in\mathbb F_2^N\). This factorization is one permutation of the complete
outer output.

The inner encoder is RandomStepConv with memory (M=22). Set
(sigma_0=0\in\mathbb F_2^M). For each position (t\in\{0,\ldots,N-1\}),
setup independently samples

\[
 H_t\gets\mathbb F_2^{(M+1)\times(M+1)}.
\]

The encoder computes

\[
 (y_t,\sigma_{t+1})=H_t(x_t,\sigma_t)
\]

and returns (y=(y_0,\ldots,y_{N-1})). The encoder discards
(sigma_N). All setup objects are sampled once and shared by every message.
Thus each setup defines one linear code.

## The spectrum event

For a matrix (G), define

\[
 A_w(G):=
 \left|\left\{u\in\mathbb F_2^K\setminus\{0\}:
                 \operatorname{wt}(uG)=w\right\}\right|.
\]

For every (w\in\{1,\ldots,B\}), let

\[
 \mu_w=(2^K-1)\binom Bw2^{-B}.
\]

The verifier constructs an integer cap (T_w\ge0). If
(mu_w\le2^{-51}), it sets (T_w=0). Otherwise it chooses (T_w) so that

\[
 \frac{\mu_w}{\mu_w+(T_w+1-\mu_w)^2}\le2^{-51}.
 \tag{2}
\]

Let (mathcal E) be the event that (G) has row rank (K) and
(A_w(G)\le T_w) for all (w). The nonzero caps have support
(42\le w\le470).

For a uniform (G), the image of each fixed nonzero (u) is uniform in
(mathbb F_2^B). Any two distinct nonzero binary messages are linearly
independent. Their images are therefore independent, and

\[
 \operatorname{Var}(A_w)\le\mu_w.
\]

Markov's inequality handles shells with (T_w=0). Cantelli's inequality and
(2) handle the other shells. A union bound also includes rank failure. The
exact rational calculation gives

\[
 \Pr_G[\neg\mathcal E]
 \le 1.896832534308249\mathbin\cdot10^{-13}
 <2^{-42.2614729204}.
 \tag{3}
\]

Event (mathcal E) is only a proof event. The sampler does not evaluate it.

## Inner transfer and monotonicity

Fix (z\in(0,1)), put (q=2^{-M}), and set (b=(1+z)/2). Classify the
inner state as zero or live. The exact tilted transfers for a fixed zero or
one input bit are

\[
 T_0(z)=
 \begin{pmatrix}
 1&0\\ qb&(1-q)b
 \end{pmatrix},
 \qquad
 T_1(z)=
 \begin{pmatrix}
 qb&(1-q)b\\ qb&(1-q)b
 \end{pmatrix}.
 \tag{4}
\]

For a Bernoulli-(p) input bit, define

\[
 T_p(z)=(1-p)T_0(z)+pT_1(z).
\]

The tilted moment is nonincreasing in every input probability (p). To see
this, let (v=(v_0,v_1)^\mathsf T) satisfy (v_0\ge v_1\ge0), and put

\[
 h=b(qv_0+(1-q)v_1).
\]

Then (h\le v_0) and

\[
 T_pv=((1-p)v_0+ph,h)^\mathsf T.
\]

Thus (T_p) preserves the cone (v_0\ge v_1\ge0), and

\[
 \frac{\partial}{\partial p}T_pv=(h-v_0,0)^\mathsf T\le0.
\]

Backward induction through a product of transfers proves the stated
coordinatewise monotonicity. This lemma lets the proof dominate both
spectrum tails by one low-density reference.

## Conditional distance bound

Fix any (G\in\mathcal E). Divide its nonzero weights into

\[
 [42,79],\qquad[80,432],\qquad[433,470].
\]

For a band (mathcal W) and (p\in(0,1)), define the exact rational
majorant

\[
 c(\mathcal W,p)=
 \max_{w\in\mathcal W}
 \frac{T_w}{\binom Bw p^w(1-p)^{B-w}}.
 \tag{5}
\]

Use (p_-=79/512) for the low band, (p_0=1/2) for the central band, and
(p_+=433/512) for the high band. Symmetry gives equal tail majorants.
Monotonicity replaces the high-density transfer by the transfer at (p_-).
Summing the low and high labels gives one defect category with majorant
(2c([42,79],p_-)). The central category uses
(c([80,432],p_0)).

Let (Z_Q) count nonzero messages with exactly (Q) nonzero outer rows and
output weight below (D). The proof bounds
(mathbb E[Z_Q\mid G]), where the expectation is over both permutation
layers and RandomStepConv.

For (Q=1), the verifier retains every exact shell cap (T_w). For
(Q=2), it retains every ordered shell pair and uses the valid count bound
(T_aT_b). These outward bounds are

\[
 \log_2\mathbb E[Z_1\mid G]\le-67.9441013828,
\]

\[
 \log_2\mathbb E[Z_2\mid G]\le-105.3577609749.
\]

For (3\le Q\le159), the verifier fixes each row category once. Independent
region permutations randomize only the order of that same category multiset.
A two-dimensional coefficient recurrence therefore avoids resampling the
category in each region. It gives

\[
 \log_2\sum_{Q=3}^{159}\mathbb E[Z_Q\mid G]
 \le-51.6796100317.
 \tag{6}
\]

For (160\le Q\le L), a categorical reference bounds the conditioning cost
inside each region. For a fixed reference, the unclipped exponent is convex
in the three type counts. Where that exponent is negative, the final log
bound is also convex. Vertex checks therefore certify each triangle in a
simplex cover.

The outward verifier reconstructs 884 accepted triangles exactly. It also
reconstructs 37 residual triangles and verifies their 300 integer lattice
points individually. Every composition has a pointwise bound of (2^{-80}).
There are 8,381,873 dense compositions, so

\[
 \log_2\sum_{Q=160}^{4096}\mathbb E[Z_Q\mid G]
 \le -57.0011587683.
 \tag{7}
\]

Combining all occupations gives

\[
 \log_2\mathbb E\left[\sum_{Q=1}^{L}Z_Q\mid G\right]
 \le-51.6439589890.
 \tag{8}
\]

Markov's inequality converts (8) to a conditional probability bound. Adding
that bound to (3) proves (1).

## Verification

The dense verifier is
`certify_single_random_constituent_dense_outward.py`. It uses 192-bit Arb
intervals and exact rational geometry. Its receipt is
`single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json`.

The sparse verifier is
`certify_single_random_constituent_sparse_outward.py`. It uses Arb for
transcendental endpoints and positive binary64 recurrences with one-ULP
upward rounding after each operation. Power-of-two scaling prevents
underflow from deleting positive mass. Its receipt is
`single_random_constituent_B512_sparse_outward_s22.json`.

The combined verifier is
`certify_single_random_constituent_combined.py`. Its receipt is
`single_random_constituent_B512_combined_outward_s22.json`.

Run, in order,

```powershell
python workstreams/finite_asymptotic_theory/certify_single_random_constituent_dense_outward.py
python workstreams/finite_asymptotic_theory/certify_single_random_constituent_sparse_outward.py
python workstreams/finite_asymptotic_theory/certify_single_random_constituent_combined.py
```

## Scope

This theorem is a finite certificate for the RandomStepConv comparison
ensemble. It does not prove the same transfer bound for RM2Sub. It also does
not claim a competitive implementation of the random inner. Explicitly
storing (N) independent (23)-by-(23) matrices requires substantial
space.

The outer sampler itself is efficient: it samples 131,072 independent bits
once. The theorem neither biases that sampler nor tests its output. A compact
structured replacement must reproduce the certified spectrum and transfer
properties through a different argument.
