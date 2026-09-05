# Finite repeated-random-outer RandomStepConv certificate

## Status

This document states a complete finite mathematical certificate for the
RandomStepConv comparison model. The construction uses one sampled
\([512,256]\) outer constituent and repeats that same constituent in every
outer row. It does not sample a different outer code per row.

The result proves a literal relative-distance target of $0.11$ at message
dimension $2^{20}$. It also supplies a bounded outer-selection procedure
whose abort probability is included in the 40-bit failure budget.

The result is not an implementation claim. The outer acceptance test and the
random inner are intentionally strong proof-model objects.

## Construction

Set

\[
 B=512,
 \qquad K=256,
 \qquad L=4140,
 \qquad L_0=4096,
\]

and hence

\[
 k=K L_0=2^{20},
 \qquad N=B L=2{,}119{,}680.
\]

The parent outer code has $K L=1{,}059{,}840$ input coordinates. The
construction shortens its final $L-L_0=44$ complete information rows to
zero.

For a binary $K\times B$ matrix $G$, define the message-multiplicity
enumerator

\[
 A_w(G)=\bigl|\{u\in\mathbb F_2^K\setminus\{0\}:
                 \operatorname{wt}(uG)=w\}\bigr|.
\]

When $G$ has full row rank, this is the ordinary nonzero weight spectrum of
its row-space code. Call $G$ admissible when it has row rank $K$ and

\[
 A_w(G)\le \frac{13}{2}\,(2^K-1)\frac{\binom Bw}{2^B}
 \tag{1}
\]

for every $1\le w\le B$.

The bounded setup procedure samples at most 28 independent uniform binary
$K\times B$ matrices. It retains the first admissible matrix and aborts if
none is admissible. The retained matrix is the single outer constituent used
in all $L$ rows.

For each outer row, sample an independent uniform permutation of its $B$
codeword coordinates. For each coordinate $j\in\{1,\ldots,B\}$, sample an
independent uniform permutation of the $L$ row positions. Concatenate the
$B$ permuted regions, each of length $L$, to obtain the $N$-bit inner
input.

The RandomStepConv inner has state $s_t\in\mathbb F_2^{30}$, initial state
$s_0=0$, and input bit $x_t$. At each position
$t\in\{0,\ldots,N-1\}$, sample an independent uniform binary
$31\times31$ matrix $M_t$ and set

\[
  (y_t,s_{t+1})=M_t(x_t,s_t).
\]

The encoder emits $y_0,\ldots,y_{N-1}$ and discards the terminal state. All
outer permutations and all matrices $M_t$ are sampled once and are shared
by every message. Consequently, each completed setup defines one linear
map from $\mathbb F_2^k$ to $\mathbb F_2^N$.

## Certified theorem

Let

\[
 D=233{,}165=\left\lceil0.11N\right\rceil.
\]

For every fixed admissible outer matrix $G$, the probability over the row
permutations, region permutations, and RandomStepConv matrices satisfies

\[
 \Pr[d_{\min}<D]<2^{-108}.
 \tag{2}
\]

One uniform outer sample is admissible with probability greater than

\[
 0.6365388449160215.
 \tag{3}
\]

For the bounded 28-attempt setup, the probability of setup abort or
$d_{\min}<D$ is less than $2^{-40}$. Outside this combined failure event,
the resulting code has

\[
 \frac{d_{\min}}N\ge
 \frac{233{,}165}{2{,}119{,}680}
 =0.11000009435386474\ldots.
\]

The probability in the combined statement includes every sampled object:
the at most 28 outer candidates, the retained constituent's row-coordinate
permutations, the region permutations, and all RandomStepConv matrices.

## Outer-selection bound

For a uniform binary $K\times B$ matrix and any fixed nonzero message, the
corresponding $B$-bit word is uniform. Therefore

\[
 \mu_w:=\mathbb E[A_w]
 =(2^K-1)\frac{\binom Bw}{2^B}.
\]

Distinct nonzero binary messages are linearly independent. Their codewords
are therefore pairwise independent, which gives

\[
 \operatorname{Var}(A_w)\le \mu_w.
\]

When $\mu_w<2/13$, failure of (1) implies $A_w\ge1$, so Markov's
inequality gives failure probability at most $\mu_w$. For the remaining
weights, Cantelli's inequality gives

\[
 \Pr\!\left[A_w>\frac{13}{2}\mu_w\right]
 \le
 \frac{1}{1+(11/2)^2\mu_w}.
\]

A union bound over all weights, together with the usual sequential row-rank
bound, gives

\[
 \Pr[G\text{ is not admissible}]
 <0.3634611550839785.
 \tag{4}
\]

The verifier evaluates the sum in (4) as an exact rational number. It then
checks directly that its 28th power, plus the conservative conditional term
$2^{-108}$, is less than $2^{-40}$.

## Distance bound

Fix an admissible $G$. After a uniform coordinate permutation, (1) bounds
the probability mass of each nonzero row word by the pointwise envelope

\[
 \Lambda\,2^{-B},
 \qquad
 \Lambda=\frac{13}{2}(2^K-1).
 \tag{5}
\]

Thus, for a message with $Q$ nonzero outer rows, its $Q$ permuted row
words are bounded by $\Lambda^Q$ times $Q$ independent uniform $B$-bit
reference words. This is a pointwise comparison, so it remains valid inside
the causal inner transfer.

For $z\in(0,1)$, put $q=2^{-30}$ and $b=(1+z)/2$. The exact tilted
RandomStepConv transfers, with states `zero` and `nonzero`, are

\[
 T_0(z)=
 \begin{pmatrix}
 1&0\\
 qb&(1-q)b
 \end{pmatrix},
 \qquad
 T_1(z)=
 \begin{pmatrix}
 qb&(1-q)b\\
 qb&(1-q)b
 \end{pmatrix}.
 \tag{6}
\]

For the uniform reference row, a candidate input position has transfer

\[
 C(z)=\frac{T_0(z)+T_1(z)}2.
\]

Let $R_Q(z)$ be the average ordered product over all placements of $Q$
candidate positions among one region's $L$ positions:

\[
 R_Q(z)=\binom LQ^{-1}
 \sum_{\substack{S\subseteq\{1,\ldots,L\}\\|S|=Q}}
 \prod_{r=1}^L
 \begin{cases}
 C(z),&r\in S,\\
 T_0(z),&r\notin S.
 \end{cases}
 \tag{7}
\]

The region permutations give exactly this uniform-subset law. The $B$
reference coordinates are independent, so the complete tilted moment is

\[
 e_0^{\mathsf T}R_Q(z)^B\mathbf 1.
\]

The calculation deliberately retains the globally all-zero reference tuple.
This only enlarges the upper bound and keeps every verifier operation
positive.

For every $1\le Q\le L_0$, define

\[
 U_Q=
 \binom{L_0}{Q}\Lambda^Q
 \inf_{s>0}
 \left(
 e^{Ds}
 e_0^{\mathsf T}R_Q(e^{-s})^B\mathbf 1
 \right).
 \tag{8}
\]

Equation (8) upper-bounds the expected number of nonzero messages with
output weight below $D$ and occupation $Q$. The outward verifier checks
fixed dyadic witnesses for all 4096 occupations and obtains

\[
 \max_{1\le Q\le4096}\log_2 U_Q
 <-120.16539169026359.
\]

The maximum occurs at $Q=1$. Hence

\[
 \sum_{Q=1}^{4096}U_Q
 \le4096\max_Q U_Q
 <2^{-108.16539169026358}.
\]

Markov's inequality proves (2).

## Reproducibility

The nearest-binary64 search is in
`evaluate_repeated_random512_randomstepconv_g1.py`. Its receipt is
`repeated_random512_randomstepconv_g1_s30_allq_d11.json`.

The accepting verifier is
`certify_repeated_random512_randomstepconv_g1_outward.py`. Its receipt is
`repeated_random512_randomstepconv_g1_s30_allq_outward_d11.json`.
`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json` hash-binds the
theorem, search, verifier, and both receipts.

The search receipt supplies only fixed tilt witnesses. The accepting verifier
recomputes every inequality. It uses Arb enclosures for transcendental
endpoints, exact rational arithmetic for the outer-selection probability,
and positive binary64 recurrences advanced upward after every elementary
operation. Per-degree power-of-two scaling prevents underflow from deleting
positive mass.

Run

```powershell
python workstreams/finite_asymptotic_theory/certify_repeated_random512_randomstepconv_g1_outward.py
```

from the repository root.

## What remains open

The following statements are not proved here.

- The exact admissibility test is efficient. A direct test may enumerate all
  $2^{256}$ row messages and is therefore not an implementation path.
- RandomStepConv has competitive storage or encoding time. Storing all
  $31\times31$ step matrices requires about 242.83 MiB before routing and
  outer-code data.
- RM2Sub-S19 satisfies the RandomStepConv transfer bound. Its fixed
  low-rank syndrome interface requires a different argument.
- A structured, efficiently testable outer family satisfies (1), or a
  comparably useful envelope, with the required setup guarantee.

The certificate answers the immediate baseline question: at these finite
dimensions, a bounded-memory random linear convolution is strong enough.
The next proof task is therefore to recover the relevant transfer property
with an implementable structured inner, rather than to strengthen this
random inner further.
