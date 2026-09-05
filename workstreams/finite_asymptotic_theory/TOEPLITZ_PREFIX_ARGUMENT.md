# The Toeplitz prefix argument

## Purpose

This note isolates the argument that made the random-convolution proof work.
The argument applies to every fixed binary linear outer code. It does not
depend on the outer weight spectrum or on independent sampling of outer
blocks.

The key property is stronger than ordinary mixing. For each fixed nonzero
outer word, its first nonzero coordinate exposes a triangular sequence of
fresh Toeplitz coefficients. The remaining output coordinates are therefore
independent fair bits.

This note states that property, derives the all-message prefix identity, and
records the limits of the argument. The concrete finite certificate is in
`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md`.

## Construction and probability space

Fix an integer \(N\ge 1\) and a binary linear code

\[
 V\subseteq\mathbb F_2^N.
\]

The code \(V\), its generator, its shortening rule, and all routing
permutations are fixed before the following sampling step.

Sample independent bits

\[
 h_1,\ldots,h_{N-1}\gets\mathbb F_2
\]

and set \(h_0:=1\). Define the lower-triangular Toeplitz map
\(T_h:\mathbb F_2^N\to\mathbb F_2^N\) by

\[
 (T_hu)_t
 :=
 \sum_{i=0}^{t}h_i u_{t-i},
 \qquad 0\le t<N.
 \tag{1}
\]

The same sampled vector \(h\) defines the map for every codeword. The
probability space contains only the \(N-1\) bits in \(h\), unless a later
statement explicitly adds outer-code randomness.

The matrix of \(T_h\) is lower triangular with unit diagonal. Hence \(T_h\)
is invertible for every realization of \(h\).

## Exact suffix law

Fix a nonzero word \(u\in\mathbb F_2^N\). Define its activation coordinate

\[
 \tau(u):=\min\{t:u_t=1\}.
\]

Write \(y:=T_hu\). The output before activation is zero, and the activation
output is one:

\[
 y_t=0\quad(t<\tau),
 \qquad
 y_\tau=1.
 \tag{2}
\]

For each \(r\ge1\) with \(\tau+r<N\), equation (1) gives

\[
 y_{\tau+r}
 =
 h_r+
 \sum_{i=0}^{r-1}h_i u_{\tau+r-i}.
 \tag{3}
\]

The second term in (3) depends only on \(u\) and
\(h_1,\ldots,h_{r-1}\). The new coefficient \(h_r\) appears with coefficient
\(u_\tau=1\).

### Lemma 1: Toeplitz suffix law

For every fixed nonzero \(u\in\mathbb F_2^N\) with activation coordinate
\(\tau\),

\[
 \bigl((T_hu)_{\tau+1},\ldots,(T_hu)_{N-1}\bigr)
 \mathrel{\overset{\mathrm d}=}
 \mathsf{Unif}\bigl(\mathbb F_2^{N-\tau-1}\bigr).
 \tag{4}
\]

Consequently,

\[
 \operatorname{wt}(T_hu)
 \mathrel{\overset{\mathrm d}=}
 1+\operatorname{Bin}(N-\tau-1,1/2).
 \tag{5}
\]

#### Proof

Fix \(u\). Consider the map

\[
 (h_1,\ldots,h_{N-\tau-1})
 \longmapsto
 (y_{\tau+1},\ldots,y_{N-1}).
\]

Equation (3) shows that its \(r\)-th output has the form

\[
 y_{\tau+r}=h_r+f_r(h_1,\ldots,h_{r-1})
\]

for a deterministic function \(f_r\). The map is triangular with unit
diagonal. It is therefore a bijection of
\(\mathbb F_2^{N-\tau-1}\). A bijection maps a uniform input to a uniform
output, which proves (4). Equation (2) then gives (5). \(\square\)

The word \(u\) can otherwise be arbitrary. After \(\tau\) is fixed, the
distribution in (5) does not depend on the remaining coordinates of \(u\).

The term *fresh coefficient* in this proof is algebraic. It means that
\(h_r\) has not appeared in an earlier suffix equation. The encoder does not
resample coefficients for different codewords or during each invocation.

## Counting all outer words by prefix ranks

Lemma 1 reduces the failure probability for a fixed word to its activation
coordinate. The outer code must therefore provide only the number of words
that activate at each coordinate.

For \(0\le t\le N\), define the shortened prefix subcode

\[
 V_t
 :=
 \{v\in V:v_0=\cdots=v_{t-1}=0\}
\]

and its dimension

\[
 \kappa_t:=\dim V_t.
\]

Set

\[
 Z_t:=|V_t\setminus\{0\}|=2^{\kappa_t}-1.
\]

The number of nonzero words with activation coordinate exactly \(t\) is

\[
 Z_t-Z_{t+1}.
 \tag{6}
\]

These counts follow from ranks of prefixes of a generator matrix. They do
not require enumeration of the \(2^{\dim V}\) codewords.

## Exact first-moment identity

Fix an integer threshold \(D\) with \(1\le D<N\). For a sampled \(h\), let

\[
 Z_{\mathrm{bad}}(h)
 :=
 \bigl|
 \{u\in V\setminus\{0\}:\operatorname{wt}(T_hu)\le D\}
 \bigr|.
\]

Define

\[
 F_t
 :=
 \Pr\bigl[
 \operatorname{Bin}(N-t-1,1/2)\le D-1
 \bigr].
\]

Lemma 1 and (6) give the exact identity

\[
 \mathbb E_h[Z_{\mathrm{bad}}]
 =
 \sum_{t=0}^{N-1}(Z_t-Z_{t+1})F_t.
 \tag{7}
\]

No independence between the bad-word events is used in (7). Linearity of
expectation applies even though all codewords share the same Toeplitz seed.

For computation, define

\[
 F_N^{\star}
 :=
 \Pr\bigl[\operatorname{Bin}(N-1,1/2)\le D-1\bigr].
\]

Summation by parts transforms (7) into

\[
 \mathbb E_h[Z_{\mathrm{bad}}]
 =
 Z_0F_N^{\star}
 +
 \sum_{t=1}^{N-D}
 Z_t
 \frac{\binom{N-t-1}{D-1}}{2^{N-t}}.
 \tag{8}
\]

Equation (8) is an all-message identity for the fixed code \(V\). Its inputs
are \(N\), \(D\), and the prefix dimensions \(\kappa_t\).

### Corollary 2: distance-failure bound

For every fixed binary linear code \(V\subseteq\mathbb F_2^N\),

\[
 \Pr_h\bigl[d_{\min}(T_h(V))\le D\bigr]
 \le
 \mathbb E_h[Z_{\mathrm{bad}}],
 \tag{9}
\]

where the expectation is given exactly by (7) or (8).

#### Proof

The event on the left of (9) is the event
\(Z_{\mathrm{bad}}\ge1\). Markov's inequality gives

\[
 \Pr_h[Z_{\mathrm{bad}}\ge1]
 \le \mathbb E_h[Z_{\mathrm{bad}}].
\]

The invertibility of \(T_h\) ensures that every nonzero input remains
nonzero. \(\square\)

If the right side of (9) is below \(2^{-\lambda}\), then one sampled
Toeplitz map fails the distance target with probability below
\(2^{-\lambda}\). If the right side is below one, the same calculation also
proves that at least one satisfactory Toeplitz map exists.

## Why the interface is unusually economical

A usual concatenated-code proof tracks outer weights, block occupations, or
local spectra. None of those statistics appears in (8). A fixed word's
post-activation suffix has the same distribution regardless of that word's
later support.

The proof therefore compresses the fixed outer code to one sequence:

\[
 (\kappa_0,\ldots,\kappa_N).
\]

This sequence can be computed by incremental Gaussian elimination on routed
generator columns. For a repeated constituent, the computation can exploit
the row decomposition. It still certifies all nonzero global messages.

The proof also explains why a typical fixed repeated BA code can perform
much better than its ensemble average. The ensemble average pays for rare
outer samples with poor late-prefix ranks. A certificate for one fixed code
uses its actual prefix ranks and pays no such averaging loss.

## Conditions that must remain explicit

The argument requires the following conditions.

1. The outer code and every routing permutation are fixed independently of
   the Toeplitz seed used in the probability claim.
2. The diagonal coefficient is fixed to \(h_0=1\).
3. The coefficients \(h_1,\ldots,h_{N-1}\) are mutually independent fair
   bits.
4. One seed may be shared by all messages. The proof does not claim that
   different messages have independent output words.
5. A theorem about a repeated random outer code needs either a fixed-code
   prefix check or a valid analysis of the required high moments. Replacing
   a high moment by a power of an expected local spectrum is invalid.

The proof remains valid for a shortened subcode. One can instead analyze a
larger parent code: a minimum-distance bound for the parent also applies to
every subcode.

## Why RM2Sub does not inherit the lemma

The Toeplitz activation condition is \(u_\tau=1\). Every nonzero word has
such a coordinate. Equation (3) then exposes one new random bit per later
output coordinate.

RM2Sub observes a 128-bit epoch \(x_e\) through a 19-bit syndrome \(Bx_e\).
A nonzero epoch in \(\ker B\) does not activate the state. An outer word can
satisfy

\[
 Bx_e=0
 \quad\text{for every epoch }e
\]

while remaining nonzero. RM2Sub acts deterministically on this silent
subcode. Raw prefix dimensions cannot record the weights of those invisible
epochs.

The missing RM2Sub analogue is therefore not another prefix-rank identity.
It is a weight-sensitive syndrome-prefix bound that separates silent input
weight from live-state mixing.

## Concrete certified instance

The finite \(k=2^{20}\) experiment fixes a routed repeated Golay--BA-3 outer
with

\[
 N=2{,}119{,}680,
 \qquad
 D=233{,}164.
\]

Its outward evaluator computes the complete prefix profile and proves

\[
 \log_2\mathbb E_h[Z_{\mathrm{bad}}]
 <-190.5216915553.
\]

The certificate claims the conservative bound

\[
 \Pr_h[d_{\min}\le233{,}164]<2^{-180}.
\]

Thus the sampled code has minimum relative distance above 11%, except with
probability below \(2^{-180}\). The certificate, manifest, evaluator, and
outward receipt are:

- `FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md`;
- `FINITE_K20_BA240_REPEATED_TOEPLITZ_MANIFEST.json`;
- `certify_ba240_repeated_toeplitz_prefix_outward.py`; and
- `ba3_B240_repeated_toeplitz_prefix_outward.json`.

## Scope

The suffix law, first-moment identity, and distance-failure reduction are
proved facts. The concrete 180-bit claim is an outward-certified instance of
those facts.

This argument does not prove that full-length Toeplitz multiplication is
linear time. It also does not transfer the certificate to RM2Sub or another
bounded-state inner. Those constructions need a separate transfer theorem.
