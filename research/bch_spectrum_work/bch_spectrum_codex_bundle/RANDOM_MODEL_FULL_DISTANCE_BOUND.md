# A conditional distance bound for the BCH outer and random inner

Updated: 2026-09-04.

The remaining occupation tail is now bounded. Combining every occupation gives
an outward-certified first-moment upper bound below 2^-40 for 2^20 message bits.
The full bound has approximately 40.04769 bits. The tail from occupations four
through 8192 is below 2^-83; it needed only to be at most 2^-45.

The statement remains conditional on three BCH shell inequalities supported by
the completed statistical tests. It concerns the ideal RandomStepConv-M22 model
defined below, not RM2Sub. No variance estimate is used.

## Construction and statement

Fix the same binary linear outer code C of length 256 and dimension 128.
It lies between Q and P, the even extensions of the primitive narrow-sense
BCH codes with designed distances 39 and 37, respectively. Thus C has minimum
distance at least 38. It contains the all-one word and is complement-symmetric.
Write A_w for the number of its words of weight w. The three additional
hypotheses are

\[
A_{38}\le3{,}827{,}351{,}840{,}403,\qquad
A_{40}\le50{,}000{,}000{,}000{,}000,\qquad
A_{42}\le5{,}000{,}000{,}000{,}000{,}000.
\tag{1}
\]

Encode each of 8192 message rows into C using a fixed linear isomorphism.
Independently permute the 256 coordinates of each encoded row. Transpose the
array into 256 regions, each containing 8192 bits, and independently permute
each region. All these permutations are uniform and belong to the setup.

The inner encoder starts with a zero 22-bit state. At position i, a fresh
uniform binary linear map sends the state and input bit to the next state and
one output bit. Equivalently, sample an independent uniform 23-by-23 binary
matrix at each position. The setup fixes these matrices and all permutations
once; every message uses that same setup. Discard the final state. Denote the
resulting linear encoder by E_theta. Its output length is N=2^21.

**Conditional distance statement.** For this fixed C, assuming (1),

\[
\Pr_\theta\left[
  \exists m\in\mathbb F_2^{2^{20}}\setminus\{0\}:
  \operatorname{wt}(E_\theta(m))\le209716
\right]
\le U<2^{-40}.
\tag{2}
\]

The exact rational U is recorded in
`generated/bch256_full_random_inner_conditional.json`. On a setup outside this
bad event, the encoder is injective and its image has minimum distance greater
than 209716. The hypotheses in (1) are inequalities about a fixed code; they
are not a conditioning event saying that sampling tests passed.

For q in {1,...,8192}, let F_q be the expected number of bad nonzero messages
with exactly q nonzero message rows. The expectation is over theta. If Z_theta
counts all bad nonzero messages, then

\[
\Pr[Z_\theta\ge1]\le\mathbb E[Z_\theta]
=\sum_{q=1}^{8192}F_q.
\tag{3}
\]

Different messages may have correlated outputs. Equation (3) requires neither
independence between those outputs nor a bound on their variance.

## Why the moment identities apply to this code

Let L be the even extension of the primitive BCH code with designed distance
59 and dimension 71. Its generator is divisible by the generator of Q, so
L is contained in Q and hence in C. This divisibility was checked directly.

Table 7 of Fujiwara and Kusaka, *The Weight Distributions of the (256,k)
Extended Binary Primitive BCH Codes with k<=71 and k>=187*, IEICE 2021,
DOI 10.1587/transfun.2020EAP1119, supplies the complete spectrum of L.
The local source is `paper/e104-a_9_1321.pdf`, printed page 1326.
Every entry of its k=71 column was checked against the retained anchor.
The exact MacWilliams transform has dual minimum distance 16.

Since C contains L, its dual is contained in the dual of L. Therefore C has
dual minimum distance at least 16. A uniform word of C has uniform projections
on every set of at most 15 coordinates: otherwise a nonzero linear dependence
among those coordinates would give a dual word of weight at most 15.
Consequently, every polynomial of degree at most 15 in its weight has the
same expectation as under W distributed as Bin(256,1/2).

This step uses a published spectrum, not a fresh enumeration of 2^71 words.
The anchor, containment, and moment checks are retained in
`generated/bch256_tail_envelope_checks.json`.

## Squared polynomials bound cumulative shell counts

The earlier envelope allowed excessive cumulative mass near weight 70.
We need cumulative bounds rather than the complete unknown spectrum.
Define T_w=sum_(v=1)^w A_v. For 38<=w<128, pair every word in this prefix
with its complement. This gives T_w<2^127 before using any moment information.

Define the binary Krawtchouk polynomials

\[
K_j(v)=\sum_i(-1)^i\binom vi\binom{256-v}{j-i}.
\]

For W distributed as Bin(256,1/2), their orthogonality identity is

\[
\mathbb E[K_i(W)K_j(W)]
=\mathbf 1_{i=j}\binom{256}{j}.
\tag{4}
\]

One proof expands K_j as the sum of all j-coordinate parity characters of
a uniform binary vector. Only identical characters survive in the expectation
of the product. All 64 identities with 0<=i,j<=7 were also checked exactly.

Choose integer coefficients c_j with indices of one parity and j<=7.
Put g(v)=sum_j c_j K_j(v) and h(v)=g(v)^2. Then h is nonnegative, symmetric
under v->256-v, and has degree at most 14. The moment identity for C gives

\[
\sum_{v=0}^{256}A_v h(v)
=2^{128}\sum_j c_j^2\binom{256}{j}.
\]

The zero and all-one words each contribute h(0). If
\(b_w:=\min\{h(v):v\text{ even},\ 38\le v\le w\}>0\), complement pairing yields

\[
T_w\le
\left\lfloor
\frac{2^{128}\sum_j c_j^2\binom{256}{j}-2h(0)}{2b_w}
\right\rfloor.
\tag{5}
\]

The coefficient choices affect sharpness, not validity. The implementation
tests individual orthogonal squares and same-parity kernel polynomials.
It combines (5) with cumulative sums of the existing deterministic shell caps,
the total mass, and bounds from later prefixes. No statistical shell cap is
used in the occupation tail.

An independent checker evaluated each retained square at every weight, checked
its binomial expectation by direct summation, and verified its prefix minimum.
All 68 retained polynomial-square witnesses passed.

## From prefix bounds to an independent-coordinate reference

For 0<rho<1, let nu_rho be Bin(256,rho), and write V_rho(w) for its cumulative
probability through w. Choose a positive integer Gamma_rho such that

\[
T_w\le\Gamma_\rho V_\rho(w)\quad(0\le w\le256).
\tag{6}
\]

For every nonnegative nonincreasing function f on {0,...,256}, summation by
parts gives

\[
\sum_w A_w^+f(w)
=T_{256}f(256)+\sum_{w=0}^{255}T_w\bigl(f(w)-f(w+1)\bigr)
\le\Gamma_\rho\,\mathbb E_{W\gets\nu_\rho}[f(W)],
\tag{7}
\]

where A_0^+=0 and A_w^+=A_w for w>0. All terms multiplying the prefix bounds
are nonnegative. The checker verified (6) exactly for 21 rational choices of
rho and all 257 prefixes: 5,397 comparisons in total.

To apply (7), we need monotonicity of the averaged inner law. For one fixed
input, let a_i indicate whether the state is nonzero. Independently sample
B_i as a fair bit and J_i as a Bernoulli(2^-22) bit at each position. The law is

\[
u_i=a_i\lor x_i,\qquad y_i=u_iB_i,\qquad
a_{i+1}=u_i(1-J_i).
\tag{8}
\]

Indeed, a uniform linear map sends each fixed nonzero state-input vector to
a uniform 23-bit vector. Its output bit and next state are independent.
Coupling two copies of (8) with the same coins preserves activity and output
order when input ones are added. Thus the output-weight moment at z in (0,1]
decreases when ones are added.

Uniform supports of different row weights can be nested using one random
ordering of the 256 coordinates. Region permutations preserve this inclusion.
Holding other rows fixed and averaging the inner setup therefore gives the
nonincreasing function required by (7). Apply (7) successively to the q rows.
The resulting reference has independent Bernoulli(rho) bits in every row and
a multiplicative factor Gamma_rho^q. Reference rows may be all-zero; excluding
them would change the valid reference distribution.

The coupling establishes each message's marginal law. It does not claim that
one sampled linear encoder is pointwise monotone on all messages.

## Transfer bound for every occupation

Fix s>0, put z=exp(-s), b=(1+z)/2, and epsilon=2^-22. The zero-input and
one-input moment matrices are

\[
Z=\begin{pmatrix}1&0\\\epsilon b&(1-\epsilon)b\end{pmatrix},\qquad
T=\begin{pmatrix}\epsilon b&(1-\epsilon)b\\\epsilon b&(1-\epsilon)b\end{pmatrix}.
\]

The states are inactive and active. Set D_rho=(1-rho)Z+rho T. In a reference
region, q candidate positions form a uniform q-subset of the 8192 positions.
A candidate uses D_rho; every other position uses Z. Let R_q denote the
average ordered product over these subsets. Independence between the 256
reference regions gives moment

\[
M_q(s,\rho)=(1,0)R_q^{256}(1,1)^\mathsf T.
\]

The Chernoff bound and the row-envelope argument yield

\[
F_q\le\binom{8192}{q}\Gamma_\rho^q e^{209716s}M_q(s,\rho).
\tag{9}
\]

Both s and rho may differ between occupations. Each choice gives a valid
upper bound for its particular F_q.

For q=4,...,86, the certificate computes R_q by a positive coefficient
recurrence. The independent check instead sums unnormalized support counts
in Arb and divides by the exact binomial coefficient at the end.

For larger q, choose 0<p<1, with p=1 also allowed when q=8192. Define

\[
P_q(p)=\binom{8192}{q}p^q(1-p)^{8192-q},\qquad
M(p)=(1-p)Z+pD_\rho.
\]

At the endpoint, use P_8192(1)=1.

Expanding the positive matrix product shows

\[
M(p)^{8192}=\sum_j P_j(p)R_j\ \ge\ P_q(p)R_q
\]

entrywise. Nonnegative matrix multiplication preserves that order, so

\[
M_q(s,\rho)\le
\frac{(1,0)M(p)^{2^{21}}(1,1)^\mathsf T}{P_q(p)^{256}}.
\tag{10}
\]

The parameter p is an arbitrary coefficient-bound witness, not necessarily
q/8192. Optimizing it avoids a substantial loss in the earlier diagnostic.

For a fast bound on the numerator, write
\(M(p)=\left(\begin{smallmatrix}a&c\\d&f\end{smallmatrix}\right)\), where

\[
\eta=\rho p,\quad d=\epsilon b,\quad f=(1-\epsilon)b,\quad
a=1-\eta+\eta d,\quad c=\eta f.
\]

Its positive eigenvalue and a positive eigenvector are

\[
\lambda=\frac{a+f+\sqrt{(a-f)^2+4cd}}2,\qquad
v=\left(1,\frac d{\lambda-f}\right)^\mathsf T.
\]

Let kappa=max{1,(lambda-f)/d}. Then the all-one vector is at most kappa v
entrywise. It follows that the numerator of (10) is at most kappa lambda^N.
The independent checker does not use this eigenvalue formula: it computes the
full-length matrix power directly at 256-bit precision and checks all 8,106
large-occupation bounds.

## Numerical certification and full sum

The initial binary64 sweeps locate witnesses only. The final calculation
treats every saved s and p as an exact dyadic rational. Arb encloses
transcendental functions. The small-occupation recurrence rounds every positive
binary64 operation upward, tracks powers of two, and repairs positive underflow.

For each q, the certificate records a rational upper bound e_q on log2(F_q).
It then uses the dyadic upper bound 2^ceil(e_q). Every occupation q=4,...,8192
is present exactly once. The exact sum of these dyadic bounds is below 2^-83;
its displayed margin is approximately 83.9999998 bits. In particular, each
bound is at most 2^-58, which alone would give a tail below 2^-45.

The earlier occupation-one through occupation-three bounds, including the
all-one shell, are unchanged. Combining their exact rational aggregate with
the new tail proves (2), under (1). The total has approximately 40.04769 bits.

| Verified component | Result |
|---|---|
| Published anchor and containment | Every Table 7 k=71 entry matched; dual minimum 16 |
| Polynomial and reference bounds | 68 squares and 5,397 exact prefix comparisons |
| Independent small-occupation arithmetic | 83 full-Arb count recurrences |
| Independent large-occupation arithmetic | 8,106 direct full-length matrix powers |
| Complete occupation coverage | q=1,...,8192; exact total below 2^-40 |

The final receipt is `generated/bch256_full_random_inner_conditional.json`.
Run the non-mutating integrity and aggregate checker with

    python -B code/verify_bch_full_random_inner.py

It passed, including 30 JSON dependencies and exact reaggregation. The linked
receipts retain the independent transfer and polynomial checks. This command
does not rerun the long shell experiments or re-enumerate the published anchor.

## What statistical acceptance does and does not establish

The three predeclared shell tests passed, with full independent replay and
audit checks. Their individual false-accept allocations are 2^-40, 2^-41,
and 2^-41. If the joint cap claim is accepted only when all tests pass, its
false-accept probability is at most 2^-40 under ideal IID sampling: for any
fixed code violating a cap, joint acceptance implies acceptance by that false
cap's test. This does not require independence between the three tests.

There is also a precise two-stage interpretation. Consider the ideal experiment
consisting of those tests followed, on acceptance, by an independent ideal
setup. The probability of accepting and producing a bad setup is at most
2^-40. If all caps hold, the event is bounded by U in (2). If a cap is false,
the event is bounded by false acceptance of the joint claim. This is a bound
on the joint experiment, not the probability of a bad setup conditioned on
the already-observed test outcomes. No repeated test run is proposed here.

The actual test tapes came from Windows CNG. The information-theoretic bounds
assume ideal IID bytes; a computational replacement requires the applicable
distinguishing-advantage qualification. Likewise, the setup in (2) is ideal.
A PRG-based or otherwise modified implementation needs a separate argument.

The occupation tail is no longer an open obligation in this random-inner model.
An unconditional deterministic theorem for the fixed BCH code still requires
deterministic proofs of (1). RM2Sub still requires its own transfer analysis.
Before incorporating the result into the paper, the next recommended step is
to audit this explicit construction and conditional statement against the
intended SPIN interface, then package the argument as a theorem and lemmas.
