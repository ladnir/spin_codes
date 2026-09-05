# Random outer blocks with the Structured SPIN inner

## First step

The first intermediate ensemble keeps the factored permutation and the
RM2Sub-S19 inner. It replaces only the local structured outer map by an
independent uniform random injection. This change isolates the structured
inner without introducing a new inner constituent family.

At the frozen dimensions, the one-active first moment remains favorable.
The binary64 diagnostic gives a 54.8791-bit margin at relative distance
\(0.11\). The modeled ParityFanout-31x33 spectrum gives 55.8647 bits under
the same inner transfer and tilt grid. Thus the modeled structured outer is
0.9855 bits better in this class.

This document gives the exact one-active reduction and the minimal scalable
construction. The numerical margin is diagnostic because the transfer uses
entrywise envelopes and binary64 optimization.

## Minimal scalable construction

Fix the frozen inner parameters

\[
  t:=128,
  \qquad
  s:=19.
\]

Let \(\mathbb F:=\mathbb F_{2^s}\). Fix the audited binary linear maps

\[
  A:\mathbb F_2^s\to\mathbb F_2^t,
  \qquad
  C:\mathbb F_2^t\to\mathbb F_2^s,
  \qquad
  C\circ A=0.
  \tag{1}
\]

The symbol \(C\) denotes the state-update map in this document. It avoids a
collision with the outer block length \(B\). The frozen materials call the
same map \(B\).

Choose an even outer length \(B\), set \(K_B:=B/2\), and choose a positive
integer \(L\) divisible by \(t\). Define

\[
  N:=LB,
  \qquad
  K:=LK_B=N/2,
  \qquad
  E:=L/t.
  \tag{2}
\]

Setup samples the following independent objects in the listed order.

1. For every \(i\in[L]\), sample a uniform linear injection

   \[
     O_i:\mathbb F_2^{K_B}\to\mathbb F_2^B.
   \]

2. For every \(i\in[L]\), sample a coordinate permutation
   \(\sigma_i\gets S_B\).
3. Transpose the resulting \(L\)-by-\(B\) matrix. For every region
   \(j\in[B]\), sample an independent permutation \(\pi_j\gets S_L\).
4. For every inner epoch \(r\in[N/t]\), sample
   \(\alpha_r\gets\mathbb F^*\).

The local coordinate permutations are distributionally redundant for the
random outer. They remain in the ensemble because the later structured outer
requires them.

Serialize the permuted regions in the order used by the inner recurrence.
Write the result as
\(X_0,\ldots,X_{N/t-1}\in\mathbb F_2^t\), and set \(Q_0:=0\).
For each epoch \(r\), compute

\[
  Y_r:=X_r+A(Q_r),
  \qquad
  Q_{r+1}:=\alpha_r Q_r+C(X_r).
  \tag{3}
\]

The frozen implementation evaluates the equivalent recurrence in reverse
storage order. Reversing the epoch labels does not change Hamming weight.
An implementation-equivalence proof remains a separate obligation.

This generalization changes fewer inner objects than the alternative
\(t=B/2\). Keeping \(t=B/2\) would require a new family of maps \(A\) and
\(C\) whenever \(B\) changes. The present model keeps the audited
\(t=128,s=19\) constituent and varies only the number of outer coordinates
and regions.

## The inner map is bijective

**Lemma 1.** For every fixed multiplier schedule and every fixed pair of
linear maps \(A,C\), recurrence (3) defines a bijection from
\((\mathbb F_2^t)^{N/t}\) to itself.

*Proof.* Fix an output sequence \(Y_0,\ldots,Y_{N/t-1}\). The initial state
\(Q_0\) is known. At epoch \(r\), compute

\[
  X_r=Y_r+A(Q_r),
\]

and then compute \(Q_{r+1}\) from (3). This procedure recovers one input
sequence from every output sequence. Therefore the inner map is bijective.
\(\square\)

Lemma 1 does not use \(C\circ A=0\). That identity supports the efficient state
analysis and implementation, not injectivity.

## Epoch-transfer interface

Fix \(z\in(0,1)\). The state class is \(0\) when \(Q=0\) and \(1\) when
\(Q\ne0\). For \(a\in\{0,1\}\), let \(W_a(z)\) be a nonnegative two-by-two
matrix. Its rows index the entering state class. Its columns index the
leaving state class.

The matrix \(W_a(z)\) is valid if every entry upper-bounds the corresponding
conditional expectation

\[
  \mathbb E\!\left[
    z^{\operatorname{wt}(Y_r)}
    \boldsymbol1\{Q_{r+1}\text{ is in the declared class}\}
    \mathrel{\big|}
    Q_r\text{ is in the entering class},
    \operatorname{wt}(X_r)=a
  \right].
  \tag{4}
\]

The bound must hold after conditioning on every reachable prior history
represented by the entering row. The expectation in (4) is over the fresh
multiplier, the uniform input support, and any declared live-state averaging
that remains after this conditioning. This uniform history condition permits
the epoch matrices to be multiplied.

For \(a=1\), the input support is uniform among the \(t\) epoch coordinates.
For a live entering state, the distribution in (4) must match the declared
uniform-nonzero or punctured-uniform envelope. A zero/live label alone does
not establish this distribution.

The existing RM2Sub receipts provide the data used to construct \(W_0\) and
\(W_1\):

- the exact spectrum of the image of \(A\);
- the exact weight-shell fractions in the kernel of \(C\); and
- the punctured-live likelihood factor \((2^s-1)/(2^s-2)\).

The current evaluator treats the resulting matrices as entrywise upper
bounds. A paper proof must derive (4) from the sampled multiplier law and the
audited fixed maps.

This two-state envelope is suitable for the frozen finite diagnostic but not
for the limit \(L\to\infty\) with fixed \(s\). Its zero-input live entry pays
\((2^s-1)/(2^s-2)\) once per epoch, so a region pays that factor to the power
\(L/t\). RM2SUB_ONE_ACTIVE_CONTINUUM.md replaces it by the exact
\(2^s\)-state epoch kernel and reduces to two states only after taking the
long-region limit. No construction parameter changes.

## One-region transfer

One transposed region contains \(L=tE\) bits. Define the inactive-region
matrix

\[
  R_0(z):=W_0(z)^E.
  \tag{5}
\]

For one active outer block, an active region contains one input bit. The
region permutation makes its position uniform among all \(L\) positions.
Its epoch is therefore uniform in \([E]\). Its position inside that epoch is
uniform in \([t]\). Define

\[
  R_1(z)
  :=
  \frac1E
  \sum_{r=0}^{E-1}
  W_0(z)^r W_1(z)W_0(z)^{E-1-r}.
  \tag{6}
\]

Equation (6) is exact once \(W_0,W_1\) are exact. It is an entrywise upper
transfer when the epoch matrices only satisfy (4).

Let \(e_0=(1,0)^{\mathsf T}\) and
\(\boldsymbol1=(1,1)^{\mathsf T}\). For \(0\le h\le B\), define

\[
  C_{B,h}(z)
  :=
  [u^h]
  e_0^{\mathsf T}
  \big(R_0(z)+uR_1(z)\big)^B
  \boldsymbol1.
  \tag{7}
\]

The coefficient in (7) sums the conditional output-weight transforms over
all \(\binom Bh\) choices of the active regions. Matrix multiplication keeps
the state boundary between consecutive regions.

## Exact random-outer reduction

For one random outer injection, define

\[
  \rho_{B,K_B}:=
  \frac{2^{K_B}-1}{2^B-1}.
  \tag{8}
\]

For every nonzero word \(v\in\mathbb F_2^B\), the expected number of nonzero
messages mapped to \(v\) is \(\rho_{B,K_B}\). This constant is independent
of the weight and support of \(v\).

Let \(Z_{D,1}\) count messages supported on exactly one outer block whose
encoded output has weight at most \(D\). The probability and expectation
below include the random outer maps, coordinate permutations, region
permutations, and multiplier schedule.

**Theorem 2 (one-active first moment).** Suppose \(W_0(z)\) and \(W_1(z)\)
satisfy (4) for every \(z\in(0,1)\). Then

\[
\begin{split}
  \mathbb E[Z_{D,1}]
  \le
  L\rho_{B,K_B}
  \sum_{h=1}^{B}
  \min\left\{
    \binom Bh,
    \inf_{0<z<1}z^{-D}C_{B,h}(z)
  \right\}.
  \tag{9}
\end{split}
\]

A weaker common-tilt form is

\[
\begin{split}
  \mathbb E[Z_{D,1}]
  \le{}&
  L\rho_{B,K_B}
  \inf_{0<z<1} z^{-D}
  e_0^{\mathsf T}
  \left(
    (R_0(z)+R_1(z))^B-R_0(z)^B
  \right)
  \boldsymbol1.
  \tag{10}
\end{split}
\]

*Proof.* Choose the active outer block in \(L\) ways. For a fixed nonzero
outer word of weight \(h\), the local coordinate permutation and region
permutations give the coefficient in (7). The Chernoff inequality

\[
  \boldsymbol1\{w\le D\}\le z^{-D}z^w
\]

gives the second term inside the minimum in (9). The first term uses the
trivial probability bound one. Equation (8) supplies the expected message
multiplicity of every nonzero outer word. Summing over \(h\) proves (9).

For (10), use one common \(z\). Expanding \((R_0+R_1)^B\) sums all binary
active-region patterns. Subtracting \(R_0^B\) removes the zero pattern.
The difference is the nonnegative sum of the remaining matrix words; it does
not subtract two independently derived probability bounds.
\(\square\)

Markov's inequality gives

\[
  \Pr[Z_{D,1}>0]\le\mathbb E[Z_{D,1}].
  \tag{11}
\]

Theorem 2 controls only occupation one. It does not yet imply a
minimum-distance theorem for the complete code.

## All-active anchor

The opposite endpoint does not require an RM2Sub state envelope. Let
\(Z_{D,L}\) count bad messages for which all \(L\) outer blocks are active,
and define the Hamming-ball volume

\[
  V_N(D):=\sum_{w=0}^{D}\binom Nw.
  \tag{12}
\]

**Theorem 3 (all-active first moment).** For every fixed factored permutation
and every fixed RM2Sub multiplier schedule, the expectation over the random
outer maps satisfies

\[
  \mathbb E[Z_{D,L}]
  \le
  \left(
    \frac{2^{K_B}-1}{1-2^{-B}}
  \right)^L
  \frac{V_N(D)}{2^N}.
  \tag{13}
\]

*Proof.* Fix an all-active message. Its random outer rows are independent and
uniform on \(\mathbb F_2^B\setminus\{0\}\). Relax them to independent uniform
binary rows. Conditioning every relaxed row to be nonzero costs
\((1-2^{-B})^{-L}\).

Before conditioning, the complete outer word is uniform on
\(\mathbb F_2^N\). The factored permutation preserves this law. Lemma 1 shows
that the fixed inner maps the uniform law to itself. Therefore the probability
of output weight at most \(D\) is \(V_N(D)/2^N\). There are
\((2^{K_B}-1)^L\) all-active messages. \(\square\)

At rate one half and fixed \(\delta<1/2\), equation (13) has normalized
exponent at most

\[
  H_2(\delta)-\frac12+o(1).
  \tag{14}
\]

Thus every \(\delta\) strictly below the rate-half Gilbert--Varshamov point
closes the all-active class exponentially. `RM2SUB_DENSE_OCCUPATION.md` now
interpolates uniformly between bounded occupation and this endpoint.

## Frozen-dimension diagnostic

The diagnostic uses

\[
  B=256,
  \quad K_B=128,
  \quad L=8192,
  \quad t=128,
  \quad s=19,
  \quad D=230686.
\]

It uses a log-surprisal grid from \(-10\) through \(-6\) with spacing
\(0.01\). The results are:

| outer spectrum | one-active margin | dominant weight |
|---|---:|---:|
| uniform random injection | 54.879111670 bits | 28 |
| modeled ParityFanout-31x33 | 55.864654977 bits | 29 |

The common-tilt random-outer form (10) gives 51.876707694 bits. Optimizing
the tilt separately for every weight yields the stronger first row.

The all-active bound (13) gives 188.379385174 bits at the same distance. It
depends only on the exact random-outer relaxation and inner bijectivity.

The structured comparison differs from the frozen receipt by approximately
\(1.3\times10^{-5}\) bits because this run uses a uniform \(0.01\) grid. The
comparison uses the same inner matrices and grid for both outer spectra.

These values remain diagnostics. The later dense-occupation proof does not
use this grid. It uses an exact four-state certificate below density
\(10^{-4}\) and 100-digit outward interval boxes above that density.

## Admissible asymptotic family

Fix a constant \(c>0\). For every positive integer \(m\), set

\[
  L_m:=tm.
\]

Let \(B_m\) be the least positive even integer satisfying

\[
  B_m\ge c\log_2(L_mB_m),
  \tag{15}
\]

and define \(N_m:=L_mB_m\). Then

\[
  B_m=\Theta(\log N_m),
  \qquad
  \frac{N_{m+1}-N_m}{N_m}=o(1).
  \tag{16}
\]

The family has rate exactly \(1/2\). The fixed inner has \(N_m/t\) epochs.
The random dense local outers cost \(O(N_mB_m)=O(N_m\log N_m)\) bit
operations under direct multiplication.

Define \(M_{B,L}(D)\) as the right side of (9) without the leading factor
\(L\). A sufficient one-active asymptotic target is

\[
  M_{B_m,L_m}(\lfloor\delta N_m\rfloor)
  \le
  N_m^a2^{-\lambda_1B_m}
  \tag{17}
\]

for constants \(a\ge0\) and \(\lambda_1>0\). Equations (11), (15), and (17)
give \(\Pr[Z_{D,1}>0]=o(1)\) whenever

\[
  c>\frac{a+1}{\lambda_1}.
  \tag{18}
\]

The one-active continuum calculation proves a positive exponent. For
\(\delta=0.11002\), it gives the optimized binary64 threshold
\(c>3.587179987\ldots\). An exact rational inequality certifies the convenient
choice \(c=18/5\). See RM2SUB_ONE_ACTIVE_CONTINUUM.md.

## Step-by-step proof ramp

1. **One-active continuum -- complete.** Set \(z=e^{-\theta/L}\), use the
   exact full-state epoch kernel, and take \(E=L/t\to\infty\). The limiting
   marked-region matrices and logarithmic-block constant are proved in
   RM2SUB_ONE_ACTIVE_CONTINUUM.md.
2. **Occupation two -- complete.** For two active outer blocks, use the
   exact zero-, one-, and two-impulse region transfers. The same-epoch
   collision has probability \((t-1)/(L-1)\). The continuum matrix and exact
   \(c=18/5\) certificate are proved in
   RM2SUB_TWO_ACTIVE_CONTINUUM.md.
3. **General fixed occupation -- complete.** The ordered-integral and
   beta-transform formulas in RM2SUB_FIXED_OCCUPATION_CONTINUUM.md give the
   exact common-tilt exponent for every fixed \(q\).
4. **One constant for every fixed occupation -- complete.** The beta-integral
   domination in RM2SUB_UNIFORM_FIXED_OCCUPATION.md proves an exact rational
   \(c=9\) certificate at \(\delta=0.11\), simultaneously for the numerical
   value of every fixed \(q\ge3\). The \(c=18/5\) receipts remain the sharper
   results for \(q=1,2\).
5. **Growing sparse occupation -- complete by a different finite lift.**
   `RM2SUB_DENSE_OCCUPATION.md` conditions an iid Bernoulli candidate law on
   exactly \(q\) marked positions. Its four-state epoch transfer proves a
   uniform \(qB\)-scale exponent for every \(q\to\infty\) with
   \(q/L\le10^{-4}\). It does not require \(q^2B/L=o(1)\).
6. **Bulk occupation -- complete.** The three-state Bernoulli transfer and
   31 outward-certified Collatz boxes cover every
   \(10^{-4}\le q/L\le1\).
7. **Structured outer.** Replace the random local spectrum by a structured
   spectrum through the comparison in
   STRUCTURED_TRANSLATION_FROM_RANDOM_BASELINE.md.

Stages 1--4 change no construction parameter beyond the scalable \(B,L\)
family. The next proof task is a finite-\(L\) transfer bound that remains
uniform through a growing sparse range.

## Status

- **Complete random-outer theorem:**
  `RANDOM_OUTER_RM2SUB_CERTIFICATE.md` proves rate one half and asymptotic
  relative distance \(0.11\) for \(B=9\log_2N+O(1)\).
- **Proved here:** the admissible family definition, inner bijectivity, the
  region-placement identity, Theorem 2 conditional on valid epoch envelopes,
  and the unconditional all-active reduction in Theorem 3. The companion
  continuum notes prove the exact exponent for every fixed occupation
  without using the finite two-state envelope.
- **Imported exact data:** the complete image spectrum of the frozen \(A\)
  map and the exact kernel-shell counts of the frozen \(C\) map.
- **Diagnostic:** the 54.8791-bit random-outer margin and the 0.9855-bit
  structured advantage at the frozen dimensions.
- **Still open:** implementation equivalence for a scalable encoder and a
  transfer from the random outer to a structured outer family with a proved
  uniform spectrum comparison. The obsolete finite two-state diagnostic is
  not needed by the completed asymptotic theorem.
