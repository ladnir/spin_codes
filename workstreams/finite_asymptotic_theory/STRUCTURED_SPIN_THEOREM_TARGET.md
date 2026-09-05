# Conditional theorem target for Structured SPIN

## What this theorem is meant to establish

Structured SPIN replaces the uniform global interleaver and random dense inner
with structured components.  The proof must therefore preserve the types that
the factored interleaver preserves.  A Hamming-weight-only transfer probability
is valid only after proving the corresponding transitivity statement.

This document gives two targets:

- a finite theorem that can certify the frozen instance; and
- a conditional asymptotic theorem for a scalable Structured SPIN family.

The frozen ParityFanout target remains open. The generic asymptotic target
is now discharged for the one-sampled Golay--BA-3/RM2Sub-S19 family in
`SINGLE_SAMPLED_BA_RM2SUB_D11.md`. Other structured families must still
verify the hypotheses below.

## Construction and probability space

At an admissible length \(N\), let \(B\mid N\) and write \(L=N/B\).  The
message space is \(\mathbb F_2^K\), where \(K=L K_B\).  The setup experiment
must specify the following objects in chronological order.

1. For each outer block \(i\in[L]\), sample the declared local outer
   randomness \(R_i^{\mathrm{out}}\), and construct
   \[
     O_i:\mathbb F_2^{K_B}\to\mathbb F_2^B.
   \]
2. Sample the coordinate, transpose, region, and lane randomness
   \(R^{\mathrm{perm}}\), and construct the factored permutation
   \(\Pi_N\).
3. Sample the declared inner schedule \(R^{\mathrm{in}}\), and construct the
   length-preserving linear map
   \[
     I_N:\mathbb F_2^N\to\mathbb F_2^N.
   \]
4. Output
   \[
     E_N:=I_N\circ\Pi_N\circ(O_1\oplus\cdots\oplus O_L).
   \]

The theorem must state which random variables are independent.  A
pseudorandom implementation can instantiate a sampled setup after correctness
is proved, but a theorem cannot infer independent uniform objects merely from
independent-looking seeds.

For the frozen implementation, the manifest records

\[
  B=256,\quad K_B=128,\quad t=128,\quad s=19,
  \quad L=8192,
\]

and \((K,N)=(2^{20},2^{21})\) in 128-bit blocks.  Each block uses the fixed
extended-BCH-based outer followed by an independently sampled disjoint
ParityFanout-31x33 map.  The interleaver uses per-block coordinate routing, a
transpose, and region permutations.  The inner uses the fixed RM2Sub-S19
constituents.  A paper theorem must replace this implementation description
with exact mathematical sampling algorithms.

## Type interface

For a nonzero message \(x=(x_1,\ldots,x_L)\), define its occupation

\[
  q(x):=|\{i:x_i\ne0\}|.
\]

Let \(\mathcal T_{N,q}\) be a finite set of types rich enough to determine the
orbit of the outer word under the sampled factored interleaver.  A type may
contain local output weights, region occupancies, lane occupancies, and any
boundary data used by the inner transfer proof.  It must not omit an invariant
that changes the transfer probability.

For a realized outer setup \(O=(O_1,\ldots,O_L)\), define

\[
  A_{O,q,\tau}
  :=|\{x\ne0:q(x)=q,\ \operatorname{type}_O(x)=\tau\}|.
\]

Define the conditional transfer average

\[
  P_{O,q,\tau}(D)
  :=\frac{1}{A_{O,q,\tau}}
  \sum_{\substack{x:q(x)=q\\\operatorname{type}_O(x)=\tau}}
  \Pr\!\left[
    \operatorname{wt}(I_N(\Pi_N(O(x))))<D
    \mid O
  \right]
\]

when \(A_{O,q,\tau}>0\).  The probability includes the remaining interleaver
and inner randomness.  Set the product to zero when the count is zero.

The exact first moment is

\[
  \mu_N(D)
  =\sum_{q=1}^{L}\sum_{\tau\in\mathcal T_{N,q}}
   \mathbb E_O[A_{O,q,\tau}P_{O,q,\tau}(D)].
\]

This identity is unconditional.  All theorem work lies in bounding its terms.

## Finite theorem target

Suppose a verifier establishes outward bounds \(U_{q,\tau}\) satisfying

\[
  \mathbb E_O[A_{O,q,\tau}P_{O,q,\tau}(D)]
  \le U_{q,\tau}
\]

for every feasible \((q,\tau)\).  Suppose it also verifies that the type cells
cover every nonzero message exactly once, or with a documented safe
overcount.  Define

\[
  U:=\sum_{q=1}^{L}\sum_{\tau\in\mathcal T_{N,q}}U_{q,\tau}.
\]

**Finite Structured SPIN theorem target.** The setup experiment satisfies

\[
  \Pr\!\left[E_N\text{ is not injective or }d_{\min}(E_N)<D\right]
  \le U.
\]

If \(U<1\), at least one deterministic setup is a binary
\([N,K,d_{\min}\ge D]\) code.  If \(U\le2^{-\lambda}\), the setup-failure
probability is at most \(2^{-\lambda}\).

The theorem is deliberately agnostic about how a cell is bounded.  Exact
enumeration, rational convexity, interval arithmetic, or analytic inequalities
are all valid if the verifier checks the complete domain.

## Weak sufficient hypotheses for an asymptotic theorem

Let

\[
  B_N=\lceil c_S\log_2N\rceil_{\mathcal B},
  \qquad L_N=N/B_N.
\]

Fix \(\delta_S>0\) and
\(D_N=\lfloor\delta_SN\rfloor+1\).  The following two weighted hypotheses are
sufficient.  They do not require a weight-only spectrum or uniform Hamming-
slice transfer law.

### S: sparse-profile product bound

There are constants \(a\ge0\), \(\lambda_S>0\), and
\(\alpha\in(0,1)\) such that, for every sufficiently large admissible \(N\)
and every \(1\le q\le\alpha L_N\),

\[
  \sum_{\tau\in\mathcal T_{N,q}}
  \mathbb E_O[A_{O,q,\tau}P_{O,q,\tau}(D_N)]
  \le
  \binom{L_N}{q}N^a2^{-\lambda_SB_Nq}.
  \tag{S}
\]

This is the weakest sparse interface needed by the first-moment proof: it
bounds the outer spectrum weighted by the actual structured transfer law.
A modular proof may discharge it using separate bounds

\[
  \mathbb E[A_{O,q,\tau}]\le\overline A_{q,\tau},
  \qquad
  P_{O,q,\tau}(D_N)\le\overline P_{q,\tau},
\]

provided the transfer bound is pointwise in the outer realization, followed by

\[
  \sum_\tau\overline A_{q,\tau}\overline P_{q,\tau}
  \le\binom{L_N}{q}N^a2^{-\lambda_SB_Nq}.
\]

Separate average bounds are insufficient if outer and transfer quantities are
dependent.

### B: many-block exponent gap

There is \(\tau_S>0\) such that, for every sufficiently large admissible
\(N\),

\[
  \sum_{q=\lceil\alpha L_N\rceil}^{L_N}
  \sum_{\tau\in\mathcal T_{N,q}}
  \mathbb E_O[A_{O,q,\tau}P_{O,q,\tau}(D_N)]
  \le N^{O(1)}2^{-\tau_SN}.
  \tag{B}
\]

An outer type exponent \(\Phi_S\) and inner type exponent \(E_S\) imply (B)
when their difference is uniformly negative and the number of types is
subexponential in \(N\).

### Conditional theorem

**Theorem target (asymptotic Structured SPIN).** Assume the setup experiment,
admissibility rules, (S), and (B).  If

\[
  c_S>\frac{a+1}{\lambda_S},
\]

then the admissible Structured SPIN family has rate
\(K_N/N\to R_S:=\lim K_{B_N}/B_N\) and

\[
  \Pr[d_{\min}(E_N)<D_N]
  \le N^{-\Omega(1)}+2^{-\Omega(N)}=o(1).
\]

*Proof target.* Sum (S) over \(q\).  With
\(z_N:=2^{-\lambda_SB_N}\),

\[
  \sum_{q=1}^{\alpha L_N}
  \binom{L_N}{q}N^az_N^q
  \le N^a\big((1+z_N)^{L_N}-1\big).
\]

The condition on \(c_S\) gives
\(N^aL_Nz_N=N^{-\Omega(1)}\).  Hypothesis (B) bounds the remaining
occupations.  The exact type first-moment theorem and Markov's inequality
finish the proof. \(\square\)

The one-block term has only \(\Theta(\log N)\) margin under (S).  Thus the
global theorem should not claim an \(\Omega(N)\) failure margin.

### Polynomial sparse alternative

A BA outer with a fixed constituent and a fixed number of accumulators can
have polynomial, rather than exponential, sparse suppression.  Replace (S)
by the following condition.  There are constants \(C>0\), \(p>0\), and
\(\alpha\in(0,1)\) such that

\[
  \sum_{\tau\in\mathcal T_{N,q}}
  \mathbb E_O[A_{O,q,\tau}P_{O,q,\tau}(D_N)]
  \le
  \binom{L_N}{q}(CB_N^{-p})^q
  \tag{S-poly}
\]

for \(1\le q\le\alpha L_N\).  Let

\[
  B_N=N^{\beta+o(1)},
  \qquad L_N=N/B_N.
\]

If \(\beta>1/(p+1)\), then

\[
  L_NCB_N^{-p}=N^{1-\beta(p+1)+o(1)}=o(1).
\]

Summing (S-poly) and applying the unchanged bulk hypothesis (B) gives

\[
  \Pr[d_{\min}(E_N)<D_N]=o(1).
\]

This alternative preserves a sublinear outer block.  It does not imply the
logarithmic-block conclusion of the exponential hypothesis (S).

## Stronger but easier-to-state specialization

If the structured permutation and inner satisfy a weight-only bound

\[
  \Pr[\operatorname{wt}(I_N(\Pi_N(u)))<D_N]
  \le N^a\rho^{\operatorname{wt}(u)}
\]

for every nonzero outer word \(u\), then (S) follows from a local weighted
outer-spectrum bound

\[
  G_B(\rho):=
  \mathbb E\!\left[
    \sum_{c\in C_B\setminus\{0\}}
    \rho^{\operatorname{wt}(c)}
  \right]
  \le2^{-\lambda_SB}.
\]

This specialization is convenient but stronger than necessary.  It should not
be imposed if the structured orbit proof naturally yields profile-dependent
transfer bounds.

## Random-reference comparison specialization

For the coordinate-permutation, transpose, and region-permutation law, define
the expected local weight spectrum \(\overline A_{B,h}\).  Let

\[
  A^{\mathrm{rnd}}_{B,h}
  :=(2^{K_B}-1)\frac{\binom Bh}{2^B-1}
\]

be the expected spectrum of a uniform random injective outer map, and set

\[
  g_B:=
  \max_{h:\,\overline A_{B,h}>0}
  \log_2\frac{\overline A_{B,h}}{A^{\mathrm{rnd}}_{B,h}}.
  \tag{C}
\]

`STRUCTURED_TRANSLATION_FROM_RANDOM_BASELINE.md` proves that every
nonnegative region-profile transfer for the structured outer is at most
\(2^{g_Bq}\) times the corresponding random-outer transfer at occupation
\(q\).  Therefore, if the random outer paired with the same structured inner
has sparse exponent \(\lambda_{\mathrm{rnd}}\), the pointwise comparison (C)
gives

\[
  \lambda_S
  <
  \lambda_{\mathrm{rnd}}
  -\limsup_{B\to\infty}\frac{g_B}{B}.
\]

In particular, \(g_B=o(B)\) preserves the limiting random-outer block
constant.  This condition is sufficient, not necessary.  A direct
transfer-weighted comparison may be sharper when a few weights have large
spectral excess but small bad-output probability.

This specialization does not replace the RM2Sub proof with the
random-convolution proof.  The reference and candidate outers must be paired
with the same inner kernel before their spectra can be compared.

The first reference ensemble is defined in
`RANDOM_OUTER_RM2SUB_RAMP.md`.  It keeps the frozen \(t=128,s=19\) inner,
lets the even outer length grow, and requires \(128\mid L\).  Its exact
conditional one-active theorem is the first subcase of (S).
`RM2SUB_ONE_ACTIVE_CONTINUUM.md` now proves the positive one-active exponent
uniformly for \(B=O(\log L)\). At \(\delta=0.11002\), it gives the optimized
binary64 threshold \(c>3.587179987\ldots\) and an exact rational certificate
at \(c=18/5\). `RM2SUB_TWO_ACTIVE_CONTINUUM.md` proves occupation two. Its
optimized threshold is \(c>3.587809926\ldots\), and the same exact
\(c=18/5\) schedule closes it.
`RM2SUB_FIXED_OCCUPATION_CONTINUUM.md` now proves the common-tilt theorem for
every fixed occupation. Its diagnostic threshold grows with occupation and
first exceeds \(18/5\) at sampled \(q=16\) for \(\delta=0.11\).
`RM2SUB_UNIFORM_FIXED_OCCUPATION.md` replaces that diagnostic extrapolation
with an exact beta-integral bound. It certifies \(c=9\) for every fixed
\(q\ge3\) at \(\delta=0.11\). The quantifier still fixes \(q\) before taking
\(N\to\infty\). The immediate obligation is a finite-\(L\) bound for a growing
sparse range, followed by the bulk bridge.

## Exact obligations for the frozen constituents

The frozen implementation and receipts identify the following proof tasks.

### Outer constituent

1. Give the exact generator of the fixed 256-to-128 extended-BCH-based map.
2. Prove the spectrum or a sufficient weighted envelope after the sampled
   disjoint ParityFanout-31x33 transform.
3. Prove that the setup samples one transform independently for each of the
   8192 outer blocks with the distribution used by the spectrum theorem.
4. Prove every claimed minimum-weight cutoff.  The current modeled receipt
   reports source minimum distance 38 and transformed minimum distance 37,
   but it is not the required spectrum proof.

### Factored interleaver

1. Define the per-block coordinate permutation, transpose, and 256 region
   permutations as mathematical maps.
2. Prove their joint distribution and independence from outer and inner setup.
3. Identify the complete orbit type.  Prove transitivity within each type if
   a representative replaces the average in the transfer definition.
4. Prove that implementation tiling and power-of-two routing realize the same
   permutation as the mathematical construction.

### RM2Sub-S19 inner

1. Authenticate the fixed \(A\)- and \(B\)-constituent matrices used by the
   transfer proof.
2. Define the distribution of each recursive state multiplier or schedule.
3. Prove the transfer inequality for every feasible type, including the first
   step, final step, zero-state transitions, and boundary types.
4. Prove that the transfer theorem uses the same \(t=128,s=19\) recurrence as
   the frozen source.

### Numerical certificate

1. Replace nearest-binary64 evaluation with independently checkable outward
   arithmetic.
2. Cover all occupation classes and all type cells.
3. Bind every spectrum table, transfer table, partition, and witness by hash.
4. Recompute the final sum with outward log-sum-exp or exact rational
   arithmetic.

## Status of the current frozen evidence

At \(N=2^{21}\) and \(D=230687\), the current diagnostic evaluates the event
of output weight at most 230686.  It partitions occupations into

\[
  q=1,\qquad 2\le q\le100,\qquad 101\le q\le8192.
\]

The reported nearest-binary64 margins are approximately 55.8646, 108.8244,
and 121.0192 bits.  The one-active class is limiting.  Lower-distance
diagnostics report approximately 66.2701 bits at relative distance 0.09 and
61.0696 bits at relative distance 0.10.

These numbers are numerical diagnostics.  They rely on a modeled outer
spectrum and nearest-binary64 arithmetic.  They do not establish the finite
theorem target.

## Certified scalable reference points

`RANDOM_OUTER_RM2SUB_CERTIFICATE.md` now proves the structured-inner
reference theorem at rate one half and distance \(0.11\). It keeps the frozen
\((t,s)=(128,19)\) inner and uses \(B=9\log_2N+O(1)\). Thus the inner
many-block exponent in the random-outer comparison is no longer an open
hypothesis.

`SINGLE_SAMPLED_BA_RM2SUB_D11.md` now proves a scalable structured theorem.
It samples one Golay--BA-3 outer, reuses that code at every outer-block
position, and keeps the fixed \((t,s)=(128,19)\) inner. A weight-dependent
concave majorant avoids charging the old \(0.36\) shell penalty at every
weight. Outward certificates cover every positive occupation, and exact
likelihood transfer covers growing sparse occupation and \(Q=1,2\). The
weight-coupled transfer covers every fixed \(Q\ge3\). The result has rate one
half, relative distance \(0.11\),
\(B=(39/4)\log_2N+O(1)\), and linear ordinary and transposed work. The
fixed-occupation improvement keeps the BA weight fugacity inside the RM2Sub
transfer.

The companion linear-time theorem at distance \(0.101\) remains a valid
more conservative reference. Its \(B=\Theta((\ln N)^2)\) schedule is no
longer needed for the one-sampled Golay--BA-3/RM2Sub-S19 family.

## Additional obligations for a scalable family

The actual frozen constituents have fixed \(B,t,s\).  An asymptotic theorem
also requires:

1. a family of structured outer maps for
   \(B_N=\Theta(\log N)\), with uniform constants in (S);
2. a family of factored permutations defined at every admissible \(N\);
3. the fixed RM2Sub family with \((t_N,s_N)=(128,19)\), or another explicit
   inner family;
4. uniform sparse-profile bounds, especially at \(q=1\);
5. a structured-outer comparison to the proved random-outer many-block gap;
   and
6. a bounded-gap admissible-length rule for the arbitrary-length wrapper.

A sufficient outer-family target for item 1 is \(g_B=o(B)\) in (C), provided
items 3--5 establish the random-outer reference bound for the same structured
inner.  If only \(g_B=\epsilon B+o(B)\) is available, the sparse exponent
loses \(\epsilon\).  If the one-active mass is merely polynomial in \(B\), the
logarithmic schedule must be replaced by the polynomial alternative above.

Proving these items defines an asymptotic Structured SPIN family.  It does not
retroactively change the frozen Structured SPIN (B=256, t=128, s=19) instance.
