# Finite 11% certificate with one repeated BCH subcode

## Certified statement

Let \(C_{124}:\mathbb F_2^{124}\to\mathbb F_2^{250}\) be the fixed
codimension-one subcode of the audited \([250,125,\ge 38]\) shortened BCH
code obtained by fixing one basis coordinate to zero. Use the same map
\(C_{124}\) in every outer row.

Set

\[
  L:=8576,\qquad
  N:=250L=2{,}144{,}000,\qquad
  K_0:=124L=1{,}063{,}424,
\]

and

\[
  D:=235{,}840=0.11N.
\]

For every row \(i\in[L]\), setup samples an independent composition \(F_i\)
of 256 ParityFanout-31x33 layers. It then samples an independent uniform
coordinate permutation \(\sigma_i\in S_{250}\). The row map is

\[
  x_i\longmapsto \sigma_i(F_i(C_{124}(x_i))).
\]

The maps \(F_i\) and \(\sigma_i\) are row-local inner wrappers. The outer
constituent \(C_{124}\) is fixed and repeated.

Each ParityFanout layer is an involution because its source and target sets
are disjoint. The fixed RM2Sub recurrence is also bijective. Hence the parent
encoder is injective for every setup realization.

After transposition, setup samples one independent uniform permutation of
the \(L\) positions in each of the 250 regions. It also samples the nonzero
RM2Sub-S19 multiplier for every 128-bit epoch. All setup objects are mutually
independent.

Let \(\mathcal C_{\mathrm{parent}}\) be the resulting random binary linear
code. Then

\[
 \Pr\!\left[
   d_{\min}(\mathcal C_{\mathrm{parent}})\le D
 \right]
 <2^{-47}.
 \tag{1}
\]

The probability in (1) is over the row-local wrappers, local coordinate
permutations, region permutations, and RM2Sub multipliers.
Setup uses no rejection event, so there is no separate setup-failure term.

Fix any injective linear map

\[
  J:\mathbb F_2^{2^{20}}\longrightarrow\mathbb F_2^{K_0}
\]

that appends 14,848 zero coordinates. Define
\(\mathcal C_{20}:=\mathcal C_{\mathrm{parent}}\circ J\). Shortening removes
codewords. Therefore

\[
 \Pr\!\left[
   d_{\min}(\mathcal C_{20})\le235{,}840
 \right]
 <2^{-47}.
 \tag{2}
\]

In particular, with probability greater than \(1-2^{-47}\),

\[
 d_{\min}(\mathcal C_{20})\ge235{,}841
 \quad\text{and}\quad
 \frac{d_{\min}(\mathcal C_{20})}{N}>0.11.
\]

The exact code rate is

\[
 \frac{2^{20}}{2{,}144{,}000}
 =0.4890746268\ldots.
\]

## Spectrum comparison

For a nonzero vector \(v\in\mathbb F_2^{250}\), let \(\mu(v)\) be the
expected multiplicity of \(v\) after one row-local fanout wrapper and its
uniform local coordinate permutation. Let

\[
  \rho:=\frac{2^{124}-1}{2^{250}-1}.
\]

The spectrum verifier constructs the one-layer weight transition over the
exact rationals. It raises that matrix to the 256th power with 384-bit Arb
arithmetic.

For every output shell, the verifier maximizes a valid upper bound over all
source spectra with total mass \(2^{124}-1\) and the authenticated packing
caps. It uses an LP-dual threshold, so floating-point ordering selects no
accepted inequality. The receipt proves

\[
  \mu(v)<2^{1/1000}\rho
  \qquad
  \text{for every }v\ne0.
 \tag{3}
\]

The worst shell is weight one. Its displayed excess is
0.000791813847 bits. The strict rational target in (3), not this display,
is used by the distance verifiers.

Independent wrappers make the occupation-\(Q\) product comparison cost at
most \(2^{Q/1000}\). Equation (3) would not justify this product bound if one
fanout composition were reused across rows.

## Occupation bounds

Let \(Z_Q\) count nonzero parent messages supported on exactly \(Q\) outer
rows whose encoded word has weight at most \(D\).

The outward small-occupation receipt proves

\[
 \sum_{Q=1}^{31}\mathbb E[Z_Q]<2^{-48}.
 \tag{4}
\]

It uses one-ULP outward positive arithmetic, exact binary64 witnesses, the
audited RM2Sub impulse recurrence, and the comparison in (3).

The outward dense receipt partitions \(32\le Q\le8576\) into 21 intervals.
For a fixed witness \((p,z,v)\), the finite logarithmic bound is a linear
function of \(Q\) plus

\[
  -(250-1)\log\binom{8576}{Q}.
\]

The bound is convex in \(Q\). Thus both integer endpoints cover the complete
interval. A 100-digit interval replay proves

\[
 \sum_{Q=32}^{8576}\mathbb E[Z_Q]<2^{-134}.
 \tag{5}
\]

Equations (4) and (5) imply

\[
 \mathbb E\!\left[\sum_{Q=1}^{8576}Z_Q\right]
 <2^{-48}+2^{-134}<2^{-47}.
\]

Markov's inequality proves (1).

## Complexity

The constituent sizes, number of fanout layers, RM2Sub epoch width, and
state size are fixed. The proof-model encoder therefore uses \(O(N)\) bit
operations.

Before circuit fusion, the fanout proxy is

\[
  256(31+33-1)=16{,}128
\]

scalar XORs per row and 138,313,728 scalar XORs for all rows. This is an
operation count, not a performance measurement.

## Receipts

- bch250_124_parityfanout31x33_l256_spectrum_outward.json proves (3).
- bch250_124_parityfanout31x33_l256_q1_31_outward_d11.json proves (4)
  conditional on (3).
- bch250_124_parityfanout31x33_l256_dense_q32_8576_outward_d11.json proves
  (5) conditional on (3).
- FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_MANIFEST.json binds the theorem,
  scripts, and receipts by hash.

## Scope and remaining work

This is a finite theorem for the stated proof-model ensemble. It does not
give an asymptotic theorem with fixed \(B=250\).

The following implementation obligations remain:

1. identify the exact basis coordinate used for \(C_{124}\) in the optimized
   BCH encoder;
2. implement the 256 independent row-local fanout compositions and the local
   coordinate permutations;
3. prove equivalence between the ordinary and transposed implementations and
   the proof-model ordering; and
4. benchmark setup, memory, ordinary encoding, and transposed encoding.

These items do not change the mathematical proof-model statement. They are
necessary before attributing the certificate to a concrete implementation.
