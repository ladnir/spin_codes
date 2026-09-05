# Finite 11% certificate with 56 row-local fanout layers

## Certified statement

Let

\[
 C_{124}:\mathbb F_2^{124}\longrightarrow\mathbb F_2^{250}
\]

be the fixed codimension-one subcode of the audited
\([250,125,\ge 38]\) shortened BCH code. The encoder uses the same map
\(C_{124}\) in every outer row.

Set

\[
 L:=8576,\qquad N:=250L=2{,}144{,}000,\qquad
 K_0:=124L=1{,}063{,}424,
\]

and set \(D:=235{,}840=0.11N\).

For each row \(i\in[L]\), setup samples an independent composition \(F_i\)
of 56 ParityFanout-31x33 layers. Setup also samples an independent uniform
coordinate permutation \(\sigma_i\in S_{250}\). The row map is

\[
 x_i\longmapsto \sigma_i(F_i(C_{124}(x_i))).
\]

After transposition, setup samples an independent uniform permutation of the
\(L\) positions in each of the 250 regions. It also samples the nonzero
RM2Sub-S19 multiplier for every 128-bit epoch. All setup objects are mutually
independent.

Let \(\mathcal C_{\mathrm{parent}}\) be the resulting random binary linear
code. Then

\[
 \Pr\!\left[d_{\min}(\mathcal C_{\mathrm{parent}})\le D\right]<2^{-40}.
 \tag{1}
\]

The probability in (1) is over the row-local wrappers, local coordinate
permutations, region permutations, and RM2Sub multipliers. Setup has no
rejection event.

Fix any injective linear map

\[
 J:\mathbb F_2^{2^{20}}\longrightarrow\mathbb F_2^{K_0}
\]

that appends 14,848 zero coordinates. Define
\(\mathcal C_{20}:=\mathcal C_{\mathrm{parent}}\circ J\). Shortening removes
codewords. Therefore (1) also holds for \(\mathcal C_{20}\). In particular,

\[
 d_{\min}(\mathcal C_{20})\ge235{,}841>0.11N
\]

except with probability less than \(2^{-40}\). The exact code rate is

\[
 \frac{2^{20}}{2{,}144{,}000}=0.4890746268\ldots.
\]

## Shell-sensitive proof

For a nonzero vector \(v\in\mathbb F_2^{250}\), let \(\mu(v)\) be its
expected multiplicity after one row-local wrapper and its local coordinate
permutation. Define

\[
 \rho:=\frac{2^{124}-1}{2^{250}-1}.
\]

The spectrum verifier constructs the exact rational transition of one
ParityFanout layer. It evaluates the 56th matrix power with 384-bit Arb
arithmetic. For each output shell, it maximizes a valid LP upper bound under
the authenticated source-shell caps and the exact mass \(2^{124}-1\).

Define the tail set

\[
 \mathcal T:=\{1,\ldots,12\}\cup\{238,\ldots,250\}.
\]

Let \(G\) be the event that no wrapped row code contains a nonzero word whose
weight lies in \(\mathcal T\). The outward spectrum receipt proves

\[
 \Pr[G^c]<2^{-45}.
 \tag{2}
\]

The same receipt proves the central pointwise comparison

\[
 \mu(v)<2^{389/1250}\rho
 \quad\text{for every }v\text{ with }13\le\operatorname{wt}(v)\le237.
 \tag{3}
\]

The displayed worst central excess is 0.3111493815 bits at weight 13. The
strict rational exponent \(389/1250=0.3112\), not the display, appears in the
verified comparison.

Let \(Z_Q\) count nonzero parent messages supported on exactly \(Q\) outer
rows whose encoded word has weight at most \(D\). The outward receipts prove

\[
 \mathbb E[Z_1]<2^{-53},
 \tag{4}
\]

and

\[
 \sum_{Q=2}^{31}\mathbb E[Z_Q]<2^{-52}.
 \tag{5}
\]

These two bounds use the complete shell-sensitive spectrum. They do not
condition on \(G\). For \(Q=1\), a matrix-coefficient recurrence averages the
ordered RM2Sub region product at every exact row weight. For
\(2\le Q\le31\), a Bernoulli majorant is selected separately for each finite
RM2Sub witness.

On \(G\), every active wrapped row has central weight. Equation (3), the
finite convex interval cover, and the outward RM2Sub replay prove

\[
 \mathbb E\!\left[
   \mathbf 1_G\sum_{Q=32}^{8576}Z_Q
 \right]<2^{-41}.
 \tag{6}
\]

The four strict bounds imply

\[
 \begin{aligned}
 \Pr\!\left[d_{\min}(\mathcal C_{\mathrm{parent}})\le D\right]
 &<2^{-45}+2^{-53}+2^{-52}+2^{-41}\\
 &<2^{-40}.
 \end{aligned}
\]

This proves (1). The displayed receipt endpoints give a combined margin of
about 41.51 bits, but the theorem claims only the integer 40-bit bound.

## Complexity

The constituent sizes, fanout count, RM2Sub epoch width, and state size are
fixed. The proof-model encoder therefore uses \(O(N)\) bit operations.

Before circuit fusion, the fanout proxy is

\[
 56(31+33-1)=3528
\]

scalar XORs per row. Across all rows, the proxy is 30,256,128 scalar XORs.
The earlier 256-layer certificate required 138,313,728 scalar XORs. The new
proof reduces this proxy by 78.125%, or a factor of \(256/56\approx4.57\).
These figures are operation counts, not performance measurements.

## Receipts

- `bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json` proves
  (2) and (3).
- `bch250_124_parityfanout31x33_l56_q1_outward_d11.json` proves (4).
- `bch250_124_parityfanout31x33_l56_q2_31_outward_d11.json` proves (5).
- `bch250_124_parityfanout31x33_l56_cutoff12_dense_q32_8576_outward_d11.json`
  proves (6) conditional on (3).
- `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json` binds the theorem,
  verifiers, witnesses, and receipts by hash.

## Scope and remaining work

This theorem applies to the stated finite proof-model ensemble. It does not
establish an asymptotic theorem with fixed outer length 250.

The implementation obligations are:

1. identify the exact basis coordinate used for \(C_{124}\) in the optimized
   BCH encoder;
2. implement the 56 independent row-local wrapper compositions and local
   coordinate permutations;
3. prove equivalence between the ordinary and transposed implementations and
   the proof-model ordering; and
4. benchmark setup, memory, ordinary encoding, and transposed encoding.

The proof optimization does not establish that 56 layers are intrinsically
necessary. It establishes that 52 layers fail the present tail-plus-central
pointwise argument. A stronger multi-band proof could reduce the count.
