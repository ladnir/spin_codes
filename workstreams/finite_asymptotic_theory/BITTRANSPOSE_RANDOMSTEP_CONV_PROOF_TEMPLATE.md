# Bit-transpose and RandomStepConv proof template

## Frozen reference

This template records the proof architecture of the certified power-of-two
EBCH128 instance. The immutable reference is commit

`b2788fdb81e30cb78d68efc2ed06eb8af2b94407`.

The reference manifest is
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json`. At the frozen
commit, its SHA-256 digest is

`a24b1853921d7410f9edcc57993d1d99bfdab69aebf59829240d80e73f9e052e`.

The reference theorem has

\[
 (k,N,D,M)=(2^{20},2^{21},228590,20)
\]

and proves

\[
 \Pr_\omega[d_{\min}(C_\omega)<D]
 <2^{-7.58436084118}.
\]

Future experiments must create new receipts and manifests. They must not
overwrite the frozen reference artifacts.

## Construction interface

Fix a binary linear constituent \(C_{\rm out}\subseteq\mathbb F_2^B\) of
dimension \(K\). Repeat this same constituent in \(L\) outer rows. The
message length and routed length are

\[
 k=KL,\qquad N=BL.
\]

The setup samples the following independent objects once:

1. For each outer row, a uniform permutation in \(S_B\).
2. For each of the \(B\) transposed regions, a uniform permutation in
   \(S_L\).
3. The random maps required by the selected inner encoder.

The encoder reuses the sampled setup for every message. The first two steps
form the region-permuted bit transpose. They do not form a uniform
permutation of all \(N\) coordinates.

For RandomStepConv-M, the state lies in \(\mathbb F_2^M\). At position \(t\),
the setup samples a uniform linear map

\[
 M_t:\mathbb F_2^{M+1}\longrightarrow\mathbb F_2^{M+1}.
\]

On input bit \(x_t\) and state \(s_t\), the encoder computes

\[
 (y_t,s_{t+1}):=M_t(x_t,s_t).
\]

The initial state is zero. The output contains \(y_t\) and omits the terminal
state.

## Required outer information

The proof needs more than the constituent dimension and minimum distance.
It requires one of the following authenticated inputs:

- the exact weight enumerator \((A_0,\ldots,A_B)\); or
- a proved nonnegative majorant for every weight shell.

The exact spectrum should be used for the first few active-row occupations.
A pointwise comparison may be used only when its total multiplicative cost
fits inside the desired failure margin.

If the proof uses parity-pivot deletion, every compared reference row must be
uniform over an even-weight shell mixture. A different parity structure
requires a new reference decomposition.

## Required inner property

Fix a Chernoff parameter \(z\in(0,1)\). The reference proof uses nonnegative
two-state transfers \(T_0(z)\) and \(T_1(z)\) for zero and nonzero input
differences. It requires the following properties:

1. Matrix products give the exact tilted output-weight moment.
2. Replacing an input one by zero cannot decrease that moment.
3. The setup randomness used at different steps has the independence stated
   by the transfer.

For RandomStepConv-M, these matrices are

\[
 T_0=
 \begin{pmatrix}
 1&0\\2^{-M}b&(1-2^{-M})b
 \end{pmatrix},
 \qquad
 T_1=
 \begin{pmatrix}
 2^{-M}b&(1-2^{-M})b\\
 2^{-M}b&(1-2^{-M})b
 \end{pmatrix},
 \qquad
 b=\frac{1+z}{2}.
\]

A replacement inner must prove its own transfer and monotonicity statements.
Numerical similarity is insufficient.

## Occupation decomposition

Let \(Q\) be the number of nonzero outer rows in a message difference. The
template bounds every \(Q\in\{1,\ldots,L\}\) separately.

### One active row

Use the exact constituent spectrum. Average each weight shell over the row
permutation and each selected position over its region permutation. This
range detects late activation by minimum-weight outer words.

### Small occupations

For each region, compute the exact uniform-subset transfer coefficient. Use
the aligned parity-pivot deletion bound when its loss is affordable. The
frozen verifier uses this method for \(2\le Q\le99\).

### Middle occupations

Choose positive coefficient fugacities for the candidate positions and the
parity-pivot loads. Positivity converts the required coefficients into matrix
evaluations. The frozen verifier uses this dispersed-pivot method for
\(100\le Q\le16319\).

### Near-full occupations

Write \(R:=L-Q\). Index each region coefficient by its noncandidate count:

\[
 S_h:=\binom Lh^{-1}[v^h](C+vT_0)^L,
 \qquad C:=\frac{T_0+T_1}{2}.
\]

If a region receives \(m\) deleted parity pivots, its noncandidate count is
\(R+m\). For every \(r>0\), use

\[
 \mathbb E[z^W\mid Q]
 \le
 \frac{Q!}{B^Q}r^{-Q}
 e_0^{\mathsf T}
 \left(\sum_{m=0}^Q\frac{r^m}{m!}S_{R+m}\right)^B
 \mathbf1.
\]

The frozen verifier applies this exact-complement bound for \(0\le R\le64\).

## Certificate condition

Let \(U_Q\) be the outward upper bound for occupation \(Q\), including the
outer multiplicity and the Chernoff factor \(z^{-D}\). A certificate with
failure margin \(\lambda\) must establish

\[
 \sum_{Q=1}^L U_Q<2^{-\lambda}.
\]

The probability is over the sampled setup. The union already covers every
nonzero message. No encoding-time randomness remains after setup.

The outward verifier must treat stored witnesses as exact values. It must
enclose all transcendental operations and round every positive matrix
operation upward.

## Substitution checklist

Changing one component creates the following proof obligations:

- **Outer constituent:** authenticate its spectrum or shell majorant and
  rerun every occupation.
- **Row permutation law:** rederive the shell-to-region distribution.
- **Region permutation law:** rederive the uniform-subset coefficients.
- **Inner encoder:** prove its state transfer, independence, and monotonicity.
- **Distance or memory:** regenerate witnesses and rerun the outward verifier.
- **Occupation split:** prove that the new ranges are exhaustive and disjoint.
- **Compact permutation family:** prove that its distribution supports the
  same bounds; implementation convenience does not imply proof equivalence.

## Known limitation and next use

For the frozen EBCH128 constituent, occupation one is dominant. Its
minimum-weight words can place every nonzero coordinate late in the ordered
regions. Increasing inner memory prevents state extinction after activation,
but it does not remove the preceding zero prefix.

The template should next evaluate one of these changes:

1. an outer constituent with a stronger low-weight spectrum;
2. information-set-aware row routing that guarantees early activation; or
3. a bidirectional inner that removes sensitivity to the first activation
   position.

The frozen proof does not cover RM2Sub, a deterministic plain transpose, or
a compact local permutation family without an additional comparison theorem.
