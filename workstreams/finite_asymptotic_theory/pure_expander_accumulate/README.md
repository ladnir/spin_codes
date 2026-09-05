# Pure sparse-mixer--accumulate route

## Status

This folder now contains a complete finite certificate for the one-stage
degree-33 sparse-mixer--accumulator outer with RandomStepConv-M22. The
certificate proves 10.9% relative distance with 41.3621 bits of overall
margin at message dimension \(2^{20}\). It defines a new SPIN proof variant;
it does not alter or certify the frozen Structured SPIN construction.

Initial kernel derivations and parameter diagnostics are recorded in
`PAIR_KERNEL_AND_FIRST_GATE.md`.

The one-stage concentration reduction is recorded in
`ONE_STAGE_VARIANCE_ROUTE.md`. The completed theorem and reproduction commands
are in `../ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`.

The route replaces the fixed BCH base and its permutation--accumulator
stages by a sparse linear map followed by an accumulator. The purpose is to
retain inexpensive mixing while making the first and second spectrum moments
accessible. The construction should use one sampled constituent repeatedly
across all SPIN outer rows. It must not resample a constituent for each row.

## Question

The finite SPIN proof needs a high-probability bound on the realized spectrum
of one repeated outer constituent. An expected spectrum is not sufficient.
For the prior BCH--permutation--accumulator route, the missing variance bound
depends on the unknown joint pair distribution of the BCH code.

The present question is whether a sparse random linear map removes that
obstruction. The full message space has an exact pair-type enumerator, and a
suitable sparse ensemble can have an explicit two-word transition law.

## Candidate architecture

The initial rate-half map has type

\[
 E_0:\mathbb F_2^K\longrightarrow\mathbb F_2^B,
 \qquad B=2K.
\]

The first candidate parameters are

\[
 (K,B)\in\{(128,256),(256,512)\}.
\]

Let \(A_B\) denote the zero-initialized binary prefix accumulator on
\(B\) coordinates. A one-stage constituent would be

\[
 C_1(u):=A_B E_0u.
\]

If one stage is insufficient, later sparse maps have type

\[
 E_i:\mathbb F_2^B\longrightarrow\mathbb F_2^B
 \qquad(i\ge1),
\]

and the tentative depth-\(t\) constituent is

\[
 C_t:=A_BE_{t-1}\cdots A_BE_1A_BE_0.
\]

This display also describes multistage candidates. The completed finite
certificate fixes one stage, degree 33, 16-attempt rank-tested setup, uniform
row and region permutations, and RandomStepConv-M22.

Setup samples the maps \(E_0,\ldots,E_{t-1}\) once. Every SPIN outer row uses
the same realized map \(C_t\). Distinct stages may use independent setup
randomness. No map is resampled as a function of the message or row index.

## First ensemble to analyze

The proof-oriented baseline gives each output coordinate a fixed degree
\(r\). Independently for every output coordinate, setup samples an
\(r\)-subset of the input coordinates and outputs the parity on that subset.
This ensemble has irregular input degrees. Its principal advantage is
independence among output neighborhoods.

For a fixed input word of weight \(h\), define

\[
 \beta_r(h):=\frac{K_r(h)}{\binom nr},
\]

where \(n\) is the input length and \(K_r\) is the binary Krawtchouk
polynomial. For an input pair \((c,d)\), one output pair has exact law

\[
 \Pr[(Ec)_j=a,(Ed)_j=b]
 =\frac14\left(
 1+(-1)^a\beta_r(\operatorname{wt}c)
  +(-1)^b\beta_r(\operatorname{wt}d)
  +(-1)^{a+b}\beta_r(\operatorname{wt}(c+d))
 \right).
\]

The output neighborhoods are independent, so this law tensorizes before the
accumulator. The common accumulator then gives a four-state pair process.
This factorization is the intended source of a tractable variance proof.

Left-regular or biregular ensembles are alternatives, not current
requirements. They may suppress weak columns and improve distance per XOR.
They also introduce dependence among output neighborhoods. The project should
add regularity only after measuring its benefit against the proof cost.

## Required outer event

For a realized constituent \(C_t\), define

\[
 A_w(C_t):=
 \left|\{u\in\mathbb F_2^K\setminus\{0\}:
              \operatorname{wt}(C_t(u))=w\}\right|.
\]

The preferred proof interface is an integer cap vector
\((T_w)_{w=1}^B\) and the event

\[
 \mathcal E_{\mathrm{out}}
 :=\{\dim C_t=K\ \text{and}\ A_w(C_t)\le T_w
      \text{ for every }w\}.
\]

The finite certificate must prove

\[
 \Pr[\neg\mathcal E_{\mathrm{out}}]
 \le \varepsilon_{\mathrm{out}}.
\]

It must also prove, uniformly for every realized constituent in
\(\mathcal E_{\mathrm{out}}\),

\[
 \Pr[d_{\min}<D\mid C_t]
 \le \varepsilon_{\mathrm{inner}}.
\]

The two failure terms must satisfy

\[
 \varepsilon_{\mathrm{out}}+\varepsilon_{\mathrm{inner}}<2^{-40}.
\]

The initial finite target is \(k=2^{20}\) and relative distance \(0.109\).
If the calculation has sufficient slack, it should also test the literal
\(0.11\) target. The exact parent length, shortening rule, and puncturing rule
must be fixed before any result is called a certificate.

## Why the pure route may simplify concentration

For messages \(u,v\in\mathbb F_2^K\), the four coordinate counts of
\((u_i,v_i)\) have an exact multinomial enumerator. The base therefore has no
unknown genus-two spectrum. A sparse-mixer pair kernel can, in principle,
give exact or outward-bounded values for

\[
 \mathbb E[A_w(C_t)]
 \quad\text{and}\quad
 \mathbb E[A_w(C_t)^2].
\]

A useful first concentration target is

\[
 \operatorname{Var}(A_w(C_t))
 \le F\,\mathbb E[A_w(C_t)]
 \qquad(1\le w\le B).
\]

For the existing \([512,256]\), 10.9% cap interface, \(F\le512\) is a
known sufficient three-stage benchmark. This value is a comparison target,
not a proved bound for the sparse-mixer ensemble. A two-stage construction
would require a substantially smaller factor under the same interface.

## Work plan and decision gates

1. **Freeze a baseline ensemble.** Specify neighborhood sampling, degree,
   rank handling, stage independence, accumulator termination, and final
   routing. Obtain approval before treating these choices as the construction.
2. **Derive exact kernels.** Derive the one-word and ordered-pair transition
   laws for one sparse-mixer--accumulator stage. Verify them exhaustively at
   small lengths.
3. **Test concentration.** Evaluate depths one and two across a small degree
   range. The main diagnostic is the worst shell ratio
   \(\operatorname{Var}(A_w)/\mathbb E[A_w]\), not only the expected
   spectrum.
4. **Select by distance per XOR.** Count the XORs in every sparse layer and
   accumulator. Retain only parameters that improve the prospective finite
   certificate enough to justify their online cost.
5. **Build a spectrum event.** Use exact or outward-bounded moments to produce
   simultaneous integer shell caps and an explicit outer failure bound.
6. **Attach the SPIN transfer.** Reuse the existing cap-conditional finite
   transfer when its routing assumptions hold. If the final uniform row
   permutation is removed, define and prove a stronger ordered-support or
   transfer-weighted interface.
7. **Certify and benchmark.** Only after the mathematical interface closes,
   build outward verifiers and compare an optimized implementation with the
   current performance reference.

The one-stage target is complete. Compiled complement-orbit receipts cover
sector zero at levels 1 through 159, sector one at levels 1 through 159, and
sector two at levels 2 through 159. The resulting diagonal relaxation closes
only the first defect shells; a certified obstruction begins at shell 65.

The final proof instead uses a dominant-character likelihood comparison. It
certifies variance factors below 245.146 on shells 42 through 79 and 433
through 470. A cancellation-free deviation form supplies every central-shell
variance bound needed for cap construction. Exact integer caps, bounded
rank-tested setup, and the existing RandomStepConv transfer give 41.3621 bits
of overall margin at 10.9% distance. The manifest is
`../ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json`.

The next design question is performance, not closure of this theorem. The
raw row-by-row one-stage count is 16,895 XORs per constituent. On one exact
sample, common-subexpression synthesis reduced this to 9,477 forward XORs
including the accumulator. The corresponding transposed count is 9,733.
These circuit identities are exact, but the counts are not wall-clock
predictions. A matched Peach benchmark gives 8.703 ms for the sparse outer
and 4.896 ms for two optimized BCH constituents at \(k=2^{20}\). The sparse
outer is 1.777 times slower.

The depth-two candidates have compatible mean spectra, but their covariance
certificate remains open. After exact synthesis, a \((17,7)\) sample costs
9,300 forward XORs including two accumulators, only 177 fewer than the
optimized one-stage sample. A degree-7 square-mixer diagnostic improves the
degree-5 positive covariance relaxation by 85.25-fold but still misses the
required factor by roughly \(10^{87}\). Invertible conditioning makes the
dense uniform reference exact, yet requires a new conditional sparse-defect
lemma. The exact probability space, cost receipts, obstruction, and restart
condition are in TWO_STAGE_SPARSE_EA_ROUTE.md.

## Choices for the next variant

The completed theorem fixes one point in the design space. The following
choices remain open for a practical successor:

- block size \(256\) or \(512\);
- row-independent, left-regular, or biregular sparse maps;
- systematic versus nonsystematic \(E_0\);
- direct rank-failure accounting versus efficient rank-tested resampling;
- degree or degree distribution;
- one stage versus a lower-XOR multistage outer;
- independent maps at different stages versus reuse of one map;
- retention or removal of the final uniform row-coordinate permutation; and
- RandomStepConv versus RM2Sub or another practical inner.

## Related work

The prior permutation--accumulator concentration attempt is documented in
`../EBCH128X4_BA_REPEATED_CONCENTRATION_STATUS.md`. The cap interface and a
complete 10.9% comparator appear in
`../FINITE_K20_PAIRWISE_SYSTEMATIC_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`.
The overall finite experiment history appears in
`../FINITE_K20_ATTEMPT_LOG.md`.
