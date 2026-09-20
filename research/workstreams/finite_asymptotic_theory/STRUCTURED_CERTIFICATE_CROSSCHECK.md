# Cross-check of the structured linear-time certificate

## Audited result

The companion linear-time workstream states a complete theorem for a
Structured SPIN family with

\[
  R=\frac12,
  \qquad
  \delta=0.101,
  \qquad
  B=\Theta((\ln N)^2).
\]

Its ordinary and transposed encoding work is \(O(N)\). The outer is a
selected Golay--BA-3 family. The route uses independent local and region
permutations. The inner is the fixed RM2Sub constituent with
\((t,s)=(128,19)\).

The theorem appears in
`../linear_time_audit/ASYMPTOTIC_DISTANCE_CERTIFICATE.md`. Its manifest is
`../linear_time_audit/CERTIFICATE_MANIFEST.json`.

Every SHA-256 entry in that manifest matched during this cross-check:

- three proof receipts;
- five verifier scripts;
- the frozen main-code document and source manifest; and
- the frozen RM2Sub kernel and image-spectrum receipts.

This check authenticates the files named by the manifest. It does not replace
the interval verifiers' mathematical arguments.

## Relation to the random-outer certificate

The two theorems use the same fixed inner but different outer costs. At the
half-Bernoulli reference point, the selected Golay--BA envelope costs

\[
  d_{1/2}<0.36
\]

natural units per active outer block. The random rate-half outer costs

\[
  d_{1/2}=\frac12\ln2=0.34657359\ldots.
\]

The difference is small on an ordinary scale but material near the
rate-half Gilbert--Varshamov point. At \(\delta=0.11\), the random-outer
all-active margin is positive but very small. Substituting the \(0.36\) bound
removes that margin.

Therefore the structured \(0.101\) theorem is a rigorous fallback, not a
transfer of the random-outer \(0.11\) theorem.

## Block-length interpretation

The companion schedule sets \(B_m=24m\) and

\[
  L_m=128\left\lceil\frac{e^{\sqrt{B_m}}}{128}\right\rceil.
\]

Hence

\[
  B_m=\Theta((\ln N_m)^2),
  \qquad
  B_m=o(N_m).
\]

This schedule satisfies the requested sublinear-block fallback. It does not
give a logarithmic structured outer block.

## Remaining \(0.11\) obligation

A structured \(0.11\) theorem needs one of the following improvements:

1. a structured outer family with half-density spectrum excess \(o(B)\)
   relative to the random rate-half outer;
2. a transfer-weighted comparison sharper than the pointwise spectrum
   envelope; or
3. a sharper certified spectrum for the selected Golay--BA family that
   restores the high-occupation margin.

Until one of these statements is proved, the theorem-safe structured result
is \(0.101\), while \(0.11\) is certified only for the random outer.
