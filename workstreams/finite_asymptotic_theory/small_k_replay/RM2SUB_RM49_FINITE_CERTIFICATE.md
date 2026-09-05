# Fixed-RM RM2Sub finite distance certificate

## Result

There is now a complete outward certificate for one finite structured-route
SPIN instance. The instance has message dimension (k=2^{16}=65{,}536) and
output length (N=2^{17}=131{,}072). Except with probability less than
(2^{-42.5779817562755}) over the setup randomness defined below, its minimum
distance satisfies

\[
  d_{\min}\ge 13{,}108=\lfloor 0.10N\rfloor+1.
\]

In particular, the failure probability is less than (2^{-40}). This is an
outward-rounded finite claim, not a binary64 diagnostic.

This instance is not the frozen Structured SPIN (B=256, t=128, s=19)
implementation. It is a separate fixed-RM proof-model instance with
(B=512), (t=64), and (s=14). The certificate makes no performance or
decoder claim.

## Encoder family and probability space

Let (C_{\mathrm{RM}}) be the fixed binary RM((4,9)) code with parameters
([512,256,32]). Split a message in (mathbb F_2^{65{,}536}) into (L=256)
rows of 256 bits. Apply the same fixed encoder (C_{\mathrm{RM}}) to every
row. The outer encoder is not resampled between rows.

Setup samples the following independent objects.

1. For each outer row, sample a uniform permutation of its 512 output
   coordinates.
2. Transpose the resulting (256)-by-(512) bit array into 512 regions of
   length 256. In each region, sample an independent uniform permutation of
   its 256 coordinates.
3. Divide every region into four RM2Sub epochs of length (t=64). For every
   epoch, sample an independent uniform nonzero scalar in the field used by
   the fixed (s=14) RM2Sub map.

The RM2Sub (A/B) maps are the fixed audited maps bound by the input
manifests. For each setup outcome (omega), all maps are binary linear and
define one binary linear encoder

\[
  \operatorname{Enc}_{\omega}:\mathbb F_2^{65{,}536}
  \longrightarrow \mathbb F_2^{131{,}072}.
\]

The only probability in the theorem is over the row permutations, region
permutations, and nonzero epoch scalars. The outer RM code and RM2Sub (A/B)
maps are fixed.

## Theorem

Let (Omega) be the setup space above, with the displayed choices mutually
independent and uniform. Then

\[
\Pr_{\omega\leftarrow\Omega}
  \left[
    \exists m\ne0:
    \operatorname{wt}(\operatorname{Enc}_{\omega}(m))\le13{,}107
  \right]
\le 1.523176871255522\times10^{-13}
<2^{-42.5779817562755}.
\]

Consequently,

\[
  \Pr_{\omega\leftarrow\Omega}
  [d_{\min}(\operatorname{Enc}_{\omega})\ge13{,}108]
  >1-2^{-40}.
\]

The machine receipt gives the certified enclosing interval

\[
\Pr[\mathrm{failure}]
\in
[1.523176871255521289789405701102375552093461859083080373215117088448001989
 \times10^{-13}\;\mathbin{+/-}\;4.29\times10^{-86}].
\]

## Proof decomposition

For a nonzero message, let (Q) be the number of nonzero 256-bit outer
message rows. Thus (1\le Q\le256). The proof partitions the bad event by
(Q) and covers each occupation exactly once.

| occupation | outward reduction | certified margin for the displayed union |
|:--|:--|--:|
| (Q=1) | exact RM spectrum and exact support average | 42.577981756608 bits |
| (Q=2) | shell-dependent Bernoulli references | 74.611094272379 bits |
| (Q=3,4) | regular/all-one separation | 83.102187705873 bits |
| (5\le Q\le29) | three-band pointwise envelopes | 81.502991282228 bits |
| (Q=30) | refined three-band bridge | 1228.637966758792 bits |
| (31\le Q\le256) | refined three-band bridge with packed witnesses | 1293.444358402287 bits |

The first-moment union is positive throughout. The exact RM spectrum supplies
the number of outer words in each shell. Uniform row and region permutations
turn each routed support into a uniform support conditional on its weight.
The fixed RM2Sub transfer then supplies a positive two-state matrix (F_u(z))
for every live count (u). Each occupation bound applies a finite Chernoff
witness and sums only nonnegative terms.

For (Q\ge30), the refined spectrum groups are

\[
  [1,95],\qquad[96,416],\qquad[417,512],
\]

with Bernoulli reference probabilities

\[
  \frac14,\qquad\frac12,\qquad
  \frac{8677722630069483}{10^{16}}.
\]

The checker evaluates the corresponding exact pointwise density envelopes.
For (31\le Q\le256), one signed byte records the selected exact-tenth
Chernoff witness for each canonical band composition. The 2,856,753 packed
witness bytes replace about 800 MB of discarded JSON discovery output.

## Arithmetic and authentication

The sparse certificates use 256-bit Arb intervals directly. The dense checker
also constructs the RM2Sub transfer with Arb. Its large positive batched
stages use binary64 only with an explicit upper enclosure:

- each binomial coefficient is generated from an Arb recurrence and converted
  to an upper binary64 endpoint;
- every positive convolution and dot product is inflated by the standard
  (gamma_n) rounding bound and an explicit subnormal budget;
- every two-by-two matrix product and sum is rounded toward (+\infty) with
  `nextafter`; and
- matrix normalization uses exact powers of two.

The dense aggregation bounds each occupation by its largest certified term
times the exact number of compositions. It then applies the same maximum-term
rule to the occupation and chunk unions. This deliberately spends at most
about 23 bits, while the dense proof has more than 1,293 bits after that loss.

The final checker verifies the hashes of all six component checkers and
receipts, verifies exact nonoverlapping coverage of (Q=1,\ldots,256), and
sums their upper probability intervals with Arb.

Run the final authentication step from the repository root with

```powershell
python workstreams/finite_asymptotic_theory/small_k_replay/certify_rm2sub_rm49_full_distance.py
```

The authoritative final files are:

- `RM2SUB_RM49_FULL_CERTIFICATE_MANIFEST.json`;
- `certify_rm2sub_rm49_full_distance.py`; and
- `rm2sub_rm49_t64_s14_full_distance_outward.json`.

Their SHA-256 hashes at closure are, respectively,

```text
dc7bd009d4eea2030fdbe8c6a660c96fc0c189611bcb591fb1c7e04165556576
f2d6aa634a366f372d86b23da4ddf0a33caa576eb22d249db20be63fd0c97f01
7efa61f29b91aede4773b382e78690743931bdbe720f20083c253470d4605432
```

## Remaining obligations

This finite theorem closes the stated proof-model instance. It does not close
the following separate tasks.

1. Prove implementation equivalence between this exact model and an optimized
   encoder, then measure its XOR count and runtime.
2. Compare this (B=512,t=64,s=14) instance against the finite BCH and
   random-outer parameter frontiers at other message lengths.
3. Instantiate the arbitrary-length wrapper and quantify its exact rate and
   distance losses.
4. Prove an asymptotic constituent schedule from authenticated spectra or
   spectrum envelopes; one finite point does not imply an asymptotic family.
5. Certify the frozen Structured SPIN (B=256, t=128, s=19) outer interface.
   The theorem here does not transfer to that different constituent or inner
   parameter set.
