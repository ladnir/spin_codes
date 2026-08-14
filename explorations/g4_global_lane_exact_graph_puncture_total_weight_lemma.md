# Global-lane puncture bound for the total-weight outer

## Scope

This note justifies the data-block factor used by the
`exact_graph_puncture_total_weight` outer witness.  The construction samples
one puncture lane for all selected band-zero tiles.  Puncture events in
different tiles are therefore correlated.  The argument below does not use
independent puncture lanes.

The result applies only under the conditions stated next.  In particular, it
does not justify multiplying arbitrary one-block puncture adjustments.

## Setup

Let \(\mathcal C\subseteq\{0,1\}^{128}\) be the binary \([128,64]\) EBCH
code.  The outer enumerator contains

\[
 B:=16384
\]

data blocks.  Each block has an independent codeword variable in
\(\mathcal C\).  Equivalently, the unnormalized outer sum ranges over the
Cartesian product \(\mathcal C^B\).

Fix \(q\in(0,1)\), and define

\[
 Z(q):=\sum_{C\in\mathcal C}q^{\operatorname{wt}(C)},
 \qquad
 \mu_q(C):=\frac{q^{\operatorname{wt}(C)}}{Z(q)}.
\]

Each block uses an independent uniform permutation of its 128 codeword
coordinates.  The permutation assigns 42 coordinates to band zero.  For a
permuted codeword \(C\), let \(J(C)\) be the number of its nonzero band-zero
coordinates.  Conditioned on \(\operatorname{wt}(C)=w\),

\[
 \Pr[J(C)=j]
 =\frac{\binom{w}{j}\binom{128-w}{42-j}}{\binom{128}{42}}.
\]

The global-lane sampler chooses a lane \(L^*\), a 128-subset \(H\) of the
256 band-zero tiles, and the selected block set

\[
 S(L^*,H):=\{(t,L^*):t\in H\}.
\]

Thus \(|S(L^*,H)|=P:=128\), and all selected blocks are distinct.  The random
variables \((L^*,H)\) are independent of the codeword choices and coordinate
permutations.  In each selected block, the construction punctures an
independent uniform band-zero coordinate.

## Lemma

**Lemma (global-lane puncture moment).**  Under the setup above, define

\[
 c:=q^{-1}-1,
 \qquad
 \lambda:=\frac{c}{5376}=\frac{c}{128\cdot42}.
\]

For every fixed realization of \((L^*,H)\), the punctured data-block sum is at
most

\[
 \left(
   \sum_{C\in\mathcal C}
   q^{\operatorname{wt}(C)}e^{\lambda J(C)}
 \right)^B,
\]

where the expression averages \(e^{\lambda J(C)}\) over the uniform coordinate
permutation.  Consequently, the same bound holds after averaging over
\((L^*,H)\).

### Proof

Fix \((L^*,H)\).  The set \(S(L^*,H)\) contains \(P\) distinct blocks.
Puncturing a zero coordinate leaves the factor \(q^{\operatorname{wt}(C)}\)
unchanged.  Puncturing a nonzero coordinate multiplies that factor by
\(q^{-1}\).  Averaging the punctured coordinate of one selected block therefore
gives the multiplier

\[
 1+c\frac{J(C)}{42}.
\]

The codeword sum ranges over \(\mathcal C^B\), and the coordinate permutations
are independent across blocks.  Hence the exact data-block sum, averaged over
the coordinate permutations and punctured coordinates, equals

\[
 Z(q)^B
 \left(\mathbb E_{C\gets\mu_q}\left[1+c\frac{J(C)}{42}\right]\right)^P.
\]

Linearity of expectation gives

\[
 \mathbb E_{C\gets\mu_q}\left[1+c\frac{J(C)}{42}\right]
 =1+\frac{c}{42}\mathbb E_{C\gets\mu_q}[J(C)].
\]

The inequality \(1+x\le e^x\), applied with \(x\ge0\), implies

\[
 \left(1+\frac{c}{42}\mathbb E[J]\right)^P
 \le
 \exp\left(\frac{Pc}{42}\mathbb E[J]\right).
\]

Because \(B=16384\), \(P=128\), and \(\lambda=c/5376\),

\[
 B\lambda=\frac{Pc}{42}.
\]

The exponential function is convex.  Jensen's inequality therefore gives

\[
 \exp\left(\lambda\mathbb E[J]\right)
 \le \mathbb E\left[e^{\lambda J}\right].
\]

Raising this inequality to the power \(B\) yields

\[
 \exp\left(\frac{Pc}{42}\mathbb E[J]\right)
 \le
 \left(\mathbb E\left[e^{\lambda J}\right]\right)^B.
\]

Multiplication by \(Z(q)^B\) proves the claim.  The bound is independent of
the identities of the blocks in \(S(L^*,H)\), so averaging over the global
lane and tile set preserves it. \(\square\)

## Spectrum form used by the evaluator

Let \(A_w\) be the number of codewords in \(\mathcal C\) with weight \(w\).
The lemma's one-block factor equals

\[
 \sum_{w=0}^{128} A_w q^w
 \sum_j
 \frac{\binom{w}{j}\binom{128-w}{42-j}}{\binom{128}{42}}
 e^{\lambda j}.
\]

The evaluator raises this factor to \(B=16384\).  It separately multiplies
the normalized exact graph factor

\[
 2^{-24}\sum_{w=0}^{128}G_wq^w,
\]

where \(G_w\) is the committed weight spectrum of the graph code.  The graph
factor does not require puncture-lane independence.  Its use requires the 128
graph coordinates to map bijectively to the 128 holes, as specified by the
global-lane construction.

The evaluator must use outward enclosures for the logarithm, exponential, and
log-sum-exp operations.  This note proves the real-valued inequality that
those enclosures must replay.

## Required proof bindings

A theorem-facing certificate that invokes this lemma must bind all of the
following inputs:

1. the global-lane sampler and its proof that the 128 selected blocks are
   distinct;
2. independence of the outer codeword variables across the \(B\) blocks;
3. independent uniform coordinate permutations for the data blocks;
4. independent uniform punctured-coordinate choices within band zero;
5. the full EBCH weight spectrum \((A_w)_{w=0}^{128}\);
6. the graph spectrum \((G_w)_{w=0}^{128}\); and
7. the graph-to-hole bijection.

The lemma does not assume that the sampled puncture lanes are independent
across tiles.  It also does not authorize the legacy product of marginal
adjustments returned by `puncture_selection_adjustment`.

For the present repository state, the manifest must bind these concrete data
and construction inputs:

| Role | Path | SHA-256 |
|---|---|---|
| global-lane specification | `explorations/g4_global_lane_construction_spec.md` | `89244f76e3bf988753ef8b4d814b141fec9b6599fc027bb5954a5f5a3c62e77b` |
| exact layout certificate | `scripts/certify_global_lane_puncture_layout.py` | `2e0f70cde90f94974f4bd86931905d975eaea8a8c1bc727ab5ba0fa1f810aaf2` |
| full EBCH spectrum | `scripts/ebch128_64_spectrum.csv` | `26164b40cee431cd24995887ba485563f02048d64fdea00dc86b3d9a4897a94e` |
| committed EBCH spectrum source | `scripts/EBCH128_64.wd` | `de6b49d292f6742488ab76d799718e1385feab42b1b375c91508ecd7b13e1cb2` |
| graph spectrum | `scripts/ebch128_graph24_spectrum.csv` | `79a3f8f34280996d46ccd076c515ef1ccdc3b1f80f389ddaf2e7a24286c4f2b3` |
| selected witness source | `out/g4_support_lower30_failed8_retuned_local.json` | `5c0d29c9e879ba21a20068e642549a5ffe4992c8807eb37c37f3dc94fc05ebd9` |

The manifest must also bind this note, the witness hardener, the outward
arithmetic module, and the full certifier after their final edits.  A digest
recorded before those edits is not a valid proof binding.

## Required evaluator metadata

The implementation of `_exact_graph_puncture_total_weight_outer` currently
describes an independent-lane Maclaurin argument.  Its numerical expression
already matches the lemma above, but its theorem-facing report must replace
that description.  The outer report must include at least these fields:

```json
{
  "type": "exact_graph_puncture_total_weight",
  "puncture_lemma": "global_lane_exchangeability_jensen_v1",
  "construction_rule": "global-lane-puncture-v1",
  "lane_events_independent": false,
  "data_blocks": 16384,
  "selected_blocks": 128,
  "band0_coordinates": 42,
  "lambda_denominator": 5376,
  "full_spectrum_sha256": "26164b40cee431cd24995887ba485563f02048d64fdea00dc86b3d9a4897a94e",
  "graph_spectrum_sha256": "79a3f8f34280996d46ccd076c515ef1ccdc3b1f80f389ddaf2e7a24286c4f2b3"
}
```

The report must also include the SHA-256 digest of this note.  The hardener
must check \(B=16384\), \(P=128\), and 42 band-zero coordinates before it
emits the lemma identifier.  The end-to-end certifier must copy the
construction-manifest digest and the witness-source digest into its final
report.  The reproduction preflight must reject a manifest that omits any
binding in the table above.

## Selected witness

The completed g=4 certificate selects exactly one source row with this outer
type:

- reference: `g4_support_lower30_failed8_retuned_local.json:0`;
- repository path: `out/g4_support_lower30_failed8_retuned_local.json`;
- row name: `s0f_anchor_0131`;
- source profile: \((393209,26217,78647,26215,0)\);
- exact pole:
  \(q=3702232328233459/9007199254740992\); and
- source-artifact SHA-256:
  `5c0d29c9e879ba21a20068e642549a5ffe4992c8807eb37c37f3dc94fc05ebd9`.

The source row binds the graph-spectrum SHA-256 digest
`79a3f8f34280996d46ccd076c515ef1ccdc3b1f80f389ddaf2e7a24286c4f2b3`.
The rerun manifest must bind this note and the same spectra before the final
certificate can treat the row as theorem-facing.
