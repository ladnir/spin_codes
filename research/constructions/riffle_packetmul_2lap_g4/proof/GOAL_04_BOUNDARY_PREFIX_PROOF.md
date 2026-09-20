# Goal 04: boundary-prefix certificate

## Result

Every nonzero four-node autonomous window has binary weight at least 36. This
bound moves the sufficient zero-prefix cutoff from 22,653 nodes to 20,976
nodes. It closes 1,677 additional first-node strata for every nonzero first-lap
terminal state.

In particular, Goal 04 closes the target strata

\[
r\in\{22{,}650,22{,}651,22{,}652\}
\]

without a probabilistic assumption about PacketMul values. The full
construction remains open.

## Conditioned experiment

Fix one authenticated support-33 outer word \(x\). Setup samples the packet
permutation \(\Pi\) uniformly. It independently samples each active multiplier

\[
A_i\gets\mathbb F_{16}^{\times}.
\]

Define \(y:=\Pi(D_Ax)\). Let \(E_r\) be the event that the first occupied node
of \(y\) is \(r\). For \(M=524{,}352\),

\[
\Pr[E_r]
=
\frac{
\binom{M-16r}{33}-\binom{M-16(r+1)}{33}
}{
\binom{M}{33}
}.
\]

Conditioned only on \(E_r\), the active multipliers remain independent and
uniform. The event \(E_r\) depends only on the placement of nonzero packets.

Let \(L(y)\) be the first-lap terminal state. Conditioning further on
\(L(y)\ne0\) generally couples \(\Pi\) and the active multipliers. The resulting
law is uniform over the feasible pairs \((\Pi,A)\) that satisfy both events. It
is not generally a product law.

The proof below does not use independence under this conditional law. It holds
for every feasible realization with \(L(y)\ne0\).

## Four-node lemma

During a zero-input prefix of the second lap, let \(o_j\in\mathbb F_2^{64}\)
denote the output at node \(j\). Consecutive outputs satisfy

\[
o_{j+1}=U(o_j),
\qquad
U(o):=\operatorname{Acc}(P(o)).
\]

The map \(U\) is invertible. Hence every \(o_j\) is nonzero when the prefix is
entered from \(L(y)\ne0\).

**Lemma.** Every consecutive four-node autonomous window has weight at least
36.

**Proof.** Suppose a four-node window has weight at most 35. At least one of
its four outputs then has weight at most eight. The primary search enumerates
every nonzero 64-bit output of weight at most eight. It checks all four window
alignments that contain each enumerated output.

The search enumerates 5,130,659,560 outputs. It finds no aligned window of
weight below 42 among the covered alignments. This contradicts the supposed
window. Therefore, every nonzero four-node window has weight at least 36.
\(\square\)

The smallest weight observed among windows containing an enumerated anchor is
42. The certificate does not claim that 42 is the global minimum. A window of
weight between 34 and 41 might have all four outputs of weight at least nine,
which this enumeration does not anchor.

## Independent verification

The primary implementation constructs the systematic map through a
carryless-polynomial inverse. It updates neighboring outputs incrementally as
the anchor combination changes.

The independent implementation uses binary Gaussian elimination. It computes
accumulation bit-by-bit and reevaluates every neighboring output through byte
tables. The independent implementation reproduces all primary counts:

| Quantity | Primary | Independent |
|---|---:|---:|
| Enumerated anchors | 5,130,659,560 | 5,130,659,560 |
| Rejected aligned windows | 0 | 0 |
| Minimum weight in covered windows | 42 | 42 |

The exact receipts are `receipts/goal04_four_node_primary.json` and
`receipts/goal04_four_node_independent.json`.

## Improved cutoff

The smallest number of complete four-node blocks whose certified weight
reaches the required distance is

\[
\left\lceil\frac{188{,}766}{36}\right\rceil=5{,}244.
\]

These blocks occupy 20,976 nodes and contribute at least

\[
5{,}244\cdot36=188{,}784.
\]

One fewer block contributes only 188,748 under this argument. Therefore,
20,976 is the minimal cutoff obtained by partitioning the prefix into complete
four-node blocks.

For every \(r\ge20{,}976\), the second lap begins with at least 20,976
zero-input nodes. If \(L(y)\ne0\), those nodes already meet the distance
requirement. The later drive word cannot reduce their weight because it
occupies disjoint output positions.

## Ledger update

The support-33 terminal-zero row from Goals 01 and 02 already applies to every
placement. The newly closed support-33 event has \(L(y)\ne0\), so its bad-event
probability is zero. Goal 04 therefore adds 1,677 closed first-node strata
without adding a probability charge.

The disjoint partial upper bound remains

\[
[2^{-40.972830672896},2^{-40.972830672895}].
\]

The remaining numerical budget below \(2^{-40}\) remains

\[
[2^{-41.027690825931},2^{-41.027690825930}].
\]

The four-node certificate also applies to supports 35 through 39 when
\(L(y)\ne0\). However, charging the entire placement region after node 20,976
would also charge its unbounded terminal-zero part. That coarse union bound,
combined with the support-33 terminal-zero row, is about

\[
2^{-39.548488921864},
\]

which exceeds the target. The active ledger therefore retains the narrower
higher-support placement charges from Goal 03. The larger terminal-zero bands
for supports 35 through 39 remain open.

## Remaining obstruction

The first open event now consists of the 26 authenticated support-33 words
with

\[
L(y)\ne0
\qquad\text{and}\qquad
0\le r\le20{,}975.
\]

For an earlier first node, the available deterministic prefix bound is

\[
B(r):=max\left\{
25\left\lfloor\frac r3\right\rfloor,
36\left\lfloor\frac r4\right\rfloor
\right\}.
\]

Let \(W_{\mathrm{tail}}\) be the output weight outside the complete blocks used
by \(B(r)\). A bad output in this event must satisfy

\[
W_{\mathrm{tail}}\le188{,}765-B(r).
\]

This tail event is the next finite proof obligation. Goal 04 does not assign it
a probability bound.

## Conclusion

Goal 04 proves the requested boundary result and a stronger cutoff. The proof
is pointwise in every nonzero terminal state, so it avoids the correlation
created by conditioning on \(L(y)\ne0\). The complete support-33 complement and
the full construction remain open.
