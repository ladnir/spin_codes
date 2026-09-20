# Finite conditioning-window revision at \(k=2^{20}\)

## Result

The first finite setup should condition each independent \(B=240\) BA row on
having no nonzero word outside

\[
  [23,217],
  \tag{1}
\]

not the earlier interval \([25,215]\). This changes only the setup acceptance
test. It does not change the online encoder, its stored representation, the
shortening rule, or the RM2Sub recurrence.

The outward BA enumerator and Markov's inequality prove

\[
  \Pr[\mathcal G_{240}^{23}]
  \ge \mathtt{0x1.fcc11f1547a88p-1}.
  \tag{2}
\]

The corresponding expected number of rejection trials is at most
\(1.00638\). The decision procedure for the event in (1) remains open.

## Outward sparse receipts

For the conditional law in (1),
`golay_ba3_rm2sub_finite_B240_q1_outward_w23_217_k20_d11.json` proves

\[
  \mathbb E[Z_{233164,1}]
  \le \mathtt{0x1.6f3c66666f370p-47}<2^{-46}.
  \tag{3}
\]

The displayed margin is 46.479 bits. The separate receipt
`golay_ba3_rm2sub_finite_B240_q2_64_outward_w23_217_k20_d11.json` proves

\[
  \mathbb E\!\left[\sum_{Q=2}^{64}Z_{233164,Q}\right]
  \le \mathtt{0x1.4dc4b8b530941p-1}\,2^{-84}<2^{-84}.
  \tag{4}
\]

Equations (2)--(4) are outward results, subject to the already recorded O6
source-interface audit. They are not a complete distance certificate.

## Why the revision is necessary

At \(Q=L=8832\), take the central BA shell of weight 120 and compare it with
iid Bernoulli-half input to the finite three-state RM2Sub transfer. Nearest
binary64 evaluation gives the following diagnostic margins:

| quantity | margin |
|---|---:|
| uniform random rate-half block outer | 180.464 bits |
| unconditioned BA central shell | 233.666 bits |
| BA central shell conditioned on \([25,215]\) by the available Markov bound | -365.245 bits |
| BA central shell conditioned on \([23,217]\) by the available Markov bound | 152.637 bits |

The old event has proved probability only \(0.9540842447\). Applying its
conditional-expectation factor independently in all 8,832 rows costs 598.911
bits at \(Q=L\). This loss exceeds the finite endpoint slack. The revised
event has proved probability \(0.9936608995\), so the same calculation costs
81.029 bits.

This comparison proves that the old conditional-spectrum *bound* cannot close
the central endpoint. It does not prove that the old conditional ensemble has
bad distance; the factor \(1/\Pr[\mathcal G]\) may be loose for central
weights.

## Window sweep

`diagnose_finite_k20_conditioning_window.py` and its JSON receipt compare all
symmetric windows with lower endpoint 18 through 25. The binary64 sweep shows
that \([23,217]\) balances the two visible bottlenecks: its one-active margin
is 46.497 bits diagnostically and its central endpoint margin is 152.637 bits.
The outward value in (3) is slightly smaller, as required.

The wider \([19,221]\) window leaves only 0.444 diagnostic bits beyond the
40-bit target in the one-active class. The narrower \([24,216]\) window leaves
only 3.533 bits beyond the target at the central endpoint. Neither is the
preferred starting point for the remaining union.

## Dense status

The full-spectrum single-reference Rényi relaxation remains insufficient. On
the revised law, its sampled margins are negative at \(Q=3072,4096,8704\),
and 8832. At \(Q=8832\), the failure is about 49,750 bits. This is an
endpoint-shell likelihood artifact: the same relaxation forces one high
Rényi reference to cover all permitted BA shells simultaneously.

The dense proof must therefore retain at least a finite band or shell label.
The next proof object is a positive multitype transfer or an equivalent
coefficient extraction that sums mixed shell compositions without charging a
single endpoint likelihood to every row. It must exploit the 152-bit central
endpoint slack sharply; a factor such as \((B+1)^Q\) is too expensive at this
fixed block size.

## Proof status

- **Proved outward:** (2), (3), and (4), within the stated mathematical
  probability space and conditional on the O6 source audit.
- **Exact in form, binary64 only:** the central-shell transfer comparison and
  the conditioning-window sweep.
- **Open:** every occupation from 65 through 8,832, the efficient test for
  (1), the RM2Sub implementation-equivalence audit, and the final sum over all
  occupations.
