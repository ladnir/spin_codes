# Riffle Outer256Model-RandomStepConv g=4 sigma=20

**Status:** `ACTIVE_MODEL_EXPLORATION`

This folder records a fixed-rate outer-code diagnostic. It does not specify a
deployable code.

The model pairs the 16384 data symbols in \(\mathbb F_{2^{64}}\). It applies
one rate-one-half \([256,128,38]\)-shaped binary constituent to each pair.
The 8192 data pairs and one parity pair produce

\[
8193\cdot64=524352
\]

four-bit packets. Thus the model preserves the packet count, overall rate,
global permutation, and one-lap RandomStepConv inner.

The provisional spectrum is complement symmetric and random-like. It has
zero multiplicity below weight 38. The evaluator also inflates its low-weight
window by 10, 20, and 40 bits.

The model must not be promoted to a construction until an explicit linear
code and an adequate spectrum bound replace the guessed spectrum.

## Initial result

The baseline 9% first-moment exponent is \(+4710.17\), compared with
\(+8097.18\) for the 128-bit constituent. This is a substantial improvement,
but it does not close the first moment. Raising only the modeled minimum
distance through 68 has almost no effect because weights near 100 drive the
bulk saddle.

The saved canonical test event is late placement without inner termination.
A representative modeled word has 32748 active packets. They all land in the
last 100477 positions, and the state remains live after activation. The
no-termination condition costs only 0.14 bits of exponent. This event does not
witness a failure: after including the relaxed modeled outer coefficient, its
expected-count exponent is \(-6901.48\).

See `REPORT.md`, `receipts/fixed_rate_outer256_model.json`, and
`receipts/fixed_rate_outer256_distance_sweep.json`. The outer multiplicity
audit is in `receipts/outer_multiplicity_audit.json`.
