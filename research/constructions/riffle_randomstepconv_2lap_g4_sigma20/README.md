# Riffle RandomStepConv-2Lap g=4 sigma=20

**Status:** `PAUSED_LOW_COMPUTE_VALUE`

This candidate is paused. A second convolution pass adds almost one full
inner traversal but improves the measured 9% exponent by only about 10 bits.
The next proposed direction is a stronger or doubled outer construction. Its
exact map has not yet been specified or named.

This candidate uses one unretained burn-in lap followed by one retained lap.
Both laps use independent random step maps. The packet order is unchanged
between laps, and the state crosses the wrap.

The candidate addresses the measured obstruction for the one-lap version.
At 9% relative output weight, that version was nearly saturated by packet
supports contained in a late suffix. The burn-in lap converts this special
boundary into an ordinary gap between active packets.

## Current proof route

For fixed packet support, decompose the retained lap into zero-input gaps.
Every active packet refreshes the state. The off-count generating function
therefore factors over the gaps. Averaging over a uniform support becomes one
univariate coefficient extraction over weak compositions of \(N-h\).

The first goal has four parts:

1. derive the exact gap generating function;
2. verify it against exhaustive small instances;
3. obtain a two-tilt bound for retained weight at fixed support \(h\); and
4. evaluate the full-size 9% saddle region from the one-lap calculation.

The outer BCH moment machinery from Goals 03 and 04 of the parent candidate
remains available. It cannot be transferred until the new inner bound has a
certified all-support envelope.

See `proof/GOAL_01_CYCLIC_LIVE_COVERAGE.md` for the active statement.

## Initial Goal 01 result

The exact identity passes exhaustive small checks. The first full-size
floating diagnostic is negative for 9%: at \(h=16384\), the upper exponent is
about \(-41480\), only 10 bits better than the one-lap value. An explicit
cyclic-cluster event gives a lower exponent near \(-41755\). See
`proof/GOAL_01_CYCLIC_LIVE_COVERAGE_REPORT.md` and
`receipts/goal01_cyclic_live_coverage.json`.

The new obstruction is one long cyclic gap followed by a state reset. Typical
supports remain well mixed, but the rare clustered-support event dominates
the proof-facing bound.
