# Goal 02: global pointwise cover

## Result of the diagnostic scan

The global outer moment retains both field parity symbols. It uses the exact
ensemble-average packet spectrum of each independently sampled data encoder.
It also uses the exact extended-BCH spectrum for both parity symbols.

The scan evaluates every total support through 512, every eighth support
through 4096, and checkpoints through 32767. Its largest sampled point is

\[
h=360,
\qquad
\log_2 A_h\le 751.9332,
\qquad
\log_2 P_{\mathrm{inner}}(h)\le -836.3347.
\]

The combined exponent is `-84.4015`.

The tilted outer profile has occupation mode three and mean occupation 3.75.
Its median occupation is four. Each active data group contributes about 94.20
packets under the tilt. The two parity blocks contribute about 6.94 packets
in total. The inner upper is dominated by eleven terminations.

## Proof target

Convert the scan into a cover of every support from 1 through 524352. It is
enough to prove that no unscanned support exceeds `-84.40`. Summing this bound
over all supports would give `-65.40`, leaving 25.40 bits beyond the target
`2^-40` first moment.

The cover must retain the exact parity moment. It must not return to a
worst-case parity charge or the old wide support intervals.

## Scope

The current outer coefficient is a Chernoff upper. The pointwise inner values
are floating diagnostics based on the verified fixed-support identity. The
sampled curve is strong evidence, but it is not yet a complete certificate.
