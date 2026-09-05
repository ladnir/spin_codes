# FrozenRandom512-P2 initial report

The random linear constituent removes RM from the structural model. Each
data-group encoder is sampled during setup and then frozen. The model uses the
exact unconditioned ensemble spectrum. A random constituent typically has
minimum distance near 58, but the model imposes no local distance floor.

The first wide-interval calculation gives `+3003.95` bits. This number is not
a current obstruction. The same interval machinery previously lost about
10779 bits by reusing the support-16384 inner tilt through support 32767.

At data support 95, the random constituent has relaxed outer coefficient
`+228.9165`. The pointwise inner upper is `-197.2281`, giving `+31.6883` if
the two determined parity symbols are discarded. This reproduces the old RM
diagnostic because RM was already random-like near this moderate support.

The exact one-active-group calculation retains both field parities. The two
parity equations have rank two over \(\mathbb F_{2^{64}}\) on a four-symbol
group. Hence their output pair is uniform as the group input varies. Combining
this fact with the exact extended-BCH \([128,64,22]\) spectrum gives

\[
\log_2 \mathbb E[Z_{\mathrm{one\ group}}]\le -91.830594.
\]

The largest individual profile has data packet support 95 and parity support
zero. Its outer coefficient is `+100.9165`; its inner upper is `-197.2281`;
the combined exponent is `-96.3117`.

Thus the previous one-group obstruction was caused by dropping the parity
correlation. The next unresolved sparse case has at least two active data
groups. Their contributions to the global parity symbols can cancel. A
parity ladder should therefore be evaluated against multi-group cancellation,
not against the already closed one-group row.

## Global pointwise retry

The global retry retains both field parity symbols and sums every outer
occupation inside one moment. The scan covers every support through 512,
every eighth support through 4096, and checkpoints through 32767.

The largest sampled point occurs at total packet support 360:

```text
outer coefficient Chernoff upper:  +751.9332
pointwise inner upper:              -836.3347
combined exponent:                   -84.4015
```

The tilted occupation mode is three, and the mean is 3.75 active data groups.
Each active group contributes about 94.20 packets. The two parity blocks
contribute about 6.94 packets in total. The inner profile uses eleven
terminations.

The old positive global result was therefore an artifact of two relaxations:
wide inner intervals and discarded parity correlations. If a complete cover
keeps every support below `-84.40`, summing all 524352 supports gives at most
`-65.40`. This would leave 25.40 bits beyond the `2^-40` target.

The present scan is not that certificate. Unscanned supports and the floating
inner evaluation still require a rigorous cover.

The dense random encoders are a proof device. This report makes no efficiency
claim for their implementation.
