# Random-convolution distance benchmark

## Result

The random convolution succeeds with state size 16. A larger state is not
needed for this benchmark.

Let the message length be \(2^{20}\), and let the encoded length be
\(N=2^{21}\). At distance

\[
D=230730,
\qquad
D/N=0.1100206375,
\]

the combined floating-point first-moment margin is 40.9510 bits. At
\(D=230731\), the margin is 38.3027 bits.

The rate-one-half Gilbert--Varshamov distance is approximately
0.1100278644. The finite-length random-code 40-bit endpoint is
\(230735/N=0.1100230217\). The benchmark endpoint is five output bits below
that random-code endpoint.

## Regime audit

At relative distance 0.11, the uniform-nonzero outer model has 42.9001 bits
of aggregate margin. One active outer block is limiting. The all-active
occupation has 174.4350 bits of pointwise margin.

At \(D=230730\), the all-active occupation becomes limiting. Its pointwise
margin is 41.3880 bits. Summing the neighboring occupations leaves 40.9510
bits.

The exact one-active calculation with the modeled even `[256,128]` spectrum
has 53.5612 bits at \(D=230730\). Its dominant outer weight is 38. Thus the
modeled BCH tail does not limit the endpoint.

The random inner removes the zero-state obstruction. When an epoch input and
state are nonzero, a terminated next state does not force a low output. The
output remains uniform after conditioning on termination.

## Scope

The full occupation calculation uses a uniform nonzero 256-bit word for
each active outer block. This is the ideal random-linear outer model. The
one-active calculation separately uses the modeled even spectrum with
minimum weight 38.

The middle and dense occupations do not yet track the even parity of each
modeled outer block across all 256 regions. Their margins are therefore
benchmark evidence, not an outward-rounded certificate for one explicit
BCH code. The dense endpoint is largely insensitive to the outer spectrum
because essentially every epoch is live and emits a uniform random block.

The two-state transfer identities are exact. The reported optimization uses
nearest binary64 arithmetic and two Chernoff grids. The sparse grid covers
log-surprisal values from -10 to 1 in increments of 0.5. The dense grid
covers 0.66 to 0.82 in increments of 0.01. Each occupation uses its better
candidate.

## Reproduction

The full occupation calculation uses
`scripts/analyze_riffle_transpose_bitshuffle.py`. The modeled-spectrum
one-active calculation uses
`scripts/analyze_riffle_spectrumperm_bitshuffle_oneblock.py`.

The next question is how to emulate the random epoch law efficiently. A
structured replacement must randomize the output even when the next state
is zero. Randomizing only the state update does not meet this requirement.
