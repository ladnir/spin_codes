# State-versus-parity report

## Current best point

Four parity symbols and 16 state bits give the lowest setup entropy found in
the fixed 512-bit, four-bit-packet family. The setup entropy is

| Setup object | Entropy in bits |
|---|---:|
| 4096 independent FieldPair512 maps | 2,097,152 |
| 524416 independent 20-bit scalar inner maps | 10,488,320 |
| Uniform global packet permutation | 9,207,527.151 |
| Four uniform parity-block permutations | 2,864.647 |
| Total | 21,795,863.798 |

The total is approximately 2.598 MiB. It saves 1,044,775.660 bits relative
to the P2, sigma-18 candidate.

The exact-saddle computation over supports 900 through 3000 gives total
exponent `-44.905837`. The largest exact point in that interval is
`-53.380379` at support 1813. The scan outside the interval is lower, but the
tails have not yet been converted into a formal cover.

## Why sigma 15 does not follow from more parity

Five parity symbols with 15 state bits have sampled maximum `+87.981540` at
support 288. Six parity symbols with 15 state bits have sampled maximum
`+79.371762` at support 320.

These are not late-start profiles with many active outer blocks. One active
512-bit data block can fill the data output and all parity outputs. That
single block represents about 256 message bits, so the family has very high
multiplicity. Increasing the parity count lengthens the packet support, but
a 15-bit inner state does not make its low-output trajectories rare enough to
offset that multiplicity.

Thus extra global parity is no longer an efficient way to buy a smaller
inner. A sigma-15 candidate needs a structural change, such as a different
outer constituent geometry, a wider packet, or a stronger inner transition
law. It should not be described as a tuning change to this candidate.

## Evidence scope

The P4 number is diagnostic. It relies on an exact outer coefficient
calculation in the central interval and a floating point inner bound. A
complete `2^-40` claim requires a support-tail cover and outward numerical
rounding.

The P5 and P6 results are refutation evidence for the present first-moment
model. Their positive margins are large enough that numerical rounding is not
the issue.
