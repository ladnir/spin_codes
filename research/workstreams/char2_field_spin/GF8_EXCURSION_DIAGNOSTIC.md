# Excursions in the unpunctured three-quarter-rate bound

The existing sparse transfer already records every transition between zero and live state.
This diagnostic asks which transitions dominate its failing upper bound for the full-Cauchy width8/state8 packet.
It does not sample codes or exhibit low-distance codewords.

For each path in the positive transfer, let B count births, R count returns, and T count packets with live incoming state.
The path begins at zero.
It therefore satisfies B-R=1 when its final state is live, and B-R=0 otherwise.
Births and returns require nonempty inputs, so B+R is at most the total raw input weight W.
Both constraints already hold term by term in the existing regional polynomial.
Marking births by exp(theta) and returns by exp(-theta) is only a diagonal similarity transform.
It changes the terminal live contribution, not the internal bound.
Consequently these counters alone do not improve the proof.

The helper differentiates positive marked path sums at the optimized output/rank marks.
Its means describe the tilted upper-cover distribution, not the distribution of an actual code.
At K20, outer [128,96], and excluded distance 10%, it gives:

| Active rows | Bound margin | Births | Returns | Live packets | Raw input weight |
|---:|---:|---:|---:|---:|---:|
| 128 | -2312.06 bits | 80.38 | 79.38 | 21819.62 | 6295.07 |
| 256 | -4795.49 bits | 94.90 | 93.90 | 22204.71 | 12225.71 |

There are 175104 packets in the full traversal.
These witnesses spend approximately one eighth of their time live.
They have approximately 79 or 130 input symbols per completed excursion, not isolated one-input cancellations.
Their mean live duration per birth is approximately 271 or 234 packets.
Thus a special treatment of very short excursions would not directly target these witnesses.

One potentially useful coupling remains unproved.
Fix a union S of packet intervals occupied by excursions.
Every active outer row must have at least 33 routed coordinates available in S, because the outer's minimum weight is 33.
At a uniform 12.5% occupation fraction, an individual row has only about 16 of its 128 routed coordinates available on average.
A row-availability tail could therefore be useful.
However, S is selected by the message and setup, and its interval-position entropy cannot be ignored.
A valid proof must combine that tail with the chronological transfer or union bound, rather than condition on a conveniently chosen S.
This is a candidate proof direction, not an additional certified factor.

`theory/gf8_excursion_screen.py` reproduces the diagnostics.
`theory/gf8_excursion_tests.py` checks path constraints, the similarity identity, and its effect after regional composition.
All three tests pass.
