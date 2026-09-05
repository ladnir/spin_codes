# Cancellation analysis

Two impulses in one accumulator lane have even terminal parity. The
checkpoint cannot revive that lane when the other 63 terminal parities are
also zero. This note measures that event under the permutations in the
construction.

## Local cancellation

An epoch contains 64 lanes with 16 positions per lane. Condition on two
uniform impulse positions in one epoch. The zero state remains zero exactly
when both positions use the same lane. This event has probability

\[
 \frac{64{16\choose2}}{{1024\choose2}}
 =\frac{15}{1023}
 \approx 2^{-6.09}.
\]

A complete region contains eight epochs. Equivalently, it contains 512
epoch-lanes with 16 positions each. Condition on two uniform positions in
the region. The state remains zero throughout the region with probability

\[
 p_2:=\frac{512{16\choose2}}{{8192\choose2}}
 \approx 2^{-9.09}.
\]

Conditioned on this event, the output weight equals the distance between the
two positions within their lane. Its conditional mean is \(17/3\). The
positions are consecutive with conditional probability \(1/8\). Therefore,
the weight-one cancellation event has unconditional probability

\[
 \frac{512\cdot15}{{8192\choose2}}
 \approx 2^{-12.09}.
\]

Odd region weight always activates the state. For the first few even
weights, the complete-region nonactivation probabilities are:

| Region weight | Probability | Negative log probability |
|---:|---:|---:|
| 2 | \(1.8313\times10^{-3}\) | 9.09 bits |
| 4 | \(1.0051\times10^{-5}\) | 16.60 bits |
| 6 | \(9.1852\times10^{-8}\) | 23.38 bits |
| 8 | \(1.1740\times10^{-9}\) | 29.67 bits |

Thus cancellation becomes rapidly less likely as the occupancy grows.

## Two minimum-weight outer words

Now condition on two outer words of weight 38. Their independently permuted
supports are independent uniform 38-subsets of the 256 regions. Their
intersection has mean 5.640625 and mode five. At the modal intersection,
the two words produce five overlap regions and 66 singleton regions.

Every singleton region activates the zero state. An overlap region fails to
activate with probability \(p_2\). Averaging over the exact hypergeometric
intersection law gives

\[
 \Pr[\text{first occupied region fails}]
 =1.4847\times10^{-4}
 \approx2^{-12.72}.
\]

Repeated failures lose about 12.7 bits each:

| Leading failed occupied regions | Probability | Negative log probability |
|---:|---:|---:|
| 1 | \(1.4847\times10^{-4}\) | 12.72 bits |
| 2 | \(2.1762\times10^{-8}\) | 25.45 bits |
| 3 | \(3.1488\times10^{-12}\) | 38.21 bits |
| 4 | \(4.4966\times10^{-16}\) | 50.98 bits |

For a stronger late-start test, fix a prefix containing \(r\) complete
regions. The exact nonactivation probability sums over the intersection size
and the number of overlap regions in the prefix. It excludes every prefix
containing a singleton and charges \(p_2\) for every overlap in the prefix.

| Prefix regions | Nonactivation exponent | Exponent without cancellation | Loss from cancellation |
|---:|---:|---:|---:|
| 64 | 34.52 bits | 34.53 bits | 0.01 bits |
| 128 | 85.27 bits | 85.33 bits | 0.06 bits |
| 192 | 183.47 bits | 184.43 bits | 0.96 bits |
| 209 | 236.32 bits | 241.87 bits | 5.55 bits |
| 210 | 240.23 bits | 246.64 bits | 6.40 bits |

At distance 188743, region 210 is the first complete-region boundary after
which a permanently live state has remaining mean output at most the
distance threshold. The modeled number of weight-38 two-block messages has
base-two logarithm 73.55. Subtracting this multiplicity from the region-210
nonactivation exponent leaves 166.68 bits of margin before bounding the
output-weight tail.

The event that the two words never activate the state is smaller still. Its
modeled union margin is 423.25 bits.

## Interpretation

Pair cancellation is a genuine local low-weight event. The outer transpose
prevents it from becoming the evident global obstruction. Two outer words
usually create many singleton regions, and the first singleton activates the
state with certainty.

This calculation does not complete the two-active distance bound. It tracks
activation, not the tilted output-weight distribution after activation. It
also conditions on two outer words of weight 38. Other outer weights remain
to be summed.

The next proof goal should compute the exact two-active Chernoff bound using
the region matrices \(R_0(z)\), \(R_1(z)\), and \(R_2(z)\). That calculation
will determine whether the local running-parity weight introduces a failure
that the activation diagnostic does not see.
