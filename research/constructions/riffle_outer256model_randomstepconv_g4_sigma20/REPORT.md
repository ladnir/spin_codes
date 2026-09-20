# Fixed-rate Outer256 model report

## Result

The fixed-rate \([256,128,38]\)-shaped model improves the 9% first-moment
diagnostic, but it does not close it. The baseline modeled spectrum gives

\[
\log_2\mathbb E[Z_{0.09}]\le 4710.17.
\]

The corresponding diagnostic for the exact \([128,64,22]\) spectrum was
\(8097.18\). Doubling the constituent block size therefore recovers about
3387 bits within the current bounding method.

The dominant packet-support interval remains 16384--32767. Its occupation
mode changes from about 1401 active 64-bit blocks to about 691 active 128-bit
blocks.

These values are floating diagnostics. A positive first-moment upper bound
does not refute the modeled construction.

## Canonical late-placement test

The baseline saddle has a concrete late-placement interpretation. It contains
about 691 active data-pair blocks. Under the saddle tilt, one active block has
mean binary weight 101.10 and mean packet support 47.39. The joint modal local
profile has binary weight 100 and packet support 47. A representative outer
word therefore has packet support

\[
h=32748.
\]

Fix any outer word with this packet support. Consider the following event over
the global packet permutation and the independent inner maps.

1. All active packets land in the final \(R=100477\) packet positions.
2. After the first active packet, the inner state remains nonzero.
3. The emitted binary weight is at most \(D=188766\), or 9% of the output.

The first condition forces 423875 leading packet positions to emit zero. The
probability of the complete event is at least \(2^{-86331.95}\) for the fixed
outer word. Its exponent has the following decomposition:

| Component | Logarithm of probability |
|---|---:|
| Late placement | -85259.39 bits |
| Low output weight in the live suffix | -1072.42 bits |
| No state termination after activation | -0.14 bits |

If all \(R\) suffix steps were live, their mean relative binary weight would
be 9.581%. Thus the output condition asks for only a modest downward deviation
from the mean.

This event requires no inner termination. The no-termination condition costs
only 0.14 bits within this test.

The outer profile comes from a real-valued spectrum model. It is not an
explicit codeword witness. The late-suffix calculation is a valid per-word
probability lower bound for any fixed word with packet support \(h\).

The relaxed modeled outer coefficient at \(h=32748\) is about
\(2^{79430.47}\). Multiplying this coefficient by the test-event probability
gives an expected count of \(2^{-6901.48}\). Thus this event does not witness
a distance failure, even in the modeled spectrum. It only supplies a concrete
late-placement mechanism against which to compare the upper bound.

## Outer multiplicity audit

The evaluator does not assign minimum weight 38 to every active block. Let
\(A_w\) denote the modeled multiplicity at binary weight \(w\), and let
\(P_{w,s}\) denote the exact probability that a weight-\(w\) word occupies
\(s\) packets after its local permutation. The data-block generating function
is

\[
F(x)=1+\sum_{w>0}\sum_s A_wP_{w,s}x^s.
\]

Expanding \(F(x)^{8192}\) sums every block-weight and packet-support
configuration. Each term includes its multinomial placement count and the
product of its local multiplicities. This summation is exact relative to the
guessed spectrum.

An unconditioned nonzero modeled block has mean binary weight 128 and mean
packet support 60. Its central 98% weight range is 110--146, and 97.0% of its
mass has weight at least 114.

Conditioning toward the dominant low-output saddle changes this distribution.
The conditioned mean weight is 101.10, and its median is 102. The conditioned
mean packet support is 47.39. About 78.7% of the conditioned mass lies at
weights 90--112. Only 0.077% lies at weight 70 or below. The saddle therefore
selects moderately lighter but still highly numerous words. It does not use
minimum-weight words as typical blocks.

The support tilt also selects about 691 active blocks. The mean is 691.39,
with standard deviation 25.16. Its 1%--99% range is 633--751 active blocks.

At representative support \(h=32748\), an FFT coefficient extraction gives
the following comparison:

| Quantity | Exponent |
|---|---:|
| Relaxed modeled outer coefficient | 79430.47 bits |
| Inner probability upper bound | -74747.57 bits |
| Pointwise first-moment upper bound | 4682.90 bits |
| Whole support-segment upper bound | 4710.17 bits |

The support-segment machinery loses only 27.26 bits relative to the pointwise
calculation. At this support, 8.85 bits come from the secant bound for the
reciprocal binomial coefficient. Another 9.37 bits come from the support-tail
Chernoff bound.

The many-block branch does make one block-level worst-case relaxation: it
discards the determined parity-pair block when the support tilt is below one.
Treating that block heuristically as an ordinary nonzero modeled block would
recover about 131 bits. Correlations prevent using this heuristic as a proof,
and one parity block cannot explain the remaining thousands of bits.

The outer bulk sum is therefore using the modeled spectrum's multiplicity.
The unresolved gap at high support is mainly between the inner probability
upper bound and concrete inner lower events. At \(h=32748\), the explicit
late-suffix event has exponent \(-86331.95\), about 11584 bits below the inner
upper bound.

## Model

The model groups the 16384 data symbols into 8192 ordered pairs. One modeled
rate-one-half binary constituent maps each 128-bit pair to 256 bits. The two
field parity symbols form one additional input pair. Thus the outer produces

\[
8193\cdot256=2097408
\]

bits, or 524352 four-bit packets. The one-lap inner is unchanged.

The provisional local spectrum has these properties:

- one zero word and one all-one word;
- complement symmetry;
- no nonzero weights below 38; and
- random-like even-weight multiplicities from weight 38 through weight 218.

The spectrum has real-valued multiplicities. It is a moment model, not the
weight distribution of an explicit code.

## Parity treatment

For one active data pair, the two field parity equations define an invertible
map from that pair to the parity pair. Cauchy's inequality bounds the product
of their local moments by the second local moment.

For two or more active data pairs, the calculation discards the parity block
when the support tilt is at most one. For larger tilts, it charges the maximum
64-packet support. This step is a valid worst-case bound relative to the
modeled local moments.

The bulk saddle contains about 691 data blocks. A sharper treatment of one
parity block cannot plausibly supply the missing 4710 bits by itself.

## Sensitivity to the modeled spectrum

The evaluator inflates the multiplicities in weights 38--70 and their
complements. It then renormalizes the nonzero spectrum mass.

| Low-weight inflation | First-moment exponent | Dominant support | Occupation mode |
|---:|---:|---:|---:|
| 0 bits | 4710.17 | 16384--32767 | 691 |
| 10 bits | 5348.48 | 16384--32767 | 773 |
| 20 bits | 17897.87 | 32768--65535 | 1819 |
| 40 bits | 99367.26 | 131072--262143 | 7622 |

The calculation is therefore sensitive to deviations from the random-like
bulk spectrum.

## Minimum-distance sweep

Increasing only the modeled minimum distance does not address the dominant
term.

| Modeled minimum distance | First-moment exponent |
|---:|---:|
| 38 | 4710.168 |
| 44 | 4710.168 |
| 50 | 4710.168 |
| 56 | 4710.167 |
| 62 | 4710.157 |
| 68 | 4710.009 |

At the dominant support tilt, the local moment has mean binary weight about
101.1. Its largest contributions come from constituent weights 90--112.
Removing additional words near the minimum distance therefore has negligible
effect on the bulk saddle.

## Interpretation

Increasing the constituent block size is a better use of computation than an
unretained second inner lap. It preserves the packet count and recovers
thousands rather than tens of exponent bits. However, the 256-bit model does
not close the proof under the current inner upper bound.

The next diagnostic should test a 512-bit rate-one-half block model. It must
report the same canonical-event fields: active-block count, local weight,
packet support, placement window, output deviation, and required termination
count. It should also identify the bulk-spectrum envelope near its optimizing
local weights. Searching only for a larger minimum distance is unlikely to
help.
