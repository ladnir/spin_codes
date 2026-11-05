# Exact balanced-update drift of the auxiliary hole process

2026-11-01. This is a diagnostic for the global rank proof. The encoder
and its routing distribution remain unchanged. In particular, actual
routing fixes a quota over sixteen updates, not separately in each update.

Let K_x be the adaptive population operator at y=2, with physical state
(a,w): a is the number of future holes and w their persistent overlap.
For this diagnostic only, define

```
P = [x^256] K_x / binom(512,256),
Q(a,w) = P 1
       = [x^256] (2+x)^w (1+x)^(512-w) / binom(512,256).
```

Normalize each row by Q(a,w). This gives a probability distribution on
the next auxiliary population. It is a moment-weighted distribution,
not the distribution of encoder states or of unweighted setup samples.

The [exact rational checker](quartet_quota_drift.py) first samples the
old colors in one quartet by the multivariate hypergeometric count.
Conditional on those colors and the local input labels, the remaining
508 positions contribute their exact quota coefficient. If u positions
in the quartet belong to the old hole set, o belong to its persistent
overlap, j current labels are outside H, and r of these j labels belong
to the old overlap, the hole probability is

```
g = (2^max(0,4+u-2j+r)-1)/(2^(4+u)-1).
```

For j=0 this is one. The expected next populations are 512 E[g] and
128 E[u g]. All local probabilities, normalization checks, and means
use exact fractions. The following decimals only display those values.

| Old (a,w) | log2 Q, display | Mean next a | Mean next w |
| --- | ---: | ---: | ---: |
| (64,8) | 4.671 | 62.201 | 8.977 |
| (128,32) | 18.562 | 67.795 | 19.399 |
| (128,48) | 27.721 | 69.226 | 20.423 |
| (160,64) | 36.797 | 72.255 | 26.226 |
| (192,96) | 54.698 | 76.825 | 33.468 |
| (256,128) | 72.261 | 82.685 | 46.250 |
| (512,512) | 256 | 160.415 | 160.415 |

The parent ran the fixed ten-point diagnostic serially; every exact
normalization equaled one. Independent mathematical and source review
passed. There are no empirical samples or fitted transition laws here.

The expensive populations in this balanced auxiliary process do not
retain their high overlap in the first moment. This rules out treating
the tested high-overlap band as approximately stationary on the strength
of its row mass alone. It is not a concentration theorem: rare paths can
still dominate later weighted moments. Nor may P replace the actual
regional operator, which permits different label counts in its sixteen
updates as long as their sum is 4096.

The next comparison should therefore retain the remaining regional quota
alongside the population. A proposed finite envelope in that quota is
being checked on the restricted three-state operator before any further
full-population control. The full distance certificate remains open.
