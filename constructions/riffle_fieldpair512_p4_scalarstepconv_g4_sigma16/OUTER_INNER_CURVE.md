# Outer-constituent size versus inner state

## Comparison rule

The comparison fixes a (2^{20})-bit message, rate-one-half independently
sampled random outer constituents, four-bit packets, one global packet
permutation, and a one-lap scalar random convolution. Changing the constituent
size changes the number of independent outer blocks but not the total number
of data-output packets.

The clean control curve fixes two field-parity symbols. P2 is surjective for
outer lengths 256, 512, 1024, and 2048 bits, so its one-active-block parity
law is available exactly.

| Outer length | Outer dimension | P2 inner result | Evidence |
|---:|---:|---:|---|
| 256 | 128 | No certified sigma | Sigma 32 exact central sum `-38.848659`; the zero-termination profile is already near its infinite-state limit |
| 512 | 256 | Sigma 18 | Exact evaluated sum `-43.598633`; sigma 17 central interval is positive |
| 1024 | 512 | Sigma 17 | Exact central sum `-44.042347`; sigma 16 sampled maximum `+13.984046` |
| 2048 | 1024 | Sigma 17 | Exact central sum `-57.715141`; sigma 16 sampled maximum `-10.273467` and therefore cannot meet `2^-40` |

These are diagnostic thresholds. The finite support intervals and numerical
rounding still require certification.

The curve is not monotone by a simple one-bit-per-doubling rule. Moving from
512 to 1024 bits buys one state bit. Moving from 1024 to 2048 bits buys margin
but no additional integer state reduction.

## Active P4 comparison

The active construction uses four parity symbols. Its currently reliable
points are:

| Outer length | Inner state | Result |
|---:|---:|---|
| 512 | 16 | Exact central sum `-44.905837` |
| 1024 | 15 | Sampled maximum `-7.394418`; fails |
| 1024 | 16 | Sampled maximum `-112.912186`; comfortably below the pointwise target |

Thus doubling the P4 outer constituent from 512 to 1024 bits gives a large
margin improvement but does not reduce the minimum observed state below 16.
Because the total FieldPair setup entropy is the outer length times the block
count, both outer sizes use the same 2,097,152 outer random bits. They also use
the same inner and permutation entropy when sigma is fixed. The 1024-bit point
therefore buys proof margin rather than setup-entropy savings, while likely
increasing field-arithmetic cost.

## Why P4 at 256 bits is unresolved

A 256-bit rate-one-half constituent has only two 64-bit input symbols. Four
parity symbols are then an MDS image of those two symbols, rather than a
surjective four-symbol output. The present safe bound retains only the minimum
packet support of each nonzero parity symbol. It produces a one-active-block
maximum of `-27.578326` at sigma 24, but that number is too loose to interpret
as the actual threshold.

Computing the P4/256 point requires a joint packet-support moment for four
correlated extended-BCH encodings of an MDS `[4,2,3]` field word. Until that
moment is derived, P2/256 is the trustworthy indication: 256-bit constituents
are at best very close to the viability boundary, even with a much larger
inner state.

## Current design conclusion

For the P4 construction, 512 bits is the smallest outer size presently shown
to support a 16-bit inner. A 1024-bit outer offers much more proof margin at
the same information-theoretic setup entropy, but it does not lower sigma.
The next useful calculation is the exact P4/256 joint-parity moment. It decides
whether a smaller and potentially cheaper outer can compensate with a modest
increase in state, or whether 512 bits is a genuine constituent-size floor.
