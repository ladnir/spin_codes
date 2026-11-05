# Independent local output shuffles: layout assessment

This note assesses a proposed new outer ensemble, not an implemented optimization.
After `D2 H* D1` in each sixteen-row RS column, apply an independent uniform
permutation of its sixteen 128-bit output records. Sample each permutation once
at setup, independently for every source group and column. The regional paired
route and inner encoder remain unchanged logically. No native code was changed,
and no build or benchmark was run for this assessment.

## What the permutation could repair in the proof

For an arity-five affine-basis input set A, the current outside-weight-four
cap assigns count 255 to every four-subset of its eleven outside positions.
There are 330 such subsets, but exactly five occur. Their total count is 1275.
The five supports are the stars of K5 on ten outside positions; the eleventh
position is absent from all five supports. This is a factor 66 overcount at
equal markers, before any regional or whole-word power.

An independent output permutation makes any fixed projected support of size w
uniform among the C(16,w) output supports. For nonnegative markers at most 1,
dropping the other output coordinates is a valid pointwise majorization.
The permutation average of this projected monomial is then `e_w(t)/C(16,w)`.
Actual outside weight counts can replace maxima per labelled support in this
comparison. This statement does not prove a distance bound, and it does not
recover the inside coordinates already dropped by the projection.

An alternative under the unchanged encoder is to retain exact low-weight rays.
The new `theory/hd_bit_syndrome_basis_rays.py` recovers the complete sixteen-output
mask on each of the five outside-weight-four rays. Its initial three labelled
examples all have full output weight nine, not four. Those examples do not
establish that conclusion for every affine-basis A.

## Producer-side placement is not free

The current `bit_syndrome_native/HdGroupedOuter.h` kernel stores sixteen ZMMs
after mixing one column. Each ZMM contains four independent source groups,
one 128-bit record per group. Each store therefore fills one 64-byte cache line
in one regional plane. The hot loop makes sixteen such stores per four-group
column batch.

Four independently sampled output permutations generally send the four records
of a ZMM to different regional planes. A dynamic change to one store address
cannot implement them. Direct placement needs sixty-four 16-byte stores, plus
record extraction, instead of sixteen 64-byte stores. A scatter instruction does
not restore the existing contiguous non-temporal store behavior.

Alternatively, the sixteen registers can be rearranged before storing. That
requires four independent sixteen-way lane permutations across registers, not
one fixed transpose or a compile-time renaming. A group-major temporary layout
still needs this rearrangement or additional transpose work; it also changes
the consumer's locality.

Sharing one permutation across the four groups would preserve full-width stores.
It would also change the proposed ensemble. Bounds using independent group
permutation averages could not be reused without a new dependence argument.

## Exact absorption into route preprocessing

The cleanest implementation experiment would leave the producer unchanged and
compose the local output permutation into the setup-time route table.
For target row r in outer column j and original group g, the physical source
record would be

\[
\texttt{scratch}+
  \bigl(16j+\sigma_{g,j}^{-1}(r)\bigr)\texttt{stride}+g.
\]

Here g is the group index after the existing paired route and fair orientation,
including its circular boundary convention. Setup can calculate complete
offsets and emit them in the consumer's lane order. Encoding need not load a
second shuffle-key table or perform runtime permutation composition.

This exact composition loses the current contiguous-pair advantage. Adjacent
groups g and g+1 choose the same physical plane with probability 1/16. The
present `innerInput4` uses two 256-bit loads and one 256-bit insertion to form a
ZMM. Generic predecoded shuffled input requires four 128-bit loads and three
insertions. Fair pair orientation can be incorporated into the offsets, removing
the separate pair-swap shuffle, but it cannot restore contiguous loads.

At K20, the current decoded route uses eighty bytes per sixty-four-record packet:
thirty-two 16-bit pair bases and sixteen mask bytes. Complete 32-bit record offsets
use 256 bytes per packet. The hot route table grows from 2.5 MiB to 8 MiB.
All offsets fit within 32 bits; ordinary 16-bit record offsets do not cover the
sixteen possible source planes.

Useful outer payload remains 32 MiB. However, one logical region now reads across
sixteen physical planes, an approximately 2 MiB address range instead of 128 KiB.
The current sequential next-plane prefetch no longer describes the access pattern.
Cache-line reuse and key traffic, not additional GFNI arithmetic, are the likely
costs. These are source-level consequences, not a measured regression.

## Decision boundary

Do not describe independent S16 output shuffles as a zero-cost producer rename.
If their proof gain warrants one native experiment, use precomputed individual
record offsets first. It preserves the optimized outer arithmetic and stores,
isolates the consumer cost, and implements the intended independent ensemble.
Keep a shared-permutation variant separate: it has different probabilistic
dependencies even if its native layout is cheaper.

Before that experiment, test the exact low-output rays under the unchanged
construction. They address both the impossible-support overcount and the dropped
inside coordinates without introducing a new route layout.
