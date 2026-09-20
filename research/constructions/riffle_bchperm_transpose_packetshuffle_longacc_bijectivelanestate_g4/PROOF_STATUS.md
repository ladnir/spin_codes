# Proof status

The construction has strong positive evidence at relative distance 0.09,
but it does not yet have a distance certificate.

## Inner reduction

Use two state classes.  Class `Z` is the zero 256-bit state.  Class `U` is a
uniform nonzero state.  A fresh nonzero field scalar refreshes class `U`
after every nonzero epoch output.

For a live epoch and fixed input, the emitted 256-bit word is an affine
translate of the full state space, with at most one state removed by the
nonzero conditioning.  If `0<z<1`, its output-weight moment is at most

\[
\frac{1}{1-2^{-256}}
\left(\frac{1+z}{2}\right)^{256}.
\]

For a zero-state epoch containing candidate dimension `d`, invertibility of
the four accumulators gives an information set of size `d`.  Its complete
output moment is at most

\[
\left(\frac{1+z}{2}\right)^d,
\]

and its exact-zero probability is `2^{-d}`.  These bounds retain the output
weight that earlier support-only recurrences discarded after a state reset.

## Numerical evidence

The current calculation uses `2^20` message bits, target distance
`floor(0.09*2^21)`, and the modeled regular spectrum of a `[256,128,38]`
outer code.  Arithmetic is nearest binary64, not outward-rounded interval
arithmetic.

For the packed rank profile `4q+r`, the pointwise minimum and summed regular
margin are both 55.951 bits.  The one-active-block row is dominant.  Some
reference rows are:

| Active outer blocks | Margin (bits) |
|---:|---:|
| 1 | 55.951 |
| 2 | 126.087 |
| 3 | 189.819 |
| 4 | 209.335 |
| 8 | 512.571 |
| 8192 | 119986.722 |

The candidate-dimension recurrence is exact for the packed placement model
under the stated state-coset envelope.  The receipt is
`receipts/profile_independent_k4_delta09.json`, field
`packed_rank_dimension_diagnostic`.

The cancellation-free live-state moment with two, three, and four state
dimensions per four-bit column gives at most about 617k, 972k, and 1,177k
bits of dense inner suppression, respectively.  The dense outer envelope is
about 1,057k bits.  Thus even the optimistic live-state calculation fails
for dimensions two and three; only the full four-dimensional state clears
the dense outer multiplicity at distance 0.09.

## Open lemma

An arbitrary active-block set produces a histogram
`(h1,h2,h3,h4)`, where `hr` counts packet groups of rank `r`.  Its total
candidate rank is

\[
a=h_1+2h_2+3h_3+4h_4.
\]

After the packet permutation, an epoch depends on a group only through its
rank.  Lane masks no longer appear because the live state spans all four
lanes.  The central target is therefore a rank-merging lemma: replacing two
partial ranks by a fuller rank, while preserving their total, must not
decrease the complete two-state moment after all 32 epochs and 256 outer
coordinates.  Repeated merging would reduce every histogram to `4q+r`.

This is a scalar rank/occupancy comparison.  It is materially simpler than
the earlier four-lane mask comparison, but it has not yet been proved.

After the rank-merging lemma, the remaining certificate work is to handle
the exceptional all-one outer word, replace or justify the modeled outer
spectrum, and use outward-rounded arithmetic.

## Bound that must not be used as evidence against the construction

The receipt's top-level `pointwise_minimum_margin_bits` comes from a
profile-independent reset union bound that discards every output bit in a
region containing a cancellation.  It fails by roughly one million bits.
That relaxation permits a cheap reset in every region and is intentionally
too coarse.  The packed candidate-dimension diagnostic retains the output
weight at each reset and is the relevant result.
