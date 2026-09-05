# Riffle LDPCSplitState g=4 t=256 s=64

**Status: paused packet branch as of 2026-08-29.** The packet route measured
about 1 ms faster than the comparable full-permutation route, but the fixed
four-block groups create a persistent shared-profile proof problem. That
gross speedup is not large enough to justify keeping this branch active. The
artifacts remain available if a substantially cleaner packet proof or a
larger measured speedup appears.

This folder explores a 64-bit randomized state with a paired sparse expansion
and compression.  The design objective is to make output cancellation and
state cancellation distinct events.

The interface now has a first fixed nested-code realization.  It uses a
64-bit systematic part, 128 layered-accumulator coordinates, and 64 sparse
parity coordinates.  `CONSTRUCTION.md` defines both the interface and this
fixed pair.  `PROOF_PLAN.md` gives the exact two-state reduction and the
coset-moment lemma used to handle adversarial epoch inputs.

The preferred fixed realization uses a degree-7 auxiliary map with a
degree-3 compressor.  It has no output-coordinate dependency of size at most
three and no compressor cancellation supported on at most three packets.
An exact degree-7 ensemble calculation separately proves that distance 40 is
available with ample all-live margin.

The current two-state occupation recurrence combines the modeled full outer
spectrum, exact packet placement, exact zero-state activation, and a
three-moment live bound.  Subject to distinct four-block packet groups and a
distance-40 constituent, it has positive 9% margin at every tested regular
occupation from 1 through 128.  Shared packet groups, the all-one outer word,
the denser occupations, and a fixed-instance distance certificate remain
open.

Shared packet groups have now been included exactly through occupation 16.
They remain positive and are easier than the distinct-group profile in that
range.  Packed and representative mixed profiles are also positive at
occupation 128.  A general packing lemma is still open because a split pair
can support a termination–reactivation path that disappears after merging.
Termination is not uniformly rare: the all-quad tilted moment concentrates
near 123 terminations.  The next proof target therefore compares complete
one-region matrices.  The maximally packed region matrix dominates every
profile in the exact low-occupation and near-packed scans.  All 270 recorded
elementary packing edges also have the required entrywise order.  This makes
a six-move region-level induction the current proof target.  The order is
not universal in the Chernoff parameters, but it is stable on a broad tested
range around the certificate's symmetric `p=1/2` point.
