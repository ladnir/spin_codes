# Two active permutation proof tracks

The independent-row experiment is an additional route, not a replacement
for the existing proof effort. Both retain BCH[256,128], K = 2^20,
N = 2^21, IMT(128,19), two updates, and four-row packets. The target remains
10% relative distance with at least 40 bits of setup-failure margin over
all nonzero messages. Neither route currently has that full certificate.

| Track | Coordinate shuffle | Measured transpose | Proof status |
|---|---|---:|---|
| Shared shuffle | One permutation per four-row group | 5.42 ms | All occupancies 1–14 covered; unrestricted 15–2048 open |
| Independent rows | Separate permutation for every BCH row | 6.50 ms | All occupancies 1–32 covered, combined margin above 43.74 bits; unrestricted 33–2048 open |

Times are the recent precomputed 128-bit-element measurements, not setup
costs. See [PACKET_SPLIT.md](PACKET_SPLIT.md).

## Shared-shuffle track: retained and active

The entry point remains `occupancy_cdf_cover.py`. Its default ensemble,
switches, numerical witnesses, and replay behavior are unchanged. The
two-update results and reproduction commands remain in
[TWO_UPDATES.md](TWO_UPDATES.md). The combined bound for messages supported
on at most fourteen groups exceeds 41.83 bits. Restricted dense results
do not cover their unrestricted complements.

The remaining challenge is the intermediate/dense contribution, including
joint pattern multiplicities and returns to zero. The latest joint
cancellation work remains in [JOINT_CANCELLATION.md](JOINT_CANCELLATION.md).
Further refinements of those bounds should continue to be tested against
this route. Its lower implementation cost remains an independent reason
to pursue it even if the new route becomes easier to prove.

## Independent-row track: isolated exploration

New proof code lives under `independent_rows/`. It reads authenticated BCH
shell bounds and, for the diagnostic, existing universal per-shape inner
operators. It does not edit or replace the shared-shuffle verifier.
The new route averages independent coordinate permutations per row;
old joint-rank/support counts and old certificates must not be transferred.

See [independent_rows/README.md](independent_rows/README.md) for the exact
counting reduction and tests, and [STRONG_INNER.md](independent_rows/STRONG_INNER.md)
for the eleven-coordinate outward results. The later
[complete sparse batch](independent_rows/LOW_OCCUPANCIES.md) independently
covers every occupancy from 1 through 32; their combined contribution is
below 6.786358e-14. No interpolation between occupancies is used.
Proof status must be recorded under the correct ensemble, even when a local
inner lemma is shared. [BRIDGE.md](independent_rows/BRIDGE.md) states the
complete reduction from conditional inner bounds to the averaged outer measure.

At occupancy 64, selected homogeneous support events remain open; retuning
and total-input-weight tilts did not close the middle supports. The next
analysis uses packet-shape averages under independent row conditioning.
These diagnostics do not replace mixed-support covers or missing occupancies.

A separate [dense-component bound](independent_rows/DENSE_UNIFORM.md) controls
part of a positive decomposition of the averaged measure. Its component
index is not active-group occupancy, so it is not added to the coverage
range in the table.

This exploration left `occupancy_cdf_cover.py`, `occupancy_memory.py`,
`TWO_UPDATES.md`, and `JOINT_CANCELLATION.md` byte-for-byte unchanged.
No historical proof output, source, or benchmark data was removed.
