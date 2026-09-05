# Proof status

Fix message length \(n=2^{20}\), output length \(N=2^{21}\), and target
distance \(d=\lfloor0.09N\rfloor\).  The calculations below use the modeled
complement-symmetric even spectrum for a binary \([256,128,38]\) outer code.
Probability is over the independent outer-coordinate permutations, the
independent packet permutations in the 256 regions, and the independent
nonzero field multipliers at checkpoint boundaries.

## One active outer block

Condition on an outer word of weight \(w\).  Its coordinate permutation
selects a uniform \(w\)-subset of the 256 regions.  In every selected region,
the packet permutation places its impulse in a uniform packet slot.  The
block's packet lane is fixed, so the impulse reaches only 16 of the 64 state
lanes.  This restriction does not change the transfer: for one impulse, the
epoch moment depends on its four possible visit numbers and not on its lane
label.  A checkpoint then maps every nonzero state to the uniform-nonzero
state class.

Consequently, the one-active packet transfer equals the bit-shuffle transfer
exactly.  Summing the modeled outer spectrum gives 81.5283 bits of aggregate
margin.  The dominant modeled outer weight is 38.

## Two active outer blocks of weight 38

Two blocks have one of three relationships.

| Relationship | Number of block pairs | Class margin |
|---|---:|---:|
| same group, different lanes | 12,288 | 178.3113 bits |
| different groups, same lane | 8,384,512 | 164.0327 bits |
| different groups, different lanes | 25,153,536 | 167.8930 bits |

The union of the three classes has 163.9366 bits of modeled margin.  The
uniform bit-shuffle baseline has 165.7100 bits, so packetization costs 1.7735
bits on this shell.

The adverse class consists of blocks in different groups but the same packet
lane.  In a region where both outer supports are active, distinct packet slots
can still address the same accumulator lane at different visits.  The two
impulses can then cancel before the checkpoint.  At the optimized tilt, the
dominant low-output term has 30 common active regions between the two
weight-38 supports.  Such a large intersection is rare; its hypergeometric
cost is included in the calculation.  Blocks in the same packet group occupy
different lanes, so they do not exhibit this cancellation mechanism.

The fixed-support epoch formula was checked against exhaustive occupancy
enumeration in a small instance to error below \(3.4\times10^{-16}\).  Every
packet-region transfer was stochastic at \(z=1\) to error below
\(6.7\times10^{-16}\).

## Open proof work

The subsequent regular-spectrum audit found a medium-density obstruction.
Activate all 2048 blocks in one fixed packet lane.  Under the fair-word
envelope, every epoch then receives candidates in only 16 of the 64 state
lanes.  The exact packed-residue calculation leaves approximately 132,900
bits of uncovered multiplicity at 9% distance.  The deficit is too large to
attribute to tilt resolution.

Thus, the favorable one- and two-active results do not extend to the complete
64-bit-state construction.  Work has moved to the distinct `g4 s256` variant.
Its checkpoint epoch visits every state lane once and removes the aligned-lane
mechanism.

## Artifacts

- `receipts/one_active_exact_transfer.json`
- `receipts/two_active_w38_packet_law.json`
- `receipts/packed_residue_full_regular_delta09.json`
- `../../scripts/analyze_riffle_packet4_fieldcheckpoint_twoactive.py`

The receipt values use optimized floating-point Chernoff tilts.  A final
certificate requires fixed tilts, outward-rounded arithmetic, and an explicit
outer-spectrum bound.
