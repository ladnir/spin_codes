# Cheaper XOR-column outer with the paired route

The matched K20 implementation takes 3.46 ms, versus 3.76 ms for the
certified paired-route control: 8.13% less full-encoding time. The candidate
has not earned a whole-code certificate. Its first correlated-input gate
now has a directed 50.89-bit class bound. Fresh 256/384-bit replays agree
on every cell, with independent 512-bit accounting of the complete union.

The retained [3.80-ms checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md) remains
the latest fully certified candidate. No production profile or paper changes.

## What changes

Keep K=2^20, rate 1/2, 8192 outer groups, 256 routed regions, and the h4
inner with width 64 and state 16. Each region independently permutes adjacent
pairs and independently swaps each pair. Matching parity alternates between
regions. The setup is fixed across encoding calls.

Replace four RS[64,32] constituents and a random invertible Toeplitz
column mixer by eight RS[32,16] constituents and D2(I+11^T)D1 over GF256.
The diagonal entries are independent uniform nonzero field elements.
The existing gauge normalization computes the same sampled outer map.
The group still encodes 128 symbols into 256 symbols, using one global route.

This removes about 19.7% of outer GFNI products and 983040 mixer XORs in
the static schedule. Hot outer keys shrink from 4.5 to 3.75 MiB. Neither
the operation counts nor the partial proof imply a full-code speed winner.

## The bounded proof gate

Fix exactly 17 active groups, each with exactly two nonzero constituent
words. The ideal-setup failure event is that some such scalar message has
output weight at most 734003. Other active-group counts and mixed
constituent counts are outside this gate.

Two geometries must remain fixed across all 32 RS columns:

- The active group support G. Let a,b count its complete pairs in the
  even and odd matchings. Exact cyclic counting gives 81 possible ordered
  pairs (a,b). Class (8,8) has 8192 supports, whereas (0,0) has about 96.7%
  of all 17-group supports.
- Each group's unordered pair of active constituent rows. Let mu be the
  maximum multiplicity of one of the 28 row pairs. The cumulative count is
  17! [x^17](sum_(j=0)^mu x^j/j!)^28; consecutive differences count mu exactly.

For a regional output subset U of G, with m occupied groups, its full-pair
count f obeys max(0,m-17+F)<=f<=min(F,floor(m/2)), where F is a or b.
Maximize the paired-route transfer matrix only over this interval.
The eight regions belonging to one RS column then use the ordered
tables (a,b,a,b,a,b,a,b), without reversing the inner chronology.

The outer's two-source column law has a 1/255 cancellation branch. The
geometry bound averages that branch before taking the state-potential
quotient. It retains the pairwise occupancy restriction controlled by mu.
The positive support polynomial is raised to power 32, with total source
support at least 34*17, before charging either geometry count. There is
one projective factor 255, not one factor per group.

Five disjoint mu ranges per (a,b) give 405 cells. Their weighted counts sum
exactly to binom(8192,17)*28^17. Each fixed cell may use its own literal
Chernoff mark and positive state potential; neither depends on setup draws.
The replay recalculates every count and bound rather than trusting saved
FLOAT estimates.

## What the iteration learned

A common potential initially left the full gate above one. Refining the
source counts alone was insufficient. The most concentrated row-pair
stratum needs different state potentials for different ordered matchings:

| Source(a,b), mu | z | log2 state potential | FLOAT log2 contribution |
|---|---:|---:|---:|
| (8,0),17 | .996598217 | -67.5 | -71.19867 |
| (0,8),17 | .996598217 | -61.5 | -70.19323 |

The complete saved screen gives -50.891481584647515. Selection used finite
literal banks, not an unbounded optimizer. Fresh 256/384-bit replays give
the same aggregate ceiling -3180717599041/62500000000, or
-50.891481584656. All 405 individual cell ceilings also agree. The largest
cell is source (4,0), mu=17, with ceiling -51.614684956866.
Independently summing the rounded cell bounds at 512-bit precision gives
the slightly looser ceiling -50.891481584655. This certifies the stated
homogeneous class only; it does not cover mixed constituent counts or
other active-group counts.

An exact common-row-pair specialization was also tested. It improves the
skew endpoint bounds by at most 0.033 bit and does not remove the need for
source-dependent potentials. Do not repeat a full sweep of that refinement.

The frozen proposal is
`theory/hd_xor_paired_joint_gate17_potential_screen.json`, SHA256
`ea67da7063d20164394f9a8aaea323bc26358adca16aea2d40770532bee699e7`.
The FLOAT proposal can be reconstructed with these bounded commands from
the theory directory. Output files must not already exist:

    python hd_xor_paired_joint_gate.py --output hd_xor_paired_joint_gate17_fixed_schedule_screen.json
    python hd_xor_paired_joint_select.py --seed hd_xor_paired_joint_gate17_fixed_schedule_screen.json --output hd_xor_paired_joint_gate17_selected_screen.json
    python hd_xor_paired_joint_select.py --seed hd_xor_paired_joint_gate17_selected_screen.json --marks 0.996598217 --potentials -67.5 -65.5 -63.5 -61.5 -59.5 --output hd_xor_paired_joint_gate17_potential_screen.json

Elapsed-time metadata varies, so regenerated proposal hashes need not match.
The literal selected witnesses are the inputs to fresh directed replay:

    python hd_xor_paired_joint_gate.py --precision 256 --witnesses hd_xor_paired_joint_gate17_potential_screen.json --output hd_xor_paired_joint_gate17_verified256.json
    python hd_xor_paired_joint_gate.py --precision 384 --witnesses hd_xor_paired_joint_gate17_potential_screen.json --output hd_xor_paired_joint_gate17_verified384.json

The helper refuses to overwrite existing receipts. Compare every cell
ceiling, not only the final union. The focused integration suite passes 41
tests, including ordered geometry, directed upper envelopes, exact counts,
literal witness loading, and rejection of incomplete/overlapping partitions.

Replay SHA256 identifiers:

- `hd_xor_paired_joint_gate.py`: `ea1916c445b3dadc3ca9a58d1739be268307d96be44571fee3e4df34a6bc1958`.
- `hd_xor_paired_joint_gate17_verified256.json`: `f1a30c19232b117b116545f66ef5cf7885e0203b7f6285684334fa501079ff35`.
- `hd_xor_paired_joint_gate17_verified384.json`: `d6ecffaf84c570c263fdadb7e20b47ad6349085f7b6faca41d17effe1073e94c`.

The source-support and constituent-pair counts sum exactly to
binom(8192,17)*28^17. Saved FLOAT counts and estimates are not inputs to
the directed arithmetic. Hashes identify these receipts; they do not
replace their mathematical scope or fresh arithmetic replay.

## Matched native measurements

All measurements are serial on Peach CPU 2, with 501 complete encodings per
run, precomputed setup, 128-bit records, and streaming output. The order was
control/candidate, candidate/control, then control/candidate.

| Build | Three run medians, ms | Median of medians, ms |
|---|---|---:|
| Toeplitz4 paired-route control | 3.762143,3.813519,3.761642 | 3.762143 |
| XOR2 paired-route candidate | 3.456162,3.484674,3.455270 | 3.456162 |

Separate phase runs give 1.864741 versus 1.628049 ms for the outer, and
1.723748 versus 1.720221 ms for the routed inner. These independently warmed
phase timings are diagnostics, not an additive decomposition of full time.
They locate the improvement in the outer and leave the inner as the next
large target. The older 3.30--3.33-ms unswapped timings remain separate.

The candidate passes 30720 packet basis directions, every pair-swap mask,
small/natural-size oracle comparisons, and an independent K20 oracle.
Its full-output checksum is 2343136256414281289. Executable SHA256:
`55bacc5d4999ffa7e00f0de70297dac626a3c4161f7347bb00527c7fd432a1e6`.
The control hash remains
`103c38e52809ab4a1030150280e55080dc1fdba844b4ee17a7ca95e791ed7312`.

The benchmark used the frozen paired-route native source snapshot, with
`HDGROUP_TOEPLITZ=1` replaced by `HDGROUP_XOR_MIXER=2`. All paired-route,
swap, h4, byte-key, and native-output flags remain enabled. The full command
uses `--bench 501 full 20 nt 8 4`; the control uses `4 4` for its last two
arguments. Source hashes in that snapshot are:

- `HdGroupedOuter.h`: `67f7ed49d635a99fc0502ccfd842d7c133f80db47575128fc665985567c81fcd`.
- `hd_grouped_forward.cpp`: `54d87c70f05db222c7abc49fd2227d92cb8e209688563c131297433ac2dced8a`.

## Next decisions

The first neighboring mixed-class gates now have useful results; the next
step is broader cohort coverage rather than another isolated corner fit.
For denser groups, the active pair may change between columns. Reusing the
fixed-pair multiplicity count would therefore be invalid. A mixed bound
must retain fixed row-set containment or use a proved density comparison.

The first isolated mixed-cohort helper, `theory/hd_xor_mixed_gate.py`, has
six exact/toy regression tests. Its four-corner screens inherit the saved
homogeneous witnesses without fitting new ones. Replacing one two-row
group by a one-row group gives a -7.20-bit bound, dominated by the (8,0)
source class at mu 9..14. Replacing it by a group with 3..8 active rows gives
a +683.02-bit bound under the initial dense comparison. Neither result
closes the mixed class, and neither is a low-distance construction witness.
These are historical initial screens, superseded for the specific classes
below by better support accounting and class-dependent witnesses.

A separate [one-group gate](HD_XOR_FIXED_MASK_FLOOR.md) now closes all
219 fixed active-row sets of sizes 3 through 8, at 164.796 bits. It replaces
the dense product comparison by the exact column-support mixture, holds
the row set fixed across all 32 columns, and retains at least 17 nonempty
columns. Fresh 256/384-bit replays agree on every contribution. This does
not close the neighboring mixed class: it supplies the sharper dense-group
law to use there. The improvement changes the proof, not the encoder or
its measured time.

The [mixed exact-mask bound](HD_XOR_MIXED_MASK_FLOOR.md) now closes one
dense group plus sixteen two-row groups for source matching counts (8,0).
It covers all 219 dense row sets and every constituent-pair assignment at
61.000916136179 directed bits. For each
fixed dense row set A and its output mask T in one column, add the
indicator of T to each region's sparse occupancy before evaluating the
paired-route cap. Average the exact mask law before the state-potential
quotient, using one potential for both empty and nonempty dense columns.
The resulting positive polynomial has one variable for sparse source
entries and another for nonempty dense columns. After 32 columns, retain
degrees at least 544 and 17 respectively. Sum A only after this extraction.
Source counting uses the 17-group support count and its 17 choices of
the dense group; it does not also charge 8192 choices. Five fixed witness
classes suffice. Full 256/384-bit replays agree, and independent 512-bit
union accounting includes all `28^16` sparse pair assignments. Other source
matching classes and multiple dense groups remain open.

The [one-single-row/sixteen-two-row cohort](HD_XOR_SPARSE_MIXED.md) has a
70.3541-bit floating-point proposal covering all 81 source classes. Its
fresh directed replay is in progress. The [coverage plan](HD_XOR_COVERAGE_PLAN.md)
also supplies a terminal argument for all cohort mixtures with at least
61,191 active constituent rows; this argument gives 165.78 directed bits.
The intermediate-density classes still need a joint bound.

In parallel, the [rank-one refresh study](HD_RANKONE_STATE_MASKS.md) and
[diagonal-XOR-diagonal native study](HD_XOR_STATE_NATIVE.md) preserve their
negative or modest-gain results. They do not justify changing the retained
inner while the cheaper outer's coverage remains incomplete.
Do not promote XOR2 or replace the retained certificate until all message
classes are covered.
