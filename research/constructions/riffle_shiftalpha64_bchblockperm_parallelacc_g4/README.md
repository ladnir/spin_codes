# Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4

**Status:** `FINITE_DISTANCE_PROVED`

This candidate shifts the geometric double-parity coefficients from
\(\gamma^i\) to \(\gamma^{64+i}\). All BCH, permutation, packetization, and
accumulator mechanisms remain unchanged.

The shift removes every weight-22 pair in the complete extended-BCH minimum
shell. An exact scan tested 3,995,074,560 source/exponent pairs and found zero
intersections. Consequently, the shifted one-data-symbol outer spectrum has
no \((22,22,22)\) profile. This is a complete minimum-shell statement, not a
complete joint-spectrum theorem.

Goal 02 computes the complete target-weight histogram for those same
weight-22 sources. The target minimum is 30, attained 1,401 times. Hence the
minimum profile with repeated source weight 22 is \((22,22,30)\), of total
weight 74.

The complete weight-24 shell contains 6,855,968 messages. Its shifted target
minimum is 24, attained 1,175 times. Combining both exact shells with the BCH
minimum proves that the global one-data-symbol minimum is 72, attained by
\((24,24,24)\). This does not yet cover inputs with several nonzero data
symbols.

Goal 03 transfers the actual minimum profile through the accumulator at
threshold \(0.09N\). The conditional profile bound is \(2^{-40.2950}\).
Charging all 1,175 profiles gives \(2^{-30.0965}\), which passes a 20-bit
target by 10.0965 bits. Relative to the complete unshifted minimum shell, the
shift gains 14.0296 probability bits.

The arithmetic cost is one fixed multiplication by \(\gamma^{64}\) per
complete outer word. The encoder retains the existing per-symbol Horner
recurrence.

Goal 04 proves the deterministic end-to-end bound `d_min >= 33`. The outer
field code has block distance three, and every nonzero BCH block has weight at
least 22. Thus the accumulator input has weight at least 66. The identity
`u_t = s_(t-1) + s_t` gives `wt(u) <= 2 wt(s)`, so the output has weight at
least 33 for every setup.

The earlier probabilistic distance-26 calculation was valid but weaker than
this deterministic inequality. The first open threshold is output weight 33:
excluding it would prove `d_min >= 34`. The current parity-dropped first
moment there is 43.4018 and is therefore trivial. Retaining at least one
parity BCH block for multi-data messages is the next improvement.

Goal 05 replaces the inadmissible uniform-profile viewpoint with the actual
low-block occupation kernel. At threshold `0.09L`, the rigorous bounds are
`2^-40.2950` for `(24,24,24)` and `2^-44.6721` for
`(22,22,22,22)`. Real-law samples show nearly uniform 16-state behavior:
typical zero returns occur at rate `H/16`, and the prefix-state weight sum is
about `2H`. The rigorous low-tail sums instead tilt zero returns toward
`H/3`. The next proof target is therefore a 16-state spectral bound on the
joint excessive-return and favorable-gap event.

Goal 06 compresses the exact return count by a zero-return tilt and then
compresses the total input-weight coordinate by a coefficient tilt. The
resulting dynamic program retains only packet support and the 16-state
accumulator state. It loses 7.8908 bits on `(24,24,24)` and 5.4740 bits on
`(22,22,22,22)` relative to the exact Goal 05 sums, but scales to 64 active
minimum-weight blocks in a few seconds. The fixed-profile bounds strengthen
from `2^-39.1982` at four blocks to `2^-1073.4895` at 64 blocks. The next
proof task is to compare that contraction with a coarse outer
block-occupation spectrum.

Goal 07 begins that occupation climb. The exact MDS shell multiplicity
outgrows both the equal-weight-22 and equal-weight-64 kernel exponents through
128 occupied symbols. A new Hölder envelope retains both field parity
equations and compresses the exact BCH spectrum to one scalar moment. That
moment is dominated by BCH weights 22 and 128, so it remains trivial at the
first two shells. The next interface must retain a small BCH-weight band
count; occupation alone is insufficient, while exact weight profiles are
unnecessary.

Goal 08 implements light, normal, and heavy BCH-weight bands. The split gains
23.051 bits at occupation three and 41.021 bits at occupation four. Both
complete-shell bounds remain worse than the trivial message-count cap. Mixed
light-heavy patterns dominate. The next task is therefore to measure whether
the two outer field equations actually permit those mixed endpoint patterns,
not to add more weight bands.

Goal 09 performs that compatibility test exactly for the dominant
occupation-three endpoint. Across all 733,141,975,040 three-position
supports, the apparent `(22,106,128)` profile never occurs. The lightest
endpoint profile is `(32,96,128)`. Using the complete compatibility
histogram improves the corresponding inner bound by about 119 bits, from
`2^159.773` to `2^40.744`, but does not yet close it. The next calculation
should exploit that the weight-128 block contributes exactly 32 fixed
`1111` packets instead of treating those packets as arbitrary.

Goal 10 preserves those 32 fixed packets in the inner operator. The finite
endpoint shell improves by another 85.297 bits to `2^-44.5527`. An exact scan
and transfer of the smaller second-parity family gives `2^-55.4174`. Thus the
complete occupation-three family containing an all-one BCH block is bounded
by `2^-44.5519`. The remaining occupation-three calculation contains no
all-one block.

Goal 11 enforces the binary weight triangle on the remaining three-block
words. The apparent two-light/one-heavy obstruction disappears, and the
minimum profile `(22,22,22)` becomes dominant. The complete minimum BCH shell
has 227,584 unordered additive triples and 352,044 multiplier ratios. None of
those ratios occurs on a finite outer support. For supports containing the
second parity position, all data-coefficient differences are distinct. This
Sidon property and the full-state inner bound the complete minimum-profile
contribution by `2^-0.6561`. Both minimum-profile subcases are closed. The
open occupation-three remainder is `2^94.9201`, dominated by
`(22,106,106)` and `(22,22,24)` on finite supports.

Goal 12 closes both adjacent finite profiles. A shared affine-orbit
enumeration gives 1,365,504 ordered `(22,106,106)` relations and 27,765,248
ordered `(22,22,24)` relations. A source-preserving multivariate inner
transfer bounds the first family by `2^-37.0377` without an outer scan. An
exact shifted-schedule ratio scan finds 9,331,089 words in the second family,
whose contribution is at most `2^-6.4985`. Removing both lowers the open
occupation-three envelope to `2^92.2729`. The new leaders form a clustered
low-shell/complement family rather than one exceptional profile.

Goal 13 compresses the family through one complement model. Complement flags
reduce every low-shell/complement triple to an XOR of the low representatives
equal to either zero or the all-one word. A weight argument makes every
odd-complement triple over low weights `{22,24,26}` impossible. For the
remaining `(22,22,26)` profile, an exact shifted-schedule scan finds
61,929,389 outer words and closes the complete contribution at `2^-4.8039`.
The scan passes a direct 12-block audit and a regression against Goal 12.
Together with the first-pass closures of `(22,104,106)` and `(24,106,106)`,
this lowers the open occupation-three envelope to `2^90.0358`. The next
frontier is the even-complement weight-24 family.

Goal 14 counts that relation core through affine orbits. There are exactly
321,121,024 ordered `L22 x L24^2` relations and 4,018,336,896 ordered
`L24^3` relations. The high-weight two-complement realizations close from
these totals and the source-preserving inner operator. More strongly, an
exact shifted-schedule scan finds no finite `(22,24,24)` outer word. The
scan covers every support and every placement of the weight-22 block, and it
passes a direct 12-block audit. The open occupation-three envelope is now
`2^87.8525`, with `(24,24,24)` as the unique leading relation-core profile.

Goal 15 quotients `(24,24,24)` by affine triple orbits and the six
permutations of each additive triple. A symmetric field invariant reduces
669,722,816 unordered BCH relations to one schedule key each. Only 17 keys
hit the shifted schedule. The exact scan finds 60,462,011 outer words and
closes the profile at `2^-6.5547`. The remaining occupation-three envelope
is `2^87.2805`.

Goal 16 changes the proof interface. Let `h_k` count input packets of Hamming
weight `k`, for `k=1,2,3,4`. Conditional on these four counts, the parallel
accumulator has an exact five-state transfer formula. The states are the
possible Hamming weights of its four-bit state. The formula reduces to the
classical scalar accumulator formula at width one. Exhaustive enumeration of
all length-five inputs, comparison with the weighted 16-state operator, and
complete one-block packetization checks all pass. The BCH weight-profile case
ladder is paused while this four-count interface is integrated with the outer
sum.

Goal 17 begins a zero-, one-, and two-parity comparison. For zero parity, the
complete expected histogram count is the coefficient sequence of
`F(t)^16384-t0^524288`, where `F` is the exact one-block BCH histogram
polynomial. The minimum weight-22 shell contains 3,995,074,560 labeled outer
words and only 136 packet histograms. Every minimum-shell coefficient is
recorded exactly. Exact accumulator composition places the expected-count-one
crossover between output weights 76,000 and 77,000. At the 9% threshold, the
minimum shell alone misses the first-moment gate by 12.4306 bits. This is not
a complete zero-parity distance proof because higher shells remain unsummed.

Goal 18 replaces exact packet-type coefficient extraction by four positive
packet tilts and one output tilt. Each evaluation uses only a five-state
matrix power. On the exact Goal 17 shell, the numerical coefficient bound
loses about 12.8 bits near the shell threshold. The method therefore looks
appropriate for the bulk, but not for sparse shells with small margins. A
crude maximum over all packet types would cost another 71.4151 bits at the
target length, so the bulk sum must use adaptive summation rather than a
worst-type multiplier.

Goal 19 implements the complete small-instance proof gym. Exact rational
spectra are available through binary length 80, and an exact clipped low tail
is available at length 112. At output weight eight, one common coefficient
tilt loses about 11 bits. Optimizing each packet histogram separately recovers
about 4.4 bits at lengths 48 and 80. This stable recovery motivates adaptive
tilt regions between the sparse exact calculation and the bulk bound.

Goal 20 compresses the typewise tilts. Eight learned regions recover 87--89%
of the available improvement at lengths 48 and 80. Sixteen recover 98--99%.
The exact accumulator identity `H <= 2W` further removes every outer type
above weight `2D`. At `D=8`, only outer weights 12 and 16 remain. Treating
weight 12 exactly still leaves 2.4--3.2 bits of loss after all fifteen
weight-16 packet types receive separate tilts. The active obstruction is now
the conditional coefficient bound for one fixed type, not the number of
packet types.

Goal 26 attempts to determine the target saddle box and finds that the target
outer bulk gate is still missing. The Goal 25 calibration range does not
contain the natural target centers: at `D=76000`, the fourth log-fugacity is
about `-6.408`. A target face audit also exposes six support graphs that did
not occur in the weight-16 proof gym; all fifteen target face graphs are
primitive. The BCH minimum distance supplies a useful rare-coordinate gate.
Conditional on any boundary word of weight `H`, the expected number of
weight-four packets is at least `0.0009973753281 H`. At `H=152000`, the
numerical Chernoff optimum is `-47.19` bits. A dyadic rational witness
certifies `Pr[h4 <= 64] <= 2^-47` for each fixed outer word. This does not
sum the outer code, but it identifies the right split before a compact
Fourier certificate.

Goal 27 uses the current exact proof-gym enumerator only as a growth
measurement. With double parity, the expected-count-one crossing is nine for
`B=1,...,8` and ten for `B=9,...,15`. The binary output length grows from 24
to 136, so the relative crossing falls from 37.5% to 7.35%. The finite
`GF(16)` family therefore shows very slow absolute growth and no evidence of
linear distance. It ends at `B=15` and is not an asymptotic model of the
target construction.

Initial artifacts:

- `CONSTRUCTION.md`;
- `GOAL_01_SHIFTED_LOW_TAIL.md`;
- `receipts/goal01_shifted_weight22_exact.json`;
- `GOAL_02_WEIGHT22_TARGET_SPECTRUM.md`;
- `receipts/goal02_weight22_target_spectrum.json`;
- `receipts/goal02_weight24_affine_enumeration.json`;
- `receipts/goal02_weight24_target_spectrum.json`;
- `GOAL_03_MINIMUM_PROFILE_TRANSFER.md`;
- `receipts/goal03_minimum_profile_transfer.json`;
- `GOAL_04_END_TO_END_DISTANCE_FLOOR.md`;
- `receipts/goal04_end_to_end_d25.json`;
- `receipts/goal04_end_to_end_d8.json`;
- `receipts/goal04_end_to_end_d33.json`;
- `receipts/goal04_frontier_adjacent.json`;
- `GOAL_05_LOW_BLOCK_OCCUPATION_KERNEL.md`;
- `receipts/goal05_lowblock_occupation_summary.json`;
- `GOAL_06_COMPRESSED_RETURN_LADDER.md`;
- `receipts/goal06_compressed_return_ladder.json`;
- `GOAL_07_OUTER_OCCUPATION_LADDER.md`;
- `receipts/goal07_outer_occupation_ladder.json`;
- `GOAL_08_THREE_BAND_OUTER_GATE.md`;
- `receipts/goal08_three_band_outer_gate.json`;
- `GOAL_09_ENDPOINT_TRIPLE_COMPATIBILITY.md`;
- `receipts/goal09_endpoint_triple_scan.json`;
- `receipts/goal09_endpoint_triple_audit.json`;
- `receipts/goal09_endpoint_triple_transfer.json`;
- `GOAL_10_FIXED_ENDPOINT_CLOSURE.md`;
- `receipts/goal10_fixed_packet_transfer.json`;
- `receipts/goal10_p1_endpoint_scan.json`;
- `receipts/goal10_p1_endpoint_audit.json`;
- `receipts/goal10_p1_endpoint_transfer.json`;
- `GOAL_11_MINIMUM_TRIPLE_GATE.md`;
- `receipts/goal11_occ3_no_endpoint.json`;
- `receipts/goal11_weight22_additive_triples.json`;
- `receipts/goal11_weight22_ratio_counts.bin`;
- `receipts/goal11_weight22_outer_triples.json`;
- `receipts/goal11_weight22_outer_triples_audit.json`;
- `receipts/goal11_data_difference_multiplicity.json`;
- `receipts/goal11_data_difference_multiplicity_audit.json`;
- `receipts/goal11_weight22_inner_exact.json`;
- `GOAL_12_ADJACENT_PROFILE_CLOSURE.md`;
- `receipts/goal12_adjacent_triple_counts.json`;
- `receipts/goal12_adjacent_triples.json`;
- `receipts/goal12_adjacent_inner.json`;
- `receipts/goal12_adjacent_closure.json`;
- `receipts/goal12_weight222224_outer_scan.json`;
- `receipts/goal12_weight222224_outer_audit.json`;
- `receipts/goal12_occ3_after_adjacent.json`;
- `GOAL_13_LOW_SHELL_COMPLEMENT_CLUSTER.md`;
- `receipts/goal13_initial_relations.json`;
- `receipts/goal13_initial_inner.json`;
- `receipts/goal13_initial_closure.json`;
- `receipts/goal13_occ3_after_initial.json`;
- `receipts/goal13_weight22_22_26_ratios_raw.bin`;
- `receipts/goal13_weight22_22_26_ratio_counts.bin`;
- `receipts/goal13_weight222226_outer_scan.json`;
- `receipts/goal13_weight222226_outer_scan_12.json`;
- `receipts/goal13_weight222226_outer_audit.json`;
- `receipts/goal12_weight222224_outer_scan_12_regression.json`;
- `receipts/goal13_schedule_closure.json`;
- `receipts/goal13_occ3_after_schedule.json`;
- `GOAL_14_EVEN_COMPLEMENT_WEIGHT24_CORE.md`;
- `receipts/goal14_weight24_relation_core.json`;
- `receipts/goal14_even_core_inner.json`;
- `receipts/goal14_relation_closure.json`;
- `receipts/goal14_weight242422_outer_scan.json`;
- `receipts/goal14_weight242422_outer_audit.json`;
- `receipts/goal14_occ3_after_mixed_schedule.json`;
- `GOAL_15_EQUAL_SHELL_S3_QUOTIENT.md`;
- `receipts/goal15_weight24_s3_scan.json`;
- `receipts/goal15_equal_shell_closure.json`;
- `receipts/goal15_occ3_after_equal_shell.json`;
- `GOAL_16_PACKET_WEIGHT_ACCUMULATOR_FORMULA.md`;
- `GOAL_17_PARITY_LADDER_ZERO_COUNT.md`;
- `GOAL_18_PACKET_TYPE_BOUND.md`;
- `GOAL_19_SMALL_EXACT_VS_TYPE_BOUND.md`;
- `GOAL_20_ADAPTIVE_TILT_REGIONS.md`;
- `GOAL_21_FIXED_TYPE_LOSS_DECOMPOSITION.md`;
- `GOAL_22_BOUNDARY_SADDLE_LOCAL_FACTOR.md`;
- `GOAL_23_ONE_DIMENSIONAL_BOUNDARY_FACES.md`;
- `GOAL_24_TWO_DIMENSIONAL_BOUNDARY_FACES.md`;
- `GOAL_25_UNIFORM_LOCAL_BOUND_CANDIDATE.md`;
- `GOAL_26_TARGET_SADDLE_BOX_ATTEMPT.md`;
- `GOAL_27_SMALL_EXACT_GROWTH.md`;
- `receipts/goal17_zero_parity_accumulator_exact.json`;
- `receipts/goal18_packet_type_bound_probe.json`;
- `receipts/goal19_small_exact_vs_type_bound_b4.json`;
- `receipts/goal19_small_exact_vs_type_bound_b8.json`;
- `receipts/goal19_small_exact_vs_type_bound_b12.json`;
- `receipts/goal19_typewise_b4_p2_d8.json`;
- `receipts/goal19_typewise_b8_p2_d8.json`;
- `receipts/goal20_clustered_b4_p2_d8.json`;
- `receipts/goal20_clustered_b8_p2_d8.json`;
- `receipts/goal20_sparse_isolation_b4_p2_d8.json`;
- `receipts/goal20_sparse_isolation_b8_p2_d8.json`;
- `receipts/goal20_exact_shell_hybrid_b4_p2_d8.json`;
- `receipts/goal20_exact_shell_hybrid_b8_p2_d8.json`;
- `receipts/goal21_one_type_loss_b4_p2_h16_d8.json`;
- `receipts/goal21_one_type_loss_b8_p2_h16_d8.json`;
- `receipts/goal22_boundary_saddle_b4_p2_h16_d8.json`;
- `receipts/goal22_boundary_saddle_b8_p2_h16_d8.json`;
- `receipts/goal22_boundary_saddle_b12_p2_h16_d8.json`;
- `receipts/goal22_boundary_saddle_b15_p2_h16_d8.json`;
- `receipts/goal22_boundary_saddle_scaling_s4.json`;
- `receipts/goal23_boundary_one_face_formula_audit.json`;
- `receipts/goal24_boundary_even_face_formula_audit.json`;
- `receipts/goal24_boundary_odd_face_formula_audit.json`;
- `receipts/goal25_boundary_uniform_local_candidate.json`;
- `receipts/goal26_target_boundary_saddle_calibration.json`;
- `receipts/goal26_boundary_h4_gate.json`;
- `receipts/goal26_target_boundary_faces.json`;
- `receipts/goal27_small_exact_growth_b1_through_b15.json`;
- `receipts/goal16_packet_weight_formula_audit.json`;
- `ebch128_weight22_messages.bin` (complete \(A_{22}\) message set);
- `ebch128_weight24_messages.bin` (complete \(A_{24}\) message set);
- parent analysis:
  `../riffle_bchblockperm_parallelacc_g4/proof/GOAL_03_ALPHA_JOINT_SPECTRUM.md`.
