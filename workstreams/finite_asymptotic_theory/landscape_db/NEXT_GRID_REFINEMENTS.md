# Next work after the complete finite grid

The active goal is now the BCH-64/128 dominance audit described in
`BCH_DOMINANCE_ANALYSIS.md`. Q1 surfaces alone do not justify a full-margin
engineering recommendation. Update the existing geometries with sparse
occupation ratios, then control the complete remaining tail at targeted
extremes and around the apparent state-size knees. Keep unresolved points
explicit; do not infer dominance from Q2..4 alone. BCH-256 stays excluded.

Both t64/s20 references at K=2^20 now have verified full bounds: BCH-64
loses 0.08133 bits from higher occupations, and BCH-128 loses 0.000002
bits. The exact-kernel lower scan finds 46 replayed first-moment
obstructions: for both blocks, t64/s7..8, t128/s8..16, and t256/s9..20.
Next target t128/s17..20 and the t64 state-size knees, then vary K at
selected boundaries. The joined engineering grid still has 82 geometries
without a useful complete-tail conclusion from the updated analysis.

The engineering-surface comparison now extends the corrected Q1 treatment
to exact RM through length 512 and random references through length 512.
Read `CONSTITUENT_ENGINEERING_SURFACES.md` and its reproduction commands.
The comparison separates the K row-count cost, state-size knees, and outer
spectrum contribution. A long-region cancellation model supplies a common
explanation across families. The random ensemble curve and the conditional
spectrum-cap curve remain separate evidence classes. Next, test selected
higher occupations near the measured knees before treating the Q1 surface
as a model of the complete failure margin.

The priorities below record the earlier BCH-focused stage and full grid.

Current priority: targeted BCH growth analysis, following the user's
direction to understand continuous trends before choosing certificate
thresholds or implementation candidates. Read `BCH_GROWTH_ANALYSIS.md`.
The new four-state Q1 study explains the message-length slope and the
BCH-64/128 state-size plateaus. The earlier fit through all four BCH sizes
mixes distinct mechanisms and should not be used as a reliable forecast.
Next, extend the uniform-refresh invariant to selected higher occupations,
handling zero syndromes explicitly, and compare matched t/s settings.
BCH-256 remains a secondary bounded-spectrum calibration point.
The implementation and RM priorities below are retained as earlier plans.

All 3,108 native tuples now have evaluated full occupation coverage.
Read `GRID_FINDINGS.md` and `complete_landscape_audit.json` for the final
findings and authenticated snapshot. Full coverage includes weak bounds;
no new outward certificate is claimed. `complete_grid_unions.csv` is the
source for the selected signed margins and their evidence labels.

The primary scope remains exact BCH through length 128, exact RM through
length 512, and reused full-rank random references through length 1024.
BCH-256 and partial-spectrum RM(5,11) remain outside the primary grid.

The next work should concentrate on implementation candidates:

1. Replay the full RM(4,9), k=2^16, t64/s18 diagnostic with outward
   arithmetic. Its present full margin is 42.9694638927 bits.
2. Refine smaller-state and larger-t candidates selectively. The s=11
   Q1 screen passes 40 bits for all three t values, but that alone is
   insufficient. The composition comparison is recorded in
   `activation_occupation_grid_rm49_boxes_tradeoff_v1`; at t64/s11,
   Q248..256 remain loose. Use these witnesses to choose the next state
   or probability search. Do not infer that the code itself fails.
3. Benchmark selected t/s alternatives sequentially, including actual
   field arithmetic, permutations, outer encoding, and shared circuits.
   The independent-XOR and state-update proxies in
   `tested_parameter_frontiers.csv` are a shortlist, not runtime results.
4. For larger messages, use the exact-family projections to select a
   modest number of candidate constituents, then tighten their full
   occupation bounds. The Q1 holdout errors are about 14 bits for BCH
   and 8 bits for RM at the median; retain substantial forecast slack.
   Random references through 1024 are computed directly.

The deterministic typed dense grid is in
`activation_occupation_grid_typed_all_v1`. Its 39 batches cover Q65..L
for all 3,108 native tuples, with raw bounds and counting fallbacks.
The exact-region composition method is documented in
`COMPOSITION_BOX_METHOD.md`. Keep all recorded producer sources immutable;
version a changed method or witness search in a new output directory.
Numerical producers and benchmarks must run sequentially. Confirm that
the preceding process has terminated before starting another job.

To rebuild the checked snapshot after registering new results, run these
commands in order and stop on any failure:

```text
python build_landscape_db.py
python query_landscape.py export landscape_export.csv
python register_complete_grid.py --require-complete
python summarize_grid_unions.py
python parameter_cost_frontiers.py
python plot_tested_frontiers.py
python landscape_snapshot.py pack
python audit_complete_landscape.py
```

Run `python -m unittest discover -v *> grid_validation.log` in PowerShell
before the final audit. Raw SQLite/CSV files are local convenience
files; verified gzip snapshots also remain local and ignored. Use
`python landscape_snapshot.py restore` to restore missing raw files.

The sections below retain the development history and earlier probes.

The Q1 grid is complete: 3,108 native parameter tuples, with 303 dominant
witnesses at a grid boundary. Every tuple now also has an adaptive Q2..16
screen, giving 46,620 new occupation observations. Boundary witnesses and
negative margins remain visible. They do not establish code failure.

The six-map composition pilot adds 90 stronger Q2/Q3/Q4 rows. At k=2^16,
exact RM(4,9), t=64 and s=12, the composition margins are 66.547, 95.594,
and 99.965 bits. The existing exact-support Q2 margin is stronger, 83.693
bits, and should remain selected in unions. Q1 is 41.348 bits. Thus the
Q1..Q4 diagnostic still has about 41.348 bits; Q5 and higher must not be
inferred from that result.

The composition producer has completed all 3,108 tuples for Q2/Q3/Q4,
giving 9,324 observations in `activation_occupation_grid_composition_all_v1`.
All 39 map batches are registered and indexed. Its directory-wide
verification/resume command is:

```text
python run_sparse_composition_grid.py --output-dir activation_occupation_grid_composition_all_v1 --occupations 2 3 4
```

The dense pilot covers Q17..L on ten tuples. Those full-range bounds are
too loose to guide parameter choice. Most are little more than a counting
fallback. `complete_grid_unions.csv` distinguishes incomplete coverage,
evaluated but unhelpful bounds, and positive diagnostic unions. None is an
outward certificate. Do not expand the one-band dense pilot merely to
replace missing coverage with uninformative counting bounds.

## Stronger random setup events

`random_spectrum_variance.py` has an exact-integer construction of tighter
simultaneous caps. For a uniform D-dimensional subspace of F_2^B, put
T=2^B-1 and M=2^D-1. A shell has n=choose(B,w) nonzero vectors. Every
distinct pair of nonzero binary vectors is linearly independent, so

    E[A_w] = n M/T,
    Var(A_w) = n M (T-M) (T-n) / (T^2 (T-1)).

The pair-inclusion probability used here is M(M-1)/(T(T-1)); no independence
assumption is made. The implementation takes the minimum of deterministic
caps, the prior Markov cap, and an integer-rounded Chebyshev cap. Each shell
receives failure budget 2^-60/B. Choosing the tighter deterministic threshold
does not add two failure charges: the selected threshold itself has a bound
at most that shell's budget.

Exact tests enumerate all 35 binary [4,2] subspaces and recover both moments.
Other tests verify the simultaneous failure charge by rational arithmetic,
check domination of the old caps, and check that the length-1024 central
shell cap is within one part per million of its mean. The balanced
occupation producer now uses this helper and records accurately named
setup-event payloads. Old receipts remain unchanged.

The new good event implies the old Markov good event whenever its caps are
entrywise smaller, for the same outer parameters. Reusing old conditional
bounds under the stronger event is valid only after recording and checking
that implication. Both the coverage ledger and the union reporter now
authenticate these events and check their containment. Incomparable events
cannot be spliced without a separate probability argument. The failure
budget is still paid once in the final union.

## Balanced probability refinement

`BALANCED_OCCUPATION_METHOD.md` describes the new fixed-shell witnesses and
the tested affine-envelope acceleration. Scaling each shell's logit toward
one half is substantially stronger than applying the same additive logit
shift to every shell. The latter also moves central shells away from one
half and can make their counting cost dominate.

`activation_occupation_grid_balanced_pilot_v1` contains 1,530 rows: every
Q2..256 for exact RM(4,9) and conditional random [512,256], at k=2^16,
t=64, and s=12,18,20. At s=12 the exact RM margins over Q5..16 are at
least 60.851 bits; over Q17..32 they are at least 48.213 bits. Larger
occupations still have weak bounds. A finer tilt/scale sweep is recorded
separately in `activation_occupation_grid_balanced_refined_v1`.

The refined sweep is complete and registered: another 1,530 observations.
Combining those rows with the earlier sparse bounds gives full-range
random diagnostics of approximately 60 bits for all three states. The
shared setup-event budget determines the displayed margin. The exact RM
cases remain open. The current reports have 12 tuples with full evaluated
coverage, 3 positive full diagnostics, and 204 Q1..Q4 unions at least 40
bits. No current outward certificate is claimed.

The batched version prepares each spectrum's affine hulls once, evaluates
all starting matrices in one matrix-power batch, and reuses the same
power-of-two region ladder across message lengths. Tests compare its
matrices, margins, selected witnesses, and producer CSV values with the
original implementation. The full Q2..64 sweep has completed all 39 map
batches, giving 195,804 observations. Its verification/resume command is:

```text
python run_balanced_occupation_grid_v2.py --output-dir activation_occupation_grid_balanced_all_v2 --maximum-occupation 64
```

Its log is `activation_occupation_balanced_all.log`. Inspect the actual
process before resuming; do not start a second numerical job. Completed
map receipts are immutable. The producer/core files are now dependencies
of those receipts and must be versioned if changed. Register completed
map batches before rebuilding and exporting the database; wait for each
dependent command to finish.

The typed-box prototype is implemented and its cover, gradient, and
coefficient inequality have tests. Its current band partition and box
budget leave negative point bounds as well as box slack. It has not been
expanded into a production grid. Improving only the box count is not a
sufficient remedy for those negative point bounds.

A second typed probe uses the historical RM bands 1..95, 96..416, and
417..512, including the all-one word in the high band. This lowers the
low-band density cost and improves pure-low-type points substantially.
For t64/s18, its Q65..256 cover with 4,095 evaluated boxes still gives
about -2,190 bits. The worst boxes have almost fixed mixed type counts,
such as (188 zero, 49 low, 17 central, 2 high). The remaining loss therefore
includes coefficient conditioning at mixed types. More box subdivision
alone will not close it. Next compare an exact region mixture for three
band compositions with the typed Cauchy bound in this middle range; retain
the coefficient-free approach for larger L where full coefficient tables
are impractical.

## Dense bounds that retain row types

The next dense route should retain a small number of outer-weight bands
instead of dominating the entire nonzero spectrum with one Bernoulli law.
The following coefficient argument extends `DENSE_RANGE_METHOD.md` and
avoids a full L-coefficient table.

Use categories g=0,...,G, where category zero is a zero row, with p_0=0
and Gamma_0=1. Every nonzero category has a band density bound
nu_g <= Gamma_g mu_(p_g). An all-one singleton has p_g=1 and Gamma_g=1.
Let c_g be the category counts, with sum c_g=L. A positive proposal vector
a_g summing to one gives theta=sum a_g p_g. The positive multivariate
coefficient bound and the independent region permutations give

    U_c <= z^(-H) e_Z M(theta)^(BL/t) 1
           prod_g Gamma_g^(c_g) a_g^(-B c_g)
           multinomial(L;c_0,...,c_G)^(1-B).

The category counts stay fixed across all regions. Only the auxiliary
coefficient evaluation uses independent types. In log form, for a fixed
witness (z,p,a), this is affine in c plus
(B-1) sum_g log Gamma(c_g+1), up to constants. It is therefore convex on
the real count simplex. This suggests a bound over integer boxes intersected
with that simplex: evaluate every feasible vertex and multiply the maximum
by an upper bound on the number of integer count vectors in the box.

A box with integer lower/upper limits and sum c_g=L has integral vertices.
They can be enumerated by fixing all but one coordinate at a lower or upper
limit, then solving for the remaining coordinate. A counting upper bound is
the smallest product of coordinate widths after omitting one coordinate.
Splitting a coordinate at an integer threshold produces disjoint children,
so an adaptive tree can preserve complete coverage while refining loose
boxes. Restricting c_0<=L-17 excludes the sparse range already handled.

Before production, implement and test the typed coefficient bound, vertex
enumeration, integer-box counts, split coverage, and witness selection on
small exhaustive examples. Retain a fixed witness within each box: taking
different optimal witnesses at its vertices does not prove an interior
bound. Also check whether coefficient-conditioning overhead is too large
at small L; those cases may need exact region coefficients instead.

For exact RM(4,9), a first partition to investigate is low weights below
B/4, central weights through 3B/4, high weights below B, and the all-one
singleton. This is a proposed diagnostic partition, not a proved optimum.
The new variance caps should reduce random-outer density slack before the
same experiment is tried for random constituents.

## Index and validation

All 48 tests passed at this checkpoint. This includes complete native Q1
coverage, the new range importer, exact finite transfer checks, sparse
composition counting, random shell moments, and setup-event union handling.

After a producer finishes, register its directory, rebuild the SQLite index,
and refresh the coverage and union reports. The source manifests remain
authoritative. `summarize_grid_unions.py` preserves gaps, charges overlaps
conservatively, and pays one explicit setup-event charge. Its tests reject
gapped coverage and splicing incompatible setup events.

The remaining goal includes useful middle/dense evaluation across the grid,
parameter/cost comparisons, and clearly labelled scaling estimates. The
old extrapolation artifact remains historical and pending activation review;
it must not be promoted to a current fit. No new full-distance certificate
or implementation benchmark has been produced.
