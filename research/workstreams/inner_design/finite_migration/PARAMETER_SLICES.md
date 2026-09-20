# IMT parameter-slice migration

## Current paper scope: Q1 slices with selected full certificates

The user refined the objective on 2026-09-16: retain broad Q1 sensitivity
curves and full-certified operating points, rather than pursuing a full
certificate at every diagnostic geometry. The paper now uses the no-constant
grid `PARAMETER_NO_CONSTANT_Q1_v1.json` and its `_VERIFIED_v1.json` receipt.
All 130 geometries and 39 maps are covered. The grid is a binary64 diagnostic,
not an outward-rounded full certificate. The five selected BCH-256 lengths
provide an exact-map-matched Q1/full comparison separately: loss is 1.774507
bits at K16 and below 0.0002 bits at K18, K20, K22, and K24.

Those BCH-256 maps are not the diagnostic BCH-64/128 maps. The selected
certificates therefore do not certify endpoints of the smaller-outer slices.
Do not resume the full diagnostic-grid search merely to finish the paper.
A future map-matched diagnostic certificate would strengthen calibration,
but is not part of this delivery's stopping condition.

Reproduce with `python -B paper/build_imt_parameter_figures.py`; append
`--check` to authenticate inputs and compare all four generated figures.
The new producer preserves the original grid producers and historical plots.

## Historical exploration record

The following records preserve the preceding full-grid objective and searches;
their open full-cover status is not an open requirement for the revised plots.

The replacement calculation retains all 130 geometries from the paper's
three parameter slices. Its one-active-row calculation is complete. The
higher-occupancy calculation remains open, so these results do not yet
replace the paper's full-margin curves.
The original diagnostic chain below includes the all-one expansion word.
A replacement chain excluding that word is now being tested; its evidence
is kept separate.

## Fixed maps and model

For each t in {64,128,256}, the expansion A uses the historical study's
ordered RM(2,log2(t)) basis. Its seed is 3390173185+t. Taking its first s
generators reproduces that study's expansion map, including its constant word.
The parameter study therefore differs from the selected BCH-256 implementation's
expansion map, which excludes the constant word.

Feedback B uses a separate ordered chain, generated with seed 2026091601+t.
At the minimum state size s0=log2(t)+1, its t columns are distinct nonzero
s0-bit vectors. Each added state bit extends these columns with a sampled
row, retaining the extension only when it increases rank. The resulting
prefix maps have full rank and distinct nonzero columns at every tested s.
The receipt records every column; no map-quality search selects these chains.

The same fixed A and B maps are used for both BCH outers and every K at
a given (t,s). These dense feedback chains are diagnostic configurations,
not the optimized sparse feedback maps used for the reported timings.

The inner uses one independently sampled transvection per epoch, with output
before update, zero initial state, persistent state across regions, and no
flush. `weight_memory.transfers` supplies the IMT Q1 envelope. Tests compare
its degree-zero and degree-one matrices with the separate-map general transfer.
No exact-refresh RM2Sub transition supplies a new numerical result.

## Coverage and results

The exact-spectrum outers are the shortened-XBCH [64,32,12] constituent and
extended BCH [128,64,22]. The cutoff is floor(N/10), with N=2K. Coverage is:

- B=64,128; (t,s)=(64,20); log2(K)=12,...,26.
- B=64,128; t=64,128,256; s=log2(t)+1,...,20; log2(K)=20.
- B=64,128; t=64; s=10,12,16; log2(K)=16,18,20,22,24.

Overlapping cells are evaluated once. `PARAMETER_Q1_v1.json` retains 130
Q1 rows and 39 exact map pairs. `PARAMETER_Q1_VERIFIED_v1.json` authenticates
the inputs, reconstructs all maps, and recomputes every row. Seven cells
also pass a separate log-domain coefficient calculation; their largest
margin difference is about 7.11e-15 bits.

| Configuration, log2(K)=20 | BCH-64 Q1 margin | BCH-128 Q1 margin |
|---|---:|---:|
| t=64, s=10 | 7.224798 | 26.989693 |
| t=64, s=16 | 10.184734 | 32.341189 |
| t=64, s=20 | 10.340377 | 32.631786 |
| t=128, s=20 | 10.291956 | 32.493523 |
| t=256, s=20 | 10.280831 | 32.462695 |

These numerical bounds concern only messages with one active outer row.
In particular, the similar Q1 values across t do not show that larger
epochs support an equally good full margin. The higher-occupancy contribution
is the missing comparison. All arithmetic here is nearest binary64, not
outward certification or a measurement of actual failure probability.

`parameter_full.py` rebuilds separate-map sparse transfers, exact feedback
kernel counts, and an IMT Fourier transfer for dense type boxes. Its first
bounded BCH-128/t64/s20/K20 pilot gives a sparse Q=2..4 bound of 23.37 bits,
but the 63-node dense search remains uninformative. The complete diagnostic
upper bound has a large negative margin; it is not evidence of poor distance.
`parameter_sparse_balanced.py` subsequently retained each Q=2,3,4 band
composition separately and scaled reference logits toward one half. Its
margins are **34.153247, 53.990145, and 72.967441 bits**, respectively.
`PARAMETER_SPARSE_BALANCED_VERIFIED_v1.json` replays all 990 composition
witnesses. These results are still binary64 diagnostics. They cover only
the three specified occupancies of the BCH-128/t64/s20/K20 pilot.

The preceding fixed-composition trial used shifted reference logits and
gave worse bounds. Its receipt is retained separately; composition retention
alone does not compensate for unsuitable reference probabilities.
The dense partition and the rest of the higher-occupancy grid remain open.

The pilot now also covers Q=5,...,64. `parameter_sparse_range.py` uses a
scaled positive matrix-polynomial product and the three balanced reference
banks. The aggregate bound for this range is **83.534893 bits**.
`PARAMETER_SPARSE_RANGE_VERIFIED_v1.json` checks every retained witness using
the separate log-domain product. Its largest log-bound discrepancy is
1.03e-12. Together with the preceding calculations, this covers Q=1,...,64
for BCH-128/t64/s20/K20, but not the full message set.

The first 1,023-node mixed dense search remains uninformative. Its selected
boxes pass exact integer-coverage and numerical witness checks, but their
union gives no useful margin. Extending the tilt grid alone also leaves a
large gap. Neither outcome refutes the distance claim: both calculations
upper-bound the failure probability, and a large upper bound is inconclusive.

## Causal moment bound for dense inputs

The generic state envelope can lose slack even for unbiased input. The
following bound recovers the exact moment in that case without changing IMT.
Fix the setup, and let the input bits be independent Bernoulli(theta).
Write N for their number and lambda>=0 for the exponential tilt.
For each epoch, condition on all preceding input epochs. The current state
is then fixed, and each output bit is a fresh input bit XOR a fixed bit.
Its conditional exponential moment is at most

```text
a(theta,lambda) = 1-p + p exp(-lambda),   p=min(theta,1-theta).
```

Conditional independence within the epoch bounds its moment by a^t.
Iterated conditioning therefore gives

```text
E[exp(-lambda wt(Y))] <= a(theta,lambda)^N.
```

At theta=1/2, the bound is an equality. This argument holds for each fixed
setup and requires output before update; it does not assume state refresh.
`parameter_causal_dense.py` takes the minimum of this complete moment bound
and the two existing IMT bounds. Exhaustive tests cover tiny causal encoders,
including nonlinear feedback, and check the equality at theta=1/2.

`parameter_three_type_dense.py` also groups every ordinary shell into one
category, leaving three types: zero, ordinary, and all-one. The ordinary
category is dominated by a Bernoulli(1/2) counting measure with density cost

```text
Gamma = max_{0<w<B} A_w 2^B / binom(B,w).
```

This coarsening increases the counting cost but reduces the partition
dimension. The all-one category remains separate. The search results must
still cover all occupancies and pass replay before they supply a full curve.

`type_box_coverage.py` counts the integer points in each selected box by
inclusion-exclusion. It checks disjointness and equality with the root-domain
count. The dense verifiers then rebuild every map and replay every selected
moment witness. Coverage is exact; the moment replay remains nearest binary64.

## Excluding the all-one expansion word

The three-type pilot isolates a problem in the original diagnostic chain.
Its worst selected types have roughly 90% all-one outer rows. Even singleton
types retain a large gap, so finer boxes cannot resolve this particular loss.
The expansion contains the all-one word, allowing a persistent state to
cancel all-one input epochs. This observation identifies a weakness in the
bound; it does not establish a code-distance counterexample.

The selected BCH-256 implementation already excludes the all-one expansion
word. `parameter_no_constant.py` tests the same restriction for the parameter
slices. It preserves the feedback chain and all geometric parameters.
Its expansion begins with the m linear functions and 20-m independent
quadratics, where t=2^m. A fixed shift adds constant coefficients to these
generators. Independence of the nonconstant parts excludes the all-one word.
The shift is chosen outside the coordinate images at the minimum state size,
so every prefix has distinct nonzero columns. These are fixed diagnostic maps,
not a change to the timed implementations.

For BCH-128/t64/s20/K20, the four previously worst singleton types now have
large positive bounds. The Q1 bound changes by less than 1e-6 bits, to
**32.631785 bits**. The new Q5..64 aggregate is **83.534892 bits**, with
independent log-domain replay. The Q2/Q3/Q4 bounds are 34.153247, 53.990145,
and 72.967440 bits; all 990 composition witnesses passed replay.
These local results motivate the new chain but do not
yet establish the full dense cover or validate all 130 configurations.
`PARAMETER_NO_CONSTANT_PROBE_v1.json` retains the exact maps and point witnesses.

The `parameter_no_constant_dense.py` and `parameter_no_constant_sparse.py`
wrappers select this model explicitly for both search and replay. They clear
the model-dependent preparation cache on entry and exit. A verifier using
the original maps must reject the new receipts rather than reinterpret them.

The first 1,023-node cover for the new chain remains inconclusive. Its
512 selected boxes pass exact coverage and numerical replay, but its union
margin is -42870.97 bits. The worst retained types now contain almost no
all-one rows. Several of those boxes contain only one integer type, so
further geometric subdivision alone cannot close the remaining gap.
The next calculation should improve the moment or counting bound at these
fixed types before attempting a larger grid. The no-constant restriction
removes the earlier extreme-type loss; it does not solve the full pilot.

Joint proposal/tilt tuning resolves the four worst ordinary-row singleton
examples from that cover. Their numerical margins become 16889.83, 17055.44,
5165.84, and 5218.05 bits. `PARAMETER_JOINT_POINT_v1.json` records these
point witnesses. The coarse tilt bank, followed by proposal-only tuning,
missed these choices. The point results do not establish a complete cover.

The first fast joint-tuned cover used only the Fourier and causal moments.
It lost useful fixed-weight bounds at smaller input densities and remained
inconclusive. `parameter_hybrid_dense.py` retains both families as complete
box bounds, recording which family supplies each selected witness. Its
verifier reconstructs that family; no entrywise minimum mixes representations.
For subsequent runs, `parameter_hybrid_safe.Dense` additionally rejects
optimizer proposals whose probability rounds to zero or one. It does not
clip those probabilities or treat an invalid proposal as a useful bound.

`PARAMETER_NO_CONSTANT_Q1_v1.json` now contains the full 130-cell grid for
the new expansion chain. `PARAMETER_NO_CONSTANT_Q1_VERIFIED_v1.json`
reconstructs all 39 map pairs and repeats every cell. Seven cells additionally
pass log-domain replay; the largest difference is 1.43e-14 bits.

The first hybrid cover finished with 512 selected boxes and a numerical
margin of -30144.78 bits. All witnesses and the exact integer partition
passed replay; the bound remains uninformative. Its worst boxes use the
fixed-weight family at about 16% active rows. Jointly retuning that family
turns the four worst singleton examples into margins of 1727.27, 1832.66,
1960.38, and 1727.38 bits, recorded in `PARAMETER_MIXED_POINT_v1.json`.
`parameter_joint_mixed_dense.py` now applies this additional optimization
only where the existing hybrid box bound is weak. Its full-cover result
requires a separate replay and aggregate check before use in a figure.

That 1,023-node search has now finished. Its 512 selected boxes passed
`parameter_joint_mixed_verify.py`, including the exact integer partition
and every selected whole-moment witness. The dense margin improves to
-19,925.44264 bits, still an uninformative bound. The worst type is the
singleton (zero, ordinary, all-one)=(12809,3575,0), about 21.8% active rows;
its selected tilt is -0.3511 and iid input probability about 0.1586.
Further subdivision of that singleton cannot help. Before another complete
cover, inspect the moment/likelihood tradeoff and try genuinely different
tilt/proposal starting points there. The completed selected BCH-256 and
quarter-rate certificates are unaffected by this diagnostic family's result.

Receipts are `PARAMETER_JOINT_MIXED_DENSE_v1.json` and
`PARAMETER_JOINT_MIXED_DENSE_VERIFIED_v1.json`; the aggregate is
`PARAMETER_FULL_JOINT_MIXED_v1.json`. All are nearest-binary64 evaluations,
not outward certificates or evidence that the actual code has poor distance.

`parameter_pilot_union.py` authenticates the Q1 grid, sparse compositions,
sparse range, and dense cover with their respective replay receipts. It
checks their common maps and contiguous Q coverage before summing the
numerical bounds. `PARAMETER_FULL_HYBRID_v1.json` retains the first complete
but uninformative aggregate. It is explicitly not an outward certificate.

## Reproduction

### Latest dense-point diagnostic

The initial worst retained type was (zero, ordinary, all-one)=(12809,3575,0)
for BCH-128, t=64, s=20, K=2^20, and target distance .10. These are the
diagnostic nested maps, not the selected timed BCH-256 maps.

`PARAMETER_FRONTIER_FACE_v1.json` records a multistart search over the
Chernoff tilt and the ordinary-row proposal probability. Five of six local
searches find a different basin: its point margin is about -1061.863 bits,
compared with -19925.443 in the accepted numerical cover. This is a point
improvement, not a replacement full cover. The all-one proposal probability
remains positive even though this type has no all-one rows.

`PARAMETER_SYNDROME_ACTIVATION_v1.json` evaluates every feedback syndrome
under the tilted Bernoulli input law. It uses the largest syndrome probability
within each expansion-weight shell to bound activation from state zero.
The best tested point improves to about -938.548 bits. Exact rational toy
enumerations check the Walsh inversion, but the production calculation uses
binary64 with a numerical guard, not certified outward arithmetic.

`PARAMETER_DIFFUSE_LAZY_v1.json` also preserves a pointwise density bound
after a lazy update from a uniform expansion-weight shell. It adds no gain
at the tested point. `test_diffuse_lazy.py` checks the nonzero-syndrome cap
with rational arithmetic and checks the implemented envelope against exhaustive
small examples, including rank-deficient feedback and zero feedback.
These tests do not validate the production floating-point calculation as a
certificate.

Those three initial diagnostics were uninformative as failure bounds. A negative
reported margin means that the evaluated union bound exceeds one; it does
not show poor distance for the code. The selected finite certificates are
unchanged.

### Retuned activation and the complete cover

Joint optimization of the sharper activation calculation resolves the first
singleton. `PARAMETER_ACTIVATION_REFINED_v1.json` records a point margin of
1981.270085 bits, with tilt -0.601272289 and ordinary-row proposal probability
0.288012527. This is a binary64 point bound, not a full distance certificate.

`parameter_activation_bank.py` reuses a moment across count boxes. For a
fixed input probability theta, the proposal is
`(1 - 2*theta + u, 2*(theta - u), u)`. Optimizing the all-one probability u
preserves theta, so the same moment remains applicable. Each proposal is
positive and held fixed over every vertex of its box. The method compares
whole bounds; it never takes entrywise minima across state representations.

The first bank improved 402 of the original 512 boxes. Adaptive subdivision
then exposed other weak singletons. `parameter_activation_frontier.py` tunes
new witnesses there and reuses them throughout the cover. Its first four
targets give the following point improvements:

| (zero, ordinary, all-one) counts | Previous margin | Retuned margin |
|---|---:|---:|
| (15359,1025,0) | -17255.83 | 2313.14 |
| (14832,1552,0) | -17254.79 | 2300.58 |
| (3071,13312,1) | -16598.73 | 5157.19 |
| (1535,14849,0) | -16434.50 | 3109.20 |

These gains establish that optimizer choices caused part of the loss.
They do not establish that optimization alone can close every remaining type.
The next broad box did not improve after retuning. Further subdivision with
the eight retained bank entries produced 1,536 boxes and a dense margin of
**-14901.086493 bits**. That complete bound is still uninformative.

`PARAMETER_ACTIVATION_PARTITION_v2.json` retains that cover;
`PARAMETER_ACTIVATION_PARTITION_VERIFIED_v2.json` reconstructs every selected
moment and checks the exact disjoint integer partition. The all-occupancy
aggregate is `PARAMETER_FULL_ACTIVATION_PARTITION_v2.json`. Replay uses the
same binary64 transfer formulas, not an independent outward calculation.
The new regression tests check constrained proposal optimization and exact
small-domain partitioning; all 68 finite-migration tests pass.

Its worst type was (16319,65,0), at the start of the dense range; another
weak singleton was (16128,256,0). That observation motivated the sparse
extension below. Singleton intersections cannot benefit from subdivision.

### Sparse prefix through Q=320

The same sparse envelope now covers Q=65..256 with an aggregate margin of
309.034529 bits. Its independent log-domain replay differs from the scaled
positive calculation by at most 3.64e-12 in logarithmic bounds. Extending
the calculation through Q=512 does not close the entire range: its aggregate
margin is -635.685800 bits. However, the Q=257..320 prefix has an aggregate
margin of 83.961541 bits. The full Q=257..512 replay passes with maximum
logarithmic error 1.10e-11, so this prefix uses verified individual terms,
not an extrapolation from Q=256.

Receipts are `PARAMETER_SPARSE_65_256_v1.json` and
`PARAMETER_SPARSE_257_512_v1.json`, with corresponding `_VERIFIED_v1.json`
records. `parameter_sparse_extension.py` supplies both search and replay
under the diagnostic no-constant maps.

`parameter_range_bridge.py` clips the old dense partition to Q>=321 and
recomputes its inherited witnesses. Exact integer coverage checks establish
that the sparse prefix and dense suffix have no gap or overlap in Q.
The combined sparse prefix Q=1..320 has margin **32.200608 bits**. This is
still a partial first moment, not a distance guarantee: Q>=321 remains.
`PARAMETER_RANGE_BRIDGE_320_v1.json` and its `_VERIFIED_v1.json` record
authenticate the assembled inputs and replay the clipped dense contribution.

The suffix-specific adapter `parameter_scoped_cover.py` accepts only an
authenticated, replayed input. It reconstructs the starting witnesses,
retunes weak boxes, and bisects a still-weak nonsingleton when useful.
Whole-moment witnesses remain fixed over each box. Its first four-step pass
resolves three further singletons but leaves a dense margin of
-13600.766098 bits over 1,433 boxes. `PARAMETER_SCOPED_321_v1.json` and
`PARAMETER_SCOPED_321_VERIFIED_v1.json` record that pass and its replay.
`PARAMETER_FULL_SCOPED_321_v1.json` joins it to the sparse prefix and confirms
that the full numerical bound is still uninformative.

Continue on the dense suffix, not on Q values already handled by the sparse
prefix. Stop repeating a local search when its improvement is negligible;
switch to another type or subdivide a nonsingleton. None of these binary64
diagnostics changes the selected outward certificates or current timings.

The subsequent twelve-step suffix pass improved the dense margin to
-12886.059788 bits. A 256-split pass retained 1,699 boxes and improved it
to -12510.052942 bits. Both passed reconstruction of their moments and exact
integer coverage checks. The split pass is retained as
`PARAMETER_SCOPED_PARTITION_321_v1.json`, with its matching verification.

An exactly fair Bernoulli input supplies an important fixed anchor. Conditional
on every previous input and the setup, each fresh input bit remains fair
after XOR with its fixed state-dependent bit. Thus the causal moment is
exactly `N * log((1 + exp(-lambda))/2)`, independently of the state maps.
The counting proposal must still have mean input probability one half.
`parameter_uniform_anchor.py` optimizes that proposal within each box and
adds this whole-moment witness at the binomial Chernoff tilt. Its reconstructed
IMT moment agrees with the causal expression; no transfer entries are mixed.

The anchor improved 151 boxes, including a previously weak singleton that
became strongly positive. The latest complete suffix is
`PARAMETER_UNIFORM_ANCHOR_321_v1.json`, replayed in
`PARAMETER_UNIFORM_ANCHOR_321_VERIFIED_v1.json`. Its dense margin is
**-12330.953193 bits**, still uninformative. The full aggregate in
`PARAMETER_FULL_UNIFORM_ANCHOR_321_v1.json` has the same displayed margin;
the sparse prefix remains 32.200608 bits.

The current worst singleton is (4095,12288,1). At its retained proposal and
tilt, the diffuse-lazy refinement gives exactly the same binary64 moment as
the activation calculation. This one-point test does not exclude gains at
other witnesses. Next, retain the fair-input anchor explicitly and test
other starting basins at the newly exposed singleton before more subdivision.
The latest run passes 75 finite-migration tests, 14 artifact tests, and six
paper evidence tests. Selected finite/asymptotic integration and the current
comparison table also pass. No encoder benchmark or commit was performed.

### Commands

Run from the repository root with the inner-design Python dependencies:

```text
python -B workstreams/inner_design/finite_migration/parameter_q1.py --output <fresh-q1.json>
python -B workstreams/inner_design/finite_migration/verify_parameter_q1.py --input <fresh-q1.json> --output <fresh-verification.json>
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p "test_parameter_imt.py" -v
python -B workstreams/inner_design/finite_migration/parameter_sparse_balanced.py --output <fresh-sparse.json>
python -B workstreams/inner_design/finite_migration/verify_parameter_sparse.py --input <fresh-sparse.json> --output <fresh-sparse-verification.json>
python -B workstreams/inner_design/finite_migration/parameter_sparse_range.py --output <fresh-range.json> --first 5 --last 64
python -B workstreams/inner_design/finite_migration/verify_parameter_range.py --input <fresh-range.json> --output <fresh-range-verification.json>
python -B workstreams/inner_design/finite_migration/parameter_three_type_dense.py --output <fresh-dense.json> --minimum 65 --nodes 1023
python -B workstreams/inner_design/finite_migration/verify_parameter_dense_v2.py --input <fresh-dense.json> --output <fresh-dense-verification.json>
```

The producer refuses to overwrite an existing receipt. Its scaled positive
coefficient calculation has one normalization per tilt and weight. This
keeps the inner products batched without giving vanished coefficients a
misleading tiny bound. Tests compare it against the retained log-domain method.
The receipt binds the actual maps, outer spectra, and source files.

For the current activation-cover adapter, use new output paths:

```text
python -B workstreams/inner_design/finite_migration/parameter_activation_frontier.py --seed <replayed-activation-cover.json> --replay <its-verification.json> --output <fresh-cover.json> --rounds 6
python -B workstreams/inner_design/finite_migration/verify_activation_bank.py --input <fresh-cover.json> --output <fresh-verification.json>
python -B workstreams/inner_design/finite_migration/parameter_activation_partition.py --seed <fresh-cover.json> --replay <fresh-verification.json> --output <fresh-partition.json> --splits 512
```

Replay the new partition with the same verifier before using it as input to
another refinement or the aggregate checker. These commands neither change
the selected maps nor perform encoder benchmarks.

The sparse extension and suffix-specific refinement use:

```text
python -B workstreams/inner_design/finite_migration/parameter_sparse_extension.py --first 65 --last 256 --output <fresh-sparse.json>
python -B workstreams/inner_design/finite_migration/parameter_sparse_extension.py --input <fresh-sparse.json> --output <fresh-sparse-replay.json>
python -B workstreams/inner_design/finite_migration/parameter_scoped_cover.py --seed <replayed-suffix.json> --replay <suffix-replay.json> --output <fresh-suffix.json> --rounds 4
python -B workstreams/inner_design/finite_migration/parameter_scoped_cover.py --verify <fresh-suffix.json> --output <fresh-suffix-replay.json>
python -B workstreams/inner_design/finite_migration/parameter_scoped_union.py --prefix <bridge.json> <bridge-replay.json> --dense <fresh-suffix.json> <fresh-suffix-replay.json> --output <fresh-full.json>
```

The scoped adapter also accepts the verified bridge as its first seed.
`parameter_scoped_partition.py` bisects a verified suffix without rerunning
the joint optimizer. `parameter_uniform_anchor.py` accepts the same
`--seed`, `--replay`, and `--output` arguments and adds the fair-input anchor.
Replay either result with `parameter_scoped_cover.py --verify` before
combining it with the sparse prefix.

Next: tighten Q>=2 under these same fixed maps, then regenerate all three
paper figures with explicit full-versus-Q1 status. The old RM2Sub curves and
status counts must not be carried into the IMT figures.
