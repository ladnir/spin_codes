# GF(16) Packet Randomization

The [proof index](../PROOF_INDEX.md) identifies all retained ensembles,
their scope, and the checksummed local archive. Shared-shuffle results
must not be combined with the independent-row certificates on this page.

This is a separate research construction, not a production default. Its
strongest complete certificate at K=2^20 gives **relative distance
greater than 10%, with 48.6501 bits of setup-failure margin**, using
four inner updates. Its precomputed transpose measures **7.484 ms**;
the paired two-update version takes **6.881 ms** and has a complete
certificate for distance greater than 9%, with 46.85 bits of margin.
The 10% hill-climb goal is closed. The complete results below concern
ideal independent setup; diagnostics and unfinished variants remain
separately labeled. Production defaults and paper claims are unchanged.

A separate [shared-shuffle exploration](SHARED_SHUFFLE.md) restores one
coordinate permutation per four-row group and adds GF16 packet mixing.
Its four-update implementation measures 6.532 ms against 7.506 ms for the
independent-row control in a matched run. It has only partial proof bounds;
the complete certificates below remain specific to independent row shuffles.
The [shared-route search](PROGRESS.md) stopped at a documented plateau after
closing occupancies 1–22. It does not establish the whole-code 10% claim.

## Complete Certificates (2026-09-29)

### Four Updates: Greater Than 10%

For K=2^20, N=2^21, and the ideal independent setup described below with
four inner updates, the probability that any nonzero message has output
weight at most 209715 is less than 2.264e-15. Outside that bad-setup
event, minimum distance is at least 209716, hence greater than 0.10 N.
The complete margin is **48.6500552955882 bits**, exceeding the target
of 40 bits by about 8.65 bits.

Fresh 256-bit checks closed all 427 dense intervals for q=49--2048.
An independent 384-bit assembly then replayed every interval and
regenerated exact q=1 placement and all support covers for q=2--48.
There are no omitted occupancies or unresolved intervals. The dense
and sparse contributions have margins 90.9171345783 and 48.6500552956
bits, respectively. Their outward sum is at most
2.263987711827718e-15. The exact scope-and-sum audit passes, and the
fresh sparse endpoint exactly matches the separate 384-bit prefix run.
That prefix was also independently regenerated at 256-bit precision.

The complete receipt is `tmp/gf16-r4-d10-assembly-p384-complete.json`,
SHA256 `944c5f41edc0d34c6d2d427f9deaa69674f926cf6eef10e0b788f8fd9f08fdf7`.
The dense input is `tmp/gf16-lowthreads-r4-d10-q49.json`, SHA256
`aaeb719ec69e5c491a5b86418ca035836bf2ca32d8d48c5106bf6499193caf10`.
The independent sparse prefix is `tmp/gf16-r4-d10-q1-48-prefix-p384.json`,
SHA256 `959f60e4b4980e3702e418abcdcde0d7ec4ba1c8cc50081f6340e9e65b8879ee`.
The complete reproduction command appears under **Fresh 10% Assembly**
below. After regeneration, run:

```powershell
python -B research/workstreams/permutation_locality/gf16_packets/audit_complete.py tmp/gf16-r4-d10-assembly-p384-complete.json tmp/gf16-lowthreads-r4-d10-q49.json --distance .1 --bits 48 --prefix tmp/gf16-r4-d10-q1-48-prefix-p384.json
```

This exact audit checks receipt consistency; the assembly command
performs the fresh numerical proof replay. All 214 Python tests pass.
Fresh normal and sanitizer checks also pass for the encoder and the
unchanged two-bit and lane-permutation controls. Source and executable
hash checks bind these implementation checks to the measured version.

### What Now Limits the Bound

The GF16 randomizer makes each fixed nonzero packet uniform over the
15 nonzero packets. The sharper proof retains feedback counts and the
expansion weights of newly activated states, and charges emitted weight
on cancellation events. Regional count bounds handle the intermediate
occupancies. Four inner updates reduce the lazy part of the nonzero-state
transition from 1/4 to 1/16, at a measured cost of about 0.603 ms over
two updates. No additional permutation is needed.

The remaining margin is controlled by the sparse bound, especially q=1.
Its contribution has 48.6551 bits of margin and accounts for about
99.65% of the full upper bound. Within its CDF-majorant calculation,
support 38 contributes about 98.46%; supports 38, 40, and 42 together
contribute about 99.9974%. These are contributions to an upper bound,
not actual BCH shell counts or observed low-weight codewords. The dense
contribution is now negligible. Further proof tightening should target
the lightest-support counting or output-tail bounds, rather than the
already closed dense range. This does not identify an intrinsic limit
on the code's distance.

The paired implementation timings are 6.880930 ms for two updates and
7.4836625 ms for four, an 8.76% increase. Both exclude setup and use
128-bit elements at K=2^20 on Peach core 15. The next engineering step
is to optimize this certified construction while keeping its distribution
fixed. The earlier two-update 9.25% assembly and 9.5% search were stopped
after the 10% closure; their checkpoints and partial logs are retained,
but neither unfinished run is claimed as a complete certificate.

### Four Updates: Greater Than 9.9%

For K=2^20, N=2^21, and the ideal independent setup described below with
four inner updates, the probability that any nonzero message has output
weight at most 207618 is less than 1.166e-15. Outside that bad-setup
event, minimum distance is at least 207619, hence greater than 0.099 N.
The full margin is **49.6084087285827 bits**. The routing, outer, field
randomizers, and inner maps are unchanged from the two-update case;
the number of inner updates is four instead of two.

The fresh 384-bit assembly covers every occupancy q=1--2048. It
regenerates exact q=1 placement and all support covers for q=2--48,
then combines their sum with a fresh replay of all 427 dense intervals
for q=49--2048. The sparse and dense contributions have margins
49.6084087285835 and 90.3779631844015 bits, respectively. No unresolved
interval or occupancy remains. The sparse contribution therefore
controls this bound; the dense contribution changes its margin by
less than 1e-12 bits.

The dense witnesses first passed 256-bit checks before this 384-bit
replay. The sparse prefix was also independently regenerated at 256
and 384 bits. Its earlier 384-bit endpoint exactly matches the sparse
endpoint regenerated by the full assembly. An exact receipt audit
checks the dense input hash, exhaustive disjoint binary-path partition,
matching setup parameters, consecutive occupancy ranges, and outward
addition of the two contributions.

The complete receipt is `tmp/gf16-r4-d099-assembly-p384-complete.json`,
SHA256 `b198aeb447bfb6b474033a144edfbdb4b603536e1e71dd954e0d75d5bde30fc8`.
The dense input is `tmp/gf16-parallel-r4-d099-q49.json`, SHA256
`ae1cba7f4b1b89b26f3b8afccaf87a546e4298e23925683d76ae48abcba72e65`.
To regenerate the complete bound:

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-parallel-r4-d099-q49.json --precision 384 --dense-workers 4 --dense-log tmp/gf16-r4-d099-assembly-p384-dense --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-target-bits 48 --sparse-small-through 8 --sparse-small-tilts .00016 .00024 .00028 .00032 .0004 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 --sparse-analytic-gradient --sparse-workers 1 --sparse-log tmp/gf16-r4-d099-assembly-p384-sparse.log --output tmp/gf16-r4-d099-assembly-p384-complete.json
```

The paired four-update timing is 7.4836625 ms, versus 6.880930 ms for
two updates: about 8.76% extra time. Existing lane-permutation and
two-bit controls are preserved. Neither production defaults nor paper
claims change automatically with this research certificate.

After regeneration, `audit_complete.py` checks the dense input hash,
occupancy coverage, disjoint exhaustive binary-path partition, exact
outward sum, and requested distance and margin. An optional independent
sparse prefix must have matching scope and an identical endpoint:

```powershell
python -B research/workstreams/permutation_locality/gf16_packets/audit_complete.py tmp/gf16-r4-d099-assembly-p384-complete.json tmp/gf16-parallel-r4-d099-q49.json --distance .099 --bits 49 --prefix tmp/gf16-r4-d099-q1-48-prefix-p384.json
```

This audit passes for the complete 9.9% receipt. Nine regression tests
cover mismatched claims, altered input hashes, missing or overlapping
intervals, incorrect sums, strict thresholds, and mismatched independent
prefixes. The complete suite passes all 214 tests; the local log is
`tmp/gf16-final-audit-full-tests.log`. This audit checks receipt consistency only; it does not replace the
fresh numerical proof replay performed by `assemble.py`.

### Two Updates: Greater Than 9%

For K=2^20, N=2^21, and the ideal independent two-update setup described
below, the probability that any nonzero message has output weight at
most 188743 is less than 7.869e-15. Outside that bad-setup event, minimum distance is at
least 188744, hence greater than 0.09 N. The complete margin is
46.8529233246 bits.

The fresh 256-bit assembly regenerates exact q=1 placement and every
sparse occupancy q=2--96. It also replays all 116 dense cells for
q=97--2048, with no unresolved cells. The sparse and dense contributions
have margins 46.8529233249 and 78.7739845636 bits, respectively.
This result uses birth-class operators and exact feedback, without the
new joint-return or lazy-density refinements. The encoder is unchanged.

The complete receipt is `tmp/gf16-birth-classes-d09-smallgrid-complete.json`,
SHA256 `5788b4b383e90774e3d4a1278134488915df48c677eecd4150d3ea2d7e373920`.
The dense input is `tmp/gf16-birth-classes-d09-tilt3over16-q97.json`, SHA256
`d378dc4f7185df56d484e21aed64b74503ece068a80bf2dd29e600e9c440a213`.
An independent full 384-bit replay has completed. It regenerates all
sparse occupancies and replays the dense partition, reproducing
46.8529233245806 bits of margin. Exact dyadic comparisons confirm that
its dense, sparse, and total upper endpoints are no larger than their
256-bit counterparts. Its receipt is
`tmp/gf16-r2-d09-independent-p384-complete.json`, SHA256
`203aa47128541ec8e3a523ad8a23f5608549bfc6f12d20fe7a05ee4dc7e5eb7f`.
The older receipt omits an explicit update-count field; both receipts
bind the same dense input by hash, and that input specifies two updates.
Neither the production parameters nor the paper claims have changed.

### Earlier 8.2% Certificate

For K=2^20, N=2^21, and the ideal independent setup described below, the
probability that any nonzero message has output weight at most 171966 is
less than 1.626e-16. Outside that bad-setup event, minimum distance is at
least 171967, hence greater than 0.082 N. The replayed margin is
52.4502533222 bits.

The complete sum contains three disjoint ranges:

- q=1: exact ordered regional placement, regenerated by `single_group.py`.
- q=2--104: freshly generated GF fixed-occupancy operators, exact feedback
  atoms through eight packets, and complete support covers for every
  occupancy. Together with q=1, aggregate log2 upper -52.4502533577.
  The operators include the `gf-rank` return bound.
- q=105--2048: 173 dense composition cells using `rank-return` and BCH
  row parity, with no unresolved cells; aggregate log2 upper -77.7295645136.

Both parts were evaluated at 256-bit precision. The dense witnesses were
replayed from their exact rational parameters; the sparse part was
regenerated without an operator cache. The stopping budgets deliberately
seek a sufficient bound, not the tightest possible margin. This is not a
measurement of the code's typical minimum distance.

The local dense cover is `tmp/gf16-parity-rank-d082-q105.json`, SHA256
`c4a33fed9e5192fe4b8d9c7157217cb4e0dd1848a00a32377a26a22f7d41cb66`.
The complete replay receipt is `tmp/gf16-parity-rank-d082-complete.json`, SHA256
`29cc3717e976d118cb8a4ed302b695f8652c958b8aac3d956cb57fb8da9ccfb6`; the commands
below regenerate and verify the result. No production code or paper
claim is changed by this research certificate.

For the earlier 8% certificate, a second run without the within-run CDF
fold cache gives identical dyadic dense, sparse, and total bounds. Its receipt is
`tmp/gf16-d08-rank-complete.json`; the JSON differs only because that
earlier program version omitted the false `single_group_exact` field.

Earlier complete replays remain available:

| Distance | Setup-failure margin | Local replay receipt |
|---|---:|---|
| >4% | 56.9744 bits | `tmp/gf16-d04-complete.json` |
| >4.5% | 56.1980 bits | `tmp/gf16-d045-complete.json` |
| >5% | 49.7438 bits | `tmp/gf16-d05-shape-complete.json` |
| >6% | 72.2350 bits | `tmp/gf16-d06-exact-complete.json` |
| >6.5% | 60.2239 bits | `tmp/gf16-d065-exact-complete.json` |
| >7% | 55.0303 bits | `tmp/gf16-d07-exact-complete.json` |
| >7.5% | 59.4211 bits | `tmp/gf16-d075-profile-complete.json` |
| >7.75% | 55.2284 bits | `tmp/gf16-d0775-conditioned-complete.json` |
| >8% | 52.4780 bits | `tmp/gf16-d08-rank-cached-complete.json` |

The margins are not monotone because the proof bounds and stopping
budgets differ. They do not estimate the code's actual distance or its
optimal failure bound. All rows concern the same two-update GF16 encoder.

## Construction and Exact Local Guarantee

Start with the independent-row four-bit route: BCH[256,128], four adjacent
rows per group, an independent uniform coordinate permutation per row,
and an independent uniform packet permutation in each of 256 regions.
Use the IMT(128,19) inner with zero initial state and no flush. The
original candidate uses two updates per step; the stronger certificate
above uses four. The update count is fixed for each construction.

Replace each packet's lane permutation with multiplication by an
independent uniform element a of GF(16) excluding zero. Interpret the four
binary coordinates in the polynomial basis modulo X^4+X+1. The multiplier
is sampled during setup, after routing, and is fixed for all messages.
All these multipliers are independent of the other setup randomness.

For any fixed nonzero packet v, the map a -> av is a bijection on the 15
nonzero packets. Thus av is uniform over those packets; zero stays zero.
Multiplication by a fixed nonzero a is binary-linear and invertible. It
preserves whether a packet is active, but does not preserve its binary
weight. This change therefore preserves rate, not an existing distance
certificate.

Fix a message before sampling setup and condition on the row and regional
permutations. The transformed active packets are independent and uniform
nonzero packets. Their weight probabilities, for weights 1 through 4, are

    (4, 6, 4, 1) / 15.

In particular, an all-one packet is no longer fixed by the packet
randomizer. These statements concern the ideal independently sampled
setup. No independence between different messages is needed when summing
their failure probabilities.

## What This Simplifies in the Proof

The existing positive BCH row-mixture envelope still dominates each
shuffled row's counting measure. Pushing both sides through the independent
GF randomizers preserves that domination. For a four-row mixture component
with packet-weight law pi, only p=1-pi(0) now matters: each packet is zero
with probability 1-p, and otherwise uniform nonzero.

Equal p values with the same active-component label can be merged by
adding their positive counting coefficients. This reduces 70 components
to 36. The labels refer to nonzero message components, not whether a
particular sample from a comparison measure happens to be nonzero.

`scalar_cover.py` applies the existing heterogeneous shuffle inequality
to the **two activity categories**, not the five packet weights. For a
positive activity tilt tau, component i has

    Z_i = 1 - p_i + tau p_i,
    f_i = tau p_i / Z_i,
    c'_i = c_i Z_i^256.

For G=2048 groups with component counts n_i, the sole cover coordinate is
x=sum_i n_i f_i/G. At the base tilt, the unnormalized reference weights
are w_0=1-x and w_1=x/tau. Dividing these weights by their sum gives an iid
activity reference; conditional nonzero packets are uniform among 15
values. Reference mass is restored once per packet, and the regional
shuffle loss is paid in each of 256 regions. These reference inputs are
proof devices, not a change to the actual input distribution.

The outer generating function sums all component assignments in a cell,
including their multiplicities. Nonnegative occupancy duals enforce the
selected lower bound on active groups. Alternate tilts use rational affine majorants
checked against every component; floating optimization only proposes
witnesses. Both the original two-state inner envelope and the new
three-coordinate refresh envelope are evaluated with Arb.

### GF-Specific Feedback and Cancellation

The current proof combines two improvements from
[the fixed-occupancy analysis](OCCUPANCY_BOUND.md). First, integer Walsh
inversion counts every feedback value through eight active packets.
For one and two packets, the largest nonzero atom after averaging their
positions is 1/480 and 5/111600. Second, cancellation must still emit
weight: with j active packets, every nonzero state emits at least
max(0,48-4j) bits in that step. The bound uses this pointwise constraint
on the cancellation event, without assuming independence.

The `gf-refresh` dense kernel averages these conditional bounds under
the iid comparison's binomial packet occupancy. This closes full 6% and
6.5% certificates without changing the encoder. A further
`weighted-return` bound retains the input weight inside the feedback
character sum; see [the refresh note](REFRESH_BOUND.md). It closed 7%.

[The profile bound](PROFILE_RETURN.md) then enforces the expansion's
minimum weight inside that character calculation. Convexity and
rearrangement give a worst-case profile with twelve full packets and
twenty zero packets. This closes the earlier full 7.5% result, without
changing the encoder or adding another inner update.

The earlier `weighted-return` dense searches, both starting at q=97,
were incomplete:
7.5% retained 74 accepted cells and 103 unresolved cells; 8% retained
37 accepted cells and 177 unresolved cells. Each used a 250-cell budget.
At 7.5%, a zero-width diagnostic near coordinate x=0.0500141 still gives
binary64 log2 upper about +2840.17. Thus subdivision alone does not resolve
this particular comparison bound. This is not a low-weight codeword or
an upper bound on the construction's distance. The profile refinement
resolves the 7.5% gap. At 8%, `profile-return` still leaves 141 cells
unresolved after 300 visited cells, alongside 80 accepted cells.

The next `shape-return` refinement retains all 19 actual expansion
histograms. Ten states have the extreme twelve-full-packet profile;
they are handled separately using the actual feedback columns at those
positions. Their signed weighted-return coefficients can be enumerated
exactly. A cheaper Cauchy bound on two character moments is sufficient
for the current default; exact enumeration remains a diagnostic. The
other 18 histogram classes use rearrangement bounds with their actual
weights and multiplicities. At a representative 7.5% comparison point,
this gains about 471 additional bits over `profile-return` in binary64.
The 8% dense search with this kernel visited 300 cells, accepting 86 and
leaving 129 unresolved. This is not a whole-code result.

The next `conditioned-return` refinement conditions on the number of
active packets before applying Holder's inequality. It uses moments of
orders 2,3,4,6,8,12,16, and also retains each expansion histogram before
taking maxima or averages. The argument is recorded in
[the return-bound note](PROFILE_RETURN.md). Exhaustive small-input and
eight-step transition checks pass. The 7.75% dense cover now has 158
accepted cells and no unresolved cells. A fresh 256-bit replay gives
aggregate log2 upper -120.1541679059 for q=97--2048. Regenerating every
occupancy in the matching sparse prefix completes the 7.75% certificate
reported above.

The last narrow 7.75% gap was numerical search loss, not a failed outward
bound. The two-variable outer-dual optimizer sometimes stopped too early.
Eliminating the binary occupancy multiplier analytically reduces the
remaining optimization to a monotone one-dimensional derivative. At
x=0.0412 this changes the proposed log2 bound from +421.44 to -299.19;
the latter was checked outward at 256-bit precision. The exact dual
inequality and the interpretation of saved witnesses are unchanged.

Two further refinements target 8%. `trimmed-return` counts outputs by
weight exactly and assigns a possible cancellation event the lightest
outcomes permitted by its integer event-count bound. `rank-return` also
bounds the number of free bits after fixing feedback, while charging for
the unchanged output outside active packet positions. Both arguments are
in [the return-bound note](PROFILE_RETURN.md). Their minimum is taken at
each fixed occupancy and expansion profile before averaging.

The new bound is available both in the dense iid comparison and in the
`gf-rank` fixed-occupancy operators. The full 8% dense search closed with
115 accepted cells and no unresolved cells. Its fresh 256-bit replay
gives aggregate log2 upper -73.0357614708 for q=97--2048. The matching
sparse prefix has now been regenerated, completing the 8% result above.
Separate full-support checks at
q=64,96,128 give margins 65.67,95.46,56.65 bits, respectively; these
selected occupancies are not a complete prefix.

At 10%, retuning and fully subdividing the earlier occupancy-one cover
still gave only about 35.85 bits. `single_group.py` removes its Bernoulli
comparison across regions by averaging ordered region-operator products
over exact-size supports. With the same encoder and local state bounds,
the resulting outward bound gives **42.1226 bits for all one-active-group
messages at 10% distance**. This is not a whole-code certificate. The
argument and its CDF counting step are recorded in
[the return-bound note](PROFILE_RETURN.md). The remaining 10% challenge
includes the intermediate and dense occupancies.

The optional [row-parity comparison](PARITY_COMPARISON.md) retains the
BCH rows' even weights after the field randomization. Its exact positive
factors reduce some dense counting losses without changing the encoder.
At selected 8.5% comparison points this gains about 373--631 bits, but
the resulting bounds still leave a gap. These floating diagnostics do
not certify 8.5%. The current complete 8% dense cover does not rely on
row parity.

The 8.1% search without parity first exhausted 400 proposals with 187
accepted cells and 27 pending cells. A rechecked resume closes its dense
partition with 222 accepted cells and no unresolved cells. Its older
full replay verified every sparse occupancy q=2--96 but left q=1
unresolved, so the assembler correctly emitted no complete receipt.
The separate exact-placement method covers q=1 through 10%; the newer
full replays use that method instead. The 8.1% run is not being repeated
while stronger replays are active. The separate 8.2% search using row parity closes the
entire dense range q=105--2048 with 173 accepted cells and no unresolved
cells. Its full replay has completed, giving the 8.2% certificate above.

`fiber_refresh_probe.py` tests another analysis option: representing lazy
feedback as a uniform density after enough active packets. Small exact
input enumerations and six-step domination checks pass, but the current
floating probe does not improve the limiting 8.5% points. It remains a
diagnostic and is not used by any certificate.

The newer [fresh-state density bound](BIRTH_DENSITY.md) instead sharpens
the zero-source row. It keeps low packet occupancies arbitrary and uses
a centered Fourier sum to bound the weighted feedback density from the
remaining inputs. At 8.5%, selected dense points improve by hundreds to
thousands of bits. A real 256-bit check at coordinate 0.04 gives log2
upper -119.0705264463, agreeing with the floating proposal. The full
`birth-refresh` dense search at q>=129 stopped after 500 proposals with
249 accepted cells and three unresolved cells. The variance refinement
below closes a larger dense range; there is no full 8.5% certificate yet.
The encoder and performance are unchanged.

The `birth_classes_probe.py` diagnostic retains the expansion
weight of low-occupancy births and uses the density envelope for the
rest. Its exact integer character census and small multistep tests pass,
but selected 9% scores remain positive. The new `birth_classes.py`
kernel adds outward evaluation and class-specific return bounds;
the dense search and assembler accept it as `birth-classes`.
No complete certificate yet uses this refinement. The derivation and a numerical
loss decomposition are recorded in [the fresh-state note](BIRTH_DENSITY.md).
At the tested 9% witness, inner returns and the shuffle comparison both
contribute substantial slack; eliminating either one alone does not
close the bound.

[The variance-sensitive shuffle bound](SHUFFLE_VARIANCE.md) bounds the
ratio between a shuffled heterogeneous
Bernoulli family and its iid reference using the variance of the active
packet count. The optional `--variance-shuffle` path checks a rational
variance minorant for every component in each cell, then bounds the
density ratio over its entire mean interval with outward arithmetic.
Deterministic types remain included; alternate tilts retain the old
comparison. At 8.5%, a selected 256-bit check gives log2 upper
-1020.3369136198. The full dense search starting at q=97 closes with
155 accepted cells and no unresolved cells. Independent assembly
replayed that partition but left q=4 unresolved under its coarse grid. A separate
full-support check at q=96 gives 61.6192 bits, but is not imported as a
proof input. The later fine-grid 9% assembly supersedes this attempt.
The fresh 256-bit dense replay gives aggregate log2 upper
-72.9420481758 before adding the sparse prefix.

The separate `variance_split_probe.py` diagnostic partitions compositions
by their Bernoulli variance and charges each part its own outer count.
Combined with the birth-class diagnostic, it lowers the selected 9%
score at coordinate 0.04 from +1806.99 to +758.49. This still does not
close the point. These are floating diagnostics, not a 9% certificate; see
[the shuffle note](SHUFFLE_VARIANCE.md) for the dual inequality and scope.

The new `--variance-bins 8` option instead applies only the variance
partition, without birth classes, and checks whole mean intervals with
outward arithmetic. At 8.6%, the interval [0.029999,0.030001] improves
from log2 bound +64.9955 to -256.7044; 256- and 384-bit evaluations agree.
The full dense cover closes with 154 accepted cells and no unresolved
cells; its full assembly was stopped after finding the same coarse-grid
q=4 gap. Separate all-support sparse checks at
q=64 and q=96 now give 79.5732 and 63.1747 bits, respectively, using the
additional output tilt .128. They do not cover the other occupancies
and are not a complete 8.6% certificate.

With `birth-classes` and eight variance bins, the complete 8.7% dense
cover closes in 156 cells, with no unresolved cells. The fixed-occupancy
version, `gf-birth-classes`, now preserves newborn expansion-weight
classes in the sparse argument as well. Exact small-instance tests
include mixed sequences of empty and active steps.

The earlier q=48 search failed with its coarse output-tilt grid. With
the expanded grid below, the rank baseline gives 63.0187 bits; the new
class bound gives 175.1443 bits. Fresh 256- and 384-bit class evaluations
agree to the displayed precision. These runs use different stopping
points, so their margins do not quantify an optimal improvement.
Separate baseline checks at q=32 and q=96 give 56.3443 and 54.6489 bits.
A full 8.7% assembly is regenerating every sparse occupancy and replaying
the dense partition. Until it completes, these remain partial results.

`frontier_probe.py` provides floating-point diagnostics, not certificates.
At an 8% comparison point x=0.04, suppressing lazy returns in the bound
changes its score from approximately +1232 to -1411. Adding a third
transvection update at the same witness instead gives about -394. These
are counterfactual comparisons, not a new construction certificate or
runtime measurement. They identify cancellation as a useful target for
the next bound, not an observed low-weight codeword.

The earlier 192-bit 8% sparse checks verified q=1,64,96 but left q=128
unresolved after five subdivisions. The rank refinement and larger work
budget now resolve that selected occupancy. Neither set of selected
checks substitutes for a complete sparse prefix.

### Earlier Hill-Climb: Retain Refresh and Split the Occupancies

[The refresh bound](REFRESH_BOUND.md) separates arbitrary nonzero-state
mass from a component dominated by the uniform nonzero distribution.
The latter uses the average expansion moment instead of the worst state.
This is a sharper analysis of the same two-update inner, not an additional
mixing operation. Exhaustive small-state tests check the separate lazy
mass and refreshed density bounds and eight successive transitions.

[The direct sparse bridge](SPARSE_TRANSFER.md) transfers the old shape-uniform
argument only with all-one penalty and input-weight tilt both equal to
one. GF multiplication preserves packet activity; conditional on the new
weight, an active packet is uniform within that weight shell. The bridge
regenerates the bounds for this ensemble and does not reuse the old
weighted lane-permutation certificates.

The dense searches establish the following partial covers. Here q is
the number of groups containing a nonzero four-row message, not the
number of active packets in an inner step.

| Distance | Dense occupancy range | Base activity tilt | Leaves | Replayed log2 upper |
|---|---|---|---:|---:|
| 3% | 65--2048 | 1/32 | 69 | -73.17538780 |
| 4% | 65--2048 | 1/8 | 30 | -74.24418257 |
| 4.5% | 73--2048 | 1/8 | 63 | -76.12591692 |
| 4.75% | 85--2048 | 1/8 | 44 | -79.01810146 |
| 5% | 90--2048 | 1/8 | 153 | -49.83553211 |
| 6% | 193--2048 | 1/8 | 28 | -210.55523123 |
| 7% | 385--2048 | 1/8 | 22 | -451.87699094 |

All listed dense covers have zero unresolved cells and were independently
replayed at 256-bit precision. This table records earlier kernels; its
6% and 7% rows are partial results, superseded where indicated above.
The 4%, 4.5%, and 5% rows have matching complete sparse coverage.
The 5% dense cover uses central mass (1001/1000) 2^129 and row
bias 21/50. The 4.75% row uses 2^129 and bias 17/40. The other listed rows use 2^130
and bias 2/5. These comparisons are checked against every BCH shell cap.

The first 5% run did not establish a whole-code
claim. Its first full replay verified q=1--65 but did not close q=66
within 100 support subdivisions (binary64 aggregate log2 upper +38.14).
That incomplete run was stopped before testing every larger occupancy.
The earlier full 3% attempt is
superseded by the complete 4% result.

The failed searches are informative. At 4%, base tilt 1/32 left 227
unresolved cells within a 300-cell budget; changing the proof coordinate
to base tilt 1/8 closed the same occupancy range. At 5%, searches from
q=65 with base tilts 1/8 and 1/4 left 417 and 243 unresolved cells,
respectively. Starting the dense part at q=129 closed it in 23 leaves.
These are changes in the proof decomposition, not changes to the code or
evidence of a low-distance word.

### Earlier Middle-Occupancy Gap

The `feedback-refresh` option requires nonzero feedback for a refreshed
state to return to zero; the derivation is in [the refresh note](REFRESH_BOUND.md).
It passes exhaustive small-instance tests, but its improvement alone is
too small to close the 5% dense cover starting at q=65.

Starting the original 5% dense cover at q=97 closed it in 68 cells. A fresh
256-bit replay gave log2 upper -74.1669826743. The sharper row comparison
closed q=90--2048, as listed above. Together with the then-verified sparse
prefix, this narrowed the interval needing coverage to **66--89**.
The later complete 5% replay checked every occupancy in that interval;
individual checks were not interpolated.

The sparse search was retuned beyond its original output-tilt grid.
The original largest tilt was .032. A screen using only smaller tilts
.0064, .0128, .0256 failed badly at q=66 even after 20 subdivisions;
tilt .064 instead verified q=66 with more than 1138 bits of margin.
Increasing the grid alone did not close q=80 or q=96 in ten subdivisions.

The new `gf-shape-tilt` option applies an exact change of measure to the
GF packet values; its proof is in [the sparse bridge](SPARSE_TRANSFER.md).
It retains the old shape-weighted state envelope and pays an explicit
normalizer per active packet. At penalty rho=3/4 and input-weight tilt
one, that normalizer is only 46/45. Outer counts remain unweighted.
At 5%, this method verified q=80 with 1076.52 bits after one subdivision,
and q=88 with 57.14 bits after thirteen subdivisions. Both used 192-bit
outward arithmetic. The q=90 screen remained unresolved after twenty
subdivisions (binary64 aggregate log2 upper +159.01). These are complete
support covers for the stated occupancies only.
The q=89 run with output tilt .08 verified all supports at 192-bit
precision, with 54.05924728 bits of margin after thirteen subdivisions.
A full 5% assembler run regenerated all sparse operators at 256-bit
precision, without memoized inputs, and closed with 49.7438455766 bits.
Its receipt is `tmp/gf16-d05-shape-complete.json`.

[The fixed-occupancy note](OCCUPANCY_BOUND.md) gives new GF-specific local
operators. They average all 15 nonzero packet values, retain a density cap,
and optionally retain expansion-weight classes through empty steps. Exact
small-instance tests pass. The first versions failed intermediate
support-cover screens. Exact feedback counts and the expansion-distance
floor subsequently made the three-coordinate version sufficient for the
complete 6% and 6.5% certificates. The additional density and weight-class
variants remain available as experiments.

The dense driver records the row comparison's central mass and bias.
Its central mass can be 2^b times a positive rational scale, avoiding
unnecessary rounding to a power of two. Saved witnesses record that scale;
older records default to one.
Every such comparison is checked against all authenticated BCH shell caps.
Reducing the central exponent from 130 to 129, or changing the bias from
2/5 to 1/3, did not close the 5% cover from q=65 within 300 cells. These
comparison parameters are proof choices, not encoder parameters.

### First Whole-Domain Attempt: Incomplete

The first bounded search used cutoff floor(0.05 N)=104857, two inner
updates, all occupancies 1 through 2048, a budget of 300 cells, and
192-bit outward arithmetic. It retained 26 accepted cells and **249
unresolved cells**. A separate replay checked the partition and the
accepted bounds. Their summed log2 upper is approximately -247.29, but
that sum excludes the unresolved cells and is not a failure margin for
the code.

The original difficulty was not only interval width. A diagnostic at the
fixed base-tilt coordinate x=1/64 gave log2 upper about +8933 at 5%.
The largest binary shuffle penalty is about 1491.44 bits across all
regions, so removing that entire penalty would not rescue this particular
bound. The later inner refinements resolve that attempt's obstruction;
the failed bound was not evidence of an actual low-weight codeword.

The saved attempt is local and ignored by Git:
`tmp/gf16-packets-d05-first.json`, SHA256
`d96085e40b5aac8f11e908f6ad41e140419631f8bfe5f2965e897c54ac9f578a`.
Old lane-permutation certificates are not imported by the verifier.

## Exact Feedback Census

`packet_census.py` authenticates the maps through the existing map loader
and examines the four feedback columns at each of the 32 packet positions
within an inner step. All 32 four-column maps are injective. Every one of
the 496 pairs of distinct positions has joint rank eight.

Consequently, one or two active packets cannot give zero feedback. Under
independent GF randomizers their largest feedback atom is respectively
1/15 or 1/225. The rank statements also hold for the old map; the new
benefit is uniformity over all nonzero packet values, regardless of the
original packet shape.

More generally, if j packet positions have joint binary rank r, any
feedback target has at most 2^(4j-r) preimages among the unrestricted
4j-bit inputs. Restricting each packet to be nonzero cannot increase that
count. Independent uniform nonzero packets therefore give the upper
2^(4j-r)/15^j on any feedback atom. Exhaustive position-set ranks give:

| Active packet positions | Rank histogram | Worst-case atom upper |
|---|---|---|
| 1 | rank 4: 32 | 1/15 |
| 2 | rank 8: 496 | 1/225 |
| 3 | rank 12: 4933; rank 11: 27 | 2/3375 |
| 4 | rank 16: 31700; rank 15: 4180; rank 14: 80 | 4/50625 |

These are unweighted local feedback bounds. They are not output-distance
bounds and cannot be multiplied by a separate output-weight moment as
though feedback and output were independent.

## Implementation and Measurement

`../Gf16Packet.h` implements the adjoint multiplier across four opaque
128-bit elements using four fixed masked AVX-512 shuffles and XORs.
It performs no GF(2^128) arithmetic. The coefficient's binary matrix is
precomputed. The hot path has no variable-length loops, lookup gathers,
allocations, or data-dependent branches. `../joint.cpp` modes 7 and 8
select GF16 streaming and cached stores, respectively. Existing modes
retain their old construction.

The full encoder is checked against materialized routing, direct
polynomial-basis matrix multiplication, and the existing BCH reference.
The SIMD adjoint is also checked exhaustively for all 15 multipliers and
all 16 binary packets. Normal-build checks passed for both new store modes,
two and three updates, seeds 1 and 17, and K=2^14, 2^18, 2^20: 24 cases,
each with three input patterns. Two unchanged-route controls also passed.
The same 26-case check passed with the research encoder compiled under
AddressSanitizer and UndefinedBehaviorSanitizer; its linked prebuilt BCH
objects and library were not instrumented in that run.
The Python tests cover field action, the binary shuffle comparison,
exact small inner transitions, local feedback counts, alternate weight
majorants, outward evaluation, refresh-density propagation, sparse and
assembler scope validation, exact fixed-occupancy moments and feedback,
density and expansion-class propagation, exact shape changes of measure,
outward region compensation, rational row comparisons, exact Walsh
inversion, joint return counts, input-weight-tilted cancellation,
profile rearrangement, the exceptional-state character census,
occupancy-conditioned Holder bounds, exact output-weight tails,
raw feedback-rank bounds, exact regional placement for one group, and
parity-filtered activity measures, variance-sensitive shuffle comparisons,
whole-cell variance minorants, grouped positive probability products,
coverage retention under a work limit and independently rechecked resumes.
Transition tests check zero-return
entries directly as well as multistep total moments.

The sparse search batches repeated Bernoulli probabilities and reuses
the vectorized CDF objective when proposing witnesses. Outward replay
evaluates the same exact Bernoulli-count product with grouped Arb
polynomials; small exact-rational tests enclose every coefficient,
including probabilities zero and one. These changes reduce proof
computation only. They do not alter the encoder, its setup distribution,
or the distance claim being checked.

Timings use the existing Peach machine, core 15, GCC -O3, 128-bit elements,
K=2^20, and precomputed setup. All benchmarks run serially under the three
shared benchmark locks. After a 31-call screen, confirmation used two
setup seeds and three 101-call runs per seed, reversing mode order in the
middle run. The table reports the median of the six run medians.

| Encoder | Precomputed transpose |
|---|---:|
| Four-bit lane permutation, two updates | 6.432 ms |
| Four-bit GF16 randomizer, two updates, streaming | 6.862 ms |
| Four-bit GF16 randomizer, cached stores | 12.688 ms |
| Two-bit padded tiled control, two updates | 7.457 ms |

The useful GF16 version costs about 6.7% more than the four-bit control
and remains about 8.0% faster than the two-bit control. The profiled screen
locates almost all the added cost in the fused inner/routing pass
(about 3.18 ms versus 2.73 ms); repacking and BCH remain about 3.7 ms.
This is not yet the approximately 5 ms target.

Remote logs are under `/tmp/spin-gf16-w2vMCK/measurements/`, with confirmation
directory `gf16-confirm-fJk2vC`. The timed executable SHA256 is
`f7bc2ee59d25752cfc5ce7a530a48682c65d19d70aeb26e2670e87d52051b677`.
Generated measurements and proof witnesses remain outside tracked source.

### Three-Update Comparison

The existing binary also supports three inner updates. A paired run on
2026-09-29 used the same machine, core, binary hash, and precomputed
GF16 streaming mode. For each of seeds 1 and 17, it measured update
counts 2,3,3,2, with 101 calls per measurement. All runs were serial
under the same three benchmark locks; their implementation checks
passed. The medians of the four run medians are:

| Inner updates | Precomputed transpose |
|---|---:|
| Two | 6.872543 ms |
| Three | 7.193953 ms |

The extra update costs about 0.3214 ms, or 4.68%, in this paired trial.
It does not require a new permutation or GF randomizer. These timings
are not accompanied by a complete three-update distance certificate yet.

### Four-Update Comparison

A second paired trial adds four updates without changing the permutation,
GF randomizer, or fixed inner maps. It uses one freshly compiled binary
for all three update counts. On Peach core 15, for each of seeds 1 and 17,
the serial order is 2,3,4,4,3,2, with 101 precomputed calls per run.
Setup, allocation, and reference checks are outside those calls.
The medians of the four run medians are:

| Inner updates | Precomputed transpose | Increase over two |
|---|---:|---:|
| Two | 6.880930 ms | -- |
| Three | 7.184282 ms | 4.41% |
| Four | 7.483663 ms | 8.76% |

Four updates cost another 0.6027 ms over two in this trial. The complete
certificate above now gives distance greater than 10% for this
four-update construction, with 48.6501 bits of setup-failure margin.
The two-update construction retains its separate greater-than-9%
certificate, with 46.85 bits of margin.

`../run_gf16_updates.sh ROOT check|sanitize|confirm` reproduces the checks
or timing order when ROOT contains `joint` and `joint-sanitize`. The
script holds all three shared benchmark locks. Regular and ASan/UBSan
checks passed for updates 2,3,4, both GF16 modes, both seeds, and the
unchanged lane-permutation and two-bit controls. All timed full-encoder
and adjoint checks also passed. The four-update dispatch is restricted
to the four-bit GF16 modes; production defaults are unchanged.

The optimized binary is `/tmp/spin-gf16-r4-Cpp6Wn/joint`, SHA256
`38d0032d8386227ed59c11baed66a1dbe57225dec9b0c7eb89db742ada0bf151`.
The local log is `tmp/gf16-r2-r3-r4-abccba.log`, SHA256
`b5c6cd30a9372b6250d708cec157ba405b7f9a2d6c4a9d026936fbcb4adbd3f2`.
No phase timers run inside the measured calls.

A read-only completion audit rechecked the remote executable hash and
all 15 normal-build and 15 sanitizer check records. They are under
`/tmp/spin-gf16-r4-Cpp6Wn/measurements/gf16-updates-check-hNQBbR/` and
`/tmp/spin-gf16-r4-Cpp6Wn/measurements/gf16-updates-sanitize-ynpBxf/`.
Each set includes the four-update larger-size case and the unchanged
lane-permutation and two-bit controls. The sanitizer executable's SHA256
is `9270a6b320eaafc298d9a7fe6fd8f27e1039bb1f01aeb15ac30ce3b1e89187cc`.
This audit did not rerun or overlap benchmarks.

The final audit also reran both checks-only modes serially. Each script
completed with exit status zero, covering all three input patterns,
the larger four-update case, and both unchanged controls. The fresh logs
are in `gf16-updates-check-RsR6YL/` and `gf16-updates-sanitize-ffeSGj/`
under the same measurements directory. Before running them, exact hashes
confirmed that the remote `joint.cpp`, `Gf16Packet.h`, `TwoBitTiled.h`,
and check script match the local sources. No timed calls were requested.

## Reproduce and Next Step

From the repository root:

```text
python -B -m unittest discover -s research/workstreams/permutation_locality/gf16_packets -p test_*.py
python -B research/workstreams/permutation_locality/gf16_packets/packet_census.py
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 1/20 --max-cells 300 --output tmp/gf16-packets-d05-first.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --replay tmp/gf16-packets-d05-first.json
```

The last two commands reproduce an incomplete proof attempt, not a full
certificate. `../run_gf16.sh ROOT check|screen|confirm` runs implementation
checks or serialized measurements when ROOT contains the built `joint`.

The refreshed-state searches and full replay use:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 1/25 --kernel refresh --base-tilt 1/8 --minimum-groups 65 --max-cells 300 --max-depth 22 --output tmp/gf16-refresh-d04-q65-base8.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-refresh-d04-q65-base8.json --precision 256 --sparse-log tmp/gf16-d04-sparse-replay.log --output tmp/gf16-d04-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 9/200 --kernel feedback-refresh --base-tilt 1/8 --minimum-groups 73 --max-cells 250 --max-depth 22 --output tmp/gf16-feedback-d045-q73.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-feedback-d045-q73.json --precision 256 --sparse-tilts .00032 .0032 .032 .064 --sparse-log tmp/gf16-d045-sparse-replay.log --output tmp/gf16-d045-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 1/20 --kernel refresh --base-tilt 1/8 --minimum-groups 129 --max-cells 300 --max-depth 22 --output tmp/gf16-refresh-d05-q129-base8.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-refresh-d05-q129-base8.json --precision 256 --sparse-log tmp/gf16-d05-sparse-replay.log --output tmp/gf16-d05-complete.json
```

`assemble.py` accepts only a complete dense partition for this GF16
ensemble with a matching two-, three-, or four-update claim. It regenerates the sparse operators without
memoized data. Its final JSON is a replay receipt, not a proof input:
verification reruns the command. Generated covers and receipts remain
ignored under `tmp/`.

The new shape-tilt experiment and the verified 5% dense range can be
reproduced separately:

```text
python -B research/workstreams/permutation_locality/gf16_packets/sparse_cover.py --inner gf-shape-tilt --occupancies 80 88 --threshold 104857 --tilts .064 .096 --shape-penalties .75 .5 --max-splits 20 --precision 192 --output tmp/gf16-d05-shape-replay.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 1/20 --kernel feedback-refresh --base-tilt 1/8 --minimum-groups 91 --central-bits 129 --row-bias 17/40 --target-bits 52 --max-cells 300 --max-depth 22 --output tmp/gf16-feedback-d05-q91-bias425-t52.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --replay tmp/gf16-feedback-d05-q91-bias425-t52.json --precision 256
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 1/20 --kernel feedback-refresh --base-tilt 1/8 --minimum-groups 90 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --target-bits 52 --max-cells 400 --max-depth 22 --output tmp/gf16-feedback-d05-q90-rational-full.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --replay tmp/gf16-feedback-d05-q90-rational-full.json --precision 256
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-feedback-d05-q90-rational-full.json --precision 256 --sparse-inner gf-shape-tilt --shape-penalties 1 .75 --sparse-tilts .00032 .0032 .032 .064 .08 .096 --sparse-log tmp/gf16-d05-shape-full-replay.log --output tmp/gf16-d05-shape-complete.json
```

The separate sparse and dense commands do not cover omitted occupancies.
The final command reproduces the complete 5% result.
`--operator-through`
may request a longer placement prefix for local memo reuse; it does not
change the requested proof domain. The existing exact prefix checks still
run. The final assembler never relies on memoized operators.

The full 6.5%, 7%, and 7.5% results are reproduced by:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 13/200 --kernel gf-refresh --base-tilt 1/8 --minimum-groups 65 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 250 --max-depth 22 --output tmp/gf16-gf-refresh-d065-q65.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-gf-refresh-d065-q65.json --precision 256 --sparse-inner gf-occupancy --exact-feedback --sparse-tilts .00032 .0032 .016 .032 .064 --sparse-log tmp/gf16-d065-exact-full-replay.log --output tmp/gf16-d065-exact-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 7/100 --kernel weighted-return --base-tilt 1/8 --minimum-groups 65 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 250 --max-depth 22 --output tmp/gf16-weighted-return-d07-q65.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-weighted-return-d07-q65.json --precision 256 --sparse-inner gf-occupancy --exact-feedback --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-d07-exact-full-replay.log --output tmp/gf16-d07-exact-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 3/40 --kernel profile-return --base-tilt 1/8 --minimum-groups 65 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 300 --max-depth 22 --output tmp/gf16-profile-return-d075-q65.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-profile-return-d075-q65.json --precision 256 --sparse-inner gf-occupancy --exact-feedback --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-d075-profile-full-replay.log --output tmp/gf16-d075-profile-complete.json
```

The two-update 9% replay and four-update 10% replay are complete.
Next, optimize the certified implementation without changing its setup
distribution. The commands below record the earlier proof searches,
not a queue of active jobs. The fixed-occupancy
`return_moment.py` census and lazy-density bound are now optional inputs
to the sparse kernel, freshly regenerated on each run. No certificate imports the old lane-permutation
or two-bit results.
Selected-point checks and partial occupancy runs are not interpolated.
Keep the complete two-bit 9.25% certificate as the fallback.

Earlier search commands are:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 31/400 --kernel conditioned-return --base-tilt 1/8 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 400 --max-depth 22 --output tmp/gf16-conditioned-return-d0775-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/sparse_cover.py --inner gf-occupancy --exact-feedback --through 96 --threshold 162529 --tilts .00032 .0032 .016 .032 .064 .096 --precision 256 --output tmp/gf16-d0775-sparse-through96.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-conditioned-return-d0775-q97-v3.json --precision 256 --sparse-inner gf-occupancy --exact-feedback --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-d0775-conditioned-full-replay.log --output tmp/gf16-d0775-conditioned-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 2/25 --kernel rank-return --base-tilt 1/8 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 400 --max-depth 22 --output tmp/gf16-rank-return-d08-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-rank-return-d08-q97.json --precision 256 --sparse-inner gf-rank --exact-feedback --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-d08-rank-full-replay.log --output tmp/gf16-d08-rank-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 41/500 --kernel rank-return --row-parity --base-tilt 1/8 --minimum-groups 105 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 400 --max-depth 22 --output tmp/gf16-parity-rank-d082-q105.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-parity-rank-d082-q105.json --precision 256 --sparse-inner gf-rank --exact-feedback --single-group-exact --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-parity-rank-d082-full-replay.log --output tmp/gf16-parity-rank-d082-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/single_group.py --output tmp/gf16-rank-d10-single-exact.json
```

The next dense refinement is reproduced with:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 17/200 --kernel birth-refresh --row-parity --base-tilt 1/8 --minimum-groups 129 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --output tmp/gf16-birth-d085-q129.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 17/200 --kernel birth-refresh --row-parity --variance-shuffle --base-tilt 1/8 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --output tmp/gf16-birth-variance-d085-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-birth-variance-d085-q97.json --precision 256 --sparse-inner gf-rank --exact-feedback --single-group-exact --sparse-tilts .00032 .0032 .016 .032 .064 .096 --sparse-log tmp/gf16-birth-variance-d085-full-replay.log --output tmp/gf16-birth-variance-d085-complete.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 43/500 --kernel birth-refresh --row-parity --variance-shuffle --variance-bins 8 --base-tilt 1/8 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --output tmp/gf16-birth-variance-partition-d086-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 87/1000 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 8 --base-tilt 1/8 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --output tmp/gf16-birth-classes-variance-d087-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/sparse_cover.py --occupancies 48 --threshold 182452 --tilts .00032 .0032 .008 .016 .024 .032 .048 .064 .096 .128 .192 --precision 384 --max-splits 150 --target-bits 52 --inner gf-birth-classes --exact-feedback --output tmp/gf16-birth-classes-d087-q48-fulltilt-p384.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-birth-classes-variance-d087-q97.json --precision 256 --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .048 .064 .096 .128 .192 --sparse-log tmp/gf16-birth-classes-d087-full-replay.log --output tmp/gf16-birth-classes-d087-complete.json
```

`base_tilt_probe.py` separately scans the proof's reference tilt at 9%,
using the birth-class and variance-partition bounds. This changes only
proof parameters. Its floating point samples do not cover all means and
are not certificates; the scan guides the next whole-domain search.

The base tilt 3/16 improves the sampled 9% frontier over the previous
1/8 choice. Four narrow mean intervals corresponding to reference
activities 1/8, 3/20, 7/40, and 1/5 have been evaluated outward at both
256 and 384 bits; their log2 bounds are approximately -642.35, -184.37,
-129.05, and -203.77. The complete dense search using that reference tilt
closes with 116 accepted cells and no unresolved cells:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 9/100 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 8 --base-tilt 3/16 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --output tmp/gf16-birth-classes-d09-tilt3over16-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-birth-classes-d09-tilt3over16-q97.json --precision 256 --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 .192 --sparse-log tmp/gf16-birth-classes-d09-full-replay.log --output tmp/gf16-birth-classes-d09-complete.json
```

The original 9% assembly's coarse grid failed at q=3 and q=4. That run
was stopped with logs preserved. A finer grid verifies all supports at
q=4 with 89.2548 bits at 384-bit precision. The following fresh assembly
uses a separate low-occupancy grid and a 48-bit per-occupancy stopping
budget. The final aggregate must still meet the same 40-bit threshold;
changing that search budget does not waive any occupancy or proof check.

The older 8.5% assembly also finished with only q=4 unresolved; all
other q=2--96 checks passed, with combined log2 upper -50.8264017078.
Its coarse grid missed the 52-bit per-occupancy stopping target there.
The 8.6% replay found the same q=4 gap and was stopped, preserving its
logs. Neither is a complete certificate. The fresh 9% run below includes
the corrected low-occupancy grid and supersedes both attempts.

```text
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-birth-classes-d09-tilt3over16-q97.json --precision 256 --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-target-bits 48 --sparse-small-through 8 --sparse-small-tilts .00016 .00032 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 .192 --sparse-log tmp/gf16-birth-classes-d09-smallgrid-full-replay.log --output tmp/gf16-birth-classes-d09-smallgrid-complete.json
```

This full assembly completed. It independently replays every dense
cell and regenerates the entire sparse prefix; it does not import
selected-occupancy receipts. Its dense contribution has log2 upper
-78.7739845636; the complete margin is 46.8529233246 bits, as stated
above. Its exact q=1 calculation gives 51.4257 bits, and its q=3 and
q=4 checks give 59.4822 and 89.2548 bits. Those individual margins are
not substitutes for the complete sum.
Separate 256-bit checks cover all supports at q=32, 48, 64, and 96,
giving 146.7917, 56.2764, 54.1585, and 71.7844 bits respectively.
The strongest complete certificate remains the one stated above.

Further reference-tilt diagnostics at 9.5% still leave positive scores.
The additional base tilts 5/32 and 11/64 do not resolve the sampled
activities 0.125--0.30. Together with the earlier 3/16--5/16 scan, this
has not found a simple retuning route to 9.5%. It does not establish an
encoder limitation or exclude other proof parameters.

The next refinement keeps the regional packet count rather than using
one worst-case shuffle-density factor for every count. Its argument and
scope are in [the shuffle note](SHUFFLE_VARIANCE.md#retaining-the-regional-packet-count).
Three complete mean intervals at 9.25% have negative outward bounds,
agreeing at 256 and 384 bits. The full dense search completed in 172
cells, with no unresolved cells:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 37/400 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 8 --regional-count --base-tilt 3/16 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-count-d0925-q97.json
```

An independent 384-bit replay of that dense cover passed, with aggregate
log2 upper bound -76.4285156748 and no unresolved cells. This covers
q=97--2048 only; it is not the whole-code margin. The replay log is
`tmp/gf16-regional-count-d0925-replay-p384.log`.
The option is stored in dense metadata and checked by full assembly.
Replay recomputes the ordered polynomial matrices and every saved MGF
majorant, with no reliance on a numerical LP objective. At the same
9.25% cutoff, selected all-support sparse checks give 52.5499 bits at
q=48 and 52.2985 bits at q=64. The initial q=96 search exhausted 200
subdivisions with floating log2 sum +64.5513, so it is not certified.
A retry adds finer output tilts around 0.096 and allows 300 subdivisions.
That finer-grid retry now covers every support at q=96, with
50.2158000153 bits at 256-bit precision. It used 247 subdivisions and
17569 retained leaves. Its receipt is `tmp/gf16-r2-d0925-q96-fine.json`,
SHA256 `bff57852a3546b65f58cb60d0d70b22e6f11bf37f06615326f2df7ef3af91548`.
This clears the tested boundary occupancy, not all earlier occupancies.
The complete sparse prefix remains a separate obligation; there is no
whole-code 9.25% certificate yet. The serial full assembly reached q=73
but is no longer running and produced no complete receipt. A fresh
four-worker assembly now regenerates the entire prefix:

```text
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-regional-count-d0925-q97.json --precision 256 --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 350 --sparse-target-bits 48 --sparse-small-through 8 --sparse-small-tilts .00016 .00024 .00028 .00032 .0004 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .088 .092 .096 .10 .104 .108 .112 .12 .128 .16 .192 --sparse-analytic-gradient --sparse-workers 4 --sparse-log tmp/gf16-r2-d0925-parallel-sparse.log --output tmp/gf16-r2-d0925-parallel-complete.json
```

A further selected-interval refinement retains anti-concentration after
the regional count's exponential tilt. At 9.5%, 32 variance parts and
a finer count-tilt grid give log2 upper bounds -107.7590 and -9.0172
near reference activities 0.15 and 0.20. The latter still misses the
40-bit target. These fresh 256-bit results are partial. The completed
9.25% dense search did not use this
extra refinement. See the shuffle note for the tilted-count identity,
variance lower bound, and verifier interpretation.
The 64-part check at the intervening reference activity 0.175 remains
positive (+50.5177). The 9.5% frontier is therefore still unresolved,
even before its complete-domain and sparse-prefix obligations.
Output-tilt retuning and a checked finite-type variance minorant improve
that interval to +37.3798, still insufficient. Reoptimizing the global
minimum-activity penalty also failed to help in a separate eight-part
diagnostic. These tests do not show an actual low-weight codeword.

The selected-interval tool permits `--updates 3` as a construction
comparison. This changes the actual inner, unlike the other proof
options. At 9.5%, the mean intervals around reference activities 0.175
and 0.20 have fresh 256-bit log2 upper bounds -1349.0381 and -1588.2668,
using 16 variance parts and the refined regional count bound. At 10%,
the intervals at activities 0.15, 0.175, and 0.20 give -43.5945,
+224.1512, and +250.0970. Thus 10% is still unresolved.

The initial three-update 9.9% dense search used:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 99/1000 --updates 3 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 16 --regional-count --base-tilt 3/16 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-count-r3-d099-q97.json
```

It was stopped after the 300-proposal checkpoint to use the stronger
direct count-mass bound. That checkpoint retains 65 accepted cells and
171 unresolved cells. The continuation rechecks the retained witnesses
before proposing direct count bounds; the original record is preserved:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-regional-count-r3-d099-q97.json --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-direct-r3-d099-q97.json
```

The direct bound already verifies two intervals of radius 1/10000 at
9.9%, centered at tilted means 7/183 and 3/67. Their 256-bit log2 upper
bounds are -70.9996 and -99.0942. Independent 384-bit evaluation agrees
to the displayed digits. These selected intervals do not establish complete coverage.
See [the direct-count argument](SHUFFLE_VARIANCE.md#bound-count-masses-directly).

The q>=97 continuation later stopped on floating regional underflow,
after its 500-proposal checkpoint retained 72 accepted and 193 unresolved
cells. Its checkpoint is preserved, but it is not complete. The proposal
routine now normalizes conditional matrices separately at every count
and placement step. Regression tests include products below binary64's
range and structurally zero counts. Outward replay is unchanged.

The full assembler takes the update count from the dense record and
regenerates its exact q=1 and remaining sparse bounds with that same
count. Three-update sparse runs require the GF rank or birth-class
operator and exact feedback; legacy engines and operator caches are
rejected. Run receipts record the update count explicitly. Tests check
forwarding through all three paths. A separate exact test enumerates the
actual transvection sampler in dimensions two through four and composes
its transition matrix. It verifies the lazy/uniform mixture for one,
two, and three updates without assuming that mixture in the test's input.
The paired implementation timing is
reported above, but no full three-update distance claim is made yet.
The exact q=1 calculation at 10% gives 46.4988994496 bits of margin;
independent 256- and 384-bit evaluations agree. All-support checks at
q=2,3,4,8 also close at 10%, with their combined log2 upper bound
-69.0763849890. These selected occupancies do not substitute for the
entire sparse prefix. The 9.9% check at q=48 now covers every support,
with 54.3053638844 bits of margin at 256-bit precision. The first q=64
attempt exhausted 200 subdivisions with a floating log2 sum of +83.2092,
so it produced no outward certificate. Its largest remaining boxes
include supports 176--202, using output tilt 0.064. A focused retry adds
tilts at spacing 0.004 around that value and allows 300 subdivisions.
That retry did not close. Its restricted grid omitted the original
low-tilt choices, and its worst remaining box contains supports 38--41
at the smallest available tilt. Its log2 score +131.3443 therefore
does not establish a limit of the combined grid. A corrected retry
uses the union of both grids and allows 400 subdivisions. It finished
without a certificate: 14870 retained leaves, with floating log2 sum
+32.8086. The largest remaining leaves include both mixed broad-support
boxes and boxes concentrated around supports 176--202. Thus the union
grid fixes the earlier omission but does not by itself close this bound.
The q=96
check also exhausted its 200-subdivision budget, with floating log2
sum +627.6369 and no outward certificate. The full dense search remains
running. Full 9.9% closure still
requires covering every omitted occupancy, either sparsely or through
an earlier dense handoff.

The q=64 retry uses the freshly checked joint-return and lazy-density
refinements, with analytic gradients for witness proposals. It closes
after one subdivision, with 65 retained leaves and outward margin
117.4341324137 bits. Its receipt SHA256 is
`45f56e3b56fc65165dfc9a825b84014a7ec55956695befa1ad22ba9383eeec45`.
This is a complete support bound for q=64, not the full code. A fresh
384-bit replay reproduces 117.4341324137 bits. The same invocation also
finishes q=96, with margin 56.0981821225 bits at 9.9% distance. This is
a fresh 384-bit result for q=96, not a replay of an earlier successful
q=96 certificate. The receipt is
`tmp/gf16-r3-d099-q64-q96-refined-p384.json`, SHA256
`9202583efbb868a5797f3ab64a21c07572fa0f8b7a5815f8f5ef692503424597`.

```text
python -B research/workstreams/permutation_locality/gf16_packets/sparse_cover.py --updates 3 --occupancies 64 --threshold 207618 --tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .052 .056 .06 .064 .068 .072 .076 .08 .088 .096 .128 .192 --precision 256 --max-splits 200 --target-bits 48 --inner gf-birth-classes --exact-feedback --joint-return-through 3 --lazy-density-through 6 --analytic-gradient --output tmp/gf16-r3-d099-q64-joint-lazy-analytic.json
```

The sparse/dense handoff is a proof parameter, not an encoder parameter.
The selected-interval tool now accepts `--minimum-groups` to test an
earlier handoff. For q>=49 and the same three-update 9.9% construction,
direct count bounds give fresh 256-bit log2 upper bounds -139.6757,
-84.0904, and -73.4606 at reference activities 0.10, 0.125, and 0.15.
These cells have radius 1/1000000 in tilted-mean coordinates. Other
sampled activities from 0.0151 through 0.08 also have negative floating
proposals, as do 0.175 and 0.20. The complete q>=49 search is now running:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 99/1000 --updates 3 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 16 --regional-count --base-tilt 3/16 --minimum-groups 49 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-direct-r3-d099-q49.json
```

Its per-cell stopping target is 60 bits. The final verifier still sums
every accepted cell and requires the whole-code aggregate to meet
40 bits; no cell or occupancy is waived. The search does not yet
discharge every remaining occupancy. A successful earlier handoff must
still be combined with a freshly regenerated sparse prefix q=1--48.
The direct-count run was stopped at its 140-proposal checkpoint, with
21 accepted and 99 unresolved cells, to enable the new lazy-density
and exact-return alternatives. The continuation retains and freshly
rechecks the old witnesses; no accepted cell is assumed correct:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-regional-direct-r3-d099-q49.json --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-joint-lazy-r3-d099-q49.json
```

That prefix completed at 256-bit precision, using exact q=1 placement,
the fine grid for q=2--8, and the full grid for q=9--48. Every occupancy
passes. Their complete sum is at most 7.668e-15, giving 46.8900721773
bits for q=1--48. The log is `tmp/gf16-r3-d099-q1-48-prefix.log`.
This prefix does not assert a whole-code bound before the dense range
is complete.

For a continuation with unchanged bounds, `--keep-unresolved-partition`
preserves previous subdivisions instead of coalescing unresolved siblings.
This avoids repeating coarse failed proposals. Every retained certificate
is still freshly checked, and the verifier reconstructs the complete
partition. Coalescing remains the default when trying a stronger bound
that may close wider intervals.
After the floating underflow fix, the main q>=49 search was resumed from
its 260-proposal checkpoint, retaining its 30 accepted and 107 unresolved
cells without coalescing:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-regional-joint-lazy-r3-d099-q49.json --keep-unresolved-partition --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-joint-lazy-stable-r3-d099-q49.json
```

A further three-update diagnostic uses 64 variance parts and varies
the output tilt at 10%. Near reference activities 0.175 and 0.20, the
best scores among scales 0.95, 0.975, 1, and 1.025 remain positive
(+162.1393 and +182.0288, respectively).
This is a floating search result, not an outward certificate or a
counterexample. It supports prioritizing complete 9.9% coverage before
another attempt to remove the last gap at 10%.

Two follow-up diagnostics narrow the remaining proof work. A conditional
birth-density bound passes exact local tests but gives at most about six
bits of iid improvement in the tested screen, so it is not integrated.
The regional counterfactual instead points to return-to-zero bounds and
the loss of distribution information after lazy updates. Exact joint
counts through three active packets improve the two selected 10% scores
from +162.1393/+182.0288 to +74.2901/+104.2342. These remain floating,
selected-interval results, not certificates. The four-packet census has
finished its independent integer checks and further improves these
scores to +63.8535/+90.3329; neither closes. All 159 local tests pass. See
[the diagnosis](BIRTH_DENSITY.md#three-update-regional-diagnosis-at-10)
and [the local return calculation](OCCUPANCY_BOUND.md#cancellation-must-still-emit-weight).

A new refinement retains the existing uniform-density coordinate through
selected lazy updates and retains expansion classes on empty lazy steps.
Combined with exact return counts through three active packets, it
verifies the two troublesome 10% mean intervals at 256-bit precision,
with log2 upper bounds -220.1615 and -68.7218. Independent 384-bit
evaluation agrees. The scalar verifier can now replay both flags
and proposes the new bound only as an alternative to the old one.
No whole-code certificate uses these refinements yet, and encoder
operations are unchanged. See
[the density argument](BIRTH_DENSITY.md#retaining-density-through-lazy-updates).

The same refinements now verify three two-update 9.5% mean intervals.
At reference activities 0.15, 0.175, and 0.20, fresh 256-bit log2 upper
bounds are -886.7217651579, -498.2153608770, and -381.6294321805. Each
interval has radius 1/1000000 in tilted-mean coordinates. They cover
q>=97, using direct count masses, 64 variance parts, and output-tilt
scale 39/40. The local receipt is
`tmp/gf16-r2-d095-joint-lazy-selected-p256.json`, SHA256
`437e8578530079c1784bfcd788ab02863349dfb44847ca8dc5221e3af7904a3e`.
An independent 384-bit evaluation agrees to the displayed digits.
Its local receipt is `tmp/gf16-r2-d095-joint-lazy-selected-p384.json`,
SHA256 `fe95b77571f3e63e7495c214cfa20d72af41d799052eb9d61d76d699d8be0229`.
These three intervals do not constitute a full dense cover. A complete
two-update 9.5% dense search has started with 16 variance parts:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 19/200 --updates 2 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 16 --regional-count --base-tilt 3/16 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-joint-lazy-r2-d095-q97.json
```

The sparse verifier now accepts `--joint-return-through 3` and
`--lazy-density-through 6` with `--inner gf-birth-classes --exact-feedback`.
It regenerates the integer return census once per invocation and checks
the one-packet density census at each output tilt. The assembler forwards
the corresponding `--sparse-joint-return-through` and
`--sparse-lazy-density-through` choices and records them in its receipt.
It still regenerates every sparse occupancy; saved run receipts are not
proof inputs. Omitting the flags preserves the original operators.
Exact small-state tests also check ordered products with both local
refinements applied together, including degenerate feedback maps.

The optional `--analytic-gradient` sparse search (or
`--sparse-analytic-gradient` during assembly) uses `cdf_gradient.py`
to differentiate grouped count probabilities, matrix powers, and CDF
folds. Quantized candidates are still evaluated with the original
floating objective, then checked outward if retained. Gradients never
serve as proof inputs. Tiny-moment, finite-difference, and end-to-end
cover tests pass. On three selected actual q=96 boxes at output tilt
0.096, both optimizers find the same values to numerical precision;
analytic derivatives reduce objective evaluations from 16/63/24 to
8/20/12. This is not an end-to-end timing measurement or a certificate.
The option defaults off, preserving the existing controls.

The remaining dense checkpoints concentrate on reference activities
about 0.25--0.47 for the three-update 9.9% search, and 0.25--0.37 for
the two-update 9.5% search. They have not reached the maximum subdivision
depth. Those ranges identify unfinished work, not a construction limit.
A narrow-cell diagnostic at activity 0.30 gives more specific evidence.
It covers q>=49 and the tilted-mean interval 9/121 +/- 1/1000000.
With three updates, 9.9% distance, joint returns through three packets,
and lazy density through six, its log2 bound is +290.1241 at output-tilt
scale 39/40. Increasing variance parts from 16 to 64 lowers this to
+248.1762 but does not close the cell. Tested output-tilt scales
9/10, 4/5, and 7/10 also fail to close it.

At the same fixed witness, retaining lazy density through four packets
gives +260.4771; extending through eight gives +981.7210. Thus extending
the existing density allocation is not automatically an improvement.
Deleting lazy contributions gives negative diagnostic values, but those
deletions are not valid bounds on the encoder. These are floating
diagnostics, not certificates or an optimized lower limit on the proof.

The updated `regional_frontier_probe.py` reconstructs direct count masses
and saved local refinements before making these comparisons. It checks
that its baseline matches the saved proposal. Alternative allocations
start from unrefined operators, avoiding a second application to a U
entry that already contains lazy mass. The ambiguous refresh-only-return
ablation is omitted when that entry includes lazy density.

[The feedback-density prototype](BIRTH_DENSITY.md#density-supplied-by-the-gf-feedback)
tests whether the GF feedback can itself spread an arbitrary entering
state. Its small exact tests pass, but its initial actual-size diagnostic
does not improve the saved bound: the best tested score is +473.5934,
versus +290.1241 before replacement. That density-only allocation is not
enabled in the certificate verifier.

Retaining both constraints instead gives a useful refinement. For each
outgoing expansion-weight class, its lazy mass is at most the smaller
of the total lazy mass and the class size times the pointwise density
cap. The two bounds concern the same class mass, so their minimum is
valid. This keeps the existing F-class invariant and comparison-matrix
interface. Exact small-state tests check each outgoing class and
multistep positive terminal costs, including zero feedback.

At activity 0.30, replacing the arbitrary lazy allocation by these class
bounds for packet counts 3--32 gives floating log2 score -398.7713.
A fresh 384-bit outward calculation verifies -398.7713020749524 over
the entire interval 9/121 +/- 1/1000000, for q>=49 at 9.9% distance.
The receipt is `tmp/gf16-r3-d099-feedback-classes-p03-p384.json`, SHA256
`a2c2dfcbee6c4a2356609701f39bc1c86e9231e32fb9d00529693667e53a4435`.
The verifier rebuilds these bounds; the scalar proposer selects them
only when better than its earlier alternatives. A separate fresh
256-bit evaluation agrees to the displayed digits. This closes one
mean interval, not the full dense domain or the whole-code proof.
At activity 0.35, the best of output-tilt scales 39/40, 9/10, and 4/5
still has floating log2 score +517.2243. An optional diagnostic extends
the class-mass bound to lazy transitions from U, using exact profile
multiplicities to average the fiber cap. It preserves any lazy branch
already allocated by the earlier density refinement. Its exact multistep
tests pass. At activity 0.35, it improves the best tested fixed-witness
score from +517.2243 to -9.1928. This is still a diagnostic and falls
short of the required 40-bit margin. Replacing complete U rows before
their earlier density allocation improves the score to -74.4654.
Fresh 384-bit outward evaluation verifies -74.4654194057170 over
21/229 +/- 1/1000000, including all 16 variance parts and q>=49 at
9.9% distance. The receipt is
`tmp/gf16-r3-d099-feedback-uniform-p035-p384.json`, SHA256
`dbc5c5ffd02895c24e3bacbceef724ce0955d7fa9514cd0e26ef4a3e061883ab`.
A fresh 256-bit check agrees to the displayed digits. These U-row options are now available
to the verifier and scalar proposer, without changing old witnesses.

The activity-0.40 diagnostic still fails (+2980.9415). The
[shared-budget argument](MASS_DENSITY.md) addresses another loss:
the existing class allocations can collectively exceed the branch's
total mass. Its nonlinear backward bound retains that budget and
passes exact multistep tests. Initial iid screens show substantial
gains, but the actual regional improvement is smaller: the activity-0.40
score becomes +2449.4363, and the activity-0.45 score +10970.7775.
Both still fail. Increasing the variance partition to 64 parts improves
these scores only to +2261.5565 and +10681.9476. Direct count distributions
for representative integer compositions also leave positive bounds.
Neither diagnostic proves an actual bad codeword. They motivate
sharpening local transitions before implementing a nonlinear regional
verifier.

The [exact-zero refinement](ZERO_TRANSITION.md) addresses a newly isolated
loss: the old zero-to-zero entry did not use its weighted Fourier count.
The new optional `--exact-zero` replaces that entry and leaves the
encoder unchanged. At three updates and 9.9% distance, the activity-0.40
interval now has outward log2 bound -177.7997549331, using output-tilt
scale 9/10 and all 16 variance parts. The fresh 384-bit receipt is
`tmp/gf16-r3-d099-exact-zero-p04-p384.json`, SHA256
`b3d93ddc0e6466e0569730ffc886d29dc87a68b5cfd9f30b4a264c50c5982993`.
It covers 1/9 +/- 1/1000000 and q>=49. An independent 256-bit check
agrees to the displayed digits. The activity-0.45 interval still fails at tested scales
4/5, 9/10, 19/20, 1, and 21/20; the best score is +2255.3303.
Combining the shared-mass refinement with the exact zero entry improves
that score to +1897.8712, still insufficient. Extra variance subdivision,
count-envelope replacement, or this shared budget alone has not closed
the remaining three-update gap. This does not establish a distance
limit for the construction.
This is not a complete 9.9% certificate.
All 192 local tests pass, including exact zero-transition counts,
feedback-class partition sums, and multistep checks. Performance is
unchanged because these edits only affect proof calculations.

```text
python -B research/workstreams/permutation_locality/gf16_packets/regional_frontier_probe.py tmp/gf16-r3-d099-high-activity-screen.json --activities .3 --precision 256 --joint-only --feedback-density-starts 2 3 4 6 8 --feedback-density-through 32 --output tmp/gf16-r3-d099-feedback-density-screen.json
python -B research/workstreams/permutation_locality/gf16_packets/regional_count.py --distance 99/1000 --updates 3 --minimum-groups 49 --activities .3 --radius 1/1000000 --variance-bins 16 --tilted-atom --fine-tilts --tilted-variance --direct-counts --tilt-scales 39/40 --lazy-density-through 6 --joint-return-through 3 --feedback-classes-from 3 --precision 384 --output tmp/gf16-r3-d099-feedback-classes-p03-p384.json
```

A bounded q=96 check tests whether these refinements can close the
two-update middle-occupancy gap at 9.5%:

```text
python -B research/workstreams/permutation_locality/gf16_packets/sparse_cover.py --updates 2 --occupancies 96 --threshold 199229 --tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 .192 --precision 256 --max-splits 200 --target-bits 48 --inner gf-birth-classes --exact-feedback --joint-return-through 3 --lazy-density-through 6 --output tmp/gf16-r2-d095-q96-joint-lazy.json
```

This run exhausted 200 subdivisions with floating log2 sum -4.3252;
it produced no outward certificate. The retry uses the finer tilt grid
from the 9.25% q=96 closure, analytic gradients, and 350 subdivisions.
It now verifies every support at q=96, with margin 50.3468317586 bits
at 9.5% distance and 256-bit precision. Its receipt is
`tmp/gf16-r2-d095-q96-refined-analytic.json`, SHA256
`0ceb24e40733341e9f361e647c8a886762396395ffd35c11e7ed5f06b02ce135`.
This is not a complete sparse prefix or a whole-code certificate.
For the same two-update 9.5% construction, the class, U-row, and exact
zero bounds now verify log2 upper bounds -2048.9354059732 and
-1979.6028057509 at activities 0.30 and 0.35, respectively. Each result
covers a mean interval of radius 1/1000000, every variance part, and
q>=97. The fresh 384-bit receipt is
`tmp/gf16-r2-d095-exact-zero-p03-p035-p384.json`, SHA256
`1cafa323cc38af87ef7cf5f19c388b927c4f77dd0f45543d97213829a2a59cda`.

### Full-Cover Searches with the Exact Zero Entry

The scalar proposer now tries three rational rescalings of the output
tilt if the strongest regional allocation still fails. It keeps all
earlier candidates. Only fresh outward evaluation can accept a cell.
A new two-update 9.5% cover starts from the entire domain, rather than
depending on the selected intervals above:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 19/200 --updates 2 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 16 --regional-count --base-tilt 3/16 --minimum-groups 97 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 500 --max-depth 22 --target-bits 60 --precision 256 --output tmp/gf16-exact-zero-r2-d095-q97.json
```

Four updates are a separate construction comparison. At 9.9%, fresh
256- and 384-bit calculations agree on log2 upper bound
-423.6688850297 for 27/203 +/- 1/1000000, all 16 variance parts, and
q>=49. This is the previously failing activity-0.45 interval. The
384-bit receipt is `tmp/gf16-r4-d099-exact-zero-p045-p384.json`, SHA256
`c1b08c9a9e07282d904468f771e41522e17d047880ffb9ca64fc8f15d5eaa6f3`.
Exact q=1 placement separately gives 49.6232469487 bits at 9.9%, at
256-bit precision. Neither result is a whole-code certificate.

The dense, sparse, and aggregate tools now forward and check four
updates explicitly. Tests enumerate the sampled transvection kernel
through eight successive updates and exercise four-update class,
return, density, and aggregate dispatch paths. Old update counts and
defaults remain unchanged. The full four-update dense search is:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --distance 99/1000 --updates 4 --kernel birth-classes --row-parity --variance-shuffle --variance-bins 16 --regional-count --base-tilt 3/16 --minimum-groups 49 --central-bits 129 --central-scale 1001/1000 --row-bias 21/50 --max-cells 800 --max-depth 22 --target-bits 60 --precision 256 --output tmp/gf16-exact-zero-r4-d099-q49.json
```

At that checkpoint, both dense searches were incomplete. The four-update
q=1--48 prefix subsequently completed at 256-bit precision. It regenerates exact q=1 placement
and all support covers for q=2--48, using fresh birth-class operators
and exact feedback. No three-update bound is reused. At cutoff 207618,
the prefix's total bad-setup contribution is at most 1.166e-15, giving
49.6084087286 bits of margin. Its local receipt is
`tmp/gf16-r4-d099-q1-48-prefix.json`, SHA256
`39af04e194fcbf9c1e104e0e26bcd1c8dd07d6f694c9d0bcab63859e97129994`.
An independent 384-bit regeneration reproduces the margin and gives
an outward upper endpoint no larger than the 256-bit result. Its
receipt is `tmp/gf16-r4-d099-q1-48-prefix-p384.json`, SHA256
`d0d4d3212122dd16a796574251910f56695907f92b0bc2ba6053718bb21e45b4`.
The dense occupancies are q=49--2048. All 427 dense intervals have since
passed fresh 256-bit checks, and the full 384-bit assembly has completed.
The complete 9.9% certificate above has 49.6084 bits of margin and
corresponds to the paired four-update timing of 7.483663 ms.

At exactly 10%, fresh 256- and 384-bit calculations agree on a bound
of **48.6550601417 bits of margin** for the complete q=1 contribution.
The 256-bit calculation is the first stage of the new sparse-prefix run.
Both calculations use four updates
and cutoff 209715, including every support of one active group. Thus
q=1 does not by itself prevent a 40-bit certificate at 10%. This is
an upper bound on that class's bad-setup contribution, not an estimate
of the code's typical distance. The receipt is
`tmp/gf16-r4-d10-q1-p384.json`, SHA256
`385eb8162d808875f7b840c3cc218b963d9bee2ad6db6e646384959f5cfed168`.
The command is:

```text
python -B research/workstreams/permutation_locality/gf16_packets/single_group.py --threshold 209715 --updates 4 --precision 384 --tilts .00016 .00024 .00028 .00032 .0004 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 .032 .064 --output tmp/gf16-r4-d10-q1-p384.json
```

The fresh q=1--48 run at 10% has completed at 256-bit precision. Its
aggregate upper bound is 2.26398771182729e-15, giving
**48.6500552956 bits of margin**. All 48 occupancies include their
complete support ranges. The receipt is
`tmp/gf16-r4-d10-q1-48-prefix.json`, SHA256
`cb925f497d86f74e4a897e65b85d0321d7807bb1e119223752fb60e47506b10a`.
An independent 384-bit regeneration reproduces the displayed margin.
Exact dyadic comparison confirms its upper endpoint is no larger than
the 256-bit endpoint, and both are below 2^-48. The second receipt is
`tmp/gf16-r4-d10-q1-48-prefix-p384.json`, SHA256
`959f60e4b4980e3702e418abcdcde0d7ec4ba1c8cc50081f6340e9e65b8879ee`.
The complete 9.9% assembly passed first. Retargeting its dense partition
and regenerating the sparse prefix subsequently closed 10%, as stated
in the complete result above.

A diagnostic retargeting of its first 80 freshly checked intervals
suggests that their existing witnesses also suffice at 10%. The three
intervals with the least predicted slack were then checked from scratch
at 384-bit precision. Their margins at cutoff 209715 are 78.3944773615,
204.6995299027, and 267.2077027414 bits. These checks cover whole saved
mean intervals for q>=49, not isolated points. Their paths are
`0000000010`, `0000000001`, and `00110101`. The local receipt is
`tmp/gf16-r4-d10-frontier-p384.json`, SHA256
`39af3196757b95ff8a1cdd1c6fd4c5e90288facbe11bf2c60c4644c8786571e8`.
These three selected checks did not complete the dense cover. All 427
intervals have since passed fresh 256-bit checks, followed by the
complete 384-bit assembly and exact audit recorded above.

The saved q=1 bound also identifies a concrete tightening target. In its
CDF fold, the support-38 term contributes about 98.4557% of the bound;
support 40 contributes another 1.5081%. Supports 38, 40, and 42 together
contribute about 99.9974%. These are contributions to the CDF-majorant
calculation, not measured shell multiplicities or observed failures.
Thus the current sparse bound is limited mainly by its treatment of
the lightest supports, rather than by large occupancies. Sharpening
their counting or output-tail bounds is a possible next step after
the dense cover closes. This attribution does not establish the code's
actual minimum distance or the limiting part of a complete certificate.

### Reuse Search Work, Not Certificate Inputs

Profiling found repeated computation of the same one-packet census and
regional count caps while testing different local envelopes. The
proposer now keeps bounded in-memory scratch for the current output
tilt and count-witness family. Count weights are reused across output
tilts only when the cell, variance partition, duals, feature family,
geometry, and working precision agree. The local census also requires
the same map instance and output tilt. Changing precision or map
instance discards the scratch.

The outward verifier does not read this numerical scratch. It rebuilds
the local bounds and count caps from the saved rational witness. Tests
check cache invalidation, unchanged proposal values, exact small-instance
domination, and fresh outward regeneration after a cached proposal.
All 192 tests pass. One repeated four-update proposal gives the same
score before and after reuse and takes about 6.3 seconds under profiling;
its cold evaluation takes about 46 seconds, including census setup.
This is a search diagnostic, not an encoder speedup or a full-search
timing claim.

The searches continue from preserved partitions with this optimization.
Every retained leaf is rechecked before further proposals. The
two-update resume starts from 36 accepted and 109 unresolved cells;
the four-update resume starts from 34 accepted and 73 unresolved cells.
Their new outputs are `tmp/gf16-exact-zero-cached-r2-d095-q97.json` and
`tmp/gf16-exact-zero-cached-r4-d099-q49.json`. The independent complete
9% replay and the 9.25% assembly continue without interruption.

The newer regional envelope does not dominate every older bound.
For example, at four updates and 9.9%, its activity-0.50 diagnostic
remains positive, but the retained variance-partition witness certifies
the broader containing mean interval with log2 upper about -808.18.
The same older argument certifies the interval containing activity
0.60 with about -770.38. Both leaves pass the fresh resume check.
These points are therefore not unresolved obstructions. The search
keeps all earlier candidates and concentrates further refinement on
its actual uncovered intervals.

### Continue from the Most Advanced Partition

The older two-update 9.5% partition contains 154 accepted cells and
213 unresolved cells. Its accepted intervals contain every interval
accepted by the newer root-started run at the checked checkpoint.
The three-update 9.9% partition contains 200 accepted cells and 227
unresolved cells. These more advanced partitions now seed the main
two- and four-update searches. All earlier checkpoints are preserved;
the redundant live searches were stopped.

`--retarget-updates 4` changes the model used to evaluate the saved
rational witnesses. It does not import three-update numerical bounds.
Every old leaf is evaluated outward for four updates; a leaf that fails
the target returns to the unresolved queue. The same procedure supports
an explicit `--retarget-distance`. Retargeting requires `--resume` and
a different output file; ordinary replay cannot change its claim.
The output records the source hash and the old claim for provenance.

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-r2-d095-advanced-start.json --keep-unresolved-partition --target-bits 60 --max-cells 800 --max-depth 22 --precision 256 --output tmp/gf16-advanced-r2-d095-q97.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-r3-d099-advanced-start.json --retarget-updates 4 --keep-unresolved-partition --target-bits 60 --max-cells 900 --max-depth 22 --precision 256 --output tmp/gf16-advanced-r4-d099-q49.json
```

The separate 9% independent replay and 9.25% full assembly continue
unchanged. Neither advanced dense search is complete yet.

### Recompute Simpler Witnesses After Changing the Inner

An unresolved cell in a checkpoint can be unvisited, rather than a
failed bound. A floating screen of eight unresolved three-update cells
found that all pass simpler four-update proposals. The screen used
whole saved intervals near activities 0.25--0.425, not isolated means.
Their best simple log2 proposals range from about -353.5 to -1332.6.
These scores are diagnostics, not outward certificates.

This suggests recomputing simple witnesses before adding further
regional refinements. The new `--repropose-all` option uses only the
saved partition. It retains no old acceptance labels or numerical
bounds; each newly proposed leaf must pass the usual outward check.
It requires `--resume --keep-unresolved-partition`. Ordinary resume
still rechecks the saved witnesses. All 195 tests pass, including
checks that this mode evaluates only fresh witnesses and preserves
coverage when a fresh proposal fails. The complete four-update search is:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-r3-d099-advanced-start.json --retarget-updates 4 --repropose-all --keep-unresolved-partition --target-bits 60 --max-cells 1200 --max-depth 22 --precision 256 --output tmp/gf16-reproposed-r4-d099-q49.json
```

The search starts with all 427 cells unresolved. Its first 20 intervals
pass fresh outward checks, leaving 407 pending and none rejected at
that checkpoint. The earlier witness replay remains live until the
new run provides a better checkpoint.
No partwise-output-tilt refinement has entered the certificate path:
the sampled cells did not require it.

The simple verifier also computed identical Fourier moments twice when
bounding a newly activated state by both its expansion-weight class and
its high-occupancy tail. These two envelopes now share one freshly
computed set of moments within the same outward call. No numerical
cache is retained between calls or read from a receipt. Regression tests
check the unchanged tail bound and require a new Fourier calculation
on the next invocation. All 196 tests pass.

For the same saved interval `001100010`, a diagnostic profile reports
16.005 seconds for the old outward check and 12.089 seconds after this
change. Both report log2 upper bound -2476.32691479969. The Fourier
calculation count falls from two to one. These are single-cell profiling
results under concurrent proof work, not an encoder benchmark or a
whole-search speedup measurement. The live searches were not restarted
for this change; subsequent fresh replays use the optimized calculation.

### Fresh Parallel Dense Checks

`parallel_repropose.py` checks the saved partition with up to four
processes. Each process reconstructs the actual requested inner and
proposes a fresh witness for each assigned interval. If that proposal
misses the target, the process may recheck the saved rational witness.
Saved scores and acceptance labels never establish a bound.

The parent checks each returned interval and its fresh dyadic upper
bound. Failed intervals remain unresolved; duplicate or missing results
prevent completion. The output preserves the full partition and can
seed the adaptive scalar search. A final fresh aggregate replay remains
required even if every interval passes. All 200 tests pass, including
the new process-result and fallback checks.

```text
python -B research/workstreams/permutation_locality/gf16_packets/parallel_repropose.py --resume tmp/gf16-r3-d099-advanced-start.json --retarget-updates 4 --workers 4 --precision 256 --target-bits 60 --output tmp/gf16-parallel-r4-d099-q49.json
```

The actual four-process run has completed all 427 intervals at 9.9%,
with no unresolved intervals. Every interval passes the 2^-60 target
at 256-bit precision. Its receipt is
`tmp/gf16-parallel-r4-d099-q49.json`, SHA256
`ae1cba7f4b1b89b26f3b8afccaf87a546e4298e23925683d76ae48abcba72e65`.
The redundant serial reproposal and old-witness replay were stopped
after the parallel checkpoint contained their accepted intervals.
Their inputs and saved checkpoints are preserved. The fresh 384-bit
full assembly has now replayed this partition and regenerated q=1--48,
closing the complete 9.9% certificate stated above.

Retargeting the first 280 accepted witnesses unchanged to 10% predicts
38 misses of the per-cell target. This does not require a new encoder:
a fresh regional witness for the worst predicted interval
`00001110001` passes at 10%, with margin 916.3486180712 bits at 384-bit
precision. Its local receipt is
`tmp/gf16-r4-d10-refined-frontier-p384.json`, SHA256
`d48774bfc2df03182ee592b472b1b692091710d7fe4603fea5809b8ced0631f3`.
This selected result covers one whole interval, not the full domain.
The complete 427-cell partition was then reproposed and checked at 10%:

```text
python -B research/workstreams/permutation_locality/gf16_packets/parallel_repropose.py --resume tmp/gf16-parallel-r4-d099-q49.json --retarget-distance 1/10 --workers 2 --precision 256 --target-bits 60 --output tmp/gf16-parallel-r4-d10-q49.json
```

The two-worker 10% run reached a saved checkpoint with 40 accepted
intervals and 387 pending intervals. Windows available commit memory
then fell to about 0.5 GB. That run was stopped after preserving its
checkpoint; the full 9.9% assembly continued uninterrupted. A temporary
single-process resume with native-library thread limits reduced memory
use substantially. After constructing the actual model, that process
used about 180 MB of private memory, versus about 1.45 GB per existing
worker. Four thread-limited workers continued the search, each using
about 183 MB while checking actual cells. The serial fallback was
stopped after this new run freshly rechecked all 40 retained intervals.
The claim and checkpoint partition are unchanged. All 205 tests also
pass under these thread limits; the log is
`tmp/gf16-lowthreads-full-tests.log`. The current command uses a
PowerShell process-local environment:

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python -B research/workstreams/permutation_locality/gf16_packets/parallel_repropose.py --resume tmp/gf16-parallel-r4-d10-q49.json --workers 4 --precision 256 --target-bits 60 --output tmp/gf16-lowthreads-r4-d10-q49.json
```

### Fresh 10% Assembly

The thread-limited run completed with all 427 intervals accepted and none unresolved.
Its claim is four updates, cutoff 209715, and q=49--2048. The dense
cover's SHA256 is
`aaeb719ec69e5c491a5b86418ca035836bf2ca32d8d48c5106bf6499193caf10`.
The full 384-bit assembly replayed that cover and regenerated q=1--48:

```powershell
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-lowthreads-r4-d10-q49.json --precision 384 --dense-workers 4 --dense-log tmp/gf16-r4-d10-assembly-p384-dense --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-target-bits 48 --sparse-small-through 8 --sparse-small-tilts .00016 .00024 .00028 .00032 .0004 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 --sparse-analytic-gradient --sparse-workers 1 --sparse-log tmp/gf16-r4-d10-assembly-p384-sparse.log --output tmp/gf16-r4-d10-assembly-p384-complete.json
```

Use the same process-local native-library thread limits shown above.
The earlier sparse receipts are comparison evidence, not inputs to this
assembly. Replay and exact audit have both passed, closing the complete
10% claim stated above.

### Fresh Parallel Assembly

`assemble.py --dense-workers 4` independently replays the dense witnesses
before regenerating the sparse ranges. Each worker authenticates the
outer data, constructs the actual inner maps, and checks its assigned
intervals at the requested precision. The parent reconstructs the full
partition and requires exactly one fresh positive dyadic result per leaf.
It sums results in the same order as serial replay. Missing, duplicate,
or mismatched results prevent a complete certificate.

The default remains serial; `--dense-log PREFIX` selects separate worker
logs. All 205 tests pass. A real two-process check using the actual
two-update maps reproduces the serial dyadic sum exactly, without loading
saved numerical bounds. Its log is `tmp/gf16-dense-parallel-smoke.log`.
That check validates the replay path, not a distance claim. The completed
four-update 9.9% assembly uses this path at 384-bit precision, with
four dense workers and serial sparse regeneration. Its output is
`tmp/gf16-r4-d099-assembly-p384-complete.json`.

`assemble.py --sparse-workers 4` distributes disjoint occupancy ranges
among four processes. The default remains serial. Each process builds
and checks fresh operators, covers every assigned support, and returns
an outward dyadic bound through that invocation's process channel.
The parent verifies the returned ranges and sums all contributions.
A missing, failed, or incomplete range prevents a complete certificate.
No run receipt or saved numerical cap is imported.

Exact q=1 placement remains separate. The finer low-occupancy tilt grid
also retains its existing scope. Worker logs append their occupancy
range to `--sparse-log`. Spawned processes select this workstream's
assembler explicitly, avoiding a namesake module in the legacy tools.

All 164 tests pass, including range coverage and fail-closed checks.
A real Windows two-process run freshly covers q=2 and q=3 at the 9.25%
cutoff. Its dyadic sum exactly matches a separate serial regeneration;
the log is `tmp/gf16-parallel-smoke-main-v2.log`.

The independent 9% replay raises precision to 384 bits and uses four
workers without changing the proof operators or search grid:

```text
python -B research/workstreams/permutation_locality/gf16_packets/assemble.py tmp/gf16-birth-classes-d09-tilt3over16-q97.json --precision 384 --sparse-inner gf-birth-classes --exact-feedback --single-group-exact --max-splits 150 --sparse-target-bits 48 --sparse-small-through 8 --sparse-small-tilts .00016 .00032 .00064 .001 .0016 .0024 .0032 .0064 .008 .016 --sparse-tilts .00032 .0032 .008 .016 .024 .032 .04 .048 .056 .064 .08 .096 .128 .192 --sparse-workers 4 --sparse-log tmp/gf16-r2-d09-independent-p384-sparse.log --output tmp/gf16-r2-d09-independent-p384-complete.json
```

The interrupted dense searches have also resumed from their last saved
partitions. Each continuation freshly rechecks retained witnesses:

```text
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-regional-joint-lazy-stable-r3-d099-q49.json --keep-unresolved-partition --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-joint-lazy-stable2-r3-d099-q49.json
python -B research/workstreams/permutation_locality/gf16_packets/scalar_cover.py --resume tmp/gf16-regional-joint-lazy-r2-d095-q97.json --keep-unresolved-partition --target-bits 60 --max-cells 500 --max-depth 22 --precision 256 --output tmp/gf16-regional-joint-lazy-stable2-r2-d095-q97.json
```

The earlier two-update `single_group.py` command covers only q=1 at 10%.
Its run-receipt SHA256 is
`a77ecaadaa389de657c58317e416ed14851f429f6e46709c809ad23ffa000195`.
The final assembler can select this exact-placement calculation with
`--single-group-exact`; it then regenerates q>=2 separately and checks
the complete aggregate. It does not import the standalone receipt.
Repeating the occupancy-one calculation at 384-bit precision reproduces
the 42.1225874867-bit margin. The dense `--row-parity` option records its
comparison choice in the witness and is reconstructed during replay.

Within each sparse replay, repeated CDF folds are cached by their count
family, interval, and exact probability numerator. This only avoids
recomputing identical outward values at the same precision; it does not
load bounds from disk or change their arithmetic order in the sum.

If a dense search exhausts its work budget, `--resume <checkpoint>` loads
the saved model parameters, reconstructs the complete partition, and
rechecks every accepted witness before continuing. `--max-cells` then
specifies additional proposals. No saved acceptance label or cell
coordinate is trusted. A complete result still needs `assemble.py`,
which freshly replays the dense cover and regenerates the sparse prefix.
