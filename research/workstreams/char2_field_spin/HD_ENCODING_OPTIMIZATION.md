# Optimizing the certified high-distance encoder

The target is precomputed forward encoding at K=2^20, rate 1/2, with
40% GF65536-symbol distance and at least 40 bits of setup-failure margin.
Each record carries eight parallel field lanes under one shared setup.
The aspiration is below 5 ms on the established Peach benchmark host.
Setup is excluded; every encoding pass, route, and output store is included.

The [final checkpoint](HD_OPTIMIZED_CHECKPOINT.md) gives the current result,
validation, and reproduction commands. The chronological entries below
retain earlier failures and intermediate next steps; those are not a list
of work still running.

The reference is the [sixteen-state certificate](HD_SUMRANK_FEEDBACK_CERTIFICATE.md).
It proves more than 59.87 bits, but has no native timing. The older 7.02-ms
four-state encoder is a performance control, not a certified alternative.
Neither reference is to be overwritten during this work.

## Candidate ledger

| Idea | Proof status | Implementation status | Decision |
|---|---|---|---|
| s16, exact transitive feedback via uniform GL16 | Original complete certificate transfers unchanged | 13.249 ms initial full baseline | Preserved control |
| Diagonal-Cauchy feedback | Complete certificate transfers with 11.542-bit loss; >48.32 bits remain | 11.156 ms first implementation | Retained |
| Compile-time FFT constants and expanded fixed coefficients | Same maps and certificate | 10.125 ms full | Retained |
| Graded local basis and quadratic-tower Moore butterflies | Same local ranks/moments; >48.32 bits with diagonal-Cauchy feedback | 9.56 ms repeated full timing | Retained |
| Direct routed scaling, without staging | Same code and certificate | 12.48 ms | Regression; preserve diagnostic only |
| Direct routed scaling with next-region streaming prefetch | Same code and certificate | 9.06--9.12 ms repeated full timing | Preserved control |
| Outer scaling fused into FFT leaf stores | Same code and certificate | 8.394 ms against 8.954/8.980-ms controls | Retained |
| Additive-coset LRS representatives and FFT state maps | Invertible row change; identical moments and >48.32 bits | 7.927 ms against 8.666/8.637-ms fused controls | Retained |
| Sequential L2 prefetch only | Same code and certificate | 7.34--7.50 ms across repeated runs | Current candidate |
| Paired 128-bit load tree | Same code and certificate | 7.361 ms against 7.352-ms nearby control | No clear win |
| Hardware gather | Same code and certificate | 9.114 ms | Regression |
| Streaming output stores | Same code and certificate | 7.541/7.522 ms; cached controls 7.340/7.483 ms | Not selected |
| Tower-field feedback multiplication, precomputed scalar | Same construction after choosing an F-linear field basis | Not implemented | Compare if feedback remains important |
| Fixed-map factorization using the period-four Moore rows | Exact algebraic factorization | Implemented with four packed state groups | Retained |
| Smaller state with tighter chronological composition | s12 diagnostic is close but incomplete | Not yet timed | Retain as proof route |
| Low-rank random feedback | Needs a new operator for possible annihilation/residual mass | Not implemented | Do not substitute an ideal refresh |
| Arbitrary binary transitive feedback | Same one-message transition; loses F-linearity and its projective counting gain | Not implemented | Not the main path |

GF65536-linearity remains a constraint of the main path. An arbitrary
binary basis change does not preserve that interface automatically.
All reported timings must identify the implemented construction and must
follow correctness checks. Benchmarks run sequentially under the host lock.

Completion requires a validated fully certified winner and a record of
the material proof and implementation alternatives tried. A small local
improvement alone does not establish that this design space is exhausted.

## First native round

The original feedback proof only uses that a fixed nonzero feedback vector
has a uniform nonzero image, independently of the other packet coins.
A uniform member of GL(16,F) has this property. Sampling an iid uniform
matrix and rejecting singular matrices gives exactly this distribution.
Thus this native baseline retains the original certificate; it is not a
claim that GL16 and field multipliers have the same joint action on messages.
The timing harness uses reproducible MT19937 streams, not an implementation
of independently uniform cryptographic setup material.

The [diagonal-Cauchy proof](HD_DIAGONAL_CAUCHY_REFRESH.md) removes dense
feedback matrix multiplication at a rigorously bounded cost in margin.
The [tower-basis proof](HD_TOWER_BASIS_INNER.md) then changes the local
polynomial basis and the LRS basis to admit smaller arithmetic circuits.
It preserves GF65536-linearity and the certificate's exact local moments.

Measurements use Peach CPU2 at the existing 4.5-GHz userspace setting,
GCC 15.2, AVX512/GFNI, three warmups and 51 timed calls (31 for the first
dense and identity controls). Generation, setup, and checksums are outside
timing. The full call includes outer encoding, staging/scaling, route reads,
the entire recurrence, and output stores. Runs acquire both benchmark locks.

The 9.711-ms run has p10/p90 9.691/9.732 ms. Its checksum is
3549692054676245512. Separate phase measurements give outer 2.305 ms and
route/inner 7.358 ms; they need not add exactly to the full median.
No production code has been promoted. The original 7.02-ms four-state
prototype remains faster but still lacks a complete distance certificate.

After removing unused setup tables and extending the validation, a
101-call A/B/A sequence returned 9.561820 ms (tower), 10.022039 ms
(old local basis), and 9.560367 ms (tower again). The two tower runs have
p10/p90 9.529/9.677 and 9.524/9.615 ms. Their checksums match the earlier
tower run. The measured executable SHA256 is
`bedf3bce8c789dabd7d6818483339b900ec41f8e40f7faac60f1df098602574a`.
This is about 28% lower latency than the first certified native baseline,
not a speedup over the older, incompletely certified four-state prototype.

Validation passes in Release and ASAN/UBSAN: every packet input/state basis
direction in both extension components, random full-lane values, nonzero
incoming state, exact aliasing, routed offsets, uint32 routes, and guards.
Full encodings at K8192 and non-power-of-two K24576 match the independent
dense inner oracle. The unchanged outer has its previous independent
validation. The exact Python and C++ F,A,C matrix digest agrees.
All 381 high-distance Python regression tests pass. The new directed
transfer helper returns a logarithmic failure upper bound below -48.3283.

The new research target is `spin_char2_hd_sumrank_forward`. Run `--check`
before timing. Its benchmark interface is `--bench 51 full 20 tower`;
`cauchy` retains the old local basis, `refresh` uses the uniform-GL reference,
and `identity` is an explicitly uncertified control. `inner` and `outer`
replace `full` for separate component timings.

## Routing and a first state-update screen

The pre-routed arithmetic control takes 5.042 ms; the standalone
staging/scaling/routing control takes 3.398 ms. These are diagnostics,
not phases whose sum is an end-to-end timing. The arithmetic control has
the same checksum as the full encoder.

Setup can conjugate each outer scalar through the route, putting its
prepared coefficient next to its routed input index. Applying that scalar
directly after each random record load removes the staging buffer pass.
The first implementation regressed to 12.483 ms: avoiding the sequential
pass made the routed reads cold. Streaming the next region into L2 while
processing the current region resolves most of this loss. Both regions
are only 128 KiB at K20 and fit in the host's L2.

A 101-call stream/staged/stream sequence gives 9.117530, 9.645145, and
9.063349 ms. All have checksum 3549692054676245512. The two stream runs
have p10/p90 9.078/9.177 and 9.026/9.108 ms. Release and ASAN/UBSAN checks
pass for every mode, including direct scaling and the streaming variant.
The measured executable SHA256 is
`7580b45ff879b9440e1374823e9ce548973941f4138fa5a96a204b5ba31f0933`.
Use `--bench 101 full 20 stream` to reproduce the current candidate;
`direct` reproduces the cold-read regression. The research plan retains
both coefficient orderings for these A/B controls; only one is read by
each timed kernel.

An initial FLOAT screen removes the separate rank-one state mixer but
keeps ideal transitive feedback. Its point-density operator replaces
`V+Gamma/q` by `V+Gamma`, and uses `Gamma` on an empty packet. The optimized
packetwise logarithmic bounds at active rows 16,128,512,1024 are
approximately +1669,+3300,+1810,-1244. This does not close the construction;
paired-profile regional extraction has not yet been tried. The default
proved operator and all pinned witnesses remain unchanged.

## Outer fusion and FFT state maps

The scaled outer preserves the frozen outer's explicit sixteen-register
FFT leaves. It applies each GF65536 scalar just before storing its output.
The inner therefore reads already-scaled values and performs no field
product after a routed load. An A/B/A run gives 8.954395 ms for `stream`,
8.394070 ms for `fused`, and 8.979683 ms for `stream`. All outputs agree.
Aligned streaming stores, unaligned cached fallback, plane padding, and
input/output offsets are checked against the old outer plus scalar scaling.

The [FFT state-map refinement](HD_FFT_STATE_MAPS.md) chooses the eight LRS
representatives on a nonzero additive coset. An invertible change of rows
then replaces the repeated weighted sums with two four-point FFTs. The
certificate's local ranks and moments do not change. The native dense maps
agree with the independent exact Python digest 17847517350586406002.
An A/B/A comparison gives 8.666217 ms (`fused`), 7.926817 ms (`fft`), and
8.637303 ms (`fused`). The new ensemble's checksum is 10512060216375813096;
it is deliberately different from the previous ensemble's output.

The remaining random-access prefetch loop costs measurable time. With
sequential next-region L2 prefetch retained, unrolling the L1 prefetches
gives 7.578 ms; shorter lookahead gives 7.589 ms. Removing the L1 prefetches
gives 7.500 ms, against 7.904/7.886-ms `fft` controls in the same sequence.
All three variants produce the same output. Paired scalar loads show no
clear benefit, whereas hardware gathers regress to 9.114 ms.

Two 301-call cached-store controls give 7.340112 and 7.483400 ms, with
p10/p90 7.315/7.377 and 7.439/7.957 ms. Intervening non-temporal-output
runs give 7.540796 and 7.521841 ms. The latter do not improve the current
candidate. Some repeats show host/timing variation despite CPU2 remaining
at the configured 4.5 GHz; report approximately 7.4 ms, not the lowest
individual sample. All benchmarks run sequentially under both host locks.

Reproduce the current candidate with `--bench 301 full 20 l2only`.
`fft`, `prefetch`, `near`, `paired-load`, `gather`, and `nt` retain the
respective controls. These are isolated research modes, not new library
APIs. The NT path checks alignment once and falls back to the cached path;
tests cover aligned and unaligned outputs, guards, and full exact aliasing.
The proof suite now contains 383 passing tests.

Final Release and ASAN/UBSAN checks pass for all thirteen research modes.
The directed transfer helper again gives -48.3283516168753 bits as the
logarithmic upper bound. This reuses the previously replayed complete
certificate through the proved map transformations; the full original
witness search/replay was not repeated in this implementation round.
The final measured executable SHA256 is
`56ba0575b0ed537b989430c9ef08e92fdc9ec9c0e603758fd110ce191d181614`.
A final 301-call comparison gives 8.960456 ms for the previous `stream`
control and 7.469022 ms for `l2only`. Preserve the wider repeated-run range
above when describing performance rather than claiming a precise 7.34-ms
constant. No production library code or public default has changed.

At this checkpoint, next proof steps were tighter composition for smaller state or less frequent
state mixing, and narrower local tiles over a larger subfield. Next native
steps are factoring the local feedthrough and eliminating layout work that
remains between these maps. None of these construction-level possibilities
has been exhausted, and the sub-5-ms goal remains active.

## Systematic feedthrough and a global RS emission map

Replacing the tile feedthrough by identity removes an input transform but
changes the local moments. The layered systematic candidate measures
7.184 ms; fixed GFNI-affine products reduce that to about 7.1 ms. Its
adaptive certificate remains incomplete at row count 2290. It is not
promoted to the certified anchor.

The identity feedthrough also makes a larger simplification possible:
only coordinate restrictions of the emission matrix need known ranks.
Replace the layered LRS maps by A=RS[64,16] evaluation in the normalized
additive-polynomial basis and C=A^T. The packet is Y=X+A Q,
Q_next=M Q+T C X. The outer, route, mixer, and diagonal-Cauchy feedback
distributions remain unchanged. All map ranks and packet moments now
follow from the ordinary RS MDS property.

The [global-MDS proof and native checkpoint](HD_SYSTEMATIC_MDS.md) closes
40% symbol distance with more than 50.06 bits of whole-code margin.
The independent 384-bit pinned replay covers every active-row count.
Native Release checks pass for all sixteen modes. The new `mds` mode
measures 6.188 ms (101 calls), then 6.234 ms (301 calls), with checksum
976174626021906176. This is now the fastest fully certified design.
All timings exclude setup and run sequentially under both host locks.

The certified construction and all previous anchors remain separate.
No production library code is changed. The active optimization goal
continues: fixed-map circuits, feedback randomization, and state/packet
parameters still offer plausible opportunities.

The follow-up `mds-affine` mode expresses the fixed byte-field products
as GFNI affine maps, including the fixed Cauchy refresh. It has identical
outputs. A 301-call A/B/A sequence gives 6.597/6.323/6.474 ms for
ordinary/affine/ordinary. A reversed sequence gives
6.261/6.343/6.297 ms for affine/ordinary/affine. The CPU remains under
the userspace governor with a 4.5-GHz ceiling, but the observed variation
precludes reporting every small difference as a robust improvement.
Use approximately 6.3 ms for the current certified implementation.

FLOAT proof screens reject several immediate work reductions under the
current operator; these are failed bounds, not counterexamples:

| Candidate | Rows 128 | Rows 256 | Rows 512 |
| --- | ---: | ---: | ---: |
| Twelve state coordinates | +509 | +964 | +1475 |
| Fourteen state coordinates | +81 | +316 | +284 |
| Fifteen state coordinates | -96 | +72 | -142 |
| No state mixer, sixteen coordinates | +632 | +281 | -414 |
| State mixer every second packet | +176 | +71 | -452 |

Entries are optimized paired-regional log2 upper screens at 40% distance,
not achieved margins. All retain ideal transitive feedback for screening.
The alternating case multiplies an unmixed packet operator by a mixed
packet operator before the fractional power. No reduced-mixing native
candidate or certificate has been claimed. Smaller state and deleting
mixers do not presently supply a straightforward speedup at the target.
The next proof opportunity is a tighter state decomposition that reduces
the separate feedback-randomizer cost, or reallocating redundancy between
the outer and inner; either would need a new full certificate.

Further FLOAT diagnostics retain a zero/flat/residual state component.
Deleting feedback refresh gives +3282/+9455/+18566 at rows 128/256/512,
even after retaining the MDS random-span bound j/64 for j<16. Refreshing
every second packet gives +585/+1245/+2515 before that span refinement.
An independent rank-one feedback randomizer gives +439/+1688/+3798 with
the current three-component envelope. Collapsing it to an almost-flat
two-component cone makes the bound worse by reintroducing residual mass
at each boundary. A dense rank-one-correction map between random diagonals
also does not close with the good-part/residual envelope. These are
proof screens only; no such kernel replaces the certified distribution.

Packet-width screens now use the appropriate regional coefficient degree
and number of packets. At width32, eight or twelve state coordinates fail
the selected middle classes; even sixteen coordinates give about +6 at
row256. Width128/state16 gives -331/-126/+132/+4012 at rows
128/256/512/1024. Thus merely changing the batch size does not yet justify
a faster construction. Larger-width and tighter composition ideas remain
separate from the complete width64/state16 certificate.

All seventeen native variants pass Release and ASAN/UBSAN checks. The
validated executable SHA256 is
`d85d0616409bd6802baf552dce0bbcb8df496c3473cb33367ea608ca0b8d4553`.
The 391-test high-distance regression suite passes; additional generalized
screen tests are recorded alongside the exact systematic-MDS map tests.
No raw timing/proof output has been added to version control, and no
production defaults have changed.

The width256/state16 screen also fails: approximately -281/+308/+4072/
+11870 at rows 128/256/512/1024. Ten focused exact-map and generalized
coefficient tests pass after adding width and three-component screening.

At this checkpoint, the next proposed experiments were sorted source indices
within each packet and fixed GFNI-affine outer circuits against the anchor.
Sorting within a packet leaves its input set unchanged; the systematic
MDS moments depend only on its potential count, so the certificate need
not require a uniformly random order inside that set. This is a proposed
proof-preserving routing optimization, not a measured win. Sorting across
packet boundaries does not have that justification and is not proposed
as an automatic substitution. The goal remains active; these native
ideas and further circuit fusion have not been exhausted.

## Routing, outer stores, and scheduling controls

Inspection confirms that the scaled outer already uses fixed GFNI-affine
products; replacing ordinary products there was not a new optimization.
Packet-local sorting is valid by the position-independent MDS argument in
HD_SYSTEMATIC_MDS.md, but does not yield a clear timing improvement.
Sequential 201-call medians are:

| Mode | Full encoding (ms) |
| --- | ---: |
| `mds-affine`, first control | 6.184 |
| `mds-sorted` | 6.191 |
| `outer32`, versus default tile64 | 6.180 |
| `outer128` | 6.203 |
| `outer-cached`, versus streaming stores | 7.566 |
| `mds-affine`, final control | 6.318 |

Sorting and tile-size differences are within the observed run variation.
Cached outer stores are materially slower. The sorted code has checksum
13129517991152346212; the other modes retain 976174626021906176.
Release checks pass for all 21 modes, including dense packet and complete
encoder comparisons, nonzero states, unaligned buffers, and exact aliasing.
An exhaustive small permutation test checks the packet-count distribution
after sorting. No candidate here displaces the approximately 6.3-ms anchor.

Assembly inspection shows register spills in the monolithic inner. Mode22,
`bounded`, puts each sixteen-output tile in a non-inlined helper, keeping
exactly the same arithmetic map. All 22 Release variants pass. A sequential
301-call A/B/A comparison gives 6.127/6.259/6.131 ms for
`mds-affine`/`bounded`/`mds-affine`. Bounding lifetimes this way does not help.
The helpers' call and state-transfer overhead outweigh any spill reduction.

The next implementation test evaluates the four RS cosets in SIMD lanes
together. This trades internal FFT shuffles for packet-boundary record
transposes. It preserves the matrices, so it needs a correctness and timing
check rather than a new distance certificate. No speedup is claimed yet.

The resulting mode23, `striped`, passes all native dense-oracle checks.
Its 301-call median is 6.187 ms between 6.182 and 6.286-ms controls.
This does not by itself establish an improvement. The implementation
evaluates four cosets together, then transposes four 128-bit records at
the packet interfaces. It does not change A or C.

## A joint refresh saves one vector of field products

[Sharing the final diagonal](HD_SHARED_DIAGONAL.md) between the rank-one
state mixer and feedback yields a new proved joint ensemble. It removes
sixteen GF65536 multiplications per 64-symbol packet. The ideal operator
comparison costs only q/(q-1) beyond the existing diagonal-Cauchy factor.
The independent 384-bit replay checks all 78 witnesses and all 8192 row
counts, then gives a full bound below 2^-49.33892764811194 at 40% distance.

The native `shared` and `shared-striped` modes explicitly set u=D2*1
and sample v conditioned on v^T u != 1. Their dense scalar oracle uses
these same joint maps, not the old independent-map construction. Both
give checksum 4678095192827670896. All 25 Release variants pass.
Exact GF4 enumeration checks every incoming state and syndrome, both
conditioned and unconditioned, against the joint density bounds. The
actual GF65536 interpolation matrix also matches the common-factor
Cauchy formula used in the proof.

A sequential 301-call comparison gives 6.211/6.232/6.086/6.167 ms for
`mds-affine`/`shared`/`shared-striped`/`mds-affine`. A separate 501-call
reversed comparison gives 6.094/6.168/6.042 ms for
`shared-striped`/`mds-affine`/`shared-striped`. The shared distribution
helps when combined with the striped circuit; the unstriped version does
not improve the control. Use approximately 6.1 ms, not the lowest sample.

Wider-packet FLOAT screens also identify a viable next construction
direction. For width128/state24, rows128/256/512/1024 give approximately
-1462/-1256/-1769/-3536. With state32 they give -2690/-2203/-2669/-5122.
These are selected row-class upper-bound screens with ideal feedback,
not full certificates. They suggest amortizing randomization over larger
packets, but do not yet establish a native speedup or complete coverage.

Validation is complete for this checkpoint: all 25 native modes pass
Release and ASAN/UBSAN checks. The high-distance regression suite passes
397 tests, followed by all four focused joint-refresh tests after adding
the full small-field weighted-operator enumeration. The final executable
SHA256 is `46cf36539f86e8d004fe88df78a998088d6429f4599c707b341e9f0bf067a288`.
The previous independent-map binary remains frozen separately.
No raw data is committed and no production-library default is changed.
The optimization goal remains active; wider-packet co-design and further
fixed-map circuit work have not been exhausted.

A narrower follow-up at width128 gives -621/-476/-600/+734 for state18
and -899/-763/-1104/-1350 for state20, at the same four row counts.
Thus twenty state coordinates are a more promising amortization target
than eighteen under this screen. A native circuit and complete directed
coverage are still needed; non-power-of-two state does not automatically
give a cheaper FFT implementation. The validated 25-mode executable is
also frozen as `shared-anchor-46cf3653` in the isolated Peach workspace.

## Wider-packet native screening

An isolated width128 driver now tests states20,24,32. Emission is the
normalized RS[128,s] map; feedback is its transpose. For non-power-of-two
state, the refresh uses the leading s-by-s submatrix of the 32-point
Cauchy interpolation map. Input padding and output truncation implement
that matrix with the existing FFT circuit. Its Cauchy minors and augmented
all-one column have the same rank property needed by the shared-diagonal
argument. Dense packet and full-encoder oracles pass for all three states.

Sequential 301-call medians are 6.089 ms for the frozen width64/state16
control, then 6.177/6.329/6.870 ms for width128/states20/24/32, then
5.876 ms for the frozen control. Host/run variation is visible, but none
of the wider implementations wins. These candidates are not promoted.
The state20 adaptive directed search also remains incomplete: it covers
rows1..1287 plus the existing dense tail, then its FLOAT proposal lacks
margin at row1288. The generalized width support preserves the old
certificate defaults; eleven focused regression tests pass.

The next structural test replaces the global MDS feedback with cheaper
block-diagonal RS feedback, retaining the global MDS emission. Plain
bucket sums fail the coarse bound. Retaining collision profiles improves
that bound but does not close it. With one/two/four local checks per
bucket, selected paired-profile screens at rows128/256/512/1024 are:

| Local checks | Local input width | Selected log2 upper screens |
| --- | --- | --- |
| 1 | 4 | +83 / +1006 / +3318 / +9787 |
| 2 | 8 | -241 / -16 / +291 / +2442 |
| 4 | 16 | -258 / -134 / -440 / -1171 |

All use total state16 and width64. These are ideal-feedback FLOAT screens,
not complete certificates or counterexamples. The four-check candidate
merits a native test; no feedback simplification has yet been promoted.

## Block-feedback native candidates

[Block-diagonal RS feedback](HD_BLOCK_FEEDBACK.md) now has two validated
native implementations. Mode26 (`bucket8`) uses two RS[32,8] transpose
maps; mode27 (`bucket4`) uses four RS[16,4] transpose maps. Both retain
the global RS[64,16] emission and shared-diagonal refresh. The striped
implementation prunes unused feedback FFT levels; no extra routing pass
or hot allocation is introduced.

All 27 native variants pass Release and ASAN/UBSAN checks. Eighteen
focused Python tests pass, including the new block-feedback small-field
operator checks and generalized width128/256 moment checks. A reversed
501-call comparison gives 5.629/6.020/6.005/5.597 ms for four-check/
reference/eight-check/four-check. The four-check candidate is consistently
about 5.6 ms, but neither block-feedback variant has a complete certificate.

The earlier dense-class obstruction was partly an aggregation artifact.
After retaining all profiles, the eight-check global-marker screens give
-2780/-8190/-32435 at rows2048/3072/4096. Four-check paired regional
screens give about -124 and -769 at rows2048 and3072; directed checks
and the intervening row classes are still required. The previous
positive screens were bounds that failed, not low-distance witnesses.

The alternating four/eight-check candidate removes the sampled middle-class
obstruction. A common-circuit implementation measures 5.745/5.752 ms around
a 5.953-ms reference. Its full directed coverage search is running; it is
not yet a certified replacement. See the block-feedback note for the
slower dispatch controls and the basis change used to share the circuit.

A smaller outer is not a free win under the current proof. Keeping the
global-MDS sixteen-state inner, but changing the outer to RS[128,64]
at the same K20 and rate1/2, gives paired-regional FLOAT upper screens
of +540/+1124/+1827/+2545 at active rows128/256/512/1024. The generalized
screen adjusts the number of rows, regions, outer measure factor, and
regional polynomial exponent. No smaller-outer native implementation or
certificate is claimed. The factor counting active-row subsets is larger
for a given active fraction; lowering outer arithmetic also changes this
part of the proof budget.

Two further same-code controls fail to improve the alternating candidate:
the natural-layout affine circuit takes 5.95--5.97 ms, and cache-line-padded
outer coefficient planes take 5.81--5.83 ms. Nearby common-circuit controls
take 5.69--5.73 ms. Their outputs agree; neither replaces the frozen binary.
Its phase controls give 2.456 ms outer, 3.174 ms routed inner, and 2.921 ms
pre-routed inner against 5.721 ms full. The next structural screen therefore
also considers cheaper emission maps, not merely routing or feedback.

The alternating four/eight-check proof search has now finished. Its 104
pinned rational witnesses cover all 8192 row counts at 256-bit directed
precision. After the full shared-refresh density loss, the setup-failure
bound is below 2^-51.04 at 40% distance. A separate 384-bit replay passes
all 8192 row counts and gives 51.04261653036409 bits. This makes the
approximately 5.7--5.8-ms common-circuit version
the latest certified research checkpoint; the 6.1-ms reference is preserved.

[Block emission](HD_BLOCK_EMISSION.md) supplies the next cost candidate.
Two RS[32,8] emission/feedback blocks do not beat the certified version.
Four RS[16,4] blocks take 5.484/5.506 ms around a 5.719-ms control, but the
current rank-one-refresh proof fails some middle classes. All 33 native
variants pass Release and ASAN/UBSAN checks. The new controls are frozen
as `block-emission-control-5a476509`, SHA256
`5a476509bdc1a36698baade8952ceea68f6a643333abdf9357d3e881b1b3d316`.

Tracking active state-block counts has not yet repaired the cheaper
emission proof. A proposed rank-two shared refresh instead reduces its
identity branch from q^-1 to q^-2. Initial ideal-operator screens at rows
256/512 give -133/-261, compared with failed bounds for the rank-one
version. The next decision requires local proof validation and a native
cost test: an extra dot product might erase the arithmetic saving. No
complete rank-two certificate or timing is claimed. The optimization goal
remains active; this tradeoff and longer state-path composition are not
exhausted, and the sub-5-ms target remains unmet.

## Carrying state through the feedback circuit

Mode34, `bucket4-carry`, changes the recurrence to
Q_next=D2 H(Q+D1 CX), retaining block emission and feedback. It removes
the random dot product and passes all 34 Release and ASAN/UBSAN variants.
Sequential 501-call medians are 5.367/5.373 ms, versus 5.696/5.718-ms
frozen certified controls and 5.517 ms for the rank-one block-emission
control. Its checksum is 2743298128861683080. The validated executable is
frozen as `carried-state-control-mode34`, SHA256
`daf6cd932491a17893b2e5591be52025662614631d6c3b1d3635bb9414b6d24b`.
The existing rank-one certificate does not transfer.

Mode35, `bucket4-common`, instead applies D2 H D1 to Q+CX. The
[common-refresh local proof](HD_COMMON_REFRESH.md) accounts for cancellation
explicitly and contains no lazy identity branch. All 35 native modes pass
Release and ASAN/UBSAN checks. Two 501-call medians are 5.318/5.328 ms,
with 5.777 ms for the frozen certified control and 5.338 ms for mode34.
Checksum: 15132924349672518052. The validated executable is frozen as
`common-refresh-control-mode35`, SHA256
`8fe2508c7b3e4d85cebe0ca288713c3613a00861636903d6d3246c5023efc65b`.
This is a cost candidate, not a certified replacement.

The first common-refresh bound fails at rows256/512/2048. A positive
zero-constraint expansion, using an exhaustive principal-nullity scan,
improves those screens to approximately +130/+298/+5185. A proposed fixed
MDS mixer between feedback blocks improves the sparse samples to about
-35.6/-34.6, but still fails at row2048. Neither closes the target. That
additional mixer has not been implemented or timed.

The active-state-count experiment also finishes without a certificate:
the paired-average screen gives about +44983 at row2048. Rank-two shared
refresh gives about +515 there. These failed bounds are not counterexamples.

The next proof step is a sharper enumerator for the common-refresh return
event. An exact zero-set-lattice implementation counts local outputs over
all potential supports of size at most four; independent small-field
enumeration tests pass. Sparse-support averages and full coverage remain
separate obligations. A feedback-coordinate transpose or small fixed block
mixer could improve return weights, but neither is a validated winner.
The fastest fully certified implementation remains approximately 5.7--5.8 ms.

The exact local-return enumeration finishes for all potential supports of
sizes 1,...,4. Its integer coefficients are pinned in the enumerator source;
replay takes about 210 seconds for size 4 on the local proof host. Using
those averages gives FLOAT screens +67.7/+4.0/+4498 at rows256/512/2048.
Thus tightening the unmodified common-refresh return is helpful but still
does not certify the candidate.

Mode36, `bucket4-scatter`, shifts each local feedback evaluation to
nodes 16,...,31 and transposes its sixteen check coordinates. The latter
removes a native record transpose. Its two 501-call medians are
5.324/5.349 ms versus 5.369 ms for mode35 and 5.797 ms for the frozen
certified control. The checksum is 13462303532037942408. All 36 Release
and ASAN/UBSAN checks pass; the dense-map digest is 1225673016892679168.
The measured executable is frozen as `scattered-feedback-control-mode36`,
SHA256 `5227dd54c14d58b6c74233f160bcc7c64e682c6f62030c10d7ab57f89ed1fd5b`.

The scattered map materially improves the one-input return: its weight is
56--63 instead of 15--16. Exact one/two-input projective enumerators and
finite check-rank tables supplement its proof screen. At rows256/512/2048,
the current FLOAT bounds are -50.1/-39.5/+23267. This is progress in the
sparse classes, not a complete certificate. The medium-density return
moment is the next proof target; the minimum-weight cap loses too much.
No extra fixed block mixer is included in this native timing.

The common-refresh note also records a deterministic rank-budget bound:
for identity feedthrough and uniform packet feedback rank f, d_min/N is
at most f/32+1/8192 at K20. Thus f<=12 cannot reach 40% distance. A cheap
four/eight-rank global feedback map is ruled out before implementation.
The active goal still has live proof/design directions and is not complete.

The focused regression run passes 21 tests across common refresh, block
emission, shared diagonals, and rank-two refresh. No benchmark or proof
process remains running at this checkpoint. The next useful step is a
joint output-zero/return bound for the scattered map at medium densities;
its sparse column improvements alone cannot establish the target.

The next return bound retains individual inactive-output zero constraints.
It uses the natural RS degree order to obtain a block-triangular rank
bound, then conditions on the syndrome to mark the remaining input
degrees of freedom. The corrected worst-labelling FLOAT screen at row512
is -133.99, improving the previous -39.5. An averaged-labelled refinement
passes five local regression tests; medium-class screens are in progress.
An earlier best-assignment prototype assumed arbitrary RS minors were
MDS. A concrete rank-zero counterexample invalidates its apparent -310.4
result. That result is discarded and the counterexample is a regression
test. No native timing or certificate claim changes.

The corrected worst-labelled row2048 screen finishes at +11182.47; the
labelled-average refinement gives -136.62/+11042.00 at rows512/2048.
The [return-bound investigation](HD_RETURN_DESIGN_SCREEN.md) records
additional local proofs and tests for a Lagrange basis, fresh scalar,
state-only refresh, and global emission. Their current bounds do not
close the target, and none has been assigned a native timing. An invalid
zero-return diagnostic isolates the return term as a useful proof target;
it is explicitly not a certificate or probability bound.

The rank-two/larger-feedback route is more promising. Alternating
four/eight-check feedback closes the sampled row2048 screen but fails
row3072. Using two eight-check feedback blocks throughout gives negative
screens at rows512/2048/3072/4096, and the row3072 witness passes the
256-bit directed verifier at an ideal log2 bound below -1247.08. Exact
rational packet generation and a complete-coverage driver are now present;
the full search is in progress, not yet a whole-code certificate.

Native modes37/38 pass all 38 Release and ASAN/UBSAN variants. Mode38
(`bucket8-rank2`) takes 5.583/5.579 ms around a 5.742-ms frozen certified
control in the stable reversed sequence. This is about a 3% candidate
improvement. Mode38 has checksum7812899659201065080 and fixed-map digest
6742237189988151056. The frozen executable is
`rank-two-eight-mode38-25334b78`, SHA256
`25334b787171b442a4500989045dca68bcd435c0cb97809ad38ff880662e555d`.
No certified checkpoint is replaced. Next, finish directed coverage and
independent replay before treating this timing as a certified improvement.

## Complementary-transition audit and arithmetic controls

The rank-two/larger-feedback grouping had a proof defect: max(K) is an
upper bound for the quiet transition, but T-max(K) is not an upper bound
for the complementary transition across different block labellings. The
corrected provider uses max(K) and T-min(K) in the respective entries.
The profile (0,0,8,8) is an exact regression counterexample to the old
choice. Entrywise checks cover all labellings of that profile and one
hundred further profiles. All earlier rank-two partial receipts are
withdrawn; the corrected full search is running. This does not affect
the fully certified alternating four/eight-check checkpoint, which uses
its own exact occupancy operator. The 39-test focused proof suite passes.

The frozen mode38 phase medians are 2.434 ms outer, 3.097 ms routed inner,
and 2.845 ms pre-routed inner, from sequential 301-call runs. A delayed
byte-swap dot product preserves the map and passes all 38 Release variants,
but loses in sequential full 501-call comparisons: 5.818/5.851 ms versus
5.760 ms for the rank-two control, and 5.751 ms versus 5.686 ms for the
certified control. It remains disabled behind HD_SUMRANK_DELAYED_DOT.
Its frozen executable is delayed-dot-control-20261007, SHA256
b5a042740dec0cd07595a1b9badf4a0c9c8818f8c0ddb59dd10029b1114c6bd6.

Cycle sampling identifies the first outer coset's source loads as a useful
diagnostic target. Two explicit prefetch controls keep the exact map:
mode39 prefetches one four-row group ahead into L1, and mode40 prefetches
sixteen rows ahead into L2. All 40 Release variants pass. Outer medians
are 2.496/2.586/2.398/2.428 ms for control/L1/L2/control. Full medians
are 5.684/5.806/6.125/5.617 ms in the same order. Neither control is
selected: an isolated outer improvement does not establish an end-to-end
improvement. An aligned streaming-output inner is the next exact-map test.

## Streaming output: new certified implementation checkpoint

Mode41 (`bucket8-rank2-nt`) preserves the rank-two construction but writes
the sequential output through aligned non-temporal stores. Mode42
(`bucket48-fused-nt`) applies the same change to the fully certified
alternating four/eight-check construction. Each full call checks output
alignment once and falls back to the corresponding cached implementation
when needed. A single fence completes the output stores before return.
There is no new allocation, no change to setup, and no additional permutation.
These timings concern encoding alone, not an immediately following consumer.

All 42 variants pass Release and ASAN/UBSAN validation. Tests cover both
scalar basis components, arbitrary incoming states, routed offsets, exact
aliasing, and full encodes at K=8192,16384,24576,32768 for the new paths.
The streaming full-encoder tests additionally exercise all four record
alignments, output guards, and aligned input/output aliasing. A setup-only
constexpr mapping ties implementation controls to their original ensembles.

Mode41 takes 5.419/5.437 ms versus cached rank-two controls of 5.721/5.688 ms.
The following sequential 501-call comparison is:

| Implementation | Full median |
| --- | ---: |
| Frozen certified mode29 | 5.675 ms |
| Certified streaming mode42 | 5.464 ms |
| Rank-two streaming mode41 | 5.463 ms |
| Certified streaming mode42 | 5.547 ms |
| Frozen certified mode29 | 5.746 ms |

Report about **5.5 ms** for the new fully certified implementation. Its
checksum remains 16994121132285933612, identical to the frozen mode29
control. It computes exactly the same linear map and retains the complete
51.04261653-bit certificate at 40% distance. The rank-two streaming checksum
is 7812899659201065080. At this checkpoint it still needed its corrected
certificate; the completed search and replay are recorded below.

The validated executable is frozen as
`streaming-output-mode41-42-fa594acc`, SHA256
`fa594accd9e34809ca8dc7fb86aedaa47a88a0e73f519a37b29bdb97234c90df`.
The old 5.7--5.8-ms executable and its witnesses remain unchanged. The
sub-5-ms target is still unmet; next finish corrected rank-two coverage
and independent replay before deciding whether its smaller emission is
a worthwhile certified alternative.

## Route-layout and outer-load controls

Modes43/44 prepare route indices in the inner FFT's striped input order.
The gathers then avoid 32 record-shuffle instructions per packet. The
original route is retained for the scalar oracle and fallback. Setup
allocates the extra table; encoding does not allocate. Output still uses
the natural order, and both checksums match the corresponding old modes.
All 44 variants pass Release and ASAN/UBSAN checks.

Two reversed 501-call sequences show no consistent full-encoding win:
the certified striped-route mode takes 5.455/5.470 ms against controls
5.420/5.512 ms. The rank-two mode takes 5.421/5.381 ms against controls
5.405/5.363 ms. The executable is frozen as
`striped-route-mode43-44-319d763b`, SHA256
`319d763b1885969712d5cbd5648bc415e7652851ed454fa20ab13f815b968209`.
It is not selected over the existing route.

`HD_SUMRANK_SHARED_OUTER_LOADS` reads each source vector once for both
outer FFT cosets, retaining the established leaf circuits. It passes all
44 Release variants, but full encoding takes 5.782 ms between controls
5.438/5.445 ms. The outer alone takes 2.763 ms versus 2.487 ms. Sharing
loads also changes the scratch and store schedule; it does not reduce
field arithmetic. The disabled control is frozen as
`shared-outer-loads-a37806e4`, SHA256
`a37806e47bef1c814cac67288a70ce97496ccb1274622114d469f70478c0b21f`.

`HD_SUMRANK_ROUTE_BYTES` precomputes 32-bit byte offsets, removing an
address-scale shift at each load but doubling the route table to 8 MiB.
All 44 Release variants pass. Certified/rank-two medians are 5.774/5.621 ms,
against 5.475/5.389-ms controls. The disabled executable is frozen as
`byte-route-3c14c78e`, SHA256
`3c14c78e4728aaafd090c94ebf96065cc8faf33d7b025d136a1d041e0a5feef4`.

`HD_SUMRANK_ROUTE_UNITS8` instead stores 16-bit offsets in eight-byte
units. This retains the 4-MiB table and permits a scaled addressing mode.
Regions above 32768 records fall back to the ordinary route. All 44
Release variants pass. Results vary across reversed runs: certified
5.755/5.861 ms versus controls 5.704/5.573 ms; rank-two 5.808/5.516 ms
versus controls 5.650/5.758 ms. These runs do not establish a repeatable
improvement. The disabled control is frozen as `word-route-412d6044`,
SHA256 `412d6044ab4200c65190a661a20973cc2ce87f47215d4ff620cb676a710cd616`.
Boundary-load regressions additionally cover the largest representable
offsets. All 44 variants also pass ASAN/UBSAN with this address control.

## Corrected rank-two closure and a lighter feedback schedule

The corrected rank-two search has completed at 256-bit precision. Its
133 pinned rational witnesses cover every active-row count. The ideal
union logarithm is below -63.25965982455041; the 12.26301241395966-bit
density loss leaves more than 50.99664741059075 bits at 40% distance.
The [separate 384-bit replay](theory/hd_rank_two_block_witnesses.py) confirms
all 8192 classes and a final logarithm below -50.99664741059075.
This result applies to modes38/41/43, not the unclosed common-map
refresh candidates. The existing 51.04-bit certificate is unchanged.

A [lighter feedback schedule](HD_RETURN_DESIGN_SCREEN.md) repeats
(4,4,4,8), retaining global emission and rank-one refresh. Full-profile
FLOAT screens are negative at the sampled middle counts, and a directed
receipt certifies counts 2560 through 2683. The subsequent 256-bit search
covers all 8192 counts using 119 pinned rational witnesses. Its ideal
union logarithm is below -68.92219588903863; the unchanged density loss
leaves more than 56.65919448196863 bits. A separate 384-bit replay confirms
all 119 witnesses and all 8192 classes. Native mode45 implements it, as
described below; sub-5-ms encoding is not claimed.

## Period-four implementation and register lifetime

Mode45 (`bucket4448-fused-nt`) repeats four/four/four/eight-check feedback.
It retains the global emission, shared rank-one refresh, and existing route.
The phase uses the global packet index, including when regions contain
an odd number of packets. The scalar oracle uses the same explicit period.
Aligned output uses streaming stores; unaligned output retains this period
in its cached fallback. All 45 variants pass Release and ASAN/UBSAN checks.
Both the default build and the compact setup control described below pass
the sanitizer suite.
The K20 checksum is 7582624837270167152.

The first comparison exposed a scheduling regression in the newly rebuilt
control: both new modes took about 5.7--5.8 ms, while the frozen control
remained near 5.5 ms. The striped-route experiment had moved output stores
past the input transpose, extending output-register lifetimes. Restoring
the original store order recovered the control. A sequential 501-call
comparison gives frozen/current mode42 at 5.484/5.486 ms, mode45 at
5.448/5.432 ms, and current/frozen mode42 at 5.492/5.503 ms. The observed
schedule gain is about 1%, not a new sub-5-ms result.

The restored executable is `restored-order-mode45`, SHA256
`f2c4f666658d329b9b3692bec0653c2a421594025cc0ae7bb1623e1a8adc9af4`.
It preserves both certified constructions and the separately identified
period-four candidate. Its first, regressed build remains available as
`lighter-feedback-first-control`, SHA256
`a9966f773323350f6a7b03358ee0bbac3cb4794f0b262598a2579cdfa0df3430`.

`HD_SUMRANK_COMPACT_STEPS` places the shared rank-one dot coefficients
and two refresh diagonals in one 192-byte, cache-line-aligned step record.
It removes the unused u coefficients from the hot setup stream, without
changing the sampled code. Full checksums match. Reversed 501-call runs
do not establish a repeatable gain: period-four candidate 5.509/5.525 ms
lies between controls 5.447/5.530 ms. The option remains disabled.
Its all-45-variant Release and ASAN/UBSAN checks pass. The executable is
`compact-steps-mode45`, SHA256
`5d0883a5bef0b6c0bfcd8c5bbb1073bdabd26f8262e829944c2c576f317f05b3`.

`HD_SUMRANK_RESTRICT_PAYLOAD` supplies a non-aliasing contract only at
the full inner's internal source/output boundary. Full encoding first
consumes the input into separate scratch, so input/output aliasing at the
encoder interface remains valid. Packet helpers retain their alias support.
All 45 Release variants pass. Reversed period-four medians are
5.454/5.455/5.433/5.437 ms for restricted/control/control/restricted.
The hint gives no repeatable gain and remains disabled. The executable is
`restrict-payload-mode45`, SHA256
`74bc8f17e7284bcd93168f38a1e10e8b769ee47ea16bd440ff9f62d9022765ac`.

An [exact circuit screen](theory/hd_cauchy_circuit_screen.py) factors H
as commuting dyadic operators. In the group algebra, write epsilon_i=I+P_i,
where P_i XOR-shifts the coordinate index by 2^i. Each epsilon_i squares
to zero. The screen factors H into fifteen terms I+c_S product_(i in S)
epsilon_i and checks every basis column against independent Lagrange
interpolation. It also checks H^2=I. The straightforward packed factor
schedule needs 32 fixed GFNI products, against 20 in the current inverse/
forward FFT pair. This does not establish a lower bound for all circuits;
it rejects this particular factor schedule as an arithmetic reduction.
No native timing is claimed for it. The focused proof regression suite
still passes all 36 tests.

A final same-executable, reversed 501-call comparison gives:

| Construction | First / reversed full median |
| --- | ---: |
| Certified global emission, alternating feedback, mode42 | 5.467 / 5.487 ms |
| Certified block emission, rank-two refresh, mode41 | 5.462 / 5.420 ms |
| Global emission, period-four feedback candidate, mode45 | 5.478 / 5.444 ms |

The rank-two implementation is the fastest independently replayed alternative
in this sequence. Differences between it and the period-four candidate are
small; no new sub-5-ms result is established. The latter's full 256-bit
certificate has now closed, and its independent 384-bit replay passes.
Its larger certified margin reflects different witnesses and bounds, not
a demonstrated monotonic benefit from reducing feedback work.

## Final target-size audit

The 384-bit period-four replay passes every witness and covers all 8192
classes, with final log2 bound below -56.65919448196863. Default Release
and ASAN/UBSAN validation pass all 45 variants. The new `--large-check`
compares every K20 output record for modes41,42,45 with the dense scalar
inner, including aligned/unaligned output, guards, and full aliasing.
All three pass and match the frozen checksums.

The final help-text rebuild passes the complete Release suite. A final
501-call reversed sequence gives mode41 5.457/5.460 ms, mode45
5.451/5.440 ms, and mode42 5.551/5.482 ms. The new validation does not
establish a timing regression. The executable is preserved as
`optimized-checkpoint-69329fef`, SHA256
`69329fefb8f3e4020488411adfcaa8c0aa694eaf43fdc8b29aaa9e88c57e46c3`.
The [consolidated checkpoint](HD_OPTIMIZED_CHECKPOINT.md) selects mode45
on comparable speed and stronger margin, preserving mode41 as its control.
The search ends around 5.5 ms, not below 5 ms or with a claim of optimality.
The final complete high-distance Python suite passes all 446 tests.
The checkpoint links a hash-verified local source archive; no raw numerical
dumps, builds, or timing logs are added to Git.
