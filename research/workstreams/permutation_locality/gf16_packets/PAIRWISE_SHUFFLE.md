# Two Independently Shuffled Row Pairs

This track tries to retain some shared-route locality without requiring a
rank-four BCH support count. Its stable ID is `pairwise4-gf16-r4`.
The independently shuffled, certified four-bit route remains the reference.
Earlier proofs and the shared-route plateau are preserved in the
[proof index](../PROOF_INDEX.md).

**Status, 2026-09-29: stopped at a documented search plateau.** The full
>10% / 40-bit certificate is open. The exact support analysis, sparse
bounds, partial dense cover, and matched implementation remain available.
No production code or paper claim changes.

## Construction and Target

Use BCH[256,128] at K=2^20, N=2^21, with 2,048 groups of four rows.
Within each group, rows 0 and 1 share one uniform coordinate permutation;
rows 2 and 3 share another. The two permutations, and those of different
groups, are independent. Pack the four bits in each column into GF16.
Retain the independent uniform regional packet permutations and independent
uniform nonzero GF16 multipliers. The IMT(128,19) inner uses four updates,
zero initial state, and no terminal flush. All randomness is sampled once
at setup and retained for every message.

The target is relative distance strictly greater than 10%, except with
probability less than 2^-40 over this setup. Equivalently, bound the
probability of any nonzero message with output weight at most 209,715.
No whole-code certificate is currently claimed. Matched implementation
timings and their scope appear below.

## Exact Conditional Support Distribution

Fix an arbitrary four-row message tuple before sampling setup. Let a and b
be the union-support sizes of its two row pairs before routing. The pair
permutations produce independent uniform subsets A and B of [256], of sizes
a and b. Bits from different rows occupy different packet coordinates;
their union has no cancellation. The active-packet support is A union B.

For n=256 and J=|A intersection B|,

\[
 \Pr[J=j]=\frac{\binom a j\binom{n-a}{b-j}}{\binom n b},
 \qquad U=a+b-J.
\]

Conditioned on U=u, simultaneous coordinate-permutation invariance makes
the support uniform among the u-subsets. Condition further on the two
routed pairs. Multiplication by independent uniform elements of GF16*
makes every nonzero packet independent and uniform in GF16*. Consequently,
given u, the existing GF16 conditional inner bound applies exactly as an
upper bound. Its state-transition estimates are not asserted to be exact.

These statements hold for each fixed message. They do not assert
independence between different messages evaluated under the same setup.
The first-moment bound does not need that independence.

## From Pair Counts to an Expected Four-Row CDF

Let F(a) count ordered BCH row pairs with union support at most a, including
the zero pair. Its total is 2^256. Let C(a) be an integer upper bound on F(a),
with C(0)=1, C(256)=2^256, and C nondecreasing.

The rank-one and rank-two subspace bounds suffice to construct C. An
h-dimensional subspace has product_(j=0..h-1)(2^g-2^j) ordered spanning
g-tuples. Divide each existing four-tuple rank cap by this number for g=4,
round down to bound the integer subspace count, then multiply by the
corresponding number for g=2. Sum h=1,2 and add the zero pair.

Define comparison masses c(0)=C(0) and c(a)=C(a)-C(a-1) for a>0.
These are not upper bounds on the actual shells. Instead, after division
by 2^256, they describe a stochastically smaller support size than the
uniformly selected actual row pair. To couple the resulting subsets,
choose a uniform ordering of [256] and take prefixes of the coupled sizes.
Use independent copies for the two pairs. Both comparison subsets are
contained in their actual counterparts, so their union is also contained.

Therefore the expected count of nonzero four-row messages with support
at most u is bounded by

\[
 D(u)=\sum_{a,b=0}^{256}c(a)c(b)
       \Pr[|A_a\cup B_b|\le u]-1.
\]

The subtraction removes precisely the all-zero four-row message. In
particular, D(0)=0 and D(256)=2^512-1. The expectation here is over the two
pair permutations; D is not a spectrum bound for every realized setup.
Independence between different groups permits products of these expected
count measures in the outer first-moment sum. The existing monotone-CDF
folding step then applies. Ceilings of D provide integer upper CDFs for
the existing implementation, without altering the endpoint or zero mass.

### Exact Evaluation

The implementation uses a quadratic-time binomial transform, not a
floating overlap approximation. For comparison shell masses c(a), let

\[
 L=\operatorname{lcm}_{a=0}^n\binom n a,\qquad
 x_a=\frac{Lc(a)}{\binom n a},\qquad
 y_j=\sum_{a\le j}\binom j a x_a.
\]

The comparison mass supported in a particular j-subset is y_j/L.
Independent union squares this quantity. Binomial inversion gives its
shell mass at u:

\[
 \frac{\binom n u}{L^2}
 \sum_{j=0}^u(-1)^{u-j}\binom u j y_j^2.
\]

All operations before the final division are integers. Tests compare the
result against every subset-size pair through n=6 and exhaustive small
binary-code counts. Separate tests check support uniformity conditioned
on union size and the direction of CDF domination.

### A Separate Dense-Range Comparison

The binary BCH dual has minimum distance at least 30; the existing
production-generator check authenticates that premise. Thus a uniformly
selected ordered pair of codewords has independent uniform two-bit symbols
on any 29 coordinates. Its union weight W has the same polynomial moments
through degree 29 as Bin(256,3/4).

For the four-ary Krawtchouk polynomials K_j, their squared norms under this
binomial distribution are 3^j binom(256,j). Define the degree-14 reproducing
kernel L(w,x)=sum_(j=0..14) K_j(w)K_j(x)/(3^j binom(256,j)). Orthogonality
gives E[L(w,W)^2]=L(w,w). Its contribution from W=w alone is
Pr[W=w] L(w,w)^2, so the actual pair shell count is at most
floor(2^256/L(w,w)). Unlike CDF differences, this bounds each shell directly.

Intersect these bounds with the rank-based pair CDF caps, then use the
existing exact positive Bernoulli-majorant routine. Two independent pair
components with activity probabilities p and q produce a union component
with activity p+q-pq. Their masses multiply. Add the two cases with one
zero pair separately, excluding only the unique pair of zero pairs.
Every other comparison component retains the active-group label, including
its artificial empty-support atom. `pairwise_moments.py` tests the exact
orthogonality identity and small-code shell inequalities. This is an
alternative counting bound, not an additional construction randomizer.

### Quaternary Enumerator Constraints

An ordered binary pair (x,y) identifies x+omega*y in the GF4 scalar
extension of the BCH code. Its symbol weight is exactly the union weight.
The scalar extension has dimension 128 over GF4; its dual is the scalar
extension of the binary dual. Let A_u and B_j be their symbol-weight counts.
The quaternary MacWilliams identity imposes

\[
 4^{128} B_j=\sum_{u=0}^{256} A_u K_j^{(4)}(u).
\]

Both zero coefficients are one, and both enumerators sum to 4^128.
`pairwise_macwilliams.py` combines these equalities with the shell and CDF
bounds above, using 446 variables and 1,030 exact inequalities. Floating
linear programs only propose multipliers. Rational verification retains
nonnegative multipliers and charges every positive residual against the
variable's known upper bound. Failed solver proposals supply no improvement.
Small-code tests verify every constraint against exhaustively counted pairs.

The same module also provides a solver-free inverse-transform bound.
For a selected primal shell, replace the dual transform's coefficients by
their nonnegative suffix maxima. The resulting sequence is decreasing,
so summation by parts bounds its functional using the dual CDF. This use
of CDF differences is a monotone comparison, not a shell-count assertion.
It is tested separately from floating feasibility and exact dual replay.

Explicitly, let D_j bound the dual pair CDF, including B_0=1, and set
h_j=max(0,max_(v>=j) K_u^(4)(v)) for j>=1. Then

\[
 A_u\le 4^{-128}\left(K_u^{(4)}(0)+
       \sum_{j=1}^{256}(D_j-D_{j-1})h_j\right).
\]

Indeed, h_j is a decreasing upper bound on K_u^(4)(j). Summation by parts
expresses its functional as a nonnegative combination of dual CDF values.
Replacing those values by D_j therefore increases the result. A second
bound sums the positive transform coefficients against actual shell caps.
Take the smaller bound, round down to an integer, and intersect with the
prior shell cap. This calculation handles the high-weight endpoints where
the floating LP fails.

### Reducing Redundant Comparison Mass

A valid mixture bounds shell u by M_u=sum_j c_j Bin(256,p_j)[u].
The greedy construction can count the same shell under several components.
`positive_prune.py` proposes multipliers x_j in [0,1] to reduce this mass.
Its linear program uses floating arithmetic only to select a candidate.

Let T_u be the candidate shell mass and S_u the required shell cap. If
T_u<S_u, the original bound M_u>=S_u permits the repair fraction
(S_u-T_u)/(M_u-T_u). Take the maximum of these fractions and zero, and
round upward to a dyadic delta. The final coefficient is
c_j(delta+(1-delta)x_j). Every required shell is then covered exactly,
and no original coefficient increases. Existing empty-atom and component
mass budgets remain valid. The verifier checks all shells after repair.

A separate joint experiment pools the activity centers from several valid
majorants before pruning. It takes the maximum coefficient at each center,
after merging equal centers within each input majorant. This pool dominates
each input majorant. Its repaired result need not dominate, or be dominated
by, an earlier selected mixture. Accordingly, each pooled-mixture probe
rebuilds and checks its own comparison model.

### Choosing the Pruning Objective After Pair Union

The shell inequalities do not determine a unique comparison mixture.
Optimizing the tilted mass of one pair ignores how two pairs combine.
For pair components (c_i,p_i), define u_ij=p_i+p_j-p_i p_j. At an activity
tilt z in (0,1], the nonzero four-row comparison has tilted mass

\[
 F_z(c)=2\sum_i c_i(1-p_i+p_i z)^{256}
       +\sum_{i,j}c_i c_j(1-u_{ij}+u_{ij}z)^{256}.
\]

The first term covers one nonzero pair; the second covers two. To choose
the pruning multipliers x_i, use the gradient of F_z(c_i x_i) at x_i=1.
Its ith coordinate, omitting the common factor two, is

\[
 c_i\left[(1-p_i+p_i z)^{256}
      +\sum_j c_j(1-u_{ij}+u_{ij}z)^{256}\right].
\]

These positive costs replace the earlier one-pair objective in the
floating proposal LP. Exact repair still verifies every shell inequality.
The objective is a search surrogate, not the SPIN failure probability.
Neither the construction nor the outward inner verifier changes.
Tests compare this gradient with exact rational evaluation on small inputs.

## Progress Ledger

Goal started 2026-09-29. Preserve prior routes and stop only at a complete
independently replayed certificate with matched performance, or a documented
plateau. A plateau requires at least three substantive experiments using
distinct approaches with no meaningful improvement, after checking obvious
numerical and search artifacts. An open proof is not a completed certificate.

| Checkpoint | Exact result / measured progress | Remaining work |
|---|---|---|
| Initial support comparison | Five exhaustive/unit tests pass. At u=112, the initial shared-route log2 CDF cap is 307.3985; pairwise is 238.9586. | This compares initial bounds only, not the best refined shared bound. No distance margin follows directly. |
| Conditional inner transfer | Given support size, routing and GF16 labels match the existing conditional model. | Fresh all-support occupancy bounds and dense coverage are required. |
| Initial distance checks | At 256-bit directed precision: q=1 gives 49.5953 bits, q=23 gives 121.1298, q=32 gives 62.3760. | Only these occupancies, not the intervening range. Initial q=48 proposal fails (log2 upper +574.8355). |
| Coupled counting and finer output tilts | Fresh q=33,40,48 bounds give 306.8478, 178.8075, 50.7569 bits, respectively. | q=48 now clears the per-occupancy target. Full sparse prefix and independent replay remain. Both count refinement and the tilt grid changed, so the gain is not attributed to either alone. |
| Initial dense probes | Outward log2 bounds at tilted means .0005, .004, .032 are -866.8003, +6281.2626, +65113.8735. | Singleton diagnostics only; the two positive values are failures, not certificates or counterexamples. |
| Pairwise implementation | Full-encoder/dense-inner/adjoint checks pass; matched streaming median 7.4071 ms versus 7.4779 ms independent control. | Only 0.95% faster in this implementation. The proof improvement has not produced the desired large speedup. |
| Pair moment shell majorants | Exact two-bit moment/kernel tests pass. Outward dense log2 probes at .0005, .004, .032 are -619.4460, +9029.6681, +63810.7169. | No new coverage: some values improve, others worsen. Do not count this as a closure or a uniform improvement. |
| Independent sparse regeneration | Fresh 384-bit regeneration completed all q=1--48. Their exact dyadic sum gives a restricted margin of 48.71746838 bits. | q=49--2048 remain open; this is not a whole-code certificate. |
| Quaternary LP, central shells | Exact-verified gains of 1.41--1.70 bits at u=128,144,160,176,192; no gain at u=112. | Shell bounds only; not yet propagated into a dense cover. |
| Quaternary LP, high shells | Exact-verified gains of 7.8111 bits at u=240 and 16.7555 bits at u=248. | Proposals fail at u=252,256. These are numerical/search failures, not evidence that the initial caps are optimal. |
| Central component plus residual | A mass-2^260 pair component with activity 3/4 is followed by an exact positive residual majorant. Dense log2 probes are -866.1140, +6541.4321, +68291.4857. | No new coverage. Singleton mean values remain .0005, .004, .032. |
| Bounded component mass | Capping each residual pair component at 2^260 changes those probes to -866.8035, +6285.2397, +64916.3149. | Eliminating large masses alone does not close the dense gap. |
| Excluding nearly deterministic centers | Restricting pair centers to at most 7/8 gives log2 bounds +6239.8477 and +56227.0972 at means .004 and .032. | A tighter failed bound, not new certified coverage. The exact shell inequalities still hold. |
| Solver-free high-shell transform | Exact reductions of 20.5904, 31.2216, and 40.4497 bits at pair weights 248,252,256. | A genuine shell-count improvement, not an end-to-end margin increase. |
| Retuned dense search | With base tilt 3/16 and 16 variance bins, singleton means .004 and .032 give log2 bounds -847.3347 and +2572.6629. | Changing the base tilt changes which compositions a mean denotes; these are not the same classes as the earlier 1/32 probes. No interval cover follows. |
| Transform propagated into dense search | Fresh exact shell reconstruction and majorant validation retain -847.3347 and +2572.6629 at those two probes. | The stronger high-shell bound alone does not improve these diagnostics at the displayed precision. |
| Combined comparison | Add the per-pair component mass cap 2^260 to the transformed, endpoint-restricted, retuned model. Outward log2 probes become -868.8207 and +1291.5995. | The harder probe improves by about 1,281 bits but still fails. This is progress in a bound, not a certified whole-code margin. |
| Exact-repaired component pruning | At base tilt 3/16 and mean .032, pruning lowers the outward log2 bound to +208.2487. The repair fraction is 11/35184372088832. | All shell inequalities are checked rationally; the remaining positive bound is still a failure. |
| Finer variance partition | Increasing from 16 to 32 variance bins gives +188.0529 at that same probe. | Only a 20.20-bit improvement. This does not explain the full remaining gap. |
| Raw-mass objective | Changing the mixture cost tilt from 1/4 to 1 gives +4814.7479 at mean .032. | Worse diagnostic; keep the original objective. Component families differ, so this is not a uniform comparison of the same classes. |
| Joint component selection | Pooling five valid majorants before pruning gives +1355.0981 at mean .032. | Worse than pruning the original majorant. No coverage is imported from this family. |
| Shell caps folded back into sparse CDF | Summing genuine shell caps and intersecting with the pair CDF closes q=49 at 50.7488 bits. | q=64 fails after 96 splits, with floating log2 bound +266.0034 and no outward certificate. The sampled CDF gains are zero through support 144, so do not attribute q=49 closure to a large counting improvement. |
| Winning-variant tilt search | Inspection found that the old search retuned the uniform-feedback variant even when another regional variant won. `pairwise_search.py` also retunes the winner. Six additional tilt factors leave the best outward log2 bound at +188.0529. | The search omission is corrected, but it does not explain this failed diagnostic. The outward verifier is unchanged. |
| Central-shell batch | Exact verification accepts improved shell bounds at 69 of the 81 supports from 144 through 224. | These are local shell bounds, not a distance certificate. Twelve supports retain their previous bound after failed numerical proposals. |
| Central-shell bounds propagated | Freshly replaying all 69 shell witnesses yields the same positive mixture as before. The outward log2 bound at mean .032 remains +188.0529. | The local counting gains do not improve this dense bound. Saved numerical shell caps are ignored; only exact-checked dual multipliers are used. |
| Limiting-shell diagnosis | The repaired pair majorant meets its cap at support 256, and nearly meets it at 66,90,38. The smallest ratios at other supports exceed one. | Central-shell tightening cannot help this chosen mixture until a limiting constraint or the mixture selection changes. This identifies bottlenecks of the comparison, not low-weight codewords. |
| Targeted limiting-shell LP | A bounded search at 18 supports near those limiting constraints yields no improvement. The raw LP succeeds without improving the cap at 66 and 90; the attempts at 38 and 256 fail numerically or reach the time limit. | This does not prove optimality of the caps. All failed or unhelpful candidates retain the existing exact bounds. |
| Half-tilt one-pair objective | Replacing the mixture cost tilt 1/4 by 1/2 gives an outward log2 bound of +554.7785 at mean .032. | Worse than the retained +188.0529 diagnostic. The comparison family changes; this is not an actual-distance measurement. |
| Union-gradient objective | Keeping the original majorant and using the union gradient at z=3/16 gives an outward log2 bound of -141.50133084 at mean .032. | A newly successful singleton diagnostic. It does not cover a nonzero-width interval or the whole dense domain. A fresh full-domain search is required. |
| Independent interval check | At 384-bit directed precision, the interval with center .032 and radius 1/65536 has outward log2 bound -138.21759150. | This verifies a nonzero-width comparison cell, not the complete dense partition or a whole-code margin. |
| Fresh assembly driver | `pairwise_assemble.py` checks the ensemble and partition, rebuilds the dense model, replays every leaf, and regenerates q=1--48. It rejects the actual incomplete checkpoint before generating any certificate. | The driver is ready, but no successful full assembly has been run. The current dense partition still has holes. |
| Parallel search validation | Eight cells produced by four numerical workers replay successfully at 384-bit precision. Each worker uses the parent's identical exact comparison components. | This checks the search plumbing, not the complete domain. The replay correctly reports 297 unresolved cells. |
| Wider interior probes | At 256-bit precision, singleton means .064, .128, .192 give outward log2 bounds -321.8867, +662.9469, +603.0121. | Two interior points still fail. Subdivision alone cannot repair those particular point bounds. This is not evidence of bad codewords. |
| Wider regional tuning | At mean .128, five regional variants and sixteen additional output-tilt trials leave the outward bound at +662.94694774. | The omitted search branch does not fix this diagnostic. The original simpler bound remains best. |
| New limiting-shell audit | The union objective makes shell 121 nearly tight. Exact LP checks at supports 116--128 improve four caps, but not shell 121. | Gains at 123,126,127,128 are 0.51--1.47 bits. Propagation gives an identical rational mixture and the same failed dense bound. |
| Pooled union objective | Pool five valid majorants before applying union-gradient pruning. Outward log2 bounds become +6717.0315 at mean .128 and +1576.1569 at .032. | No new coverage. This alternative comparison fails both selected diagnostics. |
| Finer interior variance partition | Doubling from 32 to 64 bins leaves the .128 and .192 log2 bounds at +662.9469 and +603.0121. | No meaningful improvement; finer uniform variance bins do not explain these failures. |
| Fresh obstruction replay | Regenerate the model at 384-bit precision and replay the selected mean-.128 witness without searching or trusting saved bounds. The result remains +662.9469477396. | Confirms this failed upper bound numerically, not an actual low-weight codeword or impossibility theorem. |
| Direct pairwise reads | Correctness checks pass, but matched median is 8.0306 ms versus 7.4165 ms with the initial tile copy. | Keep the copy-based path as the performance baseline. |

Sources: `pairwise_support.py`, `pairwise_sparse.py`, `pairwise_dense.py`,
`pairwise_moments.py`, `pairwise_macwilliams.py`, `positive_prune.py`,
`pairwise_search.py`, `pairwise_cover.py`, `pairwise_parallel.py`, `pairwise_retune.py`,
`pairwise_assemble.py`, and their tests.
Local numerical receipts are ignored under `tmp/pairwise-goal/`; they are
not replacements for fresh verification. The full proof-tool suite passes
291 tests after the joint-mixture, winning-variant search, shell-witness
replay, limiting-shell diagnostic, union-objective, interval-scope, assembly,
parallel-search, target-aware search-stop, width-budget, regional-shortlist,
and selected-cell replay additions (`tmp/pairwise-goal-tests-final-diagnostics.log`).
Four archive-integrity tests also pass (`tmp/pairwise-goal-tests-archive.log`).
These include the endpoint-exclusion,
inverse-transform, exact-repair, and cover scope checks. A mistakenly
overlarge local return census (degree four) was stopped before generating
any bound; the distance runs use the established degree-three census and
lazy-density cutoff six. This is an execution correction, not a proof result.

Reproduce the initial exact comparison:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_support.py --output tmp/pairwise-goal/baseline-cdf.json
python -B -m unittest discover -s research/workstreams/permutation_locality/gf16_packets -p test_pairwise_support.py
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_sparse.py --occupancies 1 23 32 48 --max-splits 64 --output tmp/pairwise-goal/baseline-sparse.json
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_sparse.py --occupancies 33 40 48 --refined-counts --coupled-counts --tilts .016 .02 .024 .028 .032 .036 .04 .048 .064 --max-splits 128 --output tmp/pairwise-goal/coupled-sparse.json
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_dense.py --minimum-groups 49 --probe .0005 .004 .032 --output tmp/pairwise-goal/baseline-dense-probes.json
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_dense.py --pair-shells --refined-counts --coupled-counts --cost-tilt 1 --zero-bits 80 --minimum-groups 49 --probe .0005 .004 .032 --output tmp/pairwise-goal/moment-dense-probes.json
```

The 384-bit prefix run passes the explicit occupancy list 1 through 48,
uses `--refined-counts --coupled-counts --max-splits 128`, and output tilts
`.00016 .00024 .00028 .00032 .0004 .00064 .001 .0016 .0032 .005 .008 .012
.016 .02 .024 .028 .032 .036 .04 .048 .064 .096`. Its log is
`tmp/pairwise-goal-prefix1-48-p384.log`. Exact upper endpoints, not floating
proposals, determine each accepted occupancy.

The completed receipt contains exactly one non-null bound for every
occupancy 1--48. An exact rational scope-and-sum audit confirms that their
sum is below 2^-40. This receipt audit is separate from the fresh 384-bit
regeneration that produced the bounds. Receipt SHA256:
`b454052532f95c80ae1503a13f65f9442732ab094901f55c1261b33deb45506c`.

The later 256-bit q=49 result extends the contiguous verified prefix to 49.
An exact dyadic sum of that result and the 384-bit prefix gives a restricted
margin of 48.40174407 bits. The q=64 entry is null and is not included.
The dense model still starts at q=49; any final assembly must account for
that overlap, for example by adding only q=1--48 to a complete dense cover.

The dense diagnostics deliberately report positive log2 bounds as failures.
They do not demonstrate low-distance codewords. At the base tilt 1/32,
the floating outer dual for the .032 probe concentrates on near-certain
activity components. For the central-plus-residual majorant, its leading
component has activity about .99847 and mass about 2^549.98. This diagnosis
motivated the mass cap and endpoint exclusion; neither experiment closes
the proof. The retained independent-row certificate also uses a different
base tilt (3/16) and 16 variance bins. Matching those search settings is a
separate experiment, not an import of that certificate into this ensemble.

New receipt SHA256 values:

| Receipt under `tmp/pairwise-goal/` | SHA256 |
|---|---|
| `transform-high-shell-probe.json` | `02230d195aedc3873d7a8bf316eca62050f29747879e61d33685cb14e65400d3` |
| `retuned-dense-probes.json` | `c335f6b1357d15835d402cf0b3b51204e75714eaf95d52f5416beb4cf99bfeb5` |
| `transformed-dense-probes.json` | `a435ba9e9cf87a0d2700397b07857f8dc2759ed4dc2011accfca592947ccf7ec` |
| `pruned-dense-probes.json` | `8172d8fa8cb06e900182d72ccac6e5e97812bd9945b0331561f5523877951db3` |
| `pruned32-dense-probes.json` | `b5b8a8cee54eeec5afbcd5a7d02acdbcbd8643fdbe01cbff552e76cf071920bc` |
| `shell-cdf-sparse-probes.json` | `6ad1bd77c059478d98bac33db1d2dce60f92035515b17df9da78dd923858e8fc` |
| `retuned-winner-dense-probes.json` | `38db585a14c68a6259c9f6de834692385471565b7a7259f4fc946f065e1c4f6e` |
| `macwilliams-central-shells.json` | `e37394fc0df1b2e415ed6477fd24c3e6f64e947bbb7a0ff2765de50e5fbf33d3` |
| `lp-shell-dense-probes.json` | `ab5bed2e774440d510a0c84df29209063e84bf53e942a8323657b4d6c7ef5629` |
| `macwilliams-limiting-shells.json` | `73cbe82baea80aae86167c53ca4b50f9403fe1f5f47d6405641001c01433174b` |
| `halfcost-dense-probes.json` | `34a85b651a1b6722394eebd146ef2932b116b2db15d200c5f7c062f160d7a4fb` |
| `unioncost-dense-probes.json` | `6082ae1d214388bb33d901ff5b649a2f2d9b257246ffa86eab451d6f077f851f` |
| `unioncost-interval-p384.json` | `edf9cf19412300d545e1ea2625cdc7332a2403a60cba015c27913f3aff099756` |
| `dense-cover-pruned32.json` | `0b9dd2fdda4dd5c60ba2020cb5422390119d2999345a81a2f1b17ca6143c237e` |

`pairwise_cover.py` searches the entire dense comparison domain, rather
than sampled means. Each checkpoint records accepted and unresolved
intervals. Resume regenerates the exact model and recomputes retained
bounds. Fresh replay checks partition coverage and sums outward bounds;
it cannot report a complete dense result while any interval is unresolved.
Even a complete dense result covers only q=49--2048. It must be combined
with the sparse prefix before claiming a whole-code certificate.

An initial bounded search requests 80 cells, depth at most 22, and 64 bits
per accepted cell. Its checkpoint is
`tmp/pairwise-goal/dense-cover-initial.json`; its log has the same stem in
`tmp/pairwise-goal-dense-cover-initial.log`. It uses the uncapped comparison
components. It finished at 80 visited cells with two outward-verified
intervals and 77 unresolved intervals. Unresolved intervals are pending
subdivisions, not counterexamples.

The next cover uses `--mass-bits 260 --prune-mixture --variance-bins 32`.
It imports only the previous partition topology through `--partition-from`:
every cell is reproposed and checked against the new model. It finished at
cumulative visit 240 with 47 accepted intervals and 145 unresolved intervals. The checkpoint
is `tmp/pairwise-goal/dense-cover-pruned32.json`. This is still incomplete.
Normal resume/replay never changes a saved model. A changed comparison
requires a separate output and complete rechecking.

The next same-model run rechecks and retains all 47 accepted intervals
before searching the unresolved partition. Its bounded work request is
160 additional cells, recorded separately as
`tmp/pairwise-goal/dense-cover-pruned32-next.json`. A running checkpoint
does not imply completion. At cumulative visit 320, all 47 intervals
remain accepted and 225 intervals are unresolved. The extra visits only
subdivided unresolved intervals at this checkpoint. After the successful
union-objective diagnostic, this old-comparison run was intentionally
stopped, not mistaken for a completed search. Its last saved checkpoint is
preserved with SHA256
`7b114d4b82341c5347b29ed4f700c5f17575040bc45c4820fec103631f2da0af`.

The union-objective cover imports only these subdivision paths and
rechecks every interval with the new comparison. It uses
`--partition-from tmp/pairwise-goal/dense-cover-pruned32-next.json
--mass-bits 260 --prune-mixture --union-cost-tilt 3/16 --variance-bins 32`,
with a bounded request of 240 additional cells, depth 24, and 64 bits per
accepted interval. Its output is `tmp/pairwise-goal/dense-cover-unioncost.json`.
No previously accepted interval is imported as already proved under the
new comparison. A separate 384-bit probe checks the interval with center
.032 and radius 1/65536. That run completed with outward log2 upper bound
-138.2175914965, after independently regenerating the shell majorant and
inner operators. Its exact interval is
[262019/8192000, 262269/8192000]. This is a separate cell check, not an
accepted leaf automatically inserted into the full cover.

At cumulative visit 440, the union-objective cover has 87 accepted leaves
and 218 unresolved leaves. All accepted bounds were freshly computed for
the new comparison. Its visit count includes the inherited partition's
320 visits; the first 120 new visits rechecked or subdivided imported cells.
This serial run was deliberately stopped at its saved checkpoint to use
four numerical workers. It did not finish the requested 240-cell budget.
The checkpoint remains unchanged with SHA256
`16d3d7bfe252c4518808cc5955bb0b22d4b01b0831ad2d397bf5cf78399f5879`.

### Parallel Continuation

`pairwise_parallel.py` regenerates the exact shell comparison in the parent.
Each worker receives the same rational components and constructs fresh
actual inner maps with its own Arb context. Component hashes, cell paths,
precision, and positive outward bounds are checked before accepting results.
The parent saves a complete partition atomically after every batch.
Previously accepted leaves are re-evaluated, not imported as proved.

A four-worker, eight-cell smoke test completed. Its receipt is
`tmp/pairwise-goal/parallel-smoke.json`, SHA256
`c61c904b309c1e14b2ec23c39f7efaef5b6948519becde6124b029ccd71cf008`.
An independent serial replay regenerated the model at 384-bit precision
and successfully checked all eight retained cells. That replay is
`tmp/pairwise-goal/parallel-smoke-replay-p384.json`, SHA256
`2fd15330f941329f288ed3daeb625a135ab4a2038fc961b57c7a1af3361e790e`.
It correctly reports `dense_complete=false` because 297 cells remain unresolved.

The real continuation starts from the preserved serial checkpoint, not
the smoke-test subset. It requests at most 700 evaluations, including
rechecks, with depth 24 and 64 bits per accepted cell:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_parallel.py tmp/pairwise-goal/dense-cover-unioncost.json --workers 4 --max-cells 700 --output tmp/pairwise-goal/dense-cover-unioncost-parallel.json
```

The worker stops proposal refinement once a floating candidate beats
66 bits; it then requires an outward bound below 2^-64. The previous
90-bit search cutoff remains the default for other callers. This change
avoids unnecessary search and does not weaken the acceptance test.
Fresh final replay does not use either proposal cutoff.

For progress tracking, the serial checkpoint covers 78.7109375% of the
comparison-parameter interval by length, versus 73.046875% before the
union objective. These percentages are not fractions of codewords or
bad-setup probabilities. Uncovered intervals remain full proof obligations.
The worker log and every checkpoint retain the distinction between a
floating proposal, an outward-verified cell, and complete domain coverage.

The parallel run rechecked all 87 inherited leaves and verified eight
more. At cumulative visit 540 it has 95 accepted leaves and 215 unresolved
leaves. Accepted cells cover 79.4921875% of the comparison interval by
length. Some unsuccessful cells were split, so the unresolved count is
not a monotone measure of progress. This checkpoint is preserved with SHA256
`dac3505e010348c3bd4456d9d58b6f0e9c4550c4181e9b69cf947b4ee098b9f2`.
The search was deliberately stopped after the interior probes below found
failed singleton bounds. Continuing uniform subdivision would not resolve
that obstruction. No running process or completed dense certificate is
implied by the checkpoint.

The optional `--aggregate-budget-bits b` allocates a single dense budget
across the binary partition. A leaf at depth d receives 2^-(b+d).
For any complete binary partition, sum_leaf 2^-d=1. Thus strict outward
acceptance of every leaf implies a dense aggregate below 2^-b, regardless
of the number of leaves. Choosing b=44 would leave room for the verified
sparse prefix and the final 40-bit target. This changes search effort,
not the construction or outward bound; the final assembler still sums
the freshly recomputed endpoints. The stopped run used the older fixed
64-bit per-cell target, not this optional allocation.

### Interior Obstruction and Regional Search

`tmp/pairwise-goal/unioncost-interior-probes.json` records the three fresh
singleton probes, with SHA256
`899c286c42c5fe9833d8f5ee7d95525eb297da1d9517b7c5328e8340e53c66ca`.
The mean-.064 witness uses the uniform-feedback regional refinement.
The failed .128 and .192 witnesses use the simpler variance-partition
bound without regional count terms. Therefore the existing winner-only
regional retuner does not explore their initially worse regional candidates.

`pairwise_retune.py` checks this search limitation. It reconstructs the
same comparison, proposes five regional variants, and retunes the two
best regional variants at eight output-tilt factors from 1/4 through 2.
The simpler baseline remains a candidate throughout. It saves exploratory
scores separately, then checks the final winner with outward arithmetic.
Its receipt covers only the stated cell, not the full dense domain.

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_retune.py tmp/pairwise-goal/dense-cover-unioncost-parallel.json --mean .128 --output tmp/pairwise-goal/regional-retune-0128.json
```

That run completed without improving the baseline. The receipt SHA256 is
`99f6a3fe547dc99b8e1ed8b64b5a92243abda1ad51920c2c3b66bfb97738ec69`.
All five initial regional bounds and all sixteen retuned candidates are
recorded; the final outward bound is +662.9469477396 in log2 units.

The union objective changes which shell inequalities constrain pruning.
Its limiting-shell audit includes support 121, which was not covered by
the earlier central-shell batch at 144--224. A new bounded LP search
checks supports 116--128. It improves shells 123,126,127,128 by
0.5109664, 1.4717944, 1.4313756, and 1.4656041 bits, respectively.
No improvement is obtained at 121; two other raw attempts hit their
time limits, so no shell optimality claim follows.
The exact-witness receipt is
`tmp/pairwise-goal/macwilliams-union-limiting-shells.json`, SHA256
`6ed1207639aa45db53d8f6e5c5a086c870ace6dec5ebd48499ee7d477e248cda`.

Three further diagnostics completed, with separate outputs:

- `unioncost-shell-refined-probe.json`: replay the new exact shell witnesses
  and propagate them into the same mean-.128 comparison. The resulting
  rational component list is identical, and the bound remains +662.9469.
- `unioncost-joint-probes.json`: apply union-gradient pruning to the pooled
  majorant family, checking means .128 and .032. This comparison changes;
  it cannot inherit the accepted interval cover. Both bounds fail, at
  +6717.0315 and +1576.1569, respectively.
- `unioncost-interior64-probes.json`: double the variance partition to
  64 bins at means .128 and .192. The bounds remain +662.9469 and +603.0121;
  the displayed differences from 32 bins are below 2e-12 bits.

All outputs are under `tmp/pairwise-goal/`. The corresponding logs use
the same stems with prefix `tmp/pairwise-goal-`. These are proof
computations, not encoder performance measurements. No benchmarks ran
concurrently. Their sessions completed successfully; these are failed
proof bounds, not interrupted computations.

| Receipt | SHA256 |
|---|---|
| `unioncost-shell-refined-probe.json` | `920235cf1d2b2e6e77227b77a8a43a0f2f9429fb7a2499fc5b5613ec2551eaec` |
| `unioncost-joint-probes.json` | `365247c84b92354f20fe67a23303889735ff80e7013b90e92fb0a446e66392ce` |
| `unioncost-interior64-probes.json` | `7d66838aef3380331672e969cf814b2aa1bc4f9812b670d0ecca8f63203732f4` |
| `regional-retune-0128-replay-p384.json` | `5c7ca9539ff9e1d2b65d252251c84d08267836470b50692fee96d20d4475d0d1` |

The last receipt rebuilds the original comparison at 384-bit precision
and re-evaluates the saved witness. Saved numerical scores and bounds
are not inputs. The outward log2 endpoint remains
+662.9469477396126725874054564383. Reproduce it with:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_retune.py tmp/pairwise-goal/regional-retune-0128.json --replay --precision 384 --output tmp/pairwise-goal/regional-retune-0128-recheck.json
```

### Plateau Decision

The earlier cover extension was progress: it verified eight additional
cells. The subsequent bounded experiments establish the requested stopping
criterion for this approach. They tested three distinct directions:
regional/output-tilt search, outer comparison selection, and exact BCH
shell constraints. None improved the failed interior bound or extended
the cover. Doubling variance resolution also gave no improvement.
Fresh high-precision evaluation, exact shell-witness replay, identical
component checks, and the search tests rule out the numerical and search
artifacts examined here. They do not prove optimality over all possible
bounds or all possible code analyses.

The retained evidence is:

- exact conditional support and expected-CDF arguments for this ensemble;
- freshly verified q=1--48 with restricted margin 48.71746838 bits;
- an additional q=49 bound, giving q=1--49 restricted margin 48.40174407 bits;
- 95 accepted dense comparison cells and 215 unresolved cells;
- a correct matched implementation at about 7.42 ms, only 1--2% faster
  than the independently shuffled reference in the recorded runs.

These do not form a complete distance certificate. The final assembly
gate still rejects the incomplete cover. The fastest shared four-row
route also remains uncertified; its timing must not be paired with the
independent-row certificate.

For a deployable proved route, retain the independently shuffled GF16
four-update construction (>10%, >48.65 bits). For another bounded
experiment, first test a fifth IMT update against the failed interior
points, then measure its cost if it helps. That would be a new ensemble
and a new certificate obligation, not a reinterpretation of this goal.
Alternatively, keep four updates and pursue a stronger joint outer/inner
counting bound; further uniform subdivision alone is not sufficient.

Reproduce the independent interval check with:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_dense.py --pair-shells --transform-shells --prune-mixture --retune-regional --refined-counts --coupled-counts --central-bits 260 --mass-bits 260 --cost-tilt 1/4 --union-cost-tilt 3/16 --zero-bits 80 --maximum-activity 7/8 --base-tilt 3/16 --variance-bins 32 --minimum-groups 49 --probe .032 --probe-radius 1/65536 --precision 384 --output tmp/pairwise-goal/unioncost-interval-recheck.json
```

### Final Verification Gate

After the dense cover has no unresolved intervals, run:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pairwise_assemble.py tmp/pairwise-goal/dense-cover-unioncost.json --precision 384 --bits 40 --output tmp/pairwise-goal/full-replay-p384.json
```

This command is not presently a completed certificate. It rejects an
incomplete or overlapping partition, the wrong ensemble, and changed
construction parameters before expensive work. It then rebuilds the
positive comparison and actual inner maps, checks the actual fixed
geometry and state size, and recomputes every dense
bound. It regenerates q=1--48 with the pairwise expected CDF, instead of
importing earlier sparse numerical endpoints. Thus q=49 belongs only to
the dense sum; the extra standalone q=49 result is not double-counted.

The final comparison with 2^-40 uses exact dyadic arithmetic. Equality is
not accepted. Missing, duplicated, null, or zero sparse bounds and changes
to the requested arithmetic precision are rejected. Both output paths
must be new, preserving existing receipts. The main receipt records the
dense input hash, fresh sparse receipt hash, scope, and exact aggregate.
Tests exercise these rejection paths, and an invocation against the real
incomplete cover produced neither a main receipt nor a sparse receipt.

The limiting-shell audit uses exact rational majorant/cap ratios to order
the constraints, then prints approximate ratios for diagnosis. Its output
is not an additional proof premise. In particular, no parity restriction
was added to the pair spectrum: production-generator row weights include
both 0 and 2 modulo four, and 4,121 generator pairs are nonorthogonal.
The doubly-even/self-orthogonal shortcut therefore does not apply here.

## Initial Matched Performance

Peach Ryzen 7950X, GCC 15.2, core 15, K=2^20, 128-bit elements. Each
sample is a median of 101 in-place precomputed transposed encodings.
Setup, allocations, and reference checks are excluded. Two seeds use
balanced serial order 7,12,13,9,9,13,12,7; all variants use four updates
and GF16 randomization. Existing benchmark locks prevent concurrent runs.

| Route | Median of four samples | Sample range |
|---|---:|---:|
| Independent rows, streaming stores | 7.4779 ms | 7.4735--7.4852 ms |
| Pairwise rows, streaming stores | 7.4071 ms | 7.3679--7.4330 ms |
| Pairwise rows, cached stores | 13.1790 ms | 13.1467--13.3083 ms |
| Shared four rows, streaming stores | 6.5159 ms | 6.4804--6.5205 ms |

The pairwise kernel keeps full 64-byte packet stores and replaces four
lane loads during BCH preparation by two aligned 32-byte loads. The
initial implementation still copies and rearranges a four-row tile.
Consequently the shared route's main memory-layout advantage is not yet
recovered. No new production default follows from this measurement.

Normal correctness checks cover streaming/cached pairwise modes, updates
2--4, two seeds, and three input patterns. They also include a K=2^20
pairwise case and independent/shared/two-bit controls. Dense-reference
agreement, inner adjoint, routing bijection, and in-place suffix preservation
pass. Address/undefined-behavior sanitizer checks also pass for this initial
implementation, including the unchanged controls. The research encoder is
instrumented; the prebuilt BCH objects and dependency library are not.
Sanitizer executable SHA256:
`07def50272ccba14d019ae06a3c59d9ee4fa706b9b0b8fbad4250611a4e7f1df`.

Research executable SHA256:
`d7d45e0cb2fca061f26657d62338e1fa12d9fe77626edc252776c18de941a396`.
Source `joint.cpp` SHA256:
`5b929be24781f0c04bd34ccf3ada657bdf09862f1661328d1b8ab0f8b79001f1`.
Remote run directory: `/tmp/spin-pairwise-gf16-asxQRV`.
Matched logs: `measurements/pairwise-gf16-confirm-UazWus`, copied locally
under `tmp/pairwise-goal/`. The runner is `../run_pairwise_gf16.sh`.

### Direct-Read Ablation

A second implementation omits the sequential 16 KiB tile copy and reads
the paired lanes directly from routed storage. Modes 14/15 select this
research-only streaming/cached ablation; modes 12/13 retain the tile copy.
Both variants preserve the code distribution and GF16 operations.
Normal checks cover both variants, updates 2--4, two seeds, three patterns,
full-size streaming cases, and independent/shared/two-bit controls.
The direct-read variant is not claimed sanitizer-validated.

The matched serial order was 7,12,14,9,9,14,12,7 for each of two seeds,
again using 101 calls per sample on core 15:

| Route | Median of four samples | Sample range |
|---|---:|---:|
| Independent rows | 7.5336 ms | 7.5132--7.5455 ms |
| Pairwise, sequential tile copy | 7.4165 ms | 7.3939--7.4443 ms |
| Pairwise, direct reads | 8.0306 ms | 8.0197--8.0419 ms |
| Shared four rows | 6.5323 ms | 6.5246--6.5413 ms |

Removing the copy makes this kernel about 8.3% slower. The copy-based
pairwise route remains within roughly 1--2% of the independent-row
reference across the two experiments; a substantial speedup is not shown.

Executable SHA256:
`2ea9cdd49c2ff9400e815bd0367a826a47081d85f9a02ef5e04370b582d6e140`.
Benchmarked source SHA256:
`8c887f04a67a63031d166fc60d01689dde3c4b10870fc0f3074ed23e44cd3b67`.
The subsequent source-only change adds a compile-time guard forbidding
direct reads without pairwise routing; it does not change these instantiations.
Remote directory: `/tmp/spin-pairwise-direct-0Gr8FD`.
Local logs: `tmp/pairwise-goal/measurements/pairwise-gf16-direct-confirm-U1GM8g`.

The bounded investigation ends at the plateau above, not at a proved
whole-code claim. Preserve the sparse prefix, the incomplete dense
partition, and all earlier proof tracks. Any resumed attempt must retain
the full target until every dense interval is covered, freshly replayed,
and combined with the sparse prefix.
