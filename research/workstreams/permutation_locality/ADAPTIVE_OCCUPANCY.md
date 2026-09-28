# Adaptive occupancy bounds beyond four groups

The state-memory bound in [STATE_MEMORY.md](STATE_MEMORY.md) closes three
and four active groups. A fixed support grid becomes expensive as occupancy
grows. This extension keeps the same ideal uniform one-column code and
replaces that grid with an adaptive partition. The input length is K=2^20,
the inner is IMT(128,19), and the bad output-weight threshold is 209715.

## Bound larger same-epoch collisions without enumerating every kernel

An epoch has 32 four-bit windows. Suppose j distinct windows are active,
with column weights a_1,...,a_j in {1,2,3,4}. Condition on all but one
group's windows and masks. Distinct window images intersect only at zero.
For any target feedback value, at most one remaining window and one mask
can complete the sum to that target. Choose the final group with largest
binom(4,a_i). Every atom of the feedback distribution is at most

    alpha = 1 / ((33-j) max_i binom(4,a_i)).

This also bounds its zero atom, but does not determine that atom. In the
state-memory transfer, use alpha to upper-bound a term proportional to the
zero probability, and use one to upper-bound its complement. Substituting
1-alpha for that complement would be invalid.

The existing exact censuses remain in use for j<=4. For j>=5,
`coarse_epoch` in `occupancy_memory.py` uses the atom bound above and the
same pointwise output-weight bounds as the nine-coordinate invariant.
It also uses the fresh component's density bound 1/32: convolution with
the new feedback has every atom at most min(alpha,1/32). This is valid for
every mixture admitted by the fresh-state invariant.

The resulting transfer depends on the input shape only through its total
weight W and max_i binom(4,a_i). Enumerating four weight multiplicities
constructs exactly these pairs, avoiding 4^j ordered shapes. The script
checks the compressed set against all ordered tuples through j=8. It also
checks zero-tilt mass at j=3,4,5,8,16,32 and compares the applicable entries
against the detailed three-/four-group envelope.

These local operators are available through j=32, the number of windows.
That availability is not a certificate for every global group occupancy.

## Partition the support domain with an exact budget

Fix an occupancy q. Each active group has support in {38,...,256}, giving
219^q ordered support vectors to cover. A node stores a sorted tuple of
support intervals and its exact number L of distinct label assignments.
Its volume is L times the product of its interval lengths.

To split an interval that appears r times, bisect it into lower and upper
halves. Create r+1 children, putting k copies in the lower half for each
k=0,...,r. Multiply the parent's label count by binom(r,k). The children's
volumes sum exactly to the parent's volume. Each labeled support vector
belongs to exactly one child.

For a node of volume V, assign failure budget

    2^-55 V / 219^q.

The Cauchy witness and authenticated outer-count method are unchanged from
STATE_MEMORY.md. An accepted node must have an outward upper contribution
strictly below its assigned budget. Otherwise it is split. A singleton
that fails is reported as unresolved; hitting a node limit reports the
remaining support volume. Neither outcome is described as a certificate.

The verifier checks exact volume conservation and only reports complete
coverage when every support vector is covered. The sum of the accepted
budgets is then 2^-55. Direct enumeration on small domains checks the
geometry, including the label multiplicities and disjointness of leaves.

The common 55-bit budget is a per-occupancy target, not a full-code claim.
Even summing 2048 such budgets would cost only eleven bits. Actual totals
are usually smaller, and previously closed occupancies retain their own
verified totals.

## Verified results

For exactly five active groups, the adaptive replay visits 5,880 nodes
and accepts 3,682 leaves. Every support and every group location is covered.
Both 192-bit and 384-bit replays give failure upper

    < 2.204277e-19,

or **61.9763 bits** of margin. The 384-bit total begins
2.20427636529944252044068362798862384636132149191143186291e-19.

Exactly six groups also pass a 192-bit outward replay: 29,082 visited
nodes, 18,469 accepted leaves, failure upper <5.537987e-19, and **60.6473
bits** of margin. Its upper begins
5.53798693845068065816951769985949289047919329809209462709e-19.
This six-group result has not yet been repeated at 384-bit precision.
Together with the earlier occupancies, the verified contribution through
six groups is below 3.182165e-15, about 48.1589 bits. Seven and more groups
are not included in that total.

Reproduce from the repository root:

    python -B research/workstreams/permutation_locality/test_occupancy_tail.py
    python -B research/workstreams/permutation_locality/occupancy_adaptive.py --groups 5 --precision 192
    python -B research/workstreams/permutation_locality/occupancy_adaptive.py --groups 5 --precision 384
    python -B research/workstreams/permutation_locality/occupancy_adaptive.py --groups 6 --precision 192

The earlier seven-tilt screen needed 6,205 leaves and gave 59.96 bits. The
current default also includes .0032,.004,.005,.0064,.008,.01, improving
some boxes and reducing the cover size. These tilts are analysis witnesses;
they do not change the encoder parameters or setup distribution.

Selected eight-group support points also pass the binary64 screen. Their
worst sampled contribution is about 2^-180.29 at support 80 in every group.
This is not full eight-group coverage; mixed support vectors must also be
covered. Larger occupancies and the 2x latency target remain open.

## Where the current extension stops

The volume-budget method stopped at its exploration limits for seven and
eight groups: respectively 50,000 and 20,000 visited nodes. No singleton
was unresolved, but most of the support domain was still pending. These
partial accepted sums do not certify either occupancy.

`occupancy_best_first.py` instead refines the largest estimated box term.
Its heap always partitions the full domain. A successful floating search
would replay every leaf outward and check the sum against the global
target. In the initial 12,000-split probes, the log2 bounds remained +13.90
for seven groups and +119.77 for eight. Neither probe reached outward
verification. A second seven-group probe with additional witnesses and the
pair-conditioned bound below remained at +86.84 after 4,000 splits.

The optional `--balanced-witness` policy proposes probabilities by balancing
the binomial masses at interval endpoints, clamped between two diagonal
point witnesses. It also retains the old midpoint witness. This affects
only search tightness: the replay verifies the chosen rational witness,
not a claim of numerical optimality.

The larger-occupancy screen shows a mathematical gap as well as a search
cost. At sixteen groups, all with support 176, the coarse point upper has
log2 contribution +483.76, despite minimizing over ten tilts from .0032
through .0256. More box refinement cannot repair that particular point
bound. This is not an observed low-weight codeword or a lower bound on the
failure probability.

## A tested but insufficient pair-conditioned improvement

For j>=5, fix j-2 groups instead of j-1. There are (34-j)(33-j) ordered
windows left for the remaining pair. The number of completions to any
nonzero target is at most the peak count in the full pair census; a zero
target has no completion because the pair map is injective. Let p2 be the
largest normalized nonzero atom over the exact ordered weight-pair census.
Then every atom of the j-group feedback is also bounded by

    p2 * 32*31 / ((34-j)(33-j)).

The optional `--pair-conditioned` setting takes the minimum of this bound
and the single-window conditional bound. It does not change the encoder.
The default remains the earlier bound so the five-/six-group commands
above reproduce their reported covers.

The pair-conditioned sixteen-group point at support 176 improves only to
log2 +479.45. At support 80, the corresponding improvement is +58.79 to
+57.94. Thus this refinement is not sufficient to close the larger
occupancy gap. The next proof experiment should separate loss in the
inner-state moment bound from loss in the outer support counts, rather
than merely growing the support grid.

## Isolate count loss and fresh-state moment loss

The follow-up keeps the encoder unchanged and examines sixteen groups with
equal support u. These are selected binary64 points, not a support cover.

`occupancy_count_gap.py` applies the exact joint-support OA shell caps from
`joint_oa_caps.py`. The required dual-distance-at-least-30 split certificate
was replayed successfully again. For rank h, an h-dimensional subspace has
|GL(h,2)| ordered bases. Multiply its ordered-h-tuple shell cap by the ratio
of spanning ordered four-tuples to ordered bases. Dependent tuples in the
OA population only enlarge this upper bound. Intersect with the existing
rank CDF cap, then sum ranks. No differences of upper CDFs are used.

At u=176, the per-group log2 count drops from 461.95185 to 439.50652;
across sixteen groups this saves 359.12524 bits. At u=80 and u=128, however,
the OA cap does not improve the existing counts.

`occupancy_fresh_moment.py` separately tightens the fresh component's
outgoing mass and refresh bounds. It enumerates the actual expansion
weights of each single-input feedback distribution. For one new input it
also enumerates the exact output-weight histogram for each old/new shape
pair. For more inputs it uses the resulting spectrum with a triangle bound.
Maximizing these averages over fresh shapes respects the mixture invariant.
The empty lazy branch retains its old pointwise multiplier: an average
cannot replace that multiplier while preserving the fresh density bound.
Cancellation and mature density entries remain unchanged.

The following table compares sixteen-group log2 contributions. Positive
entries are vacuous upper bounds. The first column uses the previously
tested pair-conditioned envelope; the next adds fresh averages; the last
also tightens the outer shell counts.

| Common support u | Pair-conditioned | Fresh averages | Also OA counts |
|---|---:|---:|---:|
| 80 | +57.943 | +51.816 | +51.816 |
| 128 | +323.932 | +312.260 | +312.260 |
| 176 | +479.447 | +462.571 | +103.445 |
| 216 | +19.935 | +2.577 | -290.045 |

Thus neither refinement closes the sampled problem points. A sensitivity
experiment then deletes the mature lazy-return entry C-to-Z. This is
**not a valid probability bound**. Even that counterfactual, with both
refinements, leaves log2 +118.43 at u=128 and +8.29 at u=96. Improving only
this return term cannot close the current envelope with the tested tilts.

The next useful experiment should retain total outer input weight together
with support size. The current entrywise envelope independently maximizes
column shapes, and can pair permissive inner transitions with large outer
counts. A joint support/weight calculation can test whether those two
pessimisms are compatible. This is a proposed analysis refinement, not
evidence that the candidate itself has bad distance.

Additional reproduction commands:

    python -B research/bch_spectrum_work/bch_spectrum_codex_bundle/code/certify_bch_shift_rank_split.py --verify
    python -B research/workstreams/permutation_locality/occupancy_count_gap.py
    python -B research/workstreams/permutation_locality/occupancy_fresh_moment.py
    python -B research/workstreams/permutation_locality/occupancy_fresh_moment.py --remove-mature-return

The last command is deliberately labeled counterfactual. None of these
screens changes the verified through-six-group total or the measured 1.9x
encoder speedup. Full-code proof coverage and the 2x timing target remain open.

The follow-up in [ALL_ONE_COLUMNS.md](ALL_ONE_COLUMNS.md) tested total
input weight and then isolated all-one columns. Total weight did not close
the difficult points. The all-one-column constraint instead verified one
sixteen-group support class at 53.20 bits; it is not full occupancy coverage.
