# Sparse-input bounds and the missing output coordinates

The 3.4-ms goal remains open. This pass separates three questions: which
constraints the bound forgot, which actual output patterns it overcounts,
and which construction changes could remove those difficulties cheaply.
No native encoder or paper changed. No benchmark ran.

## What the complete sparse test establishes

The candidate uses the fixed inverse input schedule described in
[the column-schedule note](HD_BIT_SYNDROME_COLUMN_PERMUTATIONS.md).
Each active RS[16,8] row has at least nine nonzero columns. Therefore each
active group has at least nine active columns, regardless of its arity.

The new comparison retains that condition before taking the group power.
It splits each column into zero input, nonzero quiet input, and nonquiet
input. Quiet means that the H* correction vanishes; it does not mean that
the output vanishes. The quiet polynomial's exact constant is removed
before powering. Entrywise Minkowski combines the two nonzero branches
under the same positive path measure. A first-nine-active expansion then
costs choose(16,9)^(r-1). Empty columns use that measure's S(0)^16, not a
smaller exact empty-region operator.

For r=18, the source class F<=1, and output mark 9927/10000, the homogeneous
controls are as follows. Each includes the same source, RS, and projective
normalizations and all labelled row sets of its arity.

| Active RS rows per group | log2 first-moment comparison |
|---|---:|
| 2 | -851.232 |
| 3 | -1232.521 |
| 4 | -615.894 |
| 5 | +423.627 |
| 6 | +105.531 |
| 7 | -281.520 |

These are FLOAT diagnostics, not certificates. The common atom proposal
still needs its recorded tiny dyadic deficit repaired before directed
replay. They do not cover other source weights or matching classes.
The first-nine controls are unconditional and need no bad-setup charge.

The implementation and receipts are
`theory/hd_bit_syndrome_sparse_gate.py`,
`theory/hd_sparse_activity_control.py`, and the corresponding
`hd_sparse_activity_a{2,3,4,5,6,7}_9927_float.json` files.
They reuse the authenticated common regional receipt without another
provider run. Mixed arities must be combined through their r-th roots,
not by summing the homogeneous comparisons in the table.

## A bounded constraint test that did not close the gap

The exact setup exclusion has also strengthened. Across all8192 groups,
except with probability at most2^-61.665825339754, groups of arities
five, six, seven, and eight have at least five, four, three, and two
nonquiet columns respectively. The event is counted exactly using the
rank profile of the scheduled quiet kernels. The setup probability is
added once globally; it does not permit conditioning and then reusing
independence. See `hd_bit_syndrome_nonquiet_expurgation_exact.json`.

The next test combined this condition with an input-weight marker. Quiet
input weight equals quiet output weight, so its marker is exact there.
For a nonquiet column on a allowed rows, X^a bounds its input marker.
The word receives X^(-9a) per group and the nonquiet marker receives its
separate required floor. No second first-hit combinatorial loss is added.

The prescribed fifteen marker pairs improve arity five only to +420.445
and arity six only to +105.177. Nonquiet markers worsen the comparisons.
The best restricted and unconditional results coincide. This is a stop
signal for further tuning of these markers, not evidence of bad distance.
The control and two receipts are `hd_sparse_marked_control.py` and
`hd_sparse_marked_a{5,6}_9927_float.json`.

## A concrete structural loss in the bound

The arity-five contribution has two problematic parts: the 1680 sets of
affine rank four and the 2688 affine bases of rank five. Their separate
best marked comparisons are approximately +405.059 and +399.539.
Fixing only one part cannot close the full arity-five class.

For an affine-basis input set A, exactly five outside supports have
weight four. Each contains255 vectors. The current weight-cap polynomial
instead assigns255 vectors to every one of the330 four-subsets of the
eleven outside coordinates. At equal output markers this term is66 times
too large. Moreover, the outside value determines all five inside output
coordinates, but the zero-dimensional fiber bound discards them.

The five actual supports have an explicit description. Label outside
positions by the ten edges of K5 and one additional position. The supports
are the five stars of K5; the additional position is absent. Solving five
small systems recovers the complete sixteen-coordinate outputs of these
projective rays. The helper checks their exhaustiveness against the exact
outside-image counts and checks the actual H* map.

`theory/hd_bit_syndrome_basis_rays.py` provides this exact local interface.
Its initial tests include an independent shortening calculation and all
1275 nonzero ray multiples for one input set. No full distance contraction
has yet used the recovered coordinates. The rank-four fibers must also
be retained; their nontrivial inside variation is not an affine-basis case.

A subsequent exact audit visited all2688 affine bases and all13440 rays.
Every representative has input weight five. Of the complete output masks,
13312 have weight nine and128 have weight eight. These counts concern the
minimum-outside-weight rays only, not the complete column enumerator.
The rank-four helper `theory/hd_bit_syndrome_rank4_rays.py` separately
retains the one-dimensional fibers and their inside coordinates.

This is the next proof repair to test. It changes neither the construction
nor its operation count. The correction must be made in the exact positive
kernel before taking moments, never by subtracting upper-bound matrices.

## Different construction options

An [independent local output shuffle](HD_LOCAL_OUTPUT_SHUFFLE.md) removes
the need to use worst-case labelled support caps. Its averaging argument
is valid, but it is a different ensemble with unmeasured runtime. The
layout review identifies extra SIMD and cache costs; it is not a free
reindexing of the measured encoder. The initial full-arity FLOAT screen
is saved separately. Its coarse treatment of inside fibers remains too
loose, especially at large arity.

The follow-up `theory/hd_local_shuffle_maclaurin_control.py` uses a scalar
translated-fiber bound and Maclaurin's inequality. It takes the pointwise
minimum of two valid evaluations before the group power; it does not take
a coefficientwise minimum of enumerators. The109-class FLOAT test improves
arity seven from+2107.64 to+928.07 and arity fifteen from+22564.28 to+10260.54.
The mixed comparison remains+10260.54, so this refinement does not close
the selected source class. Its receipt is
`hd_local_shuffle_maclaurin_allarity_9927_float.json`. Further tuning of
this coarse shuffled-fiber bound is parked.

Two more explicit alternatives are preserved rather than mixed into the
current certificate: a [nonbinary syndrome mixer](HD_NONBINARY_SYNDROME_PROPOSAL.md)
whose quiet code has minimum weight at least six, and an
[inner that exposes the refreshed state](HD_EXPOSED_STATE_INNER.md).
The former adds arithmetic. The latter can reuse existing state outputs,
but introduces a concrete cheaper cancellation case. Both need their
stated local stop/go tests before native or regional work.

Independent review cleared the positive activity and marker arguments.
The current67-test focused suite passes, including exact setup, schedule,
scalar-enumerator, and new local-shuffle checks. The certified fallback
manifest remains `accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.
