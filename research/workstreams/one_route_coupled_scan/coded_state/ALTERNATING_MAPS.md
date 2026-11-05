# Alternate the expansion map without enlarging the state

The fixed expansion map permits a nonzero state to retain the same 63 holes
through an arbitrarily long zero-output run. Alternating two evaluation sets
removes that particular loss, without adding full-field products. This is an
inner-only research variant, not a whole-code near-half certificate.

## Construction and setup

Keep the fixed RS[256,128] outer and the 256 independent uniform regional
permutations. Work over F=GF(2^128), with its AES subfield E=GF256, at K=2^20.
Use batches of 128 symbols, state of 64 symbols, Q_0=0, and a free terminal
state. Let e_0,...,e_255 be the existing Cantor ordering of E. Put

    Y_0={e_0,...,e_127},  Y_1={e_128,...,e_255},
    Z={e_128,...,e_191}.

For Q in F^64, let p_Q be the unique degree-<64 polynomial evaluating to Q
on Z. Define A_i Q by evaluating p_Q on Y_i in that ordering. Retain B:
interpolate a degree-<128 polynomial on Y_0, then evaluate it on Z.

For each batch t=1,...,16384, independently sample the entries of D_t from
F and the entries of R_t from E minus zero. Both matrices are diagonal.
All coins are independent of the route, sampled once, and fixed for every
message. Define epsilon_t=(t-1) mod 2 and use

    Y_t = X_t + D_t A_(epsilon_t) Q_(t-1),
    Q_t = Q_(t-1) + B R_t Y_t.

The old [column-scaled variant](QUIET_BRANCH.md) uses A_0 at every step.
Only this expansion schedule changes. No claim B A_1=I is required.
Given an output stream and R, the state is determined independently of D.
Consequently the inverse-input coordinate-cube distribution is unchanged
except for its actual feedback support.

## Two complementary observations of an unchanged state

Each A_i is an [128,64,65] MDS map over F. Their concatenation evaluates
one degree-<64 polynomial at 256 distinct points, so it is [256,64,193] MDS.
For every nonzero Q,

    wt(A_0 Q)+wt(A_1 Q) >=193.

This follows from the polynomial root bound over F; it is not a GF256-only
empirical property. A zero-output batch preserves Q. On a run of ell
consecutive feedback batches having the same nonzero incoming state, pair
consecutive observations, regardless of the starting phase. The number of
feedback holes is at most

    63 ceil(ell/2).

Thus two such batches have at most 63 holes total, instead of the separate
worst-case bound of 126. An unpaired batch still permits 63 holes.
Nonzero outputs can change the state, so this bound cannot be applied to
arbitrary adjacent live batches.

## A joint rank bound when the state changes

Fix two opposite phases i and 1-i, hole sets H_0,H_1 of sizes h_0,h_1<=63,
and an intervening output support W of size w. Fix its nonzero R entries.
Treat the old state Q and supported output y_W as variables in one system:

    A_i[H_0,:] Q = 0,
    A_(1-i)[H_1,:] (Q+B[:,W] R_W y_W) = 0.             (1)

Write h=h_0+h_1 and r=min(64,w). The rank of (1) is at least

    min(64,h) + max(0,h+r-128).                         (2)

Proof: the coefficient matrix on Q has rank min(64,h), since the two
evaluation sets are disjoint. When h>64, the row spaces of the two
observations intersect in dimension d=h-64. Eliminating Q leaves equations
L B[:,W] R_W y_W=0, where L has row rank d. The MDS property of B gives
rank B[:,W]=r. The matrix-product rank inequality gives rank at least
max(0,d+r-64), proving (2). When h<=64, the Q columns already give (2).

The exact rank after eliminating Q equals the joint rank minus min(64,h).
Only that eliminated rank can be charged directly to y_W without an
additional argument about earlier outputs. Counting the Q rank again
would double-count constraints already imposed on the incoming state.
For w>=64 the joint rank is exactly h; for w=0 it is min(64,h).
With zero output and h>=64, (1) forces Q=0 and is impossible for a live run.

This is a composable local constraint interface, not yet a summation over
all output supports and trajectory patterns.

## Tests separate rank deficiency from a live-to-live obstruction

The exact GF16 analogue uses state4 and batch8. A_1 B has fewer deficient
small minors than A_0 B, but it is not totally nonsingular: rows (0,4) and
columns (2,6) give a rank-one 2-by-2 minor. Therefore generic full-rank
substitution remains invalid. In the actual GF256 maps, every A_0 B column
has weight65 and every A_1 B column has weight128. Column weights alone do
not establish the needed multitime ranks.

The two-step GF16 test prescribes three holes in each expansion with one
intervening nonzero output symbol. Among 25088 fixtures per starting phase,
24640 have full rank5 and 448 have rank4. Importantly, all 448 deficient
phase0 fixtures force the next state to zero. All 448 deficient phase1
fixtures force the old state to zero. None admits both states nonzero.

An exhaustive projective check then fixes the nonzero output amplitude to
one, enumerates every state, and excludes activation/reset endpoints.
Across 1048544 live-to-live pairs, the minimum combined feedback weight is
11 in either phase, versus 13 with zero intervening output. Witnesses are
Q=(0,1,4,7), output coordinate1 for phase0, and Q=(0,0,2,10), coordinate0
for phase1. This is toy-field evidence, not a claimed weight bound at
state64 or over F. It motivates a sparse-update split enumerator that
excludes zero-state fibers before taking a worst rank.

The earlier half-count diagnostic also remains useful: among449 profiles
classified by left/right cardinalities and intersections, two contain
different exact ranks. Those statistics alone do not determine AB rank.

## Partial proof transfer and cost scope

The following arguments transfer analytically to this variant:

- The sparse-source and regional-episode certificates through64 outer
  blocks use per-step MDS holes, B's kernel distance, and exact uniform
  refresh. Each phase has the same MDS parameters. Every live Q gives at
  least65 random feedback coordinates, whose B R columns span F^64.
- The entirely quiet certificate is unchanged: A is multiplied by zero.
- The mixed quiet-rank argument is unchanged. Its forward state-hole
  bounds allow a predetermined MDS map at each time.
- The output-interval certificate uses only support containment and the
  per-live-batch lower bound65, so it also transfers for every fixed R.

The older numerical certificates were not recomputed for a new parameter
set: these are analytic transfers of the same inequalities and constants.
They do not cover all mixed trajectories.

The reference implementation uses an inverse64 transform on Z, zero padding,
and a forward128 transform on the selected half. A_1 Q begins with Q itself,
which may permit further pruning. Both phases use128 full-F products per
batch and the same D/R tapes as the previous variant. There is no extra
global route or random matrix. No native timing or net speedup is claimed.

## Reproduce and next decision

Run serially, with Python -B and bounded memory/time:

    python -B coded_state/test_alternating_maps.py
    python -B coded_state/test_alternating_rank.py

The first checks65535 nonzero GF16 states, both run parities through length7,
864 exact GF4 refresh fixtures, actual-size FFT/Lagrange agreement, column
weights, and16 recurrence/inverse steps. The second checks17768 square
minors,4096 joint-rank fixtures, two exact eliminated-state fiber counts,
the50176 maximal-hole fixtures, and the exhaustive live-to-live pairs above.
The final-source runs took approximately1.5 and2.1 seconds, respectively,
on Peach. These are checker runtimes, not encoder timings. No files are
generated.

Next derive the sparse-update split enumerator in the actual64-state maps,
retaining zero-state fibers and the same global message/output constraints.
Do not extrapolate the toy minimum11, or count the state as fresh randomness.
If that enumerator remains too weak, compare a targeted map change against
the specific deficient class before adding more arithmetic. The full
>=49%/40-bit certificate remains open.

The [one-half update analysis](HALF_UPDATES.md) now proves an actual-size
feedback-pair floor150/256 for every live-to-live weight-one update, with
further bounds through weight12. It also gives an exact zero-fiber
subtraction and a local dimension upper bound. Three complete restricted
one-episode classes now pass at >49.1%; their placement sums are checked
by alternating_interval_certificate.py. These do not close the global
mixed-trajectory sum. Its full-length stress test calls for retaining the
actual feedback support in the outer-shortening bound, not just its interval.
