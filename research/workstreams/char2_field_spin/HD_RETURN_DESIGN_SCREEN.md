# What limits the cheaper emission design

The certified checkpoint remains approximately 5.7--5.8 ms, with 40%
GF65536-symbol distance and 51.04 bits of setup-failure margin. The
approximately 5.3-ms common-refresh candidates remain uncertified. This
note records further proof screens. The rank-two variants at the end now
have native timings; the earlier alternatives remain untimed.

All screens use K=2^20, rate 1/2, state16, packet width64, and the existing
RS[256,128] randomized outer and regional route. A negative single-row
screen is not a complete certificate. Directed arithmetic, all 8192 row
counts, the union bound, and setup-density losses remain separate steps.

## The return moment is still the obstruction

For the timed scattered-feedback candidate, the corrected natural-degree
rank bound gives logarithmic screens -133.99 and +11182.47 at 512 and2048
active rows. Averaging labelled occupancies improves these to -136.62
and +11042.00. These are valid local upper-bound constructions evaluated
in FLOAT arithmetic; neither closes the middle-density class.

As an obstruction diagnostic only, setting the return-to-zero matrix
entry to zero gives -5505.23 at row2048, even with a global marker instead
of exact regional extraction. This artificial operator is not a valid
upper bound. It indicates that improving the return bound could matter;
it is not evidence that the actual return probability is zero or small.
The diagnostic cannot use the exact-region routine, which rejects its
structurally zero matrix coefficient.

## A different local polynomial basis

One proposal uses the four Lagrange polynomials for nodes16,...,19 in
both local RS maps, evaluated at nodes0,...,15. Keep the feedback-coordinate
transpose. Every square minor of the resulting sixteen-by-four matrix
is nonzero. To see this for k selected basis polynomials, factor out their
common 4-k roots outside the evaluation set. The remaining polynomials
form a basis of degree below k, so any k evaluation rows are independent.

Products of the local emission and feedback polynomials have degree at
most six. Their sums over the sixteen-point additive subspace vanish,
so C_sc A=0 still holds. The quiet and uniform-state moments are unchanged.
Unlike the timed degree-basis map, this new basis permits any ordering
of the input-block pivots in the inactive-zero rank argument. One may
therefore take the largest resulting lower bound.

Exact tests check all 4844 local square minors, C_sc A=0, 600 concrete
restricted ranks, and the one/two-input return caps. The row2048 FLOAT
screen is still +10651.32. This basis change does not justify a native
implementation on current proof evidence.

## A fresh scalar in the emission

Another proposal keeps the timed maps but samples an independent nonzero
scalar gamma for each packet. Compute Y=X+gamma A Q and
Q_next=T(Q+C_sc X). The scalar lies in GF256 or GF65536. Both choices
preserve GF65536-linearity. They require a new setup distribution and
additional arithmetic; neither has been timed.

At a return, the output is (I+gamma A C_sc)X. For active output-zero
positions D and inactive output-zero positions E, let c=|D| and U'=U-D.
At gamma=0, the associated constraint matrix has rank

    c + rank((A C_sc)[E,U']).

A witnessing minor is a polynomial in gamma with nonzero constant term.
Its degree is at most min(c,16), because the perturbation has rank at
most16. Thus at most min(c,16) nonzero scalars lose this rank. Outside
those exceptional scalars, subtract the quiet mass
q^(-c-rank(C_sc[:,U'])). On the exceptional scalars, discard the active
zero constraints and retain the inactive nonquiet probability.

The implementation sums the positive zero-set expansion. It uses the
natural-degree grid-rank lower bounds for the generic term and retains
the existing inactive cap for exceptions. Exact GF4 tests cover every
input support and output-zero set in a small model. A separate test
checks the one-input GF65536 return average for both scalar domains.

At the fixed row2048 witness from the corrected grid screen, GF256 and
GF65536 scalars give +11011.58 and +10757.35. Optimizing the latter gives
+10756.49. The determinant bound is insufficient; no scalar variant is
promoted.

## Applying the refresh only to the incoming state

The recurrence Y=X+A Q, Q_next=T Q+C_sc X removes the same arithmetic
as common refresh, but changes which exceptional event must be bounded.
The fresh T acts after emission. For any fixed nonzero Q, its output has
point probability at most b/(q^16-1), where b is the diagonal-Cauchy
density factor. Returns are therefore easy to bound. Births are no longer
flat and must be tracked separately.

Track zero mass Z, nonzero density coefficient B, and arbitrary nonzero
mass R. Write d=q^16-1. For a fixed potential profile, K,T_j,V retain
their previous meanings. For the arbitrary component use
L=T_j z^max(0,13-max_i l_i). Before the common density factor b, a local
comparison matrix is

    [[K,       0,   T_j-K],
     [V/d,     V,   0],
     [L/d,     L,   0]].

Both return entries are zero on an empty packet. A birth can alternatively
use the density cap B_birth=d q^-rho T_(j-rho), where
rho=sum_i min(4,l_i). The screen chooses one of the two birth envelopes
by a fixed rank threshold. This choice changes the bound, not the code.

When averaging the within-packet positions and block labels before paired
Jensen, L can use the hypergeometric puncturing moment of a weight-thirteen
word. This is a uniform bound for each fixed incoming nonzero state;
current packet positions are sampled independently of that state in the
conditional comparison experiment. It is not a pointwise bound for a
fixed support. Exact GF4 enumeration checks the pointwise operator and
birth density for every potential support and incoming state.

The fixed row2048 witness gives +25189/+29779/+53676 for pointwise birth
thresholds8/12/16. These are not optimized results. They show why moving
randomness from feedback to state is not automatically a proof improvement.

## Restoring global emission

Keeping common refresh but replacing the block emission by RS[64,16]
gives a different return map. Its I+C_sc A has exact rank16, so I+A C_sc
is invertible. For a rho-dimensional syndrome space and e inactive output
zero constraints, the combined rank is at least max(0,rho+e-16). Sixteen
zero constraints force the entire syndrome to zero. Summing these ranks
and marking the j-rho input-fiber dimensions gives a local return cap.
The global minimum-weight cap and the information-set cap also apply.

Exact small-field tests check this subspace bound and the concrete
invertibility check. The optimized row2048 screen is +13672.39. Alternating
block/global emission gives +12573.05 at the original fixed witness.
Thus the current global-emission return bound also remains insufficient.
No native timing has been attached to either changed construction.

## Rank-two refresh with larger feedback blocks

Rank-two shared refresh retains the cheaper four-block emission but pays
for one extra random dot product. It uses the setup and recurrence in the
[block-emission note](HD_BLOCK_EMISSION.md). For two RS[32,8] feedback
blocks, a sorted four-block emission profile has three possible pairings.
The corrected operator takes the largest quiet moment for its zero-to-zero
entry and the smallest quiet moment for the complementary entry T-K.
It minimizes the feedback rank separately. These choices bound every
entry for every labelling; they need not describe one common labelling.

An audit found that the first implementation incorrectly used max(K)
also in T-K. The profile (0,0,8,8) has distinct quiet moments for different
pairings, so that entry was not an upper bound. The first coverage search
was stopped and its partial coverage is withdrawn. The corrected FLOAT
and exact providers pass an entrywise regression against every labelling
of that profile and one hundred further sampled profiles. The existing
alternating-feedback certificate does not use this grouping and is unaffected.

Before this correction, alternating four/eight-check feedback gave FLOAT screens
-134.43/-284.04/-1177.27/-1158.84/+868.24 at rows256/512/1024/2048/3072.
Using eight-check feedback in every packet gives
-307.11/-2874.80/-1247.23/-22135.78 at rows512/2048/3072/4096. Thus the
alternating version still fails a sampled class, whereas the all-eight
version merited full coverage. These values are historical diagnostics,
not proof evidence; the corrected provider must recheck every class.

An exact rational provider now supplies packet matrices to the directed
verifier. The earlier row3072 witness was
z=106595/7203386, alpha=420853/6283999, p=4391671/9398723. Its pre-correction
receipt is withdrawn along with the interrupted search. The provider
retains all profiles and uses twelve-bit upward mantissa rounding.
The full density loss for rank-two refresh is 12.26301241396 bits.
The corrected 256-bit search now covers all 8192 active-row counts. It
re-evaluates the frozen alternating-check witnesses as proposals, without
importing their guarantees. Its 133 accepted rational witnesses are pinned
in [the corrected replay module](theory/hd_rank_two_block_witnesses.py).
The ideal union logarithm is below -63.25965982455041; after the
12.26301241395966-bit density loss, the failure logarithm is below
-50.99664741059075. A separate 384-bit replay confirms all 8192 classes
and the same bound. Modes38/41/43 therefore have a complete certificate.

Native mode37, `bucket48-rank2`, alternates feedback widths; mode38,
`bucket8-rank2`, uses eight checks throughout. Both use four RS[16,4]
emission blocks. The hot circuit retains unrolled GFNI operations and
reuses the old mixer's unused coefficient array for the second row of V.
No extra hot allocation or indirect dispatch is introduced. A K20 setup
test caught inappropriate reuse of rank-one validation; setup now checks
the joint two-by-two determinant only, and that seed is a regression test.
All 38 variants pass Release and ASAN/UBSAN validation.

Mode37 gives 5.597/5.690 ms in the first sequence. Its frozen binary is
`rank-two-alternating-mode37-b6a65803`, SHA256
`b6a65803c3fd24a69e7d31e00cbab0205c55823ecdeccdd2b502bcac37c71688`.
Its checksum is 3626870394749882516. The final control in that sequence
drifted to6.439 ms, so it is not used to claim a large relative win.

Mode38 gives 5.684/5.575 ms against 5.825/5.749-ms certified controls.
A reversed candidate/control/candidate sequence gives
5.583/5.742/5.579 ms. Report approximately 5.6 ms, about a 3% improvement
in the stable reversed sequence, not a sub-5-ms result. Its checksum is
7812899659201065080 and its fixed-map digest is 6742237189988151056.
The measured executable is frozen as `rank-two-eight-mode38-25334b78`,
SHA256 `25334b787171b442a4500989045dca68bcd435c0cb97809ad38ff880662e555d`.
All timings include the full encoder, exclude setup, and run sequentially
under both host benchmark locks.

## Smaller batches and omitting the second refresh row

Two further FLOAT screens use the corrected matrix grouping. Halving the
packet and state sizes gives width32/state8, two RS[16,4] emission blocks,
global RS[32,8] feedback, and rank-two shared refresh. Optimized bounds at
rows256/512/2048 are approximately +4112/+6416/+15051. These failed bounds
do not support a smaller-state implementation at the target.

Keeping width64/state16 and eight-check feedback but reverting to rank-one
refresh gives approximately +490/+606/+153 at rows256/512/2048. The two
sparse checks use regional coefficient extraction; row2048 uses the global
marker bound. Thus stronger feedback alone does not presently justify
removing the second random dot product. Neither screen is a counterexample
to the construction's actual distance, and neither has a native timing.

## Using eight-check feedback less often

This separate candidate retains the certified global RS[64,16] emission
and rank-one shared refresh. It repeats feedback ranks (4,4,4,8), instead
of (4,8). Only one quarter of its packets use eight-check feedback.

For a fixed packet support U, let P_U^(h) be the positive two-state
operator from the block-feedback proof, with local check count h. Define
the matrix polynomial B_(a,b)(x) by summing
x^(|U|+|V|) (P_U^(a) P_V^(b))^[alpha] over all support pairs (U,V).
Here [alpha] denotes entrywise powers, with 0<alpha<=1. Use the period
polynomial B_(4,4)(x) B_(4,8)(x). This is an upper bound because
(AB)^[alpha] <= A^[alpha] B^[alpha] entrywise for nonnegative matrices.
The regional coefficient extraction uses its 32nd power, since a region
has 8192 positions and each period consumes 256 positions. The outer
measure comparison and shared-refresh density factor are unchanged.

The FLOAT screen retains every occupancy profile before upward grouping.
At active-row counts 2048/2560/3072, its optimized regional logarithms are
approximately -1309.38/-1162.71/-3145.40. A coarser screen that merged
profiles above degree 20 had failed; those failures did not establish an
obstruction to this construction. Repeating B_(4,8) twice exactly matches
the original global bound, checking the composition's normalization.

The first 256-bit directed receipt uses z=139402/4423823,
alpha=459633/7739926, and marker p=309254/795265. It certifies active-row
counts 2560 through 2683 with ideal union logarithm below -80.63586432043291.
These are partial classes, not a whole-code certificate. A period with
three eight-check packets out of eight also has a negative row3072 screen.
The quarter-frequency schedule now has native mode45, measuring about
5.43--5.45 ms against nearby 5.49-ms controls; the other schedule is untimed.
The directed period evaluator
and four composition regression tests are recorded alongside the screen.
The (4,4,4,8) complete-coverage search has now closed at 256-bit precision.
It recomputes the old alternating witnesses under the new period operator,
then optimizes uncovered classes. Its 119 accepted rational witnesses
cover all 8192 active-row counts. The ideal union logarithm is below
-68.92219588903863. The unchanged 12.26300140707-bit density loss leaves
more than 56.65919448196863 bits of whole-code margin at 40% distance.
The [pinned replay](theory/hd_sparse_feedback_witnesses.py) reconstructs
every directed period polynomial without optimization. A separate 384-bit
run passes all 119 witnesses and covers all 8192 classes. Its final union
logarithm is below -56.65919448196863. Optional worker processes evaluate independent witnesses;
the final union and coverage checks consume them in fixed witness order.
This is proof arithmetic, not concurrent encoder benchmarking.
