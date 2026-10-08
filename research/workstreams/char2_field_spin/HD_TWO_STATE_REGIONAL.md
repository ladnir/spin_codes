# Regional rank contraction and the next layout question

The common-order, two-prefix, two-additive-state candidate is rejected by
the separate interval-with-holes obstruction. Its reported K20 setup-failure
lower bound exceeds 0.9985. Therefore the work here stops short of a K20
certificate search and preserves two reusable proof components: a sharper
live-pivot bound and the correct joint regional contraction.

The exact finite obstruction is implemented separately in
`theory/hd_two_state_cut.py`. It concerns the common coordinate-major row
order. This note does not reproduce that proof or claim that the obstruction
extends to unrelated row orders.

## A sharper charge for an equation with a live pivot

Fix a current source space V over a field F of size q. A proof history
retains linear equations and a list of forms required to be nonzero on its
eventual witness. After imposing an equation, discard the history if any
required live form vanishes identically on the restricted space.

Let L be one such live form, and let M be another form on V. Sample
a from F*, independently of the current history, and impose M+aL=0.
For a rank mark tau>=1, the expected surviving rank factor is at most

    g_* = max(1/tau, 1/(q-1)).

To prove the claim, first suppose M and L are independent on V. Every a
causes a rank-one restriction. The form L remains nonzero on that kernel,
so the contribution is at most 1/tau. Other live constraints may only
discard additional histories.

Otherwise M=beta*L. Every nonzero equation M+aL=0 forces L=0 and kills the
history. Only the cancelling scalar a=-beta can survive, if it is nonzero.
That choice has probability at most 1/(q-1) and causes no rank loss.

This case split improves the earlier sum of the two alternatives. They
cannot both contribute for the same fixed pair M,L. For 1<=tau<=q-1,
the sharper factor is simply 1/tau. In particular, at tau=q-1 it avoids
an approximately twofold loss for every charged zero equation.

The lemma applies to returns in both pre-prefixes. It also applies to the
post-inner's zero output and zero state equations. Expose the output first:
fresh c is its pivot for live input, and fresh f is its pivot for zero input
with both old state components live. Then expose the two state equations
with fresh d1 and d2. The full historical live-form list remains part of
the induction; four current state flags only index a uniform envelope.

These charges are derived in the independent-coefficient experiment.
For the actual conditioned post-inner, divide the final bound by Z^N,
where Z=(q^2-3q+3)/(q-1)^2. The conditioning cost is paid once, not per
equation. At K20/GF65536 it is about 46.167 bits.

## Retaining the regional occupancy vector

Fix an outer [b,k] MDS code, R source rows, and N=bR. The precode consists
of two input-weighted prefixes in common coordinate-major order. Both
prefix states carry across coordinate-region boundaries. The post-inner
is the two-additive-state map from `HD_ADDITIVE_TWO_STATE.md`; its states
also carry globally.

For one prescribed source-row support A of size w, the shortened source
dimension is max(0,w-b+k). Feasible row supports are empty or have w>b-k.
Let D be the sum of these dimensions across source rows. Source supports
receive weight tau^D, not a count of field-valued messages.

Fix a complete precode support history. Each zero return contributes g_*,
and the history fixes the middle occupancy

    h = (h_0,...,h_(b-1)),    0<=h_c<=R.

Let P_h be the sum of these rank-envelope weights over all feasible nonzero
source supports and precode histories with occupancy h. The precode starts
in state00, follows the actual global scan order, and leaves its final
state free. The all-zero source support is excluded.

For an output mark 0<y<=1, let T_0(y) and T_1(y) be the post-inner's
four-state rank-envelope matrices for zero and live input, respectively.
Each entry sums allowed output/state support branches. A live output
contributes y; a charged zero equation contributes g_*. Local branches
forbidden by the conditioned determinant are omitted before paying Z^-N.

Define the regional matrix

    F_h(y) = [x^h](T_0(y)+x*T_1(y))^R / C(R,h).

Conditional on a precode support of size h in one region, its random route
is uniform over that region's C(R,h) output masks. The coefficient sums
exactly those masks. The local rank envelope is uniform over the live
forms assigned to them; no independent resampling of those forms is used.

Consequently, the post contribution for the full occupancy vector is

    Q_h(y) = e_00^T F_(h_0)(y) ... F_(h_(b-1))(y) 1.

The matrix product retains state carry between regions. Resetting each
matrix to state00 would change the construction. The order of occupancies
can matter even when their sum is fixed.

For a forbidden output-weight cutoff e, the resulting rank-envelope bound is

    Pr[distance<=e]
      <= tau^-1 * Z^-N * y^-e * sum_h P_h Q_h(y).

The final factor tau^-1 is valid because a surviving nonzero witness requires
at least one remaining source dimension. The bound counts support histories
and ranks, not individual field messages. Its validity requires the common
live-form induction described above.

There is no factor 1/C(N,sum(h)) in this formula. That factor belongs to
one uniform permutation of all N positions. It cannot replace the regional
normalizations. A tiny exact regression also shows that profiles (2,0),
(1,1), and (0,2) can give three different post contributions.

## The bounded exact proof object

`theory/hd_two_state_regional.py` computes this envelope with rational
arithmetic for small bR. It enumerates feasible source supports, contracts
the two pre-prefixes in their actual order, and retains every regional
occupancy. It then contracts the normalized post matrices without resets.
It does not enumerate field-valued messages or claim exact failure
probabilities; the local rank transitions are upper bounds.

For [4,2], R=2, q=4, tau=3, y=1/2, and cutoff3, there are 68 contributing
occupancy profiles. The exact un-clipped probability upper bound is

    9062328938209 / 30118144,

approximately 2^18.199. Thus this tiny reference is inconclusive. It is
a regression object for the correct contraction, not a successful code.

Seven tests pass. They verify the sharp live-pivot lemma for all GF4 forms
in dimension two and several tau values, the regional kernel's uniform
mask average, the distinction from global routing, exact source shortening,
precode boundary carry, and the complete small contraction.

```text
python -B research/workstreams/char2_field_spin/theory/hd_two_state_regional_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_two_state_regional.py
```

## Changing row geometry without a second full-record pass

The failed common-order cut uses many pre-inverse columns confined to the
same short interval of source-row check groups. Different row orders for
different outer coordinates can remove that property.

Let pi_c be a fixed permutation of the R source rows for coordinate c.
Write a newly computed outer value B_c(r) directly to

    scratch[c*R + pi_c^(-1)(r)].

The two pre-prefixes can then read scratch contiguously. This placement
requires no separate read-and-rewrite permutation pass. The later sampled
regional routing remains the only sampled full routing stage.

The cost is spatial store locality. The current outer computes four
adjacent source rows together and stores one 512-bit result per coordinate.
Unrelated pi_c generally split that result into four 128-bit stores at
different locations, potentially touching four cache lines instead of one.
No latency is inferred from this operation count.

Preserving common four-row groups to retain those wide stores restores the
grouped-check obstruction from `HD_REGIONAL_THREE_PROOF.md`. With pre-order
two, even common two-row blocks have a deterministic ceiling below40%.
Changing only block order is therefore not enough.

Computing the source rows separately in each coordinate's order is not an
equivalent free schedule. It gives up the shared multi-output GRS evaluation
or repeats input reads and computation. Likewise, computing prefixes out
of order requires later carry corrections across the stored records.
Direct placement followed by contiguous scans avoids those extra passes,
but it does not avoid the scatter-store cost.

A concrete next geometry is

    pi_c(j) = a_c*j+b_c mod R,

with R=2^15 and odd fixed a_c. Its inverse address uses one fixed modular
multiplier per coordinate, rather than a full permutation table. The last
two pre-inverse columns of each region need separate boundary handling.
Every other pre-inverse column touches a three-term arithmetic progression
with difference a_c among source-row groups.

Fixed affine orders do not make the three groups independent, and the
safe-column density cannot simply be replaced by alpha^3. Offsets alone
do not change the cyclic progression multiset. The following bounded
spectral check supplies a weaker, but justified, replacement.

## A fixed stride family removes the quarter-row safe-column cut

`theory/hd_two_state_regional_strides.py` freezes 64 distinct odd strides
on R=32768 rows. They were selected from eight bounded trial families.
The accepted family and the receipt are now fixed; no random row-order
ensemble is assumed in the result below.

Let P be the symmetric Cayley averaging operator with generators +/-a_c.
Its Fourier eigenvalues are

    lambda_m = (1/64) * sum_c cos(2*pi*m*a_c/R).

A 96-bit directed-Arb computation verifies every nonconstant mode, using
conjugate symmetry, and proves

    max_(m!=0) lambda_m < 323/1000.

The largest value is approximately0.322464375702732597 at mode1082.
All generators are odd, so the parity character has eigenvalue-1.
That negative eigenvalue does not invalidate an upper bound on internal
edges. It also does not imply an independent three-coordinate law.

For any row set U of density alpha, decompose its indicator into its
constant part and its orthogonal part. The spectral bound gives

    (1/(64R)) * sum_(c,x) 1_U(x)*1_U(x+a_c)
       <= alpha^2 + (323/1000)*alpha*(1-alpha).

Every three-term progression contained in U contributes an internal edge.
Thus the same expression bounds structurally safe interior pre2 columns.
The statement holds for every U, including unions of low-bit fibers.
At most two columns per region need a separate global-boundary allowance.

For the particular quarter-row cut, a=8192 and alpha=1/4. The structural
count is at most floor(32260096/125)=258080, compared with262144 checks.
There is one coefficient exception: the middle entry of a two-prefix
inverse column vanishes when the two adjacent second-prefix weights agree.
The endpoint coefficients never vanish. Hence actual safe columns are
bounded by the structural count, plus128 boundary columns, plus the total
number C of these middle cancellations.

For independent nonzero GF65536 second-prefix weights, C has distribution
Binomial(N-1,1/65535). The adjacent-equality indicators are jointly
independent: any selected edges in the global path form a forest and impose
one independent equality per edge. Exact integer comparison proves

    Pr[C>=1024] <= C(N-1,1024)/65535^1024 < 2^-3649.

Outside that event, every row set of size8192 has at most259231 safe
columns, fewer than its262144 checks. The interval-with-holes cut requires
gain greater than the check capacity, while gain is at most the total safe
count. Its sufficient condition therefore cannot hold for this row-set
size and this placement. This disables that cut; it does not prove distance.

There is also a cheaper next candidate. With just one pre-prefix, an
interior inverse column is exactly a Cayley edge, with nonzero endpoints.
There is no middle-cancellation exception. Allowing64 boundary columns gives
at most258144 safe columns for every quarter-row set, deterministically.
The intended next topology is therefore

    GRS[64,32] with these fixed affine placements
      -> one weighted prefix
      -> independent uniform regional permutations
      -> two additive states with varying observation.

The post-inner has a dense inverse of tail rank two, so the old bound for
two bidiagonal inverses does not apply. No distance lower bound follows.
In particular, the quarter-row calculation does not eliminate cuts at all
other row-set sizes or more general check subspaces.

Six additional tests pass: affine-address bijectivity, the directed spectrum,
an exact small all-subsets edge bound, independence of path equalities,
the inverse-entry cancellation identity, and the finite quarter-cut counts.

```text
python -B research/workstreams/char2_field_spin/theory/hd_two_state_regional_strides_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_two_state_regional_strides.py --certify
```

The next turn should test remaining cuts and develop a source-support bound
for this fixed graph. Native performance remains unknown. The additional
fixed address map can be fused into outer stores without a second full-record
permutation pass, but its scatter-store cost still needs measurement.

No native implementation, build, performance benchmark, or remote command
was run. Only the new regional proof note and its new helpers were written.
