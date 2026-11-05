# Pre-H Toeplitz pairs: fixed-mask fractional control

The earlier scalar means do not test the main proof obstruction: rare
output masks are powered before they are combined with the inner law.
This control therefore composes the complete positive mass of each
labelled mask before taking its fractional power. It uses no regional
provider and does not claim a distance certificate.

## Scope and calculation

Use the construction and iid-NZ majorant J from
[the pre-H comparison note](HD_TWO_SUM_PRE_TOEPLITZ_PAIRS.md).
Fix A={0,...,7}. For each u in {0,...,4}, choose a fixed matching with
u pairs inside A, 8-2u crossing pairs, and u pairs outside A. These five
matchings are separate one-column controls. They are not independently
chosen columns of a newly sampled schedule.

For one target support T after H, define

\[
 C_{A,M}(T)=\sum_B
 \frac{L_{A,M}(B;1)}{255^{|B|}}
 N_{\rm NZ}(|B|,|B\cap T|,|T|).
\]

Here L is the raw intermediate-support mass after J on uniform input
contained in A. The exact integer N_NZ counts inputs with support exactly
B and H-output support T. Its dependence on the three sizes follows
from the established closed kernel, not from an assumed permutation
symmetry of H.

The helper forms the positive pair polynomial in |B| and |B intersect T|.
It clears each pair's rational denominators, multiplies integer
coefficient tables, and returns C as an exact rational. Only after this
composition does it evaluate

\[
 F_u(w)=\sum_{|T|=w}C_{A,M_u}(T)^{16/5}.
\]

Target multiplicities are outside the power. Grouping identical pair
types gives 165, 756, 1260, 1000, and 225 target classes for u=0,...,4.
Their exact multiplicities partition all 65536 labelled targets in each
case. The old control replaces C(T) by the old fixed-A count N_A(T).
Both zero-mask masses equal one exactly.

The full seventeen-entry vectors are the primary result. Auxiliary
scalar summaries use the three frozen marks 691/700, 172/175, 171/175,
their arity-eight activity logs ell, and the standalone mark 1/2. The
reported formula is explicitly

\[
 \sum_w F_u(w)z^w2^{16\ell\mathbf1_{w>0}}.
\]

The frozen activity log already acts after alpha, which explains 16ell
at this power. No RS floor is enforced by these summaries. Separate
totals for weights 1..4 and all nonempty masks prevent the zero mask from
hiding a change in the nonempty contribution.

This fixes A and sums over targets T. The old arity weights instead fix
T and sum over constituent sets A. A matching breaks the symmetry that
would identify those operations. Consequently, neither F_u nor its
scalar summaries can be substituted for the old W_a without a new
averaging argument. RS normalization, source counts, projective scaling,
and the inner operators are absent here.

## Completed local gate

The parent ran all five cases in 0.702656 seconds. The table gives
log2(old sum)-log2(new sum); positive values mean a smaller comparison.

| Internal pairs u | Weights1..4 | Weights5..8 | Weights9..16 | All nonempty |
|---:|---:|---:|---:|---:|
| 0 | 52.558065 | 96.342922 | -0.001511 | -0.001511 |
| 1 | 29.144530 | 72.514420 | -0.019225 | -0.019225 |
| 2 | 6.093565 | 48.615458 | -0.036728 | -0.036728 |
| 3 | 2.200945 | 24.564494 | -0.053947 | -0.053947 |
| 4 | -0.015340 | -0.071433 | -0.070877 | -0.070877 |

The old logarithmic totals for these three weight ranges are
57.256793, 153.455448, and 204.509386. High-weight masks dominate the
unrestricted total. Consequently, every tested scalar summary is slightly
worse despite the large low-weight gains. At z=1/2, the old total is
188.509386 and the five new totals are 188.510897, 188.528611,
188.546114, 188.563333, and 188.580263. The three near-one marks and
their activity variants have the same gain signs and the same displayed
six-decimal losses. The artifact retains all seventeen coefficients and
every scalar value, rather than discarding these failed totals.

This isolates a useful effect but does not settle the inner comparison.
The next calculation must weight these masks by the actual regional
operators. It must not read a 52.6-bit mask-moment gain as a 52.6-bit
whole-code failure margin.

The result artifact is `hd_two_sum_pre_toeplitz_fractional_gate_float.json`,
SHA-256 `65f670512f9f5c50d7091fbe4e547ad02c2bd8c994376fde20228fb6676bc8fb`.

At u=4, every pair is wholly inside or outside A. The exact pair maps
therefore preserve uniform input on F_256^A, so the exact column law
cannot improve. Its iid majorant can only add a comparison tax. For this
fixed A, the proposed XOR schedule (directions 1..15, then 1) contains eight
u=0 columns and eight u=4 columns. A future schedule calculation must
retain that dependence rather than selecting the best case per column.

## Correct all-row-set recovery for the next operator gate

There is an exact incidence identity that avoids interchanging those
sums without justification. Fix a matching M. Its automorphism group
permutes the eight pairs and swaps their two coordinates. That action
preserves the pair comparison and the three-size H count, hence C(A,T).
It is transitive on input sets of a given size a and internal-pair count u.
It is also transitive on targets of a given weight w and full-pair count f.

Let

\[
 n_a(u)={8\choose u}{8-u\choose a-2u}2^{a-2u},\qquad
 n_w(f)={8\choose f}{8-f\choose w-2f}2^{w-2f}.
\]

For a representative A_u, let S_u(w,f) sum C(A_u,T)^(16alpha) over
the target orbit(w,f). Counting the same incidences by inputs and by
targets gives, for each fixed T in that orbit,

\[
 \sum_{|A|=a}C_{A,M}(T)^{16\alpha}
 =\frac{\sum_u n_a(u)S_u(w,f)}{n_w(f)}.
\]

This is a fixed-target moment. It produces the correct matching-dependent
weight

\[
 W_M(T)=\left(
 \frac{\sum_{|A|=8}C_{A,M}(T)^{16\alpha}}
 {255^{64\alpha}}\right)^{1/16}.
\]

Hölder applies to the same row set A over the sixteen columns; it does
not sample A again. No extra factor binomial(16,8) belongs outside W.
Unlike the old weight, W depends on both |T| and its full-pair count in M.

The operator controller is restricted to r=18, source stratum
F=1, and the homogeneous arity-eight cohort. It retains all sixteen
XOR matchings chronologically and use the frozen regional bank at
z=691/700, alpha=1/5. The saved regional matrices already incorporate
alpha. Four checked affine regional-entry envelopes give

\[
 R_{st}(m)\le A_{st}B_{st}^m.
\]

Fix an activity marker y>=1. For a binary state path sigma across
sixteen regions, with endpoints sigma_0=s and sigma_16=t, set

\[
 G_j(\sigma)=\sum_T W_{M_j}(T)
 y^{\mathbf1_{T\ne\varnothing}-9/16}
 \prod_{i\in T}B_{\sigma_i,\sigma_{i+1}},
 \qquad
 K_j(s,t)=\sum_{\sigma:s\to t}
 \left(\prod_i A_{\sigma_i,\sigma_{i+1}}\right)G_j(\sigma)^{18}.
\]

The groups' row sets remain fixed in the underlying experiment; their
choices have already been charged by Hölder. A sequential binary tensor
transform evaluates G for all 2^17 paths using at most 2^17 working entries
and O(16*2^17) arithmetic per matching. It retains state order rather
than sorting transition markers, which would be invalid for these
matching-dependent weights.

The final single-class bound is

\[
 N(18,F=1)\,(z^{-734003}/255)^{\alpha}
 e_0 K_0K_1\cdots K_{15}{\bf1}.
\]

Here N(18,F=1) counts indexed source sets with eighteen groups whose
larger full-pair count in the two route matchings is one. This exact
source count is outside alpha. RS normalization and the nonempty
column-floor charge occur inside W and G respectively, exactly once.
The controller deliberately uses the ordinary affine envelope;
it does not retain the refined empty00 branch or exact off-diagonal
zeros. It must compare the old and new outer laws under that same
envelope and marker schedule. A successful result covers only this
homogeneous cohort and source stratum, not all cohorts or the whole code.

The parent completed this gate in 2.306020 seconds. All three fixed
activity schedules fail, although the new outer comparison improves the
same-envelope bound by about 830 bits.

| Activity schedule | Old log2 upper screen | New log2 upper screen | Improvement |
|---|---:|---:|---:|
| Zero log | 6161.631728 | 5331.172282 | 830.459445 |
| Frozen scale1 | 7294.344857 | 6463.862416 | 830.482441 |
| Frozen scale2 | 8427.061584 | 7596.579098 | 830.482486 |

These are failure-upper-bound screens, not actual failure probabilities.
Their positive values provide no distance certificate. The earlier
empty-cut comparator retained information this ordinary affine envelope
drops, so the results do not reject the construction. A cheap next check
can isolate the positive all-zero-state path after restoring that empty
branch. If that single comparator component still exceeds the target,
the full comparator cannot close at the same fixed witness.

The artifact `hd_two_sum_pre_toeplitz_schedule_gate_float.json` has SHA-256
`64b5aafb91ca7d2b11584ad9896836d9bd3e3239b5c4a36c82bb8921688d7e35`.

The proposed native product count remains 47, with no implementation
claim. Normalize the nonzero output-diagonal entry d_0, not a Toeplitz
coefficient that might be zero. Replacing every pair matrix T_j by
d_0*T_j and D_out by D_out/d_0 preserves the product through linear H.
The joint Toeplitz law is invariant under that common nonzero scale.
The resulting count is 24 variable pair products, 15 output-diagonal
products, and eight fixed H products.

## Validation and reproduction

Ten new tests pass; the combined pair-law, pre-H, and fractional suites
have 25 passing tests. Checks include labelled intermediate-mask
enumeration, old-diagonal equality, positive normalization, every target
weight's exact binomial count, C(empty)=1, unpowered multiplicities,
scalar first-moment recovery, and partial-artifact budget handling.

From `theory`:

```text
python -B -m unittest hd_two_sum_toeplitz_pairs_tests hd_two_sum_pre_toeplitz_pairs_tests hd_two_sum_pre_toeplitz_fractional_tests
python -B hd_two_sum_pre_toeplitz_fractional.py --max-seconds 30 --output NEW_OUTPUT.json
```

The CLI permits only the five fixed cases, alpha1/5, and input marker1.
It reads frozen scalar metadata without building regional operators.
It preserves completed cases if the thirty-second guard expires.

SHA-256 pins:

- Helper: `a1e2b46ce992b5ab433b2c6bbe57e9407df338b792306be548cd72539aa2ab7b`.
- Tests: `798a261e57bb5ccb077705d7fa8e8166d3eac018ef36749ad4e812e1e9cf5ea8`.
- Frozen metadata artifact: `5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf`.

The exact rational mask calculation does not make the subsequent
binary64 fractional powers directed. Any successful local diagnostic
still requires a fixed-A, fixed-schedule operator comparison and a full
class union before the distance target can be claimed.

The isolated operator helper `hd_two_sum_pre_toeplitz_schedule.py` has
seven tests. They compare incidence
weights against direct fixed-target row-set sums, recover the old
weights, and compare every small chronological tensor path with brute
force. The combined suites have 32 passing tests. The parent independently
ran the seven tests before the single operator gate.

- Schedule helper SHA-256: `2c40606f0efc615a9fc8003d839a0bcb61c4050bfb3e1a6188d2d7ddf6199305`.
- Schedule tests SHA-256: `5e6e2991a673a1e0e595af6a069d0b01edd226ff554883c97afd2189a0590871`.

## Next component check: restore empty00 without enumerating augmented paths

For the frozen slope B_00, define

\[
 c=\max_{1\le m\le18}R_{00}(m)B_{00}^{-m}\le1.
\]

The positive refinement is
R_00(m)<=B_00^m[c+(1-c)1_{m=0}]. Along an all-zero state path, a forced
empty slot requires every group's mask to omit that coordinate. For an
allowed set U, let

\[
 G_U=\sum_{T\subseteq U}W_M(T)B_{00}^{|T|}
 y^{\mathbf1_{T\ne\varnothing}-9/16}.
\]

The refined all00 column component is exactly

\[
 Q_{00}=\sum_{U\subseteq[16]}
 c^{|U|}(1-c)^{16-|U|}G_U^{18}
\]

within this positive majorant. No subtraction of upper bounds occurs.
A weighted subset-zeta transform evaluates all G_U. The cut weights
already enumerate every labelled U, so they receive no extra binomial
factor. Because W depends only on the target's matching orbit, relabelling
coordinates shows that Q_00 is identical for every perfect matching.
One transform per old/new law therefore suffices for all sixteen columns.

The source/projective prefix times Q_00^16 is one positive component of
the refined comparator. It is neither a probability lower bound nor an
upper bound on the whole failure event. If this component exceeds the
target, that particular witness and comparator cannot prove the target;
the construction itself is not thereby rejected.

The isolated helper `hd_two_sum_pre_toeplitz_zero_component.py` also
extracts the largest ordinary state-word component. It maximizes paths
within each ordered boundary-state pair, then uses a max-plus product
across columns. It does not independently select incompatible boundary
states. Both diagnostics use the existing saved W tables and the three
frozen activity marks; they build no new counts or regional operators.

Six tiny tests pass, including direct labelled group-mask enumeration,
zero/one cut probabilities, matching invariance, and boundary-state
consistency. The parent independently reran them, then completed the
component diagnostic in 1.346 seconds.

| Activity schedule | Old refined all00 component | New refined all00 component | New dominant allowed size |
|---|---:|---:|---:|
| Zero log | 1319.573022 | 1319.573022 | 0 |
| Frozen scale1 | 1575.248372 | 1440.257856 | 6 |
| Frozen scale2 | 2707.953632 | 2572.940946 | 6 |

Every entry includes the same source/projective prefix as the schedule
gate. The zero-log result is dominated by the synthetic all-empty
comparison term. Raising the activity marker suppresses that term but
exposes larger allowed sets. The new outer law reduces those nonempty
components, yet not enough for any frozen witness to close.

The largest ordinary state-word component is all00 in every column for
both outer laws and all three marks. For the new law its values are
5331.015111, 6463.705251, and 7596.421932. Thus the ordinary loss is not
mainly a large family of alternating state paths. Restoring the empty00
branch is important, but its all00 component already remains above the
target. Enumerating millions of augmented paths cannot reverse this
positive-component obstruction at the fixed witnesses.

The component artifact is
`hd_two_sum_pre_toeplitz_zero_component_gate_float.json`, SHA-256
`5431a784c3f8f1b48bb29771883c9b68dbef8f0cc41bd61f606bbe69cbcf8dc7`.

The local gains therefore do not yet overcome the loss from expanding
the positive support sum. This conclusion rejects neither the code nor
every possible proof of this outer modification. The next proof step
should test later fractionalization at column level, using the correct
source-conditioned regional law, or a stronger construction. It should
not continue an ordinary-envelope marker grid or launch the much larger
augmented-path calculation with these unchanged witnesses. No native
experiment is justified by these bounds alone.

- Component helper SHA-256: `0a2fb671c8a3869adac56254d13ae16016fd745131b74cdc40ae2efedfb842d1`.
- Component tests SHA-256: `97ecf9e744f7e31dc95cda3908994b72e5b1fff7e61a583d1272067ee35c360d`.
