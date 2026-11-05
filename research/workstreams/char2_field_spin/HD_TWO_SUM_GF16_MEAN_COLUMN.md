# Full outer masks in the zero-boundary first moment

This concerns the original G16 two-sum outer with diagonal input/output
scalings, not the experimental pre-H Toeplitz modification.
The r18, homogeneous-arity8, F1 comparison now includes every outer
moment branch along the all-zero regional boundary path. Its FLOAT
log2 bound is -88.130205, including the exact source count. This is one
positive component of a failure majorant, not a complete class or
whole-code certificate. Live boundary paths remain unresolved.

The witness is z=691/700 and alpha=1. F1 means that the larger full-pair
count of the fixed source set in the two route matchings is one.
Every regional subset therefore has at most one full pair. The source
count is

```text
154355220469939768130077702840442839490778956275433472
```

It is charged once, outside the per-source first moment. No source
geometry or constituent row set is resampled across columns.

## Positive outer decomposition and a usable input floor

Fix an active constituent set A of size eight, and let O be its
complement among sixteen regional coordinates. For B subset A, write
b=|B|. The exact closed kernel supplies three nonnegative counts:

\[
 q_b=Z_2(b),\qquad h_b=Z_1(b)-Z_2(b),\qquad
 f_b=255^b-q_b-8h_b.
\]

Here Z_1 and Z_2 count nonzero assignments satisfying one and two
moment equations, respectively. The corresponding output supports are
B, B union (O minus {j}) for each j in O, and B union O. Their counts
are q_b, h_b, and f_b. All outside-hole positions are separate labelled
branches.

The quiet branch has pre-H support weight b. A nonquiet input has
pre-H weight at most eight. Since h_0=h_1=f_0=0, every nonquiet output
has weight at least nine. Consequently

\[
 \kappa(w)=\min(w,8)
\]

equals the pre-H weight for every quiet branch and upper-bounds it for
every nonquiet branch. This upper surrogate is independent of A.
The actual sixteen-column RS word has total pre-H weight at least
9*8*18=1296. Its surrogate weight is therefore also at least 1296.
Retaining only the surrogate tail is a relaxation, not an assertion that
the two weights agree on nonquiet branches.

For formal region markers t_i and an input-floor marker X, the positive
fixed-A kernel is

\[
 \sum_{B\subseteq A}t^B\left[
 q_bX^b+X^8\left(h_b\sum_{j\in O}t^{O\setminus\{j\}}
                         +f_b t^O\right)\right].                 \tag{1}
\]

Thus no new joint input/output kernel is required for this first control.
The individual row9 floors and individual group72 floors are not
enforced by the aggregate surrogate.

## Fixed row sets and the raw Holder weight

Let N_A(T) be the exact uniform-input count for one labelled output mask
T. For any fixed weight-w target, define

\[
 W(w)=\left(\sum_{|A|=8}N_A(T)^{16}\right)^{1/16}.                 \tag{2}
\]

The established closed kernel makes this sum depend only on w. For a
fixed sixteen-column target word, Holder bounds the sum over its single
fixed A by the product of the sixteen weights in (2). Because kappa
depends only on T, its marker remains the same for every A in that sum.
The bound therefore remains coefficientwise after attaching
X raised to the sum of the sixteen surrogate weights.

Apply this argument separately to each of the eighteen groups. W
already charges their row-set choices once. There is no extra factor
binomial(16,8)^18. The implementation uses raw W, with the whole-word
RS factor 255^-1152 outside. Equivalently, one could divide each W by
255^4, but must then omit the external RS factor.

For quiet-only control, W vanishes above weight eight and equals full W
at every weight at most eight. The full comparison adds the positive
nonquiet weights nine through sixteen.

## Scalar moments and chronological state evolution

Let R(m) be the raw alpha1 regional matrix after taking the F1 maximum
over feasible pair histograms. Suppose a probability law L and a fixed
positive matrix C satisfy

\[
 R(m)\le M_m C,\qquad M_m=\mathbb E L^m,
\]

entrywise for every m from zero through eighteen. The frozen law has
Pr[L=B]=p, Pr[L=0]=1-p, with

```text
log2(B) = -8.907630419774677
log2(p) = -24.529436473994224
```

Conditional on exactly d of the sixteen independent auxiliary atoms
being nonzero, every comparison mask lies within those d positions.
Summing all labelled masks and all eighteen groups gives the column
polynomial

\[
 Q(X)=\sum_{d=0}^{16}{16\choose d}p^d(1-p)^{16-d}
 \left[\sum_{w=0}^d{d\choose w}W(w)B^wX^{\kappa(w)}\right]^{18}.   \tag{3}
\]

The binomial factors count atom subsets and target masks, respectively.
Neither is a row-set factor. The scalar moment comparison removes the
occupancy dependence from the boundary-state matrices, leaving their
ordinary chronological product C^256. No matrix commutation or
probability-row normalization is used.

With N_F1 denoting the source count above, the resulting complete
comparison, when all regional entries satisfy the envelope, is

\[
 \frac{N_{F1}z^{-734003}}{255\,255^{1152}}
 (e_0C^{256}{\bf1})
 \sum_{k\ge1296}[X^k]Q(X)^{16}.                                  \tag{4}
\]

The common-C envelope is sound but numerically poor. Its measured
state-factor log is +14699.747295. This turns the earlier quiet/all00
-432.030562 control into +14267.716734. The full-outer common-C gate
was therefore parked without another production run.

## Completed zero-boundary control

To isolate the outer remainder, the saved-bank driver retains only the
all-zero boundary path. It rechecks R_00(m)<=M_m for all nineteen
occupancies, then uses the matrix with C_00=1 and all other entries zero.
This deliberately omits live boundary paths; (4) now describes a
selected positive component, not a whole-event upper bound.

The parent ran exactly two controls in 0.182355 seconds:

| Outer branches | FLOAT log2 component |
|---|---:|
| Quiet only, raw Holder weights | -88.130205160649 |
| All quiet/hole/full branches | -88.130205160649 |

Equality at displayed FLOAT precision does not mean the nonquiet mass
is zero. Its added contribution is negligible at this witness. Fixed-A
Holder is roughly 344 bits looser here than the earlier common-incidence
quiet bound, but the component still clears -40. This result isolates
live boundary chronology, rather than nonquiet outer mass, as the next
unresolved part of this particular comparison.

Artifact: `theory/hd_two_sum_gf16_mean_zero_gate_float.json`.
SHA-256: `610eb126c6dd2c34bf8ad8a36ec447d372d5fed3321d05b216d2b54003e57cd2`.

The helpers are `hd_two_sum_gf16_mean_column.py` and
`hd_two_sum_gf16_mean_zero.py`. Eleven tiny/mocked tests pass. They check
labelled latent enumeration, full matrix products, branch separation,
fixed-A Holder, normalization, artifact scope, and timeout preservation.
The driver builds no regional bank and fits no parameters.

## Next boundary-path contraction

A positive matrix-valued moment law can avoid the single-C loss:

\[
 R(m)\le\sum_t B_t^m C_t.
\]

At m=0, every atom, including B_t=0, contributes its matrix weight.

For fixed atom counts, the scalar mask sum depends only on their
multiset. The coefficient of the corresponding monomial in
(sum_t Z_t C_t)^16 sums the ordered matrix words exactly. This permits
separate atoms for 00, cross-state, and 11 transitions without assuming
that their matrices commute. The next controlled gate is this complete
homogeneous-class comparison, not a new outer construction or a broad
witness search.

Cross-state transitions additionally satisfy R_01(0)=R_10(0)=0. The
isolated `hd_two_sum_gf16_cross_coverage.py` preserves that restriction
using a positive coverage DP. Fix one state/atom pattern with c cross
coordinates, k active00 coordinates, and l live11 coordinates. At a
scalar floor marker v, let

\[
 A_j=B_x^j\sum_{u=0}^k\sum_{s=0}^l
 {k\choose u}{l\choose s}W(u+j+s)B_0^u
 v^{\kappa(u+j+s)}.                                               \tag{5}
\]

A_j is the weight of one fixed j-subset of the cross coordinates; it
does not include binomial(c,j). If d cross positions have been covered,
the next group covers n new positions with transition weight

\[
 \sum_{h=0}^d {d\choose h}{c-d\choose n}A_{h+n}.                   \tag{6}
\]

Start at d=0, apply eighteen group transitions, and require d=c.
Formula (6) counts every labelled group-mask word once. It never
subtracts upper bounds. For c=0 it reduces to the original eighteenth
power. A caller using scalar v>=1 must also include the whole-word
floor charge v^-1296 once.

Six exact small tests pass, including labelled coverage enumeration and
the independent count (2^r-1)^c for unit subset weights. The coverage
helper has no production driver or numerical proof result yet. It
retains fixed row sets through the preceding Holder comparison; it does
not sample a fresh constituent set per group or column.

The separate ordered-atom control subsequently gave +3497.359683 with
its exact positive degree tail. Its scalar control used a marker derived
from the earlier quiet word and was weaker. The new isolated driver
`hd_two_sum_gf16_cross_coverage_control.py` instead derives one marker
from adjacent coefficients 1296 and 1297 of the full ordered word,
summing the two terminal states first. It then compares coverage and
unrestricted masks under identical path weights. The marker is fixed
across every path and column; there is no marker search.

Five small driver tests pass. They compare labelled state/mask
enumeration, reproduce the unrestricted atom operator, and check scope
and timeout handling. Independent read-only review cleared the driver.

The parent's production control exposed an unstable marker selection.
Adjacent full-word coefficients selected log2(X)=65.193451, producing
+68047.778631 for both covered and unrestricted comparisons. Odd/even
degree effects make that local ratio unsuitable for this polynomial.
The result is preserved in
`hd_two_sum_gf16_cross_coverage_gate_float.json`, SHA-256
`4286c359d475aa8d4253fac094d74eccaca65a904b5d2678ca16a07d08ea39dd`.
It gives no useful evidence about coverage gains near a good marker.

The separate wrapper `hd_two_sum_gf16_cross_coverage_convex.py` leaves
that failed driver unchanged. For the full positive word polynomial
F(X)=sum_d a_d X^d, it minimizes

\[
 \log_2 F(2^\ell)-1296\ell,\qquad \ell\ge0.
\]

The derivative is the tilted mean degree minus 1296; the second
derivative is ln(2) times the tilted degree variance. The wrapper uses
the endpoint zero when appropriate. Otherwise it brackets the stationary
point and performs 48 bisections, with an explicit finite bracket cap.
Both coverage controls share the resulting single global marker.

Five new tiny tests pass. The test polynomial
16+X^4+2^-100 X^5+X^8 demonstrates the issue: its adjacent ratio suggests
log2(X)=100, whereas its convex optimum is approximately 1/2. The wrapper
also checks the scalar matrix result against direct evaluation of the
full word polynomial. Independent read-only review cleared this control.

The parent completed its single bounded production run. The selected
log2(X) was 1.109035132306877, with tilted mean degree 1296. The matched
results were:

| Global convex marker | FLOAT log2 comparison |
|---|---:|
| Unrestricted cross coordinates | +3543.6142451091223 |
| Every cross coordinate covered | +3543.6142451091237 |

The difference is FLOAT rounding; it is not a negative coverage gain.
The exact coverage sum is a positive subsum of the unrestricted one.
The artifact is `hd_two_sum_gf16_cross_coverage_convex_gate_float.json`,
SHA-256 `0b1039c3a01ec744c41971d3ca1fbbb513d5054f5d77f6057559402908d179ae`.

The corrected marker confirms that empty cross transitions are not the
dominant obstruction at this witness. The parent's path diagnosis found
four excursions through live state, full occupancy m=18 at the cross
regions, no11 transitions, and no active00 atoms. Enforcing nonempty
cross occupancy cannot remove that path. The adjacent-ratio marker
control is discontinued; no further marker fitting is justified by this
result.

The next cycle should address the common-source route event supporting
these costly cross-state transitions, or their fractional-moment
chronology. It should preserve the passing full-outer/all00 component
and first target the complete homogeneous a8/F1 class. The current
positive full-state comparisons still provide no certificate for that
class, much less the whole code. All27 focused helper/driver tests pass;
no native work or benchmark occurred in this proof iteration.

Run the small tests from `theory`:

```text
python -B -m unittest hd_two_sum_gf16_mean_column_tests hd_two_sum_gf16_mean_zero_tests hd_two_sum_gf16_cross_coverage_tests
```
