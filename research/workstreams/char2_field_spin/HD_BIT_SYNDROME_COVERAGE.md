# Coverage plan for the bit-syndrome outer

The current selected-class replay does not cover the whole code. This note
gives a finite coverage partition and a positive composition rule for mixed
constituent arities. It records proof interfaces and scalar headroom read
from saved receipts, not newly evaluated regional or word bounds.
No regional provider, contraction, or native experiment was run for this note.

Subsequent numerical progress is recorded in the
[kernel note](HD_BIT_SYNDROME_KERNEL.md#a-first-mixed-arity-extension):
the homogeneous a8/r18/F<=1 class now has directed 256/384-bit bounds and
independent 512-bit final checks. All mixtures of a8 and a9 at r18,F<=1
also clear a fresh common-measure FLOAT bound, near -202.59 bits.
The subsequent activity-marker and exact-empty-step repairs extend this
to all mixtures of arities seven, eight, and nine, near -76.39 bits.
These mixed extensions are not directed. The remaining scheme below is
still a coverage plan, not a completed ledger.

The construction is the [bit-syndrome outer](HD_BIT_SYNDROME_OUTER.md),
with sixteen RS[16,8] rows per group and the **original one-state h4 inner**.
There are 8192 groups and 256 independently routed regions. Routing uses
alternating adjacent matchings, uniform pair permutations, and independent
fair pair swaps. The ideal outer diagonals have independent uniform nonzero
entries sampled once at setup. The native gauge requires its established
exact lift to this ideal ensemble.

The event to bound is existence of a nonzero scalar message with output
weight at most L=734003, over the randomized encoder setup. The target is
a total failure bound below 2^-40. The separately certified Toeplitz encoder
and the exploratory post-refresh inner are different constructions.

## What the selected class establishes

Fix the active group set G and write r=|G|. In each group g, let A_g be its
fixed nonempty subset of the sixteen constituent rows. Define

\[
 a_g=|A_g|,\qquad A_{\rm tot}=\sum_{g\in G}a_g.
\]

The row sets remain fixed across all sixteen RS columns. Let e(G) and o(G)
count the full source pairs in the two adjacent matchings, and put
F(G)=max(e(G),o(G)). The exact integer routine
`theory/hd_xor_group_geometry.py:matching_max_counts` counts indexed G by F.
Write N_{r,F} for that count. Zero-count strata can be omitted.
The finite range is 0<=F<=min(4096,floor(r/2)).

The selected directed result has r=18, every a_g=8, and F<=1. It includes every
labelled eight-row set in each group, not just a common A shared by groups.
The [kernel note](HD_BIT_SYNDROME_KERNEL.md) records both directed replays
and independent final checks. Their narrow scope does not establish a
whole-code certificate.

The same Fmax1 regional envelope covers F=0 and F=1. Extending the source
union costs only the exact change from N_{18,1} to N_{18,0}+N_{18,1}.
It does not cover F>=2, other r, or heterogeneous a_g.

## A mixed-arity bound without enumerating cohort histograms

Fix r, a source stratum, and one output mark 0<z<1. Use one common positive
matrix-valued atom measure that dominates every required regional matrix.
This measure must be independent of the chosen tuple (A_g). Its ordered
path functional, denoted by \(\mathcal L\), includes the initial zero state
and terminal all-ones vector. Expanding its nonnegative matrix products
makes \(\mathcal L\) an ordinary positive scalar linear functional on paths.
It need not be a probability measure, and the matrices need not commute.

For each arity a, fix an input-support marker X_a>=1. For each labelled
set A of size a, let K_A be a valid positive column majorant for its actual
H* mask law. The majorant must upper-mark true input activity. For a path
with atom labels t_{c,i}, define

\[
 g_A:=\prod_{c=1}^{16}K_A(t_{c,0},\ldots,t_{c,15};X_a),
 \qquad
 J_A:=\mathcal L[g_A^r].
\]

A positive computed upper bound on J_A may replace J_A below. In the
column-matrix implementation, J_A is bounded by
\(e_0^T Q_A(X_a)^{16}\mathbf1\), where Q_A integrates K_A^r against the
same atom measure. Any Young decomposition is applied before this
sixteen-column product, with fixed weights for that labelled A.

The per-group RS counting bound is 255^(-8a). Each of its a nonzero rows
has at least nine nonzero columns, so its input-floor factor is X_a^(-9a).
Define

\[
 b_a:=\sum_{\substack{A\subseteq\{0,\ldots,15\}\\|A|=a}}
       255^{-8a}X_a^{-9a}J_A^{1/r},
 \qquad B_r(u):=\sum_{a=1}^{16}b_a u^a.
\]

One may additionally mark nonzero columns by Y_a>=1. Every active group
has at least nine such columns, so b_a then includes Y_a^-9 and J_A
uses the marked column kernel. The [kernel note](HD_BIT_SYNDROME_KERNEL.md#enforce-the-rs-column-activity-before-combining-groups)
proves a correction that reuses the same positive atom measure.

These exponents are **per group**, not divided by r. For homogeneous
arity eight at r=18, the final eighteenth power recovers 255^-1152 and
X_8^-1296. There is no additional row-set or cohort multiplicity.

For a fixed ordered tuple of row sets, Holder's inequality gives

\[
 \mathcal L\!\left[\prod_{g=1}^{r}
      255^{-8a_g}X_{a_g}^{-9a_g}g_{A_g}\right]
 \le\prod_{g=1}^{r}
      255^{-8a_g}X_{a_g}^{-9a_g}J_{A_g}^{1/r}.
\]

Sum this inequality only over tuples with A_tot<=B. Polynomial expansion
then supplies their row-set choices and group placements exactly once.
The resulting source-stratum bound is

\[
 U_{r,F,\le B}\le
 \frac{N_{r,F}z^{-L}}{255}
 \sum_{k=r}^{\min(B,16r)}[u^k]B_r(u)^r.
\]

There is **no outer root** around this coefficient tail. The source count,
output mark, and single projective division stay outside it. With no arity
restriction, the tail becomes B_r(1)^r. Thus separately derived homogeneous
bounds under the same functional combine through a sum of r-th roots,
not through their direct sum.

The arity polynomial is much smaller than a list of all cohort histograms:
its maximum degree is 16r, or 512 at r=32. This does not make the regional
or fixed-A contractions free. Those remain the expensive proof inputs.
The exact tail can be truncated after each positive polynomial product.

Different source strata may use different witnesses. Within one use of
the formula, z and the positive functional must remain common to every A.
Type-specific X_a and fixed-A Young weights are allowed. Combining b_a
from unrelated output marks or unrelated atom measures is not justified.

## The portable high-arity comparison

There is a second bound that does not need the H* mask kernel. Fix G,
a target outer array y, and all outer output diagonals. Invert those
diagonals and H* column by column. The resulting pre-input-diagonal array
determines one constituent-row-set assignment. It does not create a union
over all assignments.

For each nonzero RS[16,8] row of support w>=9, shortened-code dimension
w-8 bounds its full-support count by 255^(w-8). Independent input diagonals
spread its nonzero amplitudes over 255^w assignments. Its expected
membership is therefore at most 255^-8. Impossible supports contribute zero.
Consequently the expected counting measure of A_tot>=A_0 is pointwise at
most 255^(-8A_0). Averaging the conditioned output diagonals preserves it.

Relative to uniform field amplitudes on the 256r potential output
coordinates, the comparison density is

\[
 \beta_{r,A_0}=256^{256r}/255^{8A_0}.
\]

There is no row-set, arity-deficit, or cohort factor. The proof uses
invertibility of H*, the ideal diagonal law, and RS[16,8], not its old
rank-two mask law. It uses a **uniform-amplitude** route provider, including
zero amplitudes, rather than the current forced-nonzero provider.

The [earlier high-arity adapter](HD_TWO_SUM_GF16_HIGH_ARITY_CERTIFICATE.md)
already supplies exact arithmetic for the following three anchors:

| r | High class | z | Fractional exponent |
|---:|---:|---:|---:|
| 18 | A_tot>=247 | 691/700 | 1/5 |
| 24 | A_tot>=317 | 172/175 | 1/5 |
| 32 | A_tot>=463 | 171/175 | 1/5 |

Its directed aggregate is -62.24877612169 bits, for the older rank-two
construction. These are practical candidate cutoffs for H*, because the
comparison proof and scalar factor are unchanged. A new H* theorem binding
and receipt are still required. The old receipts must not be relabelled.
The anchors establish neither a contiguous r band nor a useful cutoff
for every r. In particular, the sharper Toeplitz density does not transfer.

### What the saved bulk receipts suggest

The closed Toeplitz ledger provides additional inexpensive proposals.
Suppose one saved all-source weight-r bound is h_r in log2 units, its
fractional exponent is alpha, and its per-group log-density is d_old.
Changing only the outer comparison gives the formal scalar adjustment

\[
 h_r+\alpha\bigl(2048r-8A_0\log_2 255-rd_{\rm old}\bigr).
\]

For the saved `inherited-diagonal` mode, d_old=2048-128log2(255).
The adjustment is then 8alpha(16r-A_0)log2(255). For the saved
`invertible-toeplitz` mode, d_old=1280-32log2(255); one must additionally
pay alpha*r*(768-96log2(255)). Omitting this full-arity correction would
incorrectly transfer the sharper Toeplitz bound to H*.

Seven saved per-weight records give the following **FLOAT scalar proposals**
at target -60. Each saved alpha is read from its exact witness and is close
to 0.1. No regional coefficient or word contraction was rebuilt here.

| r | Saved log2 ceiling, rounded | Saved density mode | Proposed A_0 | Allowed deficit 16r-A_0 |
|---:|---:|---|---:|---:|
| 1024 | -165.87 | inherited-diagonal | 16368 | 16 |
| 1536 | -177.65 | inherited-diagonal | 24558 | 18 |
| 2048 | -1063.71 | inherited-diagonal | 32612 | 156 |
| 2496 | -2545.02 | invertible-toeplitz | 39569 | 367 |
| 3072 | -6083.57 | invertible-toeplitz | 48237 | 915 |
| 3584 | -12267.76 | invertible-toeplitz | 55466 | 1878 |
| 3936 | -18298.29 | invertible-toeplitz | 60158 | 2818 |

The source records are the `results` entries at those exact weights in
these files under `theory`, each with suffix `_directed256.json`:

- `hd_shifted_potential_band_1024`;
- `hd_shifted_potential_band_1536_1500_1571`;
- `hd_shifted_potential_band_2048_1882_2153`;
- `hd_shifted_potential_band_2496_2275_2683`;
- `hd_shifted_potential_band_3072_2828_3285_toeplitz`;
- `hd_shifted_potential_band_3584_3326_3768`;
- `hd_shifted_potential_band_3936_3769_4050`.

This is evidence for a cheap high-side adapter across bulk weights, not
for a cheap complete bulk proof. These sample thresholds remain near
full arity: the permitted deficits are about 0.1% to 4.5% of 16r. They
do not cover all intermediate or sparse arities. Nor do seven sample
records certify the surrounding bands.

A bounded next audit can extract one proposed A_0(r) from every saved
per-weight record without fitting new witnesses. The resulting complete
schedule would identify the low-side domain precisely. Then regenerate
the selected bounds with H* metadata, its established density, and fresh
directed scalar and operator arithmetic. Do not import the old receipts
into an H* ledger, or treat these rounded FLOAT cutoff estimates as
certificates. The existing providers and exact source counts remain
reusable; their old construction identifiers and densities do not.

## A finite disjoint target partition

The route-independent inverse-image argument also gives

\[
 U_{\rm terminal}(T)\le
 \frac{\operatorname{Ball}_{256}(2^{21},734003)}{255^{8T+1}}.
\]

Here one conditions on the inner setup, routing, and output diagonals,
then inverts a candidate low-weight output. No source count is added.
The existing scalar arithmetic at T=122382 has ceiling -165.775266208552.
Its applicability to H* needs the same new ensemble binding as above.

This threshold covers only the nearly full constituent tail. It is
reachable only for r>=7649; it does not cover every message at those r.
For example, r=8192 with one row per group has A_tot=8192. There is no
supported claim of a much earlier route-independent terminal cutoff.

Choose an integer A_0(r) in [r,16r+1] for each r. Choosing 16r+1 disables
the high-side class when no useful high-side witness is available. The
following target sets are disjoint and exhaust all nonzero messages:

1. A_tot>=T, handled by the global terminal comparison.
2. For each r, A_0(r)<=A_tot<min(T,16r+1), handled by the uniform comparison.
3. For each r and nonempty source stratum F,
   r<=A_tot<=min(A_0(r)-1,T-1,16r), handled by the low-arity polynomial.

An upper bound for a larger high-side event can bound the second target
set; no subtraction of rounded bounds is needed. Numerical majorants may
overlap even though the target partition is disjoint. Sum their directed
ceilings once according to the target ledger.

Exact F strata may be joined into bins. Count a bin by summing its exact
N_{r,F}, and use a regional envelope valid for every F in the bin. For
a regional subset U of G, its full-pair count is at most
min(F,floor(|U|/2)). An F bin therefore admits the upper bound from its
largest F. It does not justify a lower bound derived from an exact
matching count. Independence between regions does not resample G.

This partition is sufficient as an accounting scheme. No current set of
receipts supplies all its numerical cells. A whole-code claim requires
complete source/arity coverage and a directed total below 2^-40.
The illustrative per-weight target -60 is not an allocation for every
possible (r,F) cell; the final union must include their actual number
and individual ceilings.

## Construction-specific inputs and next gates

The new exact interface is
[`hd_bit_syndrome_cross_fibers.py`](theory/hd_bit_syndrome_cross_fibers.py).
For every labelled A it uses the actual cross map H*_(A^c,A), its image,
and its kernel. The inside differences are the image of that kernel under
H*_(A,A). They are not automatically the quiet kernel.

For fixed outside output z and allowed inside support S, it solves

\[
 H_*^{-1}[A^c,S]y_S=H_*^{-1}[A^c,A^c]z.
\]

The supported count is zero if inconsistent, and otherwise
256^(|S|-rank(H*^-1[A^c,S])). Exact subset inversion gives shifted support
counts. The helper provides local algebra, not complete all-arity kernels
or a distance certificate. Its extraction budget is explicit.

The next gates should address three distinct missing dimensions:

1. **Mixed arities at the current source weight.** Build a7 and a9 kernels
   from actual cross fibers, then combine a7/a8/a9 under one z and atom
   measure. Each neighboring arity has 11440 labelled sets. Its 240
   exceptional sets need separate fiber treatment. At a9, zero outside
   can be nonquiet. This test exercises mixed-group normalization before
   extending to all sixteen b_a and the low class A_tot<=246 at r18.
2. **The residual source geometry.** Test homogeneous a8 at r18,F=2 using
   fresh f=2 regional coefficients. The F>=2 source fraction is small,
   but its count alone does not control larger regional returns. If a
   worst-F tail is loose, split F=2 from the remaining tail before making
   the arity partition finer.
3. **Weight extension with a complete arity sum.** Once the r18 low side
   is credible, test the complete low side at r24 or r32 against its
   independently justified high-side cutoff. This is more informative
   than collecting additional homogeneous-a8 receipts.

The selected m<=18,f<=1 regional bank is mathematically reusable for
r<=18 at the same mark and same inner construction. The group exponent,
source counts, RS factors, and input floors must still be recomputed.
A result at r18 does not certify neighboring r by monotonicity or
interpolation. Extending r above18 or F above1 needs additional regional
coefficients and renewed all-entry moment domination checks.

The original singleton full-output law remains valid because every H*
entry is nonzero. The older pair law does not: H* pairs can have several
holes. Older rank-two dense, pair, and mixed-cohort receipts therefore
do not transfer. Regional geometry algorithms and uniform high-arity
comparison methods remain reusable with explicit construction bindings.
The arity-eight affine count orbits likewise do not authorize collapsing
labelled A in chronological contractions or extending their kernels to
other arities.

For large r, the established positive/alias coefficient machinery is a
possible computational method, not a ready H* certificate. The low side
still needs forced-nonzero coefficients and valid normalizer envelopes.
No present evidence supplies a universal low-r cutoff, a complete bulk
band, or an all-arity terminal range beyond the threshold stated above.

## Provenance and scope of this note

The per-group normalization and the low-arity coefficient formula were
independently reviewed. This note adds no production arithmetic or tests.
The cross-fiber implementation inspected here has SHA256
`d8658b8c3a67d2925e776b63bd62b637456ea9e22c8c75d64dda1ba958e9104e`;
its tests have SHA256
`38ffe1b53d884cc9ec4ebd2931c957519ea8d8400b382d038263686220c797e1`.
These pins identify the local interface, not a whole-code certificate.

The mixed a7/a8/a9 test is now successful after two structural proof
repairs. Next, check whether the exact-empty-step improvement supports
a simpler complete-arity bound before adding more separate arity kernels.
Retain the explicit column activity constraint in that test. The new
[rank-eight mixer](HD_BIT_SYNDROME_STAR_OUTER.md) is a separately checked
fallback, not a replacement with transferred certificates. Keep the
completed directed homogeneous-a8 result separate from wider FLOAT scope.

The subsequent broad-coverage screen rejects scalar-state collapse at
the current mark: quiet affine planes give a +1947-bit lower bound on
that proof majorant. A conditional uniform-cover calculation helps the
middle arities but fails both tails; its cover hypothesis is not yet
universal. Neither screen supplies a completed coverage range.

The next proof operation is instead the
[global quiet-setup decomposition](HD_BIT_SYNDROME_QUIET_EXPURGATION.md).
Its exact probability charge, at most 2^-50.771255305346, covers all8192
groups and all row sets of arity at most seven for the unchanged H* code.
It does not have to be paid again for each source stratum. The remaining
unconditional word count must enforce at least one nonquiet column in
each such group. The [inverse column schedule](HD_BIT_SYNDROME_COLUMN_PERMUTATIONS.md)
is a separate candidate with stronger deterministic floors at arities
four and five. Its selected a7/a8/a9 FLOAT comparison remains negative,
but no whole-code certificate or native timing transfers to that candidate.

The [subsequent sparse audit](HD_SPARSE_ACTIVITY.md) now includes every
labelled row set of arities2..7 at r18/F<=1. Four homogeneous controls
pass numerically; arities five and six still fail. This does not fill a
mixed-arity or all-source coverage cell. Exact output-fiber repairs are
the next test; further fitting of the tested input/nonquiet markers is
parked. Independent local output shuffles define a separate ensemble and
must not inherit this coverage ledger by relabelling.
