# G16 whole-code coverage: a factorized cohort gate

The G16 affine candidate still lacks a whole-code distance certificate.
The current sparse and one-/two-dense-group gates do not cover arbitrary
cohorts or active-group weights. The controls below cover every cohort and
every source matching class in one small weight band. Their scalar
common-potential majorants fail the recorded bounded tests. The remaining
proposals retain matrix states or compare a near-full constituent stratum
with the already studied uniform-output measure.

The construction here has sixteen RS[16,8] constituents in each of 8,192
groups, the checked affine mixer with labels `06 + GF16`, independently
swapped paired routing, and the unchanged width64/state16 inner. Setup
is sampled once. The objective remains a setup-failure probability below
`2^-40` for a nonzero scalar message of output weight at most 734003.
The retained Toeplitz certificate is a separate construction.

## The coverage gap is two-dimensional

Fix the active-group set `G`, let `r=|G|`, and let `A_g` be the fixed nonempty
constituent-row set in group `g`. Write `A_tot=sum_g |A_g|`. A complete
partition ranges over `1<=r<=8192` and `r<=A_tot<=16r`, including all row-set
assignments and both ordered source matching counts. Strong gates at `r=17`
with at most two dense groups are only a small portion of this domain.

The route-independent terminal argument extends to this invertible affine
outer. Condition on routing, inner setup, and the outer output diagonals.
For a candidate low-weight output, invert those maps and the fixed affine
mixer. Each nonzero prescribed RS[16,8] row has membership probability at
most `255^-8` over its still-unexposed independent input diagonal. This is
the shortened-MDS exact-support count `255^(w-8)` divided by `255^w`.
Consequently

\[
 \Pr[\text{a bad word with }A_{\rm tot}\ge A_0]
 \le \frac{\operatorname{Ball}_{256}(2^{21},734003)}{255^{8A_0+1}}.
\]

There is no extra source/cohort count in this inverse-image argument.
The native scalar gauge must be justified through the equivalent independent
ungauged diagonal ensemble, not by conditioning directly on normalized keys.
Taking `A_0=122382` gives exactly the exponent of the previously checked
XOR2 terminal threshold `61191`, since `8*122382=16*61191`. Its numerical
ceiling is `-165.775266208552`; applying it to G16 requires the stated
ensemble transfer, not relabelling a receipt.

This tail does not cover all large `r`. For example, `r=8192` with one row
per group has only 8,192 active constituents. The main remaining domain is
all intermediate constituent counts below that threshold, across essentially
the entire active-group range. A ledger based only on `r` needs a bound
that has already summed every cohort.

One narrower transfer remains available. If every active group has all
sixteen rows nonzero, its expected output counting measure is dominated by
`256^256/255^128` times the uniform 256-symbol measure. This is the legacy
outer interface of the paired Toeplitz proof. Only receipts using that
interface, or freshly corrected/replayed bounds, can transfer to this pure
full-arity class. The sharper Toeplitz atom bound does not transfer. This
reuse does not handle mixed arities.

## Near-full strata: what the uniform comparison actually transfers

The pure-full-arity comparison extends to a precisely charged mixed-arity
stratum. Fix an active-group set and, within one group, a nonempty fixed
constituent set `A` of size `a`. For prescribed nonzero row supports, the
shortened RS count and the independent nonzero input diagonal give an
expected pre-mixer counting measure at most `255^(-8a)` at every amplitude
word. This is counting-measure domination, not a likelihood statement for
one fixed codeword. Enlarge its support to every word in `F_256^(16a)`.

For a column output mask `T`, let `N_A(T)` count inputs from this coordinate
subspace whose image under the fixed invertible affine mixer has exactly
that mask. Injectivity gives `N_A(T)<=255^|T|`. The output diagonal assigns
mass `255^(-|T|)` to each prescribed amplitude word with that mask. Thus
the column pushforward has density at most one relative to counting
measure, and the complete expected group output count is at most
`255^(-8a)` at every 256-symbol word. Relative to the uniform group law,
the valid density is

\[
 \beta_a=\frac{256^{256}}{255^{8a}}.
\]

No output support or RS floor is retained in this comparison. That loss is
intentional. The diagonal law is the independent ungauged ensemble already
used in the fixed-row-set proofs; native key normalization must preserve
that scalar map distribution before applying the comparison.

### Exact normalization of the old fractional bound

The current `hd_shifted_potential_band.py` evaluates

\[
 255^{-\alpha}z^{-\alpha L}\beta^{\alpha r}\Phi_r(z,\alpha,v).
\]

The support sum `Phi_r` includes the endpoint factor `v_0/min(v)` once,
outside the fractional exponent. The packet and route operators depend on
the uniform-output comparison, not on the source RS block sizes. G16 has
the same 256 coordinates per group, 256 chronological regions, paired
route law, and inner map. Therefore `a=16` reproduces the **legacy** density
exactly:

\[
 \beta_{16}=256^{256}/255^{128}.
\]

It does not reproduce the sharper Toeplitz density
`beta_T=256^160/255^32`. Converting a sharper weight-`r` bound to this G16
comparison adds

\[
 \alpha r\log_2(\beta_{16}/\beta_T)
 =96\alpha r\log_2(256/255)
\]

bits. The unscaled difference is `0.5420700615496361` bits per active group.
The projective normalization is unchanged and must not be moved outside
`alpha` or charged a second time.

For total constituent deficit `D=16r-A_tot`, the exact number of fixed
row-set assignments is

\[
 N_r(D)=[t^D]\bigl((1+t)^{16}-t^{16}\bigr)^r
 \le \binom{16r}{D}.
\]

The subtraction excludes an empty active group. Equality with the binomial
coefficient holds for `D<16`; the upper bound remains safe for larger
deficits. Sum the counting measures for all assignments with deficit at
most `D_0` **before** the fractional-moment step. Their density is bounded by

\[
 \beta_{16}^r S_r(D_0),\qquad
 S_r(D_0)=\sum_{d=0}^{D_0}\binom{16r}{d}255^{8d}.
\]

Consequently the additional log bound is `alpha*log2(S_r(D_0))`. These
binomial row choices occur once inside the density; they are not an
additional outer union factor. A single missing row costs
`8 log2(255)=63.95482749487086` density bits before its row-choice count and
the fractional exponent. Alternatively, a fixed per-group cutoff `a>=a_0`
has density `256^256 sum_(a=a_0)^16 C(16,a)255^(-8a)` per active group.
Both statements allow arbitrary mixtures of those fixed row sets; neither
resamples a row set between columns. This summed-density argument is valid
but unnecessarily loose. The following inverse-image argument supersedes
its row-choice penalty.

### Stronger total-arity comparison: no row-choice penalty

Fix the active-group set `G` and a complete target outer word `y`. Condition
on every outer output diagonal `D2`, then invert that diagonal and the fixed
invertible affine mixer. The resulting pre-`D1` array `v` is determined.
Its nonzero constituent rows determine a **unique** row-set vector before
the independent input diagonals are exposed. For a nonzero row with support
`w`, the probability that its inverse-diagonal image lies in RS[16,8] is
at most `255^(w-8)/255^w=255^-8`; for `w<9` it is zero.

Thus the expected counting measure for all source messages with
`A_tot>=A_0` is pointwise at most `255^(-8A_0)`. There is no union over
row-set vectors at this target: the inverse array already selected one.
Average over `D2` afterward. Relative to uniform output on the `r` active
groups, the whole high-arity stratum has density

\[
 \beta(r,A_0)=\frac{256^{256r}}{255^{8A_0}}.
\]

Equivalently, the stratum `D=16r-A_tot<=D_0` has density
`beta_16^r 255^(8D_0)`, with no binomial deficit factor. This is independent
of the detailed affine mask law and preserves the uniform packet interface.
It does not improve the forced-nonzero mask kernels for the complementary
low-arity stratum. The derivation and fractional chronology are recorded
separately in `HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md`.

### Saved-receipt headroom, not new G16 receipts

A read-only FLOAT calculation checked every 384-bit potential-band receipt
against the externally pinned closed-cover manifest, SHA-256
`accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.
No packet cap, alias atlas, or witness was rebuilt. For each entire saved
band ending at `r_hi`, it conservatively added the sharper-density
correction at `r_hi`. The following values estimate which complete near-full
strata could meet a log2 band target of `-60`. The old column uses the
superseded binomial factor `S_(r_hi)(D_0)`; the improved column uses only
`255^(8D_0)` under the unique-inverse-row-set argument.
They are not directed ceilings or evidence of a new construction receipt.

| Saved weight band | Corrected full-arity log2 bound | Old `D_0` | Improved `D_0` |
| --- | ---: | ---: | ---: |
| 458--566 | -98.214 | 5 | 6 |
| 1008--1040 | -105.548 | 5 | 7 |
| 1500--1571 | -85.154 | 3 | 3 |
| 1882--2153 | -95.365 | 4 | 5 |
| 2154--2274 | -1084.772 | 139 | 160 |
| 2275--2683 | +62.990 | none | none |
| 2684--2827 | -3129.054 | 425 | 479 |
| 2828--3285 | +65.009 | none | none |
| 3286--3325 | -8249.924 | 1155 | 1280 |
| 3326--3768 | -964.438 | 121 | 141 |
| 3769--4050 | -323.715 | 34 | 41 |
| 4051--4067 | -1591.322 | 207 | 239 |

Here `none` includes deficit zero.

Thus the comparison can cover useful broad near-full strata, but at these
witnesses it cannot bridge all weights even for full arity. Those positive
entries are failures of this conservative saved-band reuse, not lower
bounds on failure probability; per-weight corrections or new witnesses
could improve them. More importantly, thousands of low-/middle-arity
constituent assignments remain far outside this near-full stratum.

Any implementation needs an explicit G16 counting-measure provider and new
construction metadata. The uniform packet arithmetic may be reused after
its law and coordinates are identified, but old Toeplitz receipts must not
be relabelled. A new directed scalar correction with verified source pins,
or fresh regional replay, must bind the new density, deficit range, witness,
and support scope. Keeping some groups exact and making only the remaining
groups uniform would require a new hybrid packet provider. It is not
necessary for the simpler partition now available: compare the **whole**
high-total-arity stratum with uniform output and retain the mask-mass bound
for its complement.

For the complementary class `A_tot<=A_0-1`, a single deterministic marker
`0<t<=1` can retain total arity in the all-cohort mask bound. Multiply the
whole-word type-`a` mass by `t^a`, equivalently its per-column type factor
by `t^(a/16)`, and multiply the final bound by `t^(-(A_0-1))`. These factors
are outside `alpha`, since this marker restricts the fixed type assignment
union after fractionalization. The same `t` and type factors must apply
to every path and column. This proposed complementary provider requires
its own bounded control and directed replay; it is not implemented by the
unrestricted controls below.

### Why the constituent-only terminal threshold is intrinsically high

For RS[16,8], a prescribed support of size nine has a one-dimensional
shortened code and exactly 255 full-support words. Its membership ratio is
`255/255^9=255^-8`, attaining the per-row terminal bound. Thus an argument
that retains only `A_tot` cannot uniformly improve that row factor. The
threshold `A_tot>=122382` covers only the nearly full constituent tail
(about 93.4 percent of all 131072 constituent rows), not the whole large-`r`
range. A substantially earlier terminal cutoff must also retain row-support
weights or constraints on the inverse inner image. The next useful broad
test should target those missing low-/middle-arity strata with the
state-retaining interface, rather than accumulating near-full receipts.

## A regional majorant that makes every cohort factorize

Fix a weight band ending at `R` and witnesses `0<z<1`, `v=(1,p)` with
`0<p<=1`. For each possible regional support size `m<=R` and full-pair
count `f`, let `R[m,f]` be the forced-nonzero-amplitude two-state upper
matrix under the actual independently swapped route. Define

\[
 c_m=\max_f\max_{j\in\{0,1\}}
       \frac{(R[m,f]v)_j}{v_j}.
\]

Choose positive `A,B` satisfying `c_m<=A B^m` for every `m=0,...,R`.
For an arbitrary fixed sequence of regional masks, chronological
multiplication then gives

\[
 e_0\Bigl(\prod_{i=0}^{255}R[m_i,f_i]\Bigr)\mathbf1
 \le p^{-1}A^{256}B^{\sum_i m_i}.
\]

The matrices are not commuted. Each acts on the same positive potential,
and `1<=v/p`. Maximizing over feasible `f` is a universal domination of
each fixed source support, not a resampling of its matching geometry.

The exponent `sum_i m_i` is the total outer output support. It is additive
over source groups. This makes the ensuing comparison count factorize by
whole groups, while keeping each group's constituent row set fixed across
all sixteen columns.

## Three per-group counts, including all dense arities

Put `q=256`. For a one-row column, the affine mixer has full support;
its occupancy factor is `g_1=B^16`. For a two-row column the exact hole law
gives

\[
 g_2=\frac{239B^{16}+16B^{15}}{255}.
\]

Let `G_1` count one-row groups under the comparison and individual RS floor:

\[
 G_1=\frac{16}{255^8}
       \sum_{w=9}^{16}\binom{16}{w}(255g_1)^w.
\]

For two-row groups, retain both labelled row floors separately:

\[
 G_2=\frac{120}{255^{16}}
 \sum_{j,k\ge9}[x^jy^k]
    \bigl(1+255g_1(x+y)+255^2g_2xy\bigr)^{16}.
\]

The factors 16 and 120 count fixed constituent choices once. The positive
polynomial counts each source support array once; its coefficients do not
assert independence between the RS supports.

For dense arity `a`, retain the existing exact cardinality Hölder weight

\[
 W_a(w)^{16}=255^{-8a}
       \sum_{|A_g|=a}N_{A_g}(T)^{16},\qquad |T|=w.
\]

Its fixed-row-set count and dense normalizer are already included. Define

\[
 S_a(B)=\sum_{w=1}^{16}\binom{16}{w}W_a(w)B^w,
 \qquad
 G_{H,a}=\sum_{j=9}^{16}\binom{16}{j}W_a(0)^{16-j}S_a(B)^j.
\]

The floor counts nonempty columns of that one dense group. Invertibility
and the RS distance imply at least nine such columns. Summing labelled
target masks now contributes `C(16,w)` because the regional action has
already been replaced by `B^w`; unlike the earlier cardinality DP, this
formula has not otherwise enumerated those masks.

Set `G_*=G_1+G_2+sum_(a=3)^16 G_(H,a)`. Crucially, arities are summed
after each group's entire sixteen-column power. The factorization therefore
avoids the large loss observed when all arities were merged inside one
column Hölder weight. Expanding `G_*^r` automatically includes every cohort
coloring and every assignment of dense arities to the active groups.
No further cohort, row-set, or arity count is needed.

For any band `r_lo<=r<=r_hi<=R`, the proposed complete-family bound is

\[
 U\le \frac{A^{256}z^{-734003}}{255p}
        \sum_{r=r_{\rm lo}}^{r_{\rm hi}}\binom{8192}{r}G_*^r.
\]

One global projective normalization and one endpoint-potential factor
appear. Different group types may occur arbitrarily often. The only
new relaxation is the regional affine envelope, including its earlier
potential quotient; its cost may be too large even when the construction
is sound.

## Bounded implementation and next decision

`theory/hd_two_sum_gf16_factorized.py` implements this FLOAT interface.
Its default requested band is `18..32`, using the existing rational `z,p`
from the G16 dense gate. It selects only the one-dimensional affine slope,
not a new mark or state potential. The intercept is the maximum over every
computed occupancy. It rebuilds the actual forced-nonzero four-block packet
caps and positive coefficient powers; it does not import the old uniform
amplitude or fractional-power matrices. The source accepts maxima through
64, not an unproved arbitrary-size extension of the old guard17 wrapper.

Seven lightweight tests pass. They exhaustively check a noncommuting toy
factorization, both individual sparse row floors, once-only cohort counts,
unit-mark normalization, complete matching-histogram scope, the fixed-grid
budget, and reuse of three regional banks for nine witnesses. A timeout
remains incomplete and cannot become a coverage receipt.
An independent read-only audit found no substantive issue in the regional
majorant, fixed-row-set factorization, once-only counts, or forced-nonzero
coefficient provider.

The first complete `18..32` control finished in 1.023 seconds, including
0.994 seconds for its regional bank. It failed with FLOAT log2 union
`+27778.02961238569`. The selected envelope touched occupancies zero and
32, with `log2(A)=1.9847331600431062e-7` and
`log2(B)=-1.1834915764267318`. The input contractions were nearly flat:
`log2(c_1)=-38.65405609`, `log2(c_16)=-38.25767192`, and
`log2(c_32)=-37.87173025`, while `c_0=1`. Hence the exponential envelope
loses 37.47 bits at occupancy one and 19.32 bits at occupancy sixteen.
The full-arity per-group bound dominates at `log2(G_(H,16))=+723.5784153`.
Every individual weight fails, from `+17531.15` at 18 through `+27778.03`
at 32. This rejects that frozen-witness control, not the construction or
every possible affine envelope.

The inherited mark and potential were selected for a small dense cohort,
not a full all-cohort band. A separate authorized diagnostic therefore uses
exactly `z in {99/100,49/50,24/25}` and `log2(p) in {-64,-96,-128}`. Its
new wrapper is `theory/hd_two_sum_gf16_factorized_grid.py`. It preserves the
first helper and failed artifact, builds one regional bank per mark, and
enforces a 20-second budget for all nine fixed pairs. It fits neither `z`
nor `p`; only the same scalar envelope slope is selected.

That grid completed in 1.588 seconds with three regional banks. Every
witness failed. The best, `z=99/100`, `log2(p)=-64`, gave complete-band
FLOAT union `+16917.38270127202` and won at each individual weight. The
best values for marks `49/50` and `24/25` were respectively `+20986.51`
and `+41329.17`; both used log potential `-64`. Every selected envelope
again touched zero and 32. At the best witness, contractions range from
`log2(c_1)=-69.73652` to `log2(c_32)=-106.56905`. The envelope still loses
66.41 bits at occupancy one and 44.61 bits at occupancy sixteen. Full
arity remains dominant, with per-group log bound `+185.01250`.

These controls reject the tested affine compression, not the G16
construction. No 64-weight extension, broader fitting, or directed replay
was launched. The retained grid artifact
`hd_two_sum_gf16_factorized_r18_32_grid9_float.json` has SHA-256
`127dedbfab246d2a831913ec074b902e1b682b68fa2823b856d1b1fa1a7e47fb`.
The wrapper SHA-256 is
`905e81e23b8dd2ce298a0bda8d26141a4a2a7a2f3264750652aa5e40351e3a97`.

The first artifact is
`hd_two_sum_gf16_factorized_r18_32_frozen_float.json`, SHA-256
`c898f203b9656fc86cb41dc7acc05fa1d5a5b4118a989b15f9933209bfa93175`.
The helper SHA-256 is
`d34f054a6a1a07cec5ae97c75acdf475763a5c5e5ccab85d8be05cdd281e5da2`.
These pins identify FLOAT evidence, not directed arithmetic.

## Proposed next majorant: retain the nonempty-region penalty

The tested envelope draws a line through the exact empty contraction and
the largest occupancy. A different majorant can keep those cases separate.
Suppose `c_0<=1` and choose `0<C<=1`, `B>0` such that

\[
 c_m\le C B^m\quad(1\le m\le R).
\]

For one outer column, let `T_g` be its complete output mask for active
group `g`. The regional product contributes at most
`B^(sum_g |T_g|) C^(|union_g T_g|)`. Introduce an auxiliary allowed set
`J`, independently containing each of the sixteen regional labels with
probability `C`. The positive identity

\[
 C^{|\bigcup_g T_g|}
 =\Pr_J[\forall g:\ T_g\subseteq J]
\]

turns the correlated union penalty into a restriction applied to every
group. After the universal regional envelope, the conditional group weight
depends only on `k=|J|`, not on the positions in `J`.

For each deterministic group type `tau` (one row, two rows, or dense arity
3 through 16), let `F_tau(k)` be its positive per-column comparison under
that restriction. Distribute its once-only source count, row normalizer,
and a support-floor Chernoff factor equally over sixteen columns. For
example, for dense arity `a` and a fixed `y_a>=1`, take

\[
 F_a(k)=y_a^{-9/16}\left[
 W_a(0)+y_a\sum_{w=1}^k\binom kw W_a(w)B^w\right].
\]

For singleton and pair types, use markers `x_1,x_2>=1` and

\[
 F_1(k)=16^{1/16}255^{-1/2}x_1^{-9/16}
       (1+255x_1B^{16}1_{k=16}),
\]
\[
 F_2(k)=120^{1/16}255^{-1}x_2^{-18/16}
 \left[1+510x_2B^{16}1_{k=16}
 +255x_2^2(239B^{16}1_{k=16}+H_kB^{15})\right],
\]

where `H_16=16`, `H_15=1`, and other `H_k=0`. The second sparse marker
retains the necessary aggregate floor 18 rather than both separate floors.
The markers are deterministic proof witnesses, not setup-dependent choices.

Put `pi_k=C(16,k)C^k(1-C)^(16-k)`. Before any relaxation of group types,
the complete fixed-type union over sixteen columns is at most

\[
 \mathbb E_{k_1,\ldots,k_{16}}
       \left(\sum_\tau\prod_{j=1}^{16}F_\tau(k_j)\right)^r.
\]

Hölder across these sixteen column positions bounds this expression by

\[
 \left[\sum_{k=0}^{16}\pi_k
       \left(\sum_\tau F_\tau(k)^{16}\right)^{r/16}\right]^{16}.
\]

Multiply by `C(8192,r) z^-734003/(255p)` and sum the desired weights.
The type choices remain fixed before applying Hölder. Their counts occur
once, and no `A^256` intercept remains. This formula has seventeen cut-size
states regardless of the number of dense groups.

### Bounded cut control and the remaining obstruction

`theory/hd_two_sum_gf16_cut_gate.py` implements this formula at the grid's
fixed `z=99/100`, `log2(p)=-64`. It uses the pinned FLOAT regional cap
vector; a certificate would require fresh directed reconstruction. Both
`c_0<=1` and `C<=1` are checked explicitly, without discarding an empty
surcharge. Five shared marker schedules use
`log2(x_1)=log2(x_2)=log2(y_a)` in `{0,1,4,8,16}`. Each schedule permits at
most 40 scalar slope-search evaluations, checks both endpoints, and
reconstructs the selected bound once. There is no new mark/potential fit.

The complete control finished in 0.999 seconds. All five schedules failed.
Their complete-band FLOAT log2 bounds were respectively
`+16580.999482`, `+16804.999482`, `+17476.999482`, `+18372.999482`, and
`+20164.999482`. The best result is about 336.38 bits below the earlier
affine bound but remains far above the target. Its envelope has

```text
log2(B) = -1.389715421093215
log2(C) = -62.09815911763466
contacts: m=5 and m=32
```

At weights 18–21 the dominant allowed-set term has `k=0` and singleton
type: the unmarked comparison admits empty source columns. At weights
22–32 the dominant term has `k=16` and full dense arity. At `r=32`, this
second term also dominates the complete union. Increasing a shared floor
marker suppresses the empty term but increases the already dominant
full-column term; it therefore does not improve this grid's union.
This is a failure of the tested scalar contraction and floor witnesses,
not a construction obstruction or a lower bound on actual setup failure.

Five further exact tests verify the auxiliary-set identity, fixed-type
Hölder expansion, heterogeneous-column Chernoff direction, once-only source
counts and normalizers, the `C=1` boundary, and the empty-surcharge guard.
Independent read-only reviews of the concept and source found no substantive
issue. No further sweep or directed replay of this failed gate was launched.

The cut artifact `hd_two_sum_gf16_cut_r18_32_fixed_float.json` has SHA-256
`929d01fefd706fc6db194ff071ef68b60ff3fe0200d2ab9b4dfb6ae516227ed1`.
Its source has SHA-256
`745f79821703c47393a09b5348a88d600a9ac852241656ced5a68be0d9bbe544`.
The three isolated suites total twelve lightweight tests.

## State-retaining follow-up and the next total-arity split

The state-retaining follow-up is now implemented and remains negative;
see `HD_TWO_SUM_GF16_STATE_ENVELOPE.md`. It majorizes the four regional
entries separately and sums 274 ordered state-path classes before the
sixteen-column Hölder step. The subsequent fractional-mask control is
documented in `HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md`. Its 25 initial witnesses
retain five fixed `(z,alpha)` controls and five deterministic type-specific
floor templates. Type markers stay fixed across all paths and columns.

The fractional control exposed a large initialization defect: the
live-to-live entry peaks near occupancy two, but its initial endpoint
secant had positive slope. At `(z,alpha)=(1/2,1/16)`, that envelope loses
460.31 bits at occupancy 32. The isolated saved-bank repair
`theory/hd_two_sum_gf16_fractional_repair.py` tested at most sixteen existing
live-state slope candidates for every frozen baseline, holding its other
three slopes, mark, exponent and all floor markers fixed. It then allowed
one four-coordinate round for the globally best repaired whole-band
witness. It rebuilt no packet operators and read saved caps only as FLOAT
proposals. Fresh dense moments were recomputed from each exponent.

The repair completed in 6.097 seconds. It reduced the best fractional-bank
bounds from `+43191.32,+72866.55,+112457.52,+143997.52` to respectively
`+21288.22,+27867.34,+38798.07,+47486.72`, but all 25 repaired controls still
failed. Alpha one remained the best bank, and its single extra coordinate
round left the complete `18..32` FLOAT union at `+16494.53668044997`.
Thus correcting the diagnosed large slack does not close this unrestricted
all-cohort control. This is not a lower bound on the construction's failure
probability or a proof that another marked partition cannot close it.

The repaired artifact `hd_two_sum_gf16_fractional_repair_r18_32_float.json`
has SHA-256
`592b49e9667c519de712b7825984931400c5620be39df9534a634482b611f774`.
The repair source has SHA-256
`83fc3f5b3320f43d1c6523c5dbacc5d6d0a4751d323a3fc9b9153a3947934c4d`.
Its six tests plus the fractional and state suites total 24 passing
lightweight tests. Independent read-only review found no substantive issue.

The next priority is the disjoint **high-/low-total-arity** split proved
above, not another unrestricted mark grid. Saved FLOAT headroom at the
inherited uniform-comparison witnesses suggests the following first anchors:

| Active groups `r` | Exact mark `z` | `alpha` | Proposed high-arity cutoff `A_0` at a `-60` target |
| ---: | ---: | ---: | ---: |
| 18 | `691/700` | `1/5` | 247 |
| 24 | `172/175` | `1/5` | 317 |
| 32 | `171/175` | `1/5` | 463 |

These are parent-reported saved-bound estimates, not fresh directed G16
receipts. The high class would use the no-row-choice uniform density
`256^(256r)/255^(8A_0)`. The low class would retain `A_tot<=A_0-1` with the
single `t<=1` type marker outside `alpha`. Both providers must be bound to
the new construction and replayed before their disjoint sum can enter a
coverage ledger. No such implementation or further sweep ran in this cycle.

Reproduce the recorded controls only with fresh output filenames:

```text
python -B -m unittest hd_two_sum_gf16_factorized_tests hd_two_sum_gf16_factorized_grid_tests hd_two_sum_gf16_cut_gate_tests
python -B hd_two_sum_gf16_factorized.py --output NEW_FACTOR.json
python -B hd_two_sum_gf16_factorized_grid.py --output NEW_GRID.json
python -B hd_two_sum_gf16_cut_gate.py --output NEW_CUT.json
```

For genuinely large `r`, direct truncated bivariate tables will eventually
be too expensive. The established paired-route positive/alias coefficient
machinery can supply regional boxes, followed by a directed affine envelope
over those boxes. That requires new forced-nonzero coefficient providers and
normalizer bounds; the old uniform-input numerical tables cannot be reused.
Together with the terminal `A_tot` tail and high-total-arity uniform
comparison, this would yield a coverage plan by broad bands rather than
dense-cohort histograms.
