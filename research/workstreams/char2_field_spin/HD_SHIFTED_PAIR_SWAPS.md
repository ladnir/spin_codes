# Alternating supports: an obstruction and a swap repair

The fixed-lane alternating-pair route has a support class that the current
density comparator cannot certify, for any Chernoff mark. Fresh independent
swaps inside routed pairs remove this particular obstruction. The statements
below distinguish the comparator, a deterministic subcode, and a probability
bound for the modified ensemble. The early local arguments are not whole-code
certificates. The completed [disjoint cover](HD_SHIFTED_CLOSED_CERTIFICATE.md)
below now establishes 35% distance with 55.09 bits of setup-failure margin.

Set $q=256$, $K=2^{20}$, $N=2^{21}$, and

\[
L=\lfloor 7N/20\rfloor=734003,\qquad
T=N/64=32768,\qquad
\beta=\frac{256^{256}}{255^{128}}.
\]

Write $e_0=(1,0)$ for the zero-state row and $\mathbf1=(1,1)^T$.
The two matrix states distinguish zero from nonzero refresh state.

There are 8192 outer groups, each with 128 input coordinates and 256 output
coordinates. For a binary cyclic support on these groups, let $r$ count
active groups. Let $a,b$ count active-active pairs in the even and odd
matchings. The class $(r,a,b)=(4096,0,0)$ contains exactly two supports:
the two alternating binary words. Indeed, an independent set of size half
the even cycle must have exactly one inactive vertex between successive
active vertices.

## Fixed-lane comparator

Fix either alternating support. Within each matching, every routed pair
contains one active group, always in the same lane. Every 64-input packet
therefore has potentially nonzero inputs in two complete 16-coordinate
blocks of the local h4 map. Individual outer coordinates may still vanish.
The occupied lanes switch between matchings, but the sorted occupancy is
always $(0,0,16,16)$.

For $0<z<1$, let $M(z)$ be the existing nonnegative two-state packet
cap for this occupancy, before its entrywise fractional power. Its
zero-state-to-zero-state entry is

\[
k(z)=\frac{W(z)^2}{256^{32}},\qquad
W(z)=\sum_{w=0}^{16} A_wz^w,
\]

where $A_w$ counts weight-$w$ words in the q-ary MDS kernel
with parameters $[16,12]$. These counts are exact:

\[
A_w=\binom{16}{w}B_w,\qquad B_0=1,
\]

with $B_w=0$ for $1\le w<5$, and

\[
B_w=255\sum_{j=0}^{w-5}(-1)^j\binom{w-1}{j}256^{w-5-j}
\quad (w\ge5).
\]

Define $d=[z^{22}]k(z)$. The exact value is

\[
d=\frac{19208186874811647000306904575190818022875}
{7237005577332262213973186563042994240829374041602535252466099000494570602496}.
\]

All coefficients of $k$ are nonnegative. Since $22T=720896\le L$,

\[
z^{-L}\,e_0M(z)^T\mathbf1
\ \ge\ z^{-L}k(z)^T
\ \ge\ d^T z^{22T-L}
\ \ge\ d^T.
\]

Consequently, the numerical comparator for this fixed support satisfies

\[
\frac{\beta^{4096}}{255}\,
z^{-L}e_0M(z)^T\mathbf1
\ \ge\ \frac{\beta^{4096}}{255}d^T.
\]

Directed 256- and 384-bit arithmetic gives the same rational lower bound
on its base-two logarithm:

\[
\log_2\!\left(\frac{\beta^{4096}}{255}d^T\right)
\ \ge\ \frac{81174087123381731}{250000000000}
=324696.348493526924.
\]

This inequality holds for every $z\in(0,1)$. Raising the entire
comparator to any positive exponent cannot make it less than one.
Entrywise fractional powers only enlarge the corresponding path-sum bound.
Exact orientation handling, larger alias grids, and higher arithmetic
precision therefore cannot repair this comparator on the fixed-lane class.
The lower bound is on an upper-bound expression, not on actual failure
probability. It does not refute relative distance $0.35$.

## A deterministic subcode of the fixed-lane construction

There is also a structural statement that does not use the comparator.
Fix all sampled maps and routes in the original construction. Restrict
input messages to one alternating half of the outer groups. This input
space has dimension $K/2$ over GF(256).

The outer map acts independently on each group. Hence all other groups
remain zero in every region. Packet routing places the possibly nonzero
inputs in two of the four local blocks. Each occupied block contributes
at most four independent feedback equations. Impose zero feedback on
every packet. These are at most $8T$ homogeneous linear constraints, so
their common kernel has dimension at least

\[
K/2-8T=K/4=262144.
\]

The initial state is zero. Zero feedback keeps the state zero under the
common refresh recurrence. On this kernel, the inner output equals the
routed outer input, up to the fixed physical output permutation.
The overall map remains injective there because the outer map is injective.
Thus the code contains a subcode of dimension at least $K/4$, supported
on a fixed set of at most $N/2$ output coordinates.

The native layout agrees with these statements:

- `HdGroupedOuter.h`, `coset` and `mixStore`, preserve each 128-bit lane as
  one independent outer group; the group mixer does not mix these lanes.
- `hd_grouped_forward.cpp`, the expanded pair-route assignment, maps a
  fixed alternating support to logical packet lanes 0 and 2, or 1 and 3.
- `HdGroupedInner.h`, `innerCoset`, applies four checks independently to
  each such lane. Its comment $C=A^T$, $CA=0$ also permits the fused
  feedback implementation to consume the output instead of the input.
- `finishStriped` maps zero state plus zero feedback to zero. The optional
  fixed output permutation does not affect support size.

This argument concerns the original width-two route without random swaps.
For the SIMD implementation, one may restrict to one GF(256) component of
each record; the same coordinate-support statement then applies to records.
The Singleton bound yields a nonzero word of weight at most

\[
N/2-K/4+1=3N/8+1=786433.
\]

This deterministic upper bound still exceeds the target $L$. A heuristic
comparison with random q-ary codes on the half-length support is not a
distance theorem and is not used in any receipt.

## Fresh independent within-pair swaps

Consider a modified ensemble. After each regional matching and pair
permutation, independently swap the two lanes of each routed pair with
probability $1/2$. All swap bits are fresh across pairs and regions,
and independent of the original routing, outer maps, and refresh maps.
The existing ideal-randomness assumptions for those maps remain unchanged.
These assumptions are not supplied merely by a reproducible native PRNG.

For either alternating support, every routed pair is a singleton.
In a packet, write $u,v\in\{0,\ldots,16\}$ for the numbers of active
first lanes in the first and second halves. Their joint probability is

\[
p_{u,v}=2^{-32}\binom{16}{u}\binom{16}{v}.
\]

The four block occupancies are $(u,16-u,v,16-v)$. Let

\[
A_\alpha(z)
=\sum_{u,v=0}^{16}p_{u,v}
\left(M_{\operatorname{sort}(u,16-u,v,16-v)}(z)
\right)^{\circ\alpha},\qquad 0<\alpha\le1,
\]

where $\circ\alpha$ denotes an entrywise power. Each matrix on the
right is the existing exact rational packet cap, rounded upward to a
60-bit dyadic mantissa before exponentiation.

The order of operations matters. First fix the swaps and apply the
fractional-moment bound. Subadditivity bounds the fractional power of
the state-path sum by products of entrywise fractional packet caps.
Then average over independent swaps. This gives $A_\alpha(z)^T$,
not an entrywise fractional power of the averaged raw cap.

Let $U_{\rm swap}$ be the probability, over the modified ensemble's setup,
that some message with either alternating active-group support has output
weight at most $L$. The existing outer density and projective normalization
give

\[
U_{\rm swap}
\ \le\ 2\left(\frac{\beta^{4096}}{255}z^{-L}\right)^\alpha
e_0A_\alpha(z)^T\mathbf1.
\]

The projective factor $1/255$ is inside the fractional power. The
support count two is outside it. At $z=1/32$, $\alpha=1/4$, both
256- and 384-bit directed replay return

\[
\log_2 U_{\rm swap}
\ \le\ -\frac{95969876597769549}{10^{12}}
=-95969.876597769549.
\]

## Inheritance of unsigned-envelope bounds

Fix $z,\alpha$. Use the same nonnegative singleton and full-pair markers
on both sides of each comparison below. Let $C_{u,v,f}$ be the existing
oriented one-step coefficient: it sums entrywise fractional packet caps
over all labeled arrangements with the indicated singleton and full-pair
counts. Conditional on the unsigned arrangement, fresh swaps assign each
of its $s=u+v$ singleton orientations independently and fairly. Thus

\[
B^{\rm swap}_{s,f}
=2^{-s}\sum_{u=0}^{s} C_{u,s-u,f}
=\sum_{u=0}^{s}\frac{\binom{s}{u}}{2^s}
  \frac{C_{u,s-u,f}}{\binom{s}{u}}
\ \le\ B^{\rm max}_{s,f},
\]

entrywise, where $B^{\rm max}_{s,f}$ is the old entrywise orientation
maximum. The weights in the middle expression sum to one. The unsigned
multinomial normalizer and the actual support counts do not change.

This comparison also preserves the existing good-pair improvement.
Write $G$ for the good one-step polynomial, with at most one singleton,
and $D^{\rm swap},D^{\rm max}$ for the dirty one-step polynomials.
Good-step packet caps are invariant under the lone singleton's orientation,
so $G$ is unchanged. Let $P$ denote the good-good polynomial, whose
ordered raw two-step matrices are composed before the entrywise alpha.
Each of its good-step matrices is orientation-invariant; hence $P$ also
remains unchanged after averaging the independent swap bits. Therefore

\[
\begin{aligned}
\mathcal T^{\rm swap}
 &=P+GD^{\rm swap}+D^{\rm swap}G+(D^{\rm swap})^2\\
 &\le P+GD^{\rm max}+D^{\rm max}G+(D^{\rm max})^2
 =\mathcal T^{\rm old}.
\end{aligned}
\]

All coefficients are nonnegative. Ordered matrix products, the regional
power $(\mathcal T)^{64}$, coefficient extraction, and positive alias sums preserve
this inequality. Existing directed unsigned-envelope bounds consequently
remain valid for the swap ensemble, with their original scalar factors.
This is an explicit inheritance argument, not an identification of the
two constructions.

Fixed-oriented receipts do not transfer through this lemma: swaps randomize
the oriented totals. Nor does the lemma justify averaging raw caps before
the fractional power. The new endpoint receipt and any inherited partial
coverage still leave a whole-code coverage obligation.

## Replay and next step

`theory/hd_shifted_endpoint.py` contains the exact coefficient calculation,
the swap law, and the two directed receipts. No FFT or floating-point
proposal is used. `theory/hd_shifted_endpoint_tests.py` checks small swap
laws exhaustively, the two-support count, the alpha chronology, the exact
quiet-polynomial identity, and precision corroboration.
The independent `theory/hd_shifted_swap_inheritance_tests.py` adds five
checks of fair-label multiplicities, good-good invariance, noncommuting
positive products, aliases, and fixed-oriented noninheritance.

From the repository root:

```text
python research/workstreams/char2_field_spin/theory/hd_shifted_endpoint_tests.py
python research/workstreams/char2_field_spin/theory/hd_shifted_endpoint.py --precision 384
python research/workstreams/char2_field_spin/theory/hd_shifted_swap_inheritance_tests.py
```

All six tests pass. The next proof step is to use the swap-averaged
unsigned operator on the unresolved support classes. Any native candidate
must preserve the fresh per-pair, per-region swap law stated above, and
must remain separate from certificates for the original fixed-lane route.

## General swap-averaged operator

`theory/hd_shifted_swap_certificate.py` implements the coefficient average
and the `SwapCachedBound` provider. It reuses the original box geometry,
count bounds, convex normalizers, positive evaluation, and interval aliases.
Only the one-step orientation aggregation changes. The ordered good-good
polynomial remains the original one.

The coefficient cache is separate from the fixed-lane cache and includes
the exact mark, exact alpha, Arb precision, and packet rounding precision.
The backend exposes an overridable coefficient-provider classmethod; its
default still calls the original provider. Swap receipts use construction
identifier `toeplitz-gf256-h4-alternating-swapped-pairs-k20-rate-half`.
Their witnesses explicitly record independent fair swaps.

The provider supports positive bounds, local alias maxima, uniform alias
maxima, and optional run-count markers. Its regional power still uses
64 ordered two-step factors, followed by 128 alternating regional pairs.
Neither the histogram normalization nor the projective normalization changes.

Two general replay controls now supplement the exact endpoint receipt:

- The positive-marker box $4064\le r\le4128$, $0\le a,b\le32$ contains
  47,377 feasible unsigned profiles. At $z=1/32$, $\alpha=1/4$, and regional
  markers $(512,1)$, its union logarithm is at most
  $-56518.096993211155$. Directed 256- and 384-bit replay agree.
- The general profile $(r,a,b)=(16,4,4)$ has a directed grid-64 alias bound
  of $-144.174076029899$ bits. Its exact rational witness is preserved in
  `hd_shifted_swap_certificate_tests.py`. This is an all-lane-balance
  profile receipt, not coverage of every weight-16 profile.

The first box is saved as `BOUNDARY_BOX`; `boundary_receipt()` replays it.
Run the general-provider checks and the saved box from the repository root:

```text
python research/workstreams/char2_field_spin/theory/hd_shifted_swap_certificate_tests.py
python research/workstreams/char2_field_spin/theory/hd_shifted_swap_certificate.py --precision 384
```

The original backend tests still pass. These controls do not fill the
remaining coverage gaps. The next step is to combine fresh swap receipts
and explicitly inherited unsigned receipts in the construction-specific
ledger, retaining unresolved profiles until they receive directed bounds.

## Normalize before taking projected maxima

A box bound can retain the cancellation between a coefficient and its
histogram normalizer. Fix a feasible unsigned profile box $\mathcal B$,
positive run markers $x,y$, and the existing $z,\alpha$. Let $R_p(r,f)$
be a certified nonnegative regional matrix for parity $p\in\{0,1\}$.
This matrix includes division by

\[
\binom{4096}{4096-r+f,\ r-2f,\ f}
\,\mu_{p,s}^{\,r-2f}\mu_{p,f}^{\,f},
\]

where $\mu_{p,s},\mu_{p,f}>0$ are that region's coefficient markers.
The numerator may come from positive evaluation or a directed alias bound.

For each parity, project $\mathcal B$ onto its integer $(r,f)$ histograms.
Including additional histograms with nonnegative entries only weakens the
bound. Define

\[
\widetilde R_p(r,f)
=R_p(r,f)\,
\beta^{\alpha r/256}
x^{-r/256}y^{-r/256+f/128},
\qquad
E_p=\max_{(r,f)\text{ in the projection}}\widetilde R_p(r,f),
\]

with an entrywise maximum. For every $(r,a,b)\in\mathcal B$, positivity
preserves the two regional comparisons through ordered composition:

\[
\frac{\beta^{\alpha r}}{x^r y^{r-a-b}}
\,e_0\bigl(R_0(r,a)R_1(r,b)\bigr)^{128}\mathbf1
\ \le\ e_0(E_0E_1)^{128}\mathbf1.
\]

The scalar identity uses 128 occurrences of each parity. Independent
entrywise maxima may select inconsistent histograms or different entries
from different histograms; these choices only enlarge the positive product.
No independence between the two support histograms is asserted.

Let $\mathcal U(\mathcal B)$ denote the indexed binary supports in the box,
and let $c(U)=r(U)-a(U)-b(U)$. Their total tilted count is bounded by

\[
\sum_{U\in\mathcal U(\mathcal B)}x^{r(U)}y^{c(U)}
\ \le\ G(x,y)
:=\operatorname{tr}
\begin{pmatrix}1&xy\\1&x\end{pmatrix}^{8192}.
\]

Thus their failure union is bounded by

\[
255^{-\alpha}z^{-\alpha L}\,
G(x,y)\,e_0(E_0E_1)^{128}\mathbf1.
\]

The exact transfer sum $G$ counts indexed supports once. An additional
three-dimensional profile factor would double-count this union. The
implementation enumerates only the two projected histogram sets, not all
triples in the box.

The unsigned domination lemma also survives this normalization: all
divisors and tilts are positive and unchanged between the two ensembles.
Consequently, an old unsigned projected bound can be explicitly inherited
by the swap construction. This does not transfer fixed-oriented receipts.
`theory/hd_shifted_projected_bound.py` implements the bound; its algebra
tests check the distributed tilt and ordered noncommuting products.

## Directed coverage of the complete 4096..4351 band

The saved plan `theory/hd_shifted_swap_bulk_plan3.json` now has independent
directed replays at 256-bit and 384-bit Arb precision. Both accept all
3,937 proposed leaves with zero failed leaves and identical ledger summaries.
The plan uses 190 witnesses across 24 distinct weight-marker/exponent pairs.
Replay reconstructs every directed bound; it does not certify the saved
FLOAT values by relabeling them.

For the ideal swapped ensemble defined above, let $U_{4096:4351}$ be the
probability that some nonzero message with $4096\le r\le4351$ has output
weight at most $L=734003$. Both replays give the exact ceiling

\[
\log_2 U_{4096:4351}
\ \le\ -\frac{48297242508857}{500000000000}
=-96.594485017714.
\]

The band contains 1,009,047,040 feasible unsigned profiles $(r,a,b)$,
including every even/odd split. All are covered. The ledger retains 902
unclaimed geometric boxes, but they contain zero feasible profiles.
Its volume and profile accounting checks pass, and
`complete_local_cover=true`. The local root is only this weight band;
`whole_code_certificate=false` remains unchanged.

The replayed plan has SHA-256

```text
401bce5f667df40cc8431ec50519cfa6af02e450c10df53dda11839360aa7339
```

From `research/workstreams/char2_field_spin/theory`, this command verifies
the hash and prints only each replay summary and failed-leaf count:

```text
python -B -c "import json; from pathlib import Path; from hashlib import sha256; from hd_shifted_swap_ledger import replay_swap_plan; data=Path('hd_shifted_swap_bulk_plan3.json').read_bytes(); assert sha256(data).hexdigest()=='401bce5f667df40cc8431ec50519cfa6af02e450c10df53dda11839360aa7339'; plan=json.loads(data)
for precision in (256,384):
 ledger,result=replay_swap_plan(plan,precision=precision); print(json.dumps(dict(precision=precision,summary=result['summary'],failed_count=len(result['failed']))),flush=True)"
```

The observed local replay times were 163.7 seconds and 177.3 seconds,
respectively. These are proof replays, not encoder benchmarks. Combining
this result with earlier partial plans requires disjoint band assignments:
use this receipt for 4096..4351 and restrict overlapping partial coverage
to its remaining domains. No combined-cover replay or whole-code claim
is supplied by this section.

## Common state potential and cyclic runs

The next reduction keeps the correlation between the two matchings but
replaces each regional matrix by a scalar bound. Fix a nonzero scalar
support weight $r$, a weight marker $0<z<1$, and $0<\alpha\le1$.
Let $R_{r,f}$ be the normalized regional operator when a matching has
$f$ full pairs. Choose one positive two-entry vector $v$, common to every
$f$ at this weight, and set

\[
c_{r,f}=\max_{i\in\{0,1\}}
  \frac{(R_{r,f}v)_i}{v_i}.
\]

For a support whose two matching counts are $a,b$, positivity gives

\[
e_0(R_{r,a}R_{r,b})^{128}\mathbf1
\le \frac{v_0}{\min_i v_i}\,c_{r,a}^{128}c_{r,b}^{128}.
\]

Using a different vector for each $f$ would invalidate this ordered-product
bound. The regional matrices need not commute.

An indexed support on the 8192-position cycle has exactly $r-a-b$ runs.
Write $C(n,r,j)$ for the number of indexed binary cycles of length $n$
with weight $r$ and $j$ runs. For a nonconstant cycle,

\[
C(n,r,j)=\frac{n}{j}\binom{r-1}{j-1}\binom{n-r-1}{j-1}
\quad (1\le j\le\min(r,n-r)),
\]

and it is zero outside that range. This counts indexed supports, not
rotation classes. It upper-bounds the number of supports with any
specified matching counts $a,b$ whose sum is $r-j$.

Define $P_r(X)=\sum_f c_{r,f}^{128}X^f$. Summing all matching profiles
therefore gives the upper bound

\[
\Phi_r=\frac{v_0}{\min_i v_i}
\sum_d C(8192,r,r-d)[X^d]P_r(X)^2.
\]

The polynomial square sums products of scalar bounds; it does not assume
independent matching counts. In particular, the impossible case $a=b=r/2$
has zero contribution for $0<r<8192$. This avoids the main loss in the
unweighted Cauchy--Schwarz bound, which uses only the two marginal counts.
With the existing outer-density factor $\beta=256^{256}/255^{128}$ and
$L=734003$, the failure probability for this entire support-weight class is
at most $255^{-\alpha}\beta^{\alpha r}z^{-\alpha L}\Phi_r$.

Implementation: `theory/hd_shifted_run_potential.py` performs the positive
Arb polynomial convolution. `theory/hd_shifted_swap_potential.py` builds
directed regional coefficient bounds. Several positive marker choices
may bound the same regional operator: normalize each bound at its own
$(r,f)$ histogram, then take the entrywise minimum. All choices retain
the same $z,\alpha$ and setup distribution. Markers and $v$ are proof
witnesses, not encoder parameters.

The exact witnesses are saved as `SAVED_RUN_WITNESSES`. Independent
256- and 384-bit replays give identical rational ceilings:

| Active groups $r$ | Feasible matching profiles | Log2 failure upper bound |
|---:|---:|---:|
| 128 | 4,224 | -360.764721228035 |
| 512 | 66,048 | -831.900235619398 |
| 1024 | 263,168 | -81.573243260603 |

These rows cover every indexed support and every lane balance at the
stated weight. They do not cover intervening weights. The $r=128$ row
overlaps the exact sparse certificate and must not be counted as a new
disjoint band. Reproduce from the theory directory with

```text
python -B hd_shifted_swap_potential.py --weights 128 512 1024 --precision 256
python -B hd_shifted_swap_potential.py --weights 128 512 1024 --precision 384
```

Exact toy tests check correlated matching counts, the common-vector
ordered-product inequality, the run convolution, the excluded zero-run
case, and normalization before taking coefficient minima. No new native
benchmark or construction change was needed for this proof reduction.

### A complete neighboring-weight band

The stronger frozen frontier witness at $r=1024$ also closes every integer
weight $1008\le r\le1040$. The witness uses
$z=0.4828132651072611$, $\alpha=0.10118724647150915$, 22 exact positive
coefficient-marker pairs, and one exact positive state potential.
These are fixed across the band. They differ from the older single-weight
1024 witness in the table above; no frozen witness was overwritten.

`theory/hd_shifted_potential_band.py` builds each directed alias table once,
normalizes its entries for every requested $(r,f)$, and applies the run
convolution separately at every integer $r$. The final union sums the
33 probability bounds directly. Both 256- and 384-bit replays give

\[
\log_2 U_{1008:1040}
\le -\frac{105547875859763}{10^{12}}
=-105.547875859763.
\]

All 33 per-weight rational ceilings agree between precisions. The band
contains 8,677,104 feasible matching profiles. Its center-weight ceiling
is $-165.868870888917$. An independent read-only audit checked chronology,
normalization, the common potential, precision contexts, and the disjoint
weight union. This is complete coverage only of the stated band.

The exact proposed witness is retained in
`theory/hd_shifted_potential_band_1024_frontier.json`, SHA-256
`23aee6ac737c2d2a2d5c3d4205775f27bcd1da91b537bdff03748f56a1949d99`.
Replays ignore all saved FLOAT bounds and reconstruct the operators:

```text
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_1024_frontier.json --precision 256
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_1024_frontier.json --precision 384
```

The construction-specific ledger accepts this proposal through
`--potential-band` and rebuilds it before acceptance. Assigned domains
must be disjoint. In particular, do not also import the older 1024
single-weight receipt or the overlapping projected checkpoint.

### Use the actual Toeplitz outer bound

The receipts above retain the older diagonal-Cauchy density bound. The
implemented outer uses four RS[64,32] constituents and independent uniform
invertible Toeplitz matrices on their four-symbol columns. This ensemble
admits a stronger bound without changing the encoder.

Fix a nonzero desired group output with $h$ nonzero columns. If $h<33$,
the output is impossible. Otherwise shorten the four constituents to
these columns and condition on the mixers at any $h-32$ of them. These
information columns determine at most one source group. Exactly 32
nonzero target columns remain. Their mixers are independent of the
conditioning. A zero source column makes the event impossible; otherwise
each mixer hits its target with probability at most

\[
\rho_T=\frac{1}{(q-1)q^3},\qquad q=256.
\]

The expected output measure at each point is therefore at most
$\rho_T^{32}$. Its density relative to uniform $\mathbb F_q^{256}$ is

\[
\beta_T=q^{256}\rho_T^{32}
       =q^{128}(1-1/q)^{-32}
       =\frac{256^{160}}{255^{32}}.
\]

Replacing the inherited $\beta$ by $\beta_T$ saves
$96\alpha r\log_2(256/255)$ bits in a weight-$r$ bound. At $r=1536$ and
$\alpha=0.102$, the saving is about 84.93 bits. The packet operators,
coefficient markers, support counts, and projective factor $1/255$ are
unchanged. In particular, this is not an improvement in measured encoding
time or a change to the setup distribution.

`hd_shifted_outer_density.py` exposes the exact two established bounds.
The old `inherited-diagonal` mode remains the default. New proposals must
explicitly select `invertible-toeplitz`; directed replay reconstructs
their operators and evaluates the sharper prefix at its working precision.
The sharper bound does not apply to the separate XOR-column ensemble.
Old artifacts are retained rather than relabeled.

### Disjoint coverage checkpoint

The following intervals cover every scalar support at every stated
integer weight. Each row has matching 256- and 384-bit directed replays.
Margins below are rounded down; they bound the union over that row, not
the whole construction. Older overlapping single-weight and narrow-band
receipts must not be added again.

| Active groups | Margin at least, bits | Proof path |
|---|---:|---|
| 1--256 | 55.36 | Exact sparse coefficients |
| 257--279 | 446.32 | Fixed 272 bank, sharper Toeplitz density |
| 280--349 | 86.24 | Fixed 320 bank |
| 350--433 | 59.11 | Fixed 400 bank, sharper Toeplitz density |
| 434--457 | 761.68 | Fixed 448 bank, sharper Toeplitz density |
| 458--566 | 98.21 | Fixed 512 bank |
| 567--633 | 320.88 | Fixed 608 bank, sharper Toeplitz density |
| 634--699 | 301.48 | Fixed 672 bank, sharper Toeplitz density |
| 700--809 | 82.29 | Fixed 768 bank |
| 810--896 | 172.76 | Fixed 864 bank, sharper Toeplitz density |
| 897--974 | 147.29 | Fixed 944 bank, sharper Toeplitz density |
| 975--1007 | 143.35 | Fixed 992 bank |
| 1008--1040 | 105.54 | Fixed 1024 frontier bank |
| 1041--1088 | 269.56 | Fixed 1064 bank, sharper Toeplitz density |
| 1089--1129 | 242.89 | Fixed 1108 bank, sharper Toeplitz density |
| 1130--1199 | 60.45 | Fixed 1168 bank, sharper Toeplitz density |
| 1200--1255 | 58.82 | Fixed 1232 bank, sharper Toeplitz density |
| 1256--1325 | 90.99 | Fixed 1280 bank, sharper Toeplitz density |
| 1326--1383 | 86.50 | Fixed 1352 bank, sharper Toeplitz density |
| 1384--1441 | 118.54 | Fixed 1412 bank, sharper Toeplitz density |
| 1442--1499 | 152.35 | Fixed 1472 bank, sharper Toeplitz density |
| 1500--1571 | 85.15 | Fixed 1536 bank |
| 1572--1645 | 243.53 | Fixed 1608 bank, sharper Toeplitz density |
| 1646--1824 | 61.99 | Fixed 1728 bank, sharper Toeplitz density |
| 1825--1881 | 785.11 | Fixed 1856 bank, sharper Toeplitz density |
| 1882--2153 | 95.36 | Fixed 2048 bank |
| 2154--2274 | 1208.03 | Fixed 2208 bank, sharper Toeplitz density |
| 2275--2683 | 82.44 | Fixed 2496 bank, sharper Toeplitz density |
| 2684--2827 | 3282.29 | Fixed 2752 bank, sharper Toeplitz density |
| 2828--3285 | 113.06 | Fixed 3072 bank, sharper Toeplitz density |
| 3286--3325 | 8430.16 | Fixed 3304 bank, sharper Toeplitz density |
| 3326--3768 | 1168.68 | Fixed 3584 bank, sharper Toeplitz density |
| 3769--4050 | 543.25 | Fixed 3936 bank, sharper Toeplitz density |
| 4051--4067 | 1811.78 | Fixed 4080 bank, sharper Toeplitz density |
| 4068--4095 | 79.94 | Frozen box bank with boundary subdivision |
| 4096--4351 | 96.59 | Complete swapped-route box cover |
| 4352--5119 | 120.80 | Complete swapped-route box cover |
| 5120--8192 | 104.53 | Inherited unsigned bounds and terminal bound |

The 38 rows are pairwise disjoint and cover all 8,192 nonzero active-group
weights and all 11,465,830,401 feasible matching profiles. The reviewed
manifest combines the two dense rows into one checkpoint, giving 37 entries.
The final directed union has log2 upper bound
`-55092738508223/1000000000000` at both aggregation precisions.
Thus the complete construction has more than 55.09 bits of setup-failure
margin at 35% symbol distance. The
[closed certificate](HD_SHIFTED_CLOSED_CERTIFICATE.md) records the setup
distribution, pinned evidence, independent sum, and matched timings.

For example, replay the new 1280 band from the theory directory with

```text
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_1280_toeplitz_1256_1325.json --max-weights 70 --precision 256
python -B hd_shifted_potential_band.py --proposal hd_shifted_potential_band_1280_toeplitz_1256_1325.json --max-weights 70 --precision 384
```

The proposal selects the stronger density explicitly. Replays rebuild
all regional bounds; they do not accept saved FLOAT margins as evidence.
The 1280 directed union ceiling is `-710866278539/7812500000`, and all
70 individual ceilings agree. This interval contains 29,219,225 feasible
matching profiles. The remaining intervals now close through the other
proof paths in the table.

### The 3072 bank: directed receipt and FLOAT edge loss

The fixed 3072 bank closes all 458 integer weights in `2828..3285`,
including every matching profile and lane split. The proposal
`theory/hd_shifted_potential_band_3072_2828_3285_toeplitz.json` explicitly
selects `invertible-toeplitz`. It retains the same encoder and setup law.
Fresh 256- and 384-bit replays agree on every per-weight rational ceiling
and on the exact log2 union ceiling
`-113060545953299/1000000000000` (-113.060545953299 bits).
The band contains 1,072,732,409 feasible matching profiles.

The saved receipts and SHA-256 hashes are:

- [256-bit receipt](theory/hd_shifted_potential_band_3072_2828_3285_toeplitz_directed256.json):
  `c1ef1ef4c26b611c4a2a7e44c89b5e9eb2169d307cf5f565708013aa68ec7f1e`.
- [384-bit receipt](theory/hd_shifted_potential_band_3072_2828_3285_toeplitz_directed384.json):
  `35bb07e89732b0ebf429f5f42e58ca5299a220d05f729d8421d303c87af6ac62`.

The FLOAT screen is unusually loose at the lower edge. At `r=2828`, its
log2 upper bound is -125.74103779126781; the directed bound is
-593.414269637913. The 467.6732318466-bit difference occurs with the same
frozen witnesses. In contrast, the center and upper endpoint differ by
only about 0.00001557 and 0.00000196 bits. The lower screening edge is
therefore not evidence of a distance obstruction.

The next screening step should reuse logarithms of directed alias upper
tables before fitting another central bank. This comparison alone does
not certify any neighboring weight. Coverage remains limited to the
disjoint intervals above; the native code and measured performance are
unchanged.
