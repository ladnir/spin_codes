# A first column-late-alpha control

This control tests whether summing masks before fractionalization helps
at the first unresolved G16 anchor. It fixes

$$
r=18,\qquad a=8,\qquad z=691/700,\qquad\alpha=1/5.
$$

Only one contribution is retained: every outer column has zero moments,
and the boundary state is zero at every regional boundary.
The other outer syndromes and state paths are omitted. Consequently the
result is not a complete distance bound or a lower bound on failure.
A large value would obstruct this particular proposed majorant, not the
construction itself.

## Fresh unpowered regional input

Fix the active source groups and one eight-row set in each group. These
sets stay fixed across all sixteen outer columns.
The affine mixer is $H=I+\mathbf1\mathbf1^T+uu^T$ over $F=\mathrm{GF}(256)$,
with the G16 labels specified in `HD_TWO_SUM_GF16_ARITY_CERTIFICATE.md`.
On the zero-moment branch, $Hx=x$.
The independent output diagonal makes the remaining nonzero amplitudes
independent and uniform, conditioned on their entire support.

Let $R_1(m,f)$ be the unpowered two-state regional upper matrix for a
fixed support with $m$ nonzero entries and $f$ full source pairs.
Its route permutation, fair swaps, amplitudes, and inner refresh keys
have been averaged. The entry $R_1(m,f)_{00}$ bounds paths that enter and
leave that region in state zero, including intermediate activations.
The helper rebuilds this table for every $m\leq18$ and feasible $f$.
It does not take a fractional table to the power $1/\alpha$.

For a fixed slope, a regional occupancy envelope would introduce a
separate intercept. This control instead fixes that intercept to one:

$$
b:=\max_{1\leq m\leq18,\ 0\leq f\leq\lfloor m/2\rfloor}
     {\log_2 R_1(m,f)_{00}\over m},\qquad B:=2^b.
$$

Zero entries contribute $-\infty$ to the maximum.
The exact empty-input transfer is $R_1(0,0)_{00}=1$.
Thus $R_1(m,f)_{00}\leq B^m$ throughout the requested range.
This is a universal envelope over fixed geometries; it does not resample
the regional subsets or their full-pair counts.

## Two orders of the same positive sum

For a fixed support of size $w$, define the integer

$$
Z_2(w):=\sum_{j=0}^w(-1)^{w-j}\binom wj256^{\max(j-2,0)}.
$$

It counts nonzero assignments whose two Vandermonde moments vanish.
These exact integers are nonnegative; the subsequent evaluations use
only positive sums. There are $\binom aw$ labelled supports of size $w$.
Define

$$
P(X):=\sum_{w=0}^a\binom aw Z_2(w)B^wX^w,
\qquad C_k:=[X^k]P(X)^r.
$$

The late-column comparison uses

$$
U_k:=\left(255^{-ar/2}C_k\right)^\alpha.
$$

The sum over masks and groups is inside the power. Jensen first replaces
the conditional column mean by its route average; the regional envelope
then bounds that unpowered average. This is the upper direction of Jensen.
It loses route fractionalization but retains the column's mask sum.

For comparison, power each individual mask count before summing:

$$
V_k:=255^{-\alpha ar/2}[X^k]
 \left(\sum_{w=0}^a\binom aw Z_2(w)^\alpha B^{\alpha w}X^w\right)^r.
$$

The binomial mask multiplicity is not powered.
Both formulas use the same fresh $R_1$ envelope, not different packet
providers. Positive subadditivity gives $U_k\leq V_k$, with equality
when $\alpha=1$.

The normalizer $255^{-ar/2}$ occurs once per column. Across sixteen
columns it becomes the prescribed-support factor $255^{-8ar}$ before
fractionalization. It is not a new row-set count.

Let $\mathcal L$ sum coefficients of degree at least $9ar=1296$.
This aggregate support floor is necessary for actual nonzero RS$[16,8]$
rows. It relaxes the individual row floors and may admit wholly empty
rows or groups in the comparison.
The two reported terms are

$$
\binom{8192}{r}\binom{16}{a}^{r}
 \left({z^{-734003}\over255}\right)^\alpha
 \mathcal L\left[\left(\sum_k U_kX^k\right)^{16}\right]
$$

and the same expression with $V_k$.
Source and row-set counts stay outside $\alpha$.
The weight-only envelope makes the same bound valid for every fixed
eight-row set, so multiplying by $\binom{16}{8}^{18}$ is valid.
The maximal whole-word polynomial degree is 2304.

## Implementation and scope

`theory/hd_two_sum_late_alpha_quiet.py` uses positive log-domain
convolutions, never FFT coefficient extraction or signed FLOAT sums.
Its fixed scope permits one fresh regional bank, no witness fitting,
and a 30-second compute budget. FLOAT arithmetic remains diagnostic;
no directed certificate is asserted.

The seven small tests check the labelled GF4 kernel, positive convolution,
the two fractional orders, the split normalizer, the exact aggregate
tail, scope guards, and complete regional envelope coverage.
They passed in 0.006 seconds. An independent reviewer checked the source,
and the root agent separately ran these seven tests with the six
common-source partition tests. All thirteen passed.

### Recorded gate

The one authorized control completed in 0.783 seconds, including its fresh
regional bank. The selected slope was
$b=-7.354023590854679$, with contacts at occupancies zero and eighteen.

| Quantity | Late column alpha | Individual mask alpha |
|---|---:|---:|
| Normalized whole-word log2 term | -1963.2140821793632 | -352.8269273556288 |
| After the same global prefix | +1203.0350360598336 | +2813.422190883568 |

The common prefix is $3166.249118239197$ bits. Delaying mask
fractionalization reduces this term by $1610.3871548237344$ bits,
but this proposed majorant still fails by a large margin.

The artifact is `theory/hd_two_sum_late_alpha_quiet_gate_float.json`,
SHA-256 `1f9bcb13d4a5120a130e3b23dfad264ffe25a3fc42f9b5c2dc4e8488732d4521`.
The frozen helper hash is
`8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3`.
The artifact also records its two direct dependency hashes. No directed
replay or distance receipt was produced.

### What the saved data already rules out

Moving the remaining exponent past the sixteen-column degree sum cannot
remove the residual gap under the same mean comparison.
There are at most $145^{16}$ column-degree allocations. For nonnegative
terms $x_i$ and $N$ summands,

$$
\left(\sum_i x_i\right)^\alpha
\leq\sum_i x_i^\alpha
\leq N^{1-\alpha}\left(\sum_i x_i\right)^\alpha.
$$

Thus this last change can save at most
$(1-\alpha)\log_2(145^{16})=91.90283635219117$ bits.
The same majorant would still be at least $+1111.1321997076425$ bits.
The statement concerns the chosen upper expression, not actual failure.

The whole-row fiber Hölder bound from `HD_TWO_SUM_SYNDROME_RANK.md`
also does not repair this envelope. Substituting the saved fresh slope
into its exact nonzero RS enumerator gives a quiet contribution of
$+1722.43201961678$ bits. That expression is conditional on the rank-good
setup event; its separate exceptional-setup charge is not included here.
It retains actual nonzero RS rows and is nevertheless worse than this gate.

The regional line itself has substantial slack: 32.29747061859421 bits
at occupancy five and 23.598112468976595 bits at occupancy ten, before
the fractional power. These values justify retaining the saved regional
curve in the next diagnostic. They do not yet separate scalar-envelope
loss from the max-over-geometries relaxation or from unconditional Jensen.

The next comparison is a conditional column mean under one placement
partition common to every mask of the entire fixed source set.
A good-placement event can eliminate a nonempty packet's 0-to-0 step,
but not a region's activation-and-return paths. Its probability model
must account for each routed pair spanning two emission blocks.

## Retaining more of the regional curve

The next two controls use only the saved unpowered curve
$R_m:=\max_fR_1(m,f)_{00}$. They build no regional provider.
An information-set bound gives $Z_2(w)\leq255^{w-2}$ for $w>0$.
Put $\eta=255^{-2}$. For each fixed row set $A$,

$$
\sum_{S\subseteq A}Z_2(|S|)\prod_{i\in S}x_i
\leq_{\rm coeff}(1-\eta)+\eta\prod_{i\in A}(1+255x_i).
$$

The constant coefficient is still exactly one. Weights one and two are
now admitted synthetically; this is a positive upper comparison, not an
exact new quiet law. Expanding the product across groups selects $h$ groups
with coefficient $\binom rh(1-\eta)^{r-h}\eta^h$.
This auxiliary expansion does not change the fixed source or row sets.

For one common row set $A$ across all groups, define

$$
F_h(X):=\sum_{m=0}^h\binom hm255^mR_mX^m.
$$

The raw column polynomial is bounded coefficientwise by

$$
\sum_{h=0}^r\binom rh(1-\eta)^{r-h}\eta^h F_h(X)^8.
$$

Its coefficients are normalized and powered only after this entire sum.
The common-row-set scope has one $\binom{16}{8}$ choice, not eighteen
independent such choices.

### A moment envelope covers arbitrary fixed row sets

Suppose $R_m\leq\mathbb E[L^m]$ for a nonnegative random variable $L$,
with the order-zero moment defined as one. Use these moments in $F_n$.
Then $F_n(X)=\mathbb E[(1+255LX)^n]$.
For $j\geq k\geq1$, let $T=1+255LX$ and let $S$ be an independent copy.
The identity

$$
F_{j+1}F_{k-1}-F_jF_k
=\tfrac12\mathbb E[(TS)^{k-1}(T-S)
                         (T^{j-k+1}-S^{j-k+1})]
$$

has nonnegative coefficients: its integrand factors into
$(255X)^2(L-L')^2$ times a polynomial with nonnegative coefficients.

For arbitrary fixed row sets, selecting $h$ groups gives regional
incidences $0\leq n_i\leq h$ with $\sum_i n_i=8h$.
Repeated coefficientwise transfers from smaller to larger incidences give
$\prod_{i=1}^{16}F_{n_i}(X)\leq_{\rm coeff}F_h(X)^8$.
Thus the same displayed column bound covers every fixed row-set vector
when its regional curve is a positive moment envelope.
The union factor is then $\binom{16}{8}^{18}$, outside $\alpha$.

The sole tested envelope places probability $p$ at $L=\widetilde B$ and
the remaining probability at zero. Its slope is the secant through the
saved occupancies five and eighteen; then
$p=\max_{m>0}R_m/\widetilde B^m$ covers every saved occupancy.
The candidate is rejected if $p>1$.
Here $\log_2\widetilde B=-4.869602774039739$ and
$\log_2p=-44.71957470266891$, with contacts at zero, five, and eighteen.
No other slope was tested.

### Recorded saved-curve controls

Both controls completed together in 0.117 seconds.

| Scope | New log2 term | Matched old term | Improvement |
|---|---:|---:|---:|
| One common eight-row set; raw curve | +771.5426255627344 | +970.9557206969962 | 199.41309513426177 |
| All fixed eight-row sets; two-atom curve | +1026.1728393386634 | +1203.0350360598331 | 176.86219672116977 |

Neither closes the selected quiet contribution. These two numbers have
different row-set scopes and must not be compared as bounds for one class.
Six exact or mocked tests check the mixture, moment convexity, every tiny
fixed row-set vector, probability guards, and scope-specific counts.
They passed in 0.019 seconds; the root independently ran them with the
seven common-source analytic tests, and all thirteen passed.

The artifact `theory/hd_two_sum_late_alpha_moment_gate_float.json` has
SHA-256 `e7d7436ede16c1738808f78119b3d01a218a69d74ae0a0042aff739eb3fcab46`.
The helper `theory/hd_two_sum_late_alpha_moment.py` has hash
`fe9a8936e99f82a4aec31f903bc55ad6ebfec7fa8a67fdb30f6e37bc4f2ca8e6`.

An important remaining relaxation is the maximum over full-pair counts.
A source set with maximum even/odd full-pair count one permits only
$f\leq1$ in every actual regional subset. The earlier curve allowed up
to nine. The next isolated control restricts that geometry and charges
its exact source count once; it does not change the witness or fit a grid.

## Restricting the fixed source to one full pair

The restricted control uses $F_{\max}(G)=1$, where $F_{\max}$ is the maximum
of the even and odd matching full-pair counts. Every subset of this fixed
source has at most one full pair in either matching. Therefore its regional
curve is $R_m^{(1)}:=\max_{f\leq\min(1,\lfloor m/2\rfloor)}R_1(m,f)_{00}$.
The exact number of such source sets is

$$
N_1=154355220469939768130077702840442839490778956275433472.
$$

This count occurs once, outside the fractional power. All eighteen row
sets remain fixed, and the control still covers all of their arity-eight
assignments. The restriction does not condition separately on each mask.

One fresh unpowered regional bank gave the following FLOAT controls in
0.805 seconds. The four exact or mocked tests passed independently.

| Comparison | Log2 term |
|---|---:|
| Per-mask fractional power; restricted regional line | +2050.5455687524573 |
| Column-coefficient fractional power; same line | +440.84731200066653 |
| Column-coefficient fractional power; one two-atom envelope | +295.0135792564174 |

The line slope is $-10.270376890552134$. In particular,
$\log_2R_{18}^{(1)}=-184.8667840299384$, whereas the unrestricted curve gave
$-132.37242463538422$. The two-atom envelope uses the same prescribed
secant through occupancies five and eighteen. Its parameters are
$\log_2B=-8.907630419774677$ and $\log_2p=-24.529436473994224$.
It contacts the curve at zero, five, and eighteen.

The artifact `theory/hd_two_sum_late_alpha_source_f1_gate_float.json` has
SHA-256 `862d9f5dd08c304d70f282563bfb5130e5ce35423e7153d3b68721337838b045`.
The helper has hash
`77d3f7da4f7439766d0bb18b53cebcd0960697094c0e30361b064884b5ee1987`.
These positive terms still do not close the selected comparator component.
They are not lower bounds on a failure probability.

## Retaining the complete row-set and word sums

The preceding controls put each row-set assignment outside the fractional
power. That relaxation is unnecessary when using unconditional Jensen.
Fix $G$, and let $Q_G(\Pi)$ be the quiet/all-regional-zero positive component
of its conditional first-moment comparator. Here $\Pi$ denotes the entire
route, and $Q_G$ already sums all fixed arity-eight row-set assignments.
All other outer and inner setup randomness has already been averaged.

Splitting the full positive comparator into $Q_G$ and its remaining
components permits a charge $\mathbb E_\Pi[Q_G(\Pi)^\alpha]$.
Jensen bounds this charge by $(\mathbb E_\Pi Q_G(\Pi))^\alpha$.
This is a bound on one comparator component; it does not independently
characterize a quiet failure event.

Let $P(X)$ be the raw two-atom column polynomial constructed above for the
restricted curve. Its coefficient bound holds for every fixed row-set
assignment. Define the whole-word tail

$$
T:=\sum_{k\geq1296}[X^k]P(X)^{16}.
$$

Then the source union charges at most

$$
N_1\left[
\frac{z^{-734003}}{255}\binom{16}{8}^{18}
255^{-1152}T
\right]^\alpha
$$

to this component. In particular, the row-set sum and the complete degree
tail remain inside the power. The source count $N_1$ stays outside.
No independence between row-set assignments is asserted or required.

The direct first-moment endpoint $\alpha=1$ uses the same $P$ and the same
unpowered regional bank. It is a distinct valid charge for this positive
component, not a new fractional packet provider. More generally,
$\min(1,\sum_i x_i)\leq\sum_i\min(1,x_i)$ permits different valid scalar
charges for different positive components.

The isolated saved-curve helper `theory/hd_two_sum_late_alpha_word_order.py`
implements these orders without a provider rebuild or fitted parameter.
Six exact or mocked tests passed in 0.006 seconds. They check the complete
word tail, both inequality directions, the common $\alpha=1$ endpoint,
source counting, component splitting, and scope guards.

### Successful selected-component control

The saved-curve gate completed in 0.093 seconds. It reproduced the old
comparison and evaluated the complete positive tail without changing the
regional curve, construction, or weight marker.

| Ordering | Log2 term, including the F1 source count |
|---|---:|
| Column-coefficient power; row-set union outside $\alpha=1/5$ | +295.0135792564174 |
| Whole-word power; row-set union outside $\alpha=1/5$ | +251.529474111078 |
| Whole-word and row-set sums inside $\alpha=1/5$ | +54.94464227432175 |
| Direct first moment, $\alpha=1$ | **−432.03056185400425** |

The raw whole-word tail has logarithm $-5340.21773984005$.
After row-set counting, marking, and normalization, the fixed-source
comparison mean has logarithm $-608.7190051604075$.
Thus the direct first moment is stronger than its $1/5$ power.
This endpoint change does not raise the already averaged provider to a
different entrywise power: every displayed control uses the same $R_1$.

This result removes the numerical obstruction from the selected
quiet/all-regional-zero comparator component. It does not cover nonquiet
outer moments or paths with live regional boundary states. It is a FLOAT
control, not a directed certificate or a lower bound on failure.

The artifact `theory/hd_two_sum_late_alpha_word_order_gate_float.json` has
SHA-256 `a48e647831b07086fd0e6f7398759f111b718680f1af3c8f5943909b1b40670c`.
The helper has hash
`1e549ed2bc7f260c8a45029edbed1727acb50f59b66b114e3401f7d8de13b95e`.
An independent read-only audit confirmed the count and power chronology.
The root also reran fifteen new helper suites: all 96 tests passed in
1.854 seconds. The existing closed-certificate manifest is unchanged.

### First interface for the live-state remainder

Let $M_m$ be the same positive moment envelope used in this successful
control, with $M_0=1$. Retain all four entries of the unpowered restricted
regional matrix $R^{(1)}(m)$, and define the nonnegative matrix

$$
D_{ij}:=\max_{0\leq m\leq18}\frac{R^{(1)}(m)_{ij}}{M_m},
\qquad i,j\in\{0,1\}.
$$

Then $R^{(1)}(m)\leq M_mD$ entrywise. For any fixed mask array, positivity
bounds its chronological product of 256 regional matrices by
$\bigl(\prod_iM_{m_i}\bigr)D^{256}$. This argument does not commute or
reorder the original regional matrices.

The same scalar polynomial $P$ therefore gives a bound on all quiet
regional-state paths by multiplying the successful first-moment comparison
by $e_0^\mathsf TD^{256}\mathbf1$, where $e_0=(1,0)^\mathsf T$ and
$\mathbf1=(1,1)^\mathsf T$. The exact cost of this particular envelope is
$\log_2(e_0^\mathsf TD^{256}\mathbf1)$ bits. No new moment slope is needed.
This is a proposed next diagnostic, not an evaluated remainder bound.

The helper `theory/hd_two_sum_late_alpha_live.py` implements exactly this
comparison and saves all four restricted regional curves for later use.
Its six exact or mocked tests passed in 0.007 seconds. They include every
length-four word from a noncommuting matrix fixture, order-zero occupancy,
the terminal mass of both states, and a one-bank scope check.
The source hash is
`e28bd5473fe60be32a984e67cbd0b64f52078a30aa4f6ff4b5b304497fd4fb59`.

### The shared-matrix control is too loose

One fresh bank completed in 0.768 seconds and reproduced the frozen
zero-state curve. The shared matrix had entrywise logarithms

$$
\log_2 D=
\begin{pmatrix}
0 & 121.37341036333908\\
-6.986180905166236 & 44.1956849541163
\end{pmatrix}.
$$

All three nonzero-state ratios attain their maxima at occupancy eighteen.
The resulting state factor costs $14699.747295373461$ bits, giving
$+14267.716733519457$ for the quiet component with all state paths.
This particular envelope therefore fails.

The saved curves explain the mismatch. The 11 entry remains between
about $2^{-141}$ and $2^{-139}$, whereas the common moment envelope falls
to $2^{-184.8668}$ at occupancy eighteen. The ratio compensates by
inflating every live step, including steps with other occupancies.
The control does not show that the actual live-path sum is large.

All four curves are retained in
`theory/hd_two_sum_late_alpha_live_gate_float.json`, whose SHA-256 is
`33da1eb9517c7d5665d22d2c8b93f5622199d4e1b70dd2a8e3aa52e28d4e04cf`.
An independent review confirmed the matrix order, all-occupancy maxima,
and terminal reference mass. No nonquiet outer contribution was included.

### A distinct-entry positive matrix measure

A next comparison can keep the distinct regional decay rates without
altering the saved bank. Write $p_0,B_0$ for the frozen 00 moment law.
Choose one cross-state slope from the 01 curve's secant at occupancies
two and eighteen. Call its positive atom $B_\times$, and define

$$
a_{01}=\max_m R^{(1)}(m)_{01}/B_\times^m,\quad
a_{10}=\max_m R^{(1)}(m)_{10}/B_\times^m,\quad
a_{11}=\max_m R^{(1)}(m)_{11}.
$$

Every maximum includes occupancies zero through eighteen. The resulting
entrywise majorant is

$$
\widehat R(m)=
\begin{pmatrix}
(1-p_0)\mathbf1_{m=0}+p_0B_0^m & a_{01}B_\times^m\\
a_{10}B_\times^m & a_{11}
\end{pmatrix}.
$$

It is a positive matrix-valued atomic measure on $\{0,B_0,B_\times,1\}$,
using the convention $0^0=1$. Its atoms have restricted matrix supports:
the first two occur only on 00, the third only off-diagonal, and the last
only on 11. This avoids enumerating all four-atom words independently.
The off-diagonal entries at occupancy zero are still relaxed, explicitly.

For a chronological sixteen-step binary state path, let $n_{ij}$ count
its transitions. Expanding its 00 factors selects $k$ active atoms with
weight $\binom{n_{00}}k p_0^k(1-p_0)^{n_{00}-k}$.
The nonzero scalar markers are $k$ copies of $B_0$,
$n_{01}+n_{10}$ copies of $B_\times$, and $n_{11}$ copies of one.
The remaining $n_{00}-k$ zero markers force every group mask to avoid
those positions. State paths sharing transition counts have the same
scalar contraction; their exact chronological multiplicity remains.

For the full homogeneous arity-eight class, let
$W_w=(\sum_{|A|=8}N_A(T)^{16})^{1/16}$ for any fixed output support $T$
of size $w$. These raw Hölder weights include row-set choices but exclude
the RS denominator. The closed kernel makes them depend only on $w$.
Every quiet output has weight at most eight and equals its pre-mixer
input. Every nonquiet output has weight at least nine, and its input
weight is at most eight. Thus $\kappa(w)=\min(w,8)$ is a valid upper
surrogate for input weight, independently of $A$.

If $b$ is the displayed list of nonzero markers and $e_w(b)$ its elementary
symmetric polynomial, the scalar contraction for that path and $k$ is

$$
\left(\sum_{w=0}^{16}W_w e_w(b)X^{\kappa(w)}\right)^{18}.
$$

Weight these polynomials by the 00 expansion, the factors
$a_{01}^{n_{01}}a_{10}^{n_{10}}a_{11}^{n_{11}}$, and the exact path counts.
Their endpoint sums define a $2\times2$ positive polynomial $Q(X)$.
The proposed complete-class charge is

$$
\frac{N_1z^{-734003}}{255^{1153}}
\sum_{k\geq1296}[X^k]e_0^\mathsf TQ(X)^{16}\mathbf1.
$$

There is no further row-set factor: it is already inside $W_w$.
Matrix multiplication across columns preserves their chronological
boundary states. The row sets remain fixed through the sixteen-column
Hölder step.

### Recorded ordered-matrix controls

The four-atom helper `theory/hd_two_sum_matrix_moment.py` implements this
interface. Seven small tests passed in 0.018 seconds, including exhaustive
labelled masks with noncommuting matrices and scalar/polynomial equality.
Its source hash is
`bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67`.
An independent audit confirmed the fixed-row-set, path, and degree order.

The complete homogeneous arity-eight/F1 control finished in 0.978 seconds.
Its fixed cross-state atom has logarithm $-3.590815685409988$.
The three scale logarithms are $1.1413086707804538$ for 01,
$-127.21828259772485$ for 10, and $-138.88133772969238$ for 11.
Every occupancy is checked against these envelopes.

| Control | Log2 upper comparison |
|---|---:|
| One scalar Chernoff marker derived from the earlier quiet word | +6467.669977230867 |
| Full positive matrix-polynomial tail | +3497.359683472526 |
| Quiet-only masks; same full matrix-polynomial tail | +3497.359683472526 |

The scalar marker has $\log_2X=4.091925496837575$.
It is fixed across every path and column; no neighboring marker was tried.
The polynomial mode needs no input-weight marker and retains the entire
degree-1296 tail of the degree-2304 word polynomial.
It sums 274 path classes and 1650 positive path/atom terms, with 834 cached
group polynomials. Off-diagonal transitions at zero occupancy are still
admitted by this envelope.

The full artifact `theory/hd_two_sum_matrix_moment_gate_float.json` has
SHA-256 `0465d355fd56c5c0cdfabac206a377ed74f62bf638b78af354bb4a147a9b5783`.
The quiet-only artifact `theory/hd_two_sum_matrix_moment_quiet_float.json`
has hash `9263a38cdf5311a8f5ddf10490c8f9472eaecee46bf8629a86f94f65ec53e823`.
Equality of the reported FLOAT totals is not an exact equality theorem.
It indicates that nonquiet masks are numerically negligible within this
particular comparison; the live-state quiet terms remain the obstruction.

The next saved-bank diagnostic decomposes the positive tail by degree,
column state word, and regional path class. It does not identify an actual
bad-codeword event or prove a lower bound on failure. A separate positive
coverage refinement can enforce that every 01 or 10 position is nonempty,
without subtracting independently bounded quantities.

### What dominates the ordered comparison

The saved-bank decomposition completed in 3.133 seconds and reproduced the
positive tail. Its peak degree is exactly 1296. One largest column word
has nine degree-144 columns and seven degree-zero columns, ending in state
one. That single column word accounts for a $2^{-12.930670458649274}$
fraction of the peak-degree coefficient. This is a comparator term, not
an actual failure event.

At column degree 144, the dominant 00 term has transition counts
$(n_{00},n_{01},n_{10},n_{11})=(8,4,4,0)$ and no active 00 atoms.
Its 495 chronological paths account for about 99.9988 percent of the
coefficient. At degree 81, the dominant counts are $(10,3,3,0)$, again
with no active 00 atoms. Thus activation-return excursions dominate this
envelope; sustained 11 paths do not.

The degree-144 term requires all eighteen groups to occupy the eight
cross-state regions. Its cross-state positions are already nonempty.
Consequently, enforcing the exact off-diagonal empty zeros cannot remove
this dominant term. A nine-column minimum-row-weight pattern is compatible
with the admitted aggregate floor; the comparison does not resample a
column support to infer that pattern.

The diagnostic helper has SHA-256
`57a4a9fa4c717a5836233c944f1c6b60fe8c21ca12f5a65f57461bca460d22f7`.
Its artifact `theory/hd_two_sum_matrix_moment_diagnostic_float.json` has
hash `b3db534a024f03a0fd390dd7c737a85a6fc92d5110f519edcce934cf67bb19c7`.
Three small tests check max-plus backtracking against every toy word,
regional decomposition against the column coefficient, and input guards.

### A degreewise outer-only cap

For quiet masks, total output-support degree $K$ equals pre-mixer input
degree. Let $C_K$ upper-bound the unweighted outer first moment at degree
$K$, and let $S_K$ upper-bound its inner-weighted first moment. Both sums
include all fixed row sets in the same scope. The bad-word contribution
at that degree is at most

$$
\min\{C_K,\;z^{-734003}S_K\}.
$$

This minimum is taken after each complete degree sum; it does not
subtract upper bounds or cap isolated state-path terms. The outer-only
comparison does not need the marked calculation's Holder relaxation.
After forgetting output positions, every fixed eight-row set has the same
column polynomial. Count its choices once per group, and take

$$
C_K=255^{-1152}\binom{16}{8}^{18}[X^K]
\left(\sum_{w=0}^8\binom{8}{w}Z_2(w)X^w\right)^{288}.
$$

The source count and projective factor are then $N_1/255$, once outside
the degree sum. The isolated helper `theory/hd_two_sum_quiet_degree_cap.py`
implements this comparison without new providers or parameters. Its
four small tests passed at the root in 0.002 seconds, including an explicit
sum over the same row set in two columns. The revised helper hash is
`4e2e7eaeeaae72921798e2badd480ae491bfca584ecde3d18553b8df2d02740c`.
The completed FLOAT control lowers the comparison from +3497.359683472526
to **+2943.1151222813605**, but still does not close it. Its dominant degree
moves from 1296 to 1830. The saved artifact
`theory/hd_two_sum_quiet_degree_cap_float.json` has SHA-256
`2c223d522278a37b6a6ffbadd7db98d2e48257f80702b70e60e2a55bc60163b1`.
This cap excludes nonquiet masks and is not a whole-code certificate.

The separate [two-stream proposal](HD_TWO_STREAM_INNER_PROPOSAL.md)
describes a construction-level alternative for the next cycle. It is
proof-first and unimplemented; no timing or certificate is transferred.
