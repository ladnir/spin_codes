# A Renyi comparison for the scheduled lambda outer

The uniform comparison leaves a gap at eighteen groups of arity twelve.
The missing information is not a scalar codimension discount: a linear
subspace and the ambient space both contain the zero vector. This note
retains that exceptional density through an exact third moment, before
applying the existing fractional uniform-amplitude bound.

This is a proof of the comparison interface, not a numerical certificate.
The first proposed scope is `r=18`, every active group having exactly
twelve nonzero RS rows, all original row-set choices, and all source
pairing profiles. No new numerical result is asserted here.

## Construction and the comparison distribution

Set `q=256`. The outer has sixteen RS[16,8] rows per group and sixteen
columns. Its column map before the final diagonal is

\[
M_j=Q_jH_\lambda P_j,\qquad
H_\lambda=I+\lambda\mathbf1\mathbf1^T+
\sum_{k=0}^{3}b_kb_k^T,
\]

where `lambda=0xbc` in the AES field and `(b_k)_i` is bit `k` of label `i`.
Source row `i` lands at
`P_j(i)=inverse_GF16(i xor j)`, with inverse zero defined as zero and
GF16 polynomial `0x13`. Mixer output `i` lands at `Q_j(i)=i xor j`.
There is no crossed-pair layer. The original one-state inner and the
paired route with fair swaps are unchanged.

The input and output diagonal entries are independent uniform nonzero
field elements sampled once at setup. The checked native normalization
has this ideal ensemble's distribution; this note does not replace that
construction-binding argument.

Fix an original row set `A`, of size `a`, across all sixteen columns.
Define the column subspace

\[
C_{A,j}=M_j\mathbb F_q^A\subseteq\mathbb F_q^{16}.
\]

Let `P_(A,j)` be the law obtained by sampling uniformly from this subspace
and applying a fresh independent nonzero diagonal. This is an auxiliary
comparison law, not the law of an actual RS column conditioned only on
`A`. Let `U_16` be uniform on the ambient field space, and let
`L_(A,j)=dP_(A,j)/dU_16`.

Write `n_(A,j)(S)` for the number of subspace vectors with exact support
`S`. Diagonal averaging makes amplitudes uniform conditional on support.
Consequently, for a target `y` with support `S`,

\[
P_{A,j}(y)=\frac{n_{A,j}(S)}{q^a(q-1)^{|S|}},
\qquad
L_{A,j}(y)=q^{16-a}\frac{n_{A,j}(S)}{(q-1)^{|S|}}.
\tag{1}
\]

In particular, `sum_y P_(A,j)(y)=1`. The likelihood at zero is
`q^(16-a)`, so it cannot be replaced by a constant close to one.

## Exact third moments from the complement

The matrix `M_j` is orthogonal: `M_j^T M_j=I`. Its permutations are
orthogonal, and `H_lambda` is symmetric and involutive. Hence a parity
check for `C_(A,j)` is

\[
B_{A,j}=M_j[:,A^c]^T.
\]

For every coordinate subset `T`, the number of code vectors supported
inside `T` is

\[
F_{A,j}(T)=q^{|T|-\operatorname{rank}(B_{A,j}[:,T])}.
\]

Integer subset Mobius inversion gives

\[
n_{A,j}(S)=\sum_{T\subseteq S}(-1)^{|S|-|T|}F_{A,j}(T).
\tag{2}
\]

These are nonnegative integers, with `n(empty)=1` and total mass `q^a`.
For arity twelve, each rank calculation has only four parity-check rows.
The entries lie in the embedded GF4, so rank over GF4 equals rank over
GF256. No enumeration of `q^12` inputs is required.

For general `b=16-a`, equation (1) yields the exact rational moment

\[
I_3(A,j):=\mathbb E_{U_{16}}L_{A,j}^3
=q^{3b-16}\sum_{S\subseteq[16]}
\frac{n_{A,j}(S)^3}{(q-1)^{2|S|}}.
\tag{3}
\]

Its value is at least one. At arity twelve the leading factor is `q^-4`.
The support counts are cubed before summation; neither a support count
nor its multiplicity is treated as a probability without normalization.

## RS counting domination

For one prescribed nonzero RS row with support `T`, shortening and an
information set give at most `(q-1)^(|T|-8)` full-support codewords.
Independent input diagonal entries assign a fixed nonzero amplitude word
mass `(q-1)^(-|T|)`. Their product is at most `(q-1)^-8`.

For `a` active rows, the expected counting measure before the column
maps is therefore at most `(q-1)^(-8a)` times counting measure on all
`16a` allowed input coordinates. This comparison drops the individual
RS support floors and permits zero rows in the enlarged space.
It does not assert independence of the actual RS columns.

Push this counting domination through the fixed column maps and average
the independent output diagonals. Define

\[
\beta_a=\frac{q^{16a}}{(q-1)^{8a}},\qquad
d_A=\left(\prod_{j=0}^{15}I_3(A,j)\right)^{1/3}.
\tag{4}
\]

The resulting word counting measure is bounded by `beta_a` times the
product comparison law over columns. This probability law has likelihood
norm exactly `d_A` relative to the uniform word law. The row set `A`
remains fixed throughout this product.

For an ordered tuple of row sets `(A_1,...,A_r)`, independent outer keys
give the product of these comparison laws. Its counting factor is
`product_g beta_(|A_g|)` and its likelihood norm is `product_g d_(A_g)`.

## Indicator Holder bound and the order of averaging

Fix an indexed source-group set `G` and its complete route realization
`Pi`, including swaps. Let `K` denote the remaining inner setup, which is
independent of the outer keys. Let `U_G` be uniform on all `256r`
potential field coordinates of these groups, including zero amplitudes.
For a fixed output-weight limit `W`, define

\[
p_\Pi=\Pr_{Y\gets U_G,K}
 [\operatorname{wt}(\operatorname{Inner}_{\Pi,K}(Y))\le W],
\qquad
F_\Pi(Z)=\mathbb E_{Y\gets U_G,K}
 Z^{\operatorname{wt}(\operatorname{Inner}_{\Pi,K}(Y))}.
\]

For each fixed row-set tuple, apply Holder with exponents `3` and `3/2`
under the joint uniform-amplitude and inner-key measure. The indicator's
`3/2` power is itself. Thus its expected bad-word count is at most

\[
\left(\prod_g\beta_{|A_g|}d_{A_g}\right)p_\Pi^{2/3}.
\tag{5}
\]

All tuples use the same `p_Pi`, so they may be summed before any
fractional power. For homogeneous arity twelve, put

\[
B_{12}=\beta_{12}\sum_{|A|=12}d_A.
\tag{6}
\]

There are 1,820 sets in this sum. The sum over ordered group assignments
in equation (5) is exactly `B_12^r`; it already includes all row-set and
cohort multiplicities.

Every nonzero bad scalar message has `q-1` projective multiples of the
same weight and source class. Markov's inequality for this actual count,
followed by the enlarged counting bound, gives

\[
\Pr[\text{bad message with source }G\mid\Pi]
\le \min\left(1,\frac{B_{12}^r}{q-1}p_\Pi^{2/3}\right).
\]

Use `min(1,x)<=x^theta` with `theta=3/10`. For a rational mark `0<Z<1`,
the ordinary output-weight Chernoff bound gives
`p_Pi<=Z^-W F_Pi(Z)`. Therefore

\[
\Pr[\text{bad message with source }G\mid\Pi]
\le
\left(\frac{B_{12}^r}{q-1}\right)^{3/10}
Z^{-W/5}F_\Pi(Z)^{1/5}.
\tag{7}
\]

This indicator proof uses no irrational mark. It also avoids moving a
fractional power through the outer support sum. Setup reuse across
messages causes no issue: the argument uses counting expectations, not
independence between messages.

## The existing uniform provider is unchanged

Equation (7) requires exactly the old uniform-amplitude functional at
mark `Z` and fractional exponent `alpha=1/5`. Let `R_(r,e)(Z,alpha)` be
the directed regional matrix from `hd_shifted_swap_sparse.py`. Let
`N_r(e,o)` be its reviewed integer count for indexed source sets with
`e` even-matching pairs and `o` odd-matching pairs. The resulting bound is

\[
\left(\frac{B_{12}^r}{255}\right)^{3/10}Z^{-W/5}
\sum_{e,o}N_r(e,o)\,
e_0^T\bigl(R_{r,e}R_{r,o}\bigr)^{128}\mathbf1.
\tag{8}
\]

All source profiles are unioned outside `theta`; for `r=18` there are
99 indexed profiles. The regional packet caps, their prescribed
composition before fractionalization, and the unpowered geometry
normalizers are unchanged. No codimension factor is multiplied onto a
finished regional matrix.

The new logarithmic prefix is

\[
\frac3{10}\bigl(r\log_2B_{12}-\log_2 255\bigr)
-\frac W5\log_2Z.
\tag{9}
\]

There is one projective factor, with exponent `3/10`, not `1/5`.
The sum in equation (6) must not receive another factor `1820`.
Fresh directed moments and the reviewed provider are still required
before equation (8) becomes a numerical receipt.

## Inclusion monotonicity and further mixing

For a linear code `C`, let `D=C^perp` and define
`psi(0)=1`, `psi(y)=-1/(q-1)` for nonzero `y`. Additive character
orthogonality gives the diagonally averaged likelihood

\[
L_C(y)=\sum_{v\in D}\prod_{i\in\operatorname{supp}v}\psi(y_i).
\]

Under uniform field symbols, the moments of `psi` through order three are

\[
m_0=1,\quad m_1=0,\quad m_2=\frac1{q-1},\quad
m_3=\frac{q-2}{(q-1)^2}.
\]

Expand exactly three copies of the likelihood. Each ordered dual triple
has nonnegative expectation, as a product of these moments. Consequently,

\[
C\subseteq C'\quad\Longrightarrow\quad I_3(C')\le I_3(C).
\tag{10}
\]

This is a property of the exact third moment, not a termwise fractional
power argument. A maximum over every arity-twelve mapped set therefore
bounds every larger-arity column: choose a twelve-element subset of its
input set. Output permutations do not change the moment under the
uniform reference law. Transferring this moment bound does not guarantee
that the final larger-arity bound closes; its message-counting factor
also grows.

Applying a fixed invertible map to a given comparison distribution
preserves its Renyi norm, since it permutes the uniform reference space.
An independent random invertible map contracts that norm by Jensen's
inequality. A fresh nonzero diagonal is one such averaging operation.
Immediately appending another diagonal gives no strict improvement to
an already diagonally averaged law. A fixed mixer followed by fresh
diagonal averaging can improve it. Changing a mixer before the existing
diagonal can also change the norm; fixed-map invariance alone does not
compare those differently averaged constructions.

## A finite mixed-arity next step

The first gate should evaluate equations (3)--(9) only for arity twelve.
It must include all original row sets and retain their actual scheduled
column products, or upper-bound each such product by a common maximum.
Neither method replaces the schedule by independent row sets in each
column.

For a later mixed class, use the same Renyi order and uniform functional,
and define `B_a=beta_a sum_(|A|=a) d_A`. If allowed tuples have total
arity in a set `T` and individual arities in `J`, replace `B_12^r` by

\[
\mathcal B_{J,T}
=\sum_{k\in T}[u^k]\left(\sum_{a\in J}B_a u^a\right)^r.
\tag{11}
\]

Apply the `3/10` power to this complete positive sum. The polynomial
already counts the ordered placements of the arities among groups.
Bounds built with different Renyi orders cannot be inserted into this
same polynomial without a new argument.

A small next partition is `J={12,...,16}`, with total arity below the
separately certified high-side cutoff. Equation (10) supplies conservative
moment bounds for the larger arities without a new support enumeration.
Tuples containing a group of arity at most eleven remain outside this
partition, even when their total arity is at least216. They need a
separate bound. This note makes no whole-code closure claim.

## Observed directed dense-group cover

The comparison has now been implemented without changing the encoder.
The exact full-labelled third-moment inventory has maximum
I3 approximately1.0078591195655953. Using that common maximum first
gave69.240061061278bits for homogeneous arity twelve at r18.

The next calculation uses both third and fourth moments and the actual
scheduled column products. Forty-eight proved orthogonal label maps
preserve the full H_lambda matrix, including its diagonal. They reduce
2,517 labelled complements of sizes0..4 to134 support calculations.
Each original complement and an explicit coordinate-permutation witness
remain in the inventory. All1,820 third moments for arity twelve agree
exactly with the earlier independent full-labelled enumeration.
The two-moment inventory took13.609675seconds locally.

For each common order p in{3,4}, define

    w_a = sum_{|A|=a} (product_j I_p(P_j(A^c)))^(1/p),
    beta = 256^16 / 255^8.

At total arity k, the comparison coefficient is

    beta^k [u^(k-216)] (sum_{a=12}^{16} w_a u^(a-12))^18.

Use theta=(1/5)*p/(p-1), retaining the existing provider's exponent1/5.
The original A is fixed throughout each sixteen-column product. The
polynomial already counts the ordered assignments of arities to groups.
Different orders are alternative bounds on the same class; the smaller
completed bound may be selected. Total arities partition the messages,
but their bad-setup events can overlap, so the final sum is a union bound.

With Z=691/700 and the unchanged35% output-weight threshold, the
fourth-moment comparison gives the following directed bounds:

| Selected r18 message class | Margin bits |
| --- | ---: |
| Every group has exactly12 active rows | 512.832514143170 |
| Every group has12..16 active rows; total arity216..245 | 40.822846778659 |
| Any individual arities; total arity246..288, using the prior comparison | 51.151086488213 |
| Union of the preceding two ranges | 40.821725029550 |

The union covers every r18 source with all group arities at least12,
and also every r18 source of total arity at least246. Every bound includes
all99 indexed pairing profiles. Fresh256/384-bit replays agree on all62
method ceilings and31 selected cumulative ceilings. The final union is
recomputed at512bits. Sixty focused regression tests pass, and independent
review checked the moment calculation, schedule, coefficients, and scope.

The remaining r18 class has total arity at most245 and a group of arity
at most11. Other r remain uncovered by this new certificate. The checked
native latency remains3.442226ms; no mixing layer or encoding work was
added. This is a partial cover, not a whole-code certificate.

Evidence in `theory/`:

- `hd_lambda_renyi_orbit_inventory34.json`: 227fd5c05f5946698d6fac8a703abe4657ec154f817f48352d0ce346f818e673.
- `hd_lambda_renyi_mixed_256.json`: 534de633bc59f905463e5260253a0b947d7c4ca70747539692a3d32b78579401.
- `hd_lambda_renyi_mixed_384.json`: bf5f2d5dcfd0d12393350779dab811b1054783f57d9f6447bde8feec21fdac74.
- `hd_lambda_renyi_dense_cover_verified.json`: ca42f73dcefd3fc252bec47e75f6b584e182385dfa7c2185675ad8b3929092ec.

Replay the exact orbit inventory with `hd_lambda_renyi_orbits.py`, then
run `hd_lambda_renyi_mixed.py` separately at256 and384bits. The frozen
`hd_lambda_renyi_cover_verify.py` checks the observed input hashes and
cross-replay agreement, then sums the mixed and high-side bounds. Hashes
bind those observed runs; they are not replacements for the proof.

## Retaining outer activity and original input weight

The next comparison retains a constraint that the preceding uniform
extension discarded. Every active RS[16,8] row has at least nine nonzero
symbols. Consequently, every active group has at least nine nonzero
columns. The invertible column map and nonzero diagonals preserve whether
a whole column is zero.

Fix an original row set A of size a and an integer p>=2. Under the uniform
reference on one column, the zero contribution to its normalized density
moment is exactly

    z_(a,p) = 256^(p*(16-a)-16).

For scheduled column j, put n_j=I_p(P_j(A))-z_(a,p). Subtract these exact
rational numbers before numerical evaluation. The density restricted to
at least nine nonzero columns has p-th moment

    T_(A,p) = sum_{k=9}^{16} [v^k] product_j (z_(a,p)+n_j*v).

The new per-group coefficient is

    B_a(p) = 256^(16a)/255^(8a) * sum_{|A|=a} T_(A,p)^(1/p).

This is a restricted, unnormalized density. There is no conditional-law
normalizer. Its norm can be below one. Every original A persists across
all sixteen columns; sorting scalar moment signatures only reuses the
commutative polynomial product.

The exact arity-eleven inventory covers all4,368 labelled complements,
using166 representatives under the same proved48-map subgroup. Its
packed binary rank engine retains the full diagonal. Independent review
checked the rank shortcuts and the coordinate-map witnesses.

Fresh256/384-bit activity replays agree on all144 method ceilings and48
selected cumulative ceilings. For r18, with every group arity in11..16,
the union over total arities198..225 has55.461496706851 margin bits.
All99 indexed source profiles are included. Total arity226 gives only
39.713395923127bits individually, and larger totals remain open under
this extension. The previously closed dense-group classes remain valid.

For3<=p<=4, Holder on the nonzero event gives

    n_p <= n_3^(4-p) * n_4^(p-3).

Applying this interpolation before the activity tail improves some lower
totals, but the tested eighth-step grid does not extend the cutoff225.
Both precision replays agree. Interpolating only the completed class
bound would not give the same comparison.

For arities1..5, a second helper retains original input weight. Let
n_A(w,S) count allowed input vectors of weight w whose mixed output has
exact support S, and define K_(A,S)(X)=sum_w n_A(w,S)X^w. The raw column
moment polynomial is

    J_(A,p)(X) = 256^(16(p-1)) * sum_S K_(A,S)(X)^p/255^((p-1)|S|).

The exact inventory retains all6,884 labelled sets and all cross terms
between input weights for the same S. Its normalization check is
J_(A,p)(1)=256^(a*p)I_p(A). In the moment expansion, the p possibly
different RS-valid terms each have total input weight at least9a.
Define a separate upper comparison coefficient by

    B_a^marked(p) := 255^(-8a) * sum_{|A|=a}
      (sum_{d>=9ap} [X^d] product_j J_(P_j(A),p)(X))^(1/p).

Both this coefficient and the earlier activity-only coefficient bound
the required RS-density norm; no ordering between them is asserted.
There is no additional256 factor. Reversing the degree and truncating nonnegative convolution
computes this tail without subtracting large nearly equal totals.

The complete marked inventory took6.753077seconds. Directed word-tail
replays took11.566307 and10.859433seconds. They certify homogeneous
arity one at r18, across all99 profiles, with975.628296489068bits.
They do not close the whole sparse class: the all1..5 comparison has a
positive log2 bound of1133.362271844945. Arity two dominates that bound.
This is a failure of this uniform-reference comparison, not a discovered
low-weight word or a contradiction of the earlier sparse-profile bounds.

## What another column-local mixing stage can and cannot fix

Every invertible column-local mixer preserves the auxiliary zero-column
probability256^(-a). Jensen's inequality gives the minimum possible
nonzero moment

    n_ideal = (1-256^(-a))^p / (1-256^(-16))^(p-1).

Insert this value and the unchanged zero contribution into the positive
activity tail. The resulting expression is a lower comparator for this
fixed proof method, not a failure-probability bound. At the current mark,
alpha1/5, and p in{2,3,4}, its all-arity r18 union through total245 is
still at least2^380.099527319243. Both precision replays agree.
Thus arbitrarily strong column-local mixing cannot close this particular
activity-only comparison. This statement does not exclude another mark,
moment order, reference distribution, or stronger RS constraints.

Two distinct next steps remain useful. A small deterministic transvection
can reduce the actual nonzero-moment spikes without another field
multiplier. Separately, a reference distribution with an explicit
whole-column zero branch can avoid charging those zero columns as though
they were sampled uniformly from256^16 vectors. Both require new checks;
neither is a completed certificate or a measured speed improvement.

Current receipts in `theory/`:

- `hd_lambda_renyi_a11_inventory234.json`: f08e68c79dea8599f2a63c8afef67dbe10eede78b56a4b3d33f0be03bc1842ba.
- `hd_lambda_renyi_activity_a11_256.json`: 77edcd1ea1f70b3557e249198912a267721b37dca82a3bba4d8abac1fe841464.
- `hd_lambda_renyi_activity_a11_384.json`: AF2B5FB312136BB37EBABC41C9523122CB8A1C8D30D8D7E8773FDDBFEE97D766.
- `hd_lambda_renyi_marked_inventory234.json`: db3f7d554cf9e6168d29d660e476610c2d0934a359a4a7edd802bbf89873ecfe.
- `hd_lambda_renyi_marked_control_256.json`: 5f2dd8a84852d7760f62a57b14479aa1e0b8457364cd4fdd5da17c1b09c97ad2.
- `hd_lambda_renyi_marked_control_384.json`: ceb3c9ec647280e997cdd16aa2015c3b63bec9bdbd79c11733a9cb7236bea563.
- `hd_lambda_renyi_limits_256.json`: f3d433f3e55d369741f93e53370d57e1423ab5b12ff033328c168e817ddf71be.
- `hd_lambda_renyi_limits_384.json`: 0a425a1acb6466880226c08abcd131ab910b4f76f3f9ed88a939a2b99644f31e.

The measured encoder and its3.442226ms timing are unchanged. These
results expand and diagnose partial proof coverage; they do not close
the whole-code goal.

The frozen activity verifier recomputes the union of the prior dense
cover, the new total198..225 class, and homogeneous arity one at512bits.
The union has40.821668515604 margin bits. Overlapping message classes
are harmless here because the calculation uses a probability union
bound, not an assertion that bad-setup events are disjoint. The receipt
is `theory/hd_lambda_renyi_activity_cover_verified.json`.

## Completed extra-stage and zero-reference controls

The first construction test adds a transvection supported on coordinates
V={0,1,2,3}. With v its binary indicator, the new mixer is
H'=H_lambda+vv^T=(I+vv^T)H_lambda. Here H_lambda v=v and v^T v=0,
so H' remains a symmetric involution. The parity of these four inputs is
already available in the circuit: applying the change needs four further
XORs per column, not another field product. This is an operation count,
not a measured full-encoder cost.

The old symmetry quotient does not apply unchanged. Its stabilizer of V
has four maps. A complete exact inventory therefore uses2,089 tables to
cover all6,885 labelled row sets of arities11..16. The inventory took
206.47seconds. All120 two-source profiles retain minimum output weight
six. Some eleven-row moments improve, but a new complement
{0,7,11,13,14} has log2(I4) approximately40.00565; the old maximum was
approximately33.0028. The full scheduled comparison is also worse:

| Tested class at r=18, all99 profiles | Original H_lambda | H' |
|---|---:|---:|
| Total active-row count225 | -55.47 | -50.94 |
| Total active-row count226 | -39.71 | -35.19 |
| Union of totals198..245 | +163.242415 | +166.260042 |

Entries are log2 upper bounds, not measured failure probabilities. Every
group has arity11..16. The old order-two control uses an order-three
upper bound at arities12..16; the new control uses exact second moments.
Thus this is a comparison of completed bounds, not an isolated estimate
of the matrix change. All144 method ceilings and48 selected ceilings
agree in fresh256/384-bit replays. No sparse certificate or timing is
transferred to H'.

The second test leaves the construction unchanged and changes the proof
reference to Q_eta=eta*delta_0+(1-eta)*U16. Here U16 is uniform on
GF(256)^16. A single latent flag selects the zero or uniform branch for
an entire column; the sixteen coordinate flags are not independent.
The density norm still retains at least nine nonzero columns per active
group. Its reference calculation also retains at least9r active flags.
It uses a positive coefficient marker and a common regional-matrix
majorant over the source profiles, not conditioning or renormalization.

The tested eta values are0 and2^-e for e in{96,104,112,120,128,136},
with moment orders3,4. Every selected total198..245 chooses eta=0.
The combined bound remains+163.242415. Both precision replays agree on
all672 method ceilings. This rejects this tested reference comparison,
not all possible zero-aware arguments. In particular, the common matrix
majorant sacrifices some information about regional source profiles.

The next test screens two or three independent cheap transvections
together. A large proportional-column class in a shortened parity matrix
provides a cheap rejection test before exact moment enumeration. A second
direction changes the local wiring across columns, since no invertible
column-local mixer can activate an entirely zero column. Neither direction
has a certificate or a new native timing at this point.

Reproduction uses the new controllers in `theory/`. They refuse to
overwrite existing evidence. The observed receipts are:

- `hd_lambda_extra_transvection_inventory234.json`: 2089c87acfb442319465a8e888cb06d5091ca5ebedde26d067dadebe5b8390e5.
- `hd_lambda_extra_transvection_control_256.json`: 4b04fa117fe7896d6ae2d2462622ad517a095fd265bec8c00477e0330865daee.
- `hd_lambda_extra_transvection_control_384.json`: 93ce9134285d5775e1cc187b52189307a91fd05744baaa14acb8eba30133e72d.
- `hd_lambda_renyi_zero_reference_final_256.json`: ae289c45b1e51159425e5645001c60c81b71ed3a1cac4a833ccc65baea033d78.
- `hd_lambda_renyi_zero_reference_final_384.json`: 96af0a052dcd6d2c16e28741e23d3e87b96451534de2535a9b4f81e60a4ec155.
- `hd_lambda_renyi_zero_control_256.json`: 8b2bebc0ece8c3a7fecf2c9fa41ba72b20cd14098916b65944758983769aaa23.
- `hd_lambda_renyi_zero_control_384.json`: 03965aa8f8117c26ec4dc0dc600c8ce6f4b4f6a3913079de58603db978094ecc.

## Several cheap stages: exact structural screens

Adding syndrome checks can also make two-input cancellation easier. The
first inexpensive two-stage example, with supports0123 and0145, raises
the syndrome rank from five to seven but lowers the local output minimum
from six to four. The three-stage examples can lower it to two. These
are local minima, not whole-code distance counterexamples.

An exhaustive screen checks all pairs of binary weight-four supports
whose labels XOR to zero. Of their9,730 unordered pairs,4,480 do not
commute and210 are dependent modulo the original five syndrome vectors.
All5,040 remaining pairs were checked on all120 two-input supports and
all three exceptional ratios, without a symmetry quotient. The result is:

| Local two-input output minimum | Number of labelled configurations |
|---|---:|
| 2 | 168 |
| 3 | 42 |
| 4 | 3,468 |
| 5 | 72 |
| 6 | 1,290 |

The exact packed-bit screen took0.405seconds. Eighteen survivors add nine
XORs under generic circuit accounting. Screening every size-five parity
complement for each, plus both controls, took2.69seconds. Fourteen have
maximum proportional-column class four, three have five, and one has
six. The original maximum is four; the single-stage maximum is five.
This screen computes exact rank-one contributions as LOWER bounds on
moments. It does not certify their full values or the scheduled bound.

The candidate V={0,1,2,3}, W={0,1,6,7} preserves the six-output minimum
and maximum class size four. A specialized circuit needs only seven
additional XORs. Write p_i=x_(2i)+x_(2i+1). The existing circuit computes
q0=p0+p1 and r=p1+p3. Compute t=p0+p3; add r to outputs0,1, q0 to2,3,
and t to6,7. The two updates on their shared coordinates cancel to r.
The total is64 XORs and one fixed field product per column, without
new randomness. Native cost is unmeasured. The next proof test should
target middle arities, where the two additional independent checks might
help. Old certificates cannot transfer to this changed matrix.

The specialized circuit is `theory/hd_lambda_double_transvection.py`
(SHA256973175e25051621470bd19f8c6932806e4a830b908a8a1eb04186e4b616430bd).
Five focused tests check field bases, deterministic examples, equality to
the new matrix and generic circuit, involution, and exact operation count.
Including the other controls and structural helpers,120 focused tests
pass. The measured native sources and certified fallback manifest retain
their previous hashes.

## Crossing columns with a staggered XOR shear

A second construction test changes the zero-column behavior itself.
Pair columns j and j'=j XOR1 after their scheduled H_lambda maps and
before D2. For an even j, XOR its mate into column j on even physical
rows; update column j' instead on odd rows. Each row undergoes an
invertible triangular shear. The stage adds128 XORs per256-symbol group.
Final D2 stays after this stage.

Fix original active rows A, with a=|A|, and let G_j=Q_j H_lambda P_j
restricted to A. Under independent uniform local comparison variables,
the marginal zero probability is256^(-rank([G_j,diag(T_j)G_j'])), where
T_j selects updated physical rows. Both columns are zero with probability
256^-2a. These are comparison probabilities, not RS-word probabilities.

The exact rank screen covers all16 singleton supports, all120 pair
supports, every column, and all15 nonzero binary normals defining the
updated half. It also tests at most four prescribed examples per larger
arity:45,120 labelled cases in1.057seconds. For the even/odd-row shear,
every singleton has marginal rank two and every pair has rank four.
Thus their zero masses change from256^-1 to256^-2 and from256^-2 to
256^-4. Larger-arity examples are not complete coverage.

A single nonzero RS row occupies at least five column pairs, hence at
least ten output columns after this shear. That fact alone is not a
distance improvement: cancellations weaken the conservative symbol-weight
floor from144 to80. The next test should retain the joint32-output law.
For two active rows this requires four GF(4) input variables and529
relation-kernel categories, so existing exact enumeration is applicable.
Use eight chronological column-pair factors, not sixteen independent
column factors or the old fresh-ratio law.

There is an additional implementation cost. The measured kernel absorbs
one D1 scalar per column into D2. A fixed shear permits only one common
scalar per coupled pair. To preserve the ideal setup distribution,
restore eight D1 products per group relative to that measured kernel.
Storage and register costs also need measurement. Neither new candidate
has a whole-code certificate or a native timing.

The exact diagnostic receipts are:

- `hd_lambda_multi_transvection_pairs_complete.json`: b6380e17da2effa4954dd8a9d15ccbd9b7a8dd9ea14834131860867211b4829c.
- `hd_lambda_multistage_control_complete.json`: e96f3331cc30998c6d69b693efcdd7a5cbe75b908fd41a3ce0715becb8504e11.
- `hd_lambda_column_shear_ranks_complete.json`: 3ffb445067b8126fe548ed72e36e604adec88dc617608f2dd2435c97afc20762.

## Ten-row extension for the unchanged encoder

The exact moment inventory now includes all 8,008 ten-row sets, represented
by 272 verified symmetry orbits. Orders two, three, and four retain the
original source set across the sixteen scheduled columns. The inventory
took 28.12 seconds; the two fresh directed controls took 3.66 and 3.55 seconds.

For r=18 and all 99 source profiles, the added class has every group arity
in 10..16 and total arity in 180..202. Its directed failure bound is
2^-56.155410727944. Total 203 does not meet 40 bits in this comparison.
All 198 method ceilings and 66 selected ceilings agree at 256 and 384 bits.
The independently reviewed verifier unions only the passing class with
the preceding partial cover, obtaining 40.821633581232 bits. This remains
a partial r=18 certificate for the unchanged H_lambda encoder. It covers
neither all source classes nor the modified mixers below.

The source and receipts are `hd_lambda_renyi_a10_control.py`,
`hd_lambda_renyi_a10_verify.py`, and:

- `hd_lambda_renyi_a10_inventory234.json`: d40ea671e1b1ff8056241cdbe3006efb3b0b92f14ed0c41f8a5e7f95b415784b.
- `hd_lambda_renyi_a10_control_256.json`: c6b56c7a25b1bf9ea88d38e6fb75b60bb6495d820ff40e2107aecda2944cfd6d.
- `hd_lambda_renyi_a10_control_384.json`: 3f12920d9a93f1982ac4c07d9f942d0aabc2eb4b3aa584a585923639b9ff11e3.
- `hd_lambda_renyi_a10_cover_verified.json`: 05301fd9a1ffedd635f7136c23a90fcdb8e607db835def5baf14ce2f96ba54e7.

## Fresh tests of the modified mixers

The staggered shear test now retains the exact joint distribution of
each pair of columns. It enumerates all 32 output forms using 2a input
variables, including their cancellations and original input weights.
Eight chronological pair factors retain the RS-implied input-weight floor.
No independence between the paired outputs is assumed.

At the pinned r=18, Fmax<=1 measure and mark 9927/10000, the complete
singleton bound is -3775.07 and the complete two-row bound is -395.23.
These are log2 FLOAT upper bounds. The corresponding unchanged-encoder
two-row result was -165.56. This improves a selected sparse class, not
the whole-code certificate. The joint controller took 6.76 seconds.
Its normalization and chronology passed independent read-only review.

The cheaper two-stage mixer also has a fresh complete two-row control:
all 120 original row sets give a combined FLOAT bound of -235.73.
Eight prescribed middle-support tests produce the following bounds under
the same selected measure. Each fixes one original row set in every
group; the entries do not include missing sets or mixed arities.

| Original row set | Arity | Log2 FLOAT upper bound |
|---|---:|---:|
| 0..5 | 6 | +529.08 |
| 0..6 | 7 | -360.60 |
| 0..7 | 8 | -1211.47 |
| {0,1,2,3,4,5,8,12} | 8 | -229.37 |
| 0..8 | 9 | -1447.13 |
| 0..9 | 10 | -1617.33 |
| 0..10 | 11 | -963.03 |
| Complement of {0,7,11,13,14} | 11 | -905.69 |

The six-row failure is not caused by inflating quiet input weights.
All six coordinates are protected in every column, so those weights are
already exact. Its all-nonquiet component is only -893.65; mixed branches
produce the loose bound. In columns 2..5 the nonzero quiet input omits
original rows 0 and 1. Tracking one of those rows separately is therefore
a targeted next test. Each active RS row must occupy at least nine columns;
the aggregate weight floor can conceal violations of this requirement.

The new controllers and fiber helpers passed independent review. Twenty-four
focused tests cover them and the ten-row controller. No new native timing
was measured, and no old certificate was transferred to a changed matrix.
The receipts are:

- `hd_lambda_column_shear_control_9927_float.json`: 2ec123db4909fcf7e09c1e05658e1e76142b79b1829a6e4690416f22a349fda3.
- `hd_lambda_double_middle_representatives_float.json`: 5f6cf543bc57646a6aba273c8f8ab8888d27106b3b489018edb23447f10e28f3.
- `hd_lambda_double_middle_pairs_float.json`: be8b7c854721b7bf71ec82c26272f3d34f05f84045fee990ebfe3e0ff743dedc.

### Resolving the six-row diagnostic with exact local counts

The single-row floor does not repair the bound: its result is +529.083489,
versus +529.083526 for the aggregate floor. A separate dynamic program
extracts one surviving term of the same proof expression with value
+380.165715. It selects seven zero columns and nine nonquiet surrogate
columns. All six surrogate rows occupy exactly nine columns in each
replica. Thus even all six individual row floors cannot remove this term.
This is a lower bound on the chosen proof expression, not on failure
probability. The surrogate's assigned input support is not an actual
input witness.

We therefore replace the nonquiet estimate by exact joint counts. For
each input x in GF(256)^6, its GF(4)-linear relation kernel determines
the input zeros and all sixteen output zeros. A kernel of codimension d
has product_(h<d)(256-4^h) preimages. Enumerating all 564,357 positive-mass
kernels gives exact integer counts n_(w,S), where w is original input
weight and S is physical output support. The total at each w is checked
against C(6,w)255^w. The zero atom remains one.

Only 63 distinct dense output forms occur across the scheduled columns,
so array-based membership tests avoid expanding every kernel into its
elements. The exact enumeration took 3.09 seconds. The resulting selected
class bound is **-3010.594151**, replacing +529.083526. The complete control
took 3.44 seconds. Multiplicities remain outside the eighteenth roots;
the RS normalization, source prefix, chronology, and input-weight floor
are unchanged. The new helper and controller passed independent review.

This closes only the former fixed-A={0,..,5} diagnostic in FLOAT at
r=18/Fmax<=1. Other six-row sets, broader source profiles, directed replay,
and whole-code coverage remain separate obligations. The next useful
extension is a reusable category table, so broader labelled sets do not
repeat the relation-kernel enumeration.

The exact helper is `hd_lambda_double_exact6.py` (SHA256
c6065ca7f6b9605aec096283d8ff2db94630692f99a4d2deed9e20b30d4dc960).
The recorded diagnostics are:

- `hd_lambda_double_rowfloor_9927_float.json`: 53d11e69d6a9375816db9eccc81babbe7ea62de1589d21e1e7e1358400a736a5.
- `hd_lambda_double_method_path_9927_float.json`: 65deda1559f57b9d3de81047e262ae9ad2445ca7169439eb59298f104737c05d.
- `hd_lambda_double_exact6_9927_float.json`: ff18bc87635c925190869c8f12cabce6e9c07c694c70f0643d6691c82000eb26.

The reusable table is now built. It merges relation kernels by total
input weight and their exact zero pattern on all 64 dense forms
lambda*1+b, b in F2^6. This reduces 564,357 kernels to 19,414 categories
without losing any information needed by the total-weight calculation.
It does not preserve individual input-row marks. A bounded cache sorts
mapped source labels only before this total-weight calculation; actual
P_j and Q_j still act separately in every scheduled column.

The catalog took 2.84 seconds to build. A cached replay matches the full
integer histogram digest and the -3010.594151 bound of the original
six-row test. Seven additional prescribed six-row sets also pass, with
fixed-set bounds ranging from -3000.89 to -3264.33. Allowing each group
to choose any of these eight sets gives -2984.441145 under the same
selected common measure. This eight-set calculation took 2.92 seconds.
The remaining 8,000 six-row sets are not covered by this result.

- `hd_lambda_double_categories_a6_q256.npz`: 4d33d1b18dccfbfb76ae250da2fabae9c9d32e04047eaf5fa20dbcc074663692.
- `hd_lambda_double_catalog_selected8_9927_float.json`: 0529a5df563a6ea0bbcb942a3c51368faaa9b13d3cd6cc68713e2d301ba4e2c5.

### Paired-column dense comparison: room, not a bound yet

For the shear, an ideal comparison leaves substantial room at arities
9..12. Under the uniform 32-coordinate reference, the zero contribution
to the pth likelihood moment is exactly q^[p(32-2a)-32]. Jensen's
inequality gives a lower bound on the nonzero contribution. Keeping at
least five active column pairs yields selected log2 method lower
comparators -748.27, -673.73, -600.91, and -529.07 for those four arities.
Orders 2, 3, and 4 were tested at the fixed full-profile reference and
mark; order 4 wins each comparison. Directed 256/384-bit endpoints agree.

These are lower comparators for a proof method, not failure upper bounds.
They justify testing actual paired-column moment upper bounds rather than
rejecting this reference in advance. No exact middle-arity pair inventory
or native timing was obtained. The receipts are:

- `hd_lambda_shear_renyi_ideal_256.json`: 25979ad175e48654e5371c90c0d9f50ffa301211b5d10d74aca615fb1ec5312f.
- `hd_lambda_shear_renyi_ideal_384.json`: ce170f6946eb01ecd32e351f0b3f8ae3d995067ff022c2565b79acf8f09ccef3.

An actual moment upper bound is available for one twelve-row set. Its
eight paired-column codes have 8-by-32 parity matrices. Each parity
matrix admits four disjoint independent sets covering all 32 columns;
exact GF(4) elimination checks each witness. For an output support S,
the corresponding shortened code has dimension |S|-rank(B_S), so its
full-support words number at most 255^[|S|-rank(B_S)]. The four witnesses
lower-bound rank(B_S) by the largest of four intersection sizes. Summing
over these four sizes gives an exact rational fourth-moment upper bound,
without enumerating all 2^32 supports.

All eight witnesses for A={0,..,11} were obtained and checked in 0.012
seconds. The parity matrices include the inverse shear, not a forward
shear applied to the old dual. This is a local moment upper bound for
one original row set. The next step is complete labelled coverage and
the full-source directed comparison.

- `hd_lambda_shear_moment_cover_a12.json`: 26ed8841dddb328003e9555ac6a76e5fb0c7c086f70bf2d956380fce15c51ffb.

The full twelve-row inventory is now complete: all 1,820 original row
sets and their 14,560 scheduled pair matrices have verified witnesses.
Of those matrices, 13,741 have a complete four-set cover, 799 leave one
column uncovered, and 20 leave two uncovered. The exact inventory took
9.77 seconds. Search failure would not establish a rank obstruction;
only the recorded, independently checked sets enter the moment bound.

For each fixed original A, multiply its eight nonzero-moment polynomials
and retain degree at least five. Write J_A for that exact rational upper
bound. The raw group weight is
256^192/255^96 times sum_A J_A^(1/4). The 1,820 row sets occur only in
this sum; there is no additional binomial factor. The directed reference
then counts all 99 source profiles outside theta and applies the
projective normalization once.

Fresh 256- and 384-bit evaluations each recheck every parity matrix,
independence witness, and rational tail. Both give the rounded upper
failure exponent **-285.081339551997**. The certified class is r=18 with
exactly twelve active RS rows in every group, allowing every original
row set and all 99 pairing profiles. This is for the **shear construction**,
not the double-transvection mixer or unchanged H_lambda. Mixed arities,
other group counts, and native performance remain open.

Independent read-only review clears the completed inventory, both source
graphs, the normalization, and the scope of this class certificate.

The controller is `hd_lambda_shear_moment_all_a12.py` (SHA256
3544a89d48894cf6c82ec478bfac4d96f3840df9f54d49b948e1c7d2f999651c).
The two directed evaluations took 3.92 and 3.87 seconds. Seventy-five
focused tests pass across the new exact-count, reference, and witness
helpers. The completed records are:

- `hd_lambda_shear_all_a12_inventory.json`: 5a9d5b2d2adc6d3f9089167bae3e6d83c4ef7810c99c0caa73a4bc4809319f41.
- `hd_lambda_shear_all_a12_256.json`: 1daa61d0e15bd13dbea1978eb68ae3f65a0a84749a66df03d1dc93cf1160b301.
- `hd_lambda_shear_all_a12_384.json`: 3fc33c4100441d47b97ae983f4044e640eb448b04e851cdb148c295fcbd7bdc9.

### Shear continuation: broader coverage and a limit of the cover bound

The following results concern the staggered adjacent-column shear only.
They do not alter the measured H_lambda encoder or certify its modified
implementation. No native code or timing changed in this round.

First, the exact three-row calculation now covers all 560 original row
sets. Its six-input relation table has 564,357 positive-mass kernels,
compressed to 66,795 categories. The category key retains total original
input weight and the zero pattern on 80 possible output forms. Those
forms include left-only, right-only, and jointly supported forms; using
only the 64 jointly supported forms would lose information.

The complete calculation takes 74.77 seconds and gives log2 bound
**-2391.10124575756** for r=18 and Fmax<=1. The same original row set is
kept across all eight chronological pairs. This is FLOAT evidence, not
a directed certificate or coverage of all 99 source profiles. A replay
also reproduces all 136 earlier one- and two-row results. The catalog
and full three-row record are:

- `hd_lambda_shear_categories_a3_q256.npz`: 4a515aca09ca731da1c171190cc9f1c8c7bed000c220d2d1e55cba8398c64116.
- `hd_lambda_shear_all_a3_9927_float.json`: 4a62c6677e60dbcef5bafe585231f45382302255b5d7ff427bf0bba07f5553db.

Second, the twelve-row directed proof extends to other active-group
counts. A fresh regional bank counts every matching profile separately
at each r, leaving the source union outside theta. The 256- and 384-bit
replays agree on these rounded upper failure exponents:

| Active groups r | Profiles | Upper log2 failure bound, all groups of arity 12 |
| ---: | ---: | ---: |
| 16 | 80 | -95.388347173201 |
| 17 | 81 | -272.359695263099 |
| 18 | 99 | -285.081339551997 |
| 19 | 100 | -56.783934935938 |
| 20 | 120 | +172.839989759933 |

The first four classes have a combined margin exceeding 54.7839 bits,
using four times their largest upper bound. The last entry does not
certify useful distance. It is a limitation of this fixed mark and
moment bound, not a low-weight codeword. The records are:

- `hd_lambda_shear_a12_r16to20_256.json`: a964c5ff1519e4d4bd185c2c47595c5c7f34d69a8786ec638d085d0683bccdfb.
- `hd_lambda_shear_a12_r16to20_384.json`: 8a65a6aad484e766c5622ffb2addf175f555c3f1edbe4ebd733f03f283894480.

Third, a complete inventory covers all 6,885 original row sets of arities
11..16, with moment orders 2, 3, and 4. It contains 55,080 actual pair
parities and 165,240 independence witnesses. The improved search finds
full covers for every twelve-row pair; the preceding inventory remains
preserved. Construction takes 55.23 seconds, and each directed replay
rechecks all witnesses in about 33 seconds.

For each fixed p, the mixed-arity calculation raises one positive
polynomial to the eighteenth power. It selects p only after each total
arity has a complete bound. Both precision replays agree on all 91
totals. Only total 198 clears 40 bits: it is the all-eleven-row class,
with margin **45.697980665791** bits. Total 199 gives -30.164884312708;
the broader mixed class does not close. This does not contradict the
separate all-twelve-row result: mixed totals also allow other row counts.

- `hd_lambda_shear_neighbors_inventory.json`: 0e05655a6d982c0b998d81359e2bbcadb6a83cfadf30369a86a92796aa54ec66.
- `hd_lambda_shear_neighbors_256.json`: fa28e45fdd1469e7e02e43e353148ef81780183e884b0dc26c1959beb213fc3f.
- `hd_lambda_shear_neighbors_384.json`: 7e044892f4e7d055224fbca7dac537675e3abff0da2cddccc176a61cf64fa98d.

Before extending the inventory to arities 4..10, we checked whether even
ideal independence witnesses could make this relaxation succeed. The
check minimizes the cover expression over every abstract set-size profile,
ignoring whether the actual parity matrix permits it. Adding a covered
coordinate cannot worsen that expression, so maximal-capacity profiles
suffice. This gives a lower comparator for the cover method itself.

Even these optimistic values miss the target at all seven arities. Their
best log2 expressions are +1781.63, +1445.41, +1106.30, +764.60, +348.52,
+356.28, and +90.46. Directed endpoints agree at both precisions. Therefore
larger inventories or better packing searches cannot repair this specific
rank-cover relaxation at the current mark, activity floor, and p=2/3/4.
This result does not reject the construction or stronger proof methods.

- `hd_lambda_shear_cover_limits_256.json`: dc1a89e54f1e4645bc9d9e792dc6aa9388a1fa2a34bc80b2e3c6720ef8f5de9a.
- `hd_lambda_shear_cover_limits_384.json`: 3961f60f9ed61700ee5b0c32c832ecaeaedc262fa4c339cddc6a139cfd86efcd.

Thirty-seven focused tests pass for the new controllers and exact catalog.
Independent review clears the source normalization and stated scopes.
The next step is a shared-scalar cancellation bound, not more of the same
middle-arity cover search. Mixed sparse/dense arities, broader group
counts, directed sparse evaluation, and native performance remain open.

### Next bound: one shared scalar, correlated cancellations

Fix the two pre-shear vectors up to nonzero scaling. Full independent D1
leaves their relative scale rho uniform in F256*. Scaling preserves each
column's input weight and support, so this conditioning is compatible
with the exact single-column histograms. It requires retaining one
relative scale per pair; the old two separately removed gauges would
invalidate the argument.

For each overlapping row, exactly one value of rho cancels its shear
target. Different cancellation values therefore partition the overlapping
rows into disjoint classes. Consider scalar output marks in [0,1]. Merging
two cancellation classes C,D increases their contribution by
the unchanged-coordinate product times
(1-product_{i in C} t_i)(1-product_{i in D} t_i), which is nonnegative.
Repeated merging gives the bound

    E_rho M <= (254/255) M_OR + (1/255) M_cancel.

Here M_OR retains every overlapping target. M_cancel deletes all those
targets but retains the unmodified donor coordinates. For even physical
rows the target is left; for odd rows it is right. This is not independent
cancellation at each coordinate, nor does it cancel both endpoints.

Apply the inequality after passing to the common positive matrix measure.
Its scalar atoms lie in [0,1], and its nonnegative state-path weights
preserve the inequality. Minkowski then bounds the r-norm by
(254/255)||M_OR||_r+(1/255)||M_cancel||_r. The coefficients remain outside
the roots. This is a bound on the marked expression, not coefficientwise
support domination or a likelihood-moment bound. It uses the existing
shear distribution without adding encoder work.

The first bounded test passes. Under the same r18/Fmax<=1 measure,
the fixed original sets {0,1,2}, {0,1,2,3}, and {0,1,4,8} give FLOAT
log2 bounds -1452.027886, -1375.053427, and -1099.295023. Allowing each
group to choose among these three sets gives -1099.294355. The complete
test takes 0.447 seconds after preparing one common measure. No input
mark fitting or regional-provider call is involved.

The three-row result is looser than its exact joint-support result,
-3077.961837. That is expected: the new argument permits an exceptional
scalar to cancel every overlap target, whether or not a particular
projective pair can realize that event. Crucially, the two four-row
tests do not require enumerating relation kernels on eight variables.
They are selected examples, not the complete class of 1,820 row sets.

Eight small tests check exact GF4 averages, true marked moment polynomials,
noncommuting state matrices, physical shear targets, and ratio normalization.
Five controller tests check chronology, original-degree floors, and timeout
scope. Independent review clears both helpers. Together with the preceding
37 tests, all 50 new focused tests pass.

- `hd_lambda_shear_shared_ratio.py`: 5b60c22fa6e120c9e4faee276de0f029dc52121c21dc53c2d7c12eafdac9b363.
- `hd_lambda_shear_shared_control.py`: 933a8a423af5b65a00c6a070d68b66eeb14171b5c00328a492eb0311fbf9b38d.
- `hd_lambda_shear_shared_selected_float.json`: 3fb5403e32ac876c51a2ab45fe0193888220fcf33ebf653986a22e6d2b679e15.

Next, test every four-row set in bounded batches, then extend exact
single-column histograms to larger row counts if these bounds retain
room. Full profile coverage and directed arithmetic remain separate work.
