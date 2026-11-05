# Fractional mask comparison for the g16 two-sum outer

This is a proposed proof route, not a numerical certificate.
It applies to the g16 two-sum outer with labels $06+\operatorname{GF}(16)$, sixteen RS$[16,8]$ constituents, the alternating independently swapped pair route, and the existing h4 DHD inner map.
The construction identifier is `two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half`.
The bounded FLOAT controls recorded below do not certify a whole-code bound.

Write $u$ for the ordered 16 labels in `LABELS` in `theory/hd_two_sum_gf16_sparse.py`.
The column mixer is $H=I+\mathbf1_{16}\mathbf1_{16}^T+uu^T$ over AES $\operatorname{GF}(256)$, where $\mathbf1_{16}$ is the all-ones vector.

The purpose is to apply a fractional moment before averaging the route.
The previously used forced-nonzero regional matrices have already averaged that randomness.
They cannot be silently reused as the finer fractional provider defined below.

## Fix the class before averaging setup

Fix an active source-group set $G$ and, for each $g\in G$, a nonzero constituent-row set $A_g\subseteq[16]$.
Each $A_g$ remains fixed across all 16 outer columns.
Define the total arity by $A_{\rm tot}=\sum_{g\in G}|A_g|$.
Let $\Pi$ denote all regional pair permutations and fair swap bits.
The outer diagonals and inner refresh setup are independent of $\Pi$ under the ideal setup law.

For one fixed row set $A$ of size $a$, define

$$
N_A(T)=\#\{x\in\operatorname{GF}(256)^A:
                 \operatorname{supp}(Hx)=T\},\qquad
c_A(T)=\frac{N_A(T)}{255^{a/2}}.
$$

The input coordinates in this counting law may be zero.
The closed kernel in `theory/hd_two_sum_closed_kernel.py` supplies the exact counts $N_A(T)$.
These mask weights are not probabilities.

The prescribed-support RS comparison contributes $255^{-8a}$ for the whole constituent array.
Distributing this factor over 16 columns gives

$$
\prod_{j=1}^{16}c_A(T_j)
   =255^{-8a}\prod_{j=1}^{16}N_A(T_j).
$$

Conditional on the complete output-mask word, the outer output diagonal supplies independent uniform nonzero amplitudes.
The comparison restricts each active group to at least nine nonempty columns.
This restriction is necessary for every actual nonzero group because RS$[16,8]$ has minimum distance nine and the column mixer is invertible.
It need not enforce every individual constituent's support floor to remain an upper bound.

Let $\mathcal T$ be the retained set of mask arrays, including this restriction, and put

$$
c_{\mathbf A}(\mathbf T)
   =\prod_{g\in G}\prod_{j=1}^{16}c_{A_g}(T_{g,j}),\qquad
K=\frac{z^{-L}}{255},\quad L=734003,\quad 0<z<1.
$$

For fixed $\Pi,\mathbf T$, let $\Phi_{\mathbf T}(\Pi)$ be the mean of
$z^{\operatorname{wt}(Y)}$ over the conditional nonzero amplitudes and all inner refresh keys.
If $B_{G,\mathbf A}$ is the event that this message class has a nonzero codeword of weight at most $L$, the conditional first moment gives

$$
\Pr[B_{G,\mathbf A}\mid\Pi]
 \leq\min\left(1,K\sum_{\mathbf T\in\mathcal T}
                 c_{\mathbf A}(\mathbf T)\Phi_{\mathbf T}(\Pi)\right).
$$

The factor $255$ removes projectively equivalent nonzero codewords.
All outer and inner setup has already been averaged in this conditional first moment; only the route remains random.

## Fractional power before route averaging

Fix $0<\alpha\leq1$.
The inequalities $\min(1,x)\leq x^\alpha$ and
$(\sum_i x_i)^\alpha\leq\sum_i x_i^\alpha$ give

$$
\Pr[B_{G,\mathbf A}]
 \leq K^\alpha\sum_{\mathbf T\in\mathcal T}
        c_{\mathbf A}(\mathbf T)^\alpha
        \mathbb E_{\Pi}[\Phi_{\mathbf T}(\Pi)^\alpha].
$$

For one packet, let $M^{\rm NZ}_z(o)$ be the existing nonnegative two-state cap for block occupancies $o=(o_0,o_1,o_2,o_3)$.
Its input law is independent uniform nonzero amplitudes at the specified positions.
Its two coordinates track zero-state mass and the established density bound over nonzero states.
It includes the usual averaged DHD refresh-density and empty-packet corrections.
Write $e_0=(1,0)^T$ for the initial zero state and $\mathbf1=(1,1)^T$ for the unrestricted terminal vector.

For fixed $\Pi,\mathbf T$, the established packet argument bounds the mean by

$$
\Phi_{\mathbf T}(\Pi)
 \leq e_0^T\left(\prod_{\ell\ \mathrm{in\ chronological\ order}}
                     M^{\rm NZ}_z(o_\ell)\right)\mathbf1.
$$

Expanding this positive product over the two-state paths before taking the fractional power yields

$$
\Phi_{\mathbf T}(\Pi)^\alpha
 \leq e_0^T\left(\prod_\ell
                     M^{\rm NZ}_z(o_\ell)^{\circ\alpha}\right)\mathbf1,
$$

where $\circ\alpha$ denotes entrywise powers, including $0^\alpha=0$.
No claim that $M^{\rm NZ}_z$ bounds each fixed refresh-key realization is needed.
This argument exploits variation of the conditional first moment across routes; it does not preserve fluctuations of the already averaged inner setup.

## New regional provider and its normalization

Fix a regional support with $m$ entries and $f$ full source pairs; write $s=m-2f$.
Let $n(s',f',o)$ count oriented arrangements with $s'$ singleton pairs, $f'$ full pairs, and occupancy profile $o$ in a 32-pair packet.
Define the matrix polynomial

$$
B_\alpha(X,Y)=\sum_{s',f',o}
 \frac{n(s',f',o)}{2^{s'}}
 M^{\rm NZ}_z(o)^{\circ\alpha}X^{s'}Y^{f'}.
$$

The fractional regional upper matrix is

$$
R_\alpha(m,f)=
\frac{[X^sY^f]B_\alpha(X,Y)^{128}}
     {\binom{4096}{s}\binom{4096-s}{f}}.
$$

The fair-swap factors, placement multiplicities, and regional denominator are **not** raised to $\alpha$.
They average the remaining route randomness after the packet entries were powered.
Independent regional permutations justify composing these regional matrices in chronological order for each fixed mask word.
The regional support sets and their full-pair counts may still be correlated across regions and columns.
An entrywise maximum over feasible $f$ is a universal upper bound, not a resampling of those sets.

The existing forced-nonzero packet cap formulas therefore suffice; the new provider must change the order of packet powers and placement averaging.
Raising an already averaged regional matrix entrywise to $\alpha$ is a separate, coarser Jensen bound, not the coefficient formula above.
The uniform-amplitude packet provider has a different input law and cannot replace this one for a general fixed output mask.

## Fixed-row-set Hölder weights

For arity $a$, the correct fractional weight is

$$
W_{a,\alpha}(w)=
\left(\frac{\sum_{|A|=a}N_A(T)^{16\alpha}}
                  {255^{8a\alpha}}\right)^{1/16},
\qquad |T|=w.
$$

The closed-kernel symmetry makes this depend only on $w$.
For any fixed mask word, Hölder over the same row set in all columns gives

$$
\sum_{|A|=a}\prod_{j=1}^{16}c_A(T_j)^\alpha
 \leq\prod_{j=1}^{16}W_{a,\alpha}(|T_j|).
$$

The row-set multiplicity is included once in this sum; there is no additional $\binom{16}{a}$.
In general, $W_{a,1}(w)^\alpha$ is smaller than $W_{a,\alpha}(w)$ and is not a valid substitute.
Coincident closed-kernel branches for arities 15 and 16 must still be combined before taking $16\alpha$ powers.

The mask-first order also preserves the amplitude normalization.
At an amplitude vector with support $T$, the corresponding point mass is $c_A(T)/255^{|T|}$.
Powering these point masses and then summing amplitudes would introduce the factor $255^{(1-\alpha)|T|}$.
The mask-first derivation instead keeps the normalized nonzero-amplitude expectation inside its fractional power and introduces no such factor.

If the stronger sparse two-row comparison is retained, its column branches must also be fractionalized before mask summation.
The singleton branches have weights $255^\alpha$, with their two labelled choices counted outside the power.
The full-output double branch has weight $(255\cdot239)^\alpha$; each specified one-hole branch has weight $255^\alpha$.
The original $239/255,1/255$ probability law cannot simply be reused unchanged at $\alpha<1$.
Whole-word RS denominators are powered; support-floor markers inserted after this step are ordinary positive counting markers.

## A stronger uniform comparison for large total arity

The mask comparison is not the only available bound.
For sufficiently large total arity, a uniform comparison avoids the row-set sum entirely.
This argument uses the ideal independent outer diagonal law before applying any native gauge normalization.

Fix $G$, with $|G|=r$, and a target outer array $y$ on its $256r$ coordinates.
Fix an integer threshold $A_0$ with $r\leq A_0\leq16r$.
Condition on all outer output diagonals $D_{\rm out}$.
For every group and column, invert the output diagonal and $H$:

$$
v_{g,j}=H^{-1}D_{{\rm out},g,j}^{-1}y_{g,j}.
$$

The resulting array $v$ uniquely determines its nonzero constituent rows.
Nonzero input-diagonal scalings $D_{\rm in}$ cannot change those row supports.
Thus only one row-set vector can produce this target for the fixed $D_{\rm out}$.

For a prescribed nonzero row vector of support size $w$, averaging its independent input diagonal gives RS membership mass $B_w/255^w$.
Here $B_w$ is the number of RS$[16,8]$ words with that exact support.
It vanishes for $1\leq w\leq8$.
For $w\geq9$, the shortened code has dimension $w-8$; projection to an information set gives $B_w\leq255^{w-8}$.
Each prescribed nonzero row therefore contributes at most $255^{-8}$.
A zero row contributes exactly one through the zero constituent word.

Write $A(v)$ for the number of nonzero rows in $v$, summed over $G$.
Independence of the input diagonals across rows gives the following bound for the expected counting measure of messages with total arity at least $A_0$:

$$
\mu_{G,\,A_{\rm tot}\geq A_0}(y\mid D_{\rm out})
 \leq \mathbf1_{\{A(v)\geq A_0\}}255^{-8A(v)}
 \leq255^{-8A_0}.
$$

The requirement that every group in $G$ be active can only reduce this measure.
Averaging $D_{\rm out}$ preserves the pointwise upper bound.
There is no factor counting row sets, cohort assignments, or arity deficits: those possibilities were uniquely determined before the input diagonals were averaged.
The same bound applies to an exact arity or an interval of arities with lower endpoint $A_0$.

Consequently,

$$
\mu_{G,\,A_{\rm tot}\geq A_0}
 \leq \beta_{r,A_0}\,\operatorname{Uniform}
                 (\operatorname{GF}(256)^{256r}),\qquad
\beta_{r,A_0}=\frac{256^{256r}}{255^{8A_0}}.
$$

This comparison includes potential zero values and can use the established **uniform-amplitude** swap operator.
Its fractional prefix is $(\beta_{r,A_0}z^{-L}/255)^\alpha$.
The union over source supports remains outside this power.
The sharper density specific to invertible Toeplitz mixers does not transfer to this construction.
At $A_0=16r$, the formula reduces to the legacy full-arity density.

For a fixed legacy uniform-provider witness, changing only this density adds

$$
8\alpha(16r-A_0)\log_2 255
$$

to its log2 bound.
Saved legacy sparse receipts suggest the following starting points for a future partition control:

| $r$ | $z$ | $\alpha$ | Legacy log2 bound | Proposed $A_0$ |
|---:|---|---|---:|---:|
| 18 | $691/700$ | $1/5$ | $-588.371637445129$ | 247 |
| 24 | $172/175$ | $1/5$ | $-921.507348929826$ | 317 |
| 32 | $171/175$ | $1/5$ | $-690.057289231772$ | 463 |

These are saved-prefix FLOAT diagnostics, not new g16 receipts or a covered interval.
The legacy source is `theory/hd_shifted_sparse_checkpoint_directed384.json`, SHA-256
`3997a8e66a6f9356faceda92fd18ab857c1e8e1ec3ebda99f2211c5dd7689d0d`.
Its scalar prefix uses the legacy density $256^{256r}/255^{128r}$.
The proposed cutoffs leave approximately 60 bits of high-arity headroom under the displayed adjustment; the low-arity complement remains to be bounded.

## Complementary low-arity bound with one marker

Partition the message classes into $A_{\rm tot}\geq A_0$ and $A_{\rm tot}\leq A_0-1$.
Use the preceding uniform comparison for the first class and the forced-nonzero mask comparison for the second.
These are two separately justified bounds; their local input laws are not interchangeable.

Let an ordered type vector $\tau=(a_1,\ldots,a_r)$ specify the fixed arity of each source group.
For $0<t<1$ and $B=A_0-1$,

$$
\mathbf1_{\{\sum_g a_g\leq B\}}
 \leq t^{\sum_g a_g-B}.
$$

This indicator acts on the union over fixed types after the fractional bound for each class.
It is therefore not an $\alpha$-powered marker.
If $F_a(P)$ is the column factor for type $a$ and regional state path $P$, replace it by

$$
F_{a,t}(P)=t^{a/16}F_a(P)
$$

before applying Hölder over fixed type vectors across the 16 columns.
For a fixed $\tau$, these factors multiply to $t^{\sum_g a_g}$ over the whole word.
The outside factor is $t^{-B}$, applied once.
The same $t$ must be used for every path, region, column, and group in this bound.

For example, the path-Minkowski expression for $r$ ordered group types changes from

$$
\left(\sum_{a=1}^{16}F_a(P)^{16}\right)^{r/16}
\quad\text{to}\quad
\left(\sum_{a=1}^{16}t^a F_a(P)^{16}\right)^{r/16}.
$$

Existing fixed-row-set multiplicities remain inside $F_a$; the marker introduces no additional row-set or cohort count.
Choosing $A_0$ and $t$ is a deterministic proof-witness choice, not an adaptation to sampled setup.
The union of the two resulting probability bounds covers all arities at that fixed source weight $r$.
Neither bound is asserted to meet the target until a numerical calculation and directed replay establish it.

## Bounded FLOAT control and initialization loss

`theory/hd_two_sum_gf16_fractional_state.py` implements the fractional mask route for all source supports and all cohort splits with $18\leq r\leq32$.
It retains the two boundary states through the fixed-type Hölder and path-Minkowski bound.
This control does not implement the high/low-total-arity partition above.

The frozen run built five regional banks and evaluated five deterministic type-specific floor templates per bank.
It then applied one coordinate-selection round only to the globally best baseline.
The 25 baselines had the following best log2 bounds within each bank:

| $z$ | $\alpha$ | Best baseline log2 bound |
|---|---:|---:|
| $99/100$ | $1$ | $+16496.0788951891$ |
| $24/25$ | $1/2$ | $+43191.3162761845$ |
| $9/10$ | $1/4$ | $+72866.5481277883$ |
| $3/4$ | $1/8$ | $+112457.520058798$ |
| $1/2$ | $1/16$ | $+143997.516136092$ |

All failed.
The selected refinement stayed at $(z,\alpha)=(99/100,1)$ and recorded 44 slope evaluations, including the initial slope vector.
Its union screen was $+16494.53668045$ bits.
The proof calculation took 2.6942441 seconds; this is not a native-encoder benchmark.

The fractional baselines were severely confounded by their initialized live-to-live envelope.
Write $\ell_{11}(m)$ for the maximum regional log2 live-to-live entry at occupancy $m$, and bound it by $a_{11}+b_{11}m$.
The initialization chose the secant between the first and last finite occupancies for $b_{11}$, then raised $a_{11}$ to cover every occupancy.
The curves at $\alpha<1$ have a deep value at $m=0,1$ followed by a sharp rise at $m=2$.
Their initialized positive slopes therefore produce very loose upper envelopes at larger occupancies.

For $(z,\alpha)=(1/2,1/16)$, the saved curve gives

$$
\ell_{11}(0)\approx-508.396834486,\qquad
\ell_{11}(2)\approx-21.6164897194,\qquad
\ell_{11}(32)\approx-84.8618180610.
$$

The initialized envelope has $a_{11}\approx-48.0874282459$ and $b_{11}\approx13.2354692633$.
At $m=32$ it gives approximately $+375.447588179$, exceeding the saved curve by **460.309406 bits in one region**.
The corresponding $m=32$ excesses for $\alpha=1/2,1/4,1/8$ are approximately 151.728, 249.733, and 372.021 bits.
These are genuine upper-envelope losses, not evidence of floating-point underflow.
None of those four fractional banks received a slope-refinement round because their initial bounds lost to the $\alpha=1$ bank.

The floor templates suppress the admitted all-empty comparison term to at most $-60$ bits at their anchor weights, before envelope-intercept losses.
They can simultaneously increase terms with nonempty columns.
Thus the failed run is evidence about this bounded initialization and selection rule, not an obstruction to fractional moments or to the construction's distance.
A saved-bank repair of the live-to-live slope can isolate this loss without changing any packet operator, weight mark, or fractional exponent.

The frozen result is `theory/hd_two_sum_gf16_fractional_state_r18_32_float.json`, SHA-256
`d1f4a8aa1d57b8c63f8bcd27a9af9dedac5274fc0b9496d3f5e35558cf07d544`.
Its source hash is
`bbeca202d77780ea325bbc70e959e73e929beec7c0341c2d80d85023c96273d8`.
The source's ten lightweight tests passed before the run.
There is no directed fractional receipt.

### Saved-bank live-to-live repair

`theory/hd_two_sum_gf16_fractional_repair.py` subsequently reused the same five saved FLOAT banks.
It re-evaluated all 25 baselines and varied only $b_{11}$ over the existing finite secant candidates in the first phase.
For every candidate it recomputed the upper intercept against all occupancies.
The other three slopes, all weight marks, all exponents, and all type-specific floor markers remained fixed.
No regional operator was rebuilt.
The best complete-band result then received one additional four-coordinate round.

The best bounds after the live-to-live-only phase were

| $\alpha$ | Best repaired log2 bound |
|---:|---:|
| $1$ | $+16496.0788951029$ |
| $1/2$ | $+21288.2184054891$ |
| $1/4$ | $+27867.3446599236$ |
| $1/8$ | $+38798.0725731014$ |
| $1/16$ | $+47486.7211703837$ |

The repair removed a large part of the fractional initialization loss, but all tested bounds still failed.
The globally selected witness remained $(z,\alpha)=(99/100,1)$.
Its final additional round again gave the union screen $+16494.53668045$ bits.
The saved-bank calculation completed in 6.0967669 seconds.
The remaining failure concerns this finite witness and envelope family; it is not a proof of insufficient code distance or an impossibility result for the fractional method.

The repair artifact is `theory/hd_two_sum_gf16_fractional_repair_r18_32_float.json`, SHA-256
`592b49e9667c519de712b7825984931400c5620be39df9534a634482b611f774`.
Its reviewed source hash is
`83fc3f5b3320f43d1c6523c5dbacc5d6d0a4751d323a3fc9b9153a3947934c4d`.
After this run, the complete lightweight `hd_two_sum*_tests.py` regression passed all 166 tests in 6.679 seconds.
No further operator banks were built in this cycle.

## Scope and next implementation check

The source supports, constituent-row sets, and cohort assignments are fixed classes before this argument.
The complete code event is bounded by their union, so their outer counting factors are not raised to $\alpha$.
Hölder may compress those fixed-class sums, but it does not make row sets independent between columns.

A directed implementation can reuse the positive bivariate regional convolution and the exact closed-kernel branch counts.
It must rebuild the fractional mask moments and forced-nonzero packet-power provider with directed arithmetic, using the prefix $(z^{-L}/255)^\alpha$.
The controls should retain the existing $\alpha=1$ checks, packet-power chronology tests, and fixed-row-set Hölder comparisons on labelled toy masks.
Any numerical proposal still requires a fresh directed replay and an explicit covered message class.
The next bounded control should combine the uniform high-arity bound with the marked low-arity complement, using the established near-$0.98$, $\alpha=1/5$ witnesses above.
Repeating a broad search on the rough fractional grid is not the next recommended step.
