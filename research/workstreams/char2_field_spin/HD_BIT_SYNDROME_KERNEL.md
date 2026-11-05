# A positive nonquiet kernel for the bit-syndrome outer

The [four-XOR proposal](HD_BIT_SYNDROME_OUTER.md) now has a complete
positive column majorant for a fixed eight-row set. The majorant includes
nonquiet words and can restore all singleton and pair inputs exactly.
It uses exact integer counts, not the old one-hole law. This helper runs
no distance contraction, regional provider, or native experiment.

The remaining relaxation is explicit: some output coordinates inside the
fixed source set are bounded without retaining their values. The helper
supports all finite nonnegative atomic markers; it does not assume that
the current moment embedding always has atoms at most one.

## Fixed row set and outside image

Fix `A subset {0,...,15}` of size eight, and let `O=A^c`. Throughout all
later columns, A is unchanged. Write `r=rank(V_A)`, which is four or five.
The previous exact checks establish

\[
 \operatorname{rank}H_{*,O,A}=\operatorname{rank}V_A
 =\operatorname{rank}V_O=r.
\]

Because `H*_(O,A)=V_O B V_A^T`, its image is exactly `im(V_O)`.
Consequently, the outside output of the counting measure on `F^A` ranges
over this affine-evaluation code. Every outside value has `q^(8-r)`
preimages, where `q=256`.

The zero outside value corresponds exactly to quiet inputs. Differences
between two preimages of the same outside value lie in

\[
 C_A:=\ker(V_A^T).
\]

The matrix H* fixes this kernel. The inside outputs of each outside fiber
therefore form a translated copy of the same `8-r` dimensional code C_A.
This description is conditional on A, not an average that selects a new
row set in each column.

Choose an information set `I subset A` for C_A. One can take the complement
in A of any basis of the columns of `V_A^T`. Projection of each translated
fiber onto I is a bijection with `F^I`.

## Exact outside counts and the positive projection

For `S subset O`, let E_O(S) count outside codewords with support exactly
S. The number supported within S is

\[
 q^{r-\operatorname{rank}(V_{O\setminus S})}.
\]

An exact Möbius inversion over the 256 subsets of O computes E_O. Define
the nonzero outside polynomial

\[
 C_O(t):=\sum_{\varnothing\ne S\subseteq O}E_O(S)\prod_{i\in S}t_i.
\]

For input marker `X>=1` and output markers `0<=t_i<=1`, the full column
count is bounded by

\[
 K_A(t,X)\le
 \underbrace{\sum_{S\subseteq A} Z_*(S)X^{|S|}\prod_{i\in S}t_i}_{Q_A(t,X)}
 +\underbrace{X^8 C_O(t)\prod_{i\in I}(1+255t_i)}_{P_A(t,X)}.
\]

The first term retains quiet inputs exactly. The second term counts each
nonzero outside fiber, retains its information coordinates, and drops
only the constrained inside output coordinates. Its input exponent eight
upper-bounds the true input weight. At `X=t_i=1`, the sum is exactly
`q^8`: there is no additional counting or normalization factor.

For arbitrary finite `t_i>=0`, multiply P_A by

\[
 D_A(t):=\prod_{i\in A\setminus I}\max(1,t_i).
\]

This bounds every omitted coordinate in either its zero or nonzero case.
The helper's `evaluate` function deliberately rejects atoms above one;
`evaluate_kernel` includes D_A and supports the general domain. The
common moment interface itself does not enforce atoms at most one. Its
frozen `691/700` example happens to satisfy that restriction, but the
general bound does not rely on it.

## Exact sparse restoration

The projected polynomial contains every nonquiet input, including its
singleton and pair inputs. Their projected multiplicities are known
exactly from the new labelled pair laws. The helper removes those integer
multiplicities from P_A before adding their true output masks and true
input exponents. It checks that every remaining coefficient is nonnegative.

This is not subtraction of two rounded upper bounds. It is deletion of
an exactly counted subset from the exact projected counting measure.
The general-marker factor D_A is applied only to the remaining projected
terms. Quiet and restored sparse terms retain their actual masks.
The restored kernel still has total mass `q^8` at unit markers.

No RS normalization, projective factor, or source/row-set multiplicity is
included in these column counts. The existing counting-measure factor
`255^(-8a)` for a fixed arity-a RS group belongs outside the complete
sixteen-column product. The comparison must also retain its required
input-support floor. For arity eight the helper's upper input activity
does not create a new floor; it only relaxes a previously justified one.

## Four count types, not four chronological source types

The evaluation code `im(V_O)` and the quiet counts are invariant under
binary affine changes of coordinate labels. The eight-point subsets
form four orbits, with sizes 30, 1,920, 840, and 10,080. Three adjacent bit
swaps, one elementary shear, and one translation generate the affine
group used by the exact orbit walk.

The rank-four orbit consists of the thirty affine hyperplanes. The other
three orbits have rank five. Taking coefficient maxima by output weight
over representatives gives the following universal nonzero envelopes:

| Weight w | Rank-four coefficient cap | Rank-five coefficient cap |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 0 | 255 |
| 2 | 0 | 255 |
| 3 | 0 | 255 |
| 4 | 255 | 65,025 |
| 5 | 0 | 65,025 |
| 6 | 64,770 | 16,581,375 |
| 7 | 16,322,040 | 4,178,895,885 |
| 8 | 4,162,573,845 | 1,065,618,450,675 |

If c_r(w) denotes a table entry, then
`C_O(t)<=sum_w c_r(w)e_w(t_O)` for every fixed O of rank r and every
nonnegative marker vector. No orbit size multiplies this universal
coefficient envelope. The exact quiet coefficient maxima through weight
eight are `(1,0,0,0,255,0,64770,16322040,4162573845)`.

These symmetries do not justify collapsing the final chronological sum
over A. The actual H* pair masks need not be affine-invariant, and regional
operators retain their labelled order. In particular, the exact sparse
restoration uses the actual fixed pair, not an orbit representative.

## An economical next contraction

The primary next gate must use the **original one-state inner**, which is
the inner retained by the native bit-syndrome candidate. The existing
one-state, unpowered Fmax1 curves in
`theory/hd_two_sum_late_alpha_live_gate_float.json` are a possible fixed-mark
input at `z=691/700`; their artifact SHA256 is
`33da1eb9517c7d5665d22d2c8b93f5622199d4e1b70dd2a8e3aa52e28d4e04cf`.
Only the regional curves are relevant, not that artifact's old outer
counting bound or its failed common-matrix factorization. Their reuse
requires checking the same field, packet map, route, occupancy range,
and independent nonzero conditional amplitudes.

The `post_refresh_*_rows` helpers supply a reusable contraction pattern,
not authority to use the post-refresh regional bank for the measured
one-state encoder. Any exploratory calculation with that eight-cone bank
must identify a different inner construction and cannot use the one-state
native timing.

The un-restored projection envelope admits a bounded version of the
existing five-atom, eight-marked-position dynamic program. Let a common
positive matrix-valued measure bound the regional matrices as

\[
 R(m)\preceq\sum_\ell W_\ell t_\ell^m,\qquad W_\ell\succeq0.
\]

Use a two-branch Young inequality for `(Q_A+D_A P_A)^18`. The nonquiet
branch has an outside polynomial raised to eighteen. The inside positions
contribute independent scalar factors to the same atom measure:

\[
 E_I=\sum_\ell W_\ell(1+255t_\ell)^{18},\qquad
 E_D=\sum_\ell W_\ell\max(1,t_\ell)^{18}.
\]

Mark the eight outside positions in the usual 495-count DP. At positions
in I use E_I; at the other inside positions use E_D. The gaps are ordered
products of those two matrices, not powers of one homogeneous gap matrix.
The input factor `X^144` is applied once to the nonquiet column branch.
The quiet branch uses the same underlying measure with its usual empty
gaps and quiet coefficient envelope.

Add the two positive column matrices before the sixteen-column power.
Whole-word group Minkowski may then sum over every labelled fixed A.
All gap matrices here integrate scalar functions against the same W_l;
none substitutes an A-dependent exact empty regional matrix.

The contraction is now implemented as a FLOAT comparison. It avoids a
coupled count DP over all sixteen positions. An alternative translated-
fiber bound averages information sets by Hölder. The exact counting helper
returns rational inclusion frequencies; the contraction engine evaluates
the resulting positive integrand, as described below.

## Validation and scope

Six new tests and eight earlier outer tests pass in 1.239 seconds.
Independent GF(4) enumeration checks the outside evaluation counts,
translated fibers, projected multiplicities, and both marker domains.
It also checks the information-set Hölder inequality on translated
repetition fibers. Production-sized exact tests check positive sparse
restoration, unit mass, orbit coverage, and coefficient-envelope samples.
An independent source review found no defect in the fiber, restoration,
marker, or count-orbit arguments.

Reproduce from `theory`:

```text
python -m unittest hd_bit_syndrome_kernel_tests hd_bit_syndrome_outer_tests
```

Kernel SHA256:
`a27e5557aceaf04273c3e67ce0a01a5d5746c6d01a2ecda8881205ba23541c77`.
Test SHA256:
`a0a7c250f40c45c6e5cb90124891f3df6a4ab29ccb9f8993e2172f47ce027c3e`.

## Full-state comparison with an information-set cover

For fixed A, give its quiet-code information sets equal weight. Let c_i
be the inclusion frequency of coordinate i. For every translated fiber,
generalized Hölder bounds its weighted count by

\[
 H_A(t)=\prod_{i:c_i>0}(1+255t_i^{1/c_i})^{c_i}
         \prod_{i:c_i=0}\max(1,t_i).
\]

This is an alternative to selecting one information set; it is not
asserted to improve that choice pointwise. At unit markers its value is
exactly the fiber size, since the c_i sum to the kernel dimension.
Each coordinate factor, raised to eighteen, is integrated against the
same common matrix-valued atom measure. Ordered products of these small
matrices fill the gaps between the eight outside positions. The dynamic
program still has 495 terminal count classes and retains all 12,870
labelled A. Affine transport prepares the exact c_i only; it does not
identify their chronological operators.

The first complete selected-class comparison uses the original one-state
regional curves at z=691/700 and X=16. The input floor is 1296, the RS
normalizer is 255^-1152, and one projective division contributes 255^-1.
The exact Fmax1 source count is included once. Independent review checks
the conditional nonzero-amplitude law, gauge normalization, initial zero
state, terminal all-ones vector, and these factors. The reused artifact's
old outer bound is not used.

The full calculation finishes in 45.11 seconds. Three predeclared Young
weight pairs, (1/2,1/2), (99/100,1/100), and (1/100,99/100), give log2 upper
bounds +9459.02, +10994.15, and **+9190.96**, respectively. Thus the bound
is vacuous. The unweighted quiet and nonquiet components are +1764.38 and
+9187.02. Rank-five row sets dominate; this is not merely the thirty
rank-four hyperplanes. These are failed upper bounds, not witnesses of
low codeword weight.

The computation and branch matrices are saved as
`theory/hd_bit_syndrome_cover_x4_float.json` and its `.npz` cache. The
engine is `hd_bit_syndrome_gate.py`, SHA256
`3f8fa49998b1aacc00a964e5c6b83611d62510424afeca98875979aea8063505`;
the pinned loader is `hd_bit_syndrome_gate_run.py`, SHA256
`4b938c7d8fbd46c81ffb4fc821143c88b95722e12d1228b442a23d3583b00a23`.
Twenty-nine outer, kernel, contraction, and loader tests pass, including
exact small noncommuting path controls and atoms above one. Independent
review covers the engine and wrapper.

The next bounded refinement separates the four exact row-set count types.
In particular, only 1,920 sets permit an outside word of weight one;
the universal rank-five envelope allows that coefficient for 12,840 sets.
Every labelled A must still be retained. A second prescribed output mark
can avoid the large nonquiet X^144 input-activity relaxation. Neither
refinement changes the implementation or supplies a certificate by itself.

## Refined selected-class result

The orbit-specific engine now uses the exact four coefficient types while
retaining all labelled operators. It needs two dynamic-program passes,
not one pass per orbit. Source
`hd_bit_syndrome_orbit_gate.py` has SHA256
`1aad0e0c7b88397f18f5c81806ba04b4a2aab64015396e740af4155ea7636023`.
Seven additional small tests and independent review cover its class
selection and noncommuting products.

Three bounded full-state controls give the following results. X is the
input-activity marker; z is the output-weight marker. Each entry is a
log2 upper bound on the same complete r18/arity8/Fmax1 stratum, in FLOAT.

| Coefficient caps | z | X | Quiet component | Nonquiet component | Best of three Young choices |
|---|---:|---:|---:|---:|---:|
| Universal rank caps | 691/700 | 16 | +1764.38 | +9187.02 | +9190.96 |
| Four exact types | 691/700 | 16 | +1673.41 | +9136.58 | +9140.52 |
| Four exact types | 497/500 | 1 | -2401.36 | +885.01 | +888.95 |
| Four exact types | 993/1000 | 1 | -1355.44 | -330.17 | **-58.17** |

The first alternative avoids the large X^144 relaxation. Its remaining
nonquiet bound is dominated by paths with live state at every later column
boundary: the `Q01 * Q11^15` contribution agrees with the full nonquiet
bound to about 2e-10 bits. This motivated the single intermediate mark
993/1000, rather than a grid search. Fresh original one-state regional
curves were generated for both new marks; no post-refresh curves were used.
The intermediate control completes in 17.87 seconds.

A subsequent cheap recombination of the same cached column matrices uses
Young weights (1/16,15/16). It gives **-251.46236328725718**. This witness
is recorded separately from the three predeclared choices in the receipt;
it requires no new bank or dynamic program. Selecting among valid upper
bounds on the same event does not add a union factor.

An independent audit uses a sixteen-step row-vector recurrence instead
of the engine's matrix squaring and an 80-digit rational prefix. It
reproduces both values within 2e-11 bits and verifies all 23 direct/nested
source pins and all labelled row sets. This checks the saved FLOAT
arithmetic; it does not replace a directed replay. Forty-five focused
tests now pass across the outer, exact kernels, contraction engines,
loaders, and singleton cosets.

The decisive receipt is `hd_bit_syndrome_orbit_993_x0_float.json`, SHA256
`4d5d2374edca033663a18de558b1eef56751e9e4fc6b8081df5ec13c0f89e369`.
Its `.npz` branch cache has SHA256
`5d68a6899c5ba9ae88836a5347dc78f5a2ab10914286a6337ceda933c745531c`.
To reproduce the additional witness, call `combine_columns` on the cached
quiet/nonquiet matrices with groups18, arity8, log2X0 and weights1/16,15/16;
raise each result to16, apply initial zero and terminal ones, and take
`18*log2(sum_A 2^(word_log2[A]/18))`. Add the pinned global prefix
`-1602.1378501149884` once.

This is the first negative full-state comparison for this new local
kernel in the selected class. It does **not** establish the whole-code
35%/40-bit claim. Fresh directed arithmetic and the remaining source and
arity classes are still required. In particular, the arity-eight outside-
image identity cannot be assumed for other arities without checking their
cross-rank and fiber dimensions. The measured encoder is unchanged.

At X=1 the factor X^-1296 is one: the valid bound does not obtain a
numerical benefit from the input floor. It does not claim to extract
the exact polynomial tail above that floor.

The same Fmax1 regional envelope also covers source sets with F=0.
Their exact count is
`8192/(8192-18) * binomial(8192-18,18)` =
`4076644197931783287105561490932293125919092974583472128`.
Adding it to the F=1 count costs 4.776672304894 bits, without rebuilding
the operators. This extends the FLOAT proposal to homogeneous arity eight
and F<=1 only; sources with F>=2 remain outside it.

## Directed replay and coverage obligations

Freeze z=993/1000, X=1 and branch weights (1/16,15/16). A fresh directed
regional replay should use the original nonzero-amplitude packet formulas
and the positive 128-packet power, restricted to the downward-closed domain
`single+2*full<=18, full<=1`. A new r18 adapter must preserve the older
r17 helper rather than bypass its guards. Use exact route normalizers and
recompute or independently verify the exact source count.

Check the frozen regional moment slopes against every directed entry.
If needed, enlarge the intercepts outward without searching new slopes.
The 00 atom probability p must be one fixed number with complementary
mass exactly1-p; inflating both independently would change the empty
moment and its interpretation.

One economical replay can generate directed dyadic upper bounds for all
primitive atom and information-cover weights using Arb, then contract
the positive matrices with verified upward binary64 addition and
multiplication and explicit power-of-two scaling. Such a hybrid must
check its IEEE/subnormal assumptions, fail closed on overflow, validate
against small pure-Arb controls, and identify its arithmetic accurately.
Saved FLOAT column matrices cannot serve as certificate inputs. Finish
the per-row-set roots, Minkowski sum, and global prefix in directed
high precision, followed by an independent final sum.

Coverage must include heterogeneous constituent-row arities. With one
common output marker and positive path measure, let J_A be the complete
word functional of `K_A^r`, for a fixed row set A of size a. Its contribution
inside the group-Minkowski sum has factor

\[
 255^{-8a}X_a^{-9a}J_A^{1/r}.
\]

The exponents here are per group, not divided by r. Raising the complete
sum to r recovers 255^-1152 and X^-1296 in the homogeneous r18/a8 case.
Different homogeneous-arity contributions under the same measure combine
as `(sum_a U_a^(1/r))^r`, not `sum_a U_a`. Thus the next useful numerical
coverage test includes neighboring arities seven, eight, and nine under
one common measure. None of the older outer's receipts transfers merely
because an arity or source count matches.

Each neighboring arity has 11,440 labelled sets. For 11,200 of them,
inside and outside ranks are both five, so the existing fiber argument
extends after recomputing its coefficient caps. The remaining240 require
separate treatment. At arity seven, the actual outside image has rank
four and is not the entire outside affine-evaluation code. At arity nine,
the zero-outside fiber has dimension five although the quiet kernel has
dimension four; a nonquiet zero-outside branch must be retained. These
are proof obligations, not completed extensions.

## Directed selected-class replay

Fresh 256- and 384-bit regional replays now close the selected class.
Both use z=993/1000, X=1 and Young weights (1/16,15/16), with the original
one-state inner. Their reported upward log2 bounds agree:

| Source class, with r=18 and every arity eight | Log2 upper bound |
|---|---:|
| F=0 | -246.73930742323176 |
| F=1 | -251.46236328727593 |
| F<=1 | -246.68569098238194 |

The calculation includes every labelled row-set assignment, both outer
branches, and all boundary-state paths. It covers neither other r nor
mixed arities nor F>=2. Thus it is a directed selected-class bound, **not
a whole-code distance certificate**.

The regional adapter recomputes all 36 required (occupancy, full-pair)
matrices with Arb. A common exact dyadic atom measure then dominates all
76 entry/occupancy values by exact rational comparisons. Arb encloses
the primitive cover and terminal weights; exact rational checks validate
their conversion to scaled binary64 upper bounds. Each positive addition,
multiplication and binary scaling rounds upward. Explicit reductions
avoid unspecified floating-point summation, and the replay rejects
flush-to-zero arithmetic. The final roots, sum and prefix use Arb again.
No saved FLOAT curve or column matrix is a certificate input.

The two serial runs took 95.80 and 98.06 seconds. These are proof-replay
times, not encoder timings. The encoder and its measured 3.41-ms native
experiment are unchanged. The focused arithmetic/provider/driver suite
passes 37 tests, including exact rational noncommuting matrix controls,
subnormal cases, and a simulated flush-to-zero rejection.

The driver is `theory/hd_bit_syndrome_directed_gate.py`, SHA256
`0ea0597506e167bc9da77c31e6e2f5679d9cbfd8e7ee65d133e56453cce00086`.
The 256-bit receipt has SHA256
`1bfa4fa760caa1ea31f7288eb523350c52af607a6031ec63a172c940cb7b17e0`;
the 384-bit receipt has SHA256
`e3c73e930bbbcae54dc82533975b969f1eba628b3c7b4994f37d4a0bb2085d36`.
Each receipt pins its source dependencies and full labelled column cache.

Independent 512-bit final checks pass for both receipts. They import the
authenticated column bounds as exact dyadics and use sixteen successive
row-vector products, rather than the driver's binary64 matrix squaring.
Every cached word upper bound dominates the independently computed upper
bound. The rechecker also counts cyclic source sets independently and
recomputes the roots, sum, probability prefix and threshold comparison.
It does not independently rebuild the earlier column dynamic program.

The rechecker source is `theory/hd_bit_syndrome_directed_recheck.py`, SHA256
`224d1bb24dd5a90f4d141c68349472265f66ec1b2270a1d8a5d3c93a3780eadd`.
Its 256-bit-input recheck receipt has SHA256
`0230fa951722461a2b2e344cc7939eb765825aae12ca1ef9bdc50a88d5edeef0`;
the 384-bit-input receipt has SHA256
`62d1e901155a9cc78453d62e1fd70e6684a76562629f9a67ff972ab410ff2cb7`.
The broader outer/kernel/arithmetic regression passes 80 tests; eight
additional rechecker tests cover independent counts, chronology and
fail-closed validation. The [coverage plan](HD_BIT_SYNDROME_COVERAGE.md)
specifies the remaining mixed-arity and source partitions.

## A first mixed-arity extension

The next numerical gate includes arity nine. For a fixed nine-set A,
write O=A^c and C=H*[O,A]. Its image is im(V_O), since the nine source
labels span all five affine syndrome coordinates. Every outside value
has an inside affine fiber with difference space

\[
 D_A=H_*[A,A]\ker C.
\]

An information-set cover of D_A bounds every translated fiber. Multiply
that bound by the exact outside-image coefficient envelope, **including
its zero word**. For 11,200 labelled sets D_A is the quiet dimension-four
code. The other 240 use the actual dimension-five fiber. The code computes
their information covers individually; it does not transport the quiet
matroid to these exceptions. All 11,440 labelled terms keep their order
through the sixteen-column product. The resulting outside-marked dynamic
program has seven marked positions and 330 count classes.

At the previous z=993/1000, this full-fiber bound fails: its homogeneous
arity-nine log2 term is +81.45930226. The ordinary rank-five sets dominate
(+81.08261110), not the 240 exceptions (-28.65908712). The path that
activates in the first column and remains live contributes +81.45930185;
the all-zero-boundary component is -406.42621758. These are components of
the numerical majorant, not demonstrations of low-weight codewords.

A single stronger marker, z=124/125, clears arity nine at -359.63 but
makes the arity-eight bound fail. The subsequent interior marker
z=9927/10000 supplies **one fresh common measure** that clears both:

| r=18, F<=1, X=1 | FLOAT log2 term |
|---|---:|
| Every group has arity eight; Young weights (1/8,7/8) | -202.69164754 |
| Every group has arity nine; full-fiber bound | -347.49712895 |
| Each group has either arity, including all mixtures | -202.59349190 |

The mixed calculation applies the per-group RS factors inside the
eighteenth-root sum, then takes its eighteenth power. It is not a sum
of the two homogeneous bounds. The source count, output marker, and
projective normalization occur once outside that expression. The run
took 22.30 seconds and changed neither construction nor implementation.

An independent calculation with successive row-vector products and an
80-digit rational prefix agrees within 3e-11 bits. It also checks both
complete labelled lists, the source pins, and all 76 moment inequalities.
This remains **FLOAT mixed-class evidence**, not a directed bound or a
whole-code certificate. Other arities, F>=2, and other r are still open.

The arity-nine kernel source is `theory/hd_bit_syndrome_arity9_gate.py`,
SHA256 `a02369609ee860cec6fd236e952e82e841ef776cb369c4735a88fbb93d186e48`;
eight small tests and independent source review pass. The
combined outer/kernel/arithmetic/rechecker/arity-nine suite passes 96 tests.
The successful control is `theory/hd_bit_syndrome_mixed_control.py`, SHA256
`e88e86cef3f0aa7f2ff616e2d3ec129899dce4653d1ec3e63502ea1cc9c870d8`.
Its receipt SHA256 is
`d94ccf2b9a9bc073b67c90256608fa67ef7be02f44efc47eede301f627c4e2c3`;
its cache SHA256 is
`0cdf0b60933ae1696c434afb4986caa9bf29f986a2f520762ba5245b1f0f6067`.
The failed 124/125 driver is preserved as
`theory/hd_bit_syndrome_mixed_control_992_snapshot.py` with its receipts.

The subsequent arity-seven test fails under this same measure at
+124.11946657 bits. Its zero-boundary paths dominate; the first-activation,
then-all-live path is already small. The repair below uses a constraint
that this relaxation had discarded, without changing the encoder.

## Enforce the RS column activity before combining groups

Every active group contains a nonzero RS[16,8] row. That row has at least
nine nonzero columns. A group cannot therefore use the zero-input term
in more than seven of its sixteen column kernels.

For fixed A, let K_A be its true all-amplitude column polynomial; its
zero-input term is one. Mark every nonzero column by Y>=1. The marked
kernel is K_Y=1+Y(K_A-1), and each whole-group bound acquires Y^-9.
This counts nonzero columns, not the total number of nonzero input symbols.

Let r be the number of active groups and let the common positive column
functional be L. Suppose Q_A=L[B_A^r], where B_A>=K_A>=1. Write S=L[1].
The scalar inequality

\[
 (1+Y(K-1))^r\le 1+Y^r(K^r-1),\qquad K,Y\ge1,
\]

gives the valid matrix majorant

\[
 Q_A(Y)=S+Y^r(Q_A-S).
\]

All entries use the same path functional. In the implementation,
S is the sixteenth power of the atom measure's constant matrix. A weaker
version subtracts only its common all-zero-atom component. After the
sixteen-column product and the r-th root, include Y^-9 in each group's
normalized term. Different arities may use fixed different Y values.

The cached r18/F<=1 control keeps z=9927/10000 and the a8/a9 bounds fixed.
Using Y=5/2 for arity seven gives:

| Selected class | FLOAT log2 bound |
|---|---:|
| All groups have arity seven, previous relaxation | +124.11946657 |
| Arity seven with the activity constraint | -71.26997479 |
| Every mixture of arities seven, eight, and nine | -71.10522351 |

This control uses no new regional coefficients, witness fit, or native
change. Four small tests pass. An independent recomputation forms S with
Arb384 arithmetic, checks every residual, and uses successive row-vector
products for the whole words. It agrees within 3.1e-11 bits. The result
remains FLOAT evidence because the input column matrices are FLOAT bounds.
Other arities, source weights, and F>=2 are not covered.

The driver is `theory/hd_bit_syndrome_activity_control.py`, SHA256
`533fe489bacad9af0505dc1d17bef94e4a0d0a7807170ad33a51e036417346b3`.
The receipt `hd_bit_syndrome_activity_9927_float.json` has SHA256
`8394caf34672aa43dd92779ab9511578953c431382302b844196340d59fbd8b3`.
For directed replay, accumulate Q_A-S as a positive excess rather than
subtracting two independently rounded upper bounds. Some cached FLOAT
entries lose roughly forty bits through this subtraction.

## Exact empty-step refresh without a new state reference

The previous regional proof also pays a refresh density penalty at every
empty inner packet. The inner refresh is D2 H16 D1, where H16 is the
superregular interpolation map. It is not the outer matrix H*.

For each nonzero input weight a, this refresh has an exact output-weight
distribution P[a,b]. Fresh nonzero diagonal entries make its output
uniform within each weight class. The same statement holds for inverse
refresh, since H16^-1 is superregular. Let e_a(z) be the exact average of
z^wt(AQ) over state vectors Q of weight a. Here A is the actual direct sum
of four RS[16,4] emission maps. Its selected coefficient columns are not
assumed to form MDS codes; their ranks are computed explicitly.

The uniform-nonzero live reference therefore admits the empty-packet cap

\[
 h_1(z)=\max_{1\le b\le16}\sum_{a=1}^{16}P[b,a]e_a(z).
\]

This bounds the output density pointwise, not only its total mass. Thus
the existing two-state regional proof can use h_1 in its empty11 entry.
No reference-distribution change or adjustment of nonempty entries is
needed. This is a proof improvement for the unchanged native encoder.

At z=9927/10000, exact rational calculation followed by Arb256 and Arb384
logarithms gives an empty-region ceiling of -86.252430843603 bits, versus
-74.68826991984223 previously. The improvement is at least
11.564160923761 bits for a completely empty region. This number cannot
simply be multiplied by the number of regions in a distance claim;
nonempty packets and all state paths must still be included.

Five exact small-field tests check every fixed input and target,
emission counts, weighted histories, and the forward/inverse density
identity. The source `theory/hd_mds_empty_weight.py` has SHA256
`6c76d7513942d67d10ebe7b5610517525fdb29f5de42e416db474b92b49e3920`.
Its scalar receipt `hd_mds_empty_weight_9927_exact.json` has SHA256
`83df60adddbd2979372db1a9f74ace1542ce6febbd5dfdb318c8148632014b5d`.
The isolated regional adapter changes only this empty11 entry. A fresh
same-mark common measure then gives these selected-class FLOAT bounds:

| r18, F<=1 | Without column marker | With Y_7=5/2 |
|---|---:|---:|
| Arity seven | +120.61114880 | -76.56842969 |
| Arity eight | -206.08246122 | unchanged |
| Arity nine | -1238.98832182 | unchanged |
| Every mixture of the three arities | +120.61123810 | -76.39182947 |

The two repairs address different losses. The activity marker removes
impossible sparse column words. Exact empty-step refresh reduces the
charge for long live-state runs, which dominated the arity-nine bound.
Neither changes native work. A second run with the complete dependency
pin set reproduces every reported bound exactly and takes 33.51 seconds.

The authoritative receipt is
`theory/hd_mds_exact_empty_mixed789_9927_pinned_float.json`, SHA256
`c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8`;
its column cache has SHA256
`f5a6e4dba34d989acb59cdc24fe5cb822d351983ee50ff44325d31cc3000485d`.
The driver SHA256 is
`33cb79621068a638408076de1dd47f9666f562f031b7d449e4900d3ca5a3726f`.
The first receipt without the emission-helper pin is preserved, but the
fully pinned second run is the reference. The numerical scope remains
partial and FLOAT; a directed whole-code certificate is still open.

An independent cached recheck verifies all 31 source pins, all 35,750
labelled row sets, the exact empty-cap receipt, and every old/new regional
comparison. Arb384 reconstruction of the constant functional and activity
residuals, followed by separate chronological products, reproduces the
mixed bound within 3.1e-11 bits. One FLOAT moment comparison is short by
about 1.78e-15 bits when its stored logarithms are treated as exact dyadics.
This is recorded, not silently accepted as directed arithmetic; a directed
replay must round the measure outward. The independent receipt is
`theory/hd_mds_exact_empty_mixed789_9927_independent_recheck.json`.
The combined old/new local, kernel, arithmetic, and rechecker suite passes
136 tests. The certified fallback manifest is unchanged.

A separate [rank-eight mixer candidate](HD_BIT_SYNDROME_STAR_OUTER.md)
uses ten additional XORs and no additional field products. Its exact
local algebra and pair laws are checked, but it has not been timed or
certified. Preserve it as a construction alternative; first exploit
these proof improvements on the already measured H* encoder.

## Exact singleton-outside fibers, reserved for tightening

The projection bound loses the translation of each inside fiber. A new
exact helper retains it when the outside output has weight one. Set that
outside value to one; scaling through the 255 nonzero values preserves all
supports and contributes that factor once.

With `c=B V_A^T x`, solve `V_O c=e_j`. The inside output satisfies

\[
 V_A^T y_A=d:=B^{-1}c+V_A^TV_Ac.
\]

For S contained in A, its supported count is zero unless d lies in the
field span of V_S. Otherwise it is `256^(|S|-rank(V_S))`. Exact subset
inversion gives every nonzero support count. This retains the actual
labelled H* map rather than assuming affine symmetry of its coshift.

All 1,920 feasible labelled fibers were checked. Of these, 1,888 have
inside minimum weight five; 32 exceptional sets have minimum four. Thus
the total output weight is at least five universally and at least six
outside those exceptions. Each fixed outside value has exactly256^3
preimages. These claims are locally exact, but this tighter kernel was
**not** used in the negative full-state result above.

`hd_bit_syndrome_singleton_cosets.py` has SHA256
`18a3ee622869a7f0450785deb5e1598b132c5fed1947673c3187dc44586c04e8`.
Six tests include the complete finite audit and independent direct-H
rank checks. This helper is available if other source classes require a
sharper nonquiet bound; it is not a new construction or a distance
certificate by itself.
