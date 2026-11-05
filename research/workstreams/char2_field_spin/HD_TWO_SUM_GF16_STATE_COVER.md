# Exact empty off-diagonal constraints: bounded FLOAT control

The positive coverage recursion is valid, but it does not close the three
tested G16 classes. Restoring `R01(0)=R10(0)=0` improves the frozen
`F=1` bounds by at most 0.2641 bits. The remaining positive bounds are
thousands of bits above the target. This is a failed proof control, not a
counterexample to the encoder and not a directed certificate.

The construction remains sixteen RS[16,8] rows, the affine two-sum mixer,
independently fair paired routes, and the unchanged h4 inner transform.
The control uses only active-group weights 18, 24, and 32. Here
`F=max(a(G),b(G))=1` is fixed for the indexed source set `G`; it is not
resampled across regions or columns. The retained source-stratum tables
come from [the source-stratum helper](theory/hd_two_sum_gf16_source_strata.py).
Every regional subset has at most one full adjacent pair in either matching.
The exact count `N(r,F=1)` is charged once, outside the fractional exponent.

## Positive coverage before the fixed-type norm

Fix a chronological inner-state path and one of its positive empty-cut
branches. A coordinate assigned a `00` delta branch must be absent from
every group mask. Coordinates with transitions `01` or `10` must instead
be covered by at least one group mask. The remaining coordinates are optional.

Let `n1,n2` be the numbers of required `01,10` coordinates. They satisfy
`n1+n2<=16`. After the common symmetric regional and type replacements,
the coverage state is just `c=(c1,c2)`, the number already covered in each
class. There are at most 81 states. This symmetry concerns the comparison
polynomial, not a new symmetry assumption about a fixed physical source.

For fixed type `a`, let `w_a(t)` be the mass of one labelled mask of weight
`t`, including its fixed row-set and floor factors. There is no binomial
mask-count factor in `w_a`. Let `B1,B2` be the regional coordinate markers
for required positions, and let `e_u` be the elementary coefficients of
the optional markers. Define

\[
 L_a(h_1,h_2)=B_1^{h_1}B_2^{h_2}
       \sum_u e_u w_a(h_1+h_2+u).
\]

The positive transition from `c` to `c+d` is

\[
 T_a(c,c+d)=
 \binom{n_1-c_1}{d_1}\binom{n_2-c_2}{d_2}
 \sum_{h_1=0}^{c_1}\sum_{h_2=0}^{c_2}
 \binom{c_1}{h_1}\binom{c_2}{h_2}
 L_a(h_1+d_1,h_2+d_2).
\]

The first binomials select newly covered coordinates. The inner sum selects
already covered coordinates again. Every group mask appears once, including
its optional-coordinate choices. For a fixed ordered type vector
`tau=(tau_1,...,tau_r)`, its covered column mass is exactly

\[
 C_\tau=(T_{\tau_1}\cdots T_{\tau_r})_{0,(n_1,n_2)}.
\]

Set `p=16` and `U(c,d)=(sum_a T_a(c,d)^p)^(1/p)`. Minkowski over coverage
histories, followed by factorization over the group indices, gives

\[
 \left(\sum_\tau C_\tau^p\right)^{1/p}
 \le (U^r)_{0,(n_1,n_2)}.
\]

The type vector remains fixed across the 16 outer columns in the preceding
Hölder step. This inequality does not sample fresh types. It uses no signed
inclusion–exclusion and never subtracts upper bounds. If no coordinate is
required, there is one coverage history and the old unrestricted type norm
is recovered exactly. Otherwise, the helper takes the smaller of the old
norm and the coverage-history bound **for each path**, before the path sum
and final 16-column matrix power. The additional Minkowski step can lose
information, so this minimum is useful even though the zero constraints
themselves are exact.

## Frozen six-control result

The regional tables, slopes, weight marks, and type-marker literals were
reused unchanged. The exponent is `alpha=1/5`. Each saved old `F=1` result
was reproduced before applying coverage. No regional provider or optimizer
was called. The 30-second control completed in 4.9735 seconds.

| r | z | type-marker template | failure log2 FLOAT screen | improvement |
|---:|---|---|---:|---:|
| 18 | 691/700 | old zero | +2757.53284290728 | 0.264062329868 |
| 18 | 691/700 | joint pre-H empty-floor1 | +3787.46868922621 | 0 |
| 24 | 172/175 | old zero | +3785.53244939005 | 0.000008476702 |
| 24 | 172/175 | joint pre-H empty-floor1 | +5448.49339048183 | 0 |
| 32 | 171/175 | old zero | +6035.76821285539 | 0 |
| 32 | 171/175 | joint pre-H empty-floor1 | +8041.83180448566 | 0 |

Each old-zero control improved eight path/cut classes, with maximum local
norm gains of 7.4063, 7.0711, and 6.7227 bits. No joint-floor path improved
materially. The all-`00`, all-forbidden path is unchanged by definition:
it has no required off-diagonal coordinate. The full bound shows no material
improvement at weights 24 and 32. Repaired path rankings were not retained,
so the artifact does not establish exact equality of dominant-path indices.

This result rules out omitted empty off-diagonal constraints as the source
of the missing thousands of bits in this particular comparison. The next
proof question was how to retain whole-group nonzero and RS-row constraints
before the mask/type relaxation. The subsequent
[group-activity control](HD_TWO_SUM_GF16_GROUP_SEEN.md) tests nonzero groups
and the aggregate nine-column floor separately. It improves the first
anchor but not the larger two; it does not enforce each group's complete
nine-column history. Further witness grids are not justified by these
controls. The encoder and all existing certificates remain unchanged.

## Reproduction and tests

From the theory directory:

```text
python -B -m unittest hd_two_sum_gf16_state_cover_tests hd_two_sum_gf16_state_cover_control_tests
python -B hd_two_sum_gf16_state_cover_control.py --max-seconds 30 --output NEW_OUTPUT.json
```

The eight helper tests include explicit labelled-mask enumeration, the
fixed-type norm inequality, exact zero propagation, unique-history equality,
and an underflow-safe log fallback. Eight driver tests cover the source-count
replacement, fixed literals, six-control limit, baseline check, scope, and
partial timeout. The 43 related tests passed in 0.250 seconds. Independent
read-only reviews cleared both the mathematical interface and the code.

Frozen SHA-256 pins:

```text
hd_two_sum_gf16_state_cover.py
789ab7ab0ad1f4ea4b3e2874a2c7393ff00aac4f016a77bdeada3cc5c8dfbda5
hd_two_sum_gf16_state_cover_control.py
0549c92786c5d2ea7a635014991d279ab3a73dd0490bea2f7809acd491739a06
hd_two_sum_gf16_source_strata_anchors_float.json
5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf
hd_two_sum_gf16_state_cover_f1_anchors_float.json
e26ae2832a871933853d2399ed80b6e4077c67f6214c481515702493dfbea118
```

The artifact is scoped to three `F=1` strata. It is neither a complete
source partition nor a continuous weight-band result.
