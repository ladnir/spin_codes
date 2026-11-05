# G16 low-total-arity anchor control

The high-total-arity uniform comparison leaves a complementary class with
fewer active constituent rows. This note records a bounded FLOAT control
for that complement. It does not certify the G16 construction or any
continuous interval of active-group weights.

The construction is the existing G16 affine candidate: sixteen RS[16,8]
constituents per group, labels `06 + GF16`, independently swapped paired
routing, and the unchanged width64/state16 inner map. The probability
space is the ideal independent setup ensemble, not one native PRNG seed.
The failure event is a nonzero scalar word of output weight at most 734003.

## The restricted event and its marker

Fix an active-group set of size `r`. Each group has a fixed nonempty
constituent-row set; let `a_g` be its cardinality and
`A_tot=sum_g a_g`. The three controls use:

| `r` | High-side cutoff `A_0` | Low-side maximum `B=A_0-1` | `z` | `alpha` |
| ---: | ---: | ---: | ---: | ---: |
| 18 | 247 | 246 | `691/700` | `1/5` |
| 24 | 317 | 316 | `172/175` | `1/5` |
| 32 | 463 | 462 | `171/175` | `1/5` |

For a deterministic marker `0<t<=1`, the inequality
`1_(A_tot<=B)<=t^(A_tot-B)` restricts the positive type union. Each type-`a`
whole-word factor acquires `t^a`; distribute this as `t^(a/16)` over its
sixteen column factors. Multiply the final bound by `t^-B` once.
These factors are outside `alpha`. The type choices are fixed across all
columns and are not resampled after the marker is introduced.

The underlying fractional packet provider is unchanged. It conditions on
the route, averages outer amplitudes and inner setup, and applies the
fractional bound before averaging route placement. Its packet cap entries
are raised to `alpha` before the geometry average. The route denominator
is unpowered. Dense weights are the fresh moments

\[
 W_{a,\alpha}(T)^{16}
 =255^{-8a\alpha}\sum_{|A|=a}N_A(T)^{16\alpha},
\]

not powers of the alpha-one weights. The derivation is recorded in
`HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md`.

For a chronological path `gamma` through the sixteen regional states,
let `A_gamma` be its regional-envelope intercept product. Let
`F_a(gamma)` be its type-`a` column mass, including the fixed support-floor
marker and once-only row-set count distributed over sixteen columns.
The existing fixed-type Hölder and path Minkowski bound gives

\[
 K_{r,t}(i,j)=\sum_{\gamma:i\to j}
 A_\gamma\left(\sum_{a=1}^{16}F_a(\gamma)^{16}t^a\right)^{r/16}.
\]

Therefore the low-side bound at the specified weight is

\[
 U_{r,\le B}\le
 \binom{8192}{r}(z^{-734003}/255)^\alpha\,
 t^{-B}e_0K_{r,t}^{16}\mathbf1.
\]

The source-support count occurs once, outside `alpha`. There is no further
cohort or row-set factor and no endpoint-potential factor. The packet
matrices remain ordered; grouping paths by their four transition counts
uses the common all-matching envelope and symmetric type mask laws.

## Finite control and result

`theory/hd_two_sum_gf16_low_arity.py` builds exactly three regional banks
through occupancy 32. Each anchor uses only its feasible prefix `m<=r`
when checking every envelope intercept. Initial slopes join each entry's
nonempty peak to its tail. Their bounded candidate sets include zero and
nonempty secants, avoiding the earlier deep-empty-endpoint initialization.

The fixed total-arity logs are `{0,-1,-4,-16,-64}`. Three floor templates
are used at each actual anchor: all-zero and the previous per-type rule
at scales one and two. All markers remain fixed across every path and
column. The approved run evaluated these 45 controls without the optional
worst-anchor refinement.

The run completed in 1.426 seconds. All controls failed:

| Anchor | Best FLOAT log2 bound | Selected marker and floor |
| ---: | ---: | --- |
| 18 | `+6209.3385828187365` | `t=1`, all-zero floor |
| 24 | `+7998.822804726728` | `t=1`, all-zero floor |
| 32 | `+10360.85023055943` | `t=1`, all-zero floor |

The saved-envelope diagnosis differs from the previous live-state defect.
The largest `K(0,0)` path term stays quiet through all sixteen regions.
Its type norm peaks at arities nine, nine and eight at the three anchors.
Their total arities are 162, 216 and 256, below the respective low-side
cutoffs. A total-arity marker cannot remove those contributions.

The all-empty-mask comparison term alone is approximately `+2762,+3683,
+4912`. Anchor-scale-one floors reduce that term below `-60`, but increase
the complete bounds to approximately `+7376,+9599,+12549`. Thus this
control still loses on positive quiet-state dense mask mass and its
regional envelope, not only on admitted empty words. These are failures
of the comparison, not lower bounds on actual setup failure.

Eight lightweight tests check the exact marker inequality, its unscaled
placement, a noncommuting path example, peak/tail initialization, source
scope, per-anchor floors and the 45-control resource limit. Root and an
independent read-only reviewer cleared the source before the run.

The artifact `hd_two_sum_gf16_low_arity_anchors_float.json` has SHA-256
`206f34cca8021de88ce8806a05c6736523261bf769ff76f147542e7cb6fdc5d5`.
The source has SHA-256
`943dd058318f0f4c5d5a371fa3a84210646e71ca8c01327399d9a60ce7b9064a`.

Reproduce only with a fresh output filename:

```text
python -B -m unittest hd_two_sum_gf16_low_arity_tests
python -B hd_two_sum_gf16_low_arity.py --max-seconds 30 --output NEW_LOW.json
```

## Paired-packet comparison at identical witnesses

The separate provider `theory/hd_two_sum_gf16_paired_nz.py` composes two
consecutive good packet caps before their entrywise fractional power.
Here a good packet has at most one singleton source pair. Its cap is
invariant under that singleton's fair orientation, so the existing
`good_geometry` multiplicities need no additional `2^-s` factor.

For each ordered good geometry pair, the new coefficient uses
`(M_g M_h)^alpha` with the geometry counts outside the power. The other
cases retain `GD+DG+DD` from packetwise fractional caps. The complete
two-packet polynomial is raised to power 64 and divided by the unchanged,
unpowered regional support count. No packet position is discarded or
commuted. Subadditivity implies that this provider is no worse than the
packetwise provider in exact arithmetic; alpha one gives equality.

Seven explicit-table tests check that comparison, alpha-one regional
equality, exact geometry counts and ordered noncommuting products. Root
and two independent reviewers cleared the provider before production use.

`theory/hd_two_sum_gf16_paired_low_control.py` then rebuilt three banks and
replayed all 45 low-side witnesses. The A/B preserved every numeric slope,
mark, fractional exponent, type-floor log and total-arity marker. Only the
provider changed; intercept maxima were recomputed over the same feasible
occupancies. The four driver tests enforce that identity and reject a
provider regression. The complete combined suite has 19 passing tests.

This A/B completed in 1.606 seconds. Every regional entry passed the
FLOAT no-worse check, but all three selected bounds were unchanged at
binary64 precision: `+6209.3385828187365`, `+7998.822804726728` and
`+10360.85023055943`. The quiet-path envelopes still dominate. No refinement
or new witness search followed the run.

The new artifact `hd_two_sum_gf16_paired_low_anchors_float.json` has SHA-256
`10dd04490a95f124767925c44954c0b375ebcc6f3d25506e52f898452c29761b`.
The provider SHA-256 is
`e6938a76e750924d03329e091a21582d5686b606e41517d43da6118c3995e048`;
the driver SHA-256 is
`7ecee98e466729c37f5c8eb8e8ffc53075ee9238f3a6ab4ebdbb8cb62c9bd9a4`.

The next structural control should preserve the empty versus nonempty
quiet-state distinction. An exact cumulative arity polynomial can
separately measure total-arity Chernoff loss, but cannot exclude the
already dominant lower-arity terms merely by enforcing these cutoffs.
These possible refinements need independent review and directed replay
before any successful proposal could become a certificate.
