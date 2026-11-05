# Fixed-source matching strata for the G16 FLOAT gate

The previous all-cohort control maximized each regional matrix entry over every feasible full-pair count. That relaxation discarded a constraint shared by all columns: every regional active set is a subset of the same source set.

This control retains that constraint without changing the encoder or its random setup. Its scope is the existing G16 affine outer construction with the paired, independently swapped route. It evaluates only the three frozen anchors below. All numerical output is a FLOAT proposal, not a directed certificate.

| Active groups r | Weight mark z | Fractional exponent |
|---:|---:|---:|
| 18 | 691/700 | 1/5 |
| 24 | 172/175 | 1/5 |
| 32 | 171/175 | 1/5 |

## Source partition and regional envelope

Fix an indexed source set G of size r in the cycle of 8192 groups. Let a(G) and b(G) count its full pairs in the even and odd adjacent matchings. Define F(G)=max(a(G),b(G)). The existing exact helper `matching_max_counts(8192,r)` returns

\[
N(r,F)=\#\{G:|G|=r,\ F(G)=F\}.
\]

These counts partition the indexed source sets; their sum is binomial(8192,r). There is no quotient by rotations.

For either matching, every regional active subset U of G satisfies

\[
0\le f(U)\le \min(F,\lfloor |U|/2\rfloor).
\]

Let R[m,f] denote the existing paired-NZ regional matrix after its prescribed fractional step and route average. For each stratum, the control uses the entrywise envelope

\[
\overline R_F(m)=\max_{0\le f\le\min(F,\lfloor m/2\rfloor)}R[m,f].
\]

The maximum is valid for every fixed source in the stratum. It does not sample a new source between regions or columns. F is only an upper bound for either matching, so the control imposes no lower bound derived from treating F as an exact count.

The control applies the existing state-preserving empty-cut bound to these restricted curves. Its original numeric slopes remain fixed; it recomputes intercepts over every occupancy 0 through r. Each of the eight inherited type templates stays fixed across all state paths, cuts, and sixteen columns. In particular, constituent row sets remain fixed when Hölder and Minkowski are applied.

For a template w, let B(r,F,w) be the resulting bound for one source set, including the fractional weight and projective factors. The source union is

\[
\sum_{F=0}^{\lfloor r/2\rfloor}N(r,F)\min_w B(r,F,w).
\]

The exact source count enters once, outside the fractional exponent. The template may depend on the deterministic class (r,F), but not on sampled setup. Selecting one bound per class introduces no extra template factor.

## Bounded implementation and status

`theory/hd_two_sum_gf16_source_strata.py` rebuilds three paired-NZ banks through occupancy32 and retains every R[m,f]. It checks each unrestricted envelope against the pinned bank. It then evaluates all40 source strata with eight cached type templates each:320 envelope evaluations, without fitting. A30-second guard leaves incomplete output explicitly uncompleted.

The first control has no total-arity restriction. It covers every cohort and every constituent row-set assignment within each requested source stratum, but not intervening active-group counts.

Seven exact-small and mocked orchestration tests passed in0.061 seconds. They check the indexed source partition, every subset of small cycle supports, the upper-only f constraint, once-only source normalization, and all320 unchanged-template calls. No production bank or gate was run as part of those tests. Root and an independent agent reviewed the implementation without finding an accounting or scope defect.

The authorized control completed all40 strata and320 evaluations in14.16149 seconds. It rebuilt exactly three regional banks. Every anchor still fails:

| r | Complete source-stratum union, log2 bound | Earlier unrestricted-source envelope |
|---:|---:|---:|
| 18 | +2697.8809382919812 | +2700.2512333446534 |
| 24 | +3785.5576260306325 | +3789.765729491666 |
| 32 | +6035.832016854997 | +6039.627029653657 |

These positive values are uninformative probability bounds, not evidence against the construction. The control imposes no low-total-arity cutoff.

At r24 and r32, F1 dominates the union. Its selected template is the inherited all-zero floor template. Its dominant regional state path alternates01 and10 eight times, and its largest type contribution has arity8. Both off-diagonal affine intercepts attain their maximum at occupancy m2 and full-pair count f1. For F0, those contacts occur at m1,f0. Thus one allowed full pair restores the small-occupancy intercept; the exponential envelope then extends that intercept to all occupancies.

The F0 contributions are +3684.53147 at r24 and +5774.19864 at r32. Even this source class remains far above the target. Fixed-source geometry alone therefore does not explain the remaining loss.

The saved artifact retains every fresh R[m,f], all320 template results, the exact source counts, and the selected union. Its name is `theory/hd_two_sum_gf16_source_strata_anchors_float.json`; SHA256 is `5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf`.

Source SHA256: `253f5bbb9d1fa1326cbeb63ad7018fce2fb784b66caac4e57374b2cfc6640b81`.

Test SHA256: `57ee01b1cc0eaed525d408b4a529229f45a1119a3ed4622313a0acf78b08514f`.

Reproduce from the theory directory after the quantitative slot is granted:

```powershell
C:/Python314/python.exe -m unittest hd_two_sum_gf16_source_strata_tests -v
C:/Python314/python.exe hd_two_sum_gf16_source_strata.py --max-seconds 30 --output hd_two_sum_gf16_source_strata_anchors_float.json
```

The next structural control should retain the off-diagonal matrices' exact zero at m0 and distinguish their small-occupancy contacts from the remaining curve. The saved banks suffice to diagnose such an envelope before rebuilding any provider. More F bins or another mark grid would not address the identified relaxation. Any new aggregation rule needs its own fixed-row-set argument and exact small tests.
