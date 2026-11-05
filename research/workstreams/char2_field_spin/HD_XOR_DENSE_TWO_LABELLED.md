# Two dense groups: a labelled occupancy control

Retaining labelled dense occupancies makes the two-group interface practical
and recovers 539.85 bits from the coarse feature envelope. The one authorized
two-group control still fails: its FLOAT log2 bound is `+316.244166693742`.
This is not a directed certificate or a failure example for the construction.

The scope, comparison measure, and counting factors are those in
[the feature-control note](HD_XOR_DENSE_TWO_FEATURE.md). Fix source matching
counts `(8,0)` among 17 active groups. The target cohort is `(0,15,2)`.
Both dense groups use the same deterministic type: the 36 arity-five row
sets assigned to the frozen `a4_mask240_fit` witness. Each row set remains
fixed across 32 columns. The sparse pairs likewise remain fixed; all
`28^15` choices are counted once. No old helper, literal, or receipt changed.

## The retained column information

For the two dense output masks `T_1,T_2`, define
`ell_i=1_{i in T_1}+1_{i in T_2}` and let `e` count nonempty masks. Each
labelled occupancy vector `ell` lies in `{0,1,2}^8`. Sum the positive weights
`W_B(T_1)W_B(T_2)` into `M[e,ell]`. There are at most 6,561 occupancy vectors
and only 65,536 mask pairs in this one-time convolution. It preserves the
whole-mask dependence used by the Hölder argument. It does not replace
the weights by independent coordinate marginals.

For each sparse occupancy interval `[l,u]=[d-k,s+d]`, the new helper computes
all labelled offsets in one suffix DP. The offset at the newly prepended
region is a base-three digit. Suffixes with the same remaining digits and
resource sum share their vector. Matrix multiplication still follows the
original labelled-region chronology.

Only final shifted sums `7s+2k` are required. Writing `w=u-l`, these are
`7w-5k`, `0<=k<=w`. A partial sum is discarded precisely when no required
final sum is reachable using the remaining coordinates in `[0,w]`. If a
path reaches a retained target, every predecessor on that path is therefore
retained. Pruning changes the work, not the resource envelope.

For each target, the helper sums `M[e,ell]` times its complete offset vector.
It then accumulates the exact binomial cancellation weights. Both positive
mixtures finish before the common-potential quotient. The same global
mark and potential apply to both dense groups and to every column.

The remaining computation is unchanged: take the positive tail `X>=510`,
`Y>=34` of the 32nd power of the scalar column polynomial. Multiply by
the source-support count, `C(17,2)=136`, and `28^15`, and divide by the
single projective factor `255`, endpoint potential, and sparse normalizer
`255^480`. The dense row-set counts and comparison factors are already
inside `W_B`; neither is multiplied again. The aggregate `Y` floor safely
relaxes the two individual floors.

## Checks and one bounded production control

Six exact/light tests pass. They cover labelled mask-pair convolution,
exhaustive reachable sums, noncommuting full-resource comparison for all
81 tiny ternary offsets, equality with the existing one-group recurrence,
complete mask/cancellation mixture, and timeout metadata without a bound.

Before the two-group run, the same implementation reproduced the frozen
one-group control. Both controls used exactly
`z=8973939342378351/9007199254740992` and
`log2(p)=-2386334719/33554432`.

| Control | FLOAT log2 bound | Elapsed | Dominant `(X,Y)` |
|---|---:|---:|---:|
| Labelled one-group reproduction | -56.177151934928 | 2.52 s | `(1019,32)` |
| Frozen one-group comparator | -56.177151934930 | existing artifact | — |
| Labelled two-group target | +316.244166693742 | 17.02 s | `(957,64)` |
| Prior coarse two-group control | +856.091648792149 | 1.80 s | `(938,64)` |

The one-group difference is `1.96e-12` bits in FLOAT arithmetic. The two-group
operator phase took 16.01 seconds. Its measured structural counts match the
prior calculation: 136 intervals and 150,969,492 two-state matrix candidates.
The largest single DP array contains 398,034 doubles, about 3.04 MiB. This
is an array-size count, not a process peak-memory measurement.

The executable enforces a 60-second compute budget, checks it during each
interval, and writes a failure record without a bound if exceeded. This
run completed within the budget. No second two-group witness, fitting,
type-histogram sweep, directed replay, or native benchmark was performed.

The dominant coefficient has `Y=64`, the maximum possible dense degree;
both dense groups are nonempty in every column. Its sparse degree `957`
is also near the maximum `960`. Tightening only the aggregate dense-column
floor therefore cannot remove this coefficient. The result does not
separate inherited-witness loss from the remaining resource-envelope loss.

The next bounded step should test a different *global* frozen mark/potential
for this same type pair before fitting or expanding to 78 type pairs. The
current interface makes that control affordable. No independent per-group
choice of the mark or potential is justified by this calculation.

## Reproduction and pins

From `theory`, use fresh output names:

```powershell
C:/Python314/python.exe -m unittest hd_xor_dense_two_labelled_tests -v
C:/Python314/python.exe hd_xor_dense_two_labelled.py --dense-groups 1 --seconds 60 --output fresh_labelled_h1.json
C:/Python314/python.exe hd_xor_dense_two_labelled.py --dense-groups 2 --seconds 60 --output fresh_labelled_h2.json
```

| File in `theory/` | SHA-256 |
|---|---|
| `hd_xor_dense_two_labelled.py` | `1a4e2d7bce9825d44cbd9b60ca71abe65223b1de1df364af7c3c721a1828fd66` |
| `hd_xor_dense_two_labelled_tests.py` | `78f89db7c27473de1d5f530d563f17150308933af1fe83617f2a6caf42ff6a7c` |
| `hd_xor_dense_two_labelled_h1_olda4.json` | `99459209bafe1232374241490d16872b7c5fe86df35fed825847714f053c55a9` |
| `hd_xor_dense_two_labelled_h2_olda4.json` | `aca631409822486aebeed5a04bd5708e12de23e6dbeba173238b0dd126685bc1` |

Hashes identify the saved FLOAT experiment. They do not validate interval
arithmetic, and no interval-arithmetic claim is made for these files.
