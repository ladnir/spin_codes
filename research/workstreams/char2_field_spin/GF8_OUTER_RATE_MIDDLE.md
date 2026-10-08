# Outer-dimension sweep: middle-weight classes

With the existing sixteen-symbol, eight-state inner held fixed, increasing
the GRS outer dimension from 87 to 88 remains compatible with the middle
proof at K18 and K20. The next dimension, 89, has a substantial unresolved
middle-weight interval. These are statements about the available bound,
not upper bounds on the realized code's minimum distance.

The new directed replay covers every weight in its stated interval.
Sparse and dense classes must be joined separately before claiming a
whole-code distance guarantee. No implementation or production parameter
was changed for this sweep.

## Fixed construction and changed parameter

The scalar field is GF256. Each outer row is a GRS[128,k] codeword, and its
128 coordinates enter separate regions. The existing uniform regional
permutations and independent nonzero input and state diagonals are
unchanged. The inner remains the symmetric Cauchy construction with
packet width 16, state dimension eight, and multiplier two.

Only the outer dimension k changes. The source is naturally padded to
`k*64` records. For requested source size K, put

\[
L=64\left\lceil\frac{K}{64k}\right\rceil,
\qquad K_{\rm pad}=kL,\qquad N=128L.
\]

The forbidden output weight is `ceil(N/10)-1`. Thus exclusion of the bad
event proves relative symbol distance at least 10%. Setup randomness is
sampled once; there is no per-message resampling assumption.

## Directed middle certificates for dimension 88

The [replay](theory/gf8_outer_rate_middle.py) uses frozen rational witnesses
and the unchanged outward-rounded two-state transfer from the existing
state-eight certificate. It does not depend on the newer rank-cap
transfer or the truncated source polynomial. The exact modal-conditioning
factor `2^L / binom(L,L/2)` replaces `L+1` for each region, as justified in
[the preceding audit](GF8_RANK_HISTORY_AUDIT.md).

| Requested size | L | N | Actual requested rate | Covered input weights | Middle failure margin |
|---|---:|---:|---:|---:|---:|
| K18 | 3008 | 385024 | 32/47 = 68.085% | 15401–92405 | 129 bits |
| K20 | 11968 | 1531904 | 128/187 = 68.449% | 45958–367656 | 183 bits |

The K18 replay covers the interval `[.04N,.24N)` using 253 directed cells.
Its upward-rounded failure bound is

`8.5785624364798674100460268364344591514194118758784119720386711573447640295289579E-40`.

The partition hash is
`2dbf084acd15c1e04cabb0a387e382e4667936e758a65f4e26d939746584a138`.

The K20 replay covers `[.03N,.24N)` using 217 directed cells. Its
upward-rounded failure bound is

`4.8007348465347809928826915906595645090305007734299212561737701874920517534159822E-56`.

The partition hash is
`1a91d49d27785379e794ae9161b1432afc5c213ab60728af5c1ac64c0bb169b5`.

Every nonzero GRS[128,88] row has at least 41 nonzero coordinates. Thus
sparse certificates through 375 active rows at K18 and 1120 at K20 would
bridge to these respective intervals. Dense certificates must start at
weights 92406 and 367657. This observation specifies the required
coverage; it is not itself a certificate for either remaining class.

At K18, the `.03N` class fails the optimized two-state/modal screen by
about 41 bits. Raising the middle cutoff to `.04N` avoids that finite-size
loss. An extra frozen witness near density `.13` is also needed for the
directed replay; simply reusing the coarser K20 witnesses loses too much
between anchors.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_middle.py --dimension 88 --k-log 20 --certify
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_middle.py --dimension 88 --k-log 18 --certify --first 1/25
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_middle_tests.py
```

All seven tests pass. They check padding, the exact modal factor, agreement
with the previous floating transfer at a fixed witness, domination of
individual-weight sums, both complete directed replays, and rejected
invalid scope or floating witnesses.

## Where the next dimension fails the current proof

The refined screen retains the incoming state-rank cap, uses the checked
weight-one and weight-two return minima nine and nine, and applies modal
conditioning. It also removes infeasible outer supports when the required
source-tilt condition holds. Each reported number below is one optimized
input-weight class, not interval coverage or a construction-failure lower
bound.

| K20 outer dimension | Input density | Refined sampled margin |
|---:|---:|---:|
| 89 | .13 | -7294.905 bits |
| 89 | .14 | -7390.680 bits |
| 90 | .03 | -2429.391 bits |
| 90 | .06 | +493.807 bits |
| 90 | .10 | -8036.817 bits |
| 90 | .14 | -19307.657 bits |
| 90 | .20 | +106996.508 bits |
| 90 | .24 | +269997.015 bits |

All these optimizations report convergence. The dimension-89 screen at
density `.13` needed 46 iterations, one beyond the initial run's limit.
The two failing nearby classes show that merely raising the sparse cutoff
would leave a substantial middle obstruction. These failures do not prove
that another witness or proof method cannot close dimension 89.

The older two-state bound at dimension 92 fails by about 25099 bits at
`.10` and 43808 bits at `.14`. Larger dimensions were not exhaustively
optimized after the nearby dimension-89 obstruction became apparent.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_middle.py --dimension 89 --density .14 --rank-history --iterations 60
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_middle.py --dimension 90 --density .03 .06 .10 .14 .20 .24 --rank-history --iterations 60
```

## Natural padding limits the useful K16 sweep

At K16, outer dimensions 87 through 93 all use 768 rows and 98304 output
records. Their requested rate remains exactly 2/3. A larger nominal rate
does not reduce routed data for those sizes.

Dimension 94 is the first improvement in requested rate: it uses 704
rows and 90112 output records, giving rate 8/11. However, the two-state
screen with modal conditioning has margins approximately -3120, -4441,
and -4757 bits at densities `.10`, `.13`, and `.14`. All three optimizations
converged. No new K16 certificate is claimed.

The useful next step is to join the dimension-88 K18 and K20 class
certificates, then determine whether the modest reduction in output size
justifies an implementation variant. Raising the outer dimension much
further requires a stronger middle proof or a changed inner, not merely
new parameter labels.
