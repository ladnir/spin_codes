# An affine two-sum dense-cohort gate

The new affine g8 mixer passes bounded FLOAT gates for one and two dense
groups, at `r=17` and source matching counts `(8,0)`. The one-group bound is
`-184.772727559437`, covering all 219 fixed dense row sets. The two-group bound
is `-2001.873508517472`, covering all `219^2=47961` ordered fixed row-set
assignments. These bounds use different global witnesses. Neither is a
directed certificate, coverage of other source classes, or a whole-code claim.

The mixer is `H=I+11^T+uu^T` with AES-field labels `u=(8,...,15)`, bracketed
by the existing nonzero diagonals. Neither the old XOR2 support counts nor
its certificates enter this calculation. The inner regional operators are
reused because the new output diagonal supplies the same conditional
independent nonzero amplitude law.

## Fixed row sets and the new mask kernel

For each fixed dense row set `A`, the new helper computes the exact count
`N_A(T)` of inputs supported inside `A` whose affine-mixer output has support
exactly `T`. Coordinates inside `A` may be zero in this comparison measure.
The kernel uses ranks of restrictions of the new affine matrix, followed
by exact subset inversion. It is not a support-weight approximation.

The shortened-RS argument remains a counting-measure domination. On an
actual nonzero row support of size `w`, the number of compatible words is
at most `255^(w-16)`. Independent input diagonals give nonzero amplitude
mass `255^-w`. Their product is at most `255^-16` for that row. Relative to
uniform inputs on a fixed `a`-row coordinate space through 32 columns,
the comparison factor is `beta_a=256^(32a)/255^(16a)`.

Partition the 219 fixed row sets into six deterministic classes by arity.
For the class `B_a={A:|A|=a}`, define

\[
 W_a(T)^{32}=\sum_{A\in B_a}\frac{N_A(T)^{32}}{255^{16a}}.
\]

Hölder bounds the sum over one fixed row set across the entire 32-column
word by the product of these column weights. The row set is not resampled
between columns. Its multiplicity and `beta_a` are already inside `W_a`.
There is no extra `C(8,a)` or dense comparison factor outside the polynomial.

## Exact joint hole mixture with a labelled dense mask

In a sparse column category, let `s` two-row groups have one active input
and `d` have two. A singleton has full output support. Each double has full
support with probability `247/255` and each of eight one-hole supports with
probability `1/255`, independent of its fixed constituent-pair label.

Fix the dense output mask `T`. For `t=s+d`, the suffix recurrence keeps the
number `k` of labelled sparse holes still to allocate. At regional position
`i`, allocating `h` holes uses

\[
 V_i^T[k]=\sum_{h=0}^k\binom kh
 R_i\bigl(t-h+1_{i\in T}\bigr)V_{i+1}^T[k-h],
 \quad V_8^T[0]=(1,p),\quad V_8^T[k>0]=0.
\]

The regional matrices are the source-conditioned upper envelopes for the
same fixed source support. This replacement is pointwise over the possible
regional subsets, so it does not require those subsets to be independent.
The recurrence preserves their labelled chronology and the joint hole law.

After the suffix, sum `C(d,k)247^(d-k)V_0^T[k]/255^d`. Then average the
complete dense masks with `W_a(T)`, separately for empty and nonempty `T`.
Only then take the common-potential quotient. There is no column-occupancy
resource maximum or constituent-pair multiplicity stratum.

The suffix recurrence depends on `t`, not separately on `d`. One bank for
each `t=0,...,16` therefore serves all categories. A binary suffix code keeps
`T` labelled while sharing its repeated suffixes.

The positive column polynomial counts sparse support degree `X` with
`C(16,s)C(16-s,d)2^s255^(s+2d)` and dense nonemptiness with `Y`. Retain
`X>=544`, `Y>=17` in its 32nd power. The common external factor is

\[
 \frac{K_{17}(8,0)\,17\,28^{16}\,z^{-734003}}
 {255\,p\,255^{512}},\qquad
 K_{17}(8,0)=15708305994709246428578365440.
\]

The source coloring, sparse pair choices, projective factor, and endpoint
potential each occur once. The dense nonempty-column floor follows from
the invertibility of `H` and the RS minimum distance.

## Two frozen controls, no fit

The initial sparse witness was
`z=996248217/1000000000`, `log2(p)=-127/2`. It closed arities 3--6 but failed
at arities 7--8. The additional control reused only the literal of an old
dense-A8 witness; it recomputed the new affine operators and mask counts:

```text
z       = 8970988817957597/9007199254740992
log2(p) = -2266954815/33554432
```

| Dense arity | Fixed row sets | Initial FLOAT log2 bound | Additional FLOAT log2 bound |
|---|---:|---:|---:|
| 3 | 56 | -325.348787 | -395.741161 |
| 4 | 70 | -296.699180 | -387.326186 |
| 5 | 56 | -175.780840 | -379.556215 |
| 6 | 28 | -48.753624 | -371.222433 |
| 7 | 8 | +77.529253 | -308.850847 |
| 8 | 1 | +202.619598 | -184.772728 |

The additional literal wins every class, so the final gate uses one global
mark and potential; it needs no per-arity witness schedule. Both evaluations
took about two seconds. The first included 0.70 seconds for exact mask
kernel/moment construction and 0.78 seconds for regional and hole operators.
There was no optimizer, broader source sweep, native work, or directed replay.

Seven exact/light tests pass. They distinguish the affine singleton/pair
kernel from XOR2, check exact moments and fixed-class Hölder, compare the
offset recurrence with independent exact noncommuting calculations, and
verify direct labelled-hole enumeration, normalization, scope, and deterministic
class selection. Source/arithmetic identities are pinned below; byte hashes
identify FLOAT evidence and do not certify numerical upper bounds.

## Two-group triage at the same literal

A subsequent isolated control retained labelled ternary occupancies for
two dense groups, each with all eight constituent rows active. Its FLOAT
log2 bound was `+581.957974095917`; the inherited one-group witness fails.
This is a witness/control failure, not a construction obstruction. The
control does not cover all `219^2` fixed dense row-set assignments.

The exact positive convolution of the two new affine mask weights is
indexed by `ell in {0,1,2}^8` and total nonempty count `e in {0,1,2}`.
For each complete labelled `ell`, the sparse-hole recurrence uses
`R_i(t-h+ell_i)`. Ternary suffix sharing computes all 6,561 occupancies in
one bank for each `t=0,...,15`. Both the sparse-hole and dense-occupancy
mixtures finish before the potential quotient. The retained floors are
`X>=510`, `Y>=34`; the second is explicitly an aggregate relaxation.

The prefix has source coloring `C(17,2)=136`, sparse pair choices `28^15`,
normalizer `255^480`, and one projective/endpoint factor. This `a8 x a8`
class has exactly one fixed row-set assignment, so no dense type multiplier
is needed. The first control did not evaluate the other arity pairs.

The control took 3.36 seconds, including 2.32 seconds for operators and
mixtures. Its retained vector bank has 1,784,592 doubles, about 13.6 MiB.
The dominant coefficient is `(X,Y)=(957,64)`, so tightening only the aggregate
dense floor would not remove it. The run remained within its 120-second
guard. That first control evaluated no second witness or remaining type-pair tail.
Six additional exact/light tests passed, including every tiny ternary offset
against an exact noncommuting hole recurrence and the one-group specialization.

## Bounded h2 fit and complete arity-pair control

A subsequent fit used only the `a8 x a8` subclass. Its authorized rectangle
was `z in [0.992,0.998]`, `log2(p) in [-100,-40]`. It made 16 optimizer
evaluations, then one fresh FLOAT evaluation of the selected rational
literal, within the 24-evaluation limit. The full trial log is preserved.
The run took 45.73 seconds. Its winner is

```text
z       = 1116892707587883/1125899906842624
log2(p) = -40
```

The winner reaches both rectangle boundaries: the lower mark and upper
log potential. It is not an interior optimum. An earlier interior trial
already gave `-82.632734441748`; the winner gave `-2001.958048806847`.
The search rectangle was not extended.

One additional frozen-witness run then evaluated all 21 unordered arity
pairs `3 <= a <= b <= 8`, without fitting. The six classes `B_a` have
sizes `56,70,56,28,8,1`. For each pair, the positive mask convolution uses
`W_a(T_1)W_b(T_2)` and retains the labelled ternary occupancy and the number
of nonempty masks. Tensorized Hölder charges each fixed row-set pair once
across all 32 columns. The source coloring is `C(17,2)=136`.

For distinct arities, multiply the resulting word bound by two, once,
after the 32-column tail. This counts which selected dense group receives
which arity. Equal arities need no such factor, since `W_a` tensor `W_a`
already counts the ordered row-set choices. In particular,

\[
 \sum_{3\le a\le b\le8}(2-1_{a=b})\binom8a\binom8b=219^2.
\]

There is no extra dense row-set factor, comparison factor, or witness-bank
factor. Every pair uses the same global mark and potential. One shared
regional/ternary operator bank serves all 21 positive tails.

The complete FLOAT log2 union bound is `-2001.873508517472`. Every pair
passes; `a8 x a8` is largest and exactly reproduces the fit result. The
smallest cell is `a3 x a3` at `-2071.112909792575`. Every cell has dominant
coefficient `(X,Y)=(510,34)`. The artifact contains all 21 cell bounds,
exact row-set memberships, moment hashes, and the once-only multiplicities.
The run took 28.10 seconds: 0.67 seconds for exact moments, 2.20 seconds for
the shared operators, and the remainder for mask mixtures and positive tails.

All 18 affine dense tests pass, including exact noncommuting tiny operators,
complete ordered assignment enumeration, multiplicity placement after the
tail, frozen-literal identity, and refusal to report coverage after timeout.
An independent read-only review found no defect in the full 21-pair wrapper's
counts, floors, chronology, witness pin, or completion metadata.
The next numerical step is an independently reviewed directed replay of
this fixed 21-cell gate; no further fitting is needed for this subclass.

## Why width64/state8 is not a viable next inner gate

The stronger grouped outer does not remove the deterministic quiet-subcode
obstruction. This applies both to g8 with eight RS[32,16] rows per group
and to g16 with sixteen RS[16,8] rows per group. Both have 8,192 independent
outer groups, each with input dimension 128 and output length 256. Fix any
setup and select 2,048 whole groups. Their message space has dimension
`2048*128=262144`, and their outer outputs occupy at most
`2048*256=524288` positions. Each affine mixer acts within its own group;
the route preserves support size.

For packet width 64 there are 32,768 packets. With state dimension eight,
the preterminal equations `C_i X_i=0` impose at most
`8*32767=262136` homogeneous constraints. The least whole-group count that
ensures a nonzero kernel is exactly `floor(262136/128)+1=2048`; the strict
dimension surplus is eight. Any nonzero message in that kernel keeps every emitted state
zero, so `Y_i=X_i` throughout. Its relative weight is at most 25%, below
the 35% target, for every setup. This is the grouped version of
[the existing quiet-subcode argument](HD_BLOCKWISE_QUIET_SUBCODE.md).

Thus four RS[16,2] maps with width64/state8 cannot be rescued by stronger
within-group mixing. Width32/state8 changes the state-to-packet ratio and
is not rejected by this argument, but it would be a different construction
and cost analysis. No lighter-inner implementation or numerical gate was run.

For a next width32/state8 gate, two RS[16,4] maps preserve the current
state-to-packet ratio. The quiet-subcode calculation then selects
`floor(8*(65536-1)/128)+1=4096` whole groups, yielding a 50% upper cap.
That necessary screen does not reject 35%; it is not a distance guarantee.
The precise regional-cap replacement would be:

- Keep the rank-four RS[16,4] single-block kernel and return bounds; compose
  two blocks instead of four.
- Replace the live-state normalization by `256^8-1`, the refresh density by
  `(256^8-1)/255^8`, and the empty-input cap by its width32/state8/h4 version.
- For a packet's 16 routed pairs, use
  `B[s,f]=sum_(a+b=s) multinomial(16;a,b,f,16-a-b-f) M_(a+f,b+f)/2^s`.
  Raise this matrix polynomial to 256 packets per region. The regional
  normalization remains `C(4096,s)C(4096-s,f)`.

The source geometry, affine mask kernels, and outer support-floor tails
could consume the resulting new `R[m,f]` interface without changing their
counting arguments. Every numerical regional cap and composed gate would
still need fresh replay under a distinct construction identity.

There is no automatic factor-two native saving. The scalar A/C work per
output is unchanged. Refresh-diagonal products and diagonal-key bytes per
output are also unchanged: half the state is refreshed twice as often.
The possible saving is in the fixed eight-point refresh instead of the
sixteen-point refresh. The current inner packs four RS blocks across four
128-bit lanes. A literal two-block adaptation would fill only two lanes;
efficient 512-bit packing would need a new layout and may add shuffles.
Packet boundaries and serial refresh calls double. A next experiment should
first build the new regional caps for the frozen sparse and dense gates,
then inspect a concrete packing before any native timing claim.

## Reproduction and pins

From `theory`, the initial command requires a fresh output filename. The
second command reads the saved initial artifact named in its source.

```powershell
C:/Python314/python.exe -m unittest hd_two_sum_affine_dense_tests hd_two_sum_affine_dense_frozen_tests -v
C:/Python314/python.exe hd_two_sum_affine_dense.py --output fresh_affine_dense_initial.json
C:/Python314/python.exe hd_two_sum_affine_dense_frozen.py --output fresh_affine_dense_two_frozen.json
```

| File in `theory/` | SHA-256 |
|---|---|
| `hd_two_sum_affine_dense.py` | `31b387b8886c0a4cf2a07468812871231bcc965dae1cb26bbb0d64b7d5ab4063` |
| `hd_two_sum_affine_dense_tests.py` | `c4e8c3a16ab2cf99c1c0f64d3276fc4f87b5d6b5567f7888bdef371d49870246` |
| `hd_two_sum_affine_dense_frozen.py` | `99f65bd1e79207f49c89dd7f7a14ef29bf5caebcb50f60701b1094a9ff2752c1` |
| `hd_two_sum_affine_dense_frozen_tests.py` | `e07d58318c82a9a1a398197fd36f18108c1bbcd3b3eba11a6a5241cf5c3a8235` |
| `hd_two_sum_affine_dense_8_0_inherited.json` | `dbded746f5522f6356ad394a8cb3b11ba5a9c9cb03bfa3bd35a30126b5b31958` |
| `hd_two_sum_affine_dense_8_0_two_frozen.json` | `922d920286a6a513af0ae1a98ba41f51a5379c6fa26de2b9c249d5d64a393b9a` |
| `hd_two_sum_affine_dense_two.py` | `9e942d31e80c04a75ff43ad92b2b1230d89237fbc88e59568d2e80fa89eaf5c3` |
| `hd_two_sum_affine_dense_two_tests.py` | `566718f847c2a4676ac3acca658913ae87d6f66da23da2312511c2e7e4392f7d` |
| `hd_two_sum_affine_dense_h2_a8_8_0_frozen.json` | `a2d482c9c2168a1c1649e70842a0685726ddc6be6c1b8a015c4ceef71190d518` |
| `hd_two_sum_affine_dense_two_fit.py` | `8565c17dea8ed62b8a00f244eea5cc2cce207f454452af1f161d8aa40fbd2767` |
| `hd_two_sum_affine_dense_two_fit_tests.py` | `35cb2dc1ed5e32bee6f8c842603654369d09439f18d001646d388ba70296b54a` |
| `hd_two_sum_affine_dense_h2_a8_8_0_fit24.json` | `358a73269e326fb50fcc845d4c2511e9cdcb71b6e9f741dbc499dec25dfbd071` |
| `hd_two_sum_affine_dense_two_types.py` | `812bd4b322a243f25650b7d1c354ac035f7bc1c9e1412bb347e3c9fd59a69bff` |
| `hd_two_sum_affine_dense_two_types_tests.py` | `853c100b2a79039684569103716dc83c7ba0efc58a37a437e5b6af3f4084a05a` |
| `hd_two_sum_affine_dense_h2_all21_8_0_frozen.json` | `4f828393417127ebbc5ab2c6ef762c1cd92b1d55d730ca3a6f0c498e4928e54a` |

Reproduce the two-group triage separately with a fresh output filename:

```powershell
C:/Python314/python.exe -m unittest hd_two_sum_affine_dense_two_tests -v
C:/Python314/python.exe hd_two_sum_affine_dense_two.py --output fresh_affine_h2_a8.json
```

Reproduce only the approved full arity-pair gate, without a new fit:

```powershell
C:/Python314/python.exe -m unittest hd_two_sum_affine_dense_two_types_tests -v
C:/Python314/python.exe hd_two_sum_affine_dense_two_types.py --output fresh_affine_h2_all21.json
```
