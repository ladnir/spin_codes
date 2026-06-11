# Dense+Dense Audit Packet

This packet is for reviewing the current dense+dense finite certificate work in
`C:\Users\peter\repo\permute_conv`.

The intended audit target is narrow:

- finite checkpoint: `N = 2^21`, relative distance `delta = 0.09`, `d = floor(delta N) = 188743`
- outer model: direct sum of `4096` copies of the binary `RM(4,9)` block code, i.e. local `[512,256,32]`
- inner model: full-codeword random-split dense recursive inner with block size `b = 64`, using the EBCH `[128,64]` weight distribution
- quantity being certified: a first-moment upper bound for the currently isolated late-window dense+dense contribution

The current finite ledger reports

```text
log2 mu_finite <= -34.767174
```

This is a finite first-moment certificate for the audited late-window split. It is not yet a polished global theorem.

## Construction Under Audit

The active inner is the full-codeword random-split recursive inner, not the older scalar fixed-tap inner.

At each block:

1. The current input block and state form a mismatch `V_i = U_i + S_{i-1}`.
2. If `V_i = 0`, the branch emits zero and the state becomes zero.
3. If `V_i != 0`, a dense scrambler makes the local EBCH codeword uniform over nonzero codewords.
4. The `128` codeword coordinates are randomly split into `64` output coordinates and `64` next-state coordinates.
5. The effective one-step termination atom used in the certificate is

```text
log2 p_term = -63.89264923803805
```

The relevant manuscript section starts in `innerDense.tex`, especially the full-split episode material and finite ledger around:

```text
innerDense.tex:2944
innerDense.tex:3155
innerDense.tex:3412
```

## Main Verification Commands

Run these first:

```powershell
python scripts\verify_fullsplit_finite_ledger.py
python scripts\check_fullsplit_T_monotonicity.py --sufficient-reduction
pdflatex -interaction=nonstopmode main_permConv.tex
```

Expected key outputs:

```text
prefix_32_500_e_le8_log2,-34.767174
prefix_32_500_e_ge9_tail_log2,-269.335258
postprefix_501_2000_eall_log2,-182.259739
interval_2001_1148736_total_log2,-588.081877
finite_ledger_total_log2,-34.767174
finite_ledger_margin_bits,34.767174
```

For the monotonicity audit, every reported middle/far bucket slack should be positive. The smallest current slack is
the far-bucket row for `725001--950000`, about `3.773620` bits.

## Ledger Pieces

The combined finite ledger has four pieces.

### 1. Tiny Prefix, Explicit Episode Grid

Range:

```text
32 <= h <= 500
1 <= r <= 64
gap buckets: 1--4000, 4001--8000, 8001--12000
episode count: e <= 8
```

Artifact:

```text
scripts/fullsplit_piecewise_h32_500_csv.csv
```

Current total:

```text
log2 mu_32_500_e_le8 = -34.767174
```

Peak:

```text
h = 32
r = 1
remaining H = 31
gap bucket = 1--4000
term_log2 = -37.38576684150053
```

This is the dominant piece of the whole finite ledger.

The per-`T` inner grids feeding this row are:

```text
scripts/fast_fullsplit_e08_T5949_H0_499_all.csv
scripts/fast_fullsplit_e08_T9949_H0_499_all.csv
scripts/fast_fullsplit_e08_T13949_H0_499_all.csv
```

The matching right-endpoint grids are:

```text
scripts/fast_fullsplit_e08_T9948_H0_499_all.csv
scripts/fast_fullsplit_e08_T13948_H0_499_all.csv
scripts/fast_fullsplit_e08_T17948_H0_499_all.csv
```

Endpoint audit:

```powershell
python scripts\check_fullsplit_exact_grid_endpoints.py
```

Current output:

```text
bucket,H_min,H_max,max_total_right_minus_left,H_total,left_total,right_total,max_column_right_minus_left,column,H_column
gap_1_4000,0,499,-43.63206913,1,-20.0117348926,-63.6438040227,0.00201643022984,e8_log2,7
gap_4001_8000,0,499,0,0,-63.892649238,-63.892649238,0.000859762181619,e8_log2,7
gap_8001_12000,0,499,0,0,-63.892649238,-63.892649238,0.000476502560048,e8_log2,7
```

Interpretation: the exact `e <= 8` totals are endpoint-safe for the full
`H = 0..499` table. The small positive column wiggles are individual fixed-`e`
columns, not the log-summed total used in the ledger. This does not yet prove
monotonicity for every interior `T` in the bucket; that remains a proof
obligation.

The row-level prefix sum can be regenerated with:

```powershell
python scripts\sum_fullsplit_piecewise_certificate.py --h-values 32:500 --inner-mode-by-gap csv,csv,csv --inner-knot-csvs "scripts\fast_fullsplit_e08_T5949_H0_499_all.csv;scripts\fast_fullsplit_e08_T9949_H0_499_all.csv;scripts\fast_fullsplit_e08_T13949_H0_499_all.csv" --require-knot-coverage --output-csv scripts\fullsplit_piecewise_h32_500_csv.csv
```

Audit priority: high. This is the narrow waist of the certificate. After the
right-endpoint check, the remaining issue is the analytic or finite certified
control of the interior `T` values.

### 2. Tiny Prefix, Crude Episode Tail

Range:

```text
32 <= h <= 500
e >= 9
```

Artifact:

```text
scripts/fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv
```

Current total:

```text
log2 mu_32_500_e_ge9_tail = -269.335258
```

Peak:

```text
h = 500
r = 1
gap bucket = 8001--12000
term_log2 = -271.14116944335876
```

This is intentionally pessimistic and far below the dominant `e <= 8` prefix.

Audit priority: medium. Check that the crude `sum_{e>=9} binom(H+1,e) p_term^e` tail really applies to all omitted
selected-gap configurations.

### 3. Post-Prefix, All-Episode Wrapper

Range:

```text
501 <= h <= 2000
```

Artifact:

```text
scripts/fullsplit_piecewise_h501_2000_eall_hsummary.csv
```

Current total:

```text
log2 mu_501_2000_eall = -182.259739
```

Peak:

```text
h = 501
r = 1
gap bucket = 1--4000
mode = cap
peak term_log2 = -184.39217327424365
```

Regenerate with:

```powershell
python scripts\sum_fullsplit_piecewise_certificate.py --h-values 501:2000 --inner-mode-by-gap cap,eallratio,eallratio --gap-sum-mode exact --inner-T-by-gap min --output-h-summary-csv scripts\fullsplit_piecewise_h501_2000_eall_hsummary.csv
```

Audit priority: medium. This is now theorem-facing and no longer relies on an `e <= 1` caveat.

### 4. Interval Certificate

Range:

```text
2001 <= h <= 1148736
```

The interval constants are encoded in:

```text
scripts/verify_fullsplit_finite_ledger.py
scripts/check_fullsplit_T_monotonicity.py
```

Current interval rows:

```text
2001--7858          (lambda,rho)=(0.02,0.01)   log2 mu=-588.081877
7859--20550         (0.05,0.03)                log2 mu=-2657.429602
20551--75000        (0.2,0.1)                  log2 mu=-6027.058821
75001--250000       (0.5,0.3)                  log2 mu=-12846.353611
250001--350000      (1.2,0.5)                  log2 mu=-28478.125315
350001--400000      (1.2,1.0)                  log2 mu=-81557.079312
400001--450000      (1.2,1.0)                  log2 mu=-91423.257760
450001--550000      (1.2,1.0)                  log2 mu=-82083.446102
550001--650000      (1.2,1.0)                  log2 mu=-1519.978235
650001--725000      (1.2,2.0)                  log2 mu=-116329.618720
725001--750000      (1.2,3.0)                  log2 mu=-169889.387695
750001--850000      (1.2,3.0)                  log2 mu=-134525.068420
850001--950000      (1.2,3.0)                  log2 mu=-117536.328304
950001--1050000     (1.2,3.0)                  log2 mu=-285054.035442
1050001--1148736    (1.2,3.0)                  log2 mu=-562862.201412
```

The combined interval total is still dominated by `2001--7858`:

```text
log2 mu_2001_1148736 = -588.081877
```

The cutoff `h = 1148736` is the feasibility edge for the current far bucket:

```text
T_max = 5949 + 12000 - 1 = 17948
64*T_max = 1148672
h <= H + r <= 1148672 + 64 = 1148736
```

Audit priority: medium. Check the endpoint placement bound, the `T_eff=max(T_min,ceil(H/b))` rule, and the sufficient
monotonicity reduction.

## Proof Objects To Inspect

Core scripts:

```text
scripts/verify_fullsplit_finite_ledger.py
scripts/sum_fullsplit_piecewise_certificate.py
scripts/check_fullsplit_T_monotonicity.py
scripts/fast_fullsplit_episode_e01.py
scripts/sum_fullsplit_turnoff_tail.py
```

Core data:

```text
scripts/rm512_256_spectrum.csv
scripts/EBCH128_64.wd
scripts/fast_fullsplit_e08_T5949_H0_499_all.csv
scripts/fast_fullsplit_e08_T9949_H0_499_all.csv
scripts/fast_fullsplit_e08_T13949_H0_499_all.csv
scripts/fullsplit_piecewise_h32_500_csv.csv
scripts/fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv
scripts/fullsplit_piecewise_h501_2000_eall_hsummary.csv
```

Manuscript locations:

```text
innerDense.tex
integration.tex
outerDense.tex
```

The active finite certificate material is in `innerDense.tex`; `integration.tex` still contains older random-banded and
scalar-inner history that should not be read as the current proof target.

## What Is Proved Versus Not Yet Proved

Currently checkable:

- the finite ledger arithmetic in `verify_fullsplit_finite_ledger.py`
- the `501..2000` all-episode wrapper row
- the high-`h` fixed-pole interval rows as numerical certificates
- the sufficient endpoint monotonicity inequalities for the listed intervals
- the finite `32..500` row-level sum from the existing CSV artifacts

Still proof debt:

- convert the `32..500`, `e <= 8` finite grid from "CSV artifact" to an audit-grade finite lemma
- verify the derivation and implementation of the effective termination atom `log2 p_term = -63.89264923803805`
- add interval-arithmetic or rational/integer safeguards for the most important numerical bounds
- state the coverage of first-active placement regimes cleanly, including what is covered by the late-window split and
  what is handled elsewhere
- make the construction definition and boundary convention crisp enough that every script is visibly evaluating the
  same object

## Do Not Trust Yet

The following are useful context but should not be treated as proof:

- old scalar fixed-tap dense-inner conclusions
- old random-banded outer `sigma` plots
- BCH/RM replacement speculation outside the RM `[512,256,32]` checkpoint
- fit-based residual laws
- any `tmp_*` CSV or PNG unless it is explicitly named above
- fixed-`h` asymptotic heuristics that are not connected to the first-moment sum

## Recommended Audit Order

1. Run the three main verification commands.
2. Check the exact formulas in `fast_fullsplit_episode_e01.py`, especially the selected-gap count and endpoint
   `e=x+1` case.
3. Check `sum_fullsplit_piecewise_certificate.py` for placement factors, outer spectrum bounds, and bucket handling.
4. Check `check_fullsplit_T_monotonicity.py` against the monotonicity argument in `innerDense.tex`.
5. Check that the `32..500` finite CSV inputs are generated from formulas that are conservative upper bounds.
6. Only then read the surrounding prose and decide what should become theorem text versus working-save-point text.
