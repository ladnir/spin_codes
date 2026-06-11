# Dense+Dense Audit Packet

This packet is for reviewing the current dense+dense finite certificate work in
`C:\Users\peter\repo\permute_conv`.

The intended audit target is narrow:

- finite checkpoint: `N = 2^21`, relative distance `delta = 0.09`, `d = floor(delta N) = 188743`
- outer model: direct sum of `4096` copies of the binary `RM(4,9)` block code, i.e. local `[512,256,32]`
- inner model: full-codeword random-split dense recursive inner with block size `b = 64`, using the EBCH `[128,64]` weight distribution
- quantity being certified: a first-moment upper bound for the currently isolated dense+dense contribution, with full first-active-position coverage for `32 <= h <= 500` and late-window coverage beyond that

The current finite ledger reports

```text
log2 mu_finite <= -34.767174
log2 mu_late_plus_window <= -34.732187
log2 mu_32_500_all_first_active_positions <= -34.732187
```

The first line is still dominated by the audited window `5949 <= T <= 17948`. The second also includes the checked
placement-only ultra-late prefix `T < 5949`, through `h <= 500`; the third records that the new early
`T > 17948`, `h <= 500` slice is now included as well. It is not yet a polished global theorem.

## Construction Under Audit

The active inner is the full-codeword random-split recursive inner, not the older scalar fixed-tap inner.

At each block:

1. The current input block and state form a mismatch `V_i = U_i + S_{i-1}`.
2. If `V_i = 0`, the branch emits zero and the state becomes zero.
3. If `V_i != 0`, a dense scrambler makes the local EBCH codeword uniform over nonzero codewords.
4. The `128` codeword coordinates are randomly split into `64` output coordinates and `64` next-state coordinates.
5. The conservative effective one-step termination atom used in the finite-prefix certificate is

```text
log2 p_term <= -63.8926492
```

This is backed by:

```powershell
python scripts\verify_fullsplit_turnoff_atom.py
```

Expected key output:

```text
self_turnoff_p0_log2,-64.00457432724919
eta_sup_H,499
eta_sup_T,5949
eta_sup_log2,-67.63641076470455
computed_pterm_log2,-63.89264922047382
certified_pterm_upper_log2,-63.89264920000000
extra_turnoff_for_certified_upper_log2,-67.63641049043137
status,PASS
```

Scope: this sharper exact-cancellation add-on is certified for the current tiny-prefix
domain `0 <= H <= 499`, `T >= 5949`. Outside this tiny prefix, use the global
Bernstein envelope below.

For theorem-facing wrappers outside the tiny prefix, use the global Bernstein
envelope:

```powershell
python scripts\verify_fullsplit_global_turnoff_envelope.py
```

Expected key output:

```text
full_state_q_log2,-63.00228535065820
finite_population_correction_log2,0.01527942647943
eta_envelope_log2,-62.98700592417877
computed_pterm_envelope_log2,-62.40787581907139
certified_global_pterm_upper_log2,-62.40787580000000
status,PASS
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
python scripts\verify_fullsplit_turnoff_atom.py
python scripts\verify_fullsplit_global_turnoff_envelope.py
python scripts\verify_fullsplit_high_interval_adjustment.py
python scripts\verify_fullsplit_finite_ledger.py
python scripts\check_fullsplit_T_monotonicity.py --sufficient-reduction
pdflatex -interaction=nonstopmode main_permConv.tex
```

Expected key outputs:

```text
prefix_interior_audit_rows,1500
prefix_interior_audit_worst,bucket=gap_8001_12000,H=75,diff=2.19205276153e-09,T=13949
late_prefix_T_lt_5949_cap_log2,-40.115431
prefix_32_500_e_le8_log2,-34.767174
prefix_32_500_e_ge9_tail_log2,-269.335258
postprefix_501_2000_eall_log2,-182.259739
interval_2001_1148736_total_log2,-586.597103
finite_ledger_total_log2,-34.767174
finite_ledger_margin_bits,34.767174
late_plus_window_total_log2,-34.732187
late_plus_window_margin_bits,34.732187
```

For the monotonicity audit, every reported middle/far bucket slack should be positive. The smallest current slack is
the far-bucket row for `725001--950000`, about `3.773620` bits.

## Ledger Pieces

The combined finite ledger now has six checked pieces.

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

The matching right-endpoint grids, used for the earlier endpoint guard, are:

```text
scripts/fast_fullsplit_e08_T9948_H0_499_all.csv
scripts/fast_fullsplit_e08_T13948_H0_499_all.csv
scripts/fast_fullsplit_e08_T17948_H0_499_all.csv
```

Full interior-`T` audit:

```powershell
python scripts\check_fullsplit_exact_grid_interior.py --h-values 0:499 --sample-step 1 --method combined --summary-only --output-csv scripts\fullsplit_exact_interior_h0_499_allT.csv
```

Current output:

```text
Full-split exact-grid interior audit
summary,rows=1500,worst_bucket=gap_8001_12000,worst_H=75,worst_diff=2.19205276153e-09,worst_T=13949,worst_repro_bucket=gap_8001_12000,worst_repro_H=185,worst_repro_delta=-5.52478240934e-09
```

Interpretation: the exact `e <= 8` totals were checked for every
`H = 0..499` and every `T` in the three placement buckets. The worst positive
value is at the left endpoint itself and is numerical roundoff, not an interior
increase. This is a finite audit, not the eventual analytic monotonicity lemma,
but it removes the dominant `T`-reuse assumption from the current finite ledger.
The standard finite-ledger verifier now checks this tracked audit artifact by
default before summing the prefix rows.

The row-level prefix sum can be regenerated with:

```powershell
python scripts\sum_fullsplit_piecewise_certificate.py --h-values 32:500 --inner-mode-by-gap csv,csv,csv --inner-knot-csvs "scripts\fast_fullsplit_e08_T5949_H0_499_all.csv;scripts\fast_fullsplit_e08_T9949_H0_499_all.csv;scripts\fast_fullsplit_e08_T13949_H0_499_all.csv" --require-knot-coverage --output-csv scripts\fullsplit_piecewise_h32_500_csv.csv
```

The `e <= 8` endpoint grids are regenerated with the certified cancellation
add-on:

```powershell
python scripts\fast_fullsplit_episode_e01.py --remaining-blocks 5949 --remaining-ones 0:499 --e-max 8 --extra-turnoff-log2 -67.63641049043137 --output-csv scripts\fast_fullsplit_e08_T5949_H0_499_all.csv
python scripts\fast_fullsplit_episode_e01.py --remaining-blocks 9949 --remaining-ones 0:499 --e-max 8 --extra-turnoff-log2 -67.63641049043137 --output-csv scripts\fast_fullsplit_e08_T9949_H0_499_all.csv
python scripts\fast_fullsplit_episode_e01.py --remaining-blocks 13949 --remaining-ones 0:499 --e-max 8 --extra-turnoff-log2 -67.63641049043137 --output-csv scripts\fast_fullsplit_e08_T13949_H0_499_all.csv
```

Audit priority: high. This is the narrow waist of the certificate. The finite
interior-`T` control is now checked; the remaining proof debt is to either keep
this as an explicit finite certificate or replace it by a clean analytic
monotonicity lemma.

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

Regenerate with:

```powershell
python scripts\sum_fullsplit_turnoff_tail.py --outer-prefix-csv scripts\block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact.csv --h-values 32:500 --first-r-values 1:64 --turnoff-log2 -63.8926492 --e-min 9 --output-csv scripts\fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv
```

Audit priority: medium. Check that the crude `sum_{e>=9} binom(H+1,e) p_term^e` tail really applies to all omitted
selected-gap configurations.

### 3. Early Small-Weight Prefix, Uniform-Survival Grid

Range:

```text
32 <= h <= 500
1 <= r <= 64
gap buckets: 12001--17000, 17001--22000, 22001--26819
episode count: e <= 16
```

Artifacts:

```text
scripts/fullsplit_early_e16_uniformsurv_T17949_22948_H0_499.csv
scripts/fullsplit_early_e16_uniformsurv_T22949_27948_H0_499.csv
scripts/fullsplit_early_e16_uniformsurv_T27949_32767_H0_499.csv
scripts/fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv
```

Current total:

```text
log2 mu_32_500_early_e_le16 = -85.403338
```

Bucket totals:

```text
12001--17000: -85.815880
17001--22000: -87.831441
22001--26819: -89.394315
```

Peak:

```text
h = 32
r = 1
remaining H = 31
gap bucket = 12001--17000
inner_log2 = -108.576588
term_log2 = -87.36339217954111
```

Regenerate the bucket-uniform inner CSVs with:

```powershell
python scripts\build_fullsplit_early_uniform_survival.py --e-max 16
```

Then regenerate the row-level ledger with:

```powershell
python scripts\sum_fullsplit_piecewise_certificate.py --h-values 32:500 --first-r-values 1:64 --gap-start 12001 --gap-stop 26819 --gap-step 5000 --inner-mode-by-gap csv,csv,csv --inner-knot-csvs "scripts\fullsplit_early_e16_uniformsurv_T17949_22948_H0_499.csv;scripts\fullsplit_early_e16_uniformsurv_T22949_27948_H0_499.csv;scripts\fullsplit_early_e16_uniformsurv_T27949_32767_H0_499.csv" --require-knot-coverage --gap-sum-mode exact --inner-T-by-gap min --turnoff-log2 -63.8926492 --output-csv scripts\fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv --output-h-summary-csv scripts\fullsplit_piecewise_early_h32_500_e16_uniformsurv_hsummary.csv
```

Interpretation: this is not an endpoint monotonicity claim. The builder uses a bucket-uniform denominator and a
bucket-uniform skipped-gap suffix while keeping the live-span Chernoff survival factor.

### 4. Early Small-Weight Prefix, Crude Episode Tail

Range:

```text
32 <= h <= 500
gap buckets: 12001--17000, 17001--22000, 22001--26819
e >= 17
```

Artifact:

```text
scripts/fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv
```

Current total:

```text
log2 mu_32_500_early_e_ge17_tail = -305.738816
```

Peak:

```text
h = 500
r = 1
gap bucket = 22001--26819
term_log2 = -306.4449663093145
```

Regenerate with:

```powershell
python scripts\sum_fullsplit_turnoff_tail.py --outer-prefix-csv scripts\block_outer_probe_k1048576_rm512_256_sig32_d009_h500_exact.csv --h-values 32:500 --first-r-values 1:64 --gap-start 12001 --gap-stop 26819 --gap-step 5000 --turnoff-log2 -63.8926492 --e-min 17 --output-csv scripts\fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv
```

Audit priority: medium. The crude tail is safe only after pushing the explicit uniform-survival grid to `e <= 16`;
starting the crude tail at `e >= 9` fails badly in the far early bucket.

### 5. Early High-Density Probe, Paired T

Range under investigation:

```text
h > 500
T > 17948
```

Important finding: the ordinary endpoint-placement interface is too pessimistic for high-density early rows because it
pairs placement at `T_max` with survival at `T_min`. For large `H`, placement inside the bucket is dominated by large
`T`, where the survival exponent is also stronger.

New probe:

```text
scripts/probe_fullsplit_early_paired_survival.py
```

Representative formerly bad row:

```powershell
python scripts\probe_fullsplit_early_paired_survival.py --h-values 1100001 --gap-start 22001 --gap-stop 26819 --gap-step 5000 --episode-slack-bits 20000 --lambdas 0.05:20:0.05
```

Current output:

```text
total_log2 = -113084.867453
peak_h = 1100001
peak_first_r = 34
peak_T = 32767
peak_lambda = 2.3
```

Sparse high-density probes with the same deliberately large `20000`-bit episode allowance:

```powershell
python scripts\probe_fullsplit_early_paired_survival.py --h-values 300001:1200000:100000,1200000 --gap-start 12001 --gap-stop 26819 --gap-step 5000 --episode-slack-bits 20000 --lambdas 0.05:20:0.05
python scripts\probe_fullsplit_early_paired_survival.py --h-values 1200001:2097152:100000,2097152 --gap-start 12001 --gap-stop 26819 --gap-step 5000 --episode-slack-bits 20000 --lambdas 0.05:20:0.05
```

Current outputs:

```text
300001..1200000 step 100000: log2 <= -112940.527337
1200001..2097152 step 100000: log2 <= -111645.563005
```

Audit priority: high but theorem-facing, not finite-artifact yet. The next proof step is a paired placement-survival
interval lemma plus a same-`T` episode multiplier. Do not use the separable endpoint-placement bound for this
high-density early branch.

Follow-up diagnostic: the near-full-density endpoint also needs the occupancy Cauchy pole on the correct side of the
polynomial. For `T=32767`, `H=bT=2097088`, `lambda=2.3`, the eall-ratio inner bound is `-893425.892127` at
`rho=10`, but improves to `-1181735.473865` at `rho=1e6`, matching the no-turnoff survival exponent. So the
full-density edge is a coefficient-envelope artifact, not a new episode obstruction. The wrapper default `rho` grid
now includes large poles through `1e4`; a later theorem version may state this as a deficit-side occupancy bound.

Outer high-weight hook: the RM512 local spectrum is complement-symmetric, so the direct-sum outer has
`A_h = A_{N-h}`. The paired probe and piecewise certificate now have `--outer-complement-symmetry` to use this fact
instead of building outer tables up to `h=N`. Smoke checks at `h=N`, `r=64`, far early bucket:

```powershell
python scripts\probe_fullsplit_early_paired_survival.py --h-values 2097152 --first-r-values 64 --gap-start 22001 --gap-stop 26819 --gap-step 5000 --outer-complement-symmetry --episode-slack-bits 0 --lambdas 2.3
python scripts\sum_fullsplit_piecewise_certificate.py --h-values 2097152 --first-r-values 64 --gap-start 22001 --gap-stop 26819 --gap-step 5000 --inner-mode-by-gap eallratio --gap-sum-mode endpoint --inner-T-by-gap feasible-min --outer-complement-symmetry --turnoff-log2 -62.4078758 --lambdas 2.3 --rhos 1000000
```

The paired smoke gives `total_log2 = -1181735.473865`; the endpoint piecewise smoke gives `-1181723.239347`, exactly
adding the `log2(4819)` bucket overcount.

Paired all-episode diagnostics: `scripts/probe_fullsplit_early_paired_survival.py` now has
`--inner-mode eallratio`, which attaches the same fixed-pole all-episode wrapper to each exact `T` placement term.
A sparse complement-symmetric survival scan over `h=501,100501,...,2000501,N` peaks at
`h=1000501`, `r=31`, `T=32767`, with `total_log2 = -136323.731114`. Local all-episode windows around the survival
peaks give:

```text
h=1000501, r=25..37: total_log2 = -133179.527574, peak_r=32
h=1150501, r=30..42: total_log2 = -132946.839848, peak_r=32
h=1200501, r=32..43: total_log2 = -132647.058774, peak_r=32
h=1250501, r=33..45: total_log2 = -132835.935679, peak_r=33
```

Endpoint/full-`r` profile at the same sampled high-density points, using `gap-start=gap-stop=26819` and `r=1..64`:

```text
h=1000501: total_log2 = -133179.355669
h=1150501: total_log2 = -132946.388557
h=1200501: total_log2 = -132646.193762, peak_r=32, peak_term=-132649.525145
h=1250501: total_log2 = -132834.786925
```

At the sampled worst `h=1200501`, the `r`-sum costs only `3.331383` bits over the peak term, while the single
`r=32` all-`T` sum is `-132649.521954`, only `0.003191` bits over the endpoint peak. This points to the final proof
shape: endpoint dominance in `T`, plus an exact first-block hypergeometric `r`-tail bound.

The theorem-facing fast helper is now `scripts/sum_fullsplit_paired_early_eall.py`. It splits the paired inner bound
into an `e=0` branch and an `e>=1` fixed-pole Cauchy branch, then sums the same-`T` placement terms directly. This
matches the endpoint/full-`r` probe exactly at displayed precision:

```text
python scripts\sum_fullsplit_paired_early_eall.py --h-values 1200501 --first-r-values 1:64 --gap-start 26819 --gap-stop 26819 --gap-step 5000 --outer-complement-symmetry --lambdas 2.3 --rhos 1
total_log2 = -132646.193762
peak_h = 1200501
peak_first_r = 32
peak_T = 32767
peak_branch = ege1
```

It also makes the full far bucket cheap:

```text
python scripts\sum_fullsplit_paired_early_eall.py --h-values 1000501,1150501,1200501,1250501 --first-r-values 1:64 --gap-start 22001 --gap-stop 26819 --gap-step 5000 --outer-complement-symmetry --lambdas 2.3 --rhos 1
total_log2 = -132646.190570
endpoint_total_log2 = -132646.193762
full_minus_endpoint_bits = 0.003192
ege1_T_endpoint_geom_bound: lambda=2.3, rho=1, log2_g=8.821856, overhead_bits=0.003199
h=1000501: -133179.352477
h=1150501: -132946.385366
h=1200501: -132646.190570
h=1250501: -132834.783733
branch_e0_log2 = -136323.731114
branch_ege1_log2 = -132646.190570
```

For `h=1200501`, the full far bucket is only `0.003192` bits above the endpoint/full-`r` value. The high-density
`T` overhead is now backed by Lemma `lem:fullsplit-paired-endpoint-dominance`: with `lambda=2.3, rho=1`,
`log2_g=8.821856`, so the geometric right-endpoint overhead bound is `0.003199` bits. The high-density proof target
is therefore sharper than before: handle the first-block hypergeometric `r`-sum or tails without a heuristic local-`r`
window, then close the remaining `h` interval/complement-spectrum envelope.

Important correction from a broader endpoint-only scan: `h=1200501` is not the sampled high-density peak. Scanning
`h=1100501:1300501:5000` moved the peak to the `1115k` neighborhood, and a step-one zoom over
`1115101:1115201:1` found the sampled endpoint peak at:

```text
h = 1115146
endpoint_h_log2 = -132529.972886
endpoint_peak_r = 32
endpoint_peak_term_log2 = -132533.304269
endpoint_total_minus_endpoint_peak_bits = 3.331383
```

Running the full far bucket at this corrected sampled peak gives:

```text
python scripts\sum_fullsplit_paired_early_eall.py --h-values 1115146 --first-r-values 1:64 --gap-start 22001 --gap-stop 26819 --gap-step 5000 --outer-complement-symmetry --lambdas 2.3 --rhos 1
total_log2 = -132529.969694
endpoint_total_log2 = -132529.972886
full_minus_endpoint_bits = 0.003192
```

This changes the honest bookkeeping but not the qualitative conclusion: the `T` and first-block `r` costs look small
and stable.

The complement-side interval closure is now handled by:

```text
python scripts\certify_fullsplit_high_density_interval.py --intervals "1048577:1148736:1;1148737:1300000:1;1300001:1500000:1;1500001:1700000:1;1700001:1900000:10;1900001:2097089:10;2097090:2097152:10" --lambda-value 2.3
```

It gives:

```text
1048577--1148736, rho=1:  -128915.315746
1148737--1300000, rho=1:  -124889.651755
1300001--1500000, rho=1:  -113710.357653
1500001--1700000, rho=1:   -35150.097366
1700001--1900000, rho=10: -524628.366199
1900001--2097089, rho=10: -894140.350713
2097090--2097152, rho=10: -893792.285010
total_log2 = -35150.097366
```

`verify_fullsplit_finite_ledger.py` now records these rows as `complement_interval_*`; the optional
`--check-complement-high-intervals all` recomputes all seven rows and passes. Current interpretation: early
high-density has ample numerical margin once `T` pairing, outer complement symmetry, and the fixed-`z` interval
endpoint bound are used. The remaining work is to polish the proof prose and remove avoidable overcount/overlap with
the older `2001--1148736` interval table.

### 6. Post-Prefix, All-Episode Wrapper

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

### 7. Interval Certificate

Range:

```text
2001 <= h <= 1148736
```

The interval constants are encoded in:

```text
scripts/verify_fullsplit_finite_ledger.py
scripts/check_fullsplit_T_monotonicity.py
```

The ledger manifest can also print the exact theorem-facing recomputation command for every high interval:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --print-high-interval-commands
```

Each printed row has the form

```powershell
python scripts\sum_fullsplit_piecewise_certificate.py --h-values A:B --inner-mode-by-gap cap,eallratio,eallratio --gap-sum-mode endpoint --inner-T-by-gap feasible-min --turnoff-log2 -62.4078758 --lambdas L --rhos R
```

To spot-check selected rows, run for example:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --check-high-intervals 2001--7858
python scripts\verify_fullsplit_finite_ledger.py --check-high-intervals 550001--650000
```

The second command is intentionally slower because it recomputes the awkward high-transition interval.

Current interval rows:

```text
2001--7858          (lambda,rho)=(0.02,0.01)   log2 mu=-586.597103
7859--20550         (0.05,0.03)                log2 mu=-2655.944828
20551--75000        (0.2,0.1)                  log2 mu=-6025.574047
75001--250000       (0.5,0.3)                  log2 mu=-12844.868837
250001--350000      (1.2,0.5)                  log2 mu=-28476.640541
350001--400000      (1.2,1.0)                  log2 mu=-81555.594538
400001--450000      (1.2,1.0)                  log2 mu=-91421.772986
450001--550000      (1.2,1.0)                  log2 mu=-82081.961328
550001--650000      (1.2,1.0)                  log2 mu=-1518.493461
650001--725000      (1.2,2.0)                  log2 mu=-116328.133946
725001--750000      (1.2,3.0)                  log2 mu=-169887.902921
750001--850000      (1.2,3.0)                  log2 mu=-134523.583646
850001--950000      (1.2,3.0)                  log2 mu=-117534.843530
950001--1050000     (1.2,3.0)                  log2 mu=-285052.550668
1050001--1148736    (1.2,3.0)                  log2 mu=-562860.716638
```

The combined interval total is still dominated by `2001--7858`:

```text
log2 mu_2001_1148736 = -586.597103
```

These rows include a checked `+1.484774` bit adjustment for switching the
high-interval all-episode wrapper from the finite-prefix atom to the global
Bernstein atom. The verifier

```powershell
python scripts\verify_fullsplit_high_interval_adjustment.py
```

shows that the required adjustment is `1.484773404861` bits, leaving about
`5.95e-7` bits of slack. A spot recomputation of the leading interval gives
`-586.597103`, matching the adjusted row. The current spot checks also verify
the transition interval:

```text
interval_2001--7858_recomputed_log2,-586.597103
interval_550001--650000_recomputed_log2,-1518.493462
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
- the high-`h` fixed-pole interval rows as numerical certificates, with printed regeneration commands and selected
  opt-in recomputation checks
- the sufficient endpoint monotonicity inequalities for the listed intervals
- the placement-only ultra-late prefix cap `T < 5949`, `h <= 500`, with log2 contribution `-40.115431`
- the finite `32..500` row-level sum from the existing CSV artifacts
- the tiny-prefix interior-`T` finite certificate now stated as Lemma `lem:fullsplit-tiny-prefix-interiorT`
- the selected-gap product/conditioning bound now stated as Lemma `lem:fullsplit-selected-gap-product`
- the finite-prefix and global termination atoms now stated as Lemma
  `lem:fullsplit-certified-termination-atoms`
- the current first-active placement partition: ultra-late prefix, audited window, and early remainder

Still proof debt:

- decide whether to keep Lemma `lem:fullsplit-tiny-prefix-interiorT` as a finite certificate or eventually replace it by
  analytic monotonicity
- decide whether the high-interval manifest should remain command-reproducible or be regenerated into checked artifacts
- add interval-arithmetic or rational/integer safeguards for the most important numerical bounds
- cover the remaining first-active regimes: the early region `T > 17948, h > 500` and the post-prefix ultra-late cap
  beyond `h > 500`
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
