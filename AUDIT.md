# Dense+Dense Audit Packet

This packet is for reviewing the current dense+dense finite certificate work in
`C:\Users\peter\repo\permute_conv`.

The intended audit target is narrow:

- finite checkpoint: `N = 2^21`, relative distance `delta = 0.09`, `d = floor(delta N) = 188743`
- outer model: direct sum of `4096` copies of the binary `RM(4,9)` block code, i.e. local `[512,256,32]`
- inner model: full-codeword random-split dense recursive inner with block size `b = 64`, using the EBCH `[128,64]` weight distribution
- quantity being certified: a first-moment upper bound for the currently isolated dense+dense contribution, with full first-active-position coverage for `32 <= h <= N`; independent rational/outward artifacts cover `h <= 500` and every `h >= 501` post-prefix family

The current finite ledger reports

```text
log2 mu_finite <= -37.383345
log2 mu_late_plus_window <= -37.2785
log2 mu_32_500_all_first_active_positions <= -37.2785
complete rational/outward mu <= 2^-37.27
```

These values use exact RM direct-sum outer coefficients for the `h <= 500`
prefix. The first line is still dominated by the audited window
`5949 <= T <= 17948`. The second also includes the checked placement-only
ultra-late prefix `T < 5949`, through `h <= 500`; the third records that the
early `T > 17948`, `h <= 500` slice is included as well. The checked
post-prefix ultra-late table covers `501 <= h <= 380736` with
`log2 <= -545.690526`, the checked early post-prefix tables cover
`501 <= h <= N/2` with `log2 <= -380.829185` to displayed precision, and the
complement-high table covers `N/2 < h <= N`.
An independent rational/outward post-prefix certificate now replaces those
floating table values as the hardened bound: it gives a diagnostic upper
logarithm `-183.800029543464` and passes the theorem-safe threshold `2^-180`.
Adding the two hardened thresholds on a common dyadic denominator gives
`(6793*2^130+1)/2^180 <= 2^-37.27`.  This is the theorem-facing whole-range
bound; the sharper overall `-37.2785` statement remains a checked-log
diagnostic.

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
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\check_fullsplit_T_monotonicity.py --sufficient-reduction
pdflatex -interaction=nonstopmode main_permConv.tex
```

Expected key outputs:

```text
prefix_interior_audit_rows,1500
prefix_interior_audit_worst,bucket=gap_8001_12000,H=75,diff=2.19205276153e-09,T=13949
late_prefix_T_lt_5949_cap_log2,-40.115431
exact_outer_local_nonzero_terms,106
exact_outer_local_min_positive,32
exact_outer_local_first_positive,32;48;56;60;64;68;72;76;80
exact_outer_coefficient_build_exponent,4096
exact_outer_coefficient_build_multiply_steps,1
exact_outer_coefficient_build_square_steps,12
exact_outer_coefficient_build_max_coefficient_bits,663
exact_outer_coefficient_build_max_coefficient_weight,500
exact_outer_support_count,115
exact_outer_support_positive_count,114
exact_outer_support_min_positive,32
exact_outer_support_next_positive_after_min,48
exact_outer_support_zero_prefix,1--31
exact_outer_support_gap_after_min,33--47
exact_outer_support_first_positive,32;48;56;60;64;68;72;76;80
exact_outer_support_split_h,32
exact_outer_support_split_coefficient_bits,38
exact_outer_support_next_coefficient_bits,52
late_prefix_T_lt_5949_exact_outer_log2,-41.113442
late_prefix_T_lt_5949_exact_outer_split_log2,-41.113442
late_prefix_T_lt_5949_exact_outer_above_split_log2,-66.416502
prefix_32_500_e_le8_log2,-34.767174
prefix_32_500_e_le8_exact_outer_log2,-37.383345
prefix_32_500_e_le8_exact_outer_split_log2,-37.383482
prefix_32_500_e_le8_exact_outer_above_split_log2,-50.738253
prefix_32_500_e_le8_dominant_row,outer_weight=32,first_r=1,remaining_ones=31,gap_min=1,gap_max=4000,bucket_T=5949,inner_mode=csv
prefix_32_500_e_le8_dominant_term_log2,-37.385767
prefix_32_500_e_le8_dominant_total_minus_peak_bits,0.002422
prefix_32_500_e_le8_dominant_total_remainder_log2,-46.602718
prefix_32_500_e_le8_dominant_split_minus_peak_bits,0.002285
prefix_32_500_e_le8_dominant_split_remainder_log2,-46.687229
prefix_32_500_e_le8_dominant_above_split_gap_bits,13.352486
prefix_32_500_e_ge9_tail_log2,-269.335258
prefix_32_500_e_ge9_tail_exact_outer_log2,-284.004805
postprefix_501_2000_eall_log2,-182.259739
high_interval_cover_status,PASS
high_interval_cover_range,2001--1148736
high_interval_count,15
high_interval_turnoff_adjustment_bits,1.484774
high_feasible_min_T_status,PASS
high_feasible_min_T_far_cutoff_h,1148736
high_feasible_min_T_gap_1_4000_cutoff_h,636736
high_feasible_min_T_gap_4001_8000_cutoff_h,892736
high_feasible_min_T_gap_8001_12000_cutoff_h,1148736
endpoint_gap_sum_status,PASS
endpoint_gap_sum_buckets,3
endpoint_gap_sum_samples,33
endpoint_gap_sum_max_bucket_loss_bits,11.965784
endpoint_gap_sum_min_sample_slack_bits,0.000000
endpoint_gap_sum_max_sample_slack_bits,11.965784
eallratio_tail_status,PASS
eallratio_tail_samples,5
eallratio_tail_max_current_adjacent_ratio_log2,-28.992971
eallratio_tail_max_sample_slack_bits,0
eallratio_branch_status,PASS
eallratio_branch_samples,6
eallratio_branch_min_formula_slack_bits,0
eallratio_branch_max_code_delta_bits,2.84217094304e-14
t_monotonicity_sufficient_status,PASS
t_monotonicity_sufficient_rows,18
t_monotonicity_sufficient_min_slack_bits,0.154987
t_monotonicity_sufficient_worst,h=7859--20550,bucket=gap_4001_8000,T=9949--13948,lambda=0.05,rho=0.03,slack=0.154987
interval_2001_1148736_total_log2,-586.597103
complement_high_cover_status,PASS
complement_high_cover_range,1048577--2097152
complement_high_interval_count,7
high_complement_overlap_range,1048577--1148736
high_complement_overlap_count,100160
complement_high_shape_status,PASS
complement_high_shape_min_convexity_growth_bits,0.276449723567
early_postprefix_501_75000_total_log2,-380.829185
early_accelerated_endpoint_75001_350000_total_log2,-1056.187188
early_accelerated_paired_350001_1048576_total_log2,-98386.449256
early_accelerated_75001_1048576_total_log2,-1056.187188
late_postprefix_cover_status,PASS
late_postprefix_cover_range,501--380736
late_postprefix_interval_count,7
late_postprefix_feasible_cutoff_h,380736
late_postprefix_z,0.39605985943459426
late_postprefix_shape_status,PASS
late_postprefix_shape_peak_h_range,501--250001
late_postprefix_shape_left_endpoint_peaks,7
late_postprefix_shape_right_endpoint_peaks,0
late_postprefix_shape_critical_peaks,0
late_postprefix_501_380736_total_log2,-545.690526
h501_plus_checked_rows_log2,-182.259739
finite_ledger_total_log2,-37.383345
finite_ledger_margin_bits,37.383345
late_plus_window_total_log2,-37.278528
late_plus_window_margin_bits,37.278528
current_checked_ledger_total_log2,-37.278528
current_checked_ledger_margin_bits,37.278528
manifest_json,C:\Users\peter\repo\permute_conv\scripts\fullsplit_finite_ledger_manifest.json
```

The finite-ledger verifier now runs the sufficient monotonicity audit by default. Every reported middle/far bucket
slack is positive; the smallest current slack is the middle-bucket row for `7859--20550`, about `0.154987` bits.

## Ledger Pieces

The combined finite ledger now has seven checked row families, plus the high-density early diagnostics below.

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

Smoothed-outer total:

```text
log2 mu_32_500_e_le8 = -34.767174
```

Exact-outer total:

```text
log2 mu_32_500_e_le8_exact_outer = -37.383345
```

Peak:

```text
h = 32
r = 1
remaining H = 31
gap bucket = 1--4000
term_log2 = -37.38576684150053
```

With exact RM direct-sum outer coefficients, this is the dominant piece of the
whole finite ledger:

```powershell
python scripts\certify_rm_outer_prefix_exact.py --ledger-csv scripts\fullsplit_piecewise_h32_500_csv.csv --ledger-csv scripts\fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv --ledger-csv scripts\fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv --ledger-csv scripts\fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv
```

Key output:

```text
local_nonzero_terms,106
local_min_positive,32
local_first_positive,32;48;56;60;64;68;72;76;80
coefficient_build_exponent,4096
coefficient_build_multiply_steps,1
coefficient_build_square_steps,12
coefficient_build_max_coefficient_bits,663
coefficient_build_max_coefficient_weight,500
support_count,115
support_positive_count,114
support_min_positive,32
support_next_positive_after_min,48
support_zero_prefix,1--31
support_gap_after_min,33--47
support_first_positive,32;48;56;60;64;68;72;76;80
support_split_h,32
support_split_coefficient_bits,38
support_next_coefficient_bits,52
late_prefix_exact_log2,-41.113442
late_prefix_exact_split_log2,-41.113442
late_prefix_exact_above_split_log2,-66.416502
ledger_exact_outer_log2,-37.383345   # fullsplit_piecewise_h32_500_csv.csv
ledger_exact_split_log2,-37.383482   # h=32 slice of fullsplit_piecewise_h32_500_csv.csv
ledger_exact_above_split_log2,-50.738253 # h>32 support remainder
ledger_dominant_row,h=32,r=1,remaining_ones=31,gap=1--4000,T=5949,inner_mode=csv
ledger_dominant_term_log2,-37.385767
ledger_dominant_total_minus_peak_bits,0.002422
ledger_dominant_total_remainder_log2,-46.602718
ledger_dominant_split_remainder_log2,-46.687229
ledger_dominant_above_split_gap_bits,13.352486
ledger_exact_outer_log2,-284.004805  # fullsplit_turnoff_tail_h32_500_r1_64_emin9.csv
ledger_exact_outer_log2,-86.910456   # fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv
ledger_exact_outer_log2,-319.977808  # fullsplit_turnoff_tail_early_h32_500_r1_64_emin17.csv
```

Interpretation: the current live finite prefix is exact-support dominated by
`h=32`, with the next possible weight only at `h=48`. The checked window split
is `-37.383482` bits from `h=32` and `-50.738253` bits from all supported
`h>32`; the checked ultra-late split is `-41.113442` and `-66.416502`.
The dominant row inside the live prefix is now explicitly checked as
`h=32,r=1,gap=1--4000,T=5949`, with term `-37.385767`; the whole exact
window prefix is only `0.002422` bits above it, and deleting it leaves
`-46.602718`.
The exact coefficient build has one multiply step and twelve squaring steps
for the `4096=2^12` direct sum, so the support gap is a finite integer
polynomial statement rather than a floating Cauchy-envelope observation.
The old adjacent `h=33,34,...` ridge below is a useful diagnostic for the
smoothed Cauchy outer envelope, but it is not the current bottleneck once the
exact RM direct-sum support is used.

### Complete rational/outward `h <= 500` certificate

The default verifier now also runs an independent certificate for every
first-active position through `h=500`.  It uses exact RM direct-sum
coefficients and exact placement sums.  At the six knots
`5949,9949,13949,17949,22949,27949`, the EBCH split law is rational, the
termination atom is certified below `2^-63`, and the fixed-pole Chernoff
envelope (`z=2333/2373`) is rounded upward at 1024 dyadic bits.  The first
three buckets retain `e=0..8`; the early three retain `e=0..16`; exact/outward
tails cover `e>=9` and `e>=17`, respectively.  The ultra-late prefix remains
an exact rational sum.

```text
fullsplit_h500_complete_rational_status,PASS
fullsplit_h500_complete_rational_log2,-37.276548513006
fullsplit_h500_complete_rational_threshold,6793/2^50
fullsplit_h500_complete_rational_threshold_log2,-37.270166861152
fullsplit_h500_complete_rational_37_27_bits_status,PASS
```

The displayed logarithms are diagnostics.  The final gate is exact cross
multiplication against `6793/2^50`, and `6793^100 <= 2^1273` proves this
threshold is at most `2^-37.27` using integers only.  The cached exact/outward
tables are SHA-256 protected in the default pass; their complete regeneration
commands are sequential:

```powershell
python scripts\certify_fullsplit_h500_rational.py --recompute-gap-sums
python scripts\certify_fullsplit_h500_rational.py --recompute-inner-bounds
python scripts\certify_fullsplit_h500_rational.py
```

The canonical artifacts are
`scripts/fullsplit_h500_gap_sums_exact.json` and
`scripts/fullsplit_h500_inner_bounds_dyadic.json`.

### Complete rational/outward `h >= 501` certificate

Every remaining first-active-position family is independently recertified from
exact rational combinatorial inputs and directed outward logarithmic intervals.
The certificate covers the critical `501..2000` rows, ultra-late placement,
prefix cap, prefix fixed-`T` and paired-`T` episode rows, early fixed-`T` and
paired-`T` rows, and the complement-high range through `h=N`.  It checks the
global turnoff envelope against `2^-62` by exact rational arithmetic and checks
56 endpoint-in-`T` monotonicity reductions exactly.

```text
fullsplit_postprefix_complete_rational_status,PASS
fullsplit_postprefix_complete_rational_log2_upper,-183.800029543464
fullsplit_postprefix_complete_rational_threshold,2^-180
```

The logarithm is a diagnostic outward upper endpoint; the theorem-facing gate
is `2^-180`.  The default verifier authenticates the cached artifact by SHA-256.
Full regeneration is sequential and takes roughly two minutes on the reference
machine:

```powershell
python scripts\certify_fullsplit_postprefix_rational.py --recompute-artifact
```

The canonical artifact is `scripts/fullsplit_postprefix_rational.json`.

### Complete rational/outward whole-range certificate

The default verifier combines the `h<=500` threshold `6793/2^50` with the
`h>=501` threshold `2^-180` exactly:

```text
fullsplit_complete_rational_status,PASS
fullsplit_complete_rational_threshold,9246152473975739929226814833136005841682433/2^180
fullsplit_complete_rational_threshold_factored,(6793*2^130+1)/2^180
fullsplit_complete_rational_37_27_bits_status,PASS
```

Writing `K=6793*2^130+1`, the final comparison is the integer inequality
`K^100 <= 2^14273`.  Hence the complete first moment over every positive outer
weight is at most `2^-37.27`; no floating logarithm is used in this final gate.

Smoothed-outer ridge decomposition:

```powershell
python scripts\analyze_fullsplit_prefix_ridge.py --expected-total-log2 -34.767174286
```

Key output:

```text
prefix_ridge_total_log2,-34.767174286
prefix_ridge_peak,h=32,r=1,gap=1-4000,term_log2=-37.385766842
prefix_ridge_peak_to_total_bits,2.618592555
gap,1-4000,log2,-34.767174286,share,1
r,1,log2,-34.769818553,share,0.998168812772
h_le,64,log2,-34.770777320,share,0.997505683424
h_le,80,log2,-34.767370670,share,0.999863886476
h_le,96,log2,-34.767184998,share,0.999992575315
max_ratio_all,0.852699299925,from,32,to,33
max_ratio_after_33,0.833794827634,from,33,to,34
prefix_ridge_first_outer_ratio,2.8086853907
prefix_ridge_first_placement_ratio,0.303593739174
prefix_ridge_first_inner_ratio,1
prefix_ridge_tail_outer_ratio,2.74641945768
prefix_ridge_tail_placement_ratio,0.303593402422
prefix_ridge_tail_inner_ratio,1
ridge_geometric_infinite_log2,-34.769785229
ridge_total_ratio_bound_log2,-34.767141024
ridge_total_ratio_bound_slack_bits,3.32625484987e-05
ridge_inner_log2_span,3.14770431942e-12
```

Interpretation for the smoothed Cauchy outer envelope: the narrow waist is a
near-boundary geometric `h` ridge, not just one bad `h=32` row and not a broad
intermediate-weight ridge. The verifier still checks this ratio skeleton as an
audit of the smoothed table: the infinite geometric ridge bound plus the exact
smoothed-table remainder gives log2 contribution `-34.767141`, only `3.4e-5`
bits above the smoothed `e <= 8` prefix table. For the dominant `r=1`,
first-gap ridge, the inner row is constant to numerical precision; the ratio is
the product of an outer-spectrum slope and the placement slope

```text
P_{h+1}/P_h = (S_h/S_{h-1}) * (h+1)/(N-h)
S_H = sum_{g=1}^{4000} binom(64*(5949+g-1), H).
```

This smoothed-ridge target has been superseded for the RM finite checkpoint by
the exact-support outer certificate above. It remains useful if we later replace
RM by a smoother outer ensemble or need a Cauchy-envelope fallback.

The placement slope is now checked independently of the row CSV by exact
integer arithmetic:

```powershell
python scripts\certify_prefix_placement_ratio.py
```

Output:

```text
placement_ratio_threshold,0.303594
placement_ratio_threshold_num,151797
placement_ratio_threshold_den,500000
placement_ratio_exact_cross_multiply_status,PASS
placement_ratio_exact_comparisons,468
placement_ratio_peak_h,32
placement_ratio_peak_to_h,33
placement_ratio_peak,0.303593738594
placement_ratio_peak_log2,-1.71978605839
placement_ratio_peak_num_bits,513
placement_ratio_peak_den_bits,515
placement_ratio_threshold_slack_bits,512
placement_ratio_threshold_gap_log2,-21.8672036595
endpoint_bound_at_h_min,0.3130655524
endpoint_bound_at_h_min_log2,-1.67546332201
endpoint_bound_slack_factor,1.03119897614
```

The verifier runs the same cross-multiplication check by default and prints:

```text
prefix_placement_ratio_threshold,0.303594
prefix_placement_ratio_threshold_num,151797
prefix_placement_ratio_threshold_den,500000
prefix_placement_ratio_exact_cross_multiply_status,PASS
prefix_placement_ratio_exact_comparisons,468
prefix_placement_ratio_peak_h,32
prefix_placement_ratio_peak_to_h,33
prefix_placement_ratio_peak,0.303593738594
prefix_placement_ratio_peak_log2,-1.71978605839
prefix_placement_ratio_peak_num_bits,513
prefix_placement_ratio_peak_den_bits,515
prefix_placement_ratio_threshold_slack_bits,512
prefix_placement_ratio_threshold_gap_log2,-21.8672036595
prefix_placement_endpoint_bound_at_h_min,0.3130655524
prefix_placement_endpoint_bound_at_h_min_log2,-1.67546332201
prefix_placement_endpoint_slack_factor,1.03119897614
```

Interpretation: endpoint domination is too loose by about `3.12%`, so the
eventual analytic proof must exploit the averaged binomial sum `S_H`; however,
the finite placement component no longer depends on floating log-binomial rows.
The threshold itself is the rational `151797/500000`, and the exact peak ratio
is below it by about `2^-21.8672`.

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
default before summing the prefix rows. In particular it verifies the finite
rectangle structure, not just the row count:

```text
prefix_interior_audit_rows,1500
prefix_interior_audit_h_range,0--499
prefix_interior_audit_buckets,gap_1_4000:5949--9948:4000;gap_4001_8000:9949--13948:4000;gap_8001_12000:13949--17948:4000
prefix_interior_audit_left_endpoint_maxima,PASS
prefix_interior_audit_tolerance_bits,1e-07
prefix_interior_audit_worst_slack_bits,9.78079472385e-08
prefix_interior_audit_reproduction_tolerance_bits,5e-06
prefix_interior_audit_worst_repro_slack_bits,4.99447521759e-06
```

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

Smoothed-outer total:

```text
log2 mu_32_500_e_ge9_tail = -269.335258
```

Exact-outer total:

```text
log2 mu_32_500_e_ge9_tail_exact_outer = -284.004805
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

Smoothed-outer total:

```text
log2 mu_32_500_early_e_le16 = -85.403338
```

Exact-outer total:

```text
log2 mu_32_500_early_e_le16_exact_outer = -86.910456
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

Smoothed-outer total:

```text
log2 mu_32_500_early_e_ge17_tail = -305.738816
```

Exact-outer total:

```text
log2 mu_32_500_early_e_ge17_tail_exact_outer = -319.977808
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

### 5. Ultra-Late Post-Prefix Placement Rows

Range:

```text
501 <= h <= 380736
T < 5949
```

Here `380736 = 64 * 5949`; for larger `h`, placing all nonzero coordinates in the final `5949` inner blocks is
impossible. The bound is placement-only:

```text
A_h^out * binom(64*5949, h) / binom(N, h)
```

For a fixed outer Cauchy pole `z`, the logarithm is concave in `h`, with adjacent ratio
`((m-h)/(N-h))/z`. Thus each interval is certified by checking endpoints and the critical point
`h* = (m - zN)/(1-z)`.

Current interval rows, all with `z = 0.39605985943459426`:

```text
501--2000:      -545.690526
2001--7858:     -2237.593414
7859--20550:    -8919.160723
20551--75000:   -23772.761808
75001--150000:  -93856.032305
150001--250000: -210536.528296
250001--380736: -417970.896784
```

Combined total:

```text
log2 mu_late_postprefix_501_380736 = -545.690526
```

The finite-ledger verifier also checks the structural cover:

```text
late_postprefix_cover_status,PASS
late_postprefix_cover_range,501--380736
late_postprefix_interval_count,7
late_postprefix_feasible_cutoff_h,380736
late_postprefix_z,0.39605985943459426
late_postprefix_shape_status,PASS
late_postprefix_shape_interval_count,7
late_postprefix_shape_peak_h_range,501--250001
late_postprefix_shape_left_endpoint_peaks,7
late_postprefix_shape_right_endpoint_peaks,0
late_postprefix_shape_critical_peaks,0
```

Print the recomputation commands with:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --print-late-postprefix-commands
```

Spot-check selected rows with:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --check-late-postprefix-intervals 501--2000
python scripts\verify_fullsplit_finite_ledger.py --check-late-postprefix-intervals 250001--380736
```

Audit priority: medium. This closes the remaining first-active placement gap. The finite-ledger verifier now checks
the fixed-pole concavity reduction directly: for the current pole all seven row maxima occur at the left endpoint.

### 6. Early Post-Prefix Rows

Range:

```text
501 <= h <= N/2
T > 17948
gap buckets: 12001--17000, 17001--22000, 22001--26819
```

Current legacy endpoint rows:

```text
501--2000:    -380.829185
2001--7858:   -852.998772
7859--20550:  -3712.365087
20551--30000: -8614.382763
30001--50000: -21332.131033
50001--75000: -34847.929043
```

Current accelerated rows:

```text
75001--90000:      -28574.869802   endpoint recurrence
90001--100000:     -22149.918612   endpoint recurrence
100001--110000:    -1056.187188    endpoint recurrence
110001--125000:    -10362.025949   endpoint recurrence
125001--160000:    -45589.980049   endpoint recurrence
160001--250000:    -60293.751685   endpoint recurrence
250001--350000:    -94682.265141   endpoint recurrence
350001--524288:    -181215.853042  paired-T recurrence
524289--750000:    -194565.849984  paired-T recurrence
750001--1048576:   -98386.449256   paired-T recurrence
```

The finite-ledger verifier now checks this accelerated table as a contiguous
cover before summing it:

```text
early_accelerated_cover_status,PASS
early_accelerated_cover_range,75001--1048576
early_accelerated_interval_count,10
early_accelerated_endpoint_cover,75001--350000
early_accelerated_endpoint_interval_count,7
early_accelerated_paired_cover,350001--1048576
early_accelerated_paired_interval_count,3
```

It also records each row's helper, pole parameters, log contribution, and
recompute command in `fullsplit_finite_ledger_manifest.json`.

Combined totals:

```text
log2 mu_early_postprefix_501_75000 = -380.829185
log2 mu_early_accelerated_75001_1048576 = -1056.187188
log2 mu_early_postprefix_501_1048576 = -380.829185  # to displayed precision
```

The verifier records the first table as `early_postprefix_interval_*` and the second as
`early_accelerated_interval_*`, then includes both log-sums in `h501_plus_checked_rows_log2`. These do not move the
global checkpoint because the late-window `501 <= h <= 2000` row at `-182.259739` is larger.

Print the recomputation commands with:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --print-early-postprefix-commands
python scripts\verify_fullsplit_finite_ledger.py --print-early-accelerated-commands
```

Spot-check selected rows with:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --check-early-postprefix-intervals 501--2000
python scripts\verify_fullsplit_finite_ledger.py --check-early-postprefix-intervals 2001--7858
python scripts\verify_fullsplit_finite_ledger.py --check-early-accelerated-intervals 100001--110000
python scripts\verify_fullsplit_finite_ledger.py --check-early-accelerated-intervals 750001--1048576
```

Audit priority: medium-high. These are finite endpoint-placement rows using the all-episode ratio wrapper. The
accelerated endpoint helper is a vectorized recurrence version of the same row bound, while the paired helper keeps the
placement length `T` attached to the survival cost. The main audit task is now to inspect these two accelerators and add
interval-arithmetic or rational safeguards where needed.

### 7. Early High-Density Probe, Paired T

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
`--check-complement-high-intervals all` recomputes all seven rows and passes. The manifest also checks:

```text
complement_high_cover_status,PASS
complement_high_cover_range,1048577--2097152
complement_high_interval_count,7
high_complement_overlap_range,1048577--1148736
high_complement_overlap_count,100160
complement_high_shape_status,PASS
complement_high_shape_min_convexity_growth_bits,0.276449723567
```

Current interpretation: early high-density has ample numerical margin once `T` pairing, outer complement symmetry,
and the fixed-`z` interval endpoint bound are used. The shape audit checks that the complement-side outer/volume/rho
term is endpoint-dominated by discrete convexity on every interval. The overlap with the older `2001--1148736` table
is deliberate union-bound overcount, not a gap in the cover.

### 8. Post-Prefix, All-Episode Wrapper

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

### 9. Interval Certificate

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

The finite-ledger verifier checks the high interval family as a contiguous cover and records the global turnoff
adjustment attached to every row:

```text
high_interval_cover_status,PASS
high_interval_cover_range,2001--1148736
high_interval_count,15
high_interval_turnoff_adjustment_bits,1.484774
high_feasible_min_T_status,PASS
high_feasible_min_T_far_cutoff_h,1148736
high_feasible_min_T_gap_1_4000_cutoff_h,636736
high_feasible_min_T_gap_4001_8000_cutoff_h,892736
high_feasible_min_T_gap_8001_12000_cutoff_h,1148736
endpoint_gap_sum_status,PASS
endpoint_gap_sum_buckets,3
endpoint_gap_sum_samples,33
endpoint_gap_sum_max_bucket_loss_bits,11.965784
endpoint_gap_sum_min_sample_slack_bits,0.000000
endpoint_gap_sum_max_sample_slack_bits,11.965784
eallratio_tail_status,PASS
eallratio_tail_samples,5
eallratio_tail_max_current_adjacent_ratio_log2,-28.992971
eallratio_tail_max_sample_slack_bits,0
eallratio_branch_status,PASS
eallratio_branch_samples,6
eallratio_branch_min_formula_slack_bits,0
eallratio_branch_max_code_delta_bits,2.84217094304e-14
t_monotonicity_sufficient_status,PASS
t_monotonicity_sufficient_rows,18
t_monotonicity_sufficient_min_slack_bits,0.154987
t_monotonicity_sufficient_worst,h=7859--20550,bucket=gap_4001_8000,T=9949--13948,lambda=0.05,rho=0.03,slack=0.154987
```

Audit priority: medium-low. The endpoint gap-sum dominance, the `T_eff=max(T_min,ceil(H/b))` cutoff, the safe
`eallratio` multiplier-tail truncation, independent finite samples of the fixed-pole branch formula, and the sufficient
monotonicity reduction for the fixed-pole all-episode wrapper are now checked in the main ledger. The manuscript now
states the finite first-moment conclusion as Theorem `thm:fullsplit-rm-finite-certificate-009`; remaining work here is
mainly independent review, interval-arithmetic hardening, and deciding which finite row helpers should be promoted into
smaller reproducible artifacts.

### `eallratio` Formula Map

The production implementation is `eall_ratio_log2` in `scripts/sum_fullsplit_piecewise_certificate.py`. For `H > 0`,
it matches Lemma `lem:fullsplit-fixed-pole-eallratio` term-by-term:

- `denom = log2_binom(b * T, H)` is the placement denominator `log2 binom(bT,H)`.
- `M = sum(p * exp(-lam*j) for j,q,p in entries if q > 0)` is `M_+(lambda)`.
- `pref = lam * distance / log(2)` is the base-2 form of `exp(lambda d)`.
- `e0 = pref + T * log2(M)` is `log2 E_0(lambda)`.
- `log_block = log2_expm1_pos(b * log1p(rho))` is `log2(((1+rho)^b)-1)`.
- `log_g = log_block + log2(M/(1-M))` is `log2 G_{lambda,rho}` for the Cauchy/suffix factor.
- `log2_arith_geom_sum(log_g, ceil(H/b), min(H,T))` is
  `log2 sum_x (x+1) G_{lambda,rho}^x`.
- `e1 = term_log2 + pref - denom - H*log2(rho) + ...` is `log2 B_1(lambda,rho)`.
- `tail_log2 = log2_ht_tail_envelope(H,T,term_log2 + log2(1-M))` is
  `log2 R_{H,T}(lambda)`.
- `tail_factor_log2 = log2_one_plus_pow2(tail_log2)` is `log2(1+R_{H,T}(lambda))`.
- `e_ge1 = e1 + tail_factor_log2` is `log2(B_1(lambda,rho)(1+R_{H,T}(lambda)))`.
- The loop computes `tail_log2` inside the fixed-`lambda` loop, so the `e=1` Cauchy envelope and the `e>=2`
  multiplier use the same lambda. The `e=0` branch is optimized separately, which is safe because the branches are
  summed after each has been independently upper-bounded.
- `return min(0, log2add(best_e0, best_tail))` is the final probability cap by `1`.

The verifier audits two implementation-sensitive pieces: `eallratio_tail_semantics` checks the positive multiplier
series and its geometric omitted-tail bound against exact finite samples, while `eallratio_branch_semantics` recomputes
the fixed-pole branch from exact nonempty-block coefficients and exact skipped-gap suffix sums on six samples, then
compares the production singleton-lambda/rho output against an independent formula evaluator.

## Proof Objects To Inspect

Core scripts:

```text
scripts/verify_fullsplit_finite_ledger.py
scripts/sum_fullsplit_piecewise_certificate.py
scripts/certify_fullsplit_late_prefix_interval.py
scripts/certify_fullsplit_early_endpoint_interval.py
scripts/certify_fullsplit_early_paired_interval.py
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
scripts/fullsplit_finite_ledger_manifest.json
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
- the canonical finite-ledger JSON snapshot `scripts/fullsplit_finite_ledger_manifest.json`, regenerated by
  `verify_fullsplit_finite_ledger.py --write-manifest-json`, with schema `fullsplit_finite_ledger.v1` and 56 row-family
  entries
- the `501..2000` all-episode wrapper row
- the high-`h` fixed-pole interval rows as numerical certificates, with printed regeneration commands and selected
  opt-in recomputation checks
- the endpoint gap-sum semantics used by the high-`h` piecewise wrapper, including boundary exact-sum samples and the
  `log2(4000)` bucket-size loss
- the `eallratio` multiplier-tail evaluation, with exact finite samples and a global decreasing-ratio bound for the
  current feasible `T <= 32767` range
- the `eallratio` fixed-pole branch formula on six exact finite selected-gap Chernoff samples, including an implementation
  comparison against an independent formula evaluator
- the sufficient endpoint monotonicity inequalities for the listed intervals
- the placement-only ultra-late prefix cap `T < 5949`, `h <= 500`, with log2 contribution `-40.115431`
- the placement-only ultra-late post-prefix intervals `T < 5949`, `501 <= h <= 380736`, with log2 contribution
  `-545.690526`
- the finite `32..500` row-level sum from the existing CSV artifacts
- the exact RM direct-sum outer prefix coefficients and exact-outer reweighting in
  `certify_rm_outer_prefix_exact.py`
- the paper-facing finite statement `lem:rm-prefix-support-split`, which packages the exact support gap and
  `h=32`/`h>32` split checked by the verifier
- the paper-facing small-weight prefix corollary `cor:fullsplit-small-prefix-all-first-active`, which combines the
  ultra-late, window, and early `h <= 500` pieces
- the paper-facing checked post-prefix row corollary `cor:fullsplit-postprefix-checked-ledger`, which combines the
  `501..2000`, ultra-late `501..380736`, early `501..N/2`, fixed-pole interval, and complement-high ledger rows
- the paper-facing current finite checkpoint `cor:fullsplit-current-finite-checkpoint`, which combines the checked
  small-prefix and post-prefix row families
- the paper-facing finite first-moment theorem `thm:fullsplit-rm-finite-certificate-009`, which states
  `E[Z_d] <= (6793*2^130+1)/2^180 <= 2^-37.27` for `N=2^21`, `delta=.09`, and
  `d=floor(.09 N)=188743`, hence the probability of noninjectivity on the outer
  code or `d_min <= d` is at most `2^-37.27`.  Therefore an actual binary
  `[2^21,2^20,d_min >= 188744]` code exists.  The verifier's
  `-37.278528` line remains a sharper checked-log diagnostic for the unrounded
  `-37.278527626...`, not the direction-sensitive theorem exponent.
- the tiny-prefix ridge decomposition in `analyze_fullsplit_prefix_ridge.py`, and its ratio skeleton in
  `verify_fullsplit_finite_ledger.py`
- the first-gap placement slope in `certify_prefix_placement_ratio.py`, checked by exact integer cross multiplication
- the tiny-prefix interior-`T` finite certificate now stated as Lemma `lem:fullsplit-tiny-prefix-interiorT`
- the selected-gap product/conditioning bound now stated as Lemma `lem:fullsplit-selected-gap-product`
- the finite-prefix and global termination atoms now stated as Lemma
  `lem:fullsplit-certified-termination-atoms`
- the current first-active placement partition: ultra-late prefix, audited window, and checked early rows
- the ultra-late post-prefix placement rows, with printed recomputation commands and selected opt-in recomputation
  checks
- the accelerated early endpoint and paired interval rows through `h = N/2`, with printed recomputation commands and
  selected opt-in recomputation checks
- the outer-mode projection driver `scripts/compare_outer_modes_fullsplit.py`, which leaves the full-split EBCH inner
  ledger fixed, reproduces the proved RM `[512,256,32]` checkpoint at `-37.278528`, and writes the generated comparison
  reports `scripts/outer_mode_comparison_delta009.csv` and `scripts/outer_mode_comparison_delta009.md`

Still proof debt:

- decide whether to keep Lemma `lem:fullsplit-tiny-prefix-interiorT` as a finite certificate or eventually replace it by
  analytic monotonicity
- replace the finite exact-support tiny-prefix decomposition by an analytic or clean finite lemma: exact `h=32` peak,
  exact RM support gap to `h=48`, and a support-weight remainder bound for `h>32`
- decide which command-reproducible rows should eventually be regenerated into smaller checked artifacts, versus kept as
  interval-helper commands inside the JSON manifest
- before any BCH-like outer row becomes theorem text, replace the random-like spectrum projection with an exact local
  spectrum or a proved low-weight envelope, and recompute the `h>2000` tail under that same outer model

## Do Not Trust Yet

The following are useful context but should not be treated as proof:

- old scalar fixed-tap dense-inner conclusions
- old random-banded outer `sigma` plots
- BCH/RM replacement speculation outside the RM `[512,256,32]` checkpoint, including the BCH rows in
  `outer_mode_comparison_delta009.*`, until a real spectrum or rigorous low-weight envelope replaces the modeled
  spectrum
- fit-based residual laws
- any `tmp_*` CSV or PNG unless it is explicitly named above
- fixed-`h` asymptotic heuristics that are not connected to the first-moment sum

## Recommended Audit Order

1. Run the three main verification commands.
2. Check the exact formulas in `fast_fullsplit_episode_e01.py`, especially the selected-gap count and endpoint
   `e=x+1` case.
3. Check Lemma `fullsplit-fixed-pole-eallratio` against `sum_fullsplit_piecewise_certificate.py`, using the finite branch
   samples as implementation checks rather than as the proof itself.
4. Check `check_fullsplit_T_monotonicity.py` against the monotonicity argument in `innerDense.tex`.
5. Check that the `32..500` finite CSV inputs are generated from formulas that are conservative upper bounds.
6. Only then read the surrounding prose and decide what should become theorem text versus working-save-point text.

## Proof Cleanup Audit (2026-08-07)

- The generic framework no longer assumes that an `n x n` inner map is
  invertible.  Its first moment bounds the joint noninjectivity-or-low-distance
  event, and the finite theorem consequently certifies dimension as well as
  distance.
- The full-split definition fixes `S_0=0`, discards unconstrained `S_B`, and
  states all independence in the random scramblers, coordinate splits, and
  interleaver.  The turnoff atom is conditioned on block occupancies and prior
  transition weights, under which the state and input supports are independent
  uniform subsets.
- The theorem-facing certificate is now factored into named complete
  `h<=500` and `h>=501` lemmas.  Legacy row diagnostics remain in this audit
  file and the manifest rather than interrupting the main proof.
- Both source spectrum tables are structurally validated and SHA-256-bound to
  the rational artifacts.  This authenticates the exact inputs used but does
  not derive either table from a code definition; the EBCH table's
  mathematical provenance in particular remains declared.
