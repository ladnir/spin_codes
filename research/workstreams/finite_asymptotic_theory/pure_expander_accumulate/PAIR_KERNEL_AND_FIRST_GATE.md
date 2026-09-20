# Pair kernel and the first finite gate

## Status

The one-stage one-word and ordered-pair laws below are proved identities for
the independent right-degree ensemble. Four small instances agree exactly
with exhaustive enumeration.

The length-512 extreme-shell results are binary64 diagnostics. They select
candidate degrees but do not certify the inequalities. No central-shell
concentration result or SPIN distance theorem is claimed.

## Ensemble under test

Let \(E:\mathbb F_2^k\to\mathbb F_2^n\). Independently for each output
coordinate \(j\), setup samples a uniform \(r\)-subset
\(S_j\subseteq\{1,\ldots,k\}\) and sets

\[
 (Ex)_j:=\bigoplus_{i\in S_j}x_i.
\]

Thus every right vertex has degree \(r\), while left degrees remain random.
The map requires exactly \(n(r-1)\) XORs when every output parity is evaluated
independently. The zero-initialized length-\(n\) accumulator requires at most
\(n-1\) additional XORs.

The expander-code paper uses the same sparse-map--recursive-map architecture.
Its Bernoulli EA ensemble samples every matrix entry independently. The
present ensemble instead fixes each output degree. This change fixes the XOR
count while retaining independence among output neighborhoods.

## Exact one-word law

Fix an input \(x\) of weight \(h\). Define

\[
 \beta_r(h):=\frac{K_r(h)}{\binom kr},
 \qquad
 q_r(h):=\frac{1-\beta_r(h)}2.
\]

The character \((-1)^{(Ex)_j}\) has expectation \(\beta_r(h)\). Therefore
the coordinates of \(Ex\) are independent Bernoulli variables with parameter
\(q_r(h)\).

Let \(Y\) mark the accumulated output weight. The accumulator transfer is

\[
 T_q(Y)=
 \begin{pmatrix}
 1-q&qY\\
 q&(1-q)Y
 \end{pmatrix}.
\]

Consequently, the expected shell count after one stage is

\[
 \mu_w
 =\sum_{h=1}^{k}\binom kh
   [Y^w]\,e_0^{\mathsf T}T_{q_r(h)}(Y)^n\mathbf1.
 \tag{1}
\]

Equation (1) includes messages in the kernel at \(w=0\).

## Exact ordered-pair law

Fix distinct nonzero inputs \(x,y\). Put

\[
 h_1=\operatorname{wt}(x),\qquad
 h_2=\operatorname{wt}(y),\qquad
 h_3=\operatorname{wt}(x+y).
\]

For \(a,b\in\mathbb F_2\), Fourier inversion gives

\[
 p_{a,b}(h_1,h_2,h_3)
 =\Pr[(Ex)_j=a,(Ey)_j=b]
\]

and

\[
 p_{a,b}
 =\frac14\left(
 1+(-1)^a\beta_r(h_1)
  +(-1)^b\beta_r(h_2)
  +(-1)^{a+b}\beta_r(h_3)
 \right).
 \tag{2}
\]

The output neighborhoods are independent. Hence the intermediate pair
symbols are independent samples from (2).

Index the common accumulated state by \(s\in\mathbb F_2^2\). Let \(X\) and
\(Y\) mark the two accumulated output weights. The exact pair transfer is

\[
 M_p(X,Y)_{s,s'}
 :=p_{s+s'}X^{s'_1}Y^{s'_2}.
 \tag{3}
\]

For a message-pair type
\(m=(m_{00},m_{01},m_{10},m_{11})\), the number of ordered input pairs is

\[
 \binom{k}{m_{00},m_{01},m_{10},m_{11}}.
\]

Exclude the types for which \(x=0\), \(y=0\), or \(x=y\). The second
factorial shell moment is then

\[
 F_w
 =\sum_m
 \binom{k}{m_{00},m_{01},m_{10},m_{11}}
 [X^wY^w]\,e_{00}^{\mathsf T}M_{p(m)}(X,Y)^n\mathbf1.
 \tag{4}
\]

Therefore

\[
 \operatorname{Var}(A_w)=\mu_w+F_w-\mu_w^2.
 \tag{5}
\]

Equations (2)--(5) are the central simplification of the pure route. The
base pair-type multiplicities are explicit multinomial coefficients. No BCH
genus-two enumerator occurs.

## Exact small checks

The script verify_pair_kernel_small.py computes (1) and (4) with rational
arithmetic. It also enumerates every sparse map and every nonzero message.
The analytic and exhaustive means and second factorial moments agree in all
four retained receipts:

| \(k\) | \(n\) | \(r\) | maps enumerated |
|---:|---:|---:|---:|
| 3 | 4 | 1 | 81 |
| 3 | 4 | 2 | 81 |
| 4 | 5 | 2 | 7,776 |
| 4 | 5 | 3 | 1,024 |

These checks validate the kernel implementation. They do not extrapolate a
variance bound to \(k=256,n=512\).

The exact multistage propagation also remains well behaved in the first
small models. For \((k,n)=(4,8)\) and degrees \((3,3,3)\), the worst positive
shell variance-to-mean ratios after stages one through three are

\[
 1.2657,\qquad1.4028,\qquad1.4758.
\]

Changing the degrees to \((3,2,3)\) gives

\[
 1.2657,\qquad1.5689,\qquad1.6934.
\]

For \((k,n)=(5,10)\) and degrees \((3,2,3)\), the corresponding values are

\[
 1.0392,\qquad1.2419,\qquad1.2996.
\]

These exact values support a moderate variance target. Their lengths and
relative degrees differ too much from the target to predict the length-512
constant.

## Extreme-shell gate at \([512,256]\)

The existing 10.9% SPIN cap interface sets every cap outside weights
42 through 470 to zero. The first diagnostic therefore bounds the expected
number of nonzero messages with output weight at most 41 or at least 471.

Setup may sample a bounded number of candidate composites and test their
rank by Gaussian elimination. For one unconditioned draw, the expected number
of nonzero kernel messages upper-bounds rank failure. If this bound is
\(\kappa<1\), then

\[
 \Pr[\text{bad nonzero extreme}\mid\text{full rank}]
 \le\frac{\mathbb E[Z_{\mathrm{extreme},\ne0}]}{1-\kappa}.
 \tag{6}
\]

With \(R\) independent attempts, setup abort is at most \(\kappa^R\). The
diagnostics below use \(R=16\). This is an efficient rank test, not a spectrum
acceptance test.

### One stage

The first odd right degree that exceeds the 42.2615-bit comparison target is

\[
 r_0=33.
\]

It has 44.6379 diagnostic bits and costs 16,895 XORs per constituent.
Degree 31 has only 41.1969 bits. The weight-one input messages dominate both
values. Its expected left degree is 66. This is consistent with the separate
rate-half regular-EA certificate in the expander paper, which requires left
degree 64 for a direct near-GV distance theorem.

### Two stages

For \(C=A E_1 A E_0\), the cheapest passing pair in the tested range is

\[
 (r_0,r_1)=(17,3).
\]

It costs 10,238 XORs per constituent. Its diagnostic quantities are

\[
 \mathbb E[|\ker C\setminus\{0\}|]
 \approx1.5453\cdot10^{-13},
\]

and

\[
 \mathbb E[Z_{\mathrm{extreme},\ne0}]
 \approx5.8830\cdot10^{-14}.
\]

Equation (6), including bounded rank-test abort, gives 43.9504 diagnostic
bits.

### Three stages

The minimum tested degree sum that passes is

\[
 (r_0,r_1,r_2)=(13,3,3).
\]

It costs 9,725 XORs and gives 42.4540 diagnostic bits. This point has little
room beyond the existing cap-event target.

A more robust point at almost the two-stage cost is

\[
 (r_0,r_1,r_2)=(15,2,3).
\]

It costs 10,237 XORs. Its kernel first moment is approximately
\(1.4137\cdot10^{-11}\). Sixteen rank attempts make the resulting abort term
negligible. The nonzero extreme-shell first moment is approximately
\(9.6731\cdot10^{-15}\), which gives 46.5549 diagnostic bits after (6).

## Interpretation

Fresh sparse mixing substantially reduces the degree needed by the
extreme-shell gate. A second stage reduces the XOR count by about 39% relative
to one stage. A third stage gives similar cost with a larger selectable
failure margin.

This result does not yet answer the main question. The proof still needs a
simultaneous high-probability cap for the central shells. Equations (2)--(5)
provide the exact second-moment interface for that work. The current proof
attack starts with the one-stage degree-33 candidate because it avoids
interstage pair-kernel composition. Its exact dual reduction is in
ONE_STAGE_VARIANCE_ROUTE.md. The target remains

\[
 \operatorname{Var}(A_w)\le512\,\mathbb E[A_w]
\]

through a scalable covariance bound. If the one-stage mechanism closes, it
should then be transferred to the three-stage candidates.

## Artifacts

- verify_pair_kernel_small.py, pair_kernel_K3_B4_r1_exact.json,
  pair_kernel_K3_B4_r2_exact.json, pair_kernel_K4_B5_r2_exact.json, and
  pair_kernel_K4_B5_r3_exact.json contain the exact small checks.
- analyze_multistage_pair_moments_small.py and the
  multistage_pair_moments_*_exact.json receipts contain exact small-chain
  variance calculations.
- sweep_one_stage_extreme_shells.py and
  one_stage_extreme_shell_sweep_B512.json contain the one-stage diagnostic.
- sweep_two_stage_extreme_shells.py and
  two_stage_extreme_shell_sweep_B512.json contain the two-stage diagnostic.
- sweep_three_stage_extreme_shells.py and
  three_stage_extreme_shell_sweep_B512.json contain the three-stage
  diagnostic.
- ONE_STAGE_VARIANCE_ROUTE.md records the exact dual covariance reduction.
- probe_one_stage_variance_scaling.py and the one_stage_variance_*_probe.json
  receipts contain moderate-size nondirected floating-point diagnostics.
- analyze_dual_walk_and_accumulator_energy.py and
  dual_walk_accumulator_energy_K256_B512_r33.json contain exact integer row-
  walk and accumulator-energy data, with binary64 logarithms.
- verify_subset_covariance_blocks_small.py and the subset_covariance_blocks_*
  receipts validate the exact row-permutation block formula at small sizes.
- probe_subset_block_norm_bound_small.py and the subset_block_norm_bound_*
  receipts measure the loss of the proposed weight-level operator bound.
- verify_radial_triple_factorization_small.py and
  radial_triple_factorization_K4_r3_m8_exact.json validate the radial
  rank-257 moment factorization.
- probe_radial_block_bound.py and the radial_block_bound_* receipts show that
  the weight-level relaxation fails at moderate sizes even though the exact
  variance remains small.
- verify_sector_projection_energy_small.py and the
  sector_projection_energy_* receipts validate the exact Johnson-projector
  energy formula and the sector-norm bound at small sizes.
