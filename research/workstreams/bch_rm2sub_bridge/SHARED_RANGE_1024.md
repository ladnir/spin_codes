# Verified coverage through 1,024 occupied rows

The fixed construction is BCH [256,128] with RM2Sub t128_s15, as specified
in `README.md`. The output-weight cutoff remains H = 209716. All
probabilities concern the fresh-independent-multiplier setup shared by every
message. Let Z_Q count bad messages with Q nonzero outer rows.

The new bounds cover every integer Q = 129,...,1024. Combined with the
retained certificates, they give

\[
\sum_{Q=129}^{1024}\mathbb E[Z_Q]<2^{-77},\qquad
\Pr\left[\sum_{Q=1}^{1024}Z_Q>0\right]<2^{-49}.
\]

The added range has diagnostic margin 77.0174168827 bits. Occupancies
1025,...,8192 remain open. Thus the full-range goal is still active.

## Separate probabilities for the weight groups

The thirteen groups are those of `occupation_three.py`. The maximum-density
argument in `GENERAL_OCCUPANCIES.md` applies to any partition, so the same
adaptive proof applies with thirteen groups and assignment factor 13^Q.
The shell caps are the exact strength-29 bounds from
`OA29_SHELL_CAPS_AND_TIGHTENING.md`.

Each group has its own Bernoulli probability. A pure-group optimization
chooses that probability as a discovery heuristic. The certificate then
uses an adaptive maximum over all thirteen groups at every entry and step.
It does not infer a bound on mixed groups from pure-group checks.

Four shared witnesses were selected at occupancies 128, 256, 512, and 1024.
For each witness, one recurrence supplies bounds at every intermediate
depth. This does not interpolate between endpoints: the matrix at depth q
depends only on region coefficients R_0,...,R_q and is explicitly computed.

Those four witnesses missed the target only at Q=339. A separate witness
at that occupancy has diagnostic margin 7972.7681 bits. The final ledger
takes the smaller certified bound at Q=339 and sums all other rows from
the shared batch. The unsuccessful shared bound at Q=339 is preserved.

## Outward calculation

`polynomial_regions.py` computes the region polynomial by binary powering
of 3-by-3 polynomial matrices. Every multiplication uses Arb enclosures.
Truncation after degree Q cannot affect any required coefficient. Taking
an upper endpoint after division by binom(8192,j) gives an entrywise upper
bound for R_j. Exact tests check every coefficient on a small instance.

`shared_range_certificate.py` converts the region and Bernoulli coefficients
to upward-rounded binary64 values. It includes all thirteen groups; the
search's floating hull pruning is not part of the certificate.
Every scalar multiplication, addition, and rescaling rounds upward.
Each depth uses a shared binary exponent to avoid growth or decay of the
largest mantissa. The exponent is an integer, not an accumulated floating
logarithm. Positive underflow rounds to a positive subnormal upper bound.

For each occupancy, Arb raises the resulting matrix to the 256th power.
It also evaluates the Chernoff correction and the exact binary scaling.
The combinatorial factors and final sums use rational arithmetic.

The producer uses 256-bit Arb; replay uses 512-bit Arb for the region,
density, and final-moment calculations and repeats the directed recurrence.
Every replayed occupancy bound is no larger than its stored bound.
This replay shares the implementation. Separate exact tests check the
directed recurrence at every depth and its handling of positive underflow.

The two receipts are:

- `generated/oa29_shared_q129_q1024_outward.json`;
- `generated/oa29_gap339_outward.json`.

Read-only replay and aggregation commands are:

```powershell
python -B workstreams/bch_rm2sub_bridge/shared_range_certificate.py --screen independent_p_anchors_screen.json --lower 129 --upper 1024 --tag oa29_shared_q129_q1024 --verify
python -B workstreams/bch_rm2sub_bridge/shared_range_certificate.py --screen independent_p_gap339_screen.json --lower 339 --upper 339 --tag oa29_gap339 --verify
python -B workstreams/bch_rm2sub_bridge/verify_coverage_1024.py
```

## Dense-range diagnostics and next steps

Plain binary64 region coefficients underflow on denser probes.
`scaled_adaptive.py` instead stores a separate binary exponent for each
degree's matrix. It aligns adjacent exponents before each positive update.
Its directed version passes exact tests spanning exponents 0 through -3000.
`screen_dense_scaled.py` uses this representation with Arb region coefficients.
No dense screen is counted as an outward certificate.

The current Q=2048 adaptive screens remain vacuous. Several pure-group
bounds are already vacuous there, including groups 50--54 and 102--128.
Thus removing only the adaptive maximization loss cannot repair those
particular witnesses. This does not prove that the construction fails.

One further valid envelope is implemented in `syndrome_activation.py`.
For a uniform weight-j input X and any syndrome b, Fourier inversion gives

\[
\Pr[BX=b]=\frac1{2^s\binom tj}
 \sum_{u\in\mathbb F_2^s}(-1)^{u\cdot b}K_j(\operatorname{wt}(Au)).
\]

Here K_j is the length-t Krawtchouk polynomial. Let a_w count nonzero
A-codewords of weight w, M = 2^s-1, and beta_j = Pr[BX=0]. Define

\[
\rho_j:=\min\left\{M(1-\beta_j),
 \frac{M}{2^s}\left(1+\sum_w a_w
          \frac{|K_j(w)|}{\binom tj}\right)\right\}.
\]

The nonzero syndrome measure is pointwise bounded by rho_j times the
uniform nonzero law. Therefore an entering-zero row may use
(beta_j z^j, 0, rho_j z^j) in place of the arbitrary-state activation row.
This is domination of measures, not a claim that activation is uniform.
The implementation retains the original activation row for j=0,1,2 and
uses this alternative for j>=3. The existing nonzero-state rows are unchanged.

Exact tests check every toy syndrome probability and repeat all 8,736
weighted-prefix comparisons. The current pure-group diagnostic improves
at Q2048 but does not close it. This variant is not used in the verified
Q1--Q1024 ledger.

Two useful next investigations are stronger outer constraints for the
remaining troublesome shells, and retaining more input-weight information
in the moment bound. The latter may matter because the Bernoulli majorant
admits input weights absent from a fixed-weight outer group. Whether that
relaxation explains the dense gap remains to be established.
