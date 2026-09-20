# Historical padded EBCH128 with RandomStepConv at finite \(k=2^{20}\)

## Status

This file records the earlier \(L=16560\) padded candidate. The current
power-of-two candidate has \(L=16384\) and no zero rows.

The fixed repeated \([128,64,22]\) extended BCH outer has a complete
11% certificate for RandomStepConv-M30. The outward verifier proves

\[
 \Pr[d_{\min}<233165]<2^{-26.1921836158}<2^{-20}.
\]

The earlier failed relaxations remain useful diagnostics. They discarded
essential cross-region input density; they were not counterexamples to the
construction. The dispersed parity-pivot transfer closes the missing dense
range.

## Candidate and probability space

Set

\[
 B=128,\qquad K=64,\qquad L=16560,\qquad L_0=16384,
\]

so that

\[
 k=KL_0=2^{20},\qquad N=BL=2{,}119{,}680,
 \qquad D=\lceil0.11N\rceil=233{,}165.
\]

Fix the concrete extended BCH constituent defined in
`FINITE_K20_EBCH128_OUTER_DEFINITION.md` and repeat it in all \(L\) outer
rows. The last \(L-L_0=176\) information rows are zero. Sample an independent
uniform coordinate permutation for every outer row and an independent uniform
row permutation in every transposed coordinate region. Sample one independent
RandomStepConv-M30 linear map at every inner position. All sampled objects are
chosen once and shared by every message.

The desired theorem bounds, over this setup randomness,

\[
 \Pr[\exists m\ne0:\operatorname{wt}(C_\omega(m))<D].
\]

The outer constituent is deterministic. It requires neither a spectrum test
nor enumeration of its \(2^{64}\) words during setup.

The distance verifier consumes the imported complete spectrum in
`scripts/EBCH128_64.wd`. It does not locally derive that spectrum from the
concrete generator polynomial. The spectrum-generator identification is an
explicit premise of the fully concrete theorem, as recorded in the exact
outer definition.

## Exact occupation-one result

For one active outer row, a constituent word of weight \(w\) activates
exactly \(w\) coordinate regions. Each active region contains one uniformly
located input difference. Summing the exact published BCH spectrum gives

\[
 \log_2 \mathbb E[Z_{D,1}]\le -26.2848527055
\]

in nearest-binary64 arithmetic. The weight-22 shell is dominant. This result
is a diagnostic, not an outward-rounded certificate.

## Rejected all-occupation relaxations

The first relaxation pointwise-majorizes every ordinary BCH word by a
uniform even row. It represents that row by 127 fair coordinates and one
parity coordinate, then deletes the parity-region inputs by monotonicity. For
occupations above 64, it also replaces an exact fixed-occupation region by a
Bernoulli-conditioned transfer. The resulting aggregate bound fails:

\[
 \log_2 \mathbb E[Z_D]\le 11624.1589
\]

at its dominant occupation \(Q=3787\). Raising the memory from 30 to 45 does
not materially change this value. The dominant loss at intermediate
occupation is the Bernoulli coefficient bound: it overweights the event that
an entire coordinate region contains no input difference.

At \(Q=L_0\), deleting the aligned parity inputs creates a 16,560-position
zero-input interval. A rare RandomStepConv state collision in that interval
can suppress most of the interval's output. The available finite-distance
budget cannot pay for that event with a 30-bit state.

A second relaxation retains only one input difference in every covered
coordinate region and uses the exact BCH spectrum to bound the number of
covered regions. It is algebraically positive but much weaker: after a state
collision, one retained input per long region permits artificial silent
gaps. Its failure therefore does not diagnose the actual dense input.

## Exact remaining obligation

A complete proof must preserve both of the following facts in one positive
transfer.

1. For each active outer row, the 128 coordinate bits obey the exact BCH
   shell law, or the uniform-even majorant with its parity dependence.
2. At intermediate and dense occupations, every coordinate region contains
   many shuffled candidate positions. The transfer may not replace the
   fixed count by a poorly tuned Bernoulli law or delete all but one input.

The closing proof chooses the parity pivot independently for each
uniform-even reference row before deleting it. The deleted positions are
therefore dispersed among the 128 regions. This representation does not
change the reference distribution. A positive two-fugacity bound retains the
resulting multinomial pivot loads. The outward verifier uses the exact BCH
spectrum at occupation one, an exact sparse recurrence through occupation 99,
and the dispersed-pivot bound thereafter.

## Evidence and reproducibility

- Exact spectrum: `scripts/EBCH128_64.wd`.
- Exact constituent, repetition, padding, and routing definition:
  `FINITE_K20_EBCH128_OUTER_DEFINITION.md`.
- Exact occupation-one evaluator:
  `evaluate_ebch128_randomstepconv_q1_exact.py`.
- Exact occupation-one receipt:
  `ebch128_randomstepconv_g1_s30_q1_exact_d11_diagnostic.json`.
- Even-majorant all-occupation evaluator:
  `evaluate_ebch128_randomstepconv_g1.py`.
- Its M30 and M45 receipts:
  `ebch128_randomstepconv_g1_s30_allq_d11_diagnostic.json` and
  `ebch128_randomstepconv_g1_s45_allq_d11_diagnostic.json`.
- Coverage-relaxation evaluator and receipt:
  `evaluate_ebch128_randomstepconv_coverage.py` and
  `ebch128_randomstepconv_g1_s30_coverage_d11_diagnostic.json`.
- Closing theorem:
  `FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_CERTIFICATE.md`.
- Closing evaluator and verifier:
  `evaluate_ebch128_randomstepconv_parity_pivot.py` and
  `certify_ebch128_randomstepconv_parity_pivot_outward.py`.
- Outward receipt and manifest:
  `ebch128_randomstepconv_g1_s30_parity_pivot_outward_d11.json` and
  `FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_MANIFEST.json`.

## Historical recommendation

The mathematical RandomStepConv comparison is closed. The next proof task is
to transfer the successful dispersed-pivot mechanism to an implementable
structured inner, preferably RM2Sub or its smallest necessary modification.
For the current unpadded candidate, first close the high-occupation gap in
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_STATUS.md`.
