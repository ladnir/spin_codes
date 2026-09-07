# What is weakening the BCH-256 endpoint bound?

The fixed t64_s20 inner remains the fully certified baseline at K=2^20.
The cheaper t64_s16 and t128_s19 maps have unresolved dense endpoints.
This investigation separates single-band behavior, band-mixture losses,
and numerical evaluation. It does not assume the unknown BCH spectrum.

All calculations below use L=8192 active outer rows, 256 outer coordinates,
and H=209715. They retain the selected maps, shell caps, and setup convention
from [CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md).
The new screens are diagnostics, not additional full certificates.

## A cheap endpoint calculation

Fix one of the existing 13 weight bands B. Consider only messages whose
nonzero outer rows all have weights in B. The existing coefficient comparison
uses an auxiliary Bernoulli support probability p and a row cost

\[
D_B(p)=\max_{w\in B}
\frac{c_w}{\binom{256}{w}p^w(1-p)^{256-w}},
\]

where c_w is the existing upper cap on the number of outer words of weight w.
The singleton band {256} uses p=1 and cost 1.
These probabilities describe the counting measure, not a change to the encoder.

Let E_j(lambda) be the frozen three-state upper transfer for an epoch with
j input ones and exponential weight exp(-lambda * output weight).
At full occupancy, the auxiliary inputs in a homogeneous band are iid.
The epoch mixture is therefore

\[
S_p(\lambda)=\sum_{j=0}^{t}\binom tj p^j(1-p)^{t-j}E_j(\lambda).
\]

The homogeneous-band bound has the form

\[
U_B=e^{\lambda H}D_B(p)^L
e_0^\mathsf{T}S_p(\lambda)^{256L/t}\mathbf1.
\]

This identity avoids degree-L region polynomials for this restricted calculation.
`endpoint_band_diagnostic.py` evaluates epoch entries with 256-bit Arb and then
uses binary64 log arithmetic for optimization and powering. It is not an
outward evaluation of U_B. Tests compare the epoch-mixture identity against
the full fixed-weight construction on a small exact-rational example.

## Single bands do not reproduce the weak mixed bound

The six-tilt screen uses log(lambda) in {-0.5, 0, 0.2, 0.4, 0.6, 0.8}.
It takes roughly eight seconds per selected map. At log(lambda)=0.6:

| Inner | Worst homogeneous-band diagnostic margin |
|---|---:|
| t64_s20 | 35652.398 bits |
| t64_s16 | 35648.355 bits |
| t128_s19 | 35652.695 bits |

In each case the worst band is {102,104,...,128}, with p close to one half.
Every homogeneous band has a positive diagnostic margin at that tilt.
The never-activated contribution is negligible within the worst homogeneous band.

The old adaptive bound also pays L*log2(13), approximately 30314 bits, for
band labels. Even subtracting that amount from the homogeneous diagnostic
leaves a positive number. However, this is only a comparison proxy:
**the maximum over homogeneous bands does not bound all mixed assignments.**
Multiplying that maximum by 13^L does not supply the missing inequality.

The pure-band optimizer selected log(lambda)=0.8 for a follow-up t64_s16 trial.
Its worst homogeneous-band diagnostic margin was 56713.942 bits.
The frozen outward adaptive evaluator nevertheless returned a weak endpoint
bound, with diagnostic margin -179199.857 bits. This trial took about 91 seconds
and was excluded from certificate coverage. No numerical replay was needed to
decide that its retained upper bound was unhelpful.

Optimizing homogeneous cases alone therefore does not optimize the current
mixed-band certificate. The current comparison allows band choices to vary
inside an entrywise-max recurrence; it is not a fixed homogeneous assignment.

## Checking the never-activated contribution

If the state never activates, the output equals the input. At full occupancy
that input has weight at least 38L=311296, above H. Such a path cannot be bad.
This observation alone does not bound paths that activate and subsequently return
to zero, nor does it justify subtracting one upper bound from another.

To diagnose this path class, write the kernel enumerator as
\(K(x)=\sum_j k_j x^j\). Its region coefficient at input weight j is

\[
\frac{[x^j]K(x)^{L/t}}{\binom Lj}e^{-\lambda j}.
\]

`endpoint_zero_log_diagnostic.py` computes the polynomial coefficients as exact
integers and evaluates the scalar adaptive recurrence in the log domain.
For the t64_s16 tilt-0.8 witness, its diagnostic margin is 71315.296 bits.
That is over 250000 bits below the contribution represented by the weak full
upper bound. Thus removing never-activated paths is not a promising repair for
this witness, even though their actual contribution to the bad event is zero.

The earlier scalar screens `endpoint_zero_diagnostic.py` and
`endpoint_zero_exact_diagnostic.py` used a scaled binary64 recurrence.
Their saved scalar results are superseded and must not be used for inference.
Exact integer polynomial coefficients alone did not fix the discrepancy.
On a single-line L=2048 test, the scaled recurrence returned natural logarithm
-449.265043 instead of -522.896501. Log-domain recurrence agreed with the
independently evaluated epoch formula. Alignment with a zero coefficient's
exponent forced a small nonzero neighbor into subnormal arithmetic.

This is evidence about the scalar diagnostic's numerical slack. It does not
by itself establish the same problem in the full certificate evaluator.
The existing outward certificates have not been edited or invalidated by
this investigation.

## Artifacts and next decision

Local outputs are under `generated/endpoint_bands_<map>_v1.json`,
`generated/endpoint_trial_t64_s16_v1/`, and
`generated/endpoint_zero_log_t64_s16_v1.json`. The scalar files without `log`
in their names are superseded as described above.
Only the source, tests, and interpretation belong in Git; outputs remain ignored.

`endpoint_log_replay_diagnostic.py` evaluated the full adaptive recurrence in
the log domain. It held the witness, 256-bit Arb region coefficients, and
outward coefficient pairs fixed. The run took approximately 86 seconds.
Its margin was -179199.856996462 bits, compared with -179199.856996465 from
the scaled outward producer: a difference of approximately 3.4e-9 bits.
Thus recurrence scaling does not explain the weak bound for this witness.
The comparison does not test slack already present in their shared Arb
coefficient enclosures.

The final envelope assigns approximately 0.9999999604 of its terminal mass
to the live-state class. These normalized entries describe the upper transfer,
not physical encoder-state probabilities. The retained coefficient envelope
uses bands beginning at 38, 102, 130, 156, and 256.
The comparison is saved in `generated/endpoint_log_t64_s16_v1.json`.

Thirty tests passed across the search and arithmetic regressions, including
the epoch-mixture identity, exact scalar kernel coefficients, long-length
log recurrence, and matrix recurrence. Existing receipt authentication tests
and repository hygiene checks also passed.

The next bounded investigation should remove selected lines from the adaptive
envelope as a diagnostic, especially the low-weight and all-one lines.
Such an ablation cannot certify the omitted band assignments. It can identify
which mixtures merit an exact split or fixed-type treatment. The other
workstream's joint fixed-type refiners are a candidate backend, but they need
our selected maps and BCH-256 caps. No hard mixture has yet been ruled out.

No new endpoint certificate or full optimized configuration was obtained.
The useful result is a sharper research target: control mixed weight types,
rather than optimize homogeneous bands or remove only never-activated paths.
