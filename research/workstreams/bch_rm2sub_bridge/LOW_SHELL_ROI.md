# Next proof target for BCH-256

Follow-up: `PDUAL_RANK_PROGRESS.md` closes case 15 and proves a distance-32
kernel subcode. Case 7 remains open; the conditional margin below is not yet
an unconditional certificate. This note records the preceding payoff study.

The best measured next target is to exclude weight-30 words from P dual.
If that statement holds, a new exact LP witness raises the selected K=2^28
distance margin from 42.510028 to 43.585033 bits. The gain is 1.075004 bits.
The additional code statement remains unproved; this is a checked implication,
not a new unconditional distance certificate.

The investigation used five serial, bounded exact-LP runs. All five produced
rational optimal solutions that passed independent primal and dual checks.
No new spectrum enumeration or rank-witness search was run.

## Fixed endpoint and counted quantity

Keep the constituent and randomized encoder of `K28_DISTANCE_CERTIFICATE.md`:
BCH [256,128], the selected t64/s20 RM2Sub map, independent permutations, and
fresh independent uniform nonzero field multipliers. The probability is over
this setup randomness; the outer code C is fixed. The event is failure to
have relative distance greater than 0.1, not a complete SPIN security game.

Let P and Q be the length-256 extensions of the designed-distance-37 and
designed-distance-39 BCH codes, of dimensions 131 and 123. Write q_w=A_w(Q)
and h_w for the common spectrum of a nonzero Q-coset in P. Existing audits
give A_w(C)=q_w+31h_w and A_w(P)=q_w+255h_w.

The saved outward coefficients c_w include the number of possible nonzero
row positions. They bound the occupation-one score by sum_w A_w(C)c_w.
All higher occupancies retain the existing certified union U_higher.

This investigation fixes these coefficients. It does not measure the gains
from changing the inner map, improving transfer bounds, or changing K.
Its probability bounds reuse the existing certificates and their premises.

## Why another solve of the old model is insufficient

`low_shell_roi.py` replays all 1163 normalized rows of the old 196-variable
model and its exact weighted primal/dual certificate. The feasible spectrum
has the following low-shell counts and current weighted contributions:

| Weight | log2 of relaxed count | Share of the paired 38/40/42 score |
|---|---:|---:|
| 38 | 44.499195 | 99.902461% |
| 40 | 42.624628 | 0.097539% |
| 42 | zero count | 0% |

This is a fractional feasible spectrum, not an enumerator of a constructed
code. Evaluating the full current objective at this feasible point lower-bounds
the maximum over the unchanged LP. Comparing that value with the current
upper bound limits fixed-coefficient reoptimization to **0.002979 bits**.
Using only its first three shells gives the slightly weaker 0.002995-bit limit.

Consequently, materially better spectrum bounds require additional facts
that exclude this witness. Excluding one witness is necessary for such an
improvement, but does not establish that the improvement is substantial.

## Payoff from hypothetical low-shell caps

The existing integer caps for A38, A40, and A42 are respectively
24,866,368,872,377; 79,419,748,142,017; and 599,029,876,560,817.
The following table tightens only A38. It preserves every other shell term
and the higher-occupancy union.

| Hypothetical A38 cap | Reduction factor | Conditional full margin | Gain over frozen certificate |
|---:|---:|---:|---:|
| 12,433,184,436,188 | 2 | 43.478876 | 0.968847 |
| 6,216,592,218,094 | 4 | 44.444238 | 1.934210 |
| 3,108,296,109,047 | 8 | 45.377362 | 2.867334 |
| 1,554,148,054,523 | 16 | 46.252248 | 3.742220 |
| 388,537,013,630 | 64 | 47.668817 | 5.158788 |
| 97,134,253,407 | 256 | 48.448176 | 5.938147 |

These are sufficient conditional bounds, not predictions of the actual A38.
The complete receipt also tightens 38/40 together and 38/40/42 together.

For each pair w,256-w, put b_w=c_w+c_(256-w). The old exact dual supplies
sum_(w=38,40,42) a_w A_w <= B, with positive rational coefficients a_w.
For proposed caps T_w, maximize sum_w b_w x_w subject to
0<=x_w<=T_w and sum_w a_w x_w<=B. Filling variables in decreasing order
of b_w/a_w solves this three-variable relaxation. The implementation checks
an exact matching dual value, not just a numerical optimizer result.
Adding the unchanged other-shell bound and U_higher gives the table.

The zero-reduction row also exploits the existing individual caps and gives
42.511844 bits, a 0.001815-bit algebraic improvement. We have not replaced
the frozen aggregate certificate with that minor variant.

Tightening A38 alone eventually approaches roughly 48.84 bits with the other
caps fixed. Removing the contributions of all three low shells would still
leave only 52.815100 bits in this aggregation. The unchanged rest term is
about 0.07904% of the current Q1 bound. Thus the roughly 20-bit gap to the
reference-spectrum model cannot be recovered through A38 alone.

## Candidate constraints: exact payoff screen

The hull of Q is HQ=Q intersect Q dual. Its dual is a different, larger code
from P dual or Q dual. The hull-shell probes below constrain its enumerator.
Each row adds only the stated premise to the old LP. The proof of that premise
is separate from the exact LP implication.

| Added premise | Conditional full margin | Gain | Fixed-coefficient gain ceiling in this LP |
|---|---:|---:|---:|
| A20(HQ dual) <= 832,972,800 | 42.513455 | 0.003426 | 0.004554 |
| A20(HQ dual) = 0 | 42.524741 | 0.014713 | 0.015849 |
| A20=A22=A24=0 in HQ dual | 42.700978 | 0.190950 | 0.192238 |
| A30(P dual) = 0 | 43.585033 | 1.075004 | 1.077399 |
| A30(Q dual) = 0 | 43.586014 | 1.075986 | 1.078382 |

Ceilings come from the new rational feasible witnesses evaluated under the
full fixed objective. They limit these particular LP refinements, not the
value of different structural information that an enumeration might uncover.
The exact LP objectives use upward-rounded normalized coefficients; the
rounding and all primal/dual inequalities are checked over rational numbers.

The partial weight-20 cap comes from a retained old orbit receipt:
133,171,200 of the known 966,144,000 supercode words were excluded from HQ dual.
This turn inspected that receipt but did not rerun its large orbit audit.
Accordingly, its new LP result remains labeled as an implication as well.
The empty-shell premises are optimistic tests, not new empirical findings.

Already A30(P dual)=0 yields almost all the gain of the stronger Q statement.
Both codes have certified dual distance at least 30 and even dual words.
Thus this P statement would raise its dual-distance lower bound to 32.
In the LP it adds the degree-30 equation K(q+255h)=0. The Q statement gives
the two degree-30 equations Kq=Kh=0, since P dual is contained in Q dual.
Here K denotes the length-256 Krawtchouk transform at the specified degree.

The old feasible witness has positive dual counts at weight 30, but zero
counts at 32, 34, and 36, for both P and Q duals. Adding only zeros at those
higher weights leaves this witness feasible and cannot produce a large gain.

## Ranked next work and stopping criteria

1. **Try the P-dual weight-30 exclusion.** The retained rank search has 16
   Fourier cases. An independent replay this turn verified every saved path
   or root-run witness. Fourteen cases already reach rank at least 31.
   The cases with pivots 7 and 15 reach ranks 30 and 28, respectively.
   These two cases need stronger witnesses or exhaustive case refinements.
   A complete rank-31 cover, together with evenness, would prove the target.
   The old search stopped below target; the two remaining cases may be hard,
   and weight-30 words have not been ruled out. Failure to find a stronger
   rank witness would not establish that such a word exists.
2. **If that route stalls, pursue a weighted low-shell inequality or a
   genuinely coupled split constraint.** Use the existing exact bad spectrum
   as the first separation test. Target the score, not the full enumerator.
   The sensitivity table supplies concrete A38 thresholds for comparison.
3. **Deprioritize hull-shell enumeration for margin alone.** Even the strong
   three-empty-shell premise above buys less than 0.2 bits with this LP.
   The existing partial weight-20 result is inexpensive to reuse but tiny.

For the next bounded proof attempt, reuse the 14 completed P cases and search
only refinements of pivots 7 and 15. Independently check every accepted leaf
and the exhaustiveness of the case partition. Stop after the chosen search
budget if no new checked witness appears; report the residual cases rather
than treating the search target as an achieved bound. If a weight-30 word is
found, verify membership and retire the zero-shell premise.

The earlier independent OA6 block correction had exact zero-cost witnesses
for a different, weaker LP. That negative result does not settle a fully
coupled split model at the current strength. Similarly, older intersection
obstructions predate the final OA29 model. Do not transfer those negative
results to stronger models without a new feasibility check.

This investigation establishes no statistical confidence interval for the
reference-spectrum curve. Its evidence concerns relaxed worst-case spectra.
Keep any subjective curve estimate separate from the certified endpoint.

## Artifacts and replay

New sources are `low_shell_roi.py`, `low_shell_constraint_probe.py`, and
`test_low_shell_roi.py`. The compact result is
`generated/low_shell_roi_v1.json`; the five exact solutions and audits are
under `generated/low_shell_constraints_v1/`. The LP inputs are regenerated
from the frozen model, so the large exported LP files need not be committed.
No frozen source or previous certificate was edited.

From this directory:

```text
python -B low_shell_roi.py --output generated/low_shell_roi_v1.json --verify
python -B low_shell_constraint_probe.py pdual30_zero --verify
python -B low_shell_constraint_probe.py qdual30_zero --verify
python -B low_shell_constraint_probe.py hqdual20_partial --verify
python -B low_shell_constraint_probe.py hqdual20_zero --verify
python -B low_shell_constraint_probe.py hqdual20_22_24_zero --verify
python -B -m unittest test_low_shell_roi test_curve_spectrum
```

The original saved basis and the local QSopt_ex runtime were restored from
the retained ba80 workspace for discovery. No system package was installed.
The P-case witnesses remain in that read-only legacy bundle under
`generated/shift_rank_p32_probe/`; that folder's name is a target, not a proof
of distance 32. Its complete saved cover alone proves only distance 28.
The separately completed Q-dual proof still supplies the current distance-30
lower bound for both duals.
