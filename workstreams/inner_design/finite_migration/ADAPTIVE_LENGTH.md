# One-round IMT: length dependence and smaller steps

## Outcome and scope

Keep one transvection per update. Smaller steps improve short-length Q1
bounds, but do not completely remove the downturn. The selected BCH-256
certificates, implementations, and timings are unchanged. None of the new
smaller-step alternatives has a full certificate or a new timing here.

The paper's Figure 1 now includes BCH-256 as its main panel, with BCH-64
and BCH-128 for comparison. It replaces the separate mixing-round figure.
The earlier `MIXING_ROUNDS.md` study remains supporting evidence, not another
parameter in the paper's proposed construction.

## Fixed inputs and reproducible comparison

All cells have rate 1/2, N=2K, cutoff floor(N/10), zero initial state,
output before feedback, persistent state across regions, and no final flush.
Setup samples the existing independent row/region permutations and one
independent transvection per epoch. Maps are fixed and shared across messages.

The 110 cells use these map families:

| Outer | Tested (t,s) | Message exponents |
|---|---|---|
| BCH-64, BCH-128 | (8,6), (16,10), (32,15), (64,20) | 12,14,...,26 |
| BCH-256 | (16,10), (32,15), (64,19), (128,19) | 14,16,...,24 |

There are 88 one-round cells and 22 full-refresh controls for the fixed-step
baselines. The BCH-256 (128,19) baseline uses exactly the implemented
expansion and weight-five feedback maps. Smaller-outer t64 maps are exactly
the prior no-constant diagnostic maps. Other cells use the explicit fixed
recipes in `adaptive_step.py`; changing t changes the maps, not just a scalar
in an otherwise unchanged encoder. The receipts store every ordered column.
Diagnostic feedback is not claimed to have the implemented sparse-map cost.

For t>=64, the earlier no-constant map recipe is reused without editing it.
For smaller powers t=2^m, the nonconstant degree-at-most-two space has
dimension m+binom(m,2). We use its available dimension, with a fixed constant
shift excluding zero columns and the all-one expansion word. The t16 and
t32 maps therefore cannot retain s=20 within this family. Feedback has
distinct nonzero columns and full rank, constructed independently.

Full refresh means replacing the transvection by a uniform invertible map
per epoch. It keeps the corresponding baseline's A, C, t, and s unchanged.
It is not a universal outer-only ceiling or an implemented timing candidate.

The figure chooses the largest tested t within 0.1 bits of the best tested
Q1 value. This figure-level policy deliberately avoids chasing negligible
differences; the producer stores every cell rather than just this selection.
It is not a search over all maps or a measured performance optimum.

| log2 K | BCH-256 fixed t128,s19 | Adaptive (t,s) | Adaptive Q1 | Fixed-map full refresh |
|---:|---:|---|---:|---:|
| 14 | 31.272991 | (16,10) | 35.099048 | 55.204760 |
| 16 | 44.361697 | (32,15) | 49.037743 | 54.119380 |
| 18 | 50.287407 | (32,15) | 50.831070 | 52.347453 |
| 20 | 50.030502 | (128,19) | 50.030502 | 50.413639 |
| 22 | 48.373392 | (128,19) | 48.373392 | 48.454145 |
| 24 | 46.453333 | (128,19) | 46.453333 | 46.472097 |

For both smaller outers the adaptive policy selects t16 at exponent 12,
t32 at 14 and 16, and t64 at 18 and above. At BCH-128/K12, Q1 improves
from 19.179149 to 24.274701 bits. No tested t8 point is selected.

A separate state-control screen holds s=15 while changing the step:

| Outer / log2 K | t32,s15 | t64,s15 | t64, larger state |
|---|---:|---:|---:|
| BCH-128 / 12 | 21.639320 | 19.259210 | 19.179149 (s20) |
| BCH-128 / 16 | 34.612094 | 33.408715 | 33.650270 (s20) |
| BCH-256 / 16 | 49.037743 | 45.245498 | 45.525780 (s19) |
| BCH-256 / 18 | 50.831070 | 50.090203 | 50.720323 (s19) |

At BCH-256/K16, reducing state from s19 to s15 at t64 loses 0.28 bits;
reducing the step from t64 to t32 at s15 gains 3.79 bits in these map
families. Thus the observed improvement is not merely an advantage of a
smaller state. The t64 state maps are nested; changing t still changes
the concrete maps. This 12-cell binary64 control screen is supplementary,
not part of the 110-cell replay or a proof of map-independent monotonicity.
Reproduce it with `adaptive_controls.py --output <fresh-controls.json>`.

The fixed-map reference and changed-map adaptive curves answer different
questions. Their vertical difference is not a controlled measurement of
mixing loss for the adaptive map. The separate exploratory screen includes
matched refresh controls for t16 and t32 and demonstrates their smaller-state
penalty. For example, BCH-128/K12 refresh gives 34.39 bits for t16,s10,
versus 38.83 for the fixed t64,s20 reference in the final study.

## Root cause 1: length changes more than the row multiplicity

Let L=2K/B and h=L/t. A single outer word of weight w produces one singleton
input in each of w active regions. The active region positions are uniformly
permuted, and the singleton's position within each active region is uniform.
There are h epochs in each region, not K/t independent opportunities between
every two inputs.
Increasing B at fixed K,t shortens each region and reduces h. Thus a larger
outer is not simply a vertical margin improvement with unchanged inner
geometry. At K16, the BCH-128 fixed curve has h=16 (t64), while the BCH-256
implemented baseline has h=4 (t128). Their dominant low outer weights also
differ; the adjacent-pair expectation is 22*21/128 versus 38*37/256.

The Q1 expression has the form L sum_w A_w Q_w(K,t,s). The factor L explains
the approximately one-bit-per-doubling loss once the conditional bounds Q_w
are nearly invariant with K. At short lengths Q_w itself changes substantially.
The K-axis is therefore not merely the union-bound multiplicity axis.

## Root cause 2: a specific short cancellation path

This calculation is exact and unweighted; it is not the entire low-output
event. Suppose state is zero before the first of two active regions and
there are g empty regions between them. Each active region contains one
input, with independent uniform positions. Let a,b in {0,...,h-1} be their
epoch indices. The first input injects C e_p. There are

    d = g h + (h-1-a) + b

empty epochs before the next input, followed by one mixer immediately
before that second feedback is added. Thus d+1 mixers act on the injected
state. Their exact nonzero-state marginal is

    2^(-(d+1)) delta_(C e_p) + (1-2^(-(d+1))) Uniform(nonzero).

Because the feedback columns are distinct and nonzero, independent uniform
input coordinates repeat the same column with probability 1/t. Averaging
the retention coefficient over the epoch positions gives

    rho_(h,g) = 2 * 2^(-g h) * (1-2^(-h))^2 / h^2,
    p_cancel = rho_(h,g)/t + (1-rho_(h,g))/(2^s-1).

`test_adaptive_step.py` checks the geometric sum exactly and independently
enumerates all actual three-bit transvections and their products to verify
the cancellation probability. No independent per-message setup is assumed.

Consequences:

- Increasing s suppresses the fresh-state term, not the retained-state term.
- For adjacent regions and large h, the retained-state term is approximately
  2/(t h^2)=2t/L^2. Shortening t helps in this regime.
- A whole empty region adds a factor 2^(-h), so close active-region pairs
  are especially relevant. For a uniformly permuted weight-w support among
  B regions, the expected number of adjacent active pairs is w(w-1)/B.
- There is no global monotonicity theorem in t: the 1/t repeat-column factor
  increases when t shrinks, h changes, and s and the maps may change too.
- Returning to zero limits how long the state can emit weight. Conversely,
  a retained state can repeatedly emit a low-weight expansion word without
  returning to zero. Cancellation is not the only relevant tail mechanism.

For L=512,s=19, the local adjacent-region cancellation probabilities are
about 0.000860 for t128 (h=4) and 0.000486 for t64 (h=8). These are local
probabilities, not claimed contributions to the full union.
For comparison, at L=64,s=20 this probability increases from about 0.00781
to 0.00879 when t falls from 64 to 32. The repeated-column factor can win
at very small h. Even the local mechanism does not justify a blanket claim
that every reduction in t improves every part of the tail.

The relevant paths are not typical gaps. Set D=d+1 for adjacent regions.
Conditional on a retained-component cancellation, the epoch endpoints are
weighted by 2^(-D), rather than uniformly. An exact geometric sum gives

    E[D | retained-component cancellation] = 3 - 2h/(2^h-1).

This tends to **three epochs**, even though the unconditioned mean is h.
The first input tends to sit near the end of its region and the second near
the start of the next. Increasing h makes these paths rarer, not longer.
For repeated input coordinate p and retained state q=C e_p, their total
emitted weight is exactly

    D * wt(A q) + 2 - 2 (A q)_p.

Thus the same mechanism both cancels the state and limits output to a short
pulse of O(t) bits. Tests check the conditional mean with exact fractions and
the output-weight identity by direct state simulation. This is a local
single-message coupling of the exact state kernel, not an assertion that
the setup samples identity matrices with probability 1/2 for all messages.
Combining many such excursions with the outer support law remains a
separate task; these identities do not prove domination of the full tail.

## Root cause 3: a measurable loss in the old bound

All final Figure 1 cells use the same injection-shell envelope from the
mixing study, including the one-round baselines. It retains the actual
expansion-weight shell of the first feedback state, rather than immediately
treating it as arbitrary. The older two state-size figures still use their
authenticated earlier envelope; the manuscript distinguishes the studies.

`ADAPTIVE_LENGTH_DIAGNOSIS_v1.json` isolates proof refinement and transfer
terms without modifying the encoder:

| Outer / exponent / fixed inner | Old envelope | Injection-shell envelope | Full refresh | Lazy-return terms removed* |
|---|---:|---:|---:|---:|
| BCH-128 / 12 / (64,20) | 16.278446 | 19.179149 | 38.830766 | 24.024018 |
| BCH-128 / 16 / (64,20) | 33.422149 | 33.650270 | 36.616338 | 34.804960 |
| BCH-256 / 16 / (128,19) | 43.582125 | 44.361697 | 54.119380 | 49.952225 |
| BCH-256 / 20 / (128,19) | 50.015358 | 50.030502 | 50.413639 | 50.145043 |

*The last column deletes positive lazy-return terms from the envelope. Its
scores are NOT valid Q1 bounds or results for a real replacement encoder.
This is a diagnostic of which transfer terms enlarge the computed bound.
It shows that lazy returns contribute substantially at short lengths,
while their removal does not reproduce the full-refresh result. The gaps
are not an additive decomposition into independently attributable causes.

The current computation is not a measurement of true code distance. We have
identified an actual local cancellation mechanism and a substantial related
term in the proof bound, but have not proved that it dominates the true
distance-failure event. Both implementation choices and proof slack matter.

## Reproduction and evidence boundary

Run from the repository root, with the existing research dependencies:

```text
python -B workstreams/inner_design/finite_migration/adaptive_length_study.py run --output <fresh-grid.json>
python -B workstreams/inner_design/finite_migration/adaptive_length_study.py verify --source <fresh-grid.json> --output <fresh-replay.json>
python -B workstreams/inner_design/finite_migration/adaptive_length_study.py diagnose --output <fresh-diagnosis.json>
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p test_adaptive_step.py
python -B paper/build_imt_length_figure.py --check
python -B paper/build_imt_parameter_figures.py
```

The producer uses scaled positive coefficient products; the verifier
reconstructs every map and reevaluates all 110 cells using log-domain vector
products. Both use nearest binary64 arithmetic. Neither is an outward proof
replay. BCH-64/128 use exact spectra; BCH-256 uses the same deterministic
weighted-spectrum inequalities in every cell, never a guessed spectrum.
Source and spectrum hashes accompany the receipts. Numerical data stay local;
source, tests, documentation, and the generated TeX figure can be shared.

All 110 cells passed the log-domain replay. The largest margin/witness
disagreement was 2.28e-13 (the largest margin discrepancy in bits is bounded
by the same value). The paper generator pins producer SHA-256
`e5fbfaace653ef5597817c4241ff31789fa82d4b57a5b4426d824424aa747de0`
and replay SHA-256
`d5d5a536e8fc1ab86b6987c713281252199cc3744d73ec638c934553b16d1d17`.

Frozen study inputs, older certificates, and proof-producing source files are
not edited by this exploration. Small differences from other retained Q1
values arise from their different tilt banks and activation envelopes.

## Recommended next experiment

Prioritize BCH-256, t32,s15 at K16 and K18. Check sparse and dense occupancies
before any performance claim, then implement and benchmark serially against
the existing t128,s19 kernel. Also calibrate s at fixed t64 using explicitly
matched map families; nominal state dimension alone does not identify a map.
The remaining very-short-length loss may require better maps or a sharper
weighted-state transfer, not simply another reduction in t.
