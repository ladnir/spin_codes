# Small-message exact-spectrum replay: current status

## 1. Question and fixed construction

The experiment asks whether authenticated constituent spectra that were too
weak or inconvenient at \(k=2^{20}\) yield a concrete 10.9% distance
certificate at smaller message lengths and smaller inner memory.

The comparison contains only three outer models: exact-spectrum BCH,
exact-spectrum RM, and one random full-rank rate-half constituent sampled once
and reused across rows. Expander--accumulate and sparse-EA constructions are
outside this experiment. A separate global random rate-half calculation is
retained only as an ideal benchmark.

For a rate-half constituent \(C_0:[K]\to[B]\), the outer code is the direct
sum of \(L=k/K\) copies of one fixed deterministic \(C_0\). The same
constituent is reused in every row. Each row receives an independent uniform
coordinate permutation. The transposed layout then receives an independent
uniform position permutation in each region. RandomStepConv uses independent
random state/output maps at every position; one sampled set of maps is shared
by all messages. The output length is \(N=2k\), and the bad event is the
existence of a nonzero message whose output weight is below

\[
  D=\lceil 0.109N\rceil.
\]

The tested constituents are extended BCH \([32,16,8]\), extended BCH
\([128,64,22]\), and RM(4,9) \([512,256,32]\). Their exact weight spectra,
not random-code substitutes, feed the calculations.

## 2. What closed diagnostically

All results in this folder use nearest binary64 arithmetic. They are candidate
selection results, not outward-rounded certificates.

The extended BCH \([32,16,8]\) constituent never clears 40 bits at occupation
one in the tested range. The extended BCH \([128,64,22]\) constituent also
never clears 40 bits; its best tested occupation-one margin is about 34.30
bits near \(k=2^{12}\), \(M=22\).

RM(4,9) has several occupation-one points above 40 bits. The meaningful point
for an all-occupation random-like proof is

\[
  k=2^{13},\qquad N=2^{14},\qquad L=32,\qquad M=22.
\]

At this point:

- the complete exact-spectrum occupation-one sum has 40.71735 bits of
  diagnostic margin;
- the exact ordinary-spectrum occupation-two sum has 84.27710 bits of
  diagnostic margin; and
- the global random rate-half first-moment benchmark has 59.9674 bits of
  margin.

The occupation-one sum is dominated by the RM weight-32 shell. The
occupation-two sum is dominated by the pair of weight-32 shells. Because the
constituent is deterministic and row permutations are independent, the
ordinary weight spectrum suffices for this occupation-two calculation; no
genus-two enumerator is assumed.

The global random benchmark explains why the apparently stronger low-occupation
points at smaller \(k\) are misleading. Its margins at \(k=2^{11}\) and
\(k=2^{12}\) are only 19.5915 and 34.2274 bits. This is not a universal
impossibility theorem for structured codes. It shows that those lengths are
incompatible with the intended near-random first-moment route to 40 bits.

## 3. What did not close

No bound in this folder covers all occupations. In particular, the work does
not prove a 40-bit failure bound or a minimum-distance theorem.

Three high-occupation reductions were tested.

1. Charging every nonzero constituent word at minimum distance 32 discards
   nearly the entire exact spectrum and is vacuous.
2. Partitioning the exact spectrum into weight bands, deleting within each
   band, and applying a Bernoulli conditioning majorant is also vacuous. The
   first receipt incorrectly charged the all-zero candidate-input atom with
   the Chernoff low-weight factor. A corrected run splits that atom exactly,
   but live compositions entirely in the lowest RM weight band still dominate:
   occupation three alone has a positive 734.161-bit log bound.
3. Keeping only the earliest retained pulse is too destructive. A one-pulse
   RandomStepConv input retains an annihilation probability on the order of
   \(2^{-M}\), which cannot pay for the message multiplicity.

All three failures identify inadequate proof reductions. None establishes a
construction failure.

## 4. Reused-random constituent comparison

For the random outer, the constituent is sampled uniformly from full-rank
binary \([B,B/2]\) maps. Rejection sampling can realize this distribution, so
there is no rank-failure event. The one sampled map is reused in every row.

At \(k=2^{13}\) and \(M=12\), the matched Q1 comparison is:

| Block size | Structured constituent | Structured margin | Random margin |
|---:|:---|---:|---:|
| 32 | extended BCH \([32,16,8]\) | -10.134 bits | -10.529 bits |
| 128 | extended BCH \([128,64,22]\) | 7.181 bits | 6.180 bits |
| 512 | RM(4,9) \([512,256,32]\) | 17.048 bits | 91.019 bits |

The 512-bit random constituent is therefore not evidence that a smaller BCH
block works. It is only a matched control showing how far RM(4,9) lies from a
random spectrum at the same block size. At the two BCH block sizes, neither
random ensemble closes Q1, and BCH is slightly better because its guaranteed
minimum distance removes the random ensemble's low-weight tail.

For the matched RM control at \(B=512\), \(k=2^{13}\), and \(M=12\), the exact Q2 law splits into equal
and distinct local messages. Equal messages share one uniform nonzero output
vector. Distinct binary messages are linearly independent and map to a
uniform ordered pair of distinct nonzero vectors. The combined Q2 diagnostic
has 180.187 bits of margin; the equal-message sector dominates.

This does not close the random-outer theorem. For \(Q\ge3\), active local
messages can satisfy nontrivial binary linear relations. The proof must count
their rank/relation types and average the corresponding uniformly embedded
output tuples. Resampling a constituent independently by row would avoid this
obligation but would analyze the wrong construction.

## 5. Remaining proof obligations

A full certificate at the selected point requires all of the following.

1. For deterministic RM, replace the failed banded high-occupation majorant
   with a bound that retains more than a coarse lower endpoint for the
   weight-32 through weight-95 band.
2. For the reused-random outer, enumerate or dominate active local-message
   tuples by rank and linear-relation type for \(Q\ge3\).
3. If either diagnostic closes, replace every binary64 calculation by outward
   interval arithmetic and record auditable witnesses for each tilt.
4. Prove every retained deletion, conditioning, and regional transfer
   inequality in the stated probability space.
5. Combine deterministic-RM occupations without spending more than the approximately
   0.717-bit slack above the requested 40-bit margin. Memory 24 raises the
   occupation-one diagnostic only to about 40.831 bits, so memory alone does
   not create comfortable certification slack.
6. Audit the final receipt independently. Only then may the result be called
   a concrete distance certificate.

The random rank-one sector for \(3\le Q\le32\) is useful only as a diagnostic
of the same-size RM gap. It does not advance the smaller-block BCH cases. The
next construction-relevant decision is whether to accept a 512-bit outer
constituent or continue searching for a stronger authenticated spectrum near
block size 128.

## 6. Supporting material

- `exact_spectra_q1_small_k_memory_phase.json`: coarse BCH/RM occupation-one
  sweep.
- `rm49_q1_boundary_phase.json`: refined RM occupation-one boundary.
- `random_rate_half_d109_margin.json`: random rate-half reference margins.
- `rm49_k13_m22_q2_diagnostic.json`: exact occupation-two calculation at the
  selected point.
- `rm49_k13_m22_q3_32_banded4_p05_diagnostic.json`: pre-split banded
  high-occupation failure.
- `rm49_k13_m22_q3_32_banded4_p05_split_diagnostic.json`: corrected banded
  high-occupation failure after the zero-atom split.
- `rm49_k13_m22_q3_32_earliest_pivot_diagnostic.json`: earliest-pivot failure.
- `reused_random_constituent_q1_phase.json`: one-map-reused random-constituent
  Q1 sweep.
- `random512_k13_m12_q2_diagnostic.json`: reuse-aware random-constituent Q2
  result.
- `evaluate_exact_spectra_q1_phase.py`,
  `evaluate_fixed_spectrum_q2.py`,
  `evaluate_random_rate_half_margin.py`,
  `evaluate_reused_random_constituent_q1.py`,
  `evaluate_reused_random_constituent_q2.py`,
  `evaluate_banded_spectrum_high_q.py`, and
  `evaluate_earliest_pivot_high_q.py`: reproducible diagnostic programs.
