# Unconfounding RM2Sub persistence and epoch length

## Question

The first 10% screen appeared to show that increasing the epoch length $t$
while decreasing the state dimension $s$ could improve the security margin at
fixed

\[
  p=s+\log_2t.
\]

That interpretation was not justified. The compared $A/B$ maps were sampled
independently, and one occupation-two screen split its tilt witnesses across
separate aggregate runs. This experiment removes both problems.

## Fixed construction

The outer constituent is one fixed RM$(4,9)$ code with parameters
$[512,256,32]$. It is repeated in every outer row. The message length is
$k=2^{16}$, so the output length is $N=2^{17}$ and every transposed region has
length $L=256$. The bad event is output weight at most

\[
  \lfloor0.10N\rfloor=13107.
\]

Routing uses independent uniform row-coordinate permutations and independent
uniform region permutations. One audited RM2Sub map is fixed for each tested
configuration. Independent nonzero epoch scalars supply the remaining inner
randomness.

All margins below use nearest binary64 arithmetic. Occupation one uses the
complete 121-point grid with spacing 0.1. Occupation two uses one complete
seven-point grid

\[
  \{-6,-5,-4,-3,-2,-1,0\}
\]

for the log-surprisal witness. Each outer-weight pair takes its own pointwise
minimum before the spectrum sum.

## Witness-aggregation correction

For each outer-weight pair $(w_1,w_2)$ and tilt $u$, let
$F(w_1,w_2;u)$ denote the logarithmic upper bound before the outer spectrum
sum. The valid finite-grid calculation is

\[
  \sum_{w_1,w_2}
  \exp\!\left(\min_{u\in\mathcal U}F(w_1,w_2;u)\right).
\]

The invalid calculation partitioned $\mathcal U$ into two subsets. It computed
one completed spectrum sum for each subset and retained the better aggregate.
In general,

\[
  \min_j\sum_{w_1,w_2}\exp\!\left(\min_{u\in\mathcal U_j}F\right)
  \;\ge\;
  \sum_{w_1,w_2}\exp\!\left(\min_{u\in\cup_j\mathcal U_j}F\right).
\]

The inequality can be large because different weight pairs require different
tilts. The split-grid result of -11.587 bits for $t=128,s=13$ was therefore a
valid but weak upper bound, not evidence that the parameterization failed.
The complete pointwise grid gives 85.947 bits.

## Nested state experiment

Fix one ordered family of self-orthogonal RM$(2,7)$ evaluation vectors. Define
$A_s$ from the first $s$ vectors and define $B_s$ as the transpose of the same
generator matrix. Then

\[
  A_s\subset A_{s+1}
  \quad\text{and}\quad
  B_sA_s=0.
\]

The chain uses seed 3390173185. It is one unselected deterministic sample.
No search or optimization occurs. Exact enumeration gives:

| $s$ | $d(A_s)$ | $d(\ker B_s)$ | weight-four kernel words | $Q=1$ margin |
|--:|--:|--:|--:|--:|
| 12 | 48 | 4 | 5,216 | 41.478 |
| 13 | 48 | 4 | 2,656 | 42.147 |
| 14 | 48 | 4 | 1,376 | 42.536 |
| 15 | 48 | 4 | 672 | 42.751 |
| 16 | 48 | 4 | 384 | 42.874 |

The complete occupation-two screen gives 85.819 bits at $s=13$ and 87.181
bits at $s=14$. The adjacent-state change is smooth. The data do not establish
a threshold at $s=14$.

This experiment isolates prefix extension within one chain. It does not prove
that every sampled chain behaves similarly.

## Matched-persistence epoch experiment

Fix persistence exponent $p=20$. The three configurations are

\[
  (t,s)\in\{(64,14),(128,13),(256,12)\}.
\]

At region length 256, these configurations have four, two, and one epoch per
region.

| $(t,s)$ | epochs per region | $Q=1$ margin | $Q=2$ screen |
|:--|--:|--:|--:|
| $(64,14)$ | 4 | 42.583 | 87.424 |
| $(128,13)$ | 2 | 42.147 | 85.947 |
| $(256,12)$ | 1 | 41.406 | 83.721 |

Both screened margins decrease as $t$ increases. The present evidence does
not support the claim that increasing $t$ while decreasing $s$ improves the
security margin at fixed persistence.

The experiment does not show that shorter epochs are globally optimal. Epoch
length also changes the exact $A$ and kernel spectra. Larger epochs require
fewer state updates and can improve encoder performance. The final design must
therefore compare complete distance bounds and measured XOR costs.

## Status and next obligations

The following items are exactly checked:

- nesting of the $A_s$ generator bases;
- full rank of every $A_s$;
- the identities $B_sA_s=0$;
- every recorded $A_s$ spectrum; and
- every kernel spectrum through an exact MacWilliams transform.

The numerical margins are diagnostics. They do not use outward rounding.

The next occupation experiment is complete. At $Q=3,\ldots,8$, the
$t=64,s=14$ map has the largest margin at every occupation, followed by
$t=128,s=13$ and then $t=256,s=12$. The $t=64,s=14$ sparse transfer remains
positive through $Q=29$, with 30.683 bits there, and becomes vacuous at
$Q=30$.
`RM2SUB_Q_LADDER.md` records the exact reduction and the binary64 receipts.

The next useful experiment is a finite dense-transfer bridge for
$30\le Q\le256$, with $Q=29$ retained as an overlap check. A multi-seed nested
replay is lower priority because the controlled epoch ordering is already
consistent through occupation eight.

## Reproduction files

- `generate_nested_rm2sub_family.py` constructs the nested family.
- `evaluate_rm2sub_nested_t128.py` evaluates its occupation screens.
- The generated `nested_rm2sub/` bundle fixes the chain and records the full
  Q1 and Q2 grids. Commit `97c7a7e` preserves that bundle; the Overleaf tip
  omits it to stay below the project file limit.
- `rm2sub_primary_frontier_rm49_t128_s13_complete_grid_d100.json` corrects the
  original selected-map result.
- `rm2sub_epoch_geometry_rm49_t64_s14_q1_d100.json` and
  `rm2sub_epoch_geometry_rm49_t64_s14_q2_d100.json` complete the matched-$p$
  epoch comparison.
- `RM2SUB_Q_LADDER.md` and the `rm2sub_q_ladder_*_d100.json` receipts extend
  the comparison beyond occupation two.
