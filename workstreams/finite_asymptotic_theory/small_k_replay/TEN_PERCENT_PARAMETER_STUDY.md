# Ten-percent parameter study for finite Structured SPIN variants

## Objective

The primary target is relative distance $0.10$. A candidate should eventually
satisfy

\[
  \Pr[d_{\min}(\mathcal C)<0.10N]<2^{-\lambda},
\]

for a requested failure margin $\lambda$. The probability is over the sampled
setup objects defined below. The default target is $\lambda=40$ bits.

A distance near $0.109$ or $0.11$ is a bonus. It does not influence parameter
selection unless the same candidate already meets the 10% target at comparable
cost.

The study has two outputs for each admissible tuple $(k,0.10,\lambda)$:

1. the fastest fully certified construction; and
2. the fastest conditional construction under one stated outer-spectrum
   assumption.

The second output is called *extrapolated*. It is not called certified.

The current occupation ladder is exploratory. It records signed margins and
does not discard a parameterization merely because its margin is below 40
bits. The 40-bit value is a later selection target, not a gate on this
baseline scan.

## Construction under study

Fix one binary linear constituent

\[
  C:\mathbb F_2^K\longrightarrow\mathbb F_2^B.
\]

The outer encoder partitions a $k$-bit message into $L=k/K$ rows. It applies
the same map $C$ to every row. The output length is $N=BL$. The study does not
sample a new constituent for each row.

Routing samples one independent uniform coordinate permutation for each outer
row. It also samples one independent uniform position permutation in each
transposed region. Setup fixes all routing permutations.

An RM2Sub epoch has input $X\in\mathbb F_2^t$, state
$Q\in\mathbb F_2^s$, output $Y\in\mathbb F_2^t$, and update

\[
  Y=X+A(Q),\qquad Q'=\alpha Q+B(X).
\]

The maps $A:\mathbb F_2^s\to\mathbb F_2^t$ and
$B:\mathbb F_2^t\to\mathbb F_2^s$ are fixed and satisfy $BA=0$. Each epoch
samples an independent $\alpha\in\mathbb F_{2^s}^{\times}$. Setup fixes these
scalars. The current selected maps are subcodes of RM$(2,\log_2t)$ and have
exactly audited $A$ and $\ker B$ spectra.

RandomStepConv remains a calibration model. It is not the primary inner for
the parameter study.

## Native lengths and wrappers

The native construction requires

\[
  K\mid k\quad\text{and}\quad t\mid L=k/K.
\]

The current proof engine uses only complete epochs. It does not pad a region,
shorten an epoch, or mix epoch sizes. A result at a non-native length therefore
requires a separate wrapper theorem. No such wrapper is active in the present
receipts.

For rate-one-half constituents and power-of-two parameters, the rule becomes

\[
  2^{\log_2 K+\log_2 t}\mid k.
\]

For example, RM$(4,9)$ has $K=256$. It first admits $t=128$ or $256$ at
$k=2^{15}$ or $2^{16}$, respectively. The first tranche uses $k=2^{16}$ so
both epoch lengths are comparable without a wrapper.

## Fully certified output

A row may enter the fully certified frontier only after all of the following
items are complete.

1. The outer constituent is fixed and its required spectrum data are proved or
   exactly authenticated.
2. The proof covers every nonzero outer occupation $1\le Q\le L$.
3. Arithmetic is outward rounded, and the finite witness grid has a verified
   interpretation.
4. The sum over occupations is below $2^{-\lambda_{\rm dist}}$.
5. Every setup rejection or selection step has a proved failure bound
   $2^{-\lambda_{\rm setup}}$.
6. The total failure probability is bounded by

   \[
     2^{-\lambda_{\rm dist}}+2^{-\lambda_{\rm setup}}
     \le 2^{-\lambda}.
   \]

7. The implementation is proved equivalent to the analyzed construction.
8. Runtime is measured for the exact candidate implementation.

Binary64 occupation screens do not satisfy these conditions.

## Extrapolated output

The extrapolated track changes only the outer-spectrum input. It uses the same
routing law, RM2Sub transfer, occupation split, outward arithmetic, and final
union bound as the certified track.

For a sampled rate-one-half constituent $C$, let $A_w(C)$ denote the number of
messages whose codeword has weight $w$. Before sampling $C$, fix integer caps
$U_0,\ldots,U_B$. The narrow outer-spectrum assumption is

\[
  A_0(C)=1
  \quad\text{and}\quad
  A_w(C)\le U_w\quad(1\le w\le B).
  \tag{OuterCap}
\]

One sampled full-rank constituent is then fixed and repeated in every row.
Conditioned on `OuterCap`, the deterministic proof engine uses the caps
$U_w$. It does not replace powers or products of realized multiplicities by
powers or products of their expectations.

The assumption is narrow because it concerns one finite table at the outer
interface. A later concentration proof can convert the conditional result into
a setup theorem. Until that proof exists, the study reports `OuterCap` as an
assumption and does not assign it a failure probability.

The same interface can describe a sampled BCH-plus-fanout constituent. The
fanout sampling procedure and the caps must be fixed before its proof-engine
run. An expected spectrum alone does not establish `OuterCap`.

## Search procedure

The search uses an adaptive funnel.

1. Evaluate occupation one on a broad parameter grid.
2. Evaluate occupation two only for candidates that survive occupation one.
3. Evaluate representative sparse, intermediate, and dense occupations for
   the remaining candidates.
4. Run the full binary64 occupation range only for survivors.
5. Outward-round the finalist calculation.
6. Benchmark only candidates that remain on the proof-cost frontier.

The primary epoch lengths are $t=128$ and $t=256$. The $t=64$ results remain
a calibration reference. Initial persistence exponents are 20, 22, 24, and
26, where

\[
  p=s+\log_2 t.
\]

Persistence is only a first-order tuning coordinate. Exact $A$ and kernel
spectra remain part of every candidate description.

The initial message lengths are $2^{14}$ and $2^{16}$. The first length screens
the BCH $[128,64,22]$ family. The second also admits RM$(4,9)$ with both primary
epoch lengths. Later scans can extend the same rules to other native lengths.

## First-tranche result

The complete 121-witness occupation-one screen gives the following RM$(4,9)$
results at $k=2^{16}$ and distance 10%. The occupation-two values are coarse
elimination screens. They use the finite witness sets recorded in the JSON
receipts.

| RM2Sub | $d(A)$ | $d(\ker B)$ | kernel words of weight four | $Q=1$ margin | $Q=2$ screen |
|:--|--:|--:|--:|--:|--:|
| $t=128,s=13$ | 56 | 4 | 2,016 | 42.147 | 85.947 |
| $t=128,s=14$ | 48 | 4 | 928 | 42.536 | 87.181 |
| $t=128,s=15$ | 48 | 4 | 416 | 42.751 | 88.062 |
| $t=128,s=17$ | 48 | 4 | 32 | 42.940 | 88.902 |
| $t=128,s=19$ | 48 | 6 | 0 | 42.991 | 89.154 |
| $t=256,s=12$ | 120 | 4 | 81,600 | 41.406 | 83.721 |
| $t=256,s=14$ | 112 | 4 | 19,136 | 42.429 | 86.883 |
| $t=256,s=16$ | 112 | 4 | 4,032 | 42.746 | 88.290 |
| $t=256,s=18$ | 96 | 4 | 576 | 42.842 | 88.737 |

Every row in the table is diagnostic. None is a distance certificate.

The original screen incorrectly reported -11.587 bits for $t=128,s=13$.
That calculation split the tilt witnesses across two runs and compared the
two aggregate margins. The proof requires a pointwise minimum over tilts for
each outer-weight pair before the spectrum sum. A complete seven-point grid
gives 85.947 bits. The earlier elimination was therefore invalid.

A nested-family replay confirms that the adjacent-state behavior is smooth.
It fixes one ordered RM$(2,7)$ basis and defines $A_s$ by prefixes. The nested
$s=13$ and $s=14$ maps give 85.819 and 87.181 occupation-two bits. Thus the
data do not establish a state threshold at $s=14$.

At fixed $k=2^{16}$ and persistence exponent 20, the corrected comparison is:

| RM2Sub | epochs per region | $Q=1$ margin | $Q=2$ screen |
|:--|--:|--:|--:|
| $t=64,s=14$ | 4 | 42.583 | 87.424 |
| $t=128,s=13$ | 2 | 42.147 | 85.947 |
| $t=256,s=12$ | 1 | 41.406 | 83.721 |

The tested margins decrease as the epoch length increases. The first tranche
therefore gives no evidence that one epoch per region improves the proof
margin. Epoch length can still improve implementation cost by reducing the
number of state updates.

The weakest surviving point is $t=256,s=12$. Its occupation-one margin is only
1.406 bits above the 40-bit target. It is therefore an aggressive performance
candidate. The matched-persistence alternatives $t=64,s=14$ and
$t=128,s=13$ retain more proof margin. Their final ordering depends on the
all-occupation sum and measured XOR cost.

The BCH $[128,64,22]$ occupation-one margin stays below 40 bits throughout the
first grid. Its best tested value is 37.967 bits at $k=2^{14}$ and
$t=128,s=19$. It remains useful as calibration data but is not a 40-bit
finalist at these lengths.

## Occupation-ladder update

A positive two-colour recurrence now handles any mixture of regular RM rows
and the unique all-one RM word without unstable finite differences. A direct
enumeration audit passes on all 15 type pairs of a six-position test.

At matched persistence 20, a 0.1-spaced complete tilt grid gives:

| $Q$ | $t=64,s=14$ | $t=128,s=13$ | $t=256,s=12$ |
|--:|--:|--:|--:|
| 3 | 82.414 | 77.438 | 69.026 |
| 4 | 110.168 | 103.795 | 90.632 |
| 5 | 135.423 | 122.339 | 95.575 |
| 6 | 166.847 | 150.952 | 109.202 |
| 7 | 174.884 | 162.595 | 111.770 |
| 8 | 217.458 | 189.310 | 83.925 |

The shorter epoch remains strongest at every tested occupation. The
$t=64,s=14$ ladder remains positive through $Q=29$, where its margin is
30.683 bits, and becomes vacuous at $Q=30$. A 0.01-spaced local tilt check
confirms that this transition is not a coarse-grid artifact. The dominant
term has no all-one row; the current obstruction is the regular-row sparse
transfer.

These results are binary64 diagnostics. `RM2SUB_Q_LADDER.md` gives the exact
change of measure, recurrence, probability space, and reproduction files.

## Open proof obligations

The diagnostic multiband bridge is now complete for $t=64,s=14$. It uses the
refined groups 1--95, 96--416, and 417--512. The calculation checks every
composition for $30\le Q\le256$. Their union has 1228.638 bits, and its
consolidation audit passes over 2,857,249 compositions.

This closes the earlier missing region as a matter of binary64 diagnostics.
It does not yet give an outward certificate.

The complete diagnostic union selects the strongest valid receipt at every
occupation from Q=1 through 256. It has 42.577982 bits and clears the 40-bit
line by 2.577982 bits. Q=1 is the bottleneck.

The following obligations remain after that replay:

1. outward-round all surviving occupation bounds;
2. prove or instantiate `OuterCap` for extrapolated outer constituents;
3. define an arbitrary-length wrapper and quantify its rate and distance loss;
4. count XORs for the selected $A$, $B$, and field-multiply circuits;
5. benchmark the exact finalist implementations, one benchmark at a time; and
6. bind the implementation to the certified construction.

## Reproduction files

- `evaluate_rm2sub_primary_tranche.py` runs the parameterized $Q=1$ and $Q=2$
  screens.
- `rm2sub_primary_tranche_d100.json` contains the complete 121-witness $Q=1$
  grid for the planned first tranche.
- `rm2sub_primary_frontier_summary_d100.json` consolidates the RM$(4,9)$
  frontier and identifies every source receipt.
- `summarize_rm2sub_primary_tranche.py` regenerates that summary.
- `audit_rm2sub_primary_tranche.py` exactly rechecks every selected map and the
  summary structure.
- `rm2sub_primary_tranche_audit.json` records the passing exact audit.
- The generated `nested_rm2sub/t128_chain0/` bundle contains the unselected
  nested family. Commit `97c7a7e` preserves it; the Overleaf tip omits it to
  stay below the project file limit.
- `RM2SUB_PERSISTENCE_UNCONFOUNDED.md` records the corrected causal
  comparison and the witness-aggregation error.
- `RM2SUB_Q_LADDER.md` records the fixed-RM occupations 3 through 32 and the
  sparse-to-bulk boundary.
- `evaluate_rm2sub_q_ladder.py` and
  `audit_rm2sub_two_colour_recurrence.py` implement and check the new
  recurrence.
- `rm2sub_calibration_constituents/` contains the selected maps and exact
  spectra used in the screen.
- `RM2SUB_TWO_BAND_BRIDGE.md` records the positive two-band reduction and the
  intermediate mixed-family result.
- `RM2SUB_REFINED_BAND_BRIDGE.md` records the complete diagnostic Q=30--256
  bridge and its remaining outward-rounding obligation.
- `rm2sub_full_occupation_q1_q256_diagnostic.json` and
  `rm2sub_full_occupation_q1_q256_audit.json` record the full 42.577982-bit
  diagnostic union.
