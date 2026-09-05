# RM2Sub calibration and rate-half family replay

## Result

The calibrated epoch length is $t=64$. At fixed first-order persistence
$s+\log_2 t=26$, the occupation-one calculation prefers $t=64,s=20$ to
$t=128,s=19$ and $t=256,s=18$ on both tested BCH constituents and at both
tested message lengths. The comparison therefore retains $t=64$.

For the message-length scan, the neutral state schedule is

\[
  s(e)=\max(7,e-4),\qquad k=2^e.
\]

For $e\ge 11$, this schedule satisfies

\[
  s(e)+\log_2 64=e+2.
\]

It therefore matches the first-order live-state persistence of the earlier
RandomStepConv schedule $M(e)=e+2$. The floor at $s=7$ is an implementation
constraint of the selected RM(2,6)-subcode family: its basis contains the
constant word and all six linear monomials.

With this schedule, RM(4,9) is the only authenticated structured constituent
that crosses the 40-bit occupation-one screen. It has at least 40 diagnostic
bits for $14\le e\le18$, with a maximum of 41.462 bits at $e=16$. Extended
BCH $[128,64,22]$ reaches 33.966 bits at $e=16$ and does not cross 40.
These statements concern the occupation-one binary64 calculation. They are
not full distance certificates.

## Code and probability space

Fix a binary linear constituent $C:\mathbb F_2^K\rightarrow\mathbb F_2^B$
of rate one half. The outer map applies the same $C$ to each of

\[
  L=k/K
\]

message rows. Thus $N=BL=2k$. The deterministic BCH and RM experiments use
one authenticated constituent. A random control samples one uniform full-rank
$[B,K]$ map once and repeats it in every row. Its occupation-one result is a
joint first moment over that sampled map; it is not a realized-spectrum claim.

Routing samples an independent uniform coordinate permutation for every outer
row and an independent uniform position permutation in every transposed
region. These permutations and the inner randomness are sampled once and
then fixed as part of the code.

One RM2Sub epoch has input $X\in\mathbb F_2^t$, state
$Q\in\mathbb F_2^s$, output $Y\in\mathbb F_2^t$, and update

\[
  Y=X+A(Q),\qquad Q'=\alpha Q+B(X).
\]

Here $A:\mathbb F_2^s\rightarrow\mathbb F_2^t$ and
$B:\mathbb F_2^t\rightarrow\mathbb F_2^s$ are fixed linear maps satisfying
$BA=0$. Every epoch samples an independent
$\alpha\in\mathbb F_{2^s}^{\times}$. The maps used by the retained 64-bit
schedule are subcodes of RM(2,6). Their spectra, the spectra of $\ker B$, and
the identity $BA=0$
are checked exactly for each state dimension.

The RandomStepConv comparison is the two-state transfer abstraction used by
the existing evaluator. At every bit, a live state returns to zero with
probability $2^{-M}$, and a live output bit is uniform. It is not an
enumeration of arbitrary nonlinear transition functions. Matching
$M=s+\log_2t$ equates the leading reset rate only; it does not make the two
transition laws identical.

## Native admissible lengths

The present RM2Sub transfer splits every transposed region into complete
64-bit epochs. A rate-half constituent is natively admissible when

\[
  K\mid k\quad\text{and}\quad 64\mid L=k/K.
\]

No padding, partial epoch, or shortening wrapper is used. For power-of-two
$k$ and $K$, the second condition is $e\ge\log_2K+6$. The scan records
every omitted case and its reason. An arbitrary-length wrapper remains an open
construction task.

## Calibration calculations

All displayed margins bound the event

\[
  \operatorname{wt}(Y)\le\lfloor0.10N\rfloor
\]

on a finite Chernoff grid. Each shell uses its own optimizing tilt.

### Epoch length at persistence exponent 26

The three configurations satisfy $s+\log_2t=26$.

| Outer | $k$ | RandomStepConv $M=26$ | RM2Sub 64/20 | RM2Sub 128/19 | RM2Sub 256/18 |
|:---|---:|---:|---:|---:|---:|
| BCH $[32,16,8]$ | $2^{13}$ | 11.504 | 11.255 | 10.532 | 9.688 |
| BCH $[32,16,8]$ | $2^{15}$ | 8.642 | 9.802 | 9.265 | 8.588 |
| BCH $[128,64,22]$ | $2^{13}$ | 39.329 | 38.776 | 38.268 | inadmissible |
| BCH $[128,64,22]$ | $2^{15}$ | 37.582 | 37.500 | 37.335 | 37.025 |

The selected $t=64,s=20$ map has $d(A)=16$ and
$d(\ker B)=8$. The $t=128,s=19$ map has $d(A)=48$ and
$d(\ker B)=6$. The selected $t=256,s=18$ map has $d(A)=96$, but
$d(\ker B)=4$ and 576 kernel words of weight four. The result shows that
live-code distance alone does not choose the epoch length.

### State dimension at $t=64$

At $k=2^{13}$, the persistence-matched state is $s=9$, corresponding to
RandomStepConv $M=15$.

| Outer | Sector | RandomStepConv $M=15$ | RM2Sub $t=64,s=9$ |
|:---|:---|---:|---:|
| BCH $[32,16,8]$ | $Q=1$ | 0.192 | 3.912 |
| BCH $[32,16,8]$ | $Q=2$ | 3.718 | 7.505 |
| BCH $[128,64,22]$ | $Q=1$ | 28.862 | 31.174 |
| BCH $[128,64,22]$ | $Q=2$ | 63.509 | 61.135 |

Thus the persistence calibration is useful but not an equality theorem.
RM2Sub differs by between -2.38 and +3.79 bits in these four tests.
The $s=9$ map has $d(A)=28$, $d(\ker B)=4$, and 2352 weight-four kernel
words. Those kernel words do not prevent these $Q=1,2$ cases from closing,
but they remain relevant to higher occupations.

A larger-state sweep at $s=15,16,18,19,20$ shows diminishing returns on
BCH $[128,64,22]$. At $k=2^{15}$, the $Q=1$ margin rises from 37.007
bits at $s=15$ to 37.500 bits at $s=20$. BCH $[32,16,8]$ is more
sensitive: its margin rises from 7.282 to 9.802 bits.

## Persistence-matched family replay

The replay covers exact BCH-derived block lengths 8, 32, 64, and 128; exact
RM block lengths 8, 32, 128, and 512; and random full-rank controls at every
power-of-two block length from 8 through 1024.

| Family | Constituent | Maximum Q1 margin | Location | First point at least 40 bits |
|:---|:---|---:|:---|:---|
| BCH-derived | extended BCH $[8,4,4]$ | 20.361 | $e=8$ | none |
| BCH-derived | extended BCH $[32,16,8]$ | 4.195 | $e=10$ | none |
| BCH-derived | shortened-XBCH $[64,32,12]$ | 13.036 | $e=15$ | none |
| BCH-derived | extended BCH $[128,64,22]$ | 33.966 | $e=16$ | none |
| RM | RM(1,3) $[8,4,4]$ | 20.361 | $e=8$ | none |
| RM | RM(2,5) $[32,16,8]$ | 4.195 | $e=10$ | none |
| RM | RM(3,7) $[128,64,16]$ | 14.195 | $e=14$ | none |
| RM | RM(4,9) $[512,256,32]$ | 41.462 | $e=16$ | $e=14$ |
| random | $[128,64]$ | 25.365 | $e=14$ | none |
| random | $[256,128]$ | 62.318 | $e=16$ | $e=13$ |
| random | $[512,256]$ | 137.647 | $e=18$ | $e=14$ |
| random | $[1024,512]$ | 289.734 | $e=20$ | $e=15$ |

The $e<11$ points use the $s=7$ floor and are stronger than an exact
persistence match to $M=e+2$. They should not be used to compare the two
inner models. This caveat affects the reported maximum for the length-8
structured curves.

At $e=20$, where both schedules use persistence exponent 22, the direct
BCH $[128,64,22]$ cross-check gives 22.102 bits for RandomStepConv and
32.338 bits for RM2Sub. The difference survives identical cutoffs and tilt
grids. It comes from the full epoch transfer, including the audited live
output spectrum; matching the reset exponent does not remove that effect.

## Status of the statements

The following facts are exact for the recorded inputs:

- the outer BCH and RM weight spectra already authenticated by the replay;
- every selected RM2Sub $A$ spectrum and $\ker B$ spectrum;
- full rank of $A$, the identity $BA=0$, and the recorded minimum
  distances for the selected maps;
- the native divisibility rule; and
- the combinatorial without-replacement placement law used for $Q=1$ and
  $Q=2$.

The numerical margins are conditional diagnostics. They use nearest binary64
arithmetic and a finite tilt grid. They are not outward-rounded certificates.
The random-outer entries are joint $Q=1$ expectations and do not assert that
one sampled constituent simultaneously realizes its expected spectrum.

The remaining proof obligations are:

1. outward-round the selected $Q=1$ and $Q=2$ calculations;
2. cover every occupation $Q\ge3$, retaining repeated-outer relation types;
3. determine whether the low-state weight-four kernels require a separate
   sparse-occupation split;
4. prove or replace the sampled-map selection rule if the setup distribution
   must produce a good $A/B$ pair with a quantified probability;
5. specify an arbitrary-length wrapper and include its rate and distance loss;
6. measure the actual XOR cost of $t=64,s(e)$, including multiplication by
   the epoch scalars; and
7. only after these steps, combine the occupation bounds into a distance
   failure probability.

## Reproduction files

- `evaluate_inner_calibration.py` compares $t=128$ state sizes with both
  stored-state and persistence-matched RandomStepConv controls.
- `generate_rm2sub_calibration_constituent.py` selects and exactly audits the
  $t=64$ and $t=256$ calibration maps.
- `evaluate_rm2sub_epoch_calibration.py` compares epoch sizes.
- `evaluate_rm2sub_t64_state_calibration.py` compares state dimensions.
- `evaluate_rm2sub_q2_calibration.py` performs the exact-spectrum $Q=2$
  comparison.
- `evaluate_rate_half_family_curves_rm2sub.py` reproduces both the fixed-$s$
  and persistence-matched family scans.
- `audit_rm2sub_calibration_replay.py` checks every locally selected $A/B$
  pair, exact MacWilliams duality, the schedule, admissibility, JSON/CSV
  consistency, and the reported numerical anchors.
- `rm2sub_calibration_replay_audit.json` is the passing audit receipt.
- `rate_half_family_k_margin_d100_rm2sub_t64_matched.json` and `.csv` are the
  principal replay receipts.
- `rm2sub-family-curves.html` is the compact curve view.

The fixed-$s=16$ replay is retained as a separate envelope in
`rate_half_family_k_margin_d100_rm2sub_t64_s16.json` and `.csv`. It must not
be confused with the persistence-matched replay.
