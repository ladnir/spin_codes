# Fixed-RM occupation ladder for RM2Sub

## Question

This experiment continues the ten-percent finite-length comparison beyond
occupations one and two. It asks two separate questions.

1. At fixed persistence exponent

   \[
     p=s+\log_2t=20,
   \]

   how do the tested epoch lengths compare at occupations
   \(Q=3,\ldots,8\)?
2. How far can the strongest sparse transfer be continued before a different
   dense-occupation argument is required?

The calculations are nearest-binary64 diagnostics. They are not outward
certificates.

This phase explores the parameter space. It does not impose a 40-bit
acceptance threshold. Every signed margin remains useful baseline data. A
positive margin gives a nonvacuous bound for that occupation. A negative
margin identifies a limitation of the current bound.

## Fixed construction and probability space

The message length is \(k=2^{16}\), and the encoded length is \(N=2^{17}\).
The outer constituent is the fixed RM\((4,9)\) code

\[
  C:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512}
\]

with minimum distance 32 and its authenticated exact spectrum. The same map
\(C\) is used in every one of the \(L=256\) outer rows. No outer code is
resampled between rows.

Setup samples independent uniform coordinate permutations for the outer rows,
independent uniform permutations in the 512 transposed regions, and one
nonzero field scalar for every RM2Sub epoch. The selected \(A/B\) maps are
fixed audited maps. The failure probability is over these setup objects.

The bad event is

\[
  \operatorname{wt}(\operatorname{Enc}(m))\le
  \lfloor0.10N\rfloor=13107
\]

for at least one nonzero message \(m\).

## Exact outer reduction

Fix an occupation \(Q\). Every active outer row contains either the unique
all-one codeword or a regular nonzero, non-all-one codeword. If \(r\) rows are
all-one and \(q=Q-r\) rows are regular, their locations have multiplicity

\[
  \binom Lr\binom{L-r}{q}.
  \tag{1}
\]

For a regular row, let \(A_w\) be the exact RM spectrum and let \(U\) be the
uniform measure on \(\mathbb F_2^{512}\). After the row-coordinate
permutation, the counting-measure density relative to \(U\) is

\[
  M(x)=
  2^{512}\frac{A_w}{\binom{512}{w}}
  \quad\text{when }\operatorname{wt}(x)=w,
  \tag{2}
\]

for \(1\le w\le511\), and it is zero at weights 0 and 512. Consequently,

\[
  \mathbb E_U[M^h]
  =\sum_{w=1}^{511}
    2^{-512}\binom{512}{w}
    \left(2^{512}\frac{A_w}{\binom{512}{w}}\right)^h
  \tag{3}
\]

is computed directly from the exact spectrum.

Let \(E\) denote the low-output event under \(q\) independent uniform
reference rows and \(r\) forced all-one rows. Independence here is over the
independent local messages and row permutations; it does not require
independent outer codes. For every \(h>1\), Holder's inequality gives

\[
  \mathbb E_U\!\left[
    \prod_{i=1}^qM(X_i)\,\mathbf1_E
  \right]
  \le
  \mathbb E_U[M^h]^{q/h}
  \Pr_U[E]^{1-1/h}.
  \tag{4}
\]

The endpoint \(h=\infty\) recovers the old worst-density envelope. The
evaluator minimizes (4) over a fixed 1,281-point grid and this endpoint.

## Positive two-colour region recurrence

The earlier one-all-one calculation recovered a forced input by a finite
difference. Repeating that subtraction is numerically unstable. The current
evaluator instead tracks two disjoint position types directly:

- a regular position whose bit is Bernoulli one-half; and
- a forced position whose bit is one.

Suppose an epoch contains \(a\) regular positions and \(b\) forced positions.
If \(K_j(z)\) is the RM2Sub epoch transfer conditioned on input weight \(j\),
the epoch matrix is

\[
  E_{a,b}(z)=
  \sum_{j=0}^a
  \binom aj2^{-a}K_{b+j}(z).
  \tag{5}
\]

Let \(R^{(u)}_{a,b}\) be the average ordered product over \(u\) complete
epochs. When the next epoch receives \((c,d)\) positions, its coefficient is

\[
  \frac{
    \binom{ut}{a-c,\,b-d,\,ut-a-b+c+d}
    \binom{t}{c,\,d,\,t-c-d}
  }{
    \binom{(u+1)t}{a,\,b,\,(u+1)t-a-b}
  }.
  \tag{6}
\]

Equations (5)--(6) use only nonnegative additions and matrix products. They
therefore preserve the entrywise RM2Sub epoch upper bounds. The direct
enumeration audit checks all 15 pairs \(a+b\le4\) in a six-position toy
instance. Its largest absolute discrepancy is
\(1.6654\times10^{-15}\).

## Matched-persistence result

The following table uses a complete log-surprisal grid from -8 through 0 in
steps of 0.1. Every entry includes the sum over all possible numbers of
all-one rows.

| \(Q\) | \((t,s)=(64,14)\) | \((128,13)\) | \((256,12)\) |
|--:|--:|--:|--:|
| 3 | 82.414 | 77.438 | 69.026 |
| 4 | 110.168 | 103.795 | 90.632 |
| 5 | 135.423 | 122.339 | 95.575 |
| 6 | 166.847 | 150.952 | 109.202 |
| 7 | 174.884 | 162.595 | 111.770 |
| 8 | 217.458 | 189.310 | 83.925 |

These are failure-margin bits. For example, 82.414 denotes an upper bound of
\(2^{-82.414}\) for that occupation. The shorter epoch is strongest at every
tested occupation. This strengthens the occupation-one and occupation-two
evidence that \(s+\log_2t\) is only a tuning coordinate, not a proof
invariant.

## Extended sparse ladder

The \((t,s)=(64,14)\) calculation was extended through occupation 32. The
full two-colour calculation is positive through \(Q=29\). Selected boundary
values are

| \(Q\) | margin bits |
|--:|--:|
| 16 | 349.830 |
| 24 | 246.835 |
| 27 | 147.629 |
| 28 | 77.159 |
| 29 | 30.683 |
| 30 | -28.886 |
| 31 | -100.726 |
| 32 | -178.577 |

The \(Q=30,31,32\) entries use a log-surprisal grid of width 0.01 around the
observed saddles. Their failure is not a half-step grid artifact.

The endpoint class is also not responsible. In every displayed result, the
dominant mixture has zero all-one rows; mixtures containing an all-one row
have margins thousands of bits larger. The obstruction is the regular-row
change of measure together with the sparse two-state inner envelope.

## What is proved and what remains open

The following facts are exact mathematical reductions:

- the fixed-constituent counting density (2);
- the Holder bound (4);
- the two-colour recurrence (5)--(6); and
- the finite union over the number of all-one rows.

The audit independently checks the implementation of (5)--(6) on a small
instance. The selected \(A/B\) maps and their spectra are covered by the
existing exact tranche audit.

The displayed margins remain diagnostic because matrix arithmetic, special
functions, and witness selection use nearest binary64 values. The calculation
also covers only the listed occupations. It is not a distance certificate.

The sparse transfer should not be extended blindly into the bulk. Its
entrywise state relaxation loses rapidly after \(Q=29\).

The first dense-bridge probes sharpen this statement. A single Holder change
of measure for the complete RM spectrum is vacuous at \(Q=29\). A seven-band
positive decomposition with jointly optimized row-bit density closes every
pure face except the two lowest bands. At \(Q=29\), the pure-band margins are

\[
\begin{array}{c|rrrrrrrrr}
\text{weight band}&1{:}47&48{:}63&64{:}95&96{:}159&160{:}191&192{:}320&321{:}352&353{:}416&417{:}512\\
\hline
\text{margin}&-765.880&-141.613&400.311&1125.948&2029.035&2955.294&9847.791&12324.213&18438.631
\end{array}
\]

in bits. The dense minimum-weight face improves with occupation. It has
-29.391 bits at \(Q=52\), 4.780 bits at \(Q=53\), and 39.057 bits at
\(Q=54\). Thus this dense pure-face bound becomes nonvacuous at occupation 53.

Conversely, the exact sparse recurrence restricted to weights 32 through 47
has 1224.573 bits at \(Q=29\), 1232.043 bits at \(Q=30\), and 2142.283 bits at
\(Q=54\). The low shell is therefore not intrinsically fatal.

A subsequent two-band calculation covers every composition supported on
weights 1--47 and one other displayed band. The union of those mixed families
over \(30\le Q\le52\) has 1172.882 diagnostic bits. Separate Bernoulli
references 0.25, 0.5, and 0.75 remove the false common-envelope obstruction.

The refined three-group bridge subsequently covers every composition for
\(30\le Q\le256\). Its binary64 union margin is 1228.638 bits. Thus the
earlier boundary was proof slack from the coarse spectrum partition.
`RM2SUB_REFINED_BAND_BRIDGE.md` gives the reduction, audit, and remaining
outward-rounding obligation.

Every bridge probe is diagnostic and uses nearest binary64 optimization. The
final bridge must be outward rounded before the construction can enter the
certified frontier.

## Reproduction files

- `evaluate_rm2sub_q_ladder.py` implements the finite fixed-occupation bound.
- `audit_rm2sub_two_colour_recurrence.py` performs direct enumeration.
- `rm2sub_two_colour_recurrence_audit.json` records the passing audit.
- `audit_rm2sub_q_ladder.py` checks receipt parameters, complete all-one
  mixtures, matched-persistence ordering, and the boundary sign;
  `rm2sub_q_ladder_audit.json` records `PASS`.
- `rm2sub_q_ladder_rm49_t{64_s14,128_s13,256_s12}_d100.json` records the
  matched-persistence occupations 3 through 8.
- `rm2sub_q_ladder_rm49_t64_s14_q09_q16_d100.json` and
  `rm2sub_q_ladder_rm49_t64_s14_q17_q32_d100.json` record the extended sparse
  ladder.
- `rm2sub_q_ladder_rm49_t64_s14_q30_q32_fine_d100.json` records the fine-grid
  boundary check.
- `evaluate_rm2sub_fixed_rm_dense.py` records the failed one-shot Holder
  bridge in `rm2sub_fixed_rm_dense_t64_s14_probe_d100.json`.
- `probe_rm2sub_fixed_rm_dense_bands.py` and the `rm2sub_fixed_rm_dense_*`
  receipts locate the pure-face dense threshold.
- `rm2sub_sparse_minweight_bridge_probe_d100.json` records the restricted
  sparse low-shell result.
- `RM2SUB_TWO_BAND_BRIDGE.md`, `probe_rm2sub_two_band_bridge.py`, and
  `rm2sub_two_band_bridge_all_partners_holder_split_reference_probe_d100.json`
  record the two-band bridge.
- `RM2SUB_REFINED_BAND_BRIDGE.md` and
  `rm2sub_refined_band_bridge_q30_q256_diagnostic.json` record the complete
  diagnostic multiband bridge.
- `audit_rm2sub_two_band_collapse.py` and
  `rm2sub_two_band_collapse_audit.json` record the passing direct audit.
