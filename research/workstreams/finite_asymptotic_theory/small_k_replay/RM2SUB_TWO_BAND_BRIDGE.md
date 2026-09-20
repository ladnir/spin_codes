# Two-band sparse bridge for fixed-RM RM2Sub

## Question

The whole-spectrum sparse bound becomes vacuous at occupation (Q=30).
Pure-band probes do not identify a bad outer-weight shell. This experiment
tests every composition supported on the minimum-weight band and one other
spectrum band for (30\le Q\le52).

The standard comparison line is 40 failure-margin bits. The exploration
retains every signed margin, including values below that line.

## Fixed construction and probability space

The message length is (k=2^{16}), and the codeword length is (N=2^{17}).
One fixed RM\((4,9)) map is repeated in (L=256) outer rows. The inner is
the fixed audited RM2Sub map with epoch length (t=64) and state dimension
(s=14).

Setup samples independent row-coordinate permutations, independent region
permutations, and one nonzero field scalar per RM2Sub epoch. The failure
probability is over these setup objects. The outer map and the selected
RM2Sub (A/B) maps remain fixed.

The bad event is

\[
  \operatorname{wt}(\operatorname{Enc}(m))\le13107
\]

for at least one nonzero message (m).

## Two-band reduction

Fix two disjoint outer-weight bands (W_1,W_2\subseteq\{1,\ldots,512\}).
Let (A_w) be the exact RM spectrum. For a reference probability
(p_j\in[0,1]), define the density of band (j) relative to the Bernoulli
product measure by

\[
  M_j(x):=
  \mathbf1\{\operatorname{wt}(x)\in W_j\}
  \frac{A_{\operatorname{wt}(x)}}
       {\binom{512}{\operatorname{wt}(x)}
        p_j^{\operatorname{wt}(x)}
        (1-p_j)^{512-\operatorname{wt}(x)}}.
  \tag{1}
\]

Endpoint probabilities are used only when their support contains the complete
band. In particular, (p_j=1) represents the all-one word exactly.

Suppose (a) active rows use (W_1), and (b) active rows use (W_2).
There are

\[
  \binom{L}{a,b,L-a-b}
  \tag{2}
\]

choices for their row locations. In one transposed region, the region
permutation maps the two row sets to disjoint uniform position sets.
Reference bits independently survive with probabilities (p_1) and (p_2).

Let (F_j(z)) be the established positive RM2Sub region transfer for a
uniform support of exactly (j) one-bits. Conditional on the total number of
surviving bits, permutation symmetry makes their union a uniform support.
Consequently, the reference transfer for composition ((a,b)) is exactly

\[
  R_{a,b}(z)=
  \sum_{i=0}^{a}\sum_{j=0}^{b}
  \binom ai p_1^i(1-p_1)^{a-i}
  \binom bj p_2^j(1-p_2)^{b-j}
  F_{i+j}(z).
  \tag{3}
\]

This collapse avoids a multitype RM2Sub state space. A direct six-position
enumerator checks (3) for all type counts of total size at most four. The
largest absolute discrepancy is (4.8850\times10^{-15}). When
(p_1=p_2), (3) also agrees with the established homogeneous recurrence to
(2.3315\times10^{-15}).

For (h>1), independence of the reference rows gives

\[
\begin{split}
 &\mathbb E\left[
   \prod_{u=1}^{a}M_1(X_u)
   \prod_{v=1}^{b}M_2(Y_v)\mathbf1_E
 \right] \\
 &\qquad\le
 \mathbb E[M_1^h]^{a/h}
 \mathbb E[M_2^h]^{b/h}
 \Pr[E]^{1-1/h}.
\end{split}
\tag{4}
\]

The calculation minimizes (4) over the 1,281-point Holder grid and its
(h=\infty) endpoint. It selects the Chernoff witness separately for every
composition before summing positive terms.

## Diagnostic result

The standardized references are

\[
  p_{1:47}=0.25,\qquad
  p_{48:416}=0.5,\qquad
  p_{417:512}=0.75.
  \tag{5}
\]

The selected inner transfer uses the exact audited spectrum of (A) in its
support-averaged live-state envelope. The calculation also verifies the
nonactivation ratios from the exact kernel shell counts.

For every partner band, the table sums all mixed splits and every occupation
(30\le Q\le52). Pure faces are excluded.

| anchor band | partner band | interval margin | gap to 40 |
|:--|:--|--:|--:|
| 1--47 | 48--63 | 1172.882 | 1132.882 |
| 1--47 | 64--95 | 1204.624 | 1164.624 |
| 1--47 | 96--159 | 1238.539 | 1198.539 |
| 1--47 | 160--191 | 1240.857 | 1200.857 |
| 1--47 | 192--320 | 1240.857 | 1200.857 |
| 1--47 | 321--352 | 1240.857 | 1200.857 |
| 1--47 | 353--416 | 1238.539 | 1198.539 |
| 1--47 | 417--512 | 1716.025 | 1676.025 |

The union of all displayed two-band families has 1172.882 diagnostic bits.
Its weakest family is the minimum-weight band paired with weights 48--63.
Within that family, the weakest occupation is (Q=30), and its dominant
mixed split has 29 minimum-band rows and one partner row.

These results strongly identify the earlier gap as proof slack. A common
envelope-minimizing reference gave a vacuous bound for some pairs. Separate
references recover more than one thousand bits without changing any code
parameter.

## Scope and next obligation

Equations (1)--(4) are exact positive reductions. The direct audit checks the
implementation of (3). The displayed margins still use nearest-binary64
matrix arithmetic and finite witness grids. They are not outward
certificates.

This intermediate result covers only compositions supported on the
minimum-weight band and one partner band. The later refined three-group
calculation covers every composition for (30\le Q\le256). Its union margin
is 1228.638 diagnostic bits. See `RM2SUB_REFINED_BAND_BRIDGE.md`.

## Reproduction files

- `probe_rm2sub_two_band_bridge.py` implements (1)--(4).
- `audit_rm2sub_two_band_collapse.py` performs both recurrence comparisons.
- `rm2sub_two_band_collapse_audit.json` records `PASS`.
- `audit_rm2sub_two_band_bridge.py` checks the complete band partition,
  reference schedule, occupation coverage, and every positive aggregation.
- `rm2sub_two_band_bridge_audit.json` records `PASS`.
- `rm2sub_two_band_bridge_all_partners_holder_split_reference_probe_d100.json`
  is the consolidated receipt.
- The `rm2sub_two_band_bridge_low_*` receipts record the failed common-reference
  probes and the successful reference separation.
- `probe_rm2sub_two_band_dense_bridge.py` records the failed dense fallback.
  Its categorical conditioning cost is too large for mixed row types.
