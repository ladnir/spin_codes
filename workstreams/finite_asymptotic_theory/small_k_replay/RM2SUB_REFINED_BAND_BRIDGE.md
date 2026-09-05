# Refined-band RM2Sub bridge

## Outward closure update

The binary64 proof template described below has now been compiled into a
complete outward certificate. For the fixed-RM(4,9), (t=64,s=14) instance
at (k=2^{16}) and (N=2^{17}), the final authenticated union proves

\[
  d_{\min}\ge13{,}108
\]

except with probability below (2^{-42.5779817562755}). Occupation one is
the bottleneck. The outward Q=31--256 union retains 1293.444358402287 bits.

The complete statement, probability space, arithmetic rules, manifests, and
remaining obligations are in `RM2SUB_RM49_FINITE_CERTIFICATE.md`. Statements
later in this note that call outward rounding open describe the earlier
diagnostic stage and are superseded by this update.

## Question

The fixed-RM occupation ladder closed through \(Q=29\). The earlier
three-group calculation closed \(30\le Q\le52\), but its dense extension
failed. This note determines whether that failure reflects the construction
or the change-of-measure bound.

The target event is

\[
  \exists m\ne0:\quad
  \operatorname{wt}(\operatorname{Enc}(m))\le13107=\lfloor0.10N\rfloor,
  \qquad N=2^{17}.
\]

The standard comparison line is 40 failure-margin bits. The calculations in
this note are nearest-binary64 diagnostics, not outward certificates.

## Fixed construction and probability space

The message length is \(k=2^{16}\). One fixed RM\((4,9)\) constituent

\[
  C:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512}
\]

is repeated in all \(L=256\) outer rows. The constituent and its exact weight
spectrum are fixed. They are not resampled between rows.

Setup samples independent uniform coordinate permutations for the outer
rows. It also samples independent uniform permutations in the 512 transposed
regions. Each RM2Sub epoch receives an independent nonzero field scalar. The
audited RM2Sub \(A/B\) maps are fixed, with epoch length \(t=64\) and state
dimension \(s=14\). The failure probability is over the sampled routing and
field scalars.

## Refined spectrum partition

The exact nonzero RM spectrum is partitioned into three groups:

| group | weights | Bernoulli reference \(p_j\) | \(\log_2\eta_j\) | maximizing weight |
|:--|:--|--:|--:|--:|
| low | 1--95 | \(0.25\) | 119.966741 | 32 |
| central | 96--416 | \(0.5\) | 260.318572 | 416 |
| high | 417--512 | \(0.8677722630069483\) | 104.761160 | 420 |

Here \(\eta_j\) is the maximum counting-density ratio relative to the
displayed Bernoulli product measure. If \(A_w\) denotes the exact constituent
spectrum, then

\[
  \eta_j=
  \max_{w\in W_j:A_w>0}
  \frac{A_w}
       {\binom{512}{w}p_j^w(1-p_j)^{512-w}}.
  \tag{1}
\]

The low group deliberately retains \(p_{\rm low}=1/4\). Minimizing only its
pointwise envelope gives a sparser reference, but that reference weakens the
inner low-output bound. Extending the low group from weights 1--47 through
weight 95 does not increase its envelope at \(p=1/4\).

The high reference is a fixed decimal selected by binary64 envelope
minimization. The central reference is exactly \(1/2\). Validity uses only
the displayed probabilities and does not require either probability to be
optimal. The fixed schedule is used for every occupation in the reported
union.

## Positive three-group reduction

Fix an occupation \(Q\) and a composition

\[
  (q_0,q_1,q_2)\in\mathbb Z_{\ge0}^3,
  \qquad q_0+q_1+q_2=Q.
\]

The row locations contribute

\[
  \binom{L}{q_0,q_1,q_2,L-Q}.
  \tag{2}
\]

In one transposed region, let \(F_u(z)\) be the established positive RM2Sub
transfer for a uniform support of \(u\) live positions. The three reference
groups independently retain their candidate positions with probabilities
\((p_0,p_1,p_2)\). Conditional on the total live weight, the region permutation
makes the union a uniform support. Therefore the reference transfer is

\[
  R_{q_0,q_1,q_2}(z)=
  \sum_{u_0=0}^{q_0}
  \sum_{u_1=0}^{q_1}
  \sum_{u_2=0}^{q_2}
  \left(\prod_{j=0}^2
    \binom{q_j}{u_j}p_j^{u_j}(1-p_j)^{q_j-u_j}
  \right)
  F_{u_0+u_1+u_2}(z).
  \tag{3}
\]

The calculation raises (3) across all 512 regions with a binary
log-semiring matrix power. It uses the exact audited spectrum of the RM2Sub
\(A\) map in the support-averaged live-state upper transfer. It also verifies
the \(B\)-kernel nonactivation probabilities from exact shell counts.

For each composition, the pointwise change of measure contributes

\[
  \eta_0^{q_0}\eta_1^{q_1}\eta_2^{q_2}.
  \tag{4}
\]

The evaluator chooses a Chernoff witness separately for every composition.
It then sums all compositions and occupations using positive log-sum-exp
operations.

The three-binomial collapse in (3) has an independent direct-enumeration
audit. Its largest absolute error is \(3.7748\times10^{-15}\). The binary
matrix-power audit has maximum logarithmic error
\(1.0658\times10^{-14}\).

## Diagnostic result

The calculation covers every composition at every occupation
\(30\le Q\le256\). It checks 2,857,249 compositions and 227 occupations.

| occupation or interval | diagnostic margin in bits |
|:--|--:|
| \(Q=30\) | 1228.638 |
| \(Q=53\) | 2134.716 |
| \(Q=128\) | 4866.548 |
| \(Q=192\) | 6297.419 |
| \(Q=232\) | 7173.271 |
| \(Q=248\) | 4983.658 |
| \(Q=256\) | 1611.249 |
| union \(30\le Q\le256\) | **1228.638** |

The union clears the 40-bit comparison line by 1188.638 bits. Occupation 30
is the weakest occupation.

The earlier dense failure was proof slack. With groups 1--47, 48--416, and
417--512 and references \((1/4,1/2,3/4)\), the Q=224 bound was \(-3665.254\)
bits. Moving weights 48--95 into the low group does not enlarge the low
envelope. It reduces the central envelope enough to give 6998.852 bits at
Q=224 and 1611.249 bits at Q=256.

As a cross-check, the exact nine-band pure-face evaluator closes every pure
face at Q=256. Its weakest pure face has 2015.626 bits. A separate split of
the old central group at weight 95 gives 2488.774 bits for the union of its
mixed Q=256 faces.

## Full occupation union

The strongest established diagnostic receipt is selected separately for each
occupation:

| occupations | selected reduction |
|:--|:--|
| 1 | exact-spectrum Q1 evaluator |
| 2 | exact fixed-spectrum Q2 evaluator |
| 3--4 | positive two-colour RM ladder |
| 5--29 | three-group bridge |
| 30--256 | refined three-group bridge |

The resulting Q=1--256 union has 42.577982 diagnostic bits. It clears the
40-bit comparison line by 2.577982 bits. Occupation one is the bottleneck,
with 42.577982 bits before the negligible addition of the other occupations.

This result is a complete binary64 proof template for the displayed
construction at ten-percent distance. It is not an outward certificate.

## Status and remaining proof obligations

The following items are proved as exact reductions:

- the spectrum partition and counting-density formula (1);
- the row-location factor (2);
- the three-binomial transfer identity (3);
- the pointwise density bound (4); and
- the positive union over all displayed compositions and occupations.

The direct audits validate the three-binomial collapse, binary matrix power,
composition ordering, occupation aggregation, and complete Q=30--256
coverage.

The numerical margin is not yet a formal certificate. The remaining items
are:

1. replace nearest-binary64 arithmetic with directed outward rounding;
2. authenticate the finite Chernoff witness table in the outward checker;
3. state the resulting end-to-end distance theorem only after the outward
   audit passes.

No additional band decomposition is currently required for Q=30--256.

## Reproduction files

- `probe_rm2sub_three_group_bridge.py` implements the fixed-composition
  reduction.
- `audit_rm2sub_three_group_bridge.py` checks the three-binomial collapse and
  binary matrix power.
- `consolidate_rm2sub_refined_band_bridge.py` audited the raw chunk receipts
  and compacted their occupation summaries.
- `rm2sub_refined_band_bridge_q30_q256_diagnostic.json` is the compact
  diagnostic receipt.
- `rm2sub_refined_band_bridge_audit.json` records `PASS` for all 2,857,249
  compositions and all 227 occupations.
- `consolidate_rm2sub_full_occupation_diagnostic.py` combines the strongest
  receipt for every occupation.
- `rm2sub_full_occupation_q1_q256_diagnostic.json` records the 42.577982-bit
  full diagnostic union.
- `rm2sub_full_occupation_q1_q256_audit.json` records `PASS` for exact
  Q=1--256 coverage.
- `rm2sub_fixed_rm_dense_pure_bands_q256_probe_d100.json` records the
  nine-band pure-face cross-check.
- `rm2sub_central_split_95_q256_probe_d100.json` records the old-central-group
  split cross-check.

The raw per-composition chunk receipts were removed after the consolidation
audit because they occupied about 800 MB. The compact receipt records their
SHA-256 hashes. The probe can regenerate them from the displayed parameters.
