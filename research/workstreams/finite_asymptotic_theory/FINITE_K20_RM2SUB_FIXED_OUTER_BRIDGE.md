# Returning the fixed \(k=2^{20}\) outer to RM2Sub-S19

## Outcome of the first bridge audit

The fixed-outer Toeplitz certificate does not transfer to RM2Sub by replacing
one suffix probability. Toeplitz convolution observes one new random
coefficient per output coordinate. RM2Sub-S19 observes only 19 syndrome bits
per 128-bit epoch.

At the current length, all RM2Sub syndromes have rank at most 314,640.
Consequently, the shortened code contains a subspace of dimension at least
733,936 that is silent in every epoch. RM2Sub is the identity on this silent
subcode for every multiplier schedule.

This fact is not a distance counterexample. A random code of the silent
subcode's rate has a binary64 random-code benchmark of 326,094 bits at the
11% cutoff. In contrast, the universal dimension-only systematic bound is
larger than one by 661,901 bits. These diagnostics show that the dimension is
compatible with excellent distance but cannot prove it. The missing object
is a weight-sensitive profile of the silent subcodes.

The silent sector has enough provisional budget for a coarse BA selection
factor. Charging \(B^2\) once in every one of the 8,832 rows costs

\[
 2L\log_2B=139{,}667.3155\text{ bits}.
\]

Subtracting this cost from the random-code benchmark leaves about 186,427
bits. This subtraction is only a budget diagnostic. It does not prove that
the structured syndrome map produces a uniformly random silent subcode.

## RM2Sub interface

Partition a routed outer word into \(E=N/128=16{,}560\) epochs:

\[
 x=(x_0,\ldots,x_{E-1}),
 \qquad x_e\in\mathbb F_2^{128}.
\]

The selected RM2Sub constituent supplies linear maps

\[
 A:\mathbb F_2^{19}\to\mathbb F_2^{128},
 \qquad
 B:\mathbb F_2^{128}\to\mathbb F_2^{19}.
\]

The source receipts prove that \(A\) is injective, \(B\) is surjective, and
\(BA=0\). Let \(\alpha_e\gets\mathbb F_{2^{19}}^*\) be independent. With
terminal state \(q_E:=0\), the reverse recurrence is

\[
 y_e:=x_e+Aq_{e+1},
 \qquad
 q_e:=\alpha_e q_{e+1}+Bx_e.
 \tag{1}
\]

Only the multipliers \(\alpha_e\) are random in this bridge. The repeated BA
outer and every routing permutation remain fixed as in the Toeplitz
certificate.

The identity \(BA=0\) makes the epoch syndrome

\[
 \sigma_e:=Bx_e
\]

independent of the entering state. This property gives RM2Sub its small
zero/live transfer. It also creates the silent subcode below.

## Silent-subcode lemma

Let \(V\subseteq\mathbb F_2^N\) be the fixed routed outer code. Define

\[
 V_{\mathrm{silent}}
 :=\{x\in V:Bx_e=0\text{ for every }e\}.
 \tag{2}
\]

For every \(x\in V_{\mathrm{silent}}\), recurrence (1) has \(q_e=0\) for all
\(e\). Therefore

\[
 \operatorname{RM2Sub}_{\boldsymbol\alpha}(x)=x
 \quad
 \text{for every multiplier schedule }\boldsymbol\alpha.
 \tag{3}
\]

Each epoch contributes at most 19 independent constraints. Hence

\[
 \dim V_{\mathrm{silent}}
 \ge \dim V-19E.
 \tag{4}
\]

For the parent and shortened codes, (4) gives

\[
 \dim V_{\mathrm{silent}}^{\mathrm{parent}}ge745{,}200,
\]

\[
 \dim V_{\mathrm{silent}}^{\mathrm{short}}ge733{,}936.
\]

Thus a finite RM2Sub certificate must prove

\[
 d(V_{\mathrm{silent}})>233{,}164.
 \tag{5}
\]

Condition (5) is necessary. Random multipliers cannot repair its failure.

## Why the Toeplitz prefix profile is insufficient

For Toeplitz convolution, a zero output prefix forces the routed input prefix
to be zero. The dimension \(\kappa_t\) therefore counts every word that can
activate at or after coordinate \(t\).

For RM2Sub, a zero syndrome suffix requires only \(Bx_e=0\). The kernel of
\(B\) has parameters \([128,109,6]\). A silent epoch can therefore contain
many nonzero input bits. The raw prefix dimension does not record their
weight.

For \(0\le r\le E\), define the suffix-silent subcode

\[
 V_r:=
 \{x\in V:Bx_{E-r}=\cdots=Bx_{E-1}=0\}.
 \tag{6}
\]

The dimensions of \(V_r\) identify how many words can delay RM2Sub
activation. They do not identify how much deterministic output weight those
words accumulate in the silent suffix. The finite proof needs a weighted
version of (6), not only its dimension.

## A useful uniform live-state bound

The live state still mixes strongly. Let

\[
 W_A(z):=\sum_{a\in\operatorname{im}A}z^{\operatorname{wt}(a)},
 \qquad 0\le z\le1.
\]

For every \(v\in\mathbb F_2^{128}\), the coset weight enumerator satisfies

\[
 \sum_{a\in\operatorname{im}A}z^{\operatorname{wt}(v+a)}
 \le W_A(z).
 \tag{7}
\]

To prove (7), expand the Hamming kernel
\(z^{\operatorname{wt}(v)}\) in the Walsh basis. Every Walsh coefficient is
nonnegative for \(0\le z\le1\). The subgroup coset sum is therefore largest
on the subgroup itself.

If \(q\) is uniform over \(\mathbb F_2^{19}\setminus\{0\}\), equation (7)
gives

\[
 \mathbb E_q[z^{\operatorname{wt}(v+Aq)}]
 \le \frac{W_A(z)}{2^{19}-1}.
 \tag{8}
\]

The exact selected-\(A\) spectrum makes the right side explicit. An
optimistic calculation that assumes an uninterrupted live state gives more
than 228,000 bits of single-word Chernoff margin across 16,559 live epochs.
It first exceeds 40 bits at 3,684 live epochs.

These figures are diagnostics. They ignore state termination and all
deterministic segments. They show that the live-state weight production is
not the apparent bottleneck.

## Exact certificate target

For \(z\in(0,1)\), define the transfer-weighted enumerator

\[
 \Phi_V(z)
 :=\sum_{x\in V\setminus\{0\}}
   \mathbb E_{\boldsymbol\alpha}
   [z^{\operatorname{wt}(\operatorname{RM2Sub}_{\boldsymbol\alpha}(x))}].
 \tag{9}
\]

Chernoff's inequality gives

\[
 \mathbb E_{\boldsymbol\alpha}[Z_{\mathrm{bad}}]
 \le z^{-D}\Phi_V(z),
 \qquad D=233{,}164.
 \tag{10}
\]

The next finite certificate should compute or upper-bound (9) by combining
three objects:

1. the weight enumerator of \(V_{\mathrm{silent}}\);
2. weight-sensitive profiles of the suffix-silent subcodes \(V_r\); and
3. the existing zero, deterministic, uniform-live, and punctured-live
   RM2Sub transfer.

The first object discharges the multiplier-independent case. The second
charges deterministic output before activation and between live intervals.
The third charges the output produced while the state is live.

A dimension-only prefix sum cannot replace these objects. It assigns weight
zero to nonzero vectors in \(\ker B\), including at least 733,936 dimensions
of shortened messages.

## Most promising next proof step

The smallest useful next experiment is a weighted syndrome-prefix evaluator
for the frozen outer. It should begin with a tractable relaxation:

1. retain exact suffix-syndrome constraints from (6);
2. attach a fugacity \(z\) to each silent input bit;
3. apply the uniform live bound (8) after activation; and
4. compare the result with the old occupation-based RM2Sub bound.

For the silent sector, the first candidate comparison should use the existing
\(B^2\) selected-BA spectrum majorant. The provisional budget above is large
enough to absorb that factor at full occupation. Sparse occupations should
continue to use the existing dedicated RM2Sub bounds.

If the relaxation closes, the proof can restore the exact four-state
termination transfer. If it fails badly, the failure will quantify the
spectrum property that the outer must supply.

## Audit artifacts and status

`audit_rm2sub_fixed_outer_bridge.py` verifies the selected \(A\) and \(B\)
maps, re-enumerates all \(2^{19}\) words of \(\operatorname{im}A\), and checks
the exact \(A\) spectrum. It also verifies \(BA=0\), the rank of \(B\), and
the silent-subcode rank budget.

The script writes `rm2sub_fixed_outer_bridge_audit.json`. Exact structural
checks and binary64 diagnostics are separated in that receipt.

The bridge is not yet a distance certificate. Equations (2)--(5) identify a
necessary condition. Equations (7)--(10) define a sufficient proof route once
the weighted syndrome-prefix profiles are bounded.
