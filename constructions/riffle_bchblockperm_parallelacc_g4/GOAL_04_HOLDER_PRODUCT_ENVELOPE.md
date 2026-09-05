# Goal 04: Hölder product-envelope gate

## Question

Can a product envelope remove the deterministic \(\alpha_i\) joint spectrum
without enumerating the relation between BCH weights?

## Proposed inequality

Suppose the conditional inner failure bound satisfies

\[
Q(h,h,h')\le r(h)^2r(h').
\]

For \(a_x=r(\operatorname{wt}(B(x)))\), Hölder's inequality gives

\[
\sum_x a_x^2a_{\alpha x}\le \sum_x a_x^3
\]

for every nonzero deterministic \(\alpha\). This removes the joint BCH
spectrum from the upper bound.

## Result

**The inequality is useful, but the present local bound fails its minimum-shell
gate.**

The exact minimum-shell mass over all 16,384 data positions is

\[
\log_2(16{,}384 A_{22})=31.8955752824.
\]

Goal 02 certifies only
\(Q(22,22,22)\le 2^{-36.4472162884}\). An envelope derived from this certified
bound must take \(r(22)^3\) at least as large as that right-hand side. It then
charges

\[
16{,}384 A_{22}Q(22,22,22)=2^{-4.5516410060}.
\]

This is a limitation of the current certificate, not a lower bound on the
actual failure probability. To prove a 20-bit failure target through this
certificate, the local profile bound must
reach \(2^{-51.8955752824}\), an improvement of 15.4484 bits. A 40-bit target
requires an improvement of 35.4484 bits.

An optimistic full-global-bit-permutation calculation gives
\(2^{-56.3416504117}\) for fixed weight 66. After the same shell charge, it
reaches \(2^{-24.4460751293}\). Thus the product route is not structurally
incapable of a 20-bit result, but the present block-permutation certificate is
too weak. Even the optimistic comparison does not reach 40 bits.

## Consequence for the shifted candidate

The generic Hölder bound discards the exact fact that the shifted candidate
has no minimum-to-minimum pair. It should not be applied unchanged to that
shell. A better combined proof uses exact shifted-spectrum exclusions for low
weights and a product envelope only above a selected cutoff.

## Reproduction

- `receipts/goal04_holder_minimum_shell_gate.json`;
- `../../scripts/analyze_riffle_bch_holder_gate.py`.
