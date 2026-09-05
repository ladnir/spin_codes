# Distance certificate for the g=4, B=1024, sigma=15 instance

This certificate concerns the random construction specified in
`../CONSTRUCTION.md`. It certifies the complete first-moment sum, not only its
sparse occupations.

## Statement

Fix

\[
n=2^{20},\qquad B=1024,\qquad g=4,\qquad \sigma=15.
\]

The construction has \(L=2048\) outer blocks, \(P=512\) fixed packet groups,
and output length

\[
N=2^{21}.
\]

Let \(C\) denote the linear encoder sampled during setup. The setup samples
the outer injections, row permutations, and inner step maps specified in the
construction. Then

\[
\Pr_C\left[
  \exists m\ne 0:\operatorname{wt}(C(m))\le 188743
\right]
\le
2^{-195.283900426038}.
\]

Consequently, except with probability at most the displayed bound,

\[
d_{\min}(C)\ge 188744>0.09N.
\]

In particular, the setup failure probability is less than \(2^{-40}\).

## Bound

Group nonzero messages by their number \(a\) of active outer blocks. For each
\(a\in\{1,\ldots,2048\}\), the checker applies these steps.

1. Replace each uniform nonzero outer word by independent uniform bits. The
   checker includes the factor \((1-2^{-B})^{-a}\).
2. Apply the packing-domination lemma from `PACKING_DOMINATION.md`.
3. Evaluate the exact packet-order transition for one transposed row.
4. Apply a Chernoff bound at a fixed recorded parameter. The checker performs
   no optimization.
5. Compose the row transition through all \(B=1024\) rows.
6. Multiply by the exact number
   \(\binom{L}{a}(2^{B/2}-1)^a\) of messages in the occupation.

Let \(U_a\) denote the resulting upper bound. The outward-rounded sum gives

\[
\log_2\left(\sum_{a=1}^{2048}U_a\right)
\le -195.283900426038041150079401119.
\]

The one-active-block term is dominant. Its certified logarithm is
\(-195.283900426038041150079401119\). The two-active-block term is at most
\(2^{-437.351063284269215239496390839}\).

## Arithmetic

The row dynamic program stores each nonnegative value as a binary64 mantissa
and an integer binary exponent. After every mantissa addition and
multiplication, `nextafter` rounds toward positive infinity. Per-entry powers
of two prevent underflow in dense occupations.

The checker uses 80-decimal-digit interval arithmetic for scalar logarithms,
outer multiplicities, Chernoff corrections, and the final sum. The fixed
Chernoff parameters come from the exploratory receipt. Their optimality is
irrelevant to validity.

Run the checker from the repository root:

```text
python scripts/certify_riffle_transpose_packetshuffle_full.py --dps 80
```

The certificate receipt is
`../receipts/g4_b1024_sigma15_full_interval.json`.

The audited files have these SHA-256 digests:

```text
78925DDCCE9CA1A4EDBB59C2933B43FF7A51ACFE2C75387CC9B3BFB387C2E4A5  scripts/certify_riffle_transpose_packetshuffle_full.py
F952B4037E1E61D6F79199B874F0FBD48D750DDFC52A82240A0BD9860E79327F  receipts/g4_b1024_sigma15_full_interval.json
```

## Scope

The statement certifies the sampled random ensemble. It does not certify a
particular setup seed. It also retains the random rate-half outer injections;
it does not yet apply to a BCH replacement or another practical outer code.
