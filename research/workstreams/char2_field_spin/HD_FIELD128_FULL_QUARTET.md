# Composing four packets before the fractional power

The pair bound splits the state path after every two packets. This split
can lose information when the proof subsequently takes a fractional power.
The quartet calculation removes half of these splits without changing the
encoder, the route, or the setup distribution.

Fix a good global gamma, a mark `0 < z < 1`, and `0 < alpha <= 1`.
Let `P_j(z)` be the existing nonnegative, two-state packet cap at occupancy
`j`, before the fractional power. The states are zero and nonzero inner
state. The translated-feedback return proof supplies these caps pointwise
in every gamma satisfying the principal-invertibility event.

## Positive coefficient construction

A packet contains 32 routed pair slots. A slot can be empty, singleton,
or full. Use `x` to mark singleton slots and `y` to mark full slots. Define

```text
g_j(x,y) = sum_f binom(32,j-2f) binom(32-j+2f,f) x^(j-2f) y^f,
max(0,j-32) <= f <= floor(j/2).
```

The singleton orientations have already been averaged over the fair swaps.
No additional factor of two belongs in `g_j`.

Write `M^alpha` for the entrywise power of a nonnegative matrix here.
The pair and quartet polynomials are

```text
T(x,y) = sum_(a,b) g_a g_b (P_a P_b)^alpha,
Q(x,y) = sum_(a,b,c,d) g_a g_b g_c g_d (P_a P_b P_c P_d)^alpha.
```

Matrix products follow the packet order. Each product includes its state
path sum before taking alpha. Histogram multiplicities stay outside alpha.
Every occupancy from zero through 64 is retained.

For nonnegative numbers `u,v`,
`2^(alpha-1)(u^alpha+v^alpha) <= (u+v)^alpha <= u^alpha+v^alpha`.
Applying this inequality to the middle state gives coefficientwise bounds

```text
2^(alpha-1) T(x,y)^2 <= Q(x,y) <= T(x,y)^2.
```

There are 32 quartets per region and 256 regions. At fixed marks, this
refinement can save at most `8192(1-alpha)` bits in the complete transfer.
This ceiling neither predicts the gain nor proves a distance obstruction.

## Computing every coefficient without a large quartet tensor

Write `g_j = x^j h_j(y/x^2)`. Each `h_j` has at most 17 consecutive
coefficients; retain its minimum exponent rather than padding from zero.
Cache the 2,145 products `h_a h_b` with `a <= b`.

For each multiset `a <= b <= c <= d`, form all 24 positional matrix
orderings. Take alpha separately for each ordered product, sum them, and
divide by the factorials of repeated occupancy multiplicities. This is the
sum over distinct orderings, not the power of their sum.

Multiply this matrix sum by the short positive polynomial
`(h_a h_b)(h_c h_d)`. If its absolute exponent is `f`, place the term at
`(s,f)`, where `s=a+b+c+d-2f`. The 814,385 multisets cover all 17,850,625
distinct ordered occupancy quartets. The coefficient table has dimensions
`4 x 129 x 129`; degree 128 must not be truncated.

For a region with `s` singleton slots, `f` full slots, and `4096-s-f`
empty slots, divide `[x^s y^f] Q^32` by the corresponding multinomial
count. The regional exponent is 32, not the pair implementation's 64.
The original active-group weight remains in this normalization and in
the outer factor.

Within one fixed `(z, alpha)` witness, take entrywise minima of the old
pair bounds and new quartet bounds. Different witnesses can be minimized
only after their own complete chronological transfer and outer factor.

## Status and arithmetic

Versions 54--56 implement numerical diagnostics, not directed enclosures.
Version 54 independently sums all ordered products at four positive marks.
Version 55 constructs the positive coefficient polynomial and uses a
power-32 residue transform. Version 56 checks those four point evaluations
and replays every midpoint profile with all prior witnesses retained.

Raw matrix products remain in logarithms until after alpha. Ordinary
positive coefficient arithmetic rejects finite-to-zero conversions at
that stage. Later rounding and Fourier cancellation still require a
separate directed implementation before any certificate can be claimed.
The setup-failure probability remains a separate, once-global term.
