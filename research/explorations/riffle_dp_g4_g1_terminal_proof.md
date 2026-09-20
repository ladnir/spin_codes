# Terminal weight-three proof for Riffle DP g=4

## Result

The terminal-placement contribution of every outer block-weight-three word
is at most

\[
2^{-40.212996059540}.
\]

This upper bound is below the target by 0.212996 bits. It closes gate G1 for
the frozen **Riffle DP g=4@g0-v1** construction.

The result concerns only the terminal-placement event. It does not cover
other inner placements or outer block weights of at least four.

## Coefficient classes

Fix an outer support containing three field symbols. The two parity equations
leave a one-dimensional solution space over \(\mathbb F_{2^{64}}\). Thus the
nonzero words on this support have the form

\[
(at,bt,ct),\qquad t\in\mathbb F_{2^{64}}^*,
\]

for fixed nonzero coefficients \(a,b,c\). Apply the local BCH encoder and
split each 128-bit word into 32 nibbles. The resulting set, including the
zero word, is an additive code of length 96 over the alphabet
\(\mathbb F_2^4\). It has size \(2^{64}=16^{16}\).

There are

\[
\binom{16{,}386}{3}=733{,}141{,}975{,}040
\]

outer supports and therefore that many coefficient classes.

## Exact classification through support 38

The local BCH code has nibble distance 11. An exact shortened-code scan gives:

| Local support | Number of words |
|---:|---:|
| 11 | 20 |
| 12 | 1,526 |
| 13 | 37,014 |
| 14 | 742,031 |
| 15 | 13,430,995 |

The scan visited 1,846,943,452 packet subsets. Its maximum shortened
dimension was five. Every emitted message passed independent re-encoding.

The global classifier treats the four outer categories separately. The two
categories with repeated local values are exact through total support 36.
The two additive categories are exact through total support 38.

For an additive triple \((u,v,w)\), the parity equation gives
\(u+v+w=0\). If its total support is at most 38, its two smallest components
have supports at most 12 and 13. The support-13 local receipt therefore
determines the third component by XOR. No support-16 enumeration is needed.

The classifier stores every projective ratio between components of each low
additive triple. A coefficient triple matches exactly when one of its ratios
appears in this table. This test is invariant under a common nonzero scalar
and under permutation of the three components.

The exact and conservative proof buckets are:

| Minimum support used in the proof | Coefficient classes | Classification |
|---:|---:|---|
| 33 | 3 | exact |
| 35 | 11 | exact |
| 36 | 32,779 | exact |
| 37 | 134,225,892 | conservative repeated-category remainder |
| 38 | 49,143 | exact additive classes |
| 39 | 733,007,667,212 | exact lower bound on the remaining classes |

The support-37 bucket places every unclassified repeated category at support
37. This choice can only increase the final upper bound.

## Delsarte bound for one class

Let \(A_i\) count words of nibble support \(i\) in one coefficient class.
For \(0\le i,j\le96\), let \(K_j(i)\) be the \(16\)-ary Krawtchouk
polynomial. The MacWilliams transform of an additive code is nonnegative, so

\[
\sum_{i\ge d} A_i\frac{K_j(i)}{K_j(0)}\ge -1
\]

when the class has minimum support at least \(d\).

For packet support \(i\), define the terminal probability

\[
P(i):=\frac{\binom{47{,}184}{i}}{\binom{524{,}352}{i}}.
\]

The verifier finds exact rational values \(u\ge0\) and \(z\) satisfying

\[
\frac{P(i)}{P(d)}
\le
z-u\frac{K_j(i)}{K_j(0)}
\qquad(d\le i\le96).
\]

Consequently,

\[
\sum_{i\ge d}A_iP(i)
\le
P(d)\bigl((2^{64}-1)z+u\bigr).
\]

The proof uses degrees 45, 43, 42, 41, 40, and 39 for minimum supports
33, 35, 36, 37, 38, and 39, respectively. Every inequality and the final
sum use exact rational arithmetic.

The contributions of the six buckets are:

| Minimum support | Logarithm of contribution |
|---:|---:|
| 33 | at most \(-57.811822709066\) |
| 35 | at most \(-62.684950871807\) |
| 36 | at most \(-54.516206338312\) |
| 37 | at most \(-45.891870591388\) |
| 38 | at most \(-60.687049654357\) |
| 39 | at most \(-40.241517984266\) |

Their exact rational sum lies in

\[
\left[2^{-40.212996059541},2^{-40.212996059540}\right].
\]

## Refutation status

The authenticated support-33 and support-35 terminal families contribute

\[
2^{-109.944337962560}
\]

in total. They remain far below the refutation threshold. The classifier
also identifies new low coefficient classes at supports 36 and 38. Those
classes should seed the nonterminal placement search in gate G2.

## Reproduction

The proof artifacts are:

- `explorations/riffle_dp_g4_construction_manifest.json`;
- `bch_g4_le15.bin`;
- `explorations/riffle_dp_g4_support38_classification.txt`;
- `explorations/riffle_dp_g4_g1_terminal_ledger.json`;
- `scripts/enumerate_bch_g4_low_support.cpp`;
- `scripts/classify_riffle_dp_g4_support38.cpp`;
- `scripts/riffle_dp_g4_weighted_triples.py`;
- `scripts/verify_riffle_dp_g4_g1.py`.

Run the independent verifier with:

```powershell
python scripts\verify_riffle_dp_g4_g1.py
```

The expected final line is:

```text
status=EXACT_G1_LEDGER_VERIFIED
```
