# Goal 02: shifted low-source target spectrum

## Question

Fix a message \(x\) whose extended-BCH encoding has weight 22. For
\(64\le e<16{,}448\), how small can the BCH weight of \(\gamma^e x\) be?

## Result

**The exact minimum target weight is 30.**

The complete scan used all 243,840 minimum-weight messages and all 16,384
shifted exponents. It evaluated 3,995,074,560 ordered pairs. The target
histogram contains

\[
N_{22}=N_{24}=N_{26}=N_{28}=0,
\qquad
N_{30}=1{,}401.
\]

Therefore, among one-data-symbol profiles whose repeated source weight is 22,
the minimum profile is

\[
(22,22,30),
\]

with total binary weight 74. In particular, the logically possible profiles
\((22,22,24)\), \((22,22,26)\), and \((22,22,28)\) do not occur.

## Comparison with uniform alpha

Independent uniform nonzero coefficients would give only 3.591 expected
weight-30 pairs. The shifted schedule has 1,401, about 390.1 times that
expectation. All weight-30 hits occur at exponents 64 through 89 and 106
through 142.

Thus the shift removes the most dangerous correlations but does not produce a
fully random low tail. The ordinary bulk histogram remains centered at weight
64.

## Scope

## Weight-24 source shell

The authenticated support-at-most-15 records generate all 6,855,968
weight-24 messages through 532 affine orbits. Every affine image was
re-encoded, and the total equals the committed coefficient \(A_{24}\).

The second exact scan evaluated

\[
6{,}855{,}968\cdot16{,}384=112{,}328{,}179{,}712
\]

ordered pairs. Its minimum target weight is 24. The histogram contains 1,175
weight-24 targets and no weight-22 target. Hence \((24,24,24)\) occurs, while
\((24,24,22)\) does not.

## Global one-data minimum

Every nonzero BCH word has even weight at least 22. A one-data-symbol profile
has the form \((h,h,h')\).

- If \(h=22\), the first scan gives \(h'\ge30\), so the total is at least 74.
- If \(h=24\), the second scan gives \(h'\ge24\), so the total is at least 72.
- If \(h\ge26\), the BCH minimum gives \(2h+h'\ge74\).

The 1,175 observed \((24,24,24)\) profiles attain the lower bound. Therefore,

\[
\boxed{d_{\mathrm{one\text{-}data}}=72}.
\]

This is the exact minimum over one-data-symbol outer inputs. It is not the
minimum distance of the complete outer code, which also contains inputs with
several nonzero data symbols.

## Reproduction

- `receipts/goal02_weight22_target_spectrum.json`;
- `receipts/goal02_weight24_affine_enumeration.json`;
- `receipts/goal02_weight24_target_spectrum.json`;
- `../../scripts/analyze_riffle_shiftalpha64_weight22_targets.cpp`;
- `../../scripts/analyze_riffle_shiftalpha64_weight22_targets.exe`.
- `../../scripts/enumerate_ebch128_weight24_affine.cpp`.
