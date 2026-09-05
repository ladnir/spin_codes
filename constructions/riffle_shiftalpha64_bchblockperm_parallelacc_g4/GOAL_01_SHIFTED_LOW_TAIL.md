# Goal 01: shifted minimum-shell exclusion

## Question

Does the exponent window \([64,16{,}447]\) remove every weight-22
correlation found in the current window \([0,16{,}383]\)?

## Result

**PASS for the complete minimum-weight shell.**

The exact extended-BCH spectrum contains 243,840 words of binary weight 22.
Fifteen complete affine orbits generate all of them. Every generated word
was re-encoded, and the resulting count equals the committed spectrum
coefficient \(A_{22}=243{,}840\).

The exact lag scan tested all
\(243{,}840\cdot16{,}384=3{,}995{,}074{,}560\) pairs and found no
\((x,i)\) satisfying

\[
\operatorname{wt}(B(x))
=\operatorname{wt}(B(\gamma^{64+i}x))
=22
\]

for \(0\le i<16{,}384\). Thus a one-data-symbol outer word in the shifted
candidate cannot have BCH-weight profile \((22,22,22)\).

This statement is exact for the minimum shell. It does not yet bound profiles
\((22,22,h')\) with \(h'>22\), or profiles with first weight above 22.

## Next target

Bound the adjacent profiles beginning with \((22,22,24)\) and
\((24,24,22)\). Use explicit shifted-spectrum information for the lowest
shells and reserve a product envelope for the higher-weight remainder.

## Reproduction

- `ebch128_weight22_messages.bin` stores the complete sorted message set.
- `receipts/goal01_shifted_weight22_exact.json` records the exact scan.
- `../../scripts/enumerate_ebch128_weight22_affine.py` regenerates the set.
- `../../scripts/analyze_riffle_shiftalpha64_weight22_exact.py` repeats the
  intersection scan.
