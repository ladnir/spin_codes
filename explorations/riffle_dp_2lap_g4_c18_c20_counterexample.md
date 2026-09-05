# Riffle DP-2Lap g=4: exact counterexample to the 12-node mixed-component lemma

## Claim under test

Let \(T\) be the 64-bit zero-input transpose recurrence. For packet value
\(v\), observe a character \(\chi\) through the 16 packet slots and the first
12 node exponents:

\[
c_{\chi,v}(s,t)
=\langle\chi,T^{t+1}(v\mathbin{\ll}4s)\rangle,
\qquad 0\le s<16,\quad 0\le t<12.
\]

This is a 192-bit word. The proposed local proof required

\[
\min\{\operatorname{wt}(c_{\chi,v}),
192-\operatorname{wt}(c_{\chi,v})\}\ge48
\]

for every \(v\ne15\) and every nonzero character. The medium-sized target
restricted this claim to the 38-dimensional subcode
\(C_{18}\oplus C_{20}\), formed by the degree-18 and degree-20 irreducible
components of \(T^\mathsf{T}\).

## Exact counterexample

For packet value

\[
v=7
\]

and character

\[
\chi=\mathtt{0x852ac8fcc6b27c67},
\]

the 192-bit observation has weight 146. Its complement has weight 46.
Therefore,

\[
\min\{146,192-146\}=46<48.
\]

The character lies in \(C_{18}\oplus C_{20}\): applying the product
annihilator gives zero. Applying either factor alone gives a nonzero result.
Thus, this is a genuinely mixed character, not a character from either
individual component.

## Search certificate

The code has five disjoint 38-coordinate information sets, covering 190 of
the 192 coordinates. If a word or its complement had weight at most 47, its
restriction to at least one information set would have weight at most

\[
\left\lfloor\frac{47}{5}\right\rfloor=9.
\]

The fixed-width C++ search enumerates the radius-nine Hamming balls around
both zero and one. It found the counterexample after 138,217,485 candidates
from the first information set. Direct evaluation, rather than the
information-set argument, is sufficient to establish the refutation once the
witness is known.

An independent Python verifier reconstructs the recurrence and BCH map,
checks mixed-component membership, and replays all 192 observed bits. It does
not use the C++ information-set machinery.

## Global replay and implication

The same character was replayed over all 32,772 nodes, or 524,352 observed
bits. Its weight is 263,123, its complement weight is 261,229, and its
character sum is \(-1,894\). Its absolute global bias is therefore

\[
\frac{1,894}{524,352}\approx0.00361208.
\]

Among the 2,731 consecutive 12-node blocks, exactly one has two-sided weight
below 48. The counterexample therefore refutes the uniform 12-node lemma. It
does not refute the desired global cap \(B_7\le1/2\), the support-33
terminal-zero bound, or Riffle DP-2Lap g=4.

The failed block is followed by a block of weight 90. Across those 24 nodes,
the observation has weight 236 out of 384 and two-sided weight 148. This
supports replacing the per-block claim with a transition-aware or amortized
claim. A 24-node certificate is the smallest natural next target, but its
threshold must be chosen from the exact global accounting rather than by
simply doubling 48. There are 1,365 complete 24-node blocks and one remaining
12-node block. A universal 24-node two-sided distance 96 contributes 131,040
to the required global distance 131,088; the remainder would still need
distance 48, which is precisely the failed claim. Distance 97 on every
24-node block would suffice even with a trivial remainder bound. A sharper
phase-dependent argument could instead retain distance 96 and control the
remainder.

## Reproduction artifacts

- Search source: `scripts/certify_riffle_dp_2lap_g4_c18_c20.cpp`
- Search receipt: `explorations/riffle_dp_2lap_g4_c18_c20_v07.json`
- Independent verifier: `scripts/verify_riffle_dp_2lap_g4_c18_c20_counterexample.py`
- Verification receipt:
  `explorations/riffle_dp_2lap_g4_c18_c20_counterexample_verification.json`

The search receipt is labeled `EXACT_COUNTEREXAMPLE`. The verification
receipt is labeled `EXACT_COUNTEREXAMPLE_VERIFICATION` and records the full
192-bit complement support and full-orbit statistics.
