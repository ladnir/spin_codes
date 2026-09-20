# Riffle DP-2Lap g=4: component-split checkpoint

## Question

The endpoint-aware support sweep does not scale to dimension 64. This
checkpoint tests whether a direct component split yields a practical exact
certificate for one full local code.

The split is exact and easy to validate. Naive pair enumeration does not
improve the exponential scaling, and generic SAT does not supply the missing
pruning.

## Exact nearest-pair formulation

Fix a packet value and a 24-node window. Split a component support into
disjoint masks \(S_L\) and \(S_R\). Their observation codes satisfy

\[
C_S=C_{S_L}\oplus C_{S_R}.
\]

Therefore, the minimum two-sided weight of \(C_S\) is

\[
\min_{(x,y)\ne(0,0)}
\min\bigl(
\operatorname{wt}(c_L(x)+c_R(y)),
384-\operatorname{wt}(c_L(x)+c_R(y))
\bigr).
\]

This identity turns component cancellation into an exact nearest-pair
problem. It does not itself reduce the number of pairs.

## Dimension-33 validation

The validation uses support

\[
\mathtt{0x53}=(1,2,10,20)
\]

and packet value 1. The balanced component split is

\[
\mathtt{0x40}\oplus\mathtt{0x13},
\qquad
20+13=33.
\]

The optimized C++ engine enumerated every nonzero pair:

\[
2^{20}2^{13}-1=2^{33}-1=8{,}589{,}934{,}591.
\]

The exact minimum two-sided weight is 120. This result strengthens the
earlier dimension-33 certificate for this one case, which proved only the
required lower bound of 97.

A separate Python audit reconstructs both component bases and all generator
words. It verifies the generated input, enumeration count, minimum witness,
source hashes, executable hash, and receipt hash.

## Full-dimension scaling

A balanced full-code split has dimensions 32 and 32. Exhaustive comparison
requires

\[
2^{32}2^{32}-1=2^{64}-1
\]

pairs for one packet value and one side-combined search.

The committed 16-thread dimension-33 run measured approximately
\(1.74\times10^9\) pairs per second. Extrapolating that rate gives about 335
years for one full-dimension packet-value case. This projection is optimistic
because the dimension-33 right list fits in cache.

The balanced right list also requires at least 192 GiB for codewords, or
224 GiB with one 64-bit index per word. An unbalanced \(44+20\) split reduces
the right list to 56 MiB, but it leaves the \(2^{64}\) pair count unchanged.

## Generic solver diagnostic

The prototype also encodes the split problem with native XOR clauses and a
sequential cardinality counter. CryptoMiniSat received ten seconds for each
side of the known dimension-33 case. Both calls returned `UNKNOWN`.

This timeout is diagnostic. It does not refute SAT-based certification.
However, generic SAT is already much slower than the existing exact
information-set certificate on the validation case.

## Conclusion

The component split is the correct interface for cancellation, but it is not
yet a scalable certificate algorithm. A viable full-dimension method needs a
rigorous filter or lower bound that removes roughly 20--25 bits of pair search
before exact comparison. Neither raw enumeration nor the tested generic SAT
encoding supplies that reduction.

The direct algorithmic rescue route therefore remains open but requires a new
pruning theorem. Continued support sweeping and naive full-code pair search
are both unsuitable as completion strategies.

The next proof effort should return to the weaker global obligation. It
should test whether total character bias over all 32,772 nodes can be bounded
without requiring every 24-node block to have distance 97. If that route also
fails, a randomized local mixer becomes the justified construction fallback.

Artifacts:

- `explorations/riffle_dp_2lap_g4_component_split_pairs_s53_v01.json`;
- `explorations/riffle_dp_2lap_g4_component_split_pairs_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_split_sat_s53_v01.json`;
- `scripts/certify_riffle_dp_2lap_g4_component_split_pairs.cpp`;
- `scripts/prepare_riffle_dp_2lap_g4_component_split_words.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_split_pairs.py`;
- `scripts/solve_riffle_dp_2lap_g4_component_split_sat.py`.
