# Riffle DP-2Lap g=4: structural anti-cancellation checkpoint

## Question

Exact support certificates cover 66 of the 127 nonempty irreducible-component
supports. Continuing the frontier is valid but remains casework. This
checkpoint asks whether the 24-node and endpoint codes admit a uniform
reduction that avoids the remaining 61 support cases.

The checkpoint produced one useful exact reparameterization. It also refuted
the two simplest anti-cancellation arguments. It did not prove or refute the
full 64-dimensional local-distance lemma.

## Observation maps

Fix a packet value \(v\in\mathbb F_2^4\setminus\{0\}\). Let \(T\) be the
invertible transpose-state recurrence. For a character \(\chi\), one node
emits the 16 bits obtained from the 16 packet slots. Consecutive nodes emit
the same observations from \(T\chi,T^2\chi,\ldots\).

For every packet value, the observation ranks of the first one, two, and
three nodes are 16, 32, and 48. Five nodes always have rank 64. Thus, five
consecutive nodes determine the character exactly.

Four nodes behave differently:

- values `1,5,6,7,8,9,11,14` have rank 64;
- values `2,3,4,12,13,15` have rank 63;
- value `10` has rank 62.

This near-systematic form looked like a possible quotient reduction. Exact
kernel classification refutes that reduction. The nonzero four-node kernels
contain component supports of dimensions 61, 63, and 64. In particular, the
value-15 kernel contains a full-support character. The kernel is not confined
to the already-certified small component tranche.

Five-node injectivity remains useful as an exact structural fact. By itself,
it supplies only a nonzero-output condition. It does not imply the required
24-node weight.

## Complementary supports are not orthogonal

The dimension-32 frontier splits into complementary pairs

\[
(\mathtt{0x2d},\mathtt{0x52})
\quad\text{and}\quad
(\mathtt{0x34},\mathtt{0x4b}).
\]

If complementary component codes were orthogonal on each local window, their
separate certificates could constrain cancellation in the full code. Exact
cross-Gram calculations refute this premise.

The calculation covers both pairs, all 15 packet values, and both window
lengths. Every cross-Gram matrix is nonzero. Their ranks range from 30 to 32,
and they contain between 463 and 540 nonzero entries. Complement dimension
therefore gives no direct duality or orthogonality argument.

## Exact coefficient-histogram lift

Packet-value linearity gives a different representation. For each local
coordinate \(p\), define \(a_p(\chi)\in\mathbb F_2^4\) by

\[
\langle a_p(\chi),v\rangle
=\langle\chi,z_v(p)\rangle
\qquad\text{for every }v\in\mathbb F_2^4.
\]

Let \(n_a(\chi)\) count coordinates \(p\) with \(a_p(\chi)=a\). For a window
of length \(L\), the character sum for packet value \(v\) is

\[
B_v(\chi)
=L-2w_v(\chi)
=\sum_{a\in\mathbb F_2^4}n_a(\chi)(-1)^{\langle a,v\rangle}.
\]

Hence all 15 packet-value words are the nontrivial Walsh coefficients of one
nonnegative 16-bin histogram. Parseval gives

\[
\sum_{v\in\mathbb F_2^4}B_v(\chi)^2
=16\sum_{a\in\mathbb F_2^4}n_a(\chi)^2.
\]

The audit reconstructs these histograms for every full-code search witness
and verifies both identities exactly.

This lift removes packet value as an independent code construction. It also
makes cancellation explicit: adding two component characters xors their
four-bit coefficient words coordinate by coordinate. However, an arbitrary
16-bin histogram is not necessarily attainable. A proof still needs a
tractable description or relaxation of the attainable histograms.

## Full-code refutation search

A bit-parallel coordinate-descent search tested the full 64-dimensional
codes. It used 10,000 pseudorandom restarts per case, plus deterministic basis
and degree-one seeds. Each reported character and weight was replayed exactly.

The search found no violation. The smallest non-15 two-sided weight was 120,
attained by the degree-one character. The value-15 searches returned the
known tight weights 72 at 24 nodes and 36 at 12 nodes. Characters supported
on 60--64 dimensions produced heavier local optima.

This search is refutation evidence only. It does not establish that the full
code has non-15 distance 120 or that no violating character exists.

## Conclusion

The proposed small-kernel and complementary-support arguments are false. The
coefficient-histogram lift is exact and genuinely combines the 15 packet
values, but it does not yet bound the feasible histograms.

The next structural goal should be narrow: derive exact constraints or a
certified relaxation for the attainable coefficient histograms. The target
constraints are

\[
|B_v|\le190\quad(v\ne15),
\qquad
|B_{15}|\le240
\]

for the 24-node histogram, and \(|B_{15}|\le120\) for the endpoint
histogram. The degree-one histogram must remain feasible because it attains
the value-15 equalities.

If no useful histogram relaxation emerges, the fallback is the dimension-33
frontier. Its four supports require 12,216,115,184 information vectors per
implementation.

Artifacts:

- `explorations/riffle_dp_2lap_g4_structural_checkpoint.json`;
- `explorations/riffle_dp_2lap_g4_full_local_distance_probe.json`;
- `scripts/analyze_riffle_dp_2lap_g4_structural_checkpoint.py`;
- `scripts/probe_riffle_dp_2lap_g4_full_local_distance.py`.
