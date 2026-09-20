# Riffle DP-2Lap g=4: histogram-transition checkpoint

## Question

The coefficient-histogram lift represents all 15 packet-value sums as Walsh
coefficients of one 16-bin histogram. This checkpoint asks whether that
histogram supports a smaller transfer argument for the 24-node distance.

The answer is negative for the histogram alone. The exact ordered-state
identity remains useful, but it shows which information the histogram loses.

## Ordered-state identity

Let \(A\) denote one zero-input state step, and let \(T:=A^{\mathsf T}\).
Fix a nonzero character \(\chi\in\mathbb F_2^{64}\). Write

\[
\chi_t:=T^t\chi
\qquad (1\le t\le24).
\]

View each \(\chi_t\) as 16 ordered nibbles. Let
\(a_{t,s}(\chi)\in\mathbb F_2^4\) be nibble \(s\) of \(\chi_t\). Then the
output bit for packet value \(v\), node \(t\), and slot \(s\) is

\[
\langle a_{t,s}(\chi),v\rangle.
\]

Therefore, the aggregate coefficient histogram has the concrete form

\[
n_a(\chi)
=\#\{(t,s):a_{t,s}(\chi)=a\}.
\]

The audit verifies this identity on every slot and node for all 64 basis
characters. Linearity extends the identity to every character.

This identity removes an artificial layer from the problem. The histogram is
the occupation count of nibble values along the ordered orbit
\(T\chi,\ldots,T^{24}\chi\).

## The histogram is not a transfer state

An exact transfer argument would need to determine the next histogram from
the current one. The following two states refute that property:

\[
x=\mathtt{0x7373d52af3615f19},
\qquad
y=\mathtt{0x953f133d5f76721a}.
\]

They have the same current nibble histogram:

\[
(0,2,1,3,0,2,1,2,0,1,1,0,0,1,0,2).
\]

Their successor histograms differ:

\[
\begin{aligned}
h(Tx)&=(0,2,1,2,0,2,0,1,1,2,1,0,2,0,1,1),\\
h(Ty)&=(1,1,1,1,2,2,0,2,1,3,0,0,0,2,0,0).
\end{aligned}
\]

Thus, nibble positions affect the next step. A transfer proof must retain
ordered information or a stronger quotient of the 64-bit state.

## First uniform local constraints

For each packet value, the observation ranks after one, two, and three nodes
are 16, 32, and 48. Hence these output blocks realize every binary pattern of
the corresponding length. No inequality based on at most three consecutive
nodes can restrict the output.

Every five-node observation code has parameters \([80,64]\). The audit
constructs its 16-dimensional dual, enumerates all dual words, and applies
the MacWilliams identity exactly. The resulting five-node two-sided
distances are only 2 or 3.

Suppose a relaxation retains only the two-sided distance of each consecutive
five-node window. Four disjoint five-node windows then give a 24-node bound
between 8 and 12. These bounds are exact for the scalar relaxation: placing
all required weight at nodes 5, 10, 15, and 20 satisfies every sliding
five-node inequality. The required 24-node bounds are 97 for values 1 through
14 and 72 for value 15.

The complete five-node spectra do not give a small enumeration either. A
non-15 word with 24-node two-sided weight at most 96 has a disjoint five-node
block with two-sided weight at most 24. Depending on the packet value, the
exact number of such five-node words is approximately

\[
2^{52.88706878}.
\]

For value 15, the corresponding threshold is 17. Its exact candidate count
is approximately \(2^{41.92674088}\). Neither set is a practical direct
transition frontier.

## Conclusion

The histogram lift is an exact description of the final Walsh coefficients,
but it is not a smaller dynamical state. Its simplest recurrence-based
relaxation loses almost all of the required distance. Exact five-node spectra
do not repair that loss.

This checkpoint does not refute the full 24-node distance lemma. It refutes
the proposed histogram-only proof route. A successful structural proof now
needs a new ordered-state quotient or an algebraic inequality that uses the
component decomposition without enumerating all supports.

The justified fallback is to resume the exact component frontier at
dimension 33. That work preserves progress toward a certificate while a
stronger ordered-state invariant remains open.

Artifacts:

- `explorations/riffle_dp_2lap_g4_histogram_transition.json`;
- `scripts/analyze_riffle_dp_2lap_g4_histogram_transition.py`.
