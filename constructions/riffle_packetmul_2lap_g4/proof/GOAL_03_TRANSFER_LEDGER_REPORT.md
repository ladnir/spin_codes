# Goal 03: transfer audit and partial probability ledger

## Result

The paused candidate contributes three proof ingredients unchanged: the
one-data outer population, the uniform packet-placement law, and the universal
three-node autonomous-prefix certificate. Goals 01 and 02 replace the paused
fixed-value character analysis and prove the support-33 terminal-zero row.

These results do not prove the construction. After a disjoint event split, the
first open row is the support-33 event with nonzero first-lap terminal state
and an occupied node before node 22,653. The exact closed partial ledger is in
`receipts/goal03_transfer_ledger.json`.

## Probability space

Fix an outer word (x). The setup samples one independent
(A_i\in\mathbb F_{16}^{\times}) for each input packet coordinate and an
independent uniform packet permutation \(\Pi\). Define

\[
y=\Pi(D_Ax).
\]

Let (L(y)) be the first-lap terminal state and let (G(y)) be the retained-
state two-lap output. An output is bad when

\[
\operatorname{wt}(G(y))<188{,}766.
\]

The multiplier (A_i) is nonzero. It therefore preserves whether packet (i)
is zero. It does not preserve the packet's nonzero value or its binary weight.
The permutation \(\Pi\), the map (G), and the outer stage are unchanged.

## Transfer audit

The following table classifies the inherited proof mechanisms relevant to the
current ledger.

| Mechanism | Class | Reason |
|---|---|---|
| One-data outer enumeration through support 39 | Transferred | PacketMul preserves packet support. The counts remain (26,36,3233,510,933,81090). |
| Uniform support-(h) placement | Transferred | \(\Pi\) is unchanged and independent of (A). Active cells form a uniform (h)-subset of 524,352 cells. |
| Three-node autonomous-prefix certificate | Transferred | The certificate is universal over every nonzero entering state and the deterministic map (G) is unchanged. |
| Final-10,119-node placement bound | Transferred | It uses only support, uniform placement, and the autonomous certificate. |
| Support-33 terminal-zero row | Transferred with substitution | Goals 01 and 02 replace fixed-value character bounds with the PacketMul reduction and (q_\chi\le5/8). |
| Paused fixed-value component program | Superseded | Its target was the terminal-zero row now proved by Goals 01 and 02. Its fixed-value character averages are not the PacketMul averages. |
| Fixed-value and binary-weight probes | Invalidated | Multiplication by a nonzero field element changes both quantities. |
| Bounded fixed-multiset turnoff exclusion | Invalidated | Independent multipliers can change the packet-value multiset before convolution. |
| Authenticated six-packet turnoff family | Transferred with substitution | The same target drive word is realized by forcing one multiplier at each of its 33 active coordinates. |
| Support-33, nonzero-terminal, nonsuffix event | Open | Neither the terminal-zero lemma nor the long zero-prefix argument applies. |
| Remaining outer population | Open | The transferred enumeration covers only the authenticated one-data words through support 39. |

The machine-readable ledger records the artifacts and their SHA-256 digests.
It also authenticates the independent verifiers for the suffix bound and the
six-packet turnoff.

## Autonomous-prefix lemma

For a zero-input portion of the second lap, let (o_j) be the output at node
(j). The inherited certificate proves that every nonzero consecutive
three-node window has weight at least 25. The autonomous recurrence is
invertible, so a prefix entered from (L(y)\ne0) remains nonzero in every
complete three-node block.

A prefix of 22,653 nodes contains 7,551 complete blocks. Its weight is at least

\[
7{,}551\cdot25=188{,}775>188{,}766.
\]

Thus, if all active packet cells lie in the final 10,119 nodes, a bad output
must satisfy (L(y)=0). The cutoff is minimal for this particular block
argument: 22,652 nodes contain only 7,550 complete blocks, whose certified
weight is 188,750.

For packet support (h), the probability of the suffix-placement event is

\[
\frac{\binom{161{,}904}{h}}{\binom{524{,}352}{h}}.
\]

The original exact aggregate over all 85,828 one-data words is at most

\[
2^{-48.293011187780}.
\]

This original aggregate is a valid transferred bound, but the support-33 term
must not be added to the new terminal-zero bound because those two charges can
overlap.

## Disjoint partial ledger

For each support-33 word, split bad outputs into the cases (L(y)=0) and
(L(y)\ne0). In the second case, split again according to whether all active
packet cells lie in the final 10,119 nodes. The suffix part is empty by the
autonomous-prefix lemma. For each word of support 35 through 39, split only by
the suffix-placement event. These choices give the following disjoint rows.

| Closed event | Exact upper-bound interval |
|---|---:|
| Bad support-33 output with (L(y)=0), all 26 words | \([2^{-40.980718082062},2^{-40.980718082061}]\) |
| Bad support-33 output with (L(y)\ne0), all packets in suffix | (0) |
| Bad suffix output, supports 35 through 39 | \([2^{-48.491771752756},2^{-48.491771752755}]\) |
| Sum of the closed rows | \([2^{-40.972830672896},2^{-40.972830672895}]\) |

The remaining numerical budget below (2^{-40}) is

\[
[2^{-41.027690825931},2^{-41.027690825930}].
\]

This number is only a budget for future rows. The open events have not been
assigned probability zero, so the displayed partial sum is not a bound for
the complete construction.

## Turnoff refutation evidence under PacketMul

The authenticated paused-candidate witness uses 33 active packets and has
two-lap weight 56. Fix one of its favorable labeled placements. Every original
packet value is nonzero, and every target value is nonzero. Exactly one of the
15 possible multipliers maps the original value to the target value at each
active coordinate. Independence therefore contributes the factor

\[
15^{-33}.
\]

Multiplying the authenticated old lower family by this factor gives a valid
PacketMul lower family in the interval

\[
[2^{-602.334330424314},2^{-602.334330424313}].
\]

Consequently, PacketMul does not give deterministic distance for every
realized setup. This lower family is far below (2^{-40}) and does not refute
the intended high-probability setup claim. The calculation is a lower bound on
one restricted family, not the total bad-setup probability.

## First open lemma

Let (r\in\{0,\ldots,22{,}652\}) be the first occupied node of a support-33
drive word. The first \(\lfloor r/3\rfloor\) complete autonomous blocks
contribute at least (25\lfloor r/3\rfloor) when (L(y)\ne0). Let
(W_{\mathrm{tail}}) be the output weight outside those complete blocks.

The next finite lemma is to prove, after summing over the 26 authenticated
support-33 words and all 22,653 values of (r), that

\[
L(y)\ne0
\quad\text{and}\quad
W_{\mathrm{tail}}
\le 188{,}765-25\lfloor r/3\rfloor
\]

has aggregate probability at most the remaining numerical budget above. This
is the first construction-level blocker because it is the lowest-support row
not covered by the terminal-zero proof or by the long autonomous prefix.

## Conclusion

Goal 03 closes the transfer audit and produces a rigorous partial ledger. It
also isolates an exact complement event. The full candidate remains open.
