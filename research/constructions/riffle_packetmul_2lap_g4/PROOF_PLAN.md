# Riffle PacketMul-2Lap g=4: proof plan

**Checkpoint status:** Goals 01 through 05 complete. Goal 04 proves a sharper
four-node autonomous certificate and moves the sufficient zero-prefix cutoff
to node 20,976. Goal 05 leaves the exact four-node distance open in [36,42]
after two bounded exact-solver attempts. See
`proof/GOAL_05_FOUR_NODE_DISTANCE_REPORT.md`.

The next bounded goal should address the immediate support-33,
nonzero-terminal boundary band whose first occupied node is 20,972 through
20,975. In that band, a bad output has tail weight at most 17.

## Proof interface

Fix a terminal-state character \(\chi\ne0\). For a packet cell \(p\), define
\(a_\chi(p)\in\mathbb F_{16}\) by

\[
\langle\chi,z_y(p)\rangle
=\operatorname{Tr}(a_\chi(p)y)
\qquad(y\in\mathbb F_{16}).
\]

Let \(q_\chi\) be the fraction of packet cells for which
\(a_\chi(p)=0\). Averaging one independently randomized nonzero packet gives

\[
\beta_\chi
=q_\chi-\frac{1-q_\chi}{15}
=\frac{16q_\chi-1}{15}.
\]

Thus the fifteen fixed-value character sums of the paused candidate collapse
to one zero-symbol statistic.

The known degree-one character has \(q_\chi=1/16\), and hence
\(\beta_\chi=0\). The primary calculation and independent replay now certify
this coefficient histogram.

## Proof tasks

1. **Complete.** Reconstruct the exact map from \((p,y)\) to terminal contribution
   \(z_y(p)\).
2. **Complete.** Count cross-value collisions in that map and derive its Parseval identity.
3. **Complete.** Compute exact \(q_\chi\) maxima on each irreducible component.
4. **Complete as diagnostic evidence.** Search the full character space for a large \(q_\chi\), with structured
   mixed-component seeds.
5. **Complete for the support-33 row.** Determine a sufficient uniform cap on \(|\beta_\chi|\) for each authenticated
   packet support.
6. **Complete conditionally.** Insert the cap and the exact collision term into the terminal-zero ledger.
7. Audit which placement-only and autonomous second-lap results transfer from
   Riffle DP-2Lap g=4.

## Refutation tasks

1. Search for characters with an unusually large zero-symbol fraction.
2. Replay known double turnoffs after independent packet multiplication.
3. Search for multiplier-invariant low-weight families.
4. Separate a bad fixed multiplier schedule from a bad probability over the
   sampled schedule.

## Decision rule

Continue this candidate if the zero-symbol statistic admits a substantial
uniform cap and closes the terminal-zero ledger with margin. Reassess the
candidate if mixed-component zero sets reproduce the prior support-case
explosion.

## Decision

Continue this candidate. The exact cap \(q_\chi\le5/8\) closes the
authenticated support-33 terminal-zero row. The proof avoids component-support
case enumeration.
