# Goal 03: joint BCH spectrum under the alpha schedule

## Question

Does the deterministic schedule \(\alpha_i=\gamma^i\) make the third BCH
weight in a one-data-symbol outer word behave independently of the first two?

## Status

**Diagnostic complete; independence is false near the start of the schedule.**

A uniform random coefficient index behaves independently at ordinary BCH
weights. However, the first coefficients preserve many authenticated
low-weight words. In particular, \(\alpha_0=1\), so the three BCH blocks are
identical at data position zero.

The exact authenticated tail contains 3,622 pairs \((x,i)\) for which both
\(B(x)\) and \(B(\gamma^i x)\) have binary weight 22. This is a lower bound on
the number of current one-data-symbol outer words with profile \((22,22,22)\).

A shifted geometric schedule \(\alpha_i=\gamma^{64+i}\) preserves the cheap
recurrence. The separately registered candidate has now been checked against
all 243,840 minimum words: its shifted window contains no weight-22 pair.
This result belongs to the distinct shifted candidate and does not change the
current construction.

See `proof/GOAL_03_ALPHA_JOINT_SPECTRUM.md` for the analysis.
