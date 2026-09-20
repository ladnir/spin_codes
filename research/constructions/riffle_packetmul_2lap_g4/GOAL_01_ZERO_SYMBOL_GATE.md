# Goal 01: establish the zero-symbol Fourier gate

**Status:** `COMPLETE_OPEN_MIXED_LEMMA` (2026-08-21)

The completed checkpoint is recorded in
`proof/GOAL_01_ZERO_SYMBOL_GATE_CHECKPOINT.md`. It ends with the explicit,
unproved mixed-character lemma \(q_\chi\le5/8\). That lemma implies
\(|\beta_\chi|\le3/5\) and closes the authenticated support-33 terminal-zero
row with approximately 0.981 bits of margin.

## Objective

Reduce the first PacketMul terminal-state question to one explicit,
quantitative lemma about mixed characters.

This goal must reconstruct the relevant probability space exactly. It must
certify the component calculations that are small enough to exhaust. It must
also search for evidence that refutes the proposed reduction or its required
numerical bound.

Completing this goal does not prove Riffle PacketMul-2Lap g=4. It establishes
whether the packet multiplier gives a tractable terminal-zero interface.

## Fixed experiment

Fix a nonzero outer packet and its source coordinate. Setup samples an
independent multiplier from \(\mathbb F_{16}^{\times}\). The unchanged packet
permutation assigns the randomized packet to one of

\[
M=524{,}352
\]

inner packet cells.

For a terminal-state character \(\chi\ne0\), define
\(a_\chi(p)\in\mathbb F_{16}\) by

\[
\langle\chi,z_y(p)\rangle
=\operatorname{Tr}(a_\chi(p)y)
\qquad(y\in\mathbb F_{16}).
\]

Define

\[
q_\chi
:=\frac{|\{p:a_\chi(p)=0\}|}{M},
\qquad
\beta_\chi:=\frac{16q_\chi-1}{15}.
\]

The proof must derive this interface from the frozen construction. It must
not assume the formula merely because it holds for a diagnostic script.

## Work packages

### 1. Algebraic reconstruction

Reconstruct \(z_y(p)\) from the frozen BCH, accumulator, packet-slot, and node
maps. Verify linearity in \(y\). Derive \(a_\chi(p)\) from the four basis
packet values.

Produce a primary implementation and an independent verifier. Each artifact
must identify Riffle PacketMul-2Lap g=4 and the active manifest.

### 2. Collision and Parseval calculation

For

\[
\Phi(p,y):=z_y(p),
\qquad
(p,y)\in[M]\times\mathbb F_{16}^{\times},
\]

compute the exact collision multiplicities of \(\Phi\). Use those
multiplicities to state the exact Parseval sum for the normalized
PacketMul character coefficients.

This package must report collisions rather than assume cross-value
injectivity.

### 3. Exact component spectra

For each pure recurrence component of degree

\[
1,\ 2,\ 4,\ 9,\ 10,\ 18,\ 20,
\]

exhaust every nonzero character. Record the maximum zero-symbol fraction,
the maximum \(|\beta_\chi|\), and an exact maximizing character.

Independently replay every reported maximizer.

### 4. Mixed-character refutation search

Search the full 64-bit character space with:

- every pure-component maximizer;
- XOR combinations of component maximizers;
- the known mixed-component witnesses from the paused candidate;
- constrained searches that cannot descend into only the degree-one and
  degree-two subspace; and
- reproducible pseudorandom restarts.

Every reported statistic must be exact for its character. Absence of a found
counterexample remains diagnostic.

### 5. Terminal-zero gate

Apply iid packet-cell sampling first. Then condition on distinct cells to
recover the ordered packet-permutation law.

Using the exact collision term, determine an explicit sufficient cap

\[
|\beta_\chi|\le B_{\mathrm{req}}
\qquad(\chi\ne0)
\]

for the authenticated support-33 terminal-zero row. Evaluate the cap against
the exact component results and the mixed-character search.

## Completion criteria

The goal is complete when all of the following hold:

1. The coefficient map and its probability space have primary and independent
   reconstructions.
2. The cross-value collision histogram and Parseval identity are exact.
3. Every pure irreducible component has an exact zero-symbol certificate.
4. The mixed-character search has a reproducible receipt and exact witness
   replay.
5. The terminal-zero calculation states \(B_{\mathrm{req}}\) numerically and
   identifies its security margin.
6. The checkpoint ends with either:
   - one explicit mixed-character lemma whose proof would close this
     terminal-zero row; or
   - an exact character that refutes the required cap or the proposed
     interface.

## Excluded work

This goal does not:

- prove the mixed-character lemma;
- close every outer-word or placement stratum;
- modify the construction;
- benchmark GF(16) multiplication;
- optimize the integrated encoder; or
- revisit the paused DP-2Lap block-moment route.

These exclusions keep the checkpoint medium-sized and make a negative result
useful.
