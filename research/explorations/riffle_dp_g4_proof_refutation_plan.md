# Concurrent proof and refutation plan for Riffle DP g=4

## Objective

Fix the current sequential **Riffle DP g=4** candidate. Determine which of
the following claims is true:

1. The declared ensemble satisfies
   \[
   \mathbb E[Z_d]\le 2^{-40}.
   \]
2. An exact family of bad messages gives
   \[
   \mathbb E[Z_d]>2^{-40}.
   \]

Here, \(Z_d\) counts nonzero messages whose encoded binary word has weight
at most \(d=188{,}766\). The expectation is over every random choice in the
declared construction.

The two claims require different evidence. A proof needs an upper bound that
covers every nonzero message. A refutation needs one exact lower-bound family
whose contribution exceeds the target.

The project will pursue both outcomes concurrently. Every milestone must
produce one proof artifact and one refutation artifact.

## Frozen candidate

Before either track advances, bind the following data into one machine-readable
construction manifest:

- \(K=2^{20}\) message bits;
- \(B=16{,}384\) symbols in \(\mathbb F_{2^{64}}\);
- two parities
  \[
  p_0=\sum_i m_i,\qquad p_1=\sum_i\alpha_i m_i;
  \]
- the coefficient schedule \(\alpha_i=x^i\);
- the binary extended BCH \([128,64,22]\) encoder;
- \(B+2\) unpunctured outer BCH words;
- sequential division into four-bit packets;
- the uniform global permutation of all \(524{,}352\) packets;
- the exact recursive inner map and all its randomness;
- failure threshold \(d=188{,}766\).

The manifest must identify the implementation that computes the same binary
matrix. Any change creates a separately named candidate.

## Evidence classes

Every output must carry one of three labels.

| Label | Meaning |
|---|---|
| `EXACT` | Integer, rational, or exhaustive computation with an independent check |
| `RIGOROUS_BOUND` | A proved inequality with directed rounding where needed |
| `DIAGNOSTIC` | Sampling, optimization, or an incomplete search |

A diagnostic can prioritize work. It cannot close a proof or refutation gate.

## Shared technical core

Both tracks need a weighted enumerator for triples of local BCH words. For
fixed nonzero field coefficients \((a,b,c)\), define

\[
H_{a,b,c}(z)
:=
\sum_{t\in\mathbb F_{2^{64}}^*}
z^{s(E(at))+s(E(bt))+s(E(ct))},
\]

where \(E\) is the BCH encoder and \(s\) counts nonzero four-bit packets.

The four block-weight-three categories supply structured coefficient triples:

1. one data symbol and both parity symbols;
2. two data symbols and only \(p_0\) active;
3. two data symbols and only \(p_1\) active;
4. three data symbols and no active parity symbol.

Build one branch-and-bound engine for these linear maps. A node fixes part of
the 64-bit parameter \(t\). Rank calculations bound the packet support of all
completions.

The engine has two modes:

- **proof mode:** return a certified weighted upper bound for every unfinished
  branch;
- **refutation mode:** visit the most dangerous branches first and emit exact
  messages attaining small support.

The first version tracks total packet support. The second version tracks the
five-class packet-weight profile

\[
(a_0,a_1,a_2,a_3,a_4).
\]

The profile version is required for the complete inner-chain proof.

## Track P: proof

### P0. Freeze the ensemble

Produce the construction manifest and an implementation-equivalence test.
The test must compare the declared matrix with the transposed implementation
on fixed basis vectors and randomized inputs.

**Pass condition:** the manifest determines one probability law and one binary
matrix ensemble.

### P1. Complete the weight-three terminal bound

Use the weighted triple engine to bound every total support from 36 through
96. Aggregate the four coefficient categories with exact support
multiplicities.

Retain the completed results:

- total support 33: exact contribution \(2^{-109.9617}\);
- total support 34: upper bound \(2^{-70.3326}\);
- total support 35: upper bound \(2^{-67.1559}\).

**Pass condition:** the complete weight-three terminal contribution is below
an explicit budget smaller than \(2^{-40}\).

This gate only completes the terminal obstruction audit. It is not a proof
for all inner outputs.

### P2. Bound every inner outcome for block weight three

Replace terminal placement with the exact inner response. For each feasible
outer packet profile, bound

\[
\Pr[\operatorname{wt}(\mathsf{Inner}(C))\le d\mid\text{outer profile}].
\]

Reuse the certified \(g=4\) inner machinery only after matching its input
interface and randomness law to Riffle DP g=4. Recompute any bound that depends
on the previous outer layout.

**Pass condition:** the aggregate contribution of every block-weight-three
message is below its assigned first-moment budget.

### P3. Cover sparse outer weights

Handle outer block weights \(4\le w\le W\), where \(W\) is selected from the
computed tail. Use the exact MDS enumerator and coefficient-category bounds.

The proof should keep low-support BCH correlations for small \(w\). It must
not replace them with independent local words without justification.

**Pass condition:** every sparse shell has a verified upper bound, and their
sum fits the sparse-shell budget.

### P4. Cover dense outer weights

Use the favorable \(g=4\) entropy slope to treat \(w>W\). Candidate tools are:

- exponential moments of the nibble profile;
- MDS shortening and projection bounds;
- the existing outward-rounded packet-profile operator;
- a sparse/dense split with an overlap-safe union bound.

**Pass condition:** one rigorous inequality covers every remaining outer
word, including all parity-support patterns.

### P5. Assemble the first-moment ledger

Partition every nonzero message into disjoint or safely union-bounded cells.
For each cell, record:

- its population;
- its probability space;
- its inner failure bound;
- its first-moment contribution;
- the lemma or exact receipt that authenticates the row.

Sum all rows with outward rounding.

**Proof success:** an independent verifier confirms a total at most
\(2^{-40}\), and the manifest matches the implementation.

## Track R: refutation

### R0. Maintain a counterexample ledger

Record each tested family by outer weight, coefficient relation, BCH packet
profile, inner placement pattern, and exact or sampled contribution. Preserve
negative results so searches do not repeat.

### R1. Finish the terminal-family search

Search the same weight-three branches used by P1 in best-first order. Then
extend the search to outer weights four through eight.

Prioritize:

- repeated local BCH words;
- short additive relations among support-11 through support-16 words;
- small gaps among coefficients \(x^i\);
- cancellation of one or both parity symbols;
- supports that maximize value entropy per occupied packet.

For each candidate family, compute an exact disjoint count and its exact
terminal-placement probability.

**Refutation condition:** one authenticated family contributes more than
\(2^{-40}\).

### R2. Search nonterminal inner placements

The terminal suffix is only one bad placement. Build an exact or certified
search over sparse packet placements in the recursive inner chain.

The search state contains:

- the current inner state;
- remaining packet values;
- remaining node capacities;
- emitted weight;
- a lower bound on future emitted weight.

Use A* or branch-and-bound to minimize output weight. Count every permutation
orbit associated with an authenticated witness.

**Output:** either an exact lower-bound family or a list of worst profiles for
the proof track.

### R3. Search coefficient structure

The deterministic schedule \(\alpha_i=x^i\) may create exceptional relations.
Search for:

- low-degree additive relations among powers of \(x\);
- repeated multiplier classes in the four weight-three categories;
- small-support relations among four or more data symbols;
- coefficient gaps that preserve the BCH low-support set.

Every discovered relation must be converted into actual BCH words and an
exact probability contribution.

### R4. Use diagnostics only for discovery

Use importance sampling, local optimization, MILP, or randomized search to
find candidate messages and placements. Replay every candidate with exact
arithmetic.

Sampling cannot establish the absence of a counterexample. A sampled family
cannot refute the construction until its population and probability are
proved.

### R5. Escalate or refute

If a family crosses \(2^{-40}\), stop the proof track and publish the exact
counterexample receipt. If no family crosses the target, feed the worst
authenticated profiles back into Track P.

## Synchronization gates

| Gate | Proof deliverable | Refutation deliverable | Decision |
|---|---|---|---|
| G0: freeze | Construction manifest | Counterexample schema | Resolve any matrix ambiguity |
| G1: terminal \(w=3\) | Complete weighted upper bound | Best exact terminal families | Stop if refuted |
| G2: full inner \(w=3\) | Complete shell bound | Exact sparse-placement search | Stop if refuted |
| G3: sparse \(w\) | Bounds for \(4\le w\le W\) | Algebraic search through \(W\) | Stop if refuted |
| G4: dense \(w\) | Dense-shell inequality | Guided search of worst bound cells | Repair bound or stop |
| G5: final | Verified first-moment ledger | Final adversarial audit | Prove, refute, or report unresolved gap |

The project does not advance past a gate merely because the refutation search
finds nothing. The proof deliverable must satisfy its stated coverage test.

## Execution order for the next work tranche

Run the following tasks as two concurrent logical streams.

### Proof stream

1. Write the construction manifest.
2. Implement the weighted triple engine in proof mode.
3. Complete total supports 36 through 96 for the terminal event.
4. Compress the coefficient triples into exact multiplier classes.
5. Produce a machine-verifiable G1 upper-bound receipt.

### Refutation stream

1. Extend the exact support-11 relation search to supports 12 and 13.
2. Search all four weight-three categories in increasing total support.
3. Search outer weights four through six for parity cancellations.
4. Prototype the sparse nonterminal inner-placement search.
5. Publish the strongest exact lower-bound family, even if it remains below
   \(2^{-40}\).

The streams exchange data after each completed support layer. The proof stream
reports its loosest coefficient class. The refutation stream searches that
class first and returns witnesses or an authenticated empty region.

## Resource and change rules

- Never run two benchmarks simultaneously.
- Run only one large exhaustive enumeration at a time on this machine.
- Preserve exact receipts separately from sampled diagnostics.
- Use directed rounding for theorem-facing floating-point bounds.
- Verify every generated BCH word by re-encoding its message.
- Do not silently change coefficients, layout, puncturing, packet width, or
  inner randomness.
- Give every modified construction a distinct name and restart G0.

## Immediate recommendation

Start G0 and G1 together. Freeze the construction manifest while building the
dual-mode weighted triple engine. In parallel, extend the algebraic witness
search beyond the support-33 minimum layer. Do not begin the full profile
certificate until G1 either passes or produces a counterexample.

## Execution status on 2026-08-18

G0 and G1 are complete. The frozen manifest and implementation checks pass.
The complete terminal block-weight-three upper bound is
`2^-40.212996059540`, so G1 passes with 0.212996 bits of margin.

G2 is active. The exact local search gives:

```text
one-nibble autonomous crossing range       5871..5926 nodes
exact one-nibble turnoffs below distance              6
best authenticated turnoff lower family      2^-124.233993
```

The turnoff family does not refute the candidate. The next proof artifact
must cover multiple active packets under the actual fixed-state recurrence.
It must not import the earlier lane-permutation or state-permutation laws.

The exact autonomous-map audit now factors the 64-dimensional zero-input map
into irreducible components of degrees 1, 2, 4, 9, 10, 18, and 20. Every
individual component cycle averages at least 28 emitted bits per node.
Exhaustive combined-subspace checks through dimension 20 find no
cancellation-driven full-support average below 31.940371456. An independent
verifier replays about 3.8 million states and confirms the receipt.

This evidence closes a natural low-period refutation route. It does not close
G2 because finite prefixes and interactions among several packet impulses
remain uncovered. The next concurrent tranche is:

1. derive a finite-window lower bound for states reached by one and two packet
   impulses under the fixed recurrence;
2. search exact three-packet turnoffs by a meet-in-the-middle episode engine;
3. convert the resulting window bound and turnoff counts into a placement
   probability bound for each support-33 through support-38 outer profile.

The finite-window and three-packet tranche is now complete. A forward search
enumerates all 5,130,659,560 nonzero blocks of weight at most eight. An
independent reverse search reproduces its return statistics. The resulting
rigorous autonomous bound is

\[
\delta_{\mathrm{aut}}(L)
\geq 9L-8-8\left\lfloor\frac{L-1}{3}\right\rfloor.
\]

It forces weight above (d) at length 29,807 and gives 207,556 bits at full
length.

The full-chain distinct-node three-packet search tests 1,887,552,000
relations. It finds exactly 27 turnoffs, all by node six. The authenticated
support-33 lower family contributes `2^-134.276933547563`, so it does not
refute the candidate.

The active G2 task is now reset aggregation. Combine the autonomous bound with
exact counts for reset-compatible packet placements. Extend the refutation
search to packet collisions within a node and four-packet episodes.

The first reset aggregation is complete. It exactly counts the geometric
placement envelope implied by the autonomous lemma. At support 33, even the
zero-reset envelope is only `2^-4.475042948481`. Thus, this deterministic
route cannot supply the roughly 141 probability bits required per outer word.

The collision refutation search is also complete for one or two packets in
the initial and reset nodes. It finds 539 turnoffs, all by node eight. The
strongest authenticated support-33 family is `2^-123.048186206738` and does
not refute the candidate.

The proof stream must now exploit reached-state distribution under the random
packet permutation. The next candidate is a fixed-map transfer or
exponential-moment operator that retains nibble values and slots. The
refutation stream should search episodes with intermediate active nodes and
at least three packets in one node.

The reachable-prefix and component-mixing tranche is complete. Exact initial
orbit enumeration gives latest crossings 5,926, 5,943, and 5,945 for first-node
packet supports one, two, and three. The strengthened support-33 zero-reset
envelope is still only `2^-6.767496031893`.

An exact autonomous witness has weight 188,730 over 5,956 nodes. Its preimage
has nibble support 13, so it limits universal window proofs without giving an
authenticated sparse counterexample.

Exact Walsh spectra for all individual irreducible components and tested
combinations through dimension 20 show strong mixing of the 14 support-33
profiles. The worst reported quotient total variation is about
`2.452821e-19`; the primitive degree-20 quotient is below `9.34e-88`.

The quotient bounds now use directed rounding and are theorem-facing within
their stated projections. They do not cover characters spanning larger
component combinations. A full-character hill search finds no bias above the
known value 0.625, but it is not exhaustive.

The refutation stream now has an exact authenticated placement witness. A
support-33 profile packed into three consecutive nodes has weight 188,505 over
the final 5,940 nodes. This is 261 below (d). The same assignment remains bad
under 5,938 later terminal translations. Its single-family probability is
only `2^-555.961662409912`, so the witness refutes deterministic distance but
does not refute the random-permutation claim.

The active proof task has two parts. First, certify a bound for all full-state
characters and their aggregate Fourier mass. Second, transfer reached-state
mixing to the probability of low total output weight. The matching refutation
task is to enlarge the exact cluster into a counted basin or find broader
sparse episodes whose total contribution can approach `2^-40`.

The single-transposition cluster basin is now exact. It contains 919
assignments and 5,419,150 disjoint placements, but contributes only
`2^-544.542826136930` across the three compatible messages.

The broader suffix route is substantially stronger. Conditional on one
support-33 profile entering the final 5,940 nodes, 79,583 of 100,000 trials
fail the distance threshold. The three-message contribution estimate is
about `2^-80.06`.

An exact enumeration now covers 85,828 one-data outer words through support
39. A separate verifier rechecks all 631,767,040 candidate pairs. Sampling
100,000 placements per support gives an aggregate contribution estimate of
`2^-76.01295`. This does not refute the construction, but it identifies the
next refutation multiplier: support-39 words across the much larger pair-data
and three-data coefficient populations.

For the proof stream, the full-character problem can be written as a binary
observability code. A 12-node block emits 192 bits from a 64-bit character.
The known degree-one character has weight 36. Exact SAT and MILP attempts have
not yet certified that no block has weight at most 35. The next certificate
should exploit the recurrence or component decomposition rather than a generic
integer formulation.
