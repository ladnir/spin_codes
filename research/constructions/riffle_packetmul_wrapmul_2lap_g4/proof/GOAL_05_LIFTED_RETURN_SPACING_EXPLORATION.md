# Goal 05 exploration: lifted return spacing

## Result

The paired target proposed in this exploration is false. The later exact
search found a nonzero state with paired nine-node weights 47 and 50. Their
sum is 97, below the required bound of 106. See
`GOAL_05_PAIRED_RETURN_REFUTATION.md` for the counterexample and independent
verification.

The remainder of this note records the reduction and the evidence that led
to the exact refutation.

Orbit continuity materially improves the proof target. The ten-node
minimum-distance condition treats every block as an independent worst-case
restart. The actual lifted state restarts only at an occupied node.

Pairing adjacent nine-node blocks gives a cleaner sufficient theorem. Define

\[
W_9(s):=\sum_{t=0}^{8}\operatorname{wt}(q_t(s)),
\]

where \(s\in\mathbb F_2^{128}\) is the entering lifted state and
\(q_t(s)\) follows the zero-input lifted recurrence. Let \(R\) denote one
lifted state transition. It suffices to prove

\[
W_9(s)+W_9(R^9s)\ge106
\qquad
\text{for every }s\ne0.
\tag{RS}
\]

Equivalently, the 18-node lifted observation code must have minimum distance
at least 106. Condition (RS) is a return-cost statement: a cheap nine-node
transient must be followed by an expensive nine-node segment.

The initial exploration found no counterexample. Its best 18-node witness had
weight 319. Exact enumeration excluded every primary-component support of
lifted dimension at most 22. Native-XOR solver runs did not decide whether a
word of weight at most 105 exists. The later low-node decomposition found the
weight-97 counterexample.

## Why the paired target suffices

After the two boundary nodes, an authenticated support-33 word leaves at
least 32,737 zero-input nodes in at most 34 gaps. Partition each gap into
18-node blocks. The number of complete blocks is at least

\[
B_{18}
:=
\left\lceil
\frac{32737-34\cdot17}{18}
\right\rceil
=1787.
\]

First split off trajectories whose lifted state reaches zero. For every fixed
drive and time, at most one wrapped state causes that event. The existing
return-state union bound charges this event below \(2^{-44}\).

On the complement, every complete 18-node block begins in a nonzero lifted
state. Condition (RS) would therefore give output weight at least

\[
1787\cdot106=189422>188765.
\]

Thus (RS), together with the existing exact-return charge, closes the
support-33 nonzero-terminal row.

## Evidence for recovery after a low transient

The exact nine-node witness from Goal 04 has node weights

\[
(0,2,0,10,0,5,0,26,1).
\]

Its first nine-node weight is 44. The next overlapping nine-node window,
shifted by one node, has weight 76. Its next nonoverlapping nine-node block
has weight 275. The complete 18-node prefix has weight 319.

A deterministic search minimized the larger of two nine-node window weights
at several offsets. It found no pair in which both windows have weight at
most 52. At offset nine, the best searched pair had weights 243 and 242.
These are exact trajectories, but the search is not exhaustive.

The result supports the proposed mechanism. Low windows appear as transients
and recover immediately. The independent-block proof loses this recovery
cost by restarting the state at every block boundary.

## Branch-number formulation

The first nine outputs determine the 128-bit lifted state. Therefore,
\(\mathcal C_{18}\) is the graph of an invertible shift map on
\(\mathcal C_9\):

\[
\mathcal C_{18}
=
\{(c,S(c)):c\in\mathcal C_9\}.
\]

Condition (RS) is the Hamming branch-number bound

\[
\operatorname{wt}(c)+\operatorname{wt}(S(c))\ge106
\qquad(c\ne0).
\]

This formulation exposes the useful structure. A proof need not determine
the complete distance spectrum of a generic binary \([1152,128]\) code. It
may instead classify the low-weight part of \(\mathcal C_9\) and prove that
the shift \(S\) sends every such word to a sufficiently heavy word.

## Comparison with the ten-node target

The earlier target \(d(\mathcal C_{10})\ge59\) had only 17 bits of diagnostic
margin: the best exact witness weighed 76. The paired target (RS) has 213 bits
of diagnostic margin: the best exact witness weighed 319.

The larger margin does not make the theorem automatic. It does make a
low-list classification or branch-number argument more plausible than a
direct ten-node minimum-distance proof.

## Exact primary-component exclusion

Let (f_1,\ldots,f_7) be the square-free irreducible factors of the
64-dimensional autonomous map (T). Their degrees are

\[
(1,2,4,9,10,18,20).
\]

The lifted map (R) contains two copies of (T) on its diagonal. Direct
kernel computation gives the lifted primary spaces

\[
V_i:=\ker(f_i(R)^2)
\]

with dimensions

\[
(2,4,8,18,20,36,40).
\]

The seven spaces are invariant under (R) and form a direct sum of the full
128-dimensional lifted state space.

An exact Gray-code enumeration covered every state whose exact primary
support has total lifted dimension at most 22. This covers 12 support types
and 9,191,400 nonzero states. No enumerated state has nine-node weight at
most 52. The smallest enumerated nine-node weight is 221. The smallest
enumerated 18-node weight is 483.

The known state with nine-node weight 44 has a nonzero coordinate in every
one of the seven primary spaces. Thus the low transient is not confined to a
small exceptional orbit. It is produced by cancellation among all primary
components. A proof based only on listing low-dimensional components cannot
close (RS); it must control cancellation between components.

## Anchored form of a counterexample

There is a second exact reduction. Suppose that a nonzero state violates
(RS). Its 18 node weights are nonnegative integers with sum at most 105.
Therefore at least one node has weight at most

\[
\left\lfloor\frac{105}{18}\right\rfloor=5.
\]

For (0\le j<18), define

\[
\mathcal A_j
:=
\left\{s\ne0:
  \sum_{t=0}^{17}\operatorname{wt}(q_t(s))\le105,
  \ \operatorname{wt}(q_j(s))\le5
\right\}.
\]

Every counterexample lies in

\[
\bigcup_{j=0}^{17}\mathcal A_j.
\]

This replaces one unstructured 128-dimensional distance decision with 18
anchored decisions. For a fixed anchor value (q_j=u), the compatible states
form a 64-dimensional affine space. The zero-anchor slice is a
64-dimensional linear space and contains the observed low transient.

The anchored reduction is exact, but it is not yet a proof. Native-XOR runs
on the unrestricted problem, the two low-half cases, the node-zero
weight-at-most-five case, and the exact zero-anchor slice all timed out
without a counterexample. A timeout supplies no lower bound.

## Superseded next goal

The original low-list target remains valid:

\[
\mathcal L_{52}
:=
\{c\in\mathcal C_9:\operatorname{wt}(c)\le52\}.
\]

The target is not necessarily a complete enumeration. It is enough to
certify

\[
\operatorname{wt}(S(c))\ge106-\operatorname{wt}(c)
\qquad
\text{for every }c\in\mathcal L_{52}.
\]

The anchored decomposition produced an exact violating state in the
zero-anchor slice. No further classification of \(\mathcal L_{52}\) is needed
to decide (RS). The next construction-level question is whether a longer
return window or a global permutation argument can replace (RS).

## Scope

This note derives a sufficient return-spacing target, excludes all lifted
primary supports of dimension at most 22, and reduces every possible
counterexample to 18 anchored cases. The separate refutation note disproves
the target. Neither note closes the support-33 row or decides the full
construction.
