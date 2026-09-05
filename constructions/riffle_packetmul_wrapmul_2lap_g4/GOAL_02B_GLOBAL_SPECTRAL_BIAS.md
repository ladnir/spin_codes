# Goal 02B: global spectral-bias closure

## Objective

Prove or refute a global response-bias bound for **Riffle
PacketMul-WrapMul-2Lap g=4**. A proof closes the remaining early-start part of
Goal 02. A refutation must identify which level of the argument fails.

This is a large proof checkpoint. It replaces further enumeration by
irreducible-component dimension with one theorem over the full response code.
It does not change the construction.

## Response code

Let

\[
J:\mathbb F_2^{64}\longrightarrow\mathbb F_2^n,
\qquad n=2{,}097{,}408,
\]

be the authenticated autonomous response map from Goal 02. For
\(u\in\mathbb F_2^{64}\), define its response bias by

\[
B(u):=n-2\operatorname{wt}(J(u)).
\]

Let \(E\) contain the four nonzero states supported only on the degree-one or
degree-two irreducible component. The exact low-component census gives

\[
|E|=4,
\]

with absolute biases at most \(196{,}632\) for the degree-one state and
\(262{,}176\) for each of the three degree-two states.

## Target theorem

The proof target is

\[
\boxed{
  |B(u)|\le106{,}802
  \quad\text{for every }u\in
  \mathbb F_2^{64}\setminus(\{0\}\cup E).
}
\]

The exact census already proves the target for every nonexceptional state
whose irreducible-component support has dimension at most 22. Its largest
authenticated absolute bias is \(104{,}882\). Thus the new mathematical
work concerns only support dimension at least 23.

## Why the theorem suffices

Write the response coordinates as

\[
J(u)=(\langle a_1,u\rangle,\ldots,\langle a_n,u\rangle).
\]

For \(x\in\mathbb F_2^{64}\), let \(r(x)\) count unordered pairs of distinct
coordinate columns whose XOR equals \(x\). Define

\[
\lambda_u:=\sum_x r(x)(-1)^{\langle u,x\rangle}
           =\frac{B(u)^2-n}{2}.
\]

The pair-sum triangle energy is

\[
T(r):=\sum_{x,y}r(x)r(y)r(x+y)
     =2^{-64}\sum_u\lambda_u^3.
\]

The exact pair census gives

\[
\sum_x r(x)=2{,}199{,}559{,}110{,}528
\]

and

\[
\sum_x r(x)^2=2{,}199{,}581{,}502{,}666.
\]

After the trivial character and the four states in \(E\) are charged
separately, the target theorem implies

\[
T(r)\le12{,}543{,}557{,}200{,}693{,}295{,}704{,}650.
\]

Equivalently, it implies the sufficient sixth-moment bound

\[
A_6^\perp\le139{,}270{,}335{,}354{,}994{,}389{,}965.
\]

Together with the authenticated value

\[
A_4^\perp=3{,}732{,}023,
\]

this gives

\[
A_{\le377{,}530}\le316{,}605.
\]

Goal 02 then bounds every early bad-state list by 316,606. The completed late
charge therefore closes the authenticated support-33 nonzero-terminal row.

## Proof track

1. Authenticate the Fourier identity, the four exceptional states, and the
   exact threshold 106,802 in a replayable receipt.
2. Express \(B(u)\) as an orbit sum over the seven irreducible components of
   the autonomous map.
3. Prove cancellation for every component support of dimension at least 23.
   The proof may use a block decomposition, a correlation inequality, or an
   exact finite certificate with subexponential coverage.
4. Avoid a frontier sweep whose work scales as \(2^d\) for support dimension
   \(d\). Any finite certificate must cover many supports through one stated
   invariant.
5. Replay the global bound independently from the construction manifest and
   the authenticated recurrence.
6. Substitute the result into the triangle-energy and sixth-moment identities.
7. Update the disjoint support-33 ledger.

## Refutation track

1. Search all nonexceptional component classes for
   \(|B(u)|\ge106{,}804\). Bias is even, so 106,804 is the first violation.
2. Confirm each candidate by direct evaluation of all 32,772 response nodes.
3. If the bias theorem fails, determine whether a larger exceptional set or a
   direct triangle-energy bound still satisfies the exact energy cap.
4. Search separately for response words of weight at most 377,530 and affine
   bad-state lists larger than 316,606.
5. Treat a construction-level refutation as complete only after proving that
   an over-cap affine center is reachable and has probability-relevant setup
   mass.

A character above the bias threshold refutes this sufficient theorem only.
It does not by itself refute Goal 02 or the construction.

## Milestones

### A. Spectral interface

Produce primary and independent receipts for the exact Fourier formulas,
energy budget, exceptional-state charge, and threshold 106,802.

**Status:** complete. Both receipts pass with current dependency hashes.

### B. High-support diagnostic

Run a bounded search for violations outside \(E\). Record the search domain,
restart law, explored states, best witnesses, and exact replays. This result
is diagnostic unless the search exhausts a declared domain.

**Status:** complete as a bounded diagnostic. All 92 high-support masks were
searched. An independent implementation replayed all 184 retained extrema.
The maximum found absolute bias is 6,482.

### C. Global cancellation lemma

Prove one uniform inequality for all component supports of dimension at least
23. The lemma must expose the structural source of cancellation. A list of
independently solved supports does not satisfy this milestone.

**Status:** open. A sufficient scalar subgoal is
\(|S_j(u)|\le1668\) for all 64 output coordinates. The best retained global
bias witness has scalar maximum 584, but no applicable incomplete-character-
sum theorem currently proves the uniform scalar bound.

### D. Row closure or precise failure

Complete the independent replay and ledger update, or record the smallest
authenticated obstruction that prevents the global lemma.

## Completion criteria

Goal 02B ends with one of the following outcomes.

1. The target theorem is proved and independently replayed, and Goal 02's
   complete support-33 nonzero-terminal row is closed.
2. An exact energy computation closes Goal 02 even though the pointwise bias
   theorem fails.
3. A probability-relevant reachable affine list refutes the support-33 row.
4. A reproducible counterexample or finite solver certificate isolates a
   precise obstruction to both the pointwise and aggregate spectral routes.

The goal does not prove higher packet supports or the complete construction.
It authorizes no performance benchmark and no construction modification.
