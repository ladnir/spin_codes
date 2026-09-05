# Goal 02: support-33 state-response list decoding

**Status:** complete under the precise-obstruction outcome. The support-33
nonzero-terminal row is not closed.

## Objective

Close or refute the authenticated support-33 event with a nonzero first-lap
terminal state. The proof must control the number of wrapped states that yield
a bad two-lap output for every relevant driven offset.

This is a large checkpoint. It combines a general tail-code theorem, exact
placement accounting, low-weight response-spectrum certification, and a
search for reachable over-cap lists.

## Active objects

Let \(N=32{,}772\) be the number of nodes and let

\[
D_{\mathrm{bad}}=188{,}765.
\]

For a fixed driven word \(y\), define

\[
Y_y(z):=F(F(y))\mathbin\oplus J(z)
\qquad
(z\in\mathbb F_2^{64}).
\]

The nonzero bad-state list is

\[
\mathcal B_y
:=\{z\ne0:\operatorname{wt}(Y_y(z))\le D_{\mathrm{bad}}\}.
\]

Goal 01 proves that WrapMul samples \(z\) uniformly from the \(2^{64}-1\)
nonzero states whenever the first-lap terminal state is nonzero.

## Part I: tail response codes

Let \(r\) be the first occupied node and set

\[
q:=\left\lfloor\frac r4\right\rfloor.
\]

The first \(4q\) nodes contain \(q\) complete autonomous blocks and contribute
at least \(36q\). Remove those nodes and define

\[
n_r:=64(N-4q),
\qquad
t_r:=D_{\mathrm{bad}}-36q.
\]

Two outputs with the same drive and different wrapped states differ by an
autonomous trajectory. Partitioning the remaining nodes into four-node blocks
should give the tail response code the minimum-distance bound

\[
d_r\ge36(8193-q).
\]

Goal 02 must prove this statement from the exact recurrence. The proof must
not assume that the common drive is zero; it must cancel the common drive
before applying the autonomous certificate.

For a binary code of length \(n\), minimum distance \(d\), and a list inside a
ball of radius \(t\), the Johnson counting argument gives

\[
L\le
\left\lfloor
\frac{d}{d-2t(1-t/n)}
\right\rfloor
\]

when the denominator is positive.

Goal 02 must authenticate the following cutoffs.

- For every \(r\ge7060\), the list has at most 13,749 states.
- At \(r=7060\), the parameters are
  \[
  (n_r,d_r,t_r)=(1{,}645{,}568,231{,}408,125{,}225).
  \]
- For every \(r\ge9176\), the inequality \(2t_r<d_r\) makes the list unique.

The new charged band is \(7060\le r\le20975\). The parent certificate already
makes \(r\ge20976\) empty. Its exact placement probability is

\[
\Pr[7060\le r\le20975]
=
\frac{
\binom{524{,}352-16\cdot7060}{33}
-\binom{524{,}352-16\cdot20976}{33}
}
{\binom{524{,}352}{33}}.
\]

Using the maximum list size 13,749, the aggregate charge over 26 words should
be approximately \(2^{-57.1036423524}\). Goal 02 must reproduce this value
with exact rational arithmetic.

## Part II: early response spectrum

The remaining region is \(r\le7059\). Define the low-weight response-state
count

\[
A_{\le w}
:=
\#\{u\in\mathbb F_2^{64}\setminus\{0\}:
\operatorname{wt}(J(u))\le w\}.
\]

Fix a nonempty bad-state list and a member \(z_0\). For every other
\(z\in\mathcal B_y\),

\[
J(z+z_0)=Y_y(z)\mathbin\oplus Y_y(z_0),
\]

so

\[
\operatorname{wt}(J(z+z_0))\le2D_{\mathrm{bad}}=377{,}530.
\]

The map \(z\mapsto z+z_0\) is injective. Therefore,

\[
|\mathcal B_y|\le1+A_{\le377{,}530}.
\]

After charging the \(r\ge7060\) region, exact arithmetic should permit at most
316,606 bad states per fixed early drive. The sufficient spectrum target is

\[
A_{\le377{,}530}\le316{,}605.
\]

This target is independent of the driven offset. If proved, it closes every
early support-33 drive simultaneously.

## Proof track

1. Reconstruct the state-response map \(J\) from the authenticated recurrence.
2. Prove the tail difference-code lemma.
3. Certify the Johnson and unique-decoding cutoffs independently.
4. Charge the late placement region with exact rational arithmetic.
5. Bound or determine \(A_{\le377{,}530}\) using orbit decomposition,
   meet-in-the-middle enumeration, or another exact method.
6. Replay every low-weight state or compressed certificate independently.
7. Update the disjoint support-33 ledger without assuming independence across
   outer words.

The implementation must preserve explicit fixed-width operations and
predictable memory access in any enumeration hot path.

## Refutation track

1. Search the response code for unexpectedly many states below weight 377,530.
2. Search affine centers for lists larger than the active early cap.
3. Test whether an over-cap center equals \(F(F(y))\) for a reachable
   support-33 drive.
4. If it is reachable, compute a rigorous lower bound on its setup
   probability.

An over-cap spectrum count refutes only the sufficient spectrum route. An
over-cap affine list refutes only the uniform list lemma. A construction-level
refutation additionally requires enough probability under the setup
experiment.

## Completion criteria

Goal 02 ends with one of the following authenticated outcomes.

1. The complete support-33 nonzero-terminal row is charged within the remaining
   ledger budget.
2. A reachable affine list gives a probability-relevant counterexample.
3. The late region is closed, and an exact spectrum or solver certificate
   isolates a precise coding-theoretic obstruction in the early region.

Every exact computation requires a primary certificate and an independent
replay. A bounded solver result must record its exact predicate and resource
limit.

## Current status

Part I is complete. The exact late-band charge leaves an early-drive list cap
of 316,606.

For Part II, an exhaustive primary certificate and an independent replay
exclude every autonomous state whose exact irreducible-component support has
total dimension at most 22. These 35 supports contain 15,650,667 nonzero
states, and their minimum response weight is 917,616.

The remaining spectrum target has also been reduced to the sparse dual
enumerators of the full response code. Let (A_j^\perp) count its weight-(j)
dual words. It now suffices to prove

\[
2^{64}\bigl((360n-960)A_4^\perp+720A_6^\perp\bigr)
\le
1849740699995876325368244902659070749106176,
\qquad n=2097408.
\]

The coordinate columns are all nonzero and distinct. They nevertheless have
an authenticated weight-4 relation, so six-wise independence is false. The
primary and independent normalized-pair orbit censuses determine the complete
weight-4 dual count as

\[
A_4^\perp=3732023.
\]

After charging this exact term, it suffices to prove

\[
A_6^\perp\le139270335354994389965.
\]

The active proof task is to bound the weight-6 relations. An alternative is
to prove that the eight-node autonomous response code has distance at least
93; the current exact generic-solver encoding remains `UNKNOWN`, as do the
earlier C12 predicates. See
`proof/GOAL_02_TAIL_AND_SPECTRUM_PROGRESS.md`.

The large continuation checkpoint is
`GOAL_02B_GLOBAL_SPECTRAL_BIAS.md`. It reduces the weight-6 task to the
pointwise target

\[
|2{,}097{,}408-2\operatorname{wt}(J(u))|\le106{,}802
\]

outside four explicitly charged short-component states. The exact
low-component census already proves this target through component-support
dimension 22, except for those four states. The remaining theorem begins at
dimension 23 and must use global orbit structure rather than a frontier
enumeration.

## Excluded claims

Goal 02 does not close higher packet supports or prove the full construction.
It does not benchmark the integrated candidate.

## Outcome

Part I is proved and independently replayed. For every \(r\ge7060\), the bad
list has size at most 13,749. It is unique from \(r=9176\), and the exact
late-stratum charge has logarithm in

\[
[-57.103642352420,-57.103642352419].
\]

Part II ends at a precise authenticated coding-theoretic obstruction. Exact
enumeration excludes every state of component-support dimension at most 22.
The exact sparse-dual census gives

\[
A_4^\perp=3732023,
\]

and reduces closure to

\[
A_6^\perp\le139270335354994389965.
\]

Equivalently, the actual pair-sum function must satisfy

\[
T(r)\le12543557200693295704650.
\]

A sufficient pointwise form is

\[
|2{,}097{,}408-2\operatorname{wt}(J(u))|\le106802
\]

outside four charged degree-one and degree-two states. The remaining scope is
exactly the 92 component supports of dimension at least 23.

The authenticated pair histogram cannot prove the energy target. An exact
countermodel has the same histogram, zero value at the origin, the same
\(L_1\) and \(L_2\) data, and maximum multiplicity five, but its triangle
energy exceeds the target by a factor of 385.512913. Thus any continuation
must use the actual orbit locations or another response-code invariant.

The concurrent bounded search performs 149,744 exact high-support response
evaluations. Independent reconstruction replays all 184 retained extrema.
The largest found absolute bias is 6,482, and the smallest replayed
high-support weight is 1,045,463. No explored difference can seed a two-state
bad affine list, so no reachable over-cap center or setup-mass candidate was
found. This search is diagnostic and does not close the row.

The primary obstruction audit and its independent manifest replay are:

- `receipts/goal02_coding_obstruction_audit.json`;
- `receipts/goal02_coding_obstruction_independent.json`.

This outcome satisfies completion criterion 3. It is not a proof of the
support-33 row, a refutation of the construction, or a full-construction
claim.
