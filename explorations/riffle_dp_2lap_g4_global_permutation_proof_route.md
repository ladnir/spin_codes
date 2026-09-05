# Riffle DP-2Lap g=4: global-permutation proof route

## Terminal-state event

Fix an authenticated outer word with packet values
\(v_1,\ldots,v_h\in\mathbb F_2^4\setminus\{0\}\). Let \(M=524{,}352\)
be the number of packet cells. For a cell \(p\), let

\[
z_v(p)\in\mathbb F_2^{64}
\]

denote the contribution of value \(v\) in cell \(p\) to the first-lap
terminal state. Linearity gives

\[
L=\bigoplus_{i=1}^h z_{v_i}(P_i),
\]

where \((P_1,\ldots,P_h)\) is a uniform ordered injection into the \(M\)
cells. For the final-10,119-node stratum, the three-node autonomous
certificate reduces the bad event to \(L=0\).

## Correct reduction from the packet permutation

First sample \(Q_1,\ldots,Q_h\) independently and uniformly from the packet
cells. Let \(D\) be the event that these cells are distinct. Conditioned on
\(D\), the tuple \((Q_1,\ldots,Q_h)\) has the same law as
\((P_1,\ldots,P_h)\). Therefore, for every event \(E\),

\[
\Pr_{\mathrm{perm}}[E]
=\Pr_{\mathrm{iid}}[E\mid D]
\le \frac{\Pr_{\mathrm{iid}}[E]}{\Pr_{\mathrm{iid}}[D]},
\qquad
\Pr_{\mathrm{iid}}[D]=\frac{(M)_h}{M^h}.
\]

An earlier receipt multiplied upper bounds on adaptive conditional character
biases. That multiplication is invalid. For a balanced sign function, two
draws without replacement have expectation \(-1/(M-1)\), not a value bounded
by \(1/(M-1)^2\). The artifact
`riffle_dp_g4_g2_support33_component_mixing_probe.json` is invalidated.

## Fourier and Parseval gate

For \(\chi\in\mathbb F_2^{64}\), define the normalized character sum

\[
b_v(\chi):=\frac1M\sum_p(-1)^{\langle\chi,z_v(p)\rangle}.
\]

Under iid cell sampling, Fourier inversion gives

\[
\Pr[L=0]
=2^{-64}\sum_\chi\prod_{v=1}^{15}b_v(\chi)^{c_v},
\]

where \(c_v\) is the multiplicity of packet value \(v\). The map
\(p\mapsto z_v(p)\) is injective for every nonzero \(v\). Hence Parseval gives

\[
\sum_\chi b_v(\chi)^2=\frac{2^{64}}M.
\]

Suppose \(|b_v(\chi)|\le B_v\) for every nonzero \(\chi\). If \(c_u\ge2\),
the two copies of value \(u\) can pay the exact Parseval sum. The remaining
copies pay their character caps. Thus,

\[
\Pr_{\mathrm{perm}}[L=0]
\le
\frac{2^{-64}}{\Pr[D]}
\left(
1+\left(\frac{2^{64}}M-1\right)
B_u^{c_u-2}\prod_{v\ne u}B_v^{c_v}
\right).
\]

The best diagnostic character caps give an aggregate value near
\(2^{-59.2956}\) across the 26 authenticated support-33 words. These values
are proof targets because the full-character hill search is not exhaustive.

Exact maxima are unnecessary. The coarser caps

\[
B_{15}\le\frac58,
\qquad
B_v\le\frac12\quad(v\ne15)
\]

already give an aggregate bound in

\[
[2^{-44.464429949134},2^{-44.464429949133}].
\]

## Refuted finite local certificate

The node count factors as \(32{,}772=12\cdot2{,}731\). For value \(v\), let
\(C_v\) be the 64-dimensional binary code that observes one character state
through 16 packet slots and 12 consecutive node exponents. Each codeword has
length 192.

The following local statements would imply the coarse character caps:

\[
\min\{\operatorname{wt}(c),192-\operatorname{wt}(c)\}\ge36
\quad(c\in C_{15}\setminus\{0\}),
\]

and

\[
\min\{\operatorname{wt}(c),192-\operatorname{wt}(c)\}\ge48
\quad(c\in C_v\setminus\{0\},\ v\ne15).
\]

The transpose recurrence is invertible. Therefore, if both statements held,
every 12-node block would start from a nonzero character state and summing the
local bounds over 2,731 blocks would give the required global character caps.
The known degree-one character `0xa685aac60e5acfb3` has value-15 block weight
36, so the first proposed local bound is tight.

An exact partial certificate exhausts every character whose irreducible-
component support has total dimension at most 22. It performs 274,291,320
character/value evaluations, including overlaps between the maximal
subspaces. Every covered character satisfies the required two-sided distance.
The individual degree-18 and degree-20 components each have minimum two-sided
distance 58 across all values.

The minimal uncovered supports are

\[
\begin{gathered}
(1,2,20),\ (1,4,18),\ (4,9,10),\ (2,4,18),\ (4,20),\\
(9,18),\ (10,18),\ (9,20),\ (10,20),\ (18,20).
\end{gathered}
\]

The last support has dimension 38. Its 192 coordinates contain five disjoint
information sets. A Brouwer--Zimmermann search for weight at most 47 therefore
needs only information vectors of weight at most nine. Each information set
has

\[
\sum_{j=0}^{9}\binom{38}{j}=227{,}881{,}004
\]

such vectors. An optimized fixed-width C++ enumeration searched the balls
around zero and one. It found the exact counterexample

\[
v=7,\qquad \chi=\mathtt{0x852ac8fcc6b27c67}.
\]

The corresponding 192-bit word has weight 146 and complement weight 46.
The character is annihilated by the product of the degree-18 and degree-20
factors, but not by either factor alone. Hence cancellation between these two
components refutes the required value-7 distance 48.

An independent verifier reconstructs the maps and replays the witness. Over
the full 32,772-node orbit, the same character has absolute bias
\(1894/524352\), approximately \(0.00361208\). Only one of the 2,731
12-node blocks has two-sided weight below 48. The witness therefore refutes
the blockwise lemma, not the global character cap.

The next proof route must amortize across block transitions. The smallest
natural experiment is a 24-node observability code or an equivalent weighted
transition bound. Its threshold must include the odd final 12-node block:
2,731 is odd. Distance 96 across each of the 1,365 complete pairs contributes
131,040 to the required global distance 131,088, leaving distance 48 for the
last 12-node block. Distance 97 across each pair would instead suffice with a
trivial remainder bound. A phase-aware argument may permit a sharper target.

## Exact 24-node result for \(C_{18}\oplus C_{20}\)

The 24-node mixed-component certificate proves two-sided distance 73 for
value 15 and distance 97 for every other nonzero value. Each 384-bit,
38-dimensional code contains ten disjoint information sets. Exact
Brouwer--Zimmermann enumeration therefore uses radius seven for value 15 and
radius nine otherwise.

The primary and independent implementations each check 32,062,999,280
information vectors. The primary implementation uses split-19 lookup tables.
The independent implementation uses distinct partitions and recursive
combination enumeration. A third implementation reconstructs the intended
code and audits all information sets and receipts.

This result proves the desired full-orbit caps for characters in
\(C_{18}\oplus C_{20}\). The 1,365 complete 24-node blocks contribute at
least 132,405 for values other than 15 and at least 99,645 for value 15.
These totals exceed the respective global requirements 131,088 and 98,316.
The final 12-node block requires no lower bound.

The result removes the degree-18 plus degree-20 obstruction. The other
minimal uncovered component supports remain open. The complete proof still
needs either corresponding 24-node certificates for those supports or a
single full-code argument.

## Endpoint-aware certificate through dimension 30

The corrected asymmetric target is now certified for every nonempty
component support of total dimension at most 30. Eight maximal subcodes cover
all 57 such supports. Exact primary and independent searches prove 24-node
two-sided distance 97 for values 1 through 14, 24-node distance 72 for value
15, and 12-node endpoint distance 36 for value 15.

Each implementation passes 128 certificate cases and checks
6,425,612,528 information vectors. A separate audit reconstructs the linear
maps, checks all 256 receipts and information-set ranks, and verifies the
support inclusion cover.

For values 1 through 14, the complete blocks contribute

\[
1{,}365\cdot97=132{,}405>131{,}088.
\]

For value 15, the complete blocks and endpoint contribute

\[
1{,}365\cdot72+36=98{,}316.
\]

The pure degree-one character `0xa685aac60e5acfb3` attains both value-15
local bounds and has full-orbit bias exactly \(5/8\). Thus the corrected
bound is tight, and the final 12-node term is necessary.

The dimension-31 frontier consists of masks `0x2c`, `0x33`, `0x4a`, and
`0x51`. Exact primary and independent searches also pass all corrected local
bounds on these four supports. Each implementation checks 7,716,616,224
information vectors across 64 cases. A separate reconstruction audits all
128 full receipts and 16 obstruction-probe receipts.

The dimension-32 frontier consists of masks `0x2d`, `0x34`, `0x4b`, and
`0x52`. Exact primary and independent searches pass all 64 corrected local
cases. Each implementation checks 10,119,775,656 information vectors. A
separate reconstruction audits all 128 full receipts and 16 obstruction
probes.

The endpoint-aware tranche now covers all 65 supports of total dimension at
most 32. Together with the earlier \(C_{18}\oplus C_{20}\) certificate, the
result proves the desired full-orbit caps for 66 of the 127 nonempty
component supports. The remaining support-level obligation consists of 61
supports.

## Structural anti-cancellation checkpoint

A bounded attempt to replace the remaining support cases by a uniform local
lemma has produced a sharper formulation, but not a proof. Five consecutive
nodes determine every 64-bit character. Four-node kernels, however, include
component supports of dimensions 61, 63, and 64, so the small-kernel quotient
route is false. The complementary dimension-32 component codes are also not
orthogonal on either relevant window; their cross-Gram matrices have ranks
between 30 and 32.

The useful exact reduction comes from packet-value linearity. For each local
coordinate \(p\) and character \(\chi\), define
\(a_p(\chi)\in\mathbb F_2^4\) by

\[
\langle a_p(\chi),v\rangle=\langle\chi,z_v(p)\rangle.
\]

If \(n_a(\chi)\) counts the coordinates with coefficient word \(a\), then
all packet-value character sums are Walsh coefficients of the same histogram:

\[
B_v(\chi)=\sum_{a\in\mathbb F_2^4}n_a(\chi)
(-1)^{\langle a,v\rangle}.
\]

Thus the next structural obligation is to characterize, or safely relax, the
attainable histograms strongly enough to prove
\(|B_v|\le190\) for \(v\ne15\) and \(|B_{15}|\le240\) on 24 nodes, together
with \(|B_{15}|\le120\) on the 12-node endpoint. A full-code coordinate
descent search with 10,000 pseudorandom restarts per case found no violation,
but this remains refutation evidence only.

The subsequent histogram-transition checkpoint refutes the histogram-only
version of this route. The coefficient at node \(t\) and slot \(s\) is
exactly nibble \(s\) of \(T^t\chi\). Two states can have the same current
nibble histogram and different successor histograms, so the histogram does
not define a transfer state.

Exact MacWilliams spectra give a second limitation. Every five-node
observation code has parameters \([80,64]\), but its two-sided distance is
only 2 or 3. All scalar sliding-window constraints therefore imply a
24-node bound of only 8--12. At the non-15 violation threshold, the exact
five-node candidate population is approximately \(2^{52.88706878}\). A
direct low-window transition enumeration is not practical.

This result leaves the full local-distance lemma open. It removes the
histogram-only transfer argument from the active proof routes. The current
certificate route should resume at dimension 33 unless a stronger quotient
retains enough ordered state to exploit the recurrence.

## Dimension-33 frontier

Exact primary and independent searches certify the four dimension-33 masks
`0x2e`, `0x35`, `0x4c`, and `0x53`. Each implementation passes all 64 local
cases and checks 12,216,115,184 information vectors. A separate Python audit
reconstructs all local codes and verifies 128 full receipts and 16 preflight
receipts.

The endpoint-aware tranche now covers all 69 component supports of dimension
at most 33. Together with the separate \(C_{18}\oplus C_{20}\) certificate,
the merged result covers 70 of 127 nonempty supports. The other 57 supports
remain open.

The next frontier is dimension 34. Its four supports require a predicted
15,745,416,320 vectors per implementation. Dimension 35 increases the
information radius, so dimension 34 is the next finite checkpoint before a
substantial cost jump.

## Component-split checkpoint

The exact split formulation writes the local code as
\(C_S=C_{S_L}\oplus C_{S_R}\). It converts cancellation into a nearest-pair
problem between two component codes.

For support `0x53` and packet value 1, an optimized exhaustive engine checks
all \(2^{33}-1\) nonzero pairs under the \(20+13\) split. It proves exact
two-sided minimum 120. Thus, the split is correct and the known margin is
real for this case.

Raw pair enumeration does not rescue full dimension. A \(32+32\) split has
\(2^{64}-1\) pairs. The measured rate projects to about 335 years for one
packet-value case, before independent verification. An unbalanced split
reduces memory but preserves the pair count. A generic native-XOR SAT
encoding also fails to close the known dimension-33 case within its bounded
diagnostic run.

The component split becomes useful only if a new lower bound removes about
20--25 search bits before exact comparison. Until such a pruning theorem is
available, the active deterministic route should target total bias across
the complete 32,772-node interval instead of minimum distance on every
24-node block.

## Global block-moment checkpoint

The first cumulative target used the orbit-average square of the signed
16-slot sum at each node. Exact enumeration proves the proposed bound on
every pure irreducible component. A mixed degree-1 plus degree-2 character
refutes it: for value 7 its node sums repeat `(-10,10,2)`, giving second
moment 68 against the proposed cap 64. Its actual character bias is only
`1/24`. The failed statistic discards cancellation at the small component
periods.

Grouping 12 nodes repairs this specific defect. Let (R_j(\chi,v)) be the
signed sum over all 192 cells in block (j). The cumulative bounds

\[
\sum_jR_j(\chi,v)^2\le2731\cdot96^2\quad(v\ne15)
\]

and

\[
\sum_jR_j(\chi,15)^2\le2731\cdot120^2
\]

would imply the required full-character caps by Cauchy--Schwarz. The bounds
are equivalent to one-sided Fourier caps (47/191) and (74/191) on the
multiset of contribution differences between pairs of cells in the same
12-node block.

A structured full-state search finds no violation. Its largest non-15 block
mean square is 5,184 against 9,216; the degree-one value-15 character attains
14,400 exactly. The search includes the known mixed-component
counterexamples and the high-dimensional witnesses from the 10,000-restart
local-distance probe. This remains diagnostic evidence. The exact statement,
the refuted precursor, and the next proof/refutation checkpoint are recorded
in `riffle_dp_2lap_g4_global_block_moment_route.md`.

## Scope

The iid conditioning reduction, contribution injectivity, Parseval identity,
20-bit projected bound, and dimension-at-most-22 certificate are rigorous.
The proposed uniform 12-node reduction is false. The exact counterexample and
its independent replay are recorded in
`riffle_dp_2lap_g4_c18_c20_counterexample.md`. The 24-node certificate repairs
that failure for \(C_{18}\oplus C_{20}\), and the endpoint-aware certificate
repairs it for every component support of total dimension at most 30. The
stronger value-15 distance-73 target is false because it discards a necessary
endpoint contribution. A full global character-cap proof remains open for 61
component supports. Even a complete character-cap proof would close only the support-33
terminal-zero row. Other placement strata and larger outer populations still
require separate ledger rows.
