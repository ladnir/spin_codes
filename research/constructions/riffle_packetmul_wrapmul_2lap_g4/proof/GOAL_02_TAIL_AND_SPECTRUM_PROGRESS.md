# Goal 02: tail closure and the remaining spectrum interface

**Outcome:** Part I is closed. Part II terminates at the independently
authenticated location-sensitive coding obstruction stated below. The
support-33 row remains open.

## Scope

This note proves the late-start part of Goal 02 and records two exact
reductions for the early-start part. It does not close the complete
support-33 row. The remaining claim is an upper bound on two sparse dual
weight enumerators.

## Response map

Let (P:mathbb F_2^{64}	omathbb F_2^{64}) be the systematic BCH parity
map, and let (operatorname{Acc}) be the 64-bit prefix accumulator. For a
stored state (s) and a common drive (x), one node emits

\[
O(s,x)=\operatorname{Acc}(s+x)
\]

and stores

\[
S(s,x)=P(O(s,x)).
\]

For two stored states (s_1,s_2) with the same drive (x), linearity gives

\[
O(s_1,x)+O(s_2,x)=\operatorname{Acc}(s_1+s_2)
\]

and

\[
S(s_1,x)+S(s_2,x)=P(\operatorname{Acc}(s_1+s_2)).
\]

Thus the common drive cancels. The stored-state difference follows

\[
T=P\circ\operatorname{Acc},
\]

whereas the emitted-output recurrence is the conjugate map

\[
U=\operatorname{Acc}\circ P.
\]

Both maps have rank 64. The primary and independent tail receipts reconstruct
these maps separately and check cancellation on all 4,096 pairs of basis
state and drive directions.

## Late-start theorem

Let (r) be the first occupied node and let (q=\lfloor r/4\rfloor). Remove
the first (4q) nodes. The remaining binary response code has length

\[
n_r=256(8193-q).
\]

Every nonzero difference trajectory contributes at least 36 bits in each
complete four-node block. Therefore its minimum distance satisfies

\[
d_r\ge 36(8193-q).
\]

The removed prefix consumes at least (36q) units of the bad-weight budget,
so the remaining decoding radius is

\[
t_r=188765-36q.
\]

The binary Johnson counting bound is positive for the first time at
(r=7060). At that node,

\[
(n_r,d_r,t_r)=(1645568,231408,125225)
\]

and the exact integer list bound is 13,749. The bound does not increase for
(7060\le r\le20975). At (r=9176), the inequality (2t_r<d_r) first holds,
so the list is unique from that node onward.

The parent certificate makes (r\ge20976) empty. Hence the newly charged
band is exactly (7060\le r\le20975), with placement probability

\[
\frac{
\binom{524352-16\cdot7060}{33}
-\binom{524352-16\cdot20976}{33}
}{\binom{524352}{33}}.
\]

Charging 26 outer words and a maximum list of 13,749 states gives an aggregate
charge whose base-two logarithm lies in

\[
[-57.103642352420,-57.103642352419].
\]

After this disjoint charge, the largest uniform bad-state list allowed for
the early region is 316,606. Therefore it suffices to prove

\[
A_{\le377530}\le316605.
\]

## Exact low-component census

The autonomous map (T) has seven squarefree irreducible components of
degrees

\[
1,2,4,9,10,18,20.
\]

For each exact nonempty component support of total dimension at most 22, the
primary certificate enumerates every initial state. It follows each complete
orbit and evaluates all 32,772 cyclic response windows by rolling sums. The
35 supports contain 15,650,667 states in total. None has response weight at
most 377,530. The smallest response weight in this census is 917,616.

An independent Python implementation reconstructs (P) by binary Gaussian
elimination, implements the accumulator by bit-by-bit prefix parity, and
repeats all 35 orbit censuses. It reproduces every support count, cycle count,
minimum weight, and minimum witness state.

Consequently, every state counted by (A_{\le377530}), if any, has exact
irreducible-component support dimension at least 23. This is a restriction on
the remaining search space, not a bound on its cardinality.

## Sixth-moment reduction

Write the response code as

\[
J(u)=(\langle a_1,u\rangle,\ldots,\langle a_n,u\rangle),
\qquad n=2097408.
\]

The coordinate census proves that all (a_i) are nonzero and pairwise
distinct. For uniform (u\in\mathbb F_2^{64}), let

\[
X=\operatorname{wt}(J(u)),
\qquad
S=n-2X.
\]

Let (A_j^\perp) denote the number of weight-(j) dual words. Expanding
(S^6) counts ordered sextuples of coordinate columns whose XOR is zero.
Because the columns are nonzero and distinct, the odd-multiplicity support of
such a sextuple has size 0, 4, or 6. Exact multiplicity counting gives

\[
\mathbb E[S^6]
=15n^3-30n^2+16n
+(360n-960)A_4^\perp
+720A_6^\perp.
\]

For (X\le377530), the distance below the mean (n/2) is at least

\[
a=671174.
\]

Markov's inequality applied to ((X-n/2)^6=S^6/64) therefore bounds the
number of low codewords, including the zero state. The required nonzero
spectrum target follows from the exact sufficient inequality

\[
2^{64}\bigl((360n-960)A_4^\perp+720A_6^\perp\bigr)
\le
1849740699995876325368244902659070749106176.
\]

The baseline term alone bounds the total low-codeword count by 436. The
displayed sparse-dual budget is correspondingly large: if (A_6^\perp=0), it
allows (A_4^\perp\le132802503221863); if (A_4^\perp=0), it allows
(A_6^\perp\le139270339268776862126). The proof needs the joint inequality,
not either one-variable illustration.

The response code is not six-wise independent. The coordinate columns at

\[
(0,2),(1,62),(1,63),(2,1)
\]

XOR to zero and give an authenticated dual word of weight 4. The independent
receipt reconstructs the four columns and obtains zero again. Thus the next
step must count or bound sparse dual relations; it cannot assume their
absence.

## Exact weight-4 dual census

Every unordered pair of response coordinates has the form

\[
\{W^t e_i,W^{t+d}e_j\},
\qquad W=U^\mathsf T.
\]

For (d>0), normalize the pair by removing the common factor (W^t). For
(d=0), retain only (i<j). The finite response interval therefore has
134,232,032 normalized pair types. A type with gap (d) contributes a cyclic
phase interval of length (32772-d) in the orbit of
(e_i+W^d e_j).

The primary certificate decomposes each pair sum into the seven irreducible
components of (W). Component discrete logarithms give a collision-free
orbit key and a phase modulo the exact orbit period. It then counts equal
pair sums by sweeping the cyclic phase intervals within each orbit key.

The independent certificate changes three implementation choices. It orders
the components in the opposite direction, chooses a different Krylov seed in
every component, and replaces the sweep line with direct pairwise cyclic-arc
intersection. Both executions obtain

\[
\sum_x \binom{r_x}{2}=11196069,
\]

where (r_x) is the number of unordered coordinate pairs with XOR (x).
Every weight-4 dual support has three pair partitions, so

\[
A_4^\perp=11196069/3=3732023.
\]

The 134,232,032 normalized types occupy 98,013,821 orbit keys. No occupied
key contains more than ten normalized types.

Substitution into the sixth-moment budget leaves the exact sufficient target

\[
A_6^\perp\le139270335354994389965.
\]

Thus the weight-4 term is closed. Only the weight-6 dual count remains in the
moment route.

## Pair-sum energy interface

Let \(r(x)\) count unordered pairs of distinct response-coordinate columns
whose XOR equals \(x\). The exact pair census gives the multiplicity
histogram

\[
\begin{array}{c|rrrrr}
r(x)&1&2&3&4&5\\ \hline
\#x&2199542387068&6149917&1114128&229358&32762.
\end{array}
\]

Consequently,

\[
\sum_xr(x)=2199559110528,
\qquad
\sum_xr(x)^2=2199581502666.
\]

Define the triangle energy

\[
T(r):=\sum_{x,y}r(x)r(y)r(x+y).
\]

Direct edge-multiplicity counting gives

\[
T(r)=6\binom n3+(36n-120)A_4^\perp+90A_6^\perp.
\]

Small finite models independently replay this identity. After substituting
the exact value of \(A_4^\perp\), it suffices to prove

\[
T(r)\le12543557200693295704650.
\]

The generic Young--Cauchy estimate obtained from the authenticated
multiplicity data is

\[
T(r)\le4838109533537868620667648.
\]

It exceeds the sufficient target by a factor of approximately
\(385.704745\). Thus the pair-sum histogram alone cannot close Goal 02; a
successful argument must use the locations of the pair sums.

This limitation is exact. Place multiplicity one on every nonzero point of a
41-dimensional subspace \(K\subset\mathbb F_2^{42}\). The rest of the
authenticated histogram fits on distinct points outside \(K\). The resulting
function has \(r(0)=0\), the same complete histogram, the same \(L_1\) and
\(L_2\) values, and maximum multiplicity five. The terms with distinct
nonzero \(x,y\in K\) alone give

\[
T(r)\ge(2^{41}-1)(2^{41}-2)
=4835703278451919629058050.
\]

This is \(385.512913\) times the required cap. Primary and independent
arithmetic receipts replay the construction. The countermodel is not the
actual response pair-sum function; it proves that location data are
mathematically necessary.

## Spectral-bias reduction

For \(u\in\mathbb F_2^{64}\), define

\[
B(u):=n-2\operatorname{wt}(J(u)).
\]

Fourier transformation of the pair-sum function gives

\[
\lambda_u:=\widehat r(u)=\frac{B(u)^2-n}{2},
\]

and

\[
T(r)=2^{-64}\sum_u\lambda_u^3.
\]

The low-component census isolates four nonzero exceptions. The pure
degree-one state has absolute bias \(196632\). The three pure degree-two
states have absolute bias at most \(262176\). Charging these states and the
trivial character separately leaves the sufficient pointwise target

\[
|B(u)|\le106802
\]

for every other state. This is the largest even bias admitted by the
conservative Parseval charge: replacing it by \(106804\) makes that charge
exceed the triangle-energy budget.

Primary and independent receipts replay the Fourier identities, the four
exceptional-state bounds, all dependency hashes, and the exact threshold.
The low-component census already proves the target for every nonexceptional
component support of dimension at most 22. Its largest such bias is
\(104882\). The remaining pointwise theorem concerns exactly the
nonexceptional supports of dimension at least 23.

## High-support refutation search

A bounded fixed-width search covers all 92 exact component-support masks of
dimension at least 23. For each mask, it uses two deterministic starts and
six pseudorandom starts in each weight direction. Each move is the strict
best coordinate flip that preserves the exact support mask.

The search performs 149,744 exact full-response evaluations. Its largest
absolute bias is 6,482, attained on support mask `0x52`. An independent
implementation reconstructs the BCH map, accumulator, component kernels,
and physical state. It directly replays both retained extrema for every
mask, for 184 exact witness replays in total.

The first violating even bias is 106,804, so the best diagnostic witness has
margin 100,322. This result is strong refutation evidence against a
pointwise counterexample. It is not an exhaustive upper-bound certificate.

The minimum replayed high-support response weight is 1,045,463. The exact
low-component minimum is 917,616. Both exceed the difference threshold
377,530. Hence the explored states do not produce even two members of one bad
affine list. No over-cap center was found, so the reachability and setup-mass
stages of the refutation track were not triggered. This is a bounded search
outcome, not an exhaustive exclusion of an affine counterexample.

For a possible structural proof, write

\[
B(u)=\sum_{j=0}^{63}S_j(u),
\qquad
S_j(u):=\sum_{t=0}^{32771}(-1)^{\ell_j(T^tu)}.
\]

The scalar bounds \(|S_j(u)|\le1668\) for every output bit would imply
\(|B(u)|\le106752\). The best high-support bias witness has scalar maximum
584. However, a proof still needs a uniform incomplete character-sum bound
for the fixed degree-18 and degree-20 recurrence components. Generic
completion bounds do not reach 1,668, so this scalar formulation is a proof
interface rather than a theorem.

## Eight-node alternative

An independent sufficient route is an autonomous eight-node distance bound

\[
d_8\ge93.
\]

Indeed, (32772=4096\cdot8+4), so this bound and the authenticated four-node
distance 36 would imply

\[
d(J)\ge4096\cdot93+36=380964>377530.
\]

The exact decision predicate (d_8\le92) was encoded with native XOR
constraints, the two known four-node lower bounds, and totalizer cardinality
constraints. The bounded run exhausted its time limit with result `UNKNOWN`.
This receipt establishes the predicate and the generic-solver obstruction
only; it provides no distance claim. The earlier C12 solver predicates also
remain `UNKNOWN`.

## Remaining proof obligation

It is now enough to prove any one of the following statements.

1. \(A_6^\perp\le139270335354994389965\).
2. \(T(r)\le12543557200693295704650\).
3. Outside the four charged states, \(|B(u)|\le106802\).
4. The eight-node response code has distance at least 93.

The first two statements are equivalent through the authenticated exact
identity. The third is a sufficient location-sensitive spectral bound. The
fourth makes the low spectrum empty.

The precise remaining coding problem is therefore a high-support
discrepancy bound for the length-2,097,408 response code, outside one
degree-one and three degree-two states. Exact enumeration closes every
component support through dimension 22. Bounded generic solvers do not close
the eight-node predicate, and pair multiplicities without locations miss the
energy target by a factor of 385.704745. No diagnostic counterexample is
currently near the spectral threshold.

The refutation track remains active. A violation of the joint inequality
would refute this sufficient moment route, but it would not by itself refute
the construction. A construction-level refutation still requires an
over-cap affine list at a reachable driven center and a probability-relevant
setup mass.

## Authenticated artifacts

- `receipts/goal02_tail_list_primary.json`
- `receipts/goal02_tail_list_independent.json`
- `receipts/goal02_low_components_primary.json`
- `receipts/goal02_low_components_independent.json`
- `receipts/goal02_sixth_moment_primary.json`
- `receipts/goal02_sixth_moment_independent.json`
- `receipts/goal02_dual_a4_primary_raw.json`
- `receipts/goal02_dual_a4_independent_raw.json`
- `receipts/goal02_dual_a4_audit.json`
- `receipts/goal02_a6_energy_interface.json`
- `receipts/goal02_pair_histogram_countermodel.json`
- `receipts/goal02_pair_histogram_countermodel_independent.json`
- `receipts/goal02_spectral_bias_interface.json`
- `receipts/goal02_spectral_bias_independent.json`
- `receipts/goal02_high_support_bias_search_raw.json`
- `receipts/goal02_high_support_bias_search_audit.json`
- `receipts/goal02_c8_native_xor.json`
- `receipts/goal02_coding_obstruction_audit.json`
- `receipts/goal02_coding_obstruction_independent.json`
