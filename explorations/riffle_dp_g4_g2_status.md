# Riffle DP g=4: G2 proof and refutation status

## Question

G2 asks whether every outer block-weight-three word has inner-output failure
probability small enough for the first-moment budget. The probability is over
the uniform permutation of the 524,352 four-bit packets.

The terminal event passed G1. G2 must also cover nonterminal placements under
the fixed inner recurrence.

## Fixed recurrence

For inner node (i), let (U_i\in\mathbb F_2^{64}) contain the permuted input
packets. Let (S_{i-1}\in\mathbb F_2^{64}) be the incoming state. The forward
recurrence is

\[
V_i=\operatorname{Acc}(U_i+S_{i-1}),
\qquad
S_i=P V_i.
\]

Here, (P) is the systematic right half of the extended BCH encoder. The map
\(\operatorname{Acc}\) is the fixed prefix-XOR transform.

When (U_i=0), the state follows the invertible linear map

\[
T(S):=P\operatorname{Acc}(S).
\]

The emitted block is \(\operatorname{Acc}(S)\). Therefore, the BCH distance
does not directly lower-bound two consecutive emitted blocks. The BCH word is
\((V_i,P V_i)\), while the next emitted block is
\(\operatorname{Acc}(P V_i)\). Accumulation is not weight-preserving.

## Exact proof evidence

The minimal polynomial of (T) is

\[
x^{64}+x^{63}+x^{61}+x^{60}+x^{58}+x^{57}+x^{56}+x^{55}
+x^{54}+x^{53}+x^{51}+x^{50}+x^{49}+x^{47}+x^{45}+x^{44}
+x^{43}+x^{42}+x^{41}+x^{40}+x^{39}+x^{31}+x^{24}+x^{22}
+x^{19}+x^{18}+x^{17}+x^{16}+x^{15}+x^{13}+x^{12}+x^{11}
+x^{10}+x^6+x^5+x^4+x^2+1,
\]

whose hexadecimal representation is `0x1b7eebf80814fbc75`. It factors over
\(\mathbb F_2\) into irreducible factors of degrees

\[
1,2,4,9,10,18,20.
\]

The corresponding root orders are

\[
1,3,15,511,1023,37449,1048575.
\]

Thus, the matrix order is 91,625,532,075. The map is invertible, reducible,
and not primitive.

An exhaustive cycle enumeration gives these minimum average emitted weights:

| Component degree | Minimum average weight per node |
|---:|---:|
| 1 | 29.000000000 |
| 2 | 28.000000000 |
| 4 | 33.600000000 |
| 9 | 31.561643836 |
| 10 | 32.031280547 |
| 18 | 31.984619082 |
| 20 | 32.000030518 |

The audit also exhausts six maximal combinations of components whose total
dimension is at most 20. For cycles that activate every component in the
listed combination, the minimum average is:

| Component degrees | Dimension | Minimum full-support average |
|---|---:|---:|
| 20 | 20 | 32.000030518 |
| 18, 2 | 20 | 31.958877407 |
| 18, 1 | 19 | 31.964244706 |
| 10, 9, 1 | 20 | 32.000483976 |
| 10, 4, 2, 1 | 17 | 31.940371457 |
| 9, 4, 2, 1 | 16 | 32.031441618 |

This enumeration covers about 3.8 million states. No combined-component cycle
in the covered subspaces has a low average caused by cancellation.

The receipt is
`explorations/riffle_dp_g4_g2_autonomous_map.json`. Its SHA-256 digest is
`0458e349ccdb1bd773db256e7ca87b7c3bbc6be128fdbc482125a109f1764376`.
The independent verifier reconstructs the matrix, factors, orders, kernels,
cycles, and rational averages. It reports
`EXACT_G2_AUTONOMOUS_MAP_VERIFIED`.

### Finite-window lemma

Full-cycle averages do not control a finite prefix. A separate exhaustive
search handles finite autonomous windows.

Call an emitted block low when its weight is at most eight. There are exactly

\[
\sum_{j=1}^{8}\binom{64}{j}=5,130,659,560
\]

nonzero low blocks. The forward search enumerates every low block. It finds no
one-step return to the low set. It finds 11,039,336 two-step returns. Every
two-step first-return segment has weight at least 19.

An independent reverse search enumerates every low endpoint. It applies
\(Q^{-1}\) and \(Q^{-2}\), where \(Q=\operatorname{Acc}P\). The reverse search
reproduces the zero one-step returns, the two-step count, and minimum weight
19.

These facts give a deterministic bound for every nonzero autonomous segment
of length (L\geq 1):

\[
\delta_{\mathrm{aut}}(L)
\geq 9L-8-8\left\lfloor\frac{L-1}{3}\right\rfloor.
\]

To see the bound, partition a segment at successive low blocks. A short
first-return segment has length two and weight at least (19\geq18). A longer
first-return segment of length (g\geq3) has weight at least (9g-8). The
initial high prefix costs at least nine bits per node. The final segment loses
at most eight bits. At most \(\lfloor(L-1)/3\rfloor\) complete segments incur
the eight-bit loss.

Consequently, every nonzero autonomous segment of length 29,807 has weight
above (d). At the full chain length, the bound is

\[
\delta_{\mathrm{aut}}(32,772)\geq207,556=d+18,790.
\]

The proof ledger is
`explorations/riffle_dp_g4_g2_autonomous_window_bound.json`. Its SHA-256
digest is
`2a3b34da4c00b4180dfc007437aabd14b3e69ae238ca1848e27f8ef5a686eeb8`.

### Reset-pattern placement envelope

The finite-window lemma gives the linear relaxation

\[
\delta_{\mathrm{aut}}(L)\geq\frac{19L-16}{3}.
\]

Suppose an occupied node is live when its emitted block is nonzero. A reset
node emits zero. The first occupied node is live, and two reset nodes cannot
be consecutive.

For each packet support from 33 through 38, an exact generating-function
calculation counts every node placement and every abstract reset pattern that
can satisfy the relaxed distance condition. The calculation includes packet
collisions inside nodes.

For support 33, the zero-reset placement probability is at most

\[
2^{-4.475042948481}.
\]

The corresponding bounds decrease only to
\(2^{-5.143956675803}\) at support 38. Summing abstract reset patterns exceeds
one because several hypothetical patterns can label the same node placement.

This envelope is rigorous but far too loose for the first-moment target. It
shows that a proof cannot use the deterministic autonomous bound and reset
rarity alone. The proof must exploit the actual packet values, nibble slots,
or distribution of reached states.

### Reachable first segments

The arbitrary-state window bound is pessimistic for states reached directly
from zero. Exact enumeration now covers every initial node of packet support
one, two, or three:

| Packets in first node | Concrete drives | Latest crossing of (d) |
|---:|---:|---:|
| 1 | 240 | 5,926 |
| 2 | 27,000 | 5,943 |
| 3 | 1,890,000 | 5,945 |

Combining these exact prefix tables with the autonomous relaxation improves
the support-33 zero-reset envelope to
`2^-6.767496031893`. First-node support at least four contributes only
`2^-35.172380259029`, but the one-packet row still dominates. The same
combined envelope is `2^-7.433819992146` at support 38.

Thus, exact first-segment reachability gains about 2.3 bits. It does not
approach the full G2 budget.

### Observability diagnostic

The fixed autonomous map defines a binary linear observability code of
dimension 64. A deterministic hill search found the exact witness

```text
state   0x6331abf1b618efab
length  5956 nodes
weight  188730
```

The witness is 36 bits below (d). Its one-node preimage has nibble support
13 and does not match any authenticated support-33 value multiset. Therefore,
the witness refutes a universal 5,956-node autonomous bound, but it does not
refute Riffle DP g=4.

### Component mixing diagnostic

For a fixed nonzero nibble value, sample one of 16 slots and one of 32,772
node exponents. Exact Walsh transforms give the bias of every character in
each irreducible component and every tested combination through dimension 20.

The largest single-packet biases in the individual components are:

| Component degree | Maximum absolute bias |
|---:|---:|
| 1 | 0.625000000 |
| 2 | 0.250000000 |
| 4 | 0.133364610 |
| 9 | 0.060753082 |
| 10 | 0.001121384 |
| 18 | 0.004359667 |
| 20 | 0.005244569 |

The slow bias is confined to the low-degree factors. A former calculation for
the 14 authenticated support-33 profiles reported these values:

```text
degree 1 projection                 2.452821e-19
degree 2 projection                 1.384094e-27
degree 20 projection                9.335884e-88
all tested combined projections     <= 2.452821e-19
```

The Walsh spectra are exact, but the without-replacement reduction was wrong.
It multiplied separate bounds on adaptive conditional biases. A balanced
two-draw example disproves that multiplication. Therefore, the displayed
total-variation values are not rigorous bounds. The corrected iid-plus-
distinctness calculation gives a rigorous degree-20 projection-zero bound,
but full-state characters remain the central proof gap.

A full-character hill search performs 1,605 exact bias evaluations after
local optimization. It finds no character stronger than the degree-one
character. The largest reported bias is exactly 0.625 for nibble value 15.
This maximum is diagnostic because hill climbing does not exhaust all
\(2^{64}-1\) nonzero characters.

## Exact refutation evidence

The one-packet orbit search covers all 240 nonzero nibble impulses. Without a
later input, every orbit exceeds (d=188,766) after 5,871 to 5,926 nodes.

Only six two-packet turnoffs occur before the threshold. Every turnoff places
the second packet in the immediately following node. Four authenticated outer
support-33 words contain a compatible value pair.

The resulting exact disjoint lower-bound family contributes

\[
2^{-124.233993125178}.
\]

This family does not refute the candidate.

The component audit also searches a natural algebraic refutation route. It
finds no low-average cycle in any individual component or covered combination.

The full-chain distinct-node three-packet search tests 1,887,552,000 exact
relations. It finds 27 algebraic turnoffs. All 27 occur by node six; extending
the search from node 6,000 through node 32,771 produces no additional turnoff.

Sixteen authenticated support-33 outer words contain a compatible value
triple. The resulting disjoint exact lower family contributes

\[
2^{-134.276933547563}.
\]

This family is weaker than the two-packet family and does not refute the
candidate. The verifier replays each recurrence and recomputes every rational
placement probability.

The collision search covers initial and reset nodes containing one or two
packets. Across the full chain, it finds:

| Initial packets | Reset packets | Exact turnoffs |
|---:|---:|---:|
| 1 | 1 | 6 |
| 1 | 2 | 24 |
| 2 | 1 | 64 |
| 2 | 2 | 445 |

All 539 turnoffs occur by node eight. Every authenticated support-33 word
contains at least one compatible collision turnoff. The strongest disjoint
collision family contributes

\[
2^{-123.048186206738}.
\]

Allowing repeated two-packet episodes gives
\(2^{-123.049568554040}\). Both results remain far below (2^{-40}).

An authenticated support-33 profile gives a three-node terminal cluster of
weight 188,505 over 5,940 nodes. The weight is 261 below (d). An independent
verifier replays the recurrence and checks the packet-value multiset.

Every later terminal translation truncates the same output trajectory. This
gives 5,938 exact bad placements for each compatible outer word. Three outer
families have the required value profile. For one family, the counted event
has probability

\[
2^{-555.961662409912}.
\]

The cluster witness disproves deterministic distance for every permutation.
Its counted probability is far too small to refute the random-permutation
claim.

The exact single-transposition basin contains 919 distinct assignments and
5,419,150 disjoint physical placements. Its aggregate first-moment
contribution across the three compatible outer families is

\[
2^{-544.542826136930}.
\]

Thus, enlarging only the local cluster cannot approach the refutation target.

### Broad suffix diagnostic

Condition on all 33 packets of one authenticated profile occupying uniformly
random distinct cells in the final 5,940 nodes. Among 100,000 trials, 79,583
have output weight at most (d). The estimated conditional failure probability
is 0.79583. The corresponding three-family contribution estimate is
`2^-80.0607`.

Scans through 8,000 nodes keep the unconditional estimate near `2^-80`.
Longer-window failures mainly place their packets inside an effective shorter
suffix. These Monte Carlo results are diagnostic, not probability bounds.

An exact one-data enumeration checks all 631,767,040 pairs \((i,t)\) for which
both \(t\) and \(x^i t\) have local support at most 13. An independent verifier
confirms 85,828 outer words:

| Total packet support | Exact words |
|---:|---:|
| 33 | 26 |
| 35 | 36 |
| 36 | 3,233 |
| 37 | 510 |
| 38 | 933 |
| 39 | 81,090 |

A 100,000-trial suffix experiment for each nonempty support gives conditional
failure estimates from 0.76465 through 0.79766. The estimated aggregate
contribution of this exact outer population is `2^-76.01295`. The population
is exact, but the conditional inner probabilities remain diagnostic.

## Remaining gap

The exact evidence does not close G2. The cluster witness rules out a
universal deterministic distance proof. A full proof must bound the measure
of bad packet permutations, including reset-compatible placements.

The three-packet search covers one nonzero nibble in each of three distinct
nodes. It does not cover two packets in one node or resets that require four or
more packets. The new collision search covers up to two packets in the initial
and reset nodes. Episodes that use intermediate active nodes or at least three
packets in one node remain open.

The old profile certificates still cannot apply because they assume
independent lane or state permutations that this candidate does not sample.
The reset-pattern envelope demonstrates that deterministic gap aggregation is
insufficient. The next proof artifact must control the distribution of reached
states under random packet locations and fixed nibble slots. The component
Walsh calculations support this route. The next artifact must bound every
full-state Fourier character, including characters that span more than 20
component dimensions. It must then connect state mixing to total output weight.

The broad suffix evidence sharpens the refutation target. The terminal suffix
used by G1 is not the only dangerous region. A useful lower family must combine
broad suffix behavior with the large population of support-39 coefficient
classes. The exact one-data slice alone remains about 36 bits below the
refutation threshold.
