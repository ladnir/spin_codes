# Riffle DP-2Lap g=4: initial diagnostic

## Candidate

This document concerns **Riffle DP-2Lap g=4**. It does not modify the frozen
candidate **Riffle DP g=4@g0-v1**.

Let (F(x)) denote one zero-state inner lap. Let (L(x)) denote its terminal
state. Let (J(s)) denote the output on an all-zero word with initial state
(s). Riffle DP-2Lap g=4 computes

\[
y:=F(x),\qquad G(x):=F(y)\mathbin\oplus J(L(x)).
\]

Equivalently, the implementation retains the first-lap terminal state, wraps
to the first logical node, and processes the overwritten first-lap word. It
does not reset the state or sample another permutation.

The identity

\[
G=F^2\mathbin\oplus JL
\]

follows from linearity. The diagnostic checks the identity on 448 pseudorandom
chains with lengths from 1 through 127.

## Exact cluster replay

The authenticated three-node cluster has first-lap weight 188,505 over the
final 5,940 nodes. This reproduces the frozen-candidate receipt exactly. The
first lap terminates in state `0x4e023d608e282928`.

The retained state produces 858,256 bits in the preceding 26,832-node prefix.
The second-lap suffix contributes another 190,229 bits. The final two-lap
weight is therefore 1,048,485.

These values are exact evaluations of the scalar recurrence.

## Broad suffix diagnostic

The experiment conditions on the 33 packets of the authenticated profile
occupying distinct cells in the final 5,940 nodes. It samples the packet cells
and the assignment of the fixed nibble multiset pseudorandomly.

Among 100,000 trials:

- the first lap has 79,535 outputs below the distance threshold;
- no first-lap terminal state is zero;
- no final two-lap output is below or equal to the threshold;
- the minimum final weight is 1,045,691;
- the mean final weight is 1,048,701.71035.

The second-lap prefix weight ranges from 855,760 through 861,702. The
second-lap suffix weight ranges from 188,817 through 191,423. Thus, every
sampled suffix alone exceeds the required distance of 188,766.

The zero-failure 95% binomial upper confidence value for the conditional
failure probability is about (2.996\cdot10^{-5}). This statistical value is
not a rigorous upper bound for the construction.

## Receipts and scope

- Probe source: `scripts/probe_riffle_dp_2lap_g4_support33_suffix.cpp`
- Probe source SHA-256:
  `6aca094692f3a4a282eb29924b2417f73f69b45de328739bdecd88a3a2ad8bfd`
- Result artifact:
  `explorations/riffle_dp_2lap_g4_support33_suffix_probe.json`
- Result artifact SHA-256:
  `c6f6a11f3d1ad27043a6fc5332a818de1c83605aae17810b43f844349fe5c0ea`

The recurrence, decomposition checks, and cluster replay are exact. The broad
suffix sample supplies refutation evidence only. It does not prove the random
permutation claim.

## Exact one-data population diagnostic

The exact outer enumeration contains 85,828 one-data words with packet
supports 33, 35, 36, 37, 38, and 39. The DP-2Lap experiment samples 100,000
suffix placements independently within each nonempty support class.

Across all 600,000 trials:

- the first lap has 467,884 outputs below the distance threshold;
- no first-lap terminal state is zero;
- every wrapped prefix reaches the distance threshold;
- no two-lap output is below the distance threshold.

The latest prefix crossing occurs at node 5,944 of the available 26,832-node
prefix. The Bonferroni simultaneous 95% value gives an aggregate contribution
upper value of about `2^-90.0222` for these six sampled strata. This value is
a statistical confidence bound, not a theorem bound.

- Probe source: `scripts/probe_riffle_dp_2lap_g4_one_data_suffix.cpp`
- Probe source SHA-256:
  `a1cfc1fad90431f1d060af45cdbeda1f9f2287e18d19d941c8ef7b5b9f818e59`
- Result artifact:
  `explorations/riffle_dp_2lap_g4_one_data_suffix_probe.json`
- Result artifact SHA-256:
  `b12fb7d155d589cd6cb4ca2d86b2a36ca364eda7c9767456ca30d796e2ad6cb8`

## Exact terminal-turnoff replay

The retained state does not repair a first lap whose terminal state is zero.
An exact replay shifts each known isolated turnoff so its reset node is the
final node.

All 27 distinct-node three-packet episodes remain below the distance after two
laps. Their two-lap weights range from 2 through 185. All 539 collision
episodes also remain below the distance. Their weights range from 1 through
250.

These episodes do not yet give authenticated counterexamples. An authenticated
outer word has at least 33 active packets, whereas each replay contains at
most four. The old lower-family probabilities cannot be reused because their
remaining packets generally make the first-lap terminal state nonzero.

- Verifier source: `scripts/verify_riffle_dp_2lap_g4_terminal_turnoffs.py`
- Verifier source SHA-256:
  `29f4c3efcbed013abe50452d93e64621e50cd5b2c08aca569a242bc2b35bdfad`
- Result artifact:
  `explorations/riffle_dp_2lap_g4_terminal_turnoff_replay.json`
- Result artifact SHA-256:
  `ccf9d9ff2a898bad9b5a74d52f9ae12279ede05913b0b225e5750c453b9119bd`

The evidence now isolates the proof obligation. Ordinary suffix states spread
far beyond the required distance. Exact algebraic turnoffs remain dangerous.
The proof stream should bound the probability of an authenticated placement
whose first-lap terminal state is zero or belongs to a comparably exceptional
low-growth set.

## Exact nonzero-terminal certificate

Let \(o_i\in\mathbb F_2^{64}\) denote the output at one node of the autonomous
wrapped prefix. Consecutive outputs satisfy

\[
o_{i+1}=U(o_i),\qquad U(o):=\operatorname{Acc}(P(o)).
\]

The map \(U\) is invertible. Suppose a three-node window has total weight at
most 24. At least one of its outputs then has weight at most 8. An exhaustive
enumeration checks every nonzero 64-bit output of weight at most 8 and all
three window alignments containing that output.

The enumeration covers exactly 5,130,659,560 outputs. It finds no window of
weight at most 24. Therefore, every nonzero autonomous three-node window has
weight at least 25.

The old 5,940-node suffix leaves a 26,832-node wrapped prefix. Partition that
prefix into 8,944 consecutive three-node blocks. If the first-lap terminal
state is nonzero, invertibility keeps every block state nonzero. Hence

\[
\operatorname{wt}(\text{wrapped prefix})
\ge 8{,}944\cdot25
=223{,}600
>188{,}766.
\]

The deterministic margin is 34,834 bits. The suffix input and its correlation
with the terminal state no longer matter. More generally, any zero prefix of
at least 22,653 nodes is sufficient. This covers inputs supported within the
final 10,119 nodes whenever the first-lap terminal state is nonzero.

An independent verifier recomputes both neighboring outputs directly for all
5,130,659,560 enumerated words. It reproduces zero rejected windows and the same
minimum enumerated window weight of 25.

- Certificate source:
  `scripts/certify_riffle_dp_2lap_g4_three_node_block.cpp`
- Certificate source SHA-256:
  `a6fc6a1c0febe89d6d9e88d26d09a7586135a63e96bafdf1821ae3b38067504c`
- Independent verifier:
  `scripts/verify_riffle_dp_2lap_g4_three_node_block.cpp`
- Independent verifier SHA-256:
  `65d4ed4287b484587602a9b09ed0b53505daf521d336dcf7b1add3667bd09e5b`
- Result artifact:
  `explorations/riffle_dp_2lap_g4_three_node_block_certificate_w8.json`
- Result artifact SHA-256:
  `578cd95dd1c8798940f1cedf5543954f34ba84f1402ac0d0a9151cd23dc90b90`

Thus, the low-growth terminal set for this 26,832-node prefix is exactly the
zero state. The active proof obligation is now to bound

\[
\Pr[L(X)=0]
\]

under the existing random packet permutation and the authenticated outer-word
distribution relevant to each proof stratum.

## Bounded authenticated double-turnoff search

A primitive double turnoff starts both laps in state zero and returns both
laps to state zero. Such episodes can be concatenated without interaction.
Their first-lap and second-lap weights add across episode boundaries.

The known catalog contains 17 three-packet double turnoffs and 342 collision
double turnoffs. After deduplication by packet-value multiset, the catalog has
88 episode types. The authenticated support-33 receipts have 14 distinct
packet-value multisets across 26 outer words.

An exact subtraction search found no partition of an authenticated multiset
into known episode types. A separate verifier built the additive closure of
the allowed episode types below each target. The closure calculation also
found zero reachable targets. The largest closure contained 1,210 states.

This negative result eliminates one direct refutation mechanism: no known
two-to-four-packet double turnoffs can tile an authenticated support-33 word.
It does not exclude a primitive double turnoff with five or more packets. It
also does not exclude cancellation between components that are not themselves
turnoffs.

- Search source:
  `scripts/search_riffle_dp_2lap_g4_authenticated_double_turnoffs.py`
- Search source SHA-256:
  `1ac037218d9ae295e5f11314f6978f54abf9b10243b27d18622d8eb02dcf9724`
- Search artifact:
  `explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_search.json`
- Search artifact SHA-256:
  `b6ee5bdb026a5a3ca3c10e6170d4cb64a57eb9dd59aa1c52bab438039403f4a1`
- Independent verifier:
  `scripts/verify_riffle_dp_2lap_g4_authenticated_double_turnoffs.py`
- Independent verifier SHA-256:
  `f308482fc7a7d72eefb8916e2083290ca0d47560cbd60a683333bbf0822d49ef`
- Verification artifact:
  `explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_verification.json`
- Verification artifact SHA-256:
  `e4624d5d04654674c1162a2cc502749badeca872838432465f5035223b5d843d`

The additive closure also identifies the smallest missing primitive episode.
No primitive episode of five or fewer packets can complete an authenticated
profile when combined with the known catalog. The minimum residual size is
six. Exactly two targets attain it, with residual value multisets
\(10^3 13^3\) and \(5^3 13^3\).

An exact bounded search covers every two-node split of both residuals. It
checks 360,743,168 orbit steps and every reset offset through node 32,771.
The profile \(5^3 13^3\) has two double turnoffs in the four-plus-two split.
Both have two-lap weight 29.

One turnoff combines with eight known episodes to authenticate the outer word
`terminal-w3-s33-one-i0-tfbfc4e34c1f4398c`. The resulting 33-packet input
occupies 39 nodes and has two-lap weight 56. Thus, deterministic distance is
false for Riffle DP-2Lap g=4.

A disjoint lower family inserts arbitrary zero-node gaps before, between, and
after the nine episodes. Its exact probability is between
\(2^{-473.406940769233}\) and \(2^{-473.406940769232}\). This family does not
refute the random-permutation claim.

- Six-packet search source:
  `scripts/search_riffle_dp_2lap_g4_six_packet_two_node_turnoffs.cpp`
- Six-packet search source SHA-256:
  `3cafe085e6545e14bfdcf428e02da1dee42437ec605443796054f5b03799d680`
- Six-packet search artifact:
  `explorations/riffle_dp_2lap_g4_six_packet_two_node_turnoffs.json`
- Six-packet search artifact SHA-256:
  `573c8eccff4c57544922be2283ec7f641f1f47e6b5e2bd87df9b2b7b51431d57`
- Authenticated-family source:
  `scripts/construct_riffle_dp_2lap_g4_authenticated_six_packet_turnoff.py`
- Authenticated-family source SHA-256:
  `3c9a46ba6bc311f11741147b4e0ac89e8c368acf4e4d2fe4abbe01faa52d8f03`
- Authenticated-family artifact:
  `explorations/riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json`
- Authenticated-family artifact SHA-256:
  `c0e55110712d48689f47416301cf9b5f7a286cc1fc8c5c7673e7879f9507d312`
- Independent verifier:
  `scripts/verify_riffle_dp_2lap_g4_authenticated_six_packet_turnoff.py`
- Independent verifier SHA-256:
  `927058e929e068d1e2639e4b736fa8479b70ddc57515fd8d9a2b820183173644`
- Verification artifact:
  `explorations/riffle_dp_2lap_g4_authenticated_six_packet_turnoff_verification.json`
- Verification artifact SHA-256:
  `f807ab8cd6871613783df2d9ba479e4a5f915da74964aa60f1210af14893bfa8`

The remaining six-packet refutation geometries use at least three active
nodes. Broader cancellation families need not decompose into turnoffs.

## Global-permutation Fourier gate

The packet permutation is now used directly. Sample the labeled packet cells
independently, then condition on their being distinct. This conditioning gives
the exact ordered-injection law of the global permutation. It also gives the
valid inequality

\[
\Pr_{\mathrm{perm}}[E]
\le \Pr_{\mathrm{iid}}[E]/\Pr[\mathrm{distinct}].
\]

The earlier component-mixing receipt used a different reduction. It
multiplied bounds on adaptive conditional biases. That step is invalid. An
exact two-draw counterexample underestimates the character expectation by a
factor of 524,351. The old receipt is now labeled `INVALIDATED_BOUND`.

For iid cells, Fourier coefficients factor exactly. Each nonzero packet value
also maps the 524,352 cells injectively into terminal states. Parseval can
therefore pay for two repeated packet values exactly.

If the diagnostic full-character maxima were certified, this method would
bound the aggregate terminal-zero probability of all 26 support-33 words by a
value in

\[
[2^{-59.295584244770},2^{-59.295584244769}].
\]

The exact maxima are not required. The caps \(B_{15}\le5/8\) and
\(B_v\le1/2\) for \(v\ne15\) already give an aggregate bound in

\[
[2^{-44.464429949134},2^{-44.464429949133}].
\]

Because \(32{,}772=12\cdot2{,}731\), these caps would follow from two finite
observability-code statements. The value-15 code would need two-sided
distance 36 in each 12-node block. Every other value would need two-sided
distance 48. The known degree-one character attains the value-15 bound.

An exact 20-bit quotient calculation independently gives the rigorous, but
insufficient, aggregate bound `2^-15.298107519745`. It validates the corrected
iid conditioning route without assuming any full-character maximum.

An exact local component certificate covers every character whose component
support has total dimension at most 22. All 274,291,320 character/value
evaluations pass the proposed two-sided distance thresholds.

The first targeted larger support refutes the uniform 12-node lemma. For
value 7, the mixed \(C_{18}\oplus C_{20}\) character
`0x852ac8fcc6b27c67` produces a 192-bit word of weight 146, hence two-sided
weight 46 rather than 48. An independent verifier confirms the character's
component membership and observed word. Its full-orbit bias is only
`1894/524352`, and only one of 2,731 blocks violates the local threshold.
Thus, the counterexample invalidates the blockwise proof route, not the
global character cap or the construction. The next proof target must
amortize over transitions, naturally beginning with 24-node blocks.

That 24-node target is now certified for \(C_{18}\oplus C_{20}\). Exact
primary and independent searches prove two-sided distance 73 for value 15
and 97 for every other nonzero value. The 1,365 complete 24-node blocks alone
therefore prove the desired full-orbit character caps within this subcode.
The remaining gap consists of the other uncovered component supports and the
full 64-dimensional code.

The corrected dimension-at-most-30 extension is now complete. Eight maximal
subcodes cover all 57 nonempty component supports in that range. Exact
primary and independent searches prove 24-node distance 97 for values 1
through 14, 24-node distance 72 for value 15, and 12-node endpoint distance
36 for value 15. Each search checks 6,425,612,528 information vectors. A
separate audit checks all 256 receipts, information-set ranks, support
inclusions, and executable hashes.

The pure degree-one character `0xa685aac60e5acfb3` attains weight 72 in the
24-node code and 36 in the endpoint code. Its full-orbit bias is exactly
`5/8`, so both corrected value-15 bounds are tight.

The dimension-31 and dimension-32 frontiers also pass exact primary and
independent certification. The complete endpoint-aware tranche now covers 65
supports through dimension 32. Together with the earlier
\(C_{18}\oplus C_{20}\) result, exact certificates cover 66 of 127 nonempty
component supports. The other 61 supports remain open.

A subsequent structural checkpoint refutes two simple ways of disposing of
those supports. Four-node observation kernels contain high-dimensional and
full component supports, and complementary dimension-32 component codes are
not orthogonal on the local windows. The checkpoint instead represents the
15 packet-value character sums as the nontrivial Walsh coefficients of one
attainable 16-bin coefficient histogram. A 10,000-restart full-code search
found no local-distance violation; its smallest non-15 two-sided weight was
120, while value 15 reproduced the tight weights 72 and 36. These search
results are diagnostic, not a proof.

The follow-up histogram-transition audit identifies the coefficient words as
the ordered nibbles of \(T\chi,\ldots,T^{24}\chi\). It also gives an exact
counterexample to histogram-only propagation: two states have the same nibble
histogram and different successor histograms. Exact five-node spectra have
two-sided distance only 2 or 3, and their low-weight candidate populations
are too large for direct transition enumeration. This closes the simple
histogram-transfer route, but it does not refute the 24-node distance lemma.

The dimension-33 frontier subsequently passes exact primary and independent
certification. The endpoint-aware tranche now covers all 69 supports through
dimension 33. The separate \(C_{18}\oplus C_{20}\) result raises merged
coverage to 70 of 127 supports, leaving 57 supports open.

An exact component-split prototype proves minimum two-sided weight 120 for
support `0x53`, packet value 1, by enumerating all \(2^{33}-1\) pairs. The
same raw method requires \(2^{64}-1\) pairs at full dimension and is not a
practical completion strategy. A generic SAT encoding also fails to close
the known case within its diagnostic budget. The split remains a useful
cancellation interface only if a stronger pruning theorem is found.

- Proof-route note:
  `explorations/riffle_dp_2lap_g4_global_permutation_proof_route.md`
- Parseval-gate source:
  `scripts/analyze_riffle_dp_2lap_g4_terminal_zero_parseval_gate.py`
- Parseval-gate artifact:
  `explorations/riffle_dp_2lap_g4_terminal_zero_parseval_gate.json`
- Rigorous quotient source:
  `scripts/bound_riffle_dp_2lap_g4_support33_projection.py`
- Rigorous quotient artifact:
  `explorations/riffle_dp_2lap_g4_support33_projection_bound.json`
- Invalid-step audit:
  `explorations/riffle_dp_g4_without_replacement_bias_audit.json`
- Local component certificate:
  `explorations/riffle_dp_2lap_g4_local_character_components.json`
- Mixed-component counterexample:
  `explorations/riffle_dp_2lap_g4_c18_c20_counterexample.md`
- Counterexample search receipt:
  `explorations/riffle_dp_2lap_g4_c18_c20_v07.json`
- Independent verification receipt:
  `explorations/riffle_dp_2lap_g4_c18_c20_counterexample_verification.json`
- 24-node mixed-component certificate:
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_certificate.md`
- 24-node aggregate audit:
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_audit.json`
- Dimension-30 tranche counterexample:
  `explorations/riffle_dp_2lap_g4_component24_dimension30_counterexample.md`
- Audited dimension-30 ledger:
  `explorations/riffle_dp_2lap_g4_component24_dimension30_counterexample_ledger.json`
- Corrected endpoint-aware dimension-30 certificate:
  `explorations/riffle_dp_2lap_g4_component_endpoint_dimension30_certificate.md`
- Corrected endpoint-aware audit and support ledger:
  `explorations/riffle_dp_2lap_g4_component_endpoint_dimension30_audit.json`
- Dimension-31 endpoint-aware certificate:
  `explorations/riffle_dp_2lap_g4_component_dimension31_certificate.md`
- Dimension-31 independent audit:
  `explorations/riffle_dp_2lap_g4_component_dimension31_audit.json`
- Dimension-32 endpoint-aware certificate:
  `explorations/riffle_dp_2lap_g4_component_dimension32_certificate.md`
- Dimension-32 independent audit:
  `explorations/riffle_dp_2lap_g4_component_dimension32_audit.json`
- Structural anti-cancellation report:
  `explorations/riffle_dp_2lap_g4_structural_anticancellation_checkpoint.md`
- Structural exact audit:
  `explorations/riffle_dp_2lap_g4_structural_checkpoint.json`
- Full-code refutation probe:
  `explorations/riffle_dp_2lap_g4_full_local_distance_probe.json`
- Histogram-transition report:
  `explorations/riffle_dp_2lap_g4_histogram_transition_checkpoint.md`
- Histogram-transition exact audit:
  `explorations/riffle_dp_2lap_g4_histogram_transition.json`
- Dimension-33 endpoint-aware certificate:
  `explorations/riffle_dp_2lap_g4_component_dimension33_certificate.md`
- Dimension-33 independent audit:
  `explorations/riffle_dp_2lap_g4_component_dimension33_audit.json`
- Component-split checkpoint:
  `explorations/riffle_dp_2lap_g4_component_split_checkpoint.md`
- Component-split exact audit:
  `explorations/riffle_dp_2lap_g4_component_split_pairs_audit.json`

## Rigorous one-data suffix row

The three-node certificate also closes one placement stratum without a bound
on the terminal-state distribution. A zero-input prefix of 22,653 nodes
contains 7,551 complete three-node blocks. For every nonzero terminal state,
these blocks contribute at least

\[
7{,}551\cdot25=188{,}775>188{,}766.
\]

Therefore, a bad word whose packets all occupy the final 10,119 nodes must
have first-lap terminal state zero. We can ignore that additional condition
and charge the complete suffix-placement event.

For packet support \(h\), the exact placement probability is

\[
\frac{\binom{161{,}904}{h}}{\binom{524{,}352}{h}}.
\]

The exact one-data population has 85,828 words through support 39. Summing
this probability with the exact support multiplicities gives an aggregate
upper bound in the interval

\[
[2^{-48.293011187781},2^{-48.293011187780}].
\]

This row is below \(2^{-40}\) by more than 8 bits. It covers only placements
fully contained in the final 10,119 nodes. It does not bound the remaining
placements of these words.

- Bound source: `scripts/bound_riffle_dp_2lap_g4_one_data_suffix.py`
- Bound source SHA-256:
  `6c62ee2f32b718dceceb5be5831e64eb4c523d66323717706e4d55db4a14b2ae`
- Bound artifact:
  `explorations/riffle_dp_2lap_g4_one_data_suffix_bound.json`
- Bound artifact SHA-256:
  `a8f2404d8a4be90e38ce573af7617f117f44c68adf92faf82efa0bae51101c40`
- Independent verifier:
  `scripts/verify_riffle_dp_2lap_g4_one_data_suffix.py`
- Independent verifier SHA-256:
  `2afad43808f53b60b9b3e786642383812741e52520ca581585d41ef1c85aa860`
- Verification artifact:
  `explorations/riffle_dp_2lap_g4_one_data_suffix_bound_verification.json`
- Verification artifact SHA-256:
  `50bb5a9054ccdb112cbadbc234fd5230838fa1ff9e66ca4b415e5747a0b2f1f1`
