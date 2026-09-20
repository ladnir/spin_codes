# Continue the packet-width-eight Riffle co-design

You are continuing a long-running cryptographic proof and systems co-design project in:

```text
C:\Users\peter\repo\permute_conv
```

The associated libOTe implementation worktree is:

```text
C:\Users\peter\.codex\worktrees\permute-conv-code\libOTe
```

The final goal is a construction close to the original Riffle encoder, with packet width eight, that satisfies both requirements:

1. a complete end-to-end first-moment certificate with at least 40 bits of security;
2. an efficient transposed implementation for the PCG application.

Much progress has been made, but packet width eight is not proved. Do not describe a diagnostic bound as a certificate.

## Start from the current checkpoint

The main repository is on branch `codex/g8-proof`. The current committed checkpoint is:

```text
db7ab5e  Freeze repaired g4 global-lane certificate
```

Run these read-only checks first:

```powershell
git status --short
python scripts\verify_g4_savepoint.py
```

The worktree is intentionally dirty. Preserve all existing changes and untracked files. Do not reset, clean, or overwrite them.

Read these files before proposing construction or proof changes:

- `PROOF_STATUS.md`
- `G8_START_HERE.md`
- `G4_GLOBAL_LANE_CERTIFICATE_MANIFEST.json`
- `explorations/riffle_group_chain_proof.md`
- `explorations/conditioned_row_sloped_layout_audit.md`
- `explorations/g4_global_lane_construction_spec.md`
- `explorations/g4_construction_proof_cost_audit.md`
- `explorations/riffle_current_g4_performance_damage.md`
- `explorations/riffle_packet8_implementation_lessons.md`
- `explorations/riffle_packet8_vs_tungsten.md`

The final two files are currently untracked. Treat them as active project records.

## Fixed application and code parameters

The PCG evaluates the transposed binary map on field elements:

```text
G * c
```

where `c` has length `N` and each entry is one element of `GF(2^128)`. One logical element occupies 16 bytes. Cost therefore means field-element XORs, address generation, cache traffic, and memory passes. Packed-bit operations and popcount are irrelevant to the online encoder.

Keep these parameters unless the user explicitly approves a construction change:

```text
N = 2^21
K = 2^20
d = floor(0.09 N) = 188743
inner code = binary systematic extended BCH [128,64,22]
outer code = binary systematic extended BCH [128,64,22]
target packet width g = 8
```

Never replace either BCH code with a smaller code merely because the packet width is eight. Packet width changes the permutation atom, not the BCH dimensions.

For `g=8`, there are `M=N/8=262144` packets and nine profile classes. The complete profile domain has

```text
binom(262144+8,8)-2099
=553169839211945865258921061892182603726
```

feasible profiles. The proof target is

```text
E[Z_d] <= 2^-40,
```

where `Z_d` counts nonzero messages whose encoded word has weight at most `d`.

## Certified results and current limits

The current certified checkpoints are:

- `g=2`: outward bound at most `2^-62.7194713852`;
- `g=4`: construction-bound global-lane outward bound at most `2^-61.7881515553`.

The `g=4` certificate has 21.788 bits of margin beyond the 40-bit target. It uses one globally sampled puncture lane. The global-lane rule is not optional metadata: it repairs a real soundness defect in the earlier independent-lane conditioned-row argument.

The latest rigorous `g=8` factorized receipt checks 197 full-support leaves, 126 independently hardened components, and 8,687 vertex inequalities. Its expanded upper endpoint is approximately:

```text
+580922.1203803732 bits.
```

Thus the rigorous ledger misses the `-40` target by about 580,962 bits. Interval width is below `0.0001` bit, so outward rounding is not the problem.

Additional three-band BL2 and total-spectrum rows reduce later binary64 diagnostics to approximately `+486266` bits. They do not close the proof. More outer-column generation, geometry-only prefix refinement, and the tested nonclumping/polymer route all failed their continuation gates.

## The most important prior mistake

Proof work gradually added lane maps, state permutations, coordinate maps, outer bands, puncture rules, and graph handling. Those additions changed the construction while the project still referred to it as frozen.

Do not repeat that process.

For every proposed change:

1. state the exact construction delta;
2. identify the changed binary matrix or randomness law;
3. state which proof obstruction the change addresses;
4. estimate its online PCG cost;
5. implement and benchmark it as a candidate;
6. obtain user approval before treating it as the construction;
7. bind the approved construction into the certificate manifest.

A proof for one ensemble does not certify a nearby encoder. An implementation optimization may preserve the exact matrix, but a new permutation or mixer usually does not.

## The conditioned-row layout defect

The one-conditioned-row outer was decisive for `g=4`, but its original use assumed a global matching across all three sloped bands. Independent puncture-lane choices can produce punctured blocks that collide in another band. No matching then contains all punctures.

The repaired `g=4` construction samples one global lane and punctures that lane in 128 selected tiles. This produces a perfect matching in every band and preserves each block's puncture marginal. It creates lane-wise correlation, so arguments that need independent per-tile lanes do not transfer.

Any packet-width-eight conditioned-row argument must do one of the following:

- use the certified global-lane rule;
- prove a different valid matching mechanism;
- use a replacement bound that does not assume the missing matching.

Do not reuse the old exact-graph conditioned-row receipt under independent sloped punctures.

## Performance baseline and lessons

The intended fast packet-width-eight Riffle path retains both `[128,64,22]` codes. Its best standalone checkpoint is approximately:

```text
inner:             7.768 ms
outer:             3.636 ms
complete Riffle:  11.670 ms
Tungsten:          9.767 ms
```

The current proof-oriented `g=4` implementation is much heavier:

```text
complete proof-oriented g=4 encoder: 26.001 ms
historical group-chain encoder:      11.085 ms
```

The proof-oriented result shows existence for its declared ensemble. It is not an acceptable deployment-cost target.

The key systems insight is that memory cost is nonlinear in address distance. A regular streaming pass can be cheaper than a smaller irregular gather. A random gather among elements already resident in L1 or L2 can be inexpensive. A strided gather across the full 32 MiB vector can be expensive because it defeats prefetching and increases TLB pressure.

Therefore, do not count “permutations” or “passes” as abstract operations. For each design, record:

- locality scope;
- access order;
- transferred bytes;
- schedule bytes read online;
- reuse distance;
- expected cache level;
- instruction and XOR count;
- end-to-end paired timing.

Tungsten's main lesson is to stream the primary index and keep irregular dependencies inside a small resident window. The optimized Riffle code similarly uses explicit packet pointers, write-intent prefetch, direct outer reads, direct output writes, and graph fusion.

## Local-randomization measurements

The current local-randomization benchmark is:

```text
C:\Users\peter\.codex\worktrees\permute-conv-code\libOTe\libOTe_Tests\RifflePacket8RandomLocal_Bench.cpp
```

The packet-chain setup folding is in:

```text
C:\Users\peter\.codex\worktrees\permute-conv-code\libOTe\libOTe\Tools\RiffleCode\RifflePacketInnerChain.h
```

Current 31-trial paired deltas are:

| Candidate transform | Added complete time |
| --- | ---: |
| random packet `S8`, folded into the existing setup order | about 0.02 ms |
| fused fixed 8-by-8 transpose | about 0.47 ms |
| one independent invertible shear per packet | about 0.24 ms |
| three independent invertible shears per packet | about 0.71 ms |
| independent within-packet `S8` | about 1.21 ms |
| structured node `S8 x S8` | about 1.30 ms |
| full local `S64` | about 2.21 ms |

An elementary shear updates one packet coordinate as

```text
x[d] ^= x[s],  d != s.
```

Every shear is an invertible binary matrix. The independent per-packet variants cost essentially the same as sharing one matrix across the node. The independent variants therefore offer more proof randomness at almost no extra online cost.

Other implementation lessons are:

- a 1 KiB scratch gather beats a data-dependent swap chain;
- explicit streamed `S8` rows beat a compressed `AGL(3,2)` lookup table;
- packet relabeling is free only when setup composes it into the existing order table;
- setup samples permutations with rejection sampling, not biased modulo reduction;
- the benchmark must pin a suitable performance core because this CPU has fast and slow logical cores.

An early benchmark harness passed a 32 MiB vector by value to candidate lambdas. That copy falsely added 6--8 ms to every candidate. The current harness passes a span, restores one physical work buffer outside each timed interval, alternates baseline and candidate order, and checks an identity control. Do not cite results from before this repair.

Never run two benchmarks simultaneously.

## Recommended next experiment

Do not begin with another large proof search. First measure combined low-cost candidates in the real optimized hot path.

At minimum, implement and compare:

1. folded packet `S8` plus one independent shear per packet;
2. folded packet `S8` plus three independent shears per packet;
3. folded packet `S8`, one shear, and the fused 8-by-8 transpose;
4. the structured `S8 x S8` node permutation as the stronger reference.

Use a specialized implementation for each candidate. Do not route all designs through one generic permutation abstraction. Run the identity gate before the measurements. Then run the same candidates at the actual PCG/libOTe call site if available.

The performance output should include:

- exact transform and sampling law;
- whether setup absorbs any schedule;
- persistent and streamed metadata bytes;
- inner and complete median time;
- paired delta from the frozen packet-width-eight path;
- exactness test when the candidate claims to preserve the matrix;
- a clear statement when the candidate defines a new matrix.

## Recommended proof strategy

After the combined cost curve is stable, return to the mathematical obstruction.

First, try to prove the fast packet-width-eight Riffle construction as it currently exists. Record exactly where the proof fails. Do not infer a construction change merely because the present proof technique is weak.

Second, evaluate each lightweight candidate against the failed proof obligation. The main concern has been clumping: the proof can place too much mass in packets, bands, recursive state coordinates, or punctured rows because the available randomness does not destroy the relevant alignment.

For each candidate, answer:

1. Which clumping relation does the transform randomize?
2. Which relations remain invariant?
3. Does the transform restore a clean orbit, matching, occupancy, or transfer lemma?
4. Can that lemma be evaluated first on the worst frozen `g=8` leaves?
5. Is the expected gain remotely comparable to the remaining proof gap?

Prefer a small theorem-facing diagnostic over immediate global coverage. A useful first output is a table over the current worst leaves showing the old bound, the new bound, and the exact source of improvement.

Only after one candidate shows decisive mathematical leverage should you build a complete support-shard or BSP certificate. The final certificate must independently reconstruct the approved construction, harden every selected witness, verify exact ownership and profile counts, and perform outward aggregation.

## Working discipline

- Performance is a first-class constraint.
- Never run two benchmarks or long verification jobs simultaneously.
- Preserve fixed-width and explicitly unrolled hot kernels.
- Avoid dynamic polymorphism, hidden allocation, and generic callbacks in hot paths.
- Use exact arithmetic or outward intervals for theorem-facing claims.
- Use binary64 only for discovery and label it as diagnostic.
- Do not silently edit the construction to make a proof close.
- Do not promote a local witness, sampled profile, or incomplete region to global coverage.
- Keep the old implementation, candidate implementation, and certified construction distinct.
- Provide short status updates during long work.
- At the end of each turn, recommend the next bounded action.

## Expected deliverable from the next session

Produce a concrete co-design decision, not merely another collection of experiments.

The desired deliverable contains:

1. a validated benchmark table for the combined lightweight candidates;
2. exact pseudocode and sampling laws for the candidates that remain viable;
3. a proof audit of fast Riffle packet width eight as-is;
4. a mapping from each candidate change to the obstruction it neutralizes;
5. a recommendation for the cheapest construction worth proving;
6. a bounded proof experiment on the current worst `g=8` leaves;
7. an explicit statement of what remains before an end-to-end `2^-40` certificate.

The target is not an arbitrarily randomized construction that happens to be provable. The target is a fast encoder close to Riffle, with packet width eight, whose exact deployed ensemble is the ensemble certified by the proof.

## 2026-08-17 Double-Parity investigation

The fast sequential construction is now called **Riffle g=8**.
The equal-length structured candidate is called **Riffle Double-Parity g=8**.
The candidate uses two field-valued parity blocks and has block distance three.

The paired implementation benchmark was favorable: Double-Parity reduced complete runtime by 3.5% to 4.9%.
The proof audit was negative.

For the sequential outer layout, each BCH word occupies at most 16 eight-bit packets.
The Double-Parity outer is an MDS $[16{,}386,16{,}384,3]$ code over $\mathbb F_{2^{64}}$.
Its exact block-weight enumerator and the terminal-placement event give

```text
log2 lower bound from block weight 3  = -63.407484138361
log2 lower bound from block weight 4  = -43.043708373439
log2 lower bound from block weight 5  = -23.016236846907
log2 lower bound from block weight 6  =  -3.266186470799
```

Thus Riffle Double-Parity g=8 cannot satisfy `E[Z_d] <= 2^-40` in the benchmarked sequential layout.
The obstruction is independent of the chosen distinct parity coefficients and the detailed BCH spectrum.

The same audit found an exact equal-pair obstruction for the fast Graph-24 outer.
For the current power-of-two length, that family contributes at least

```text
2^-35.189306074384
```

to the first moment.
This result supersedes the earlier one-block value near `2^-48.17`, which omitted equal two-block messages.

The exact derivation and scope are in:

- `explorations/riffle_double_parity_proof_audit.md`;
- `scripts/analyze_riffle_double_parity.py`.

Do not start a complete Double-Parity certificate.
The next candidate must change the positive entropy slope of rate-one-half blocks with contiguous eight-bit packets.
Plausible directions are a lower effective outer rate, a costed striped layout, or a second global mixing stage without a common terminal region.
Adding a constant number of parity blocks is insufficient.

## 2026-08-17 Riffle DP g=4 first pass

The sequential packet-width-four candidate is called **Riffle DP g=4**.
It retains the unpunctured Double-Parity outer and changes only the packet
width from eight to four.

The exact maximum-packet MDS shell obstruction now has a favorable slope:

```text
block weight 3   -230.231606973875
block weight 4   -265.504413728495
block weight 5   -301.127813132726
block weight 6   -337.042932694766
```

Thus the terminal-placement family that rejects Riffle DP g=8 does not
reject Riffle DP g=4. This result is not an upper-bound certificate.

An exact enumeration of all 12,587 BCH words on at most eight bytes found
four-bit packet supports from 11 through 16. It proves that every nonzero
local BCH word occupies at least nine four-bit packets. It does not exclude
support nine or ten from words on more than eight bytes.

Sampling 300,000 words in each MDS block-weight-three category produced
terminal estimates between `2^-175.38` and `2^-190.99`. Rare low-support
words remain the proof risk.

The exact local-tail gate subsequently passed. A direct parity-check scan
visited 809,785,132 packet subsets and established:

```text
local nibble support 1..10     0 words
local nibble support 11       20 words
local nibble support 12     1526 words
local nibble support 13    37014 words
```

Thus the exact local nibble distance is 11. The exact minimum outer layer
contains 26 block-weight-three words, all in the one-data category. Its
terminal contribution is `2^-109.961733072748`. The other three categories
contain no words at total support 33.

Exact local counts give rigorous upper bounds of `2^-70.332635254404` for
total support 34 and `2^-67.155898463807` for total support 35. These bounds
do not cover support 36 or greater.

The next gate is a weighted upper bound for the complete block-weight-three
shell. A tested byte-to-nibble envelope gives the useless upper exponent
`+55.19`; the next bound must retain more nibble information. Do not continue
brute-force packet-subset enumeration by default, and do not start a complete
Riffle DP g=4 certificate before this gate passes.

The derivation and script are:

- `explorations/riffle_double_parity_g4_first_pass.md`;
- `scripts/analyze_riffle_double_parity_g4.py`.
- `scripts/enumerate_bch_g4_low_support.cpp`.

## 2026-08-18 concurrent proof and refutation plan

The active Riffle DP g=4 program now has two concurrent tracks. The proof
track builds complete upper bounds. The refutation track searches for exact
lower-bound families and adversarial inner placements. Both tracks use the
same weighted enumerator for triples of local BCH words.

The synchronized gates are:

1. freeze the construction manifest;
2. complete the weight-three terminal audit;
3. bound and attack every inner outcome at outer weight three;
4. cover and attack sparse outer weights;
5. cover and attack dense outer weights;
6. assemble and independently verify the complete first-moment ledger.

A proof requires complete coverage with total at most `2^-40`. A refutation
requires one exact family with contribution above `2^-40`. Failure to find a
counterexample is not a proof.

The detailed plan is:

- `explorations/riffle_dp_g4_proof_refutation_plan.md`.

## 2026-08-18 Riffle DP g=4 G0 and partial G1 receipts

G0 is complete. The frozen proof ensemble samples one uniform bijection of
the 524,352 four-bit packets. It has no lane bijections, local packet
relabeling, striped layout, puncturing, or additional mixer.

The finite-seed `init(seed)` implementation is not the theorem-facing
uniform distribution. The new `initOrder(packet_order)` interface realizes
the matrix for any specified bijection. Symbolic basis checks verify the
local BCH kernel, the packet-four inner kernel, and every coefficient
recurrence. A full-size explicit-order integration check also passes.

The frozen files are:

- `explorations/riffle_dp_g4_construction_manifest.json`;
- `scripts/verify_riffle_dp_g4_manifest.py`;
- `scripts/test_riffle_dp_g4_manifest.cpp`.

The support-34 and support-35 bounds from the first pass were loose. A new
exact classifier gives the complete result through support 35:

```text
support 33   26 words   3 new dirty coefficient classes
support 34    0 words   0 new dirty coefficient classes
support 35   39 words  11 new dirty coefficient classes
```

The 39 support-35 words comprise 36 one-data words and three words with two
data symbols and only the weighted parity active. The other categories are
empty through support 35.

The exact terminal lower family through support 35 contributes
`2^-109.94433796256`. It does not refute the candidate.

Exact one-Krawtchouk Delsarte certificates give an unconditional upper bound
of about `2^-30.101388902225` for the complete block-weight-three terminal
contribution. This bound is 9.8986111 bits above the `2^-40` target, so G1
remains unresolved.

The next proof gate is narrower. If at most 118,378,995 currently clean
coefficient classes contain a word of support 36, 37, or 38, the terminal
weight-three bound falls below `2^-40`. This threshold is rigorous. Search
these three layers directly or prove a stronger coefficient-specific bound.

The current G1 artifacts are:

- `scripts/riffle_dp_g4_weighted_triples.py`;
- `scripts/classify_riffle_dp_g4_low_triples.cpp`;
- `explorations/riffle_dp_g4_g1_low_classification.txt`;
- `explorations/riffle_dp_g4_g1_terminal_ledger.json`;
- `explorations/riffle_dp_g4_g1_support33_refutation.json`;
- `explorations/riffle_dp_g4_g1_delsarte_d36_j42.json`.

## 2026-08-18 Riffle DP g=4 terminal G1 closed

The support-38 classifier closes the terminal block-weight-three gate. An
exact local scan extended the nibble enumerator to:

```text
support 14       742031 words
support 15     13430995 words
```

The exact coefficient-class results are:

```text
minimum support 33                 3 classes
minimum support 35                11 classes
minimum support 36             32779 classes
additive minimum support 37        0 classes
additive minimum support 38    49143 classes
```

The repeated categories are exact through support 36. The proof assigns
every remaining repeated class the conservative support-37 bound. The
additive categories are exact through support 38, so every other additive
class has minimum support at least 39.

Exact Delsarte certificates and exact rational aggregation give

```text
complete terminal weight-three upper bound  2^-40.212996059540
margin below 2^-40                       0.212996059540 bits
```

Thus G1 passes. This is not an end-to-end proof. It covers only the event in
which every active outer packet enters the final 2,949 inner nodes.

The next gate is G2: bound every inner outcome at outer block weight three.
The refutation stream should search nonterminal placements for the newly
identified support-36 and support-38 coefficient classes.

The closing artifacts are:

- `explorations/riffle_dp_g4_g1_terminal_proof.md`;
- `explorations/riffle_dp_g4_g1_terminal_ledger.json`;
- `explorations/riffle_dp_g4_support38_classification.txt`;
- `scripts/classify_riffle_dp_g4_support38.cpp`;
- `scripts/verify_riffle_dp_g4_g1.py`.

## 2026-08-18 Riffle DP g=4 G2 autonomous audit

G2 remains active. The exact fixed recurrence is

```text
V_i = Acc(U_i + S_{i-1})
S_i = P V_i
```

The zero-input state map `T = P Acc` is invertible. Its minimal polynomial is
`0x1b7eebf80814fbc75`. Its irreducible component degrees are
`1,2,4,9,10,18,20`, and its matrix order is `91625532075`.

The minimum exact cycle averages for the individual components range from
28.0 to 33.6 emitted bits per node. An exhaustive audit of maximal combined
component subspaces through dimension 20 replays about 3.8 million states.
For cycles that activate every component in a tested combination, the worst
average is `31.9403714565` bits per node. No low-average algebraic cycle was
found in the covered subspaces.

The independent verifier reports:

```text
status=EXACT_G2_AUTONOMOUS_MAP_VERIFIED
```

The exact one-packet search still has only six below-threshold turnoffs. The
best authenticated lower family contributes `2^-124.233993125178`, so it does
not refute the candidate.

Do not use the tempting two-output BCH-distance lemma. The BCH pair is
`(V_i, P V_i)`, but the next emitted block is `Acc(P V_i)`. Accumulation is not
weight-preserving.

Current G2 artifacts are:

- `explorations/riffle_dp_g4_g2_status.md`;
- `explorations/riffle_dp_g4_g2_autonomous_map.json`;
- `scripts/analyze_riffle_dp_g4_autonomous_map.py`;
- `scripts/verify_riffle_dp_g4_g2_autonomous_map.py`;
- `explorations/riffle_dp_g4_g2_single_packet_orbits.json`;
- `explorations/riffle_dp_g4_g2_turnoff_lower_family.json`.

Next, build a finite-window bound for sparse reachable states and an exact
three-packet turnoff search. Then aggregate both over support-33 through
support-38 placement profiles.

## 2026-08-18 Riffle DP g=4 G2 finite-window tranche

The finite-window proof artifact is complete. Exhaustive forward enumeration
covers all 5,130,659,560 nonzero 64-bit outputs of weight at most eight. An
independent reverse enumeration with `Q^-1` and `Q^-2` confirms:

```text
one-step low returns                         0
two-step low returns                11,039,336
minimum two-step segment weight             19
```

For every nonzero zero-input autonomous segment of length `L`, the rigorous
bound is

```text
delta_aut(L) >= 9*L - 8 - 8*floor((L-1)/3)
```

This exceeds `d=188766` at `L=29807`. At the full chain length 32,772, the
bound is 207,556, with margin 18,790.

The full-chain distinct-node three-packet search tests 1,887,552,000 exact
relations. It finds exactly 27 turnoffs, all by node six. The exact support-33
lower family contributes `2^-134.276933547563`, so it does not refute the
candidate.

New artifacts are:

- `explorations/riffle_dp_g4_g2_autonomous_window_bound.json`;
- `explorations/riffle_dp_g4_g2_low_output_returns_w8_gap2.json`;
- `scripts/analyze_riffle_dp_g4_low_output_returns.cpp`;
- `scripts/verify_riffle_dp_g4_low_output_returns.cpp`;
- `explorations/riffle_dp_g4_g2_three_packet_turnoffs_full.json`;
- `explorations/riffle_dp_g4_g2_three_packet_lower_family.json`;
- `scripts/search_riffle_dp_g4_three_packet_turnoffs.cpp`;
- `scripts/verify_riffle_dp_g4_g2_three_packet.py`.

G2 remains open. The next proof task is to aggregate autonomous gap costs while
charging every zero-state reset to a compatible packet placement. The next
refutation task is to cover same-node packet collisions and four-packet
turnoffs.

## 2026-08-18 Riffle DP g=4 G2 reset-envelope result

The exact geometric reset envelope is now available. It combines the
finite-window lemma with every abstract live/reset pattern and all packet
collisions within nodes. Its zero-reset upper bounds are:

```text
support 33   2^-4.475042948481
support 34   2^-4.609038327823
support 35   2^-4.742912425689
support 36   2^-4.876690069363
support 37   2^-5.010371482619
support 38   2^-5.143956675803
```

These bounds are rigorous but far too loose for G2. Even perfect algebraic
charging of resets cannot repair the zero-reset row. A proof must exploit the
distribution of reached states, packet values, or nibble slots under the
uniform packet permutation.

The full-chain collision search covers initial and reset nodes of packet
support one or two. It finds:

```text
1 -> 1       6 turnoffs
1 -> 2      24 turnoffs
2 -> 1      64 turnoffs
2 -> 2     445 turnoffs
```

All 539 turnoffs occur by node eight. The strongest authenticated collision
lower family is `2^-123.048186206738`. Repeated two-packet episodes give
`2^-123.049568554040`. Neither result refutes the candidate.

Additional artifacts are:

- `explorations/riffle_dp_g4_g2_reset_pattern_placements.json`;
- `scripts/bound_riffle_dp_g4_reset_pattern_placements.py`;
- `explorations/riffle_dp_g4_g2_multi_packet_node_turnoffs.json`;
- `scripts/search_riffle_dp_g4_multi_packet_node_turnoffs.cpp`;
- `explorations/riffle_dp_g4_g2_collision_turnoff_lower_family.json`;
- `scripts/analyze_riffle_dp_g4_collision_turnoff_family.py`;
- `scripts/verify_riffle_dp_g4_g2_collision_turnoffs.py`;
- `explorations/riffle_dp_g4_g2_multi_turnoff_lower_family.json`.

Next proof direction: build a fixed-map transfer or exponential-moment bound
that retains nibble slots and values. Next refutation direction: search reset
episodes with intermediate active nodes or at least three packets in one node.

## 2026-08-18 Riffle DP g=4 reached-state mixing tranche

Exact first-node orbit enumeration now covers packet supports one through
three:

```text
support 1    240 drives      latest crossing 5926
support 2  27000 drives      latest crossing 5943
support 3 1890000 drives     latest crossing 5945
```

Using those exact prefixes improves the support-33 zero-reset envelope from
`2^-4.4750` to `2^-6.767496031893`. The one-packet first node still dominates,
so deterministic gap bounds remain insufficient.

The observability probe found an exact arbitrary-state witness:

```text
state  = 0x6331abf1b618efab
length = 5956
weight = 188730 = d - 36
```

Its preimage has nibble support 13 and is incompatible with all authenticated
support-33 value multisets. It blocks an overly strong universal window lemma
but does not refute the construction.

Exact component Walsh spectra are encouraging. Individual maximum biases are
`0.625, 0.25, 0.133365, 0.060753, 0.001121, 0.004360, 0.005245` for component
degrees `1,2,4,9,10,18,20`. The spectra are exact iid diagnostics. The former
without-replacement product bound is invalid and must not be cited.

A full-character hill search finds no bias above the known degree-one value
`0.625`. This is diagnostic rather than exhaustive.

An exact authenticated three-node cluster now gives a stronger refutation
fact. Its weight is 188,505 over the final 5,940 nodes, or 261 below `d`.
Every later terminal translation is also bad, giving 5,938 exact placements.
The independently counted single-family probability is
`2^-555.961662409912`. Thus, deterministic distance is false, but the
random-permutation claim is not refuted.

New artifacts are:

- `explorations/riffle_dp_g4_g2_initial_collision_orbits_s2.json`;
- `explorations/riffle_dp_g4_g2_initial_collision_orbits_s3.json`;
- `scripts/analyze_riffle_dp_g4_initial_collision_orbits.cpp`;
- `explorations/riffle_dp_g4_g2_first_impulse_envelope.json`;
- `scripts/bound_riffle_dp_g4_first_impulse_envelope.py`;
- `explorations/riffle_dp_g4_g2_observability_probe.json`;
- `scripts/probe_riffle_dp_g4_observability_distance.py`;
- `explorations/riffle_dp_g4_g2_component_mixing_probe.json`;
- `scripts/probe_riffle_dp_g4_component_mixing.py`;
- `explorations/riffle_dp_g4_g2_support33_component_mixing_probe.json`;
- `scripts/probe_riffle_dp_g4_support33_component_mixing.py`;
- `explorations/riffle_dp_g4_g2_full_character_bias_probe.json`;
- `scripts/probe_riffle_dp_g4_full_character_bias.py`;
- `explorations/riffle_dp_g4_g2_support33_cluster_probe_l5940.json`;
- `explorations/riffle_dp_g4_g2_support33_cluster_family.json`;
- `scripts/probe_riffle_dp_g4_support33_cluster.py`;
- `scripts/verify_riffle_dp_g4_support33_cluster.py`.

The next proof step is to certify all full-state Fourier characters and their
aggregate mass, then connect reached-state mixing to output weight. The
matching refutation step is to count the basin around the exact cluster and
search broader sparse episodes.

## 2026-08-18 Riffle DP g=4 broad-suffix tranche

The exact single-swap basin around the three-node cluster contains 919
assignments and 5,419,150 disjoint placements. Its aggregate contribution is
only `2^-544.542826136930`.

A broader diagnostic conditions on all packets occupying the final 5,940
nodes. For the authenticated support-33 profile, 79,583 of 100,000 placements
have weight at most `d`. The estimated three-family contribution is
`2^-80.0607`. Scans through 8,000 nodes produce similar unconditional
estimates.

An exact one-data enumeration checks 631,767,040 pairs `(i,t)` and produces
85,828 outer words through support 39. The exact support counts are
`26, 36, 3233, 510, 933, 81090` for supports
`33, 35, 36, 37, 38, 39`. An independent verifier reproduces the record
stream. Monte Carlo over this exact population estimates an aggregate suffix
contribution of `2^-76.01295`.

New artifacts are:

- `explorations/riffle_dp_g4_g2_support33_cluster_basin.json`;
- `scripts/analyze_riffle_dp_g4_support33_cluster_basin.py`;
- `explorations/riffle_dp_g4_g2_support33_suffix_monte_carlo.json`;
- `scripts/probe_riffle_dp_g4_support33_suffix_monte_carlo.cpp`;
- `explorations/riffle_dp_g4_one_data_support39.json`;
- `explorations/riffle_dp_g4_one_data_support39.bin`;
- `scripts/enumerate_riffle_dp_g4_one_data_support39.cpp`;
- `scripts/verify_riffle_dp_g4_one_data_support39.cpp`;
- `explorations/riffle_dp_g4_g2_one_data_suffix_monte_carlo.json`;
- `scripts/probe_riffle_dp_g4_one_data_suffix_monte_carlo.cpp`.

The next refutation task is to extend the exact population count to pair-data
and three-data coefficient classes. The next proof task is a recurrence-aware
certificate for the 12-node full-character observability code; generic MILP
and XOR-SAT formulations did not close its 35-bit cutoff.

## 2026-08-19 Riffle DP-2Lap g=4 tranche

**Riffle DP-2Lap g=4** is a separately named candidate. The frozen
**Riffle DP g=4** construction remains unchanged. The new candidate retains
the terminal state, wraps to the first logical node, and applies the same
deterministic convolution to the overwritten first-lap word. It samples no
new permutation. Its exact linear map is

\[
G(x)=F(F(x))\mathbin\oplus J(L(x)).
\]

An exhaustive certificate proves that every nonzero three-node autonomous
block on the second lap has weight at least 25. The certificate enumerates all
5,130,659,560 nonzero outputs of weight at most eight. An independent verifier
repeats the enumeration. Therefore, every nonzero first-lap terminal state
contributes at least 188,775 bits in a zero-input prefix of 22,653 nodes.

The retained-state problem has consequently collapsed for this geometry.
For inputs contained in the final 10,119 nodes, the only low-growth terminal
state is zero. The remaining proof obligation in this stratum is
\(\Pr[L(X)=0]\), although the complete suffix-placement event can sometimes
be budgeted without estimating that probability.

For the exact one-data population through support 39, charging the complete
final-10,119-node placement event gives the rigorous aggregate upper bound

```text
2^-48.293011187781 <= bound <= 2^-48.293011187780.
```

This closes that placement row for all 85,828 enumerated one-data words. It
does not cover their other placements.

On the refutation side, 17 known three-packet episodes and 342 known collision
episodes reset both laps. They reduce to 88 packet-value multisets. An exact
search found no partition of any of the 14 authenticated support-33 packet
multisets into these episode types. An independent additive-closure verifier
confirmed the negative result. The smallest missing primitive has six packets,
not five. The only minimum residual profiles are `10^3 13^3` and `5^3 13^3`.

An exact 360,743,168-step search tested every two-node split for both
six-packet profiles. The profile `5^3 13^3` has two double turnoffs in the 4+2
split. Each has two-lap weight 29.

One six-packet turnoff combines with eight known episodes to form an
authenticated support-33 input of weight 56. Hence deterministic distance is
false. A disjoint gap-placement family has probability about `2^-473.407`, so
this exact family does not refute the random-permutation claim.

The current DP-2Lap artifacts are:

- `explorations/riffle_dp_2lap_g4_initial_diagnostic.md`;
- `explorations/riffle_dp_2lap_g4_three_node_block_certificate_w8.json`;
- `scripts/certify_riffle_dp_2lap_g4_three_node_block.cpp`;
- `scripts/verify_riffle_dp_2lap_g4_three_node_block.cpp`;
- `explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_search.json`;
- `scripts/search_riffle_dp_2lap_g4_authenticated_double_turnoffs.py`;
- `explorations/riffle_dp_2lap_g4_authenticated_double_turnoff_verification.json`;
- `scripts/verify_riffle_dp_2lap_g4_authenticated_double_turnoffs.py`;
- `explorations/riffle_dp_2lap_g4_one_data_suffix_bound.json`;
- `scripts/bound_riffle_dp_2lap_g4_one_data_suffix.py`;
- `explorations/riffle_dp_2lap_g4_one_data_suffix_bound_verification.json`;
- `scripts/verify_riffle_dp_2lap_g4_one_data_suffix.py`;
- `explorations/riffle_dp_2lap_g4_six_packet_two_node_turnoffs.json`;
- `scripts/search_riffle_dp_2lap_g4_six_packet_two_node_turnoffs.cpp`;
- `explorations/riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json`;
- `scripts/construct_riffle_dp_2lap_g4_authenticated_six_packet_turnoff.py`;
- `explorations/riffle_dp_2lap_g4_authenticated_six_packet_turnoff_verification.json`;
- `scripts/verify_riffle_dp_2lap_g4_authenticated_six_packet_turnoff.py`.

Next proof direction: cover the complement of the final-10,119-node stratum.
Use terminal-zero mixing only where the geometric placement probability does
not already fit the ledger. Next refutation direction: search primitive
six-packet double turnoffs across three or more active nodes, then search
cancellations between components that do not reset separately.

## 2026-08-20 global-permutation Fourier tranche

The terminal-zero proof now uses the global packet permutation explicitly.
Sample labeled cells iid and condition on distinctness. This produces the
uniform ordered injection and permits exact iid Fourier factorization.

An audit invalidated the earlier adaptive conditional-bias multiplication.
For an exactly balanced character, the true two-draw expectation has magnitude
`1/524351`; the retired product gives `1/274943971201`. The old support-33
component-mixing receipt is labeled `INVALIDATED_BOUND`.

Each packet-value contribution map from 524,352 cells into 64-bit terminal
states is injective. Parseval therefore pays for two repeated packet values.
Certifying all diagnostic full-character maxima would give an aggregate
support-33 terminal-zero value near `2^-59.2956`.

A much weaker certificate suffices:

```text
value 15: full-character bias <= 5/8
values 1..14: full-character bias <= 1/2
```

These caps already give a conditional target between
`2^-44.464429949134` and `2^-44.464429949133` across all 26 support-33 words.

Since `32772 = 12 * 2731`, the caps would follow from two finite code-distance
claims. The 192-bit, 64-dimensional 12-node character code would need
two-sided distance 36 for value 15 and 48 for every other value. The value-15
bound is tight.
The degree-one character `0xa685aac60e5acfb3` has block weight 36.

A corrected, fully rigorous degree-20 quotient calculation gives only
`2^-15.298107519745` after aggregating the 26 words. It validates the
conditioning argument but cannot close the 40-bit target. Generic MILP and
plain CNF encodings do not close the tight local distance cutoff.

An exact partial certificate exhausts every character whose irreducible-
component support has total dimension at most 22. All 274,291,320
character/value evaluations satisfy the proposed two-sided distance bounds.
The individual degree-18 and degree-20 components both have minimum
two-sided distance 58. The unresolved characters span component sets above
dimension 22.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_global_permutation_proof_route.md`;
- `explorations/riffle_dp_2lap_g4_terminal_zero_parseval_gate.json`;
- `scripts/analyze_riffle_dp_2lap_g4_terminal_zero_parseval_gate.py`;
- `explorations/riffle_dp_2lap_g4_support33_projection_bound.json`;
- `scripts/bound_riffle_dp_2lap_g4_support33_projection.py`;
- `explorations/riffle_dp_g4_without_replacement_bias_audit.json`;
- `scripts/audit_riffle_dp_g4_without_replacement_bias.py`;
- `explorations/riffle_dp_2lap_g4_local_character_components.json`;
- `scripts/certify_riffle_dp_2lap_g4_local_character_components.py`.

The targeted 38-dimensional \(C_{18}\oplus C_{20}\) search refutes the
uniform 12-node lemma. For value 7, character
`0x852ac8fcc6b27c67` gives a 192-bit word of weight 146 and therefore
two-sided weight 46, below the required 48. An independent implementation
confirms that this is a genuinely mixed degree-18/degree-20 character and
replays the exact word.

The same witness has full-orbit bias only `1894/524352`, and exactly one of
2,731 consecutive 12-node blocks violates the local threshold. It therefore
refutes the per-block proof device, not Riffle DP-2Lap g=4 or the desired
global character cap. Its first two blocks have combined two-sided weight
148 out of 384. The next proof task is a transition-aware or 24-node
amortized certificate with exact accounting for the odd final 12-node block.
Continue the complementary placement ledger and broader cancellation search
in parallel.

New exact artifacts are:

- `explorations/riffle_dp_2lap_g4_c18_c20_counterexample.md`;
- `explorations/riffle_dp_2lap_g4_c18_c20_v07.json`;
- `scripts/certify_riffle_dp_2lap_g4_c18_c20.cpp`;
- `explorations/riffle_dp_2lap_g4_c18_c20_counterexample_verification.json`;
- `scripts/verify_riffle_dp_2lap_g4_c18_c20_counterexample.py`.

The proposed 24-node repair is now exact for
\(C_{18}\oplus C_{20}\). Ten disjoint information sets reduce the required
distance searches to radius seven for value 15 and radius nine for every
other value. The primary and independent implementations each exhaust
32,062,999,280 information vectors. Both pass all 15 packet values.

The certificate proves 24-node two-sided distance 73 for value 15 and 97 for
all other nonzero values. Across 1,365 complete pairs, these bounds contribute
99,645 and 132,405, respectively. They exceed the global requirements 98,316
and 131,088 without using the last 12 nodes. Thus, the desired full-orbit
character caps are proved for this mixed component subcode.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_c18_c20_24node_certificate.md`;
- `explorations/riffle_dp_2lap_g4_c18_c20_24node_audit.json`;
- `scripts/certify_riffle_dp_2lap_g4_c18_c20_24node.cpp`;
- `scripts/verify_riffle_dp_2lap_g4_c18_c20_24node.cpp`;
- `scripts/audit_riffle_dp_2lap_g4_c18_c20_24node.py`;
- the 15 primary and 15 independent per-value JSON receipts.

An ambitious extension targeted all 57 component supports of dimension at
most 30 through eight maximal subcodes. The stronger value-15 distance-73
claim is false. In support `(1,9,20)`, the known pure degree-one character
`0xa685aac60e5acfb3` has 24-node weight 72.

Independent replay shows that every 12-node block has weight 36 and every
complete 24-node block has weight 72. The full-orbit weight is
`1365*72+36 = 98316`, so the global bias is exactly `5/8`. This refutes only
the attempt to discard the endpoint. It does not refute the desired global
cap.

The corrected asymmetric proof target is:

```text
value 15: 24-node distance >= 72 and final-12 distance >= 36
values 1..14: 24-node distance >= 97
```

The run checked 2,917,959,325 information vectors before finding the witness.
Support `(10,20)` passed every value. Support `(1,9,20)` passed values 1
through 14. The audited ledger records 30 checked receipts, three fully
certified strong-threshold supports, 28 value-15 refutations by containment,
and the unrun remainder.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component24_dimension30_counterexample.md`;
- `explorations/riffle_dp_2lap_g4_component24_dimension30_counterexample_ledger.json`;
- `explorations/riffle_dp_2lap_g4_component24_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component24_s49_v15_counterexample_verification.json`;
- `scripts/certify_riffle_dp_2lap_g4_component_24node.cpp`;
- `scripts/run_riffle_dp_2lap_g4_component_24node.py`;
- `scripts/verify_riffle_dp_2lap_g4_component24_s49_v15_counterexample.py`;
- `scripts/audit_riffle_dp_2lap_g4_component24_dimension30_counterexample.py`.

The corrected asymmetric rerun is complete. Eight maximal subcodes cover all
57 nonempty component supports of total dimension at most 30. For every such
support, exact primary and independent searches prove:

```text
value 15: 24-node distance >= 72 and final-12 distance >= 36
values 1..14: 24-node distance >= 97
```

Each implementation passes all 128 cases and checks 6,425,612,528
information vectors. The primary path uses fixed-dimension split-half tables.
The independent path changes the partitions and uses direct recursive
fixed-weight enumeration. A Python audit reconstructs the code columns,
checks all 256 receipts and information-set ranks, proves the eight-subcode
cover of the 57 targets, and replays the tight degree-one witness.

The local result implies the required full-orbit character caps on every
covered support. Values 1 through 14 obtain weight at least
`1365*97 = 132405 > 131088`. Value 15 obtains exactly the required lower
bound `1365*72+36 = 98316`. The degree-one witness attains this equality and
global bias `5/8`.

Together with the earlier exact `(18,20)` certificate, the current finite
results cover 58 of the 127 nonempty irreducible-component supports. The
remaining 69 supports are the next local proof obligation.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component_endpoint_dimension30_certificate.md`;
- `explorations/riffle_dp_2lap_g4_component_endpoint_dimension30_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_endpoint_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_endpoint_independent_run.json`;
- the 128 primary and 128 independent per-case JSON receipts;
- `scripts/audit_riffle_dp_2lap_g4_component_endpoint_dimension30.py`;
- generalized source `scripts/certify_riffle_dp_2lap_g4_component_24node.cpp`;
- sequential runner `scripts/run_riffle_dp_2lap_g4_component_24node.py`.

The four dimension-31 frontier supports are now certified. Masks `0x2c`,
`0x33`, `0x4a`, and `0x51` pass all corrected endpoint-aware bounds. The
primary and independent implementations each pass 64 cases and check exactly
7,716,616,224 information vectors. The independent implementation uses
different partitions and direct recursive enumeration.

A separate Python reconstruction audits all 128 full receipts and the 16
primary obstruction-probe receipts. It checks the component kernels,
observation columns, information-set ranks, receipt hashes, vector counts,
and the exact four-support frontier. The complete tranche now contains all 61
supports through dimension 31. The separate `(18,20)` result raises merged
coverage to 62 of 127 supports, leaving 65 open.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component_dimension31_certificate.md`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_independent_run.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension31.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension31.py`.

The dimension-32 frontier is now certified. Masks `0x2d`, `0x34`, `0x4b`,
and `0x52` pass every corrected endpoint-aware bound. The primary and
independent implementations each pass 64 cases and check exactly
10,119,775,656 information vectors. The independent path uses different
partitions and recursive enumeration.

The Python reconstruction audits all 128 full receipts and 16 primary probe
receipts. It checks the component kernels, observation columns,
information-set ranks, receipt and executable hashes, vector counts, and the
exact dimension-32 frontier. The complete tranche now contains all 65
supports through dimension 32. The separate `(18,20)` result raises merged
coverage to 66 of 127 supports, leaving 61 open.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component_dimension32_certificate.md`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_independent_run.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension32.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension32.py`.

The next finite frontier is dimension 33: `0x2e=(2,4,9,18)`,
`0x35=(1,4,10,18)`, `0x4c=(4,9,20)`, and `0x53=(1,2,10,20)`. Its predicted
workload is 12,216,115,184 vectors per implementation.

The bounded structural anti-cancellation checkpoint is complete. It gives
two exact negative results and one exact reparameterization:

- Five consecutive nodes determine the 64-bit character, but the kernels of
  four-node observation maps are not confined to small certified components.
  They include component supports of dimensions 61, 63, and 64.
- The complementary dimension-32 support pairs are not orthogonal on either
  the 12-node or 24-node window. Every tested cross-Gram matrix is nonzero,
  with rank between 30 and 32.
- For a fixed character, all 15 packet-value character sums are the
  nontrivial Walsh coefficients of one nonnegative 16-bin coefficient
  histogram. The audit verifies this representation and its Parseval
  identity exactly.

A full 64-dimensional refutation search used 10,000 pseudorandom restarts per
case, together with deterministic basis and degree-one seeds. It found no
violation. The smallest non-15 two-sided weight found was 120. Value 15
returned the known tight weights 72 on 24 nodes and 36 on 12 nodes. This is
diagnostic evidence, not a certificate.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_structural_anticancellation_checkpoint.md`;
- `explorations/riffle_dp_2lap_g4_structural_checkpoint.json`;
- `explorations/riffle_dp_2lap_g4_full_local_distance_probe.json`;
- `scripts/analyze_riffle_dp_2lap_g4_structural_checkpoint.py`;
- `scripts/probe_riffle_dp_2lap_g4_full_local_distance.py`.

Recommended next goal: characterize or safely relax the attainable 16-bin
coefficient histograms. The required Walsh bounds are
\(|B_v|\le190\) for \(v\ne15\) and \(|B_{15}|\le240\) on 24 nodes, plus
\(|B_{15}|\le120\) on the 12-node endpoint. Preserve the degree-one
histogram, which attains both value-15 equalities. If this bounded histogram
step yields no useful restriction, resume exact certification with the four
dimension-33 supports.

The histogram-transition checkpoint is now complete. It gives the exact
ordered-state identity behind the coefficient histogram. If \(T\) is the
transpose recurrence, then the coefficient word at node \(t\) and slot \(s\)
is exactly nibble \(s\) of \(T^t\chi\). Thus, the 16-bin histogram counts
nibble values along 24 ordered transpose states.

The histogram alone is not a transfer state. The exact witness states
`0x7373d52af3615f19` and `0x953f133d5f76721a` have the same current nibble
histogram and different successor histograms. Any transfer proof must retain
ordered information or a stronger state quotient.

The first uniform local redundancy occurs at five nodes. Exact MacWilliams
spectra for all 15 `[80,64]` five-node codes show two-sided distance only 2
or 3. A relaxation that keeps only all sliding five-node distance constraints
proves a 24-node two-sided bound of only 8--12. Moreover, the exact number of
five-node words with two-sided weight at most 24 is approximately
`2^52.88706878`. Direct enumeration of the low-window frontier is therefore
not practical. For value 15 at threshold 17, the count is still approximately
`2^41.92674088`.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_histogram_transition_checkpoint.md`;
- `explorations/riffle_dp_2lap_g4_histogram_transition.json`;
- `scripts/analyze_riffle_dp_2lap_g4_histogram_transition.py`.

Recommended next goal: resume exact certification with the dimension-33
frontier. Treat a structural replacement as open research until a quotient
retains enough nibble ordering to improve materially on the five-node bound.

The dimension-33 frontier is now certified. Masks `0x2e`, `0x35`, `0x4c`,
and `0x53` pass every corrected endpoint-aware bound. The primary and
independent implementations each pass 64 cases and check exactly
12,216,115,184 information vectors. The independent path uses different
partitions and recursive enumeration.

The Python reconstruction audits all 128 full receipts and 16 primary
preflight receipts. It checks the component kernels, observation columns,
information-set ranks, receipt and executable hashes, vector counts, and the
exact dimension-33 frontier. The endpoint-aware tranche now contains all 69
supports through dimension 33. The separate `(18,20)` result raises merged
coverage to 70 of 127 supports, leaving 57 open.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component_dimension33_certificate.md`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_independent_run.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension33.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension33.py`.

The next finite frontier is dimension 34: `0x2f=(1,2,4,9,18)`,
`0x36=(2,4,10,18)`, `0x4d=(1,4,9,20)`, and `0x54=(4,10,20)`. Its predicted
workload is 15,745,416,320 vectors per implementation. Dimension 35 raises
the information radius and is substantially more expensive. Recommended
next goal: certify dimension 34, then reassess before the dimension-35 cost
jump.

The completion-strategy reassessment rejects a continued frontier sweep.
The remaining exact information-set workload grows to about 60 quadrillion
vectors per implementation for full support `0x7f`. Certifying every
remaining frontier extrapolates to centuries at the dimension-33 rates.

An exact component-split prototype then tested the direct full-code route.
For support `0x53`, packet value 1, the split `0x40+0x13` has dimensions
20+13. The optimized 16-thread engine enumerated all `2^33-1` nonzero pairs
and proved exact two-sided minimum 120. A Python reconstruction verifies the
generator words, pair count, witness, executable, and hashes.

The same raw method does not scale. A balanced 32+32 full-code split requires
`2^64-1` pair comparisons. The committed dimension-33 rate extrapolates to
about 335 years for one full-code packet-value case. The balanced right list
also needs 192 GiB for codewords. A 44+20 split fixes memory but not the pair
count. A generic native-XOR SAT encoding returned `UNKNOWN` on both sides of
the known dimension-33 case after ten seconds per side.

New artifacts are:

- `explorations/riffle_dp_2lap_g4_component_split_checkpoint.md`;
- `explorations/riffle_dp_2lap_g4_component_split_pairs_s53_v01.json`;
- `explorations/riffle_dp_2lap_g4_component_split_pairs_audit.json`;
- `explorations/riffle_dp_2lap_g4_component_split_sat_s53_v01.json`;
- `scripts/certify_riffle_dp_2lap_g4_component_split_pairs.cpp`;
- `scripts/prepare_riffle_dp_2lap_g4_component_split_words.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_split_pairs.py`;
- `scripts/solve_riffle_dp_2lap_g4_component_split_sat.py`.

Recommended next goal: stop local-distance casework and attack the weaker
global character-bias obligation directly. Determine whether low 24-node
intervals amortize over all 32,772 nodes. Return to local randomization only
if that deterministic global route also fails.

## 2026-08-22 Riffle DenseOuter-ParallelAcc g=4 warmup

The old `outerDense.tex` and `innerAcc.tex` argument has now been packetized.
The distinct candidate is **Riffle DenseOuter-ParallelAcc g=4**. It uses a
terminated random sliding dense outer, one uniform permutation of four-bit
packets, and four lane-parallel accumulators.

A one-lane projection gives a complete asymptotic theorem. For a fixed outer
word of binary weight `w`, some lane has weight at least `w/4`. The marginal
position set of that lane is uniform after the packet permutation. The scalar
accumulator tail therefore gives

```text
Pr[total output weight <= D]
    <= (4 e D / N)^(w/8),
```

where `N` is the number of four-bit packets. Composing this inequality with
the dense-outer generating function proves linear distance for every

```text
delta < (sqrt(2)-1)^8 / (16 e)
      = 1.992416147...e-5.
```

At binary length 2,097,408 and threshold 40, the resulting generating-function
parameter is `z=0.4119542066 < sqrt(2)-1`. The explicit first-moment bound has
20.46 failure bits at outer memory 90 and 40.05 failure bits at outer memory
129.

The old `outerDense.tex` proof has two defects. It emits parity only through
time `k`, so the claimed termination tail is absent. It also claims that every
message span activates `span+M` parity windows. A zero gap longer than `M`
contradicts that claim. Goal 01 repairs both defects: the construction emits
all `k+M` parity symbols, and the generating-function proof decomposes the
message support into `M`-connected clusters.

Goal 02 now improves the one-lane projection without changing the encoder. It
sums the free leading gap exactly and relaxes every attainable accumulator
state path using three facts: transition work is at most twice state work,
zero states cannot be consecutive, and the terminal-state weight has the same
parity as total input weight. A complete sweep of all 2,480 profiles with
input weight at most 80 gives

```text
Pr[output weight <= 40] <= rho^w,
rho = 0.3576340714375845.
```

The worst relaxed profile is attainable: two identical weight-four packets
produce state weights `(4,0)`. Composing the new rho with the repaired dense
outer gives 20.25 failure bits at memory 72 and 40.37 failure bits at memory
108. The one-lane proof required memories 90 and 129, respectively.

The exact construction and proof are in:

- `constructions/riffle_denseouter_parallelacc_g4/CONSTRUCTION.md`;
- `constructions/riffle_denseouter_parallelacc_g4/proof/GOAL_01_PACKETIZED_WARMUP.md`;
- `constructions/riffle_denseouter_parallelacc_g4/proof/GOAL_02_JOINT_LANE_CONTRACTION.md`.

Recommended next goal: couple the dense outer directly to the small repeated
weight-four return profiles. Determine whether the outer packet law suppresses
the `(4,0)` profile enough to reduce memory below 72 without weakening the
uniform bounds for the remaining profiles.

## 2026-08-21 Riffle PacketMul-2Lap g=4 registration

**Riffle PacketMul-2Lap g=4** is the active exploration candidate. The prior
**Riffle DP-2Lap g=4** candidate is paused and remains unchanged.

For every four-bit input packet, setup samples an independent multiplier from
`GF(16)^*`. The encoder multiplies each packet once before the existing packet
permutation and retained-state two-lap map. It samples no fresh multiplier or
permutation between laps.

The construction registry and exact candidate specification are in:

- `constructions/README.md`;
- `constructions/riffle_packetmul_2lap_g4/CONSTRUCTION.md`;
- `constructions/riffle_packetmul_2lap_g4/PROOF_PLAN.md`.

The first proof target is the zero-symbol fraction of the terminal
contribution code. The first refutation target is a mixed-component character
with an unusually large zero-symbol fraction.

## 2026-08-22 Riffle DenseOuter-ParallelAcc g=4 distance frontier

Goal 03 now tracks the finite certified distance separately from the proof
method's ceiling and the unknown true distance. At binary length 2,097,408,
the original uniform moment relaxation composes through distance 124 and
stops at 125 because of two identical `0xf` packets. The exact probability of
that two-packet event is much smaller than its moment bound.

The replay analyzer now counts all one- and two-packet supports exactly and
uses the existing joint-lane moment relaxation for larger supports. A complete
sweep of 62,832 profiles at distance 204 gives

```text
Pr[output weight <= 204] <= rho^w,
rho = 0.4139892039907315 < sqrt(2)-1.
```

Composing this value with the corrected dense outer certifies minimum distance
greater than 204 with 20.44 failure bits at outer memory 97, or 40.45 failure
bits at memory 137. At distance 205, the current relaxation crosses the
composition boundary on four identical `0xf` packets, whose state weights are
`(4,0,4,0)`. This is a proof-method ceiling, not a codeword or a construction
upper bound.

New artifacts are:

- `constructions/riffle_denseouter_parallelacc_g4/GOAL_03_DISTANCE_FRONTIER.md`;
- `constructions/riffle_denseouter_parallelacc_g4/proof/GOAL_03_DISTANCE_FRONTIER.md`;
- `constructions/riffle_denseouter_parallelacc_g4/receipts/goal03_distance204_summary.json`;
- `scripts/analyze_riffle_denseouter_parallelacc_g4_goal02.py`.

Recommended next goal: count four-packet value orders and support gaps exactly,
then resweep with exact bounds for supports of size at most four. Keep the
generic moment bound for larger supports.

## 2026-08-23 Riffle BCHBlockPerm-ParallelAcc g=4

The new proof-gym candidate is **Riffle BCHBlockPerm-ParallelAcc g=4**. It
retains the original GF(2^64) double-parity outer and the extended BCH
`[128,64,22]` local code. Setup samples an independent full 128-bit
permutation for every BCH block. The encoder then forms four-bit packets,
applies one global packet permutation, and runs four lane-parallel
accumulators.

For a fixed BCH word of weight `h`, its local permutation produces a uniform
weight-`h` binary slice. At minimum BCH weight 22, the permuted block occupies
17.0988 packets on average. The maximally packed six-packet outcome has
probability only `8.3565e-17`.

An exact 16-state dynamic program averages the joint-lane gap moment over that
slice. At global length 2,097,408 and threshold 40, it proves

```text
Pr[output weight <= 40 | one BCH block of weight 22]
    <= 2^-95.07548.
```

The old one-lane bound was only `2^-28.14777`; a full global bit permutation
has exact probability `2^-155.21709`. Thus block permutation recovers 66.93
probability bits without dividing the input exponent by four.

At the intended threshold 188,766, the current averaged moment becomes
trivial. A full global bit permutation still has exact one-block probability
`2^-20.14723`. This is a moment obstruction, not refutation evidence. Exact
slice enumeration shows that at least one zero prefix state occurs with
probability 0.68142, the mean number of zero returns is 1.13017, and the
terminal state is zero with probability essentially 1/8. Zero returns are a
bulk feature, not a rare exceptional profile.

Summing the leading gap and all zero-state gaps exactly repairs the target
scale. Conditional on packet support `H` and `m0` zero prefix states, the
valid placement bound is

```text
binom(D, H-m0) * binom(N-H+m0, m0) / binom(N,H).
```

Averaging this expression over the exact weight-22 block slice gives
`2^-11.72084` at threshold 188,766. The dominant terms have 15--18 active
packets and 4--6 zero returns. The bound still treats every positive prefix
state as weight one.

New artifacts are:

- `constructions/riffle_bchblockperm_parallelacc_g4/CONSTRUCTION.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/proof/GOAL_01_BLOCK_SLICE_CONTRACTION.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/GOAL_02_THREE_BLOCK_ZERO_GAP_TRANSFER.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/receipts/goal01_block_slice_summary.json`;
- `scripts/analyze_riffle_bchblockperm_parallelacc_g4_goal01.py`.

Recommended next goal: implement Goal 02 for three independently permuted
weight-22 BCH blocks under their global packet interleaving. Retain the exact
zero-gap count and incorporate positive state weights before summing the field
double-parity spectrum.

## 2026-08-23 Goal 02 three-block transfer

Goal 02 passes for the minimum double-parity profile `(22,22,22)`. The real
three-block law is a uniform 384-bit weight-66 slice conditioned on each
labeled 128-bit part having weight 22. The conditioning event has probability
`2^-6.06145835`, so any superblock event transfers with a 6.061-bit loss.

An exact dynamic program tracks packet support and zero prefix states. It
gives the superblock bound `2^-30.89941485`, hence the real three-block bound
`2^-24.83795650`. A second dynamic program retains the actual Hamming weights
of positive accumulator states. At scaled cost parameter 33.7, it improves
the transferred result to

```text
Pr[output weight <= 188766 | BCH weights (22,22,22)]
    <= 2^-36.44721629.
```

A 20,000-sample search under the real balanced law found a mean of 3.2534
zero prefix states and a maximum of 13. Zero returns remain common, but the
search found no distinct cross-block obstruction. The rigorous calculation
sums every zero-state gap separately.

New artifacts are:

- `constructions/riffle_bchblockperm_parallelacc_g4/proof/GOAL_02_THREE_BLOCK_ZERO_GAP_TRANSFER.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/receipts/goal02_three_block_summary.json`;
- `scripts/analyze_riffle_bchblockperm_parallelacc_g4_goal02.py`.

Recommended next goal: generalize the superblock transfer from `(22,22,22)`
to arbitrary attainable BCH weights `(h1,h2,h3)`. Determine a monotone or
piecewise-uniform envelope before combining the result with the field-outer
weight spectrum.

## 2026-08-23 Goal 03 alpha joint spectrum

The current double-parity schedule is deterministic:
`alpha_i = gamma^i` for `0 <= i < 16384`. A one-data-symbol outer word has
field values `(x,x,alpha_i*x)`. The relevant BCH statistic is therefore a
lagged weight correlation, not independent randomization.

A one-million-pair sample with uniform nonzero `x` and uniform coefficient
index found Pearson correlation `-0.0004265`. At thresholds 48 and 52, the
joint low-weight counts were 1.03 and 0.97 times the sampled independent
predictions. Thus a random lag mixes ordinary BCH weights well.

The initial lags are exceptional. With 100,000 samples per fixed lag, weight
correlation was 1.0 at lag 0, 0.5278 at lag 1, 0.2550 at lag 2, and 0.1225 at
lag 3. Bulk averaging hides this prefix.

The authenticated BCH list through packet support 13 contains 38,560 words,
including 1,420 words of binary weight 22. An exact scan found 3,622 pairs
`(x,i)` for which both `B(x)` and `B(gamma^i*x)` have binary weight 22. Lag 0
alone contributes all 1,420 messages because `alpha_0=1`. All authenticated
weight-22 hits occur by lag 43.

This proves that the current `(22,22,22)` analysis is necessary. Applying the
Goal 02 bound to the 3,622 authenticated inputs gives the partial-family union
bound `2^-24.6246`.

A shifted schedule `alpha_i = gamma^(64+i)` preserves the geometric recurrence
and adds only one fixed multiplication by `gamma^64` per complete outer word.
Its exponent window contains no weight-22 pair from the authenticated
support-at-most-13 tail. This is a distinct, unapproved construction
alternative; the current candidate remains unchanged.

There is also a construction-independent proof route. If the inner bound has
product form `Q(h,h,h') <= r(h)^2 r(h')`, Holder gives

```text
sum_x r(w(x))^2 r(w(alpha*x)) <= sum_x r(w(x))^3.
```

This removes the joint spectrum for every deterministic nonzero alpha. The
required product envelope is not yet proved.

New artifacts are:

- `constructions/riffle_bchblockperm_parallelacc_g4/GOAL_03_ALPHA_JOINT_SPECTRUM.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/proof/GOAL_03_ALPHA_JOINT_SPECTRUM.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/receipts/goal03_alpha_joint_spectrum_summary.json`;
- `scripts/analyze_riffle_bch_alpha_joint_spectrum.py`.

Recommended next decision: either keep the present schedule and target the
Holder-compatible product envelope, or register a distinctly named shifted
alpha candidate and rerun the three-block proof for that ensemble.

## 2026-08-23 shifted minimum shell and Holder gate

The distinct candidate **Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4**
uses `alpha_i = gamma^(64+i)`. The existing Horner recurrence first computes
the unshifted parity value; 64 sparse `xtime` steps then apply the shift once
per complete outer word.

The complete extended-BCH minimum shell was generated from 15 affine orbits.
Each orbit has 16,256 words, all affine images were re-encoded, and the total
equals the committed coefficient `A_22 = 243840`. An exact scan then tested

```text
243840 * 16384 = 3995074560
```

source/exponent pairs over exponents `[64,16448)`. It found zero
weight-22-to-weight-22 intersections. Therefore a one-data-symbol outer word
in the shifted candidate cannot have profile `(22,22,22)`.

The Holder product-envelope route was assessed separately. If
`Q(h,h,h') <= r(h)^2 r(h')`, Holder removes the joint spectrum. However, an
envelope justified by the current certified
`Q(22,22,22) <= 2^-36.44721629` produces a minimum-shell bound of only
`2^-4.55164101`. A 20-bit target needs a 15.448-bit stronger local
certificate; a 40-bit target needs 35.448 bits. A full-global-bit-permutation
comparison would reach `2^-24.4461` after the same shell charge, so 20 bits is
not structurally excluded.

The two approaches should be combined rather than substituted blindly. Use
exact shifted-spectrum exclusions for the lowest weights, where Holder loses
the useful anti-correlation, and use a product envelope for the high-weight
remainder.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal01_shifted_weight22_exact.json`;
- `constructions/riffle_bchblockperm_parallelacc_g4/GOAL_04_HOLDER_PRODUCT_ENVELOPE.md`;
- `constructions/riffle_bchblockperm_parallelacc_g4/receipts/goal04_holder_minimum_shell_gate.json`;
- `scripts/enumerate_ebch128_weight22_affine.py`;
- `scripts/analyze_riffle_shiftalpha64_weight22_exact.py`;
- `scripts/analyze_riffle_bch_holder_gate.py`.

Recommended next goal: compute or bound the adjacent shifted shells beginning
with `(22,22,24)` and `(24,24,22)`, while deriving a monotone Holder-compatible
envelope above a selected BCH-weight cutoff.

## 2026-08-23 ShiftAlpha64 complete weight-22 target spectrum

An optimized exact scanner advances each 128-bit BCH codeword directly under
field multiplication by `gamma`. It evaluated all 3,995,074,560 pairs from
the 243,840 weight-22 messages and exponents `[64,16448)`.

The exact target minimum is 30:

```text
target weight 22: 0
target weight 24: 0
target weight 26: 0
target weight 28: 0
target weight 30: 1401
```

Thus the first profile with repeated source weight 22 is `(22,22,30)`, with
total binary weight 74. The statement does not yet prove that 74 is the
global one-data-symbol minimum, because `(24,24,22)` and `(24,24,24)` remain
unclassified.

The 1,401 weight-30 hits are concentrated at exponents 64--89 and 106--142.
A uniform-alpha model predicts only 3.591 such hits. ShiftAlpha64 therefore
removes weights through 28 exactly, but its low tail is still structured.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_02_WEIGHT22_TARGET_SPECTRUM.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal02_weight22_target_spectrum.json`;
- `scripts/analyze_riffle_shiftalpha64_weight22_targets.cpp`;
- `scripts/analyze_riffle_shiftalpha64_weight22_targets.exe`.

Recommended next goal: classify the weight-24 source shell sufficiently to
decide whether `(24,24,22)` or `(24,24,24)` occurs. This will distinguish a
global one-data minimum of 70, 72, or 74.

## 2026-08-23 ShiftAlpha64 exact one-data minimum

The support-at-most-15 authenticated records generate the complete weight-24
shell: 6,855,968 messages in 532 affine orbits. Every affine image was
re-encoded, and the count equals the committed `A_24`.

The exact scan of 112,328,179,712 weight-24 source/exponent pairs found:

```text
target weight 22: 0
target weight 24: 1175
target weight 26: 0
target weight 28: 1499
target weight 30: 113
```

Together with the weight-22 scan, this closes the one-data minimum:

```text
source h=22: target h' >= 30, total >= 74
source h=24: target h' >= 24, total >= 72, attained 1175 times
source h>=26: BCH distance gives total >= 74
```

Therefore the exact one-data-symbol outer minimum is 72, attained by
`(24,24,24)`. This is not yet a full outer-code minimum or a final distance
proof because several nonzero data symbols may interact.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal02_weight24_affine_enumeration.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal02_weight24_target_spectrum.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/ebch128_weight24_messages.bin`;
- `scripts/enumerate_ebch128_weight24_affine.cpp`;
- `scripts/enumerate_ebch128_weight24_affine.exe`.

Recommended next goal: extend the three-block accumulator transfer from
`(22,22,22)` to the actual shifted minimum profile `(24,24,24)`, then charge
its exact multiplicity 1,175 at distance `0.09N`.

## 2026-08-23 ShiftAlpha64 minimum-profile transfer

The positive-state-weight superblock certificate was generalized to equal
part weight 24. At output threshold `188766`, it proves

```text
Pr[W <= 188766 | (24,24,24)] <= 2^-40.2949878173.
```

The exact shifted spectrum contains 1,175 minimum profiles, so their complete
contribution is `2^-30.0965427758`. This passes a 20-bit target by 10.0965
bits.

For comparison, the unshifted schedule contains 1,364,834 profiles
`(22,22,22)`. Their existing bound contributes only `2^-16.0669222272`.
ShiftAlpha64 therefore gains 14.0296 probability bits on the complete
one-data minimum shell and removes its obstruction at `0.09N`.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_03_MINIMUM_PROFILE_TRANSFER.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal03_minimum_profile_transfer.json`.

Recommended next goal: sum the remaining one-data profiles above total weight
72 using the two exact source-shell histograms first, then a Holder-compatible
tail envelope. After that, move to inputs with two or more nonzero data
symbols.

## 2026-08-23 ShiftAlpha64 end-to-end finite distance floor

Correction: the first-moment distances reported here originally were below a
simple deterministic bound. The outer field code has block distance three,
and every nonzero extended-BCH block has weight at least 22. Thus every
nonzero accumulator input has binary weight at least 66. For accumulator input
`u_t` and state/output `s_t`,

```text
u_t = s_(t-1) + s_t
wt(u) <= 2 wt(s).
```

Consequently, every setup satisfies `d_min >= 33`. The failure probability is
zero. This subsumes both the earlier probabilistic distance-26 statement and
the distance-9 statement at a 20-bit target.

The first open bad-output threshold is 33; excluding it would prove
`d_min >= 34`. The current fixed-word contraction and parity-dropped outer sum
give

```text
sum_(m != 0) Pr[W(m) <= 33] <= 43.401774584603615.
```

This exceeds one, so the probabilistic method does not yet improve the
deterministic distance. The multi-data term dominates because the present sum
discards both parity BCH blocks.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_04_END_TO_END_DISTANCE_FLOOR.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal04_end_to_end_d25.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal04_end_to_end_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal04_end_to_end_d33.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal04_frontier_adjacent.json`;
- `scripts/analyze_riffle_shiftalpha64_end_to_end_floor.py`.

Recommended next goal: retain one parity BCH block in the `r>=2` generating
function. First target a bound below one at bad-output threshold 33, which
would prove `d_min >= 34` with positive setup probability.

## 2026-08-24 ShiftAlpha64 low-block occupation kernel

The previous fixed-word contraction compressed the inner process too early.
For a BCH block-weight profile `h`, the correct inner object is

```text
K_D(h) = Pr[sum_t wt(S_t) <= D | h].
```

The exact first moment is the outer joint profile enumerator weighted by this
kernel. A new wrapper generalizes the pooled-slice 16-state dynamic program to
equal-weight three- and four-block profiles and separately samples the real
independent-block law.

At threshold `D=188766`, the rigorous positive-state bounds are:

```text
(24,24,24):    2^-40.2949878173
(22,22,22,22): 2^-44.6721034989
```

The real-law diagnostics show that the active-state walk is already close to
uniform on `GF(2)^4`. Typical zero returns occur at rate `H/16`, and the
prefix-state weight sum is about `2H`. A four-block weight-64 sample gives the
same pattern.

The rigorous exponential sums have a different saddle. Packet support remains
near its typical value, but zero returns tilt to about `H/3`: 17--20 returns
for the three-block profile and 21--25 for the four-block profile. The current
low-block mechanism is therefore a joint excessive-return and favorable-gap
event, not packet packing.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_05_LOW_BLOCK_OCCUPATION_KERNEL.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal05_lowblock_occupation_summary.json`;
- `scripts/analyze_riffle_shiftalpha64_lowblock_kernel.py`.

Recommended next goal: construct the tilted 16-state return operator. First
reproduce the three- and four-block bounds without tracking the zero-return
count as a separate DP dimension. Then determine whether the operator composes
across a variable number of BCH blocks before summing the outer spectrum.

## 2026-08-24 ShiftAlpha64 compressed return ladder

The Goal 05 kernel now has two rigorous coefficient compressions. The first
uses

```text
C(N-H+m,m) <= theta^(-m) (1-theta)^(-(N-H+1))
```

to remove the zero-return coordinate `m`. The second uses

```text
[x^W] F(x) <= x^(-W) F(x)
```

to remove total binary weight `W`. The remaining dynamic-program coordinates
are packet support `H` and the 16-state accumulator state.

The return-only compression reproduces the Goal 05 reference profiles within
3.5049 bits for `(24,24,24)` and 0.8889 bits for `(22,22,22,22)`. Removing
weight costs another 4.3859 and 4.5850 bits. The double-compressed bounds are
therefore `2^-32.4042` and `2^-39.1982`.

At fixed weight tilt `x=0.334`, the minimum-weight ladder gives:

```text
s=8:   2^-82.1280
s=16:  2^-188.9353
s=24:  2^-317.6716
s=32:  2^-458.2976
s=64:  2^-1073.4895
```

The 64-block state is only 512.3 KiB and runs in a few seconds. There is no
observed equal-minimum-weight inner-kernel obstruction through this range.

After weight compression, the pooled normalization and balance penalty cancel
exactly to `C(128,22)^(-s)`. The displayed balance cost is therefore
bookkeeping, not a separate loss.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_06_COMPRESSED_RETURN_LADDER.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal06_compressed_return_ladder.json`;
- `scripts/analyze_riffle_shiftalpha64_compressed_return.py`.

Recommended next goal: join the fixed-profile kernel to a dyadic outer
block-occupation spectrum. At each occupation step, compare outer multiplicity
growth against the measured kernel exponent. If the exact support coordinate
becomes expensive, compress `H` through the beta-integral representation of
`1/C(N,H)` rather than a worst-case support maximization.

## 2026-08-24 ShiftAlpha64 outer occupation ladder

The first dyadic outer climb now reaches 128 occupied field symbols. The exact
outer field code is `[16386,16384,3]` MDS over `GF(2^64)`, so its symbol-weight
enumerator is explicit.

Equal-weight diagnostics do not close the outer shells. At `s=128`, the MDS
shell has log multiplicity `9139.1433` bits. The compressed inner exponents
are `2389.2939` bits for BCH weight 22 and `4440.9873` bits for BCH weight 64.
These are diagnostics, not complete outer bounds.

A rigorous outer compression now uses both MDS equations. For a fixed support,
write every coordinate as a surjective linear form `L_i(z)` in `s-2` field
variables. Holder gives

```text
sum_z prod_i a(L_i(z)) <= q^(s-3) sum_x a(x)^s,
a(x) = tilt^(-wt(BCH(x))) / C(128,wt(BCH(x))).
```

The exact BCH spectrum evaluates the right side. The resulting shell bounds
are still trivial: `2^183.223` at occupation 3 and `2^302.586` at occupation
4. At occupation 4, the scalar moment is dominated almost equally by BCH
weights 22 and 128. A single coefficient tilt therefore collapses the BCH
spectrum too aggressively.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_07_OUTER_OCCUPATION_LADDER.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal07_outer_occupation_ladder.json`;
- `scripts/analyze_riffle_shiftalpha64_outer_occupation.py`.

Recommended next goal: implement a three-band BCH-weight envelope separating
the low shells, the central spectrum, and the high endpoint. Retain only the
number of occupied symbols in each band. Start with outer occupations 3 and 4,
then resume the dyadic climb.

## 2026-08-24 ShiftAlpha64 three-band outer gate

The outer calculation now separates BCH weights into three bands: light
`22..50`, normal `52..76`, and heavy `78..128`. For each band pattern, two
coordinates are solved by the outer MDS equations. The other coordinates use
the exact BCH-spectrum mass inside their bands.

The split improves the scalar bounds:

```text
occupation 3: 2^183.223 -> 2^160.172  (23.051 bits gained)
occupation 4: 2^272.163 -> 2^231.142  (41.021 bits gained)
```

Both bounds remain worse than the trivial message-count caps `2^103.415` and
`2^179.415`. The dominant patterns are one-light/two-heavy at occupation 3
and two-light/two-heavy at occupation 4. The all-normal patterns are smaller,
but still give `2^18.588` and `2^63.043`.

The present loss is now localized. The calculation maximizes the two solved
values inside their bands. It does not measure whether the actual outer field
equations permit the claimed mixed light-heavy values.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_08_THREE_BAND_OUTER_GATE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal08_three_band_outer_gate.json`;
- `scripts/analyze_riffle_shiftalpha64_threeband.py`.

Recommended next goal: build a compatibility table for the dominant mixed
patterns under the actual outer coefficients. Start with one light and two
heavy values at occupation three, then test two light and two heavy values at
occupation four.

## 2026-08-24 ShiftAlpha64 endpoint triple compatibility

The dominant occupation-three endpoint has now been scanned exactly.
Normalize one outer coordinate to the unique field value whose BCH encoding
is all ones. The first outer equation makes the remaining BCH words
complements. The second equation selects their shared field ratio.

Across all `733141975040` three-position supports, the profile
`(22,106,128)` never occurs. The smallest candidate profile is
`(32,96,128)`. A direct 12-data-block implementation reproduces the optimized
scanner's complete candidate-weight histogram.

The scan contains `2199023247360` normalized finite-support endpoint words.
Transferring their exact histogram through the compressed inner kernel gives
`2^40.7438`. This improves the previous dominant three-band pattern by about
119 bits, but remains trivial.

The remaining endpoint loss is concrete: the inner coefficient envelope
treats the all-one block's packets as arbitrary. They are exactly 32 fixed
`1111` packets.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_09_ENDPOINT_TRIPLE_COMPATIBILITY.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal09_endpoint_triple_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal09_endpoint_triple_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal09_endpoint_triple_transfer.json`;
- `scripts/scan_riffle_shiftalpha64_endpoint_triples.cpp`;
- `scripts/audit_riffle_shiftalpha64_endpoint_triples.py`;
- `scripts/analyze_riffle_shiftalpha64_endpoint_triple_transfer.py`.

Recommended next goal: build an endpoint-aware 16-state inner operator that
retains the 32 fixed `1111` packets and coefficient-bounds only the remaining
64 packets. Recompute the finite endpoint shell before returning to other
occupation-three profiles.

## 2026-08-24 ShiftAlpha64 fixed endpoint closure

The endpoint-aware inner operator now retains the exact all-one BCH block.
After the global packet permutation, it interleaves 32 fixed `1111` packets
with 64 variable packets. Only the variable packets receive the binary-weight
coefficient envelope. The operator passes a direct five-packet enumeration
and the full neutral-mass identity `C(96,32) * 16^64`.

The finite-support endpoint contribution improves from `2^40.7438` to
`2^-44.5527`. Preserving the fixed packets gains 85.297 bits and closes that
family.

Supports containing the second parity position require two normalizations.
An exact scan covers all `134225920` such supports. Both side-block weight
histograms start at 32; neither contains weight 22 or 128. The two transfers
give `2^-70.5135` and `2^-55.4175`, for a combined second-parity bound of
`2^-55.4174`.

The complete occupation-three family containing an all-one BCH block is now
bounded by

```text
2^-44.5519.
```

This is a closed subcase. Occupation-three words without an all-one block
remain open.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_10_FIXED_ENDPOINT_CLOSURE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal10_fixed_packet_transfer.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal10_p1_endpoint_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal10_p1_endpoint_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal10_p1_endpoint_transfer.json`;
- `scripts/analyze_riffle_shiftalpha64_endpoint_fixed_packets.py`;
- `scripts/scan_riffle_shiftalpha64_p1_endpoints.rs`;
- `scripts/audit_riffle_shiftalpha64_p1_endpoints.py`;
- `scripts/analyze_riffle_shiftalpha64_p1_endpoint.py`.

Recommended next goal: subtract the closed all-one family from the
occupation-three analysis and recompute its dominant remaining weight
pattern. Use that result to choose between one more exact small-support table
and a bulk spectral bound.

## 2026-08-24 ShiftAlpha64 minimum triple gate

After removing the all-one family, the old three-band bound is `2^97.464`.
Its dominant realization assigns weights `(22,22,106)`, which violates the
Hamming triangle forced by the first outer equation. An exact weight-profile
envelope enforcing all triangle inequalities gives `2^96.221` and localizes
both support families to `(22,22,22)`.

The complete weight-22 BCH shell contains `1365504` ordered additive triples,
or `227584` unordered triples. Affine symmetry reduces their enumeration to
15 representatives. Their `352044` multiplier ratios were compared with all
finite outer supports. The exact match count is zero.

Supports containing the second parity position use pair correlations rather
than additive triples. The `16384` shifted data coefficients have
`134209536` pairwise differences, all distinct. Therefore their complete
minimum-profile outer count is at most

```text
|L22|^2 = 59457945600 = 2^35.7912.
```

The full-state three-block inner bound is `2^-36.4472`. The second-parity
minimum-profile contribution is consequently at most `2^-0.6561`. Thus both
occupation-three `(22,22,22)` subcases are closed.

After both exclusions, the open occupation-three remainder is `2^94.9201`.
Finite supports dominate. The leading profiles are `(22,106,106)` and
`(22,22,24)`.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_11_MINIMUM_TRIPLE_GATE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_occ3_no_endpoint.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_additive_triples.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_ratio_counts.bin`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_outer_triples.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_data_difference_multiplicity.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_inner_exact.json`;
- `scripts/analyze_riffle_shiftalpha64_occ3_no_endpoint.py`;
- `scripts/analyze_ebch128_weight22_triples.py`;
- `scripts/scan_riffle_shiftalpha64_weight22_outer_triples.rs`;
- `scripts/certify_shiftalpha64_data_difference_multiplicity.rs`.

Recommended next goal: treat `(22,106,106)` and `(22,22,24)` through one
affine-orbit low-shell enumerator. Complementing both weight-106 blocks turns
the first profile into a low-weight additive relation. Compute global relation
counts first; reuse the Sidon budget before building more outer-ratio tables.

## 2026-08-24 ShiftAlpha64 adjacent finite-profile closure

The joint low-shell enumerator gives exact global relation counts:

```text
(22,106,106):  1,365,504 ordered relations
(22,22,24):   27,765,248 ordered relations
```

The first count is exactly the Goal 11 minimum-shell additive count after
complementing both heavy words. The second count comes from testing the XOR
of every minimum-shell affine representative with the complete weight-22
shell and retaining encoded weight 24.

The generic pooled inner transfer is useless for `(22,106,106)` because it
forgets the three source weights. A new multivariate positive-coefficient
operator retains each source's packet count and binary weight. It proves the
one-word bound `2^-85.4188`. Since a fixed relation ratio can occur on at most
`16385*16384` labeled finite supports, the complete profile contributes at
most `2^-37.0377`; no ratio scan is needed.

For `(22,22,24)`, the stronger pooled one-word bound is `2^-29.6521`. Its
27,765,248 relations occupy 7,375,088 ratios. An exact scan of all finite
supports and all three locations of the weight-24 block finds 9,331,089 outer
words. The resulting profile contribution is at most `2^-6.4985`. A direct
12-data-block audit matches the optimized scan exactly.

Together, the two adjacent profiles contribute less than `2^-6.4984`.
Removing them lowers the open occupation-three envelope from `2^94.9201` to
`2^92.2729`. The leading finite profiles are now a cluster:
`(106,106,106)`, `(24,24,24)`, permutations of `(22,22,26)` and
`(22,24,24)`, and `(22,104,106)`.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_12_ADJACENT_PROFILE_CLOSURE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_adjacent_triple_counts.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_adjacent_triples.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_adjacent_inner.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_adjacent_closure.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight222224_outer_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight222224_outer_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_occ3_after_adjacent.json`;
- `scripts/analyze_ebch128_adjacent_triples.py`;
- `scripts/analyze_riffle_shiftalpha64_adjacent_inner.py`;
- `scripts/analyze_riffle_shiftalpha64_adjacent_closure.py`;
- `scripts/scan_riffle_shiftalpha64_weight222224_outer.rs`;
- `scripts/audit_riffle_shiftalpha64_weight222224_outer.py`.

Recommended next goal: build one low-shell/complement occupation-three
envelope for `{22,24,26}` and `{106,104,102}`. Count the additive relations
through one ratio interface and pair them with the source-preserving inner
operator. The objective is to close the whole clustered family at once or
expose a genuinely different obstruction, rather than resume one-profile
casework.

## 2026-08-24 paused BlockAccumulateOnce variant

`BlockAccumulateOnce` is now a named, separate construction. For a binary
outer block code `B_N`, setup samples one uniform global bit permutation and
then applies one unterminated binary accumulator. The outer code may combine
scaled BCH constituents with scaled double parity.

For a fixed outer input weight `h`, the exact number of accumulator inputs
producing output weight `w` is

```text
C(w-1,ceil(h/2)-1) C(N-w,floor(h/2)).
```

The recorded transfer theorem accepts any outer weight-profile envelope
`Abar_h` and bounds the bad-permutation probability by

```text
sum_h Abar_h / C(N,h)
      * sum_{w<=D} C(w-1,ceil(h/2)-1) C(N-w,floor(h/2)).
```

The formula passes exhaustive verification through length 12. The variant is
paused: no scaled BCH–double-parity family is yet proved to satisfy the
profile condition at relative distance `0.09`.

Artifacts are:

- `constructions/block_accumulate_once/README.md`;
- `constructions/block_accumulate_once/CONSTRUCTION.md`;
- `constructions/block_accumulate_once/TRANSFER_THEOREM.md`;
- `constructions/block_accumulate_once/manifest.json`;
- `constructions/block_accumulate_once/transfer_audit.json`;
- `scripts/verify_block_accumulate_once_transfer.py`.

This record does not change the active Riffle construction.

## 2026-08-24 ShiftAlpha64 low-shell/complement cluster

Riffle work resumed with a family-level complement interface. Write each
heavy BCH message as `e+x`, where `e` encodes to the all-one word and `x` is
in a low shell. The first outer check reduces every three-block profile to
one of two equations:

```text
x+y+z = 0  for an even number of complements;
x+y+z = e  for an odd number of complements.
```

The first exact affine-orbit pass gives:

```text
(106,106,106): 0 global relations
(22,22,26):    168,184,576 ordered global relations
```

The `(22,104,106)` and `(24,106,106)` profiles both reuse the Goal 12
`(22,22,24)` relation family. Source-preserving inner bounds close them at
`2^-25.9584` and `2^-34.1042`, respectively. No new ratio scan is required.

After these exclusions, the open occupation-three envelope is `2^90.9275`.
The leading finite profiles are `(24,24,24)`, `(22,22,26)`,
`(106,104,106)`, and `(22,24,24)`. They are instances of the two shared
low-shell equations, not four independent cases.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_13_LOW_SHELL_COMPLEMENT_CLUSTER.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_initial_relations.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_initial_inner.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_initial_closure.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_occ3_after_initial.json`;
- `scripts/analyze_ebch128_goal13_relations.py`;
- `scripts/analyze_riffle_shiftalpha64_goal13_inner.py`;
- `scripts/analyze_riffle_shiftalpha64_goal13_initial_closure.py`.

Recommended next goal: build the shared relation table keyed by three low
shells, one twist bit, and one field ratio. Start with `(22,22,26)` and the
twisted `(22,22,24)` family behind `(106,104,106)`. Test the shifted outer
schedule before undertaking the much larger complete `L24^3` enumeration.

## 2026-08-24 Goal 13 schedule closure

The proposed shared-family step is complete.

For `(22,22,26)`, the affine enumerator emits `168,184,576` ordered
relations. These occupy `45,987,504` distinct nonzero field ratios, with
maximum multiplicity `292`. The generalized shifted-schedule scanner covers
all finite supports and all three possible roles of the weight-26 word. It
finds exactly:

```text
three data coordinates:             30,958,878
first parity plus two data coords:  30,970,511
total:                              61,929,389
```

The pooled source-preserving inner bound is `2^-30.6880351`, so the complete
profile contributes less than `2^-4.8039141`. A direct 12-data-block audit
reproduces all optimized counts, and a separate regression run reproduces
the Goal 12 `(22,22,24)` results after the scanner generalization.

The apparent twisted family admits a stronger closure. If an odd number of
the three BCH words are complements, the three low representatives must XOR
to the weight-128 all-one word. For low weights in `{22,24,26}`, their total
weight is at most `78`. Hamming subadditivity therefore makes every such
relation impossible. This proves zero contribution for the whole
odd-complement low-shell cluster, including `(106,104,106)`, without a ratio
table or schedule scan.

After all Goal 13 exclusions, the rigorous open occupation-three envelope is
`2^90.0358`: `2^90.0279` from finite supports and `2^82.5140` from supports
containing the second parity coordinate. The leading finite families are now
`(24,24,24)`, the orderings of `(22,24,24)`, `(24,104,104)`, and
`(22,102,106)`. Each has even complement parity.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratios_raw.bin`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratio_counts.bin`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight222226_outer_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight222226_outer_scan_12.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight222226_outer_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight222224_outer_scan_12_regression.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_schedule_closure.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_occ3_after_schedule.json`;
- `scripts/scan_riffle_shiftalpha64_equalpair_outer.exe`;
- `scripts/analyze_riffle_shiftalpha64_goal13_schedule_closure.py`.

Recommended next goal: attack the even-complement weight-24 core through one
shared description of `L24^3` and `L22 x L24^2`. Before expanding the entire
weight-24 shell, test whether affine-orbit correlations or a Fourier
convolution can count both relation families and feed the same shifted-
schedule interface.

## 2026-08-24 Goal 14 even-complement weight-24 core

The affine-orbit route succeeded. For an orbit `O` with representative `x`,
the number of partners `y` in a selected shell for which `x+y` lies in a
third shell is constant on `O`. Scanning one representative from each of the
532 weight-24 orbits gives the exact counts

```text
L22 x L24^2:   321,121,024 ordered relations
L24^3:       4,018,336,896 ordered relations
```

The second count is `669,722,816` unordered triples. The enumerator
reproduces the old `1,365,504` count for `L22^3`, obtains the mixed count
from both shell projections, and checks the parity and sixfold symmetry of
the pure weight-24 count.

These relation totals and the source-preserving inner operator close every
selected two-complement realization:

```text
(24,104,104):  2^-27.2989
(22,104,104):  2^-29.5505
(24,104,106):  2^-23.0548
(22,102,106):  2^-22.3250
(26,106,106):  2^-32.8686
```

The low mixed family has a stronger result. Its 321,121,024 ratios compress
to 87,806,530 distinct values, with maximum multiplicity 43. None occurs in
the shifted outer schedule. Thus every finite `(22,24,24)` profile is absent,
including all three placements of the weight-22 block. A direct 12-block
audit reproduces zero, and the generalized scanner still reproduces Goal 12.

The rigorous open occupation-three envelope is now `2^87.8525`, comprising
`2^87.8164` from finite supports and `2^82.5140` from second-parity supports.
The unique leading relation-core profile is `(24,24,24)`.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_14_EVEN_COMPLEMENT_WEIGHT24_CORE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight24_relation_core.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_even_core_inner.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_relation_closure.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight242422_outer_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight242422_outer_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_occ3_after_mixed_schedule.json`;
- `scripts/analyze_ebch128_weight24_relation_core.py`;
- `scripts/write_ebch128_goal14_mixed_ratio_stream.py`;
- `scripts/convert_riffle_ratio_products.rs`.

Recommended next goal: close `(24,24,24)` without producing its 32 GB
ordered-ratio stream. Quotient the 4,018,336,896 ordered relations into
669,722,816 unordered triples. Canonicalize each ratio under the six
transformations induced by permuting a triple, build the canonical schedule
table, and stream affine triple representatives directly into that table.

## 2026-08-24 Goal 15 equal-shell S3 quotient

The proposed quotient succeeded. For a nonzero triple `x+y+z=0`, the
symmetric field value

```text
J(x,y,z) = (x^2+xy+y^2)^3 / (xyz)^2
```

is a complete invariant under common nonzero scaling and permutation of the
three entries. This replaces the six labeled ratios with one projective
unordered key. The only exceptional orbit has `J=0` and projective
automorphism multiplicity three; no weight-24 relation lies in that orbit.

The production schedule has `134,201,345` distinct invariant keys. Its
multiplicity totals reproduce `C(16384,3)` data supports and `C(16384,2)`
supports containing `p0`. An independent implementation checks the algebra
and reproduces a complete 12-block table.

The `669,722,816` unordered `L24^3` relations reduce to 43,113 affine triple
orbits. Explicit stabilizer quotients expand these back to exactly the
committed relation total. The resulting relation set has `188,359,355`
distinct invariant keys. Only 17 keys hit the shifted schedule, producing

```text
three data coordinates:             30,226,386
first parity plus two data coords:  30,235,625
total:                              60,462,011
```

The weight-22 regression reproduces `1,365,504` ordered relations,
`227,584` unordered triples, 14 affine triple orbits, and the prior exact
outer count zero.

Combining the exact outer count with the rigorous Goal 14 inner bound
`2^-32.4042261` closes `(24,24,24)` at `2^-6.5547004`, or about `0.010638`.
The remaining occupation-three envelope is `2^87.28050`. Its leading finite
profiles are `(22,22,28)`, `(24,24,26)`, and `(22,24,26)`.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_15_EQUAL_SHELL_S3_QUOTIENT.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal15_weight24_s3_scan.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal15_weight22_s3_regression.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal15_equal_shell_closure.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal15_occ3_after_equal_shell.json`;
- `scripts/build_riffle_shiftalpha64_s3_schedule.rs`;
- `scripts/scan_riffle_shiftalpha64_equal_shell_s3.rs`;
- `scripts/verify_riffle_shiftalpha64_s3_ratio_quotient.py`;
- `scripts/analyze_riffle_shiftalpha64_goal15_closure.py`.

Recommended next goal: close `(22,22,28)` by projecting onto the two
weight-22 coordinates. Generate the weight-28 coordinate as their XOR; do
not materialize the `1,479,751,168`-word weight-28 shell. Use the existing
equal-pair shifted-schedule interface and keep the Goal 15 invariant quotient
available for later profiles with three equal shells.

## 2026-08-24 Goal 16 packet-weight accumulator formula

The case ladder is paused. The parallel accumulator has a smaller exact
interface than an ordered BCH weight profile. For `k=1,2,3,4`, define `h_k`
as the number of input packets with Hamming weight `k`. Conditional on these
four counts, the randomized packet sequence is uniform and the original BCH
block identities are irrelevant to the inner code.

For packet width `g`, the accumulator state can be quotiented by its Hamming
weight. If the old and new state weights are `a` and `b`, and their supports
overlap in `i` coordinates, then the input-packet weight is

```text
k = a+b-2i
```

and the transition multiplicity is

```text
C(a,i) C(g-a,b-i).
```

This gives a `(g+1)`-state multivariate transfer matrix. At `g=4`, it is a
five-state exact formula retaining `(h1,h2,h3,h4)`. At `g=1`, coefficient
extraction reproduces the classical accumulator formula.

The audit checks:

- all 125 four-bit transition multiplicities;
- the five-state quotient against the weighted 16-state operator;
- every scalar coefficient through length 20;
- all `1,048,576` four-bit input sequences of length five;
- all 58,905 packet histograms of a 128-bit block and every binary weight
  slice from 0 through 128.

All checks pass exactly.

New artifacts are:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_16_PACKET_WEIGHT_ACCUMULATOR_FORMULA.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal16_packet_weight_formula_audit.json`;
- `scripts/verify_riffle_parallelacc_packet_weight_formula.py`.

Recommended next goal: construct a rigorous four-tilt envelope for the exact
conditional kernel `K(h1,h2,h3,h4)`. Combine that envelope with the one-block
packetization polynomial before applying the outer-code sum. Test the result
against the exact Goal 15 occupation-three benchmark. Do not resume the BCH
weight-profile case ladder unless this aggregate interface loses too much.

## 2026-08-24 Goal 17 zero-parity histogram count

The outer count is now being developed as a parity ladder. Stage zero removes
both field-valued parity blocks but retains the BCH block permutations, the
global four-bit packet permutation, and the parallel accumulator. With
`B=16384` data blocks, define the exact one-block expected histogram
polynomial

```text
F(t)=sum_r A_r/binom(128,r)
     * sum_{|h|=32, ||h||=r} P_r(h)t^h.
```

The complete zero-parity count for nonzero messages is exactly

```text
C_0(h) = [t^h](F(t)^16384 - t0^524288).
```

The zero-parity minimum shell has one active weight-22 BCH block. It contains
`16384*243840 = 3,995,074,560` labeled outer words. Only 136 local packet
histograms occur. Their active packet support ranges from 6 to 22. The most
likely histogram is `(14,14,4,0,0)`, with active support 18 and probability
about `0.16768187`. Every coefficient is stored as an exact rational ensemble
expectation.

The weight-22 shell has also been composed exactly with the accumulator. A
period-12 gap quasipolynomial avoids iteration over all 524,288 packet
positions. The exact results are:

```text
D=76000:  expected bad shell words = 2^-0.144159
D=77000:  expected bad shell words = 2^ 0.034384
D=188743: expected bad shell words = 2^12.430624
```

Thus the zero-parity minimum shell crosses expected count one between 76,000
and 77,000 output bits. At the 9% target, the shell alone fails the
first-moment gate by 12.4306 bits. The calculation does not prove a complete
zero-parity distance at 76,000 because higher BCH weights and occupations
remain unsummed.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_17_PARITY_LADDER_ZERO_COUNT.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal17_zero_parity_histogram_count.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal17_zero_parity_accumulator_exact.json`;
- `scripts/analyze_riffle_parity_ladder_zero.py`;
- `scripts/analyze_riffle_parity_ladder_zero_accumulator.py`.

Recommended next goal: compute and compose the one-parity minimum shell in
the same exact format. Its outer values form a repeated-value pair, so its
local conditional polynomial is squared. Keep the complete-code sum as a
separate proof obligation. Add double parity only after the shell comparison.

## 2026-08-24 Goal 18 packet-type coefficient bound

The four-count histogram does not require an explicit `N^4` table. For
positive packet fugacities `x` and `0<z<1`, nonnegative coefficient
domination gives

```text
Pr[wt(y)<=D | h]
 <= z^-D x^-h e0^T M(x,z)^N 1 / T_N(h).
```

One scale in `x` is redundant. Each bound evaluation therefore uses a
five-state matrix power and five effective scalar parameters. The matrix
power retains the random order from the global packet permutation. A
per-packet product envelope was rejected because its zero-packet factor is
one and it permits all inactive packets to occur while the state is zero.

The coefficient bound was tested on the complete zero-parity weight-22
occupation-one shell:

```text
D=76000:  exact log expected=-0.1442, bound=12.6403, loss=12.7844 bits
D=77000:  exact log expected= 0.0344, bound=12.8144, loss=12.7800 bits
D=188743: exact log expected=12.4306, bound=24.8610, loss=12.4304 bits
```

Thus coefficient tilting removes the table but is too lossy for a sparse
shell with almost no margin. It remains plausible for the bulk. Replacing the
bulk sum by its maximum times all formal types would cost 71.4151 bits at
`N=524288`, which is also unacceptable. The current proof shape is hybrid:
exact sparse occupations, coefficient tilting in the bulk, and adaptive
summation or certified boxes between them.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_18_PACKET_TYPE_BOUND.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal18_packet_type_bound_probe.json`;
- `scripts/explore_riffle_parallelacc_type_bound.py`.

Recommended next goal: compute the complete zero-, one-, and two-parity
spectrum for the existing extended BCH `[8,4,4]` construction over `GF(16)`.
Use data-block counts `4,8,12,15`. Compare the exact whole spectrum against
the coefficient bound at every occupation, and use the result to select the
sparse-to-bulk cutoff before returning to the target parameters.

## 2026-08-25 Goal 19 small exact-versus-bound validation

The libOTe enumeration branch was inspected at `origin/enumeration` and
`codex/enumeration-work`. The new proof-gym oracle adopts its separation of
exact rational enumeration, floating-point approximation, exact low-tail
clipping, and reachability pruning. No libOTe source was copied.

The scaled construction uses extended BCH `[8,4,4]`, `GF(16)` double parity
with coefficients `gamma^(4+i)`, independent eight-bit block permutations,
one global four-bit packet permutation, and the four-lane accumulator.

Complete exact rational spectra were computed for `B=4` and `B=8`. An exact
low tail through output weight 16 was computed for `B=12`. The expected-count
crossings are:

```text
B=4:  p=0 n=32 crossing=2; p=1 n=40 crossing=6; p=2 n=48 crossing=9
B=8:  p=0 n=64 crossing=2; p=1 n=72 crossing=6; p=2 n=80 crossing=9
B=12: p=0 n=96 crossing=2; p=1 n=104 crossing=6; p=2 n=112 crossing=10
```

At double parity and output weight eight, the common coefficient tilt loses
`10.8906`, `11.1551`, and `11.2735` bits for `B=4,8,12`, respectively.
Giving each packet histogram its own tilt recovers `4.4014` bits at `B=4`
and `4.4498` bits at `B=8`. The residual typewise losses are `6.4892` and
`6.7053` bits.

This separates two effects. About 4.4 bits are lost by forcing unlike packet
types to share a tilt. Another 6.5--6.7 bits are lost by coefficient/Chernoff
domination itself at this checkpoint. The small family is a bound-validation
gym, not an estimator of the target construction's asymptotic distance.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_19_SMALL_EXACT_VS_TYPE_BOUND.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal19_small_exact_vs_type_bound_b4.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal19_small_exact_vs_type_bound_b8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal19_small_exact_vs_type_bound_b12.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal19_typewise_b4_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal19_typewise_b8_p2_d8.json`;
- `scripts/analyze_riffle_small_exact_vs_type_bound.py`;
- `scripts/probe_riffle_small_typewise_bound.py`.

Recommended next goal: cluster packet histograms into a small number of tilt
regions. Measure how many regions recover most of the stable 4.4-bit
typewise improvement. This is the concrete version of adaptive bulk boxes.

## 2026-08-25 Goal 20 adaptive tilt regions

The clustered-tilt experiment is complete at double parity and `D=8`.
Weighted clustering uses the per-type optimal five-parameter tilts, then
re-optimizes one valid common tilt per region.

```text
regions       1       2       4       8       16      typewise
n=48 log   10.201   9.014   7.478   6.299    5.829      5.800
recovery       0%     27.0%   61.9%   88.6%    99.3%    100%
n=80 log   10.075   9.178   7.325   6.189    5.692      5.625
recovery       0%     20.2%   61.8%   87.3%    98.5%    100%
```

Thus a constant-size region family captures most of the typewise gain at the
two exact sizes. The learned partition is diagnostic because it depends on
the per-type optima.

A new deterministic trim is also active. If `H` is total accumulator input
weight and `W` is output weight, then

```text
H <= 2W - wt(q_N) <= 2W.
```

Hence a bad output with `W<=D` requires `H<=2D`. At `D=8`, only outer binary
weights 12 and 16 can contribute. The exact weight-12 shell supplies 83.1%
of the complete bad tail at `n=48` and 75.6% at `n=80`.

After inserting that shell exactly, all fifteen eligible weight-16 packet
types were separately tilted. The hybrid bound still loses 2.415 bits at
`n=48` and 3.188 bits at `n=80`. Therefore the region count is no longer the
active obstruction. The remaining loss lies inside the conditional bound for
one fixed packet histogram.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_20_ADAPTIVE_TILT_REGIONS.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_clustered_b4_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_clustered_b8_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_sparse_isolation_b4_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_sparse_isolation_b8_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_exact_shell_hybrid_b4_p2_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal20_exact_shell_hybrid_b8_p2_d8.json`;
- `scripts/probe_riffle_small_clustered_tilts.py`;
- `scripts/probe_riffle_small_sparse_isolation.py`;
- `scripts/probe_riffle_small_exact_shell_hybrid.py`.

Recommended next goal: for each eligible weight-16 type, compare the exact
conditional tail, the exact packet coefficient followed by scalar Chernoff,
and the current multivariate coefficient bound. This separates output-tail
Chernoff loss from packet-coefficient loss and determines the next repair.

## 2026-08-25 Goal 21 fixed-type loss decomposition

The weight-16 remainder at `D=8` has now been decomposed. For `H=2D`, the
identity

```text
H = 2W - wt(q_N) - 2 sum_t |supp(q_{t-1}) intersect supp(q_t)|
```

forces every bad path to have `W=D`, terminal state zero, and disjoint
consecutive state supports. Therefore exact packet-coefficient extraction
followed by scalar Chernoff is exact; the scalar-tail loss is zero.

All fifteen eligible packet histograms attain output weight eight. The
measured decomposition is:

```text
                 exact H=16    scalar-Chernoff    boundary bound    loss
n=48                -3.2558        -3.2558            1.4594        4.7152 bits
n=80                -3.1172        -3.1172            1.9709        5.0881 bits
```

Thus the active obstruction is specifically multivariate packet-coefficient
domination. It is not the scalar output tail and not a collection of
exact-impossible packet types.

On the boundary, the exact transfer matrix simplifies to

```text
L(x)[a,b] = binom(4-a,b) x[a+b]  if a+b <= 4,
            0                     otherwise,
```

and the desired count is `[x^h] e0^T L(x)^N e0`. With both packet count and
binary weight fixed, the histogram has three free coordinates. The direct
boundary evaluation still ranges over binary weights and therefore has four
effective tilts. Extracting or conditioning on total binary weight first
would expose the desired three-dimensional local correction.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_21_FIXED_TYPE_LOSS_DECOMPOSITION.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal21_one_type_loss_b4_p2_h16_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal21_one_type_loss_b8_p2_h16_d8.json`;
- `scripts/probe_riffle_small_one_type_loss.py`.

Recommended next goal: compare the exact coefficient in the boundary matrix
with a multivariate saddle-point estimate for every weight-16 type, then turn
the observed three-dimensional local factor into a certified finite-size
upper bound.

## 2026-08-25 Goal 22 boundary saddle factor

The face-aware saddle experiment is complete. For an active packet-weight
face `S`, positive-fugacity evaluation equals the exact coefficient divided
by one tilted point probability. The local Gaussian approximation restores
that probability using the covariance determinant and the lattice index of
reachable histogram differences.

The complete `H=16,D=8` shell gives:

```text
binary length             48       80       112      136
raw coefficient loss    4.715    5.088     5.152    5.160 bits
corrected aggregate     0.250    0.257     0.265    0.270 bits high
```

The measured lattice indices are two and four. A universal parity factor is
not sufficient.

A proportional-scaling test used three base histograms with dimensions one,
three, and four. At scale four, their raw losses were `2.250`, `7.986`, and
`9.852` bits. After the Gaussian lattice correction, their errors were only
`0.011`, `0.060`, and `0.130` bits. From scale two to four, the raw loss grew
by `0.491`, `1.492`, and `2.046` bits, matching the predicted `d/2` bits.

This is strong numerical evidence that the boundary saddle calculation can
give tight full-size estimates. It is not a proof: the Gaussian factor is not
yet a certified upper bound. The exact coefficients use integer arithmetic;
the saddle and covariance calculations use floating point.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_22_BOUNDARY_SADDLE_LOCAL_FACTOR.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal22_boundary_saddle_b4_p2_h16_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal22_boundary_saddle_b8_p2_h16_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal22_boundary_saddle_b12_p2_h16_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal22_boundary_saddle_b15_p2_h16_d8.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal22_boundary_saddle_scaling_s4.json`;
- `scripts/probe_riffle_small_boundary_saddle.py`;
- `scripts/probe_riffle_boundary_saddle_scaling.py`.

Recommended next goal: prove a finite-size local upper bound on the
one-dimensional face `(4s,0,8s,0,0)` using Fourier inversion. Then lift the
same spectral argument to the three- and four-dimensional faces.

## 2026-08-25 Goal 23 exact one-dimensional faces

The proposed one-dimensional Fourier step simplified to an exact
combinatorial formula. Fix packet weight `k` and a histogram supported on
weights zero and `k`. If `h_k` is odd, its boundary coefficient is zero. If
`h_k=2r`, then

```text
A(h,kr) = binom(4,k)^r * binom(h0+r,r).
```

The corresponding exact conditional probability is

```text
Pr[W <= kr | h]
  = binom(h0+r,r) / (binom(h0+2r,2r) * binom(4,k)^r).
```

The proof decomposes every path into forced excursions `0 -> k -> 0`.
The first transition has `binom(4,k)` support choices. Zero packets are
distributed among the `r+1` zero-state gaps.

An exact integer audit passed all 8,580 cases through packet length 64 for
all four values of `k`. At `H=16`, these exact faces contribute 19.0%, 27.9%,
32.4%, and 34.7% of the boundary shell at binary lengths 48, 80, 112, and
136. The weight-two face is the largest individual type at the last three
sizes.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_23_ONE_DIMENSIONAL_BOUNDARY_FACES.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal23_boundary_one_face_formula_audit.json`;
- `scripts/verify_riffle_boundary_one_face_formula.py`.

Recommended next goal: derive the excursion formula on the two-dimensional
even face supported on packet weights zero, two, and four. That face has a
simple three-state boundary graph and should close exactly before the proof
needs general Fourier bounds.

## 2026-08-25 Goal 24 exact two-dimensional faces

Both two-dimensional faces that occur in the weight-16 shell now have closed
exact formulas.

For the even face `(h0,0,h2,0,h4)`, every positive excursion is either
`(4,4)` or `(2,4^j,2)`. Writing `h2=2b`, the coefficient is

```text
6^b * sum_a binom(a+b,a)
            * WeakComp(h4-2a,b)
            * binom(h0+a+b,a+b),
```

where `a=0,...,floor(h4/2)`. The coefficient is zero for odd `h2`.

For the odd face `(h0,h1,0,h3,0)`, every positive excursion is either
`(3,3)` or `(1,3^(2t),1)`. Writing `h1=2b` and `h3=2c`, the coefficient is

```text
sum_a 4^(a+b) * 6^(c-a)
      * binom(a+b,a)
      * WeakComp(c-a,b)
      * binom(h0+a+b,a+b).
```

It is zero unless both `h1` and `h3` are even.

Separate exact integer recurrences checked 20,825 histograms per face through
length 48. All 41,650 comparisons passed.

Goals 23 and 24 together close six of the fifteen boundary types. Exact
coverage rises from 43.0% at binary length 48 to 64.8% at binary length 136.
Using the numerical saddle only on the remaining types leaves aggregate
errors of `0.192`, `0.164`, `0.148`, and `0.140` bits across the four tested
sizes.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_24_TWO_DIMENSIONAL_BOUNDARY_FACES.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal24_boundary_even_face_formula_audit.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal24_boundary_odd_face_formula_audit.json`;
- `scripts/verify_riffle_boundary_even_face_formula.py`;
- `scripts/verify_riffle_boundary_odd_face_formula.py`.

Recommended next goal: certify a local Fourier upper bound on the
three-dimensional face supported on packet weights zero, one, two, and three.
That face contains the largest remaining type at binary length 48.

## 2026-08-25 Goal 25 uniform local-bound candidate

The remaining coefficient problem is now expressed as one Markov-bridge
point probability. A Perron transform converts the tilted boundary matrix
into a finite Markov chain. Conditioning that chain to return to state zero
at time `N` gives exactly the tilted path law used in coefficient extraction.

The candidate bound is

```text
Pr[K=h | q0=qN=0]
 <= 4 * lattice_index / ((2*pi)^(d/2) * sqrt(det(covariance))).
```

It is not yet proved.

After the exact low-dimensional faces are removed, only three support graphs
remain: `{0,1,2,3}`, `{0,1,3,4}`, and `{0,1,2,3,4}`. Their dimensions are
three, three, and four. Exact graph audits show that all are primitive and
all have lattice index two. Thus the proof no longer branches by histogram.

The factor-four candidate passes 44 exact coefficient comparisons. The
largest Gaussian deficit is 1.114 bits, leaving 0.886 bits of slack after the
two-bit safety factor. Combined with the exact faces, its losses over the
complete shell are `1.714`, `1.483`, `1.362`, and `1.300` bits at binary
lengths 48, 80, 112, and 136.

The literature confirms the general Fourier/spectral route for uniform local
limit theorems of finite Markov-additive processes. It does not directly give
the explicit finite factor required for this lattice bridge.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_25_UNIFORM_LOCAL_BOUND_CANDIDATE.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal25_boundary_uniform_local_candidate.json`;
- `scripts/audit_riffle_boundary_uniform_local_candidate.py`.

Recommended next goal: determine the compact saddle-parameter box required
by the target outer sum. Then certify the central Fourier remainder and the
off-center spectral gap uniformly on that box.

## 2026-08-25 Goal 26 target saddle-box attempt

The target box cannot yet be chosen rigorously because no target-length outer
bulk sum identifies which packet histograms retain nonnegligible first-moment
contribution. Numerical Perron calibration also shows that the Goal 25 range
was too narrow. At `N=524352` packets and `D=76000`, a uniform-weight-subset
center has counts approximately
`(388091,121291,14215,740,14.5)` and saddle parameters
`(-1.602,-3.204,-4.806,-6.408)`. Forcing one nonzero packet count to one
produces coordinates as low as `-12.863` over the tested distances. Thus a
compact theorem requires an explicit rare-coordinate split.

The uniform-weight-subset center understates `h4` for the actual block
construction. For a permuted BCH block of weight `r`, the expected number of
full four-bit packets is exactly `32 (r)_4/(128)_4`. Since every nonzero BCH
block has weight at least 22, every fixed outer word of total weight `H`
satisfies

```text
E[h4] >= 0.0009973753281 H.
```

At the `H=152000` boundary this is `151.601`. Exact one-block counts and a
numerically optimized uniform Chernoff envelope give `-47.19` bits. A dyadic
rational witness verifies every BCH-weight inequality and certifies the
per-word statement `Pr[h4 <= 64] <= 2^-47`. This conditional probability
cannot by itself discard the band after summing all outer words.

The target face audit also corrects the proof-gym scope. Across target-length
boundary weights, all 15 nonempty packet-support faces can occur. Six have
exact formulas from Goals 23--24, three were audited as open faces in Goal 25,
and six are newly exposed at the target. Every face graph is primitive, so a
face-parameterized theorem remains plausible.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_26_TARGET_SADDLE_BOX_ATTEMPT.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal26_target_boundary_saddle_calibration.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal26_boundary_h4_gate.json`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal26_target_boundary_faces.json`;
- `scripts/probe_riffle_target_boundary_saddle_box.py`;
- `scripts/probe_riffle_boundary_h4_gate.py`;
- `scripts/audit_riffle_target_boundary_faces.py`.

Recommended next goal: derive a target-length outer boundary envelope indexed
by an `h4` band. Compose the exact BCH packetization polynomial with an outer
occupation bound before optimizing the accumulator saddle. This envelope must
choose the rare-coordinate threshold and the compact Fourier box jointly.

## 2026-08-25 Goal 27 exact proof-gym growth

The proof effort is paused. The current exact enumerator was instead used to
measure the complete legal `GF(16)` ladder with double parity. It clips output
weight above 16, which leaves every reported low-tail coefficient exact.

```text
B=1,...,8:   expected-count-one crossing = 9
B=9,...,15:  expected-count-one crossing = 10
```

The binary output length grows from 24 to 136. The relative crossing falls
from 37.5% to 7.35%. At `B=15`, the cumulative expected-count logs through
weights 8, 9, and 10 are `-1.294`, `-0.361`, and `1.311`, respectively.
Thus the weight-nine tail has fallen below one, but the weight-ten tail
remains about `2^1.31`.

The finite family shows very slow absolute growth and no evidence of linear
distance. This is not a refutation of the target construction. The easy
constituent uses extended BCH `[8,4,4]`, and distinct shifted `GF(16)`
coefficients restrict it to `B<=15`.

Artifacts:

- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/GOAL_27_SMALL_EXACT_GROWTH.md`;
- `constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal27_small_exact_growth_b1_through_b15.json`;
- `scripts/analyze_riffle_small_exact_growth.py`.

Recommended next diagnostic, if desired: decompose the exact `B=15` tail by
outer occupation and binary input weight. This should identify the shell that
keeps the crossing at ten without reopening the proof program.

## 2026-08-26 fixed-rate Outer256 canonical test and multiplicity audit

The active model exploration is now
`Riffle Outer256Model-RandomStepConv g=4 sigma=20`. It preserves the total
packet count and rate. Each modeled rate-one-half constituent maps 128 input
bits to 256 output bits. The baseline random-like, complement-symmetric
`[256,128,38]` spectrum gives a 9% first-moment diagnostic of `+4710.17` bits.
The exact `[128,64,22]` diagnostic was `+8097.18` bits.

The experiment now emits a canonical late-placement test. The baseline
saddle has about 691 active data-pair blocks. An active block has mean binary
weight 101.10 and mean packet support 47.39. The joint modal local profile is weight 100 with
packet support 47. A representative total support is `h=32748`.

For any fixed word with this support, consider the event that all active
packets lie in the final `R=100477` positions. Require the state to remain
nonzero after the first active packet. Also require output weight at most
`D=188766`. This forces a zero prefix of 423875 packet positions. The event
has per-word probability at least `2^-86331.95`, decomposed as follows:

```text
late placement:                         -85259.39 bits
output deviation in the live suffix:    -1072.42 bits
no termination after activation:           -0.14 bits
```

This test needs zero post-activation terminations. However, it does not
witness a distance failure. The relaxed modeled outer coefficient at this
support is `2^79430.47`. Combining it with the event gives expected-count
exponent `-6901.48`.

The saddle profile is model evidence, not an explicit codeword. The suffix
event is a valid probability lower bound for every fixed outer word with
support `h`.

The multiplicity audit confirms that the evaluator does not assign minimum
weight to each active block. It expands

```text
F(x)^8192, where
F(x) = 1 + sum_(w>0,s) A_w Pr[S=s | w] x^s.
```

This exactly sums every data-block weight/support configuration relative to
the guessed spectrum. An unconditioned active block has mean weight 128 and
mean support 60; 97.0% of its mass has weight at least 114. The bad-event
tilt changes the mean weight to 101.10 and mean support to 47.39. It places
78.7% of its mass at weights 90--112 and only 0.077% at weights at most 70.

At `h=32748`, exact FFT coefficient extraction gives a pointwise relaxed
first-moment upper bound of `+4682.90`. The full segment bound is `+4710.17`,
so outer interval aggregation loses only 27.26 bits. The many-block parity
treatment is worst-case; treating the parity pair heuristically as an
ordinary nonzero block recovers about 131 bits, not thousands.

The large remaining uncertainty is the inner probability bound at high
support. At `h=32748`, its upper exponent is `-74747.57`, whereas the
explicit suffix event gives `-86331.95`, a gap of about 11584 bits.

Artifacts:

- `constructions/riffle_outer256model_randomstepconv_g4_sigma20/REPORT.md`;
- `constructions/riffle_outer256model_randomstepconv_g4_sigma20/receipts/fixed_rate_outer256_model.json`;
- `constructions/riffle_outer256model_randomstepconv_g4_sigma20/receipts/outer_multiplicity_audit.json`;
- `scripts/analyze_riffle_outer256model_randomstepconv_g4_sigma20.py`;
- `scripts/audit_riffle_outer256model_multiplicity.py`.

Recommended next goal: test the fixed-rate 512-bit constituent model. Every
reported distance result must include the canonical profile: active-block
count, local weight, packet support, placement window, output deviation, and
the number of required inner terminations.

## 2026-08-26 exact RM512 fixed-rate diagnostic

The 512-bit test uses the exact RM(4,9) `[512,256,32]` spectrum. Four GF(2^64)
data symbols form each 256-bit RM input, giving 4096 data groups. The two
field-parity symbols retain one 256-bit rate-half output block. The total
remains 524352 four-bit packets.

The committed RM spectrum has integer mass `2^256`, complement symmetry, and
minimum distance 32. It matches the Type-II Gleason formula exactly. Packet
support conditional on each binary weight is counted exactly.

The 9% floating first-moment upper bound is `+3003.952367`. The earlier
Outer256 model gave `+4710.167631`, so RM512 recovers 1706.22 bits without
closing the bound. The dominant support interval remains `16384--32767`,
with occupation mode 343.

The full-spectrum profile is:

```text
                          raw active block    bad-event conditioned
mean binary weight              256.00                203.51
mean packet support             120.00                 95.40
central 98% weight range       228--284               172--236
```

Thus the evaluator uses the full RM multiplicities rather than minimum
weight 32. The conditioned number of active groups has mean 343.48, standard
deviation 17.74, and 1%--99% interval 303--385.

At relaxed data support `h=32767`, the exact FFT coefficient audit gives:

```text
relaxed RM data coefficient:       +77769.24 bits
inner probability upper bound:     -74777.34 bits
pointwise first-moment upper:        +2991.90 bits
whole support-segment upper:         +3003.95 bits
```

The segment machinery loses only 12.05 bits. The outer sum is not the source
of the positive result.

The concrete lower event charges all 64 possible parity packets, for full
support 32831. All packets land in the final 100477 positions, the state
never terminates after activation, and output weight is at most 188766. Its
per-word exponent is `-86569.36`. Combining it with the exact relaxed RM
coefficient gives expected-count exponent `-8800.12`. It is not a refutation
witness.

Artifacts:

- `constructions/riffle_rm512_randomstepconv_g4_sigma20/README.md`;
- `constructions/riffle_rm512_randomstepconv_g4_sigma20/REPORT.md`;
- `constructions/riffle_rm512_randomstepconv_g4_sigma20/manifest.json`;
- `constructions/riffle_rm512_randomstepconv_g4_sigma20/receipts/rm512_randomstepconv_distance09.json`;
- `scripts/analyze_riffle_rm512_randomstepconv_g4_sigma20.py`.

Recommended next goal: audit the one-lap inner probability bound near packet
support 32768. The current inner upper exponent is about 11792 bits above the
concrete parity-uniform suffix event. Determine whether a richer bad event
narrows that gap before interpreting another positive upper bound as a design
failure.

## 2026-08-26 one-lap episode-structure audit

The apparent RM512 bulk obstruction was caused by the compressed inner
envelope. At support 32767, the end-to-end calculation reused parameters
optimized at support 16384. Pointwise optimization changes the inner upper
exponent from `-74777.34` to `-85557.15`, a recovery of 10779.31 bits. The
exact relaxed RM outer coefficient is `+77769.24`, so the pointwise bulk
contribution is `-7787.91`.

Termination separation itself recovers only about 0.50 bits. The main loss
was the wide support interval, not the common-radius formula.

After pointwise optimization, the largest sampled contribution moves to data
support 95:

```text
exact relaxed RM coefficient:       +228.9165 bits
fixed-support inner upper:          -197.2281 bits
pointwise sum:                       +31.6883 bits
dominant termination count:                 5
mean termination count:                  5.66
```

One active RM group supplies the full support-95 coefficient to displayed
precision. Two active groups are 40.89 bits smaller. The obstruction is
therefore a one-group parity-kernel problem.

A concrete lower family was derived. For full support H, termination count
e, and total live length L, count e terminated episodes followed by one live
episode that survives to the final position. The trajectory count is

```text
binom(L-1,e) * binom(N-L+e,e) * binom(L-e-1,H-e-1).
```

An exhaustive state-path enumerator checks this identity in all 210 cases
through seven packet positions.

For data support 95 and zero parity support, summing e=0,...,15 gives an
explicit expected-count exponent of `+8.2824`. The dominant rows have five
or six terminations and total live length about 94161. This refutes the
parity-relaxed model.

The same explicit family crosses below one after only four parity packets:

```text
parity packet support:       0      1      2      3       4       64
expected-count exponent:  8.28   5.90   3.51   1.13   -1.26  -145.02
```

Thus late activation remains the live-length mechanism, but it is split
across about six episodes. The proof/refutation boundary is now the joint
spectrum of one RM data group and its two field parities. In particular, the
next analysis must bound data packet support near 95 conditioned on parity
packet support at most three. The parity-zero kernel is the first case.

Artifacts:

- `constructions/riffle_rm512_randomstepconv_g4_sigma20/REPORT.md`;
- `constructions/riffle_rm512_randomstepconv_g4_sigma20/receipts/inner_episode_structure_audit.json`;
- `scripts/audit_riffle_randomstepconv_episode_structure.py`.

Recommended next goal: specify the RM message basis and the exact parity
encoder, then enumerate or bound the packet-support spectrum of the
128-dimensional parity-zero RM subcode. Without this joint object, the
RM512 model is not sufficiently specified to decide the explicit low-support
failure family.

## 2026-08-26 full-distance random-constituent pivot

The RM512 exploration is paused. Its exact spectrum made it useful as a
diagnostic, but its minimum distance 32 is deliberately weaker than the
constituent desired for the structural proof gym.

The active model is now **Riffle
FrozenRandom512-P2-RandomStepConv g=4 sigma=20**. Each of the 4096 data-group
positions receives an independent uniform linear injection from 256 to 512
bits. Setup freezes these maps. The groupwise independence is a proof device
that makes the ensemble first moment factor exactly. It is not an
implementation recommendation.

The model uses the exact unconditioned random-linear ensemble spectrum. It
does not impose the earlier distance-56 floor. A random constituent typically
has minimum distance near 58. The baseline retains two field parities encoded
by extended BCH `[128,64,22]`.

The parity-zero part of the old one-group obstruction does not require a
joint RM spectrum. Two independent field parity equations leave at most 128
input dimensions. Including the 4096 group positions gives at most `2^140`
candidates, versus the existing support-95 inner upper exponent `-197.228`.
This class therefore has at most `-57.228` bits before further parity gains.

Artifacts:

- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/README.md`;
- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/GOAL_01_RANDOM_SPECTRUM_BASELINE.md`;
- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/manifest.json`;
- `scripts/analyze_riffle_frozenrandom512_p2_randomstepconv_g4_sigma20.py`.

The initial random-spectrum diagnostic gives `+3003.95` in the old wide
interval and `+31.6883` at the pointwise support-95 row when parity is
dropped. The random support-95 coefficient is `+228.9165`, essentially the
same as RM. Thus RM's low minimum distance was not the cause of that moderate-
support row.

The exact one-active-group calculation retains both field parities. Their
rank-two map is surjective on one four-symbol group, so the parity pair ranges
uniformly. Composing this distribution with the exact extended-BCH
`[128,64,22]` spectrum gives total one-group exponent `-91.830594`. The
dominant profile has data support 95, zero parity support, and exponent
`-96.311670`.

The random ensemble is now unconditioned. Its expected binary spectrum is
`((2^256-1)/(2^512-1))*binom(512,w)`, and its expected four-bit packet
spectrum is `((2^256-1)/(2^512-1))*binom(128,s)*15^s`. A conjectured BCH
spectrum may later replace this baseline without changing the inner analysis.

Artifacts added:

- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/REPORT.md`;
- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/receipts/frozenrandom512_p2_distance09.json`;
- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/receipts/onegroup_exact_p2.json`;
- `scripts/audit_riffle_frozenrandom512_p2_onegroup.py`.

Recommended next goal: analyze two and then several active data groups while
retaining the global parity sums. Compare two through six field parities only
after identifying the dominant multi-group cancellation profile.

## 2026-08-26 random-outer global pointwise retry

The global parity sum can be retained efficiently under independent random
data encoders. For one active group, the exact parity-pair distribution is
used. For at least two active groups, the pair is uniform up to a relative
correction below `2^-512`. The exact extended-BCH packet spectrum supplies the
parity moment.

The pointwise diagnostic evaluates every support through 512, every eighth
support through 4096, and checkpoints through 32767. Its sampled maximum is

```text
total packet support:                 360
outer coefficient Chernoff upper:  +751.933161
pointwise inner upper:              -836.334675
combined exponent:                   -84.401514
```

The canonical tilted profile has occupation mode three and mean occupation
3.748. Its 10%, 50%, and 90% occupation quantiles are two, four, and six.
Each active data group contributes 94.20 packets on average. Both parity
blocks together contribute 6.94 packets. The dominant inner term has eleven
terminations.

If a complete cover keeps every support below `-84.4015`, summing all 524352
supports gives `-65.4013`. This leaves 25.40 bits beyond the `2^-40` target.
The current scan is evidence, not a certificate, because it does not yet
bound every unscanned support.

Artifacts:

- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/GOAL_02_GLOBAL_POINTWISE.md`;
- `constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/receipts/global_pointwise_p2.json`;
- `scripts/audit_riffle_random512_p2_global_pointwise.py`.

Recommended next goal: turn the sampled curve into an adaptive interval cover
of every support. Keep the exact global parity moment and pointwise inner
parameters. Do not add more parity unless that cover reveals a new positive
region.

## 2026-08-26 setup-entropy minimization

The active entropy candidate is **Riffle
FieldPair512-P2-ScalarStepConv g=4 sigma=18**. It preserves the exact marginal
laws used by the random-code first moment while replacing dense matrices.

For data group `i`, setup samples nonzero `(a_i,b_i)` in `GF(2^256)^2` and
encodes `m_i` as `(a_i*m_i,b_i*m_i)`. This uses 512 bits per block instead of
131072. For inner step `j`, setup samples `c_j` in `GF(2^(g+sigma))` and maps
the input-state vector by field multiplication. At sigma 18 this uses 22 bits
per step instead of a 22-by-22 random matrix. Both replacements give a uniform
image for every fixed nonzero input and are entropy-optimal for that exact
law.

The setup-entropy ledger is:

```text
4096 outer field-pair maps:             2,097,152 bits
524352 inner scalar maps:              11,535,744 bits
uniform global packet permutation:      9,206,311.134 bits
two parity-block permutations:              1,432.323 bits
total:                                 22,840,639.457 bits
```

The dense-matrix baseline used 848,105,407.457 bits. The exact-law reduction
is 37.13-fold.

Sigma 18 is the current boundary. Exact-saddle interval sums are `-53.130967`
on supports 50--399, `-43.600759` on 400--1400, and `-56.600958` on
1401--3200. Their combined exponent is `-43.598633`. The unsummed tails are
small in the Chernoff scan but are not yet certified. Sigma 17 fails the
current method: its central interval alone has upper exponent `+2.273900`.

Artifacts:

- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/CONSTRUCTION.md`;
- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/GOAL_01_ENTROPY_LEDGER.md`;
- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/REPORT.md`;
- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/manifest.json`;
- `scripts/audit_riffle_random512_p2_exact_saddle.py`.

Recommended next goal: certify the sigma-18 support tails and rounding. After
that, attack the 9.21-million-bit uniform permutation using a structured
interleaver or a low-entropy permutation family. The inner scalar stream is
already entropy-minimal under the current transfer law.

## 2026-08-26 smaller-inner parity trade

The active entropy point is now **Riffle
FieldPair512-P4-ScalarStepConv g=4 sigma=16**. Four field-parity symbols allow
the retained inner state to fall from 18 to 16 bits. Its setup entropy is
`21,795,863.798` bits, approximately 2.598 MiB. This saves
`1,044,775.660` bits relative to P2 sigma 18.

The exact-saddle sum on supports 900--3000 is `-44.905837`; the largest exact
point is `-53.380379` at support 1813. This remains diagnostic until the
support tails and floating-point rounding are certified.

Sigma 15 cannot be obtained by simply appending more global parity. P5 sigma
15 has sampled maximum `+87.981540` at support 288. P6 sigma 15 has sampled
maximum `+79.371762` at support 320. The canonical obstruction is one active
512-bit data block whose data and parity outputs occupy nearly all available
packets. That family contains about `2^256` messages, and a 15-bit inner state
does not suppress its low-output trajectories enough. More parity lengthens
this one-block family without creating independently active data blocks.

Artifacts:

- `constructions/riffle_fieldpair512_p4_scalarstepconv_g4_sigma16/README.md`;
- `constructions/riffle_fieldpair512_p4_scalarstepconv_g4_sigma16/REPORT.md`;
- `constructions/riffle_fieldpair512_p4_scalarstepconv_g4_sigma16/manifest.json`;
- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/receipts/global_pointwise_p5_sigma15.json`;
- `constructions/riffle_fieldpair512_p2_scalarstepconv_g4_sigma18/receipts/global_pointwise_p6_sigma15.json`.

Recommended next goal: test a structural route to sigma 15 rather than more
parity. The two clean candidates are a wider packet or a different outer
constituent geometry. State the allowed rate and constituent-size constraints
before claiming a global entropy minimum.

## 2026-08-26 outer-size versus inner-state curve

The generalized random-outer enumerator now accepts rate-one-half outer
lengths from 128 through 2048 bits. The clean P2 control curve is:

```text
outer 256:  sigma 32 central sum -38.848659; no target-reaching sigma seen
outer 512:  sigma 18 evaluated sum -43.598633; sigma 17 fails
outer 1024: sigma 17 central sum -44.042347; sigma 16 fails at +13.984046
outer 2048: sigma 17 central sum -57.715141; sigma 16 fails at -10.273467
```

The 256-bit obstruction approaches a zero-termination inner profile, so more
state has almost no remaining effect. The clean curve therefore shows one
state bit gained from 512 to 1024 bits and no further integer gain from 1024
to 2048 bits.

For the active P4 family, outer 512 with sigma 16 has exact central sum
`-44.905837`. Outer 1024 with sigma 15 fails at sampled maximum `-7.394418`,
whereas sigma 16 has sampled maximum `-112.912186`. Thus the larger outer
buys margin but not another state bit. Both sizes use the same total outer
setup entropy because block length times block count is constant.

P4 at outer 256 remains unresolved. Its two 64-bit input symbols generate four
correlated parity symbols through an MDS `[4,2,3]` map. The present minimum-
support bound is too loose. A reliable point needs their joint extended-BCH
packet-support moment.

Artifact:

- `constructions/riffle_fieldpair512_p4_scalarstepconv_g4_sigma16/OUTER_INNER_CURVE.md`.

Recommended next goal: derive or tightly upper-bound the P4/256 joint parity
moment. This is the decisive calculation for whether the 512-bit outer can be
made smaller without abandoning the 16-bit inner neighborhood.

## 2026-08-26 random-outer landscape engine and packet-width contours

The parity-free landscape model fixes a `2^20`-bit message, independent random
rate-half `[B,B/2]` outer blocks, a global permutation of `g`-bit packets, one
lap of RandomStepConv with state size sigma, and a 9% distance target. It
reports

```text
lambda(B,g,sigma) = -log2 E[number of nonzero bad codewords].
```

The target contour is lambda 40. Negative values are retained because they
show the magnitude and mechanism of failure.

The evaluator now accepts `g` in `{1,2,4,8}` and arbitrary compatible `B`.
It separates sparse outer occupations before using a global saddle. This
avoids the large false coefficient produced when a tilted distribution is
mostly the all-zero message. The inner transfer table is cached by
`(g,sigma,support)` and reused across outer sizes.

The exact-moment matrix refinement gives these `g=4` results:

```text
B=256:  zero-termination floor lambda ~= -80.78
B=512:  zero-termination floor lambda ~= -10.51
B=1024: zero-termination floor lambda ~=   3.86
B=2048: zero-termination floor lambda ~=  32.62
B=4096: sigma 19 gives 27.18; sigma 20 gives 47.25; floor ~= 91.53
```

The first observed 40-bit crossing is therefore `B=4096`, between sigma 19
and 20.

The first cross-packet-width contours are:

```text
g=1: B=256,  sigma 23 gives 38.26, sigma 24 gives 40.16
g=2: B=512,  sigma 18 gives 29.97, sigma 19 gives 43.12
g=4: B=4096, sigma 19 gives 27.18, sigma 20 gives 47.25
g=8: B=4096, sigma 40 gives -43958.12
```

For `g=8`, the dominant profile is a bulk late-start family with about 47,600
active packets and 99 active outer blocks. Increasing `B` to 8192 barely moves
the floor, so no nearby 40-bit contour was observed.

Artifacts:

- `scripts/analyze_riffle_random_outer_landscape.py`;
- `constructions/riffle_random_outer_landscape/README.md`;
- `constructions/riffle_random_outer_landscape/REPORT_G4_INITIAL.md`;
- `constructions/riffle_random_outer_landscape/REPORT_PACKET_WIDTH.md`;
- `constructions/riffle_random_outer_landscape/receipts/g1_floors_b128_2048.json`;
- `constructions/riffle_random_outer_landscape/receipts/g1_b256_sigma8_32.json`;
- `constructions/riffle_random_outer_landscape/receipts/g1_b256_sigma23.json`;
- `constructions/riffle_random_outer_landscape/receipts/g2_floors_b256_4096.json`;
- `constructions/riffle_random_outer_landscape/receipts/g2_b512_sigma18_19.json`;
- `constructions/riffle_random_outer_landscape/receipts/g4_matrix_floors_b256_4096.json`;
- `constructions/riffle_random_outer_landscape/receipts/g4_b4096_matrix_sigma18_20.json`;
- `constructions/riffle_random_outer_landscape/receipts/g8_b4096_sigma40_fullsupport.json`.

Recommended next goal: fill intermediate sigma values around the three viable
contours and certify the six cells immediately below and above lambda 40.

## 2026-08-26 structured full-stripe pivot

The new family is **Riffle S-Stripe RandomStepConv g=8**. Arrange outer
packets as a block-by-coordinate matrix. Divide packet coordinates into S
equal classes, independently permute the packets in each class, and emit the
classes as consecutive inner regions. The full-stripe instance uses one
packet coordinate per region.

The exact full-stripe evaluator conditions on the active outer-block count.
For every output tilt, a log-domain dynamic program computes all four entries
of every regional candidate-count coefficient exactly. It then composes the
two-state regional matrix across all stripes and sums every occupation.

At B=4096 and sigma 40, full striping improves lambda from about -43959 to
-36132. The remaining tilted path alternates late activity in one stripe with
early activity and termination in the next. It uses about 255 resets across
512 stripes. Thus striping replaces the bulk late suffix with an
alternating-boundary obstruction.

The sigma-40 outer curve is:

```text
B=4096:  lambda -36132.37
B=8192:  lambda -27570.36
B=16384: lambda -10260.14
B=32768: lambda +10278.20
```

The first refined contours are:

```text
B=8192:  sigma 95 still lambda -903.15; one-active-block floor in this bound
B=16384: sigma 50 gives -18.71; sigma 51 gives +1005.29
B=32768: sigma 27 gives -1325.98; sigma 28 gives +725.78
```

Artifacts:

- `constructions/riffle_striped_randomstepconv_g8/CONSTRUCTION.md`;
- `constructions/riffle_striped_randomstepconv_g8/REPORT.md`;
- `scripts/analyze_riffle_striped_random_outer.py`;
- receipts under `constructions/riffle_striped_randomstepconv_g8/receipts/`.

Recommended next goal: certify either B=32768, sigma 28 or B=16384, sigma 51
with outward-rounded interval arithmetic. The B=32768 cell is the cleaner
low-memory target; B=16384 is the smaller-outer comparison point.

## 2026-08-26 full-stripe g=4 comparison

The distinct packet-width-four family is **Riffle S-Stripe RandomStepConv
g=4**. The exact full-stripe evaluator uses a wider output-tilt grid for the
sparse one-block regime.

The useful trade points are:

```text
B=512,  sigma=20: lambda +118.52
B=1024, sigma=14: lambda  +35.37
B=1024, sigma=15: lambda +110.03
B=2048, sigma=12: lambda -470.79
B=2048, sigma=13: lambda  +84.32
```

The old global-permutation g=4 contour was B=4096, sigma 20, lambda 47.25.
Full striping therefore improves both the usable outer length and the inner
memory. The balanced current point is B=1024, sigma 15. The B=512 memory
boundary remains open.

Artifacts:

- `constructions/riffle_striped_randomstepconv_g4/CONSTRUCTION.md`;
- `constructions/riffle_striped_randomstepconv_g4/REPORT.md`;
- `constructions/riffle_striped_randomstepconv_g4/receipts/summary.json`.

Recommended next goal: locate the B=512 memory boundary, then compare the
operation cost of B=512/sigma-boundary, B=1024/sigma15, and B=2048/sigma13.

## 2026-08-26 transpose-first variants

Record these as two different construction families.

1. **Riffle TransposePacketShuffle-RandomStepConv g** partitions the outer
   blocks into fixed groups of size \(g\), transposes the outer-block matrix,
   forms one packet from each fixed group in every row, and independently
   shuffles the \(L/g\) packets in each row. This is the active version.
2. **Riffle TransposeBitShuffle-RandomStepConv g** transposes first,
   independently shuffles all \(L\) individual bits in each row, and then
   forms consecutive \(g\)-bit inner inputs. This stronger but more expensive
   version is recorded and paused.

For one active outer block, the two versions induce the same inner input
support distribution. They differ at higher occupation because the packet
version retains fixed collisions while the bit version redraws collisions in
every row.

The packet evaluator conditions on \(a\) active outer blocks and packs them
into \(\lfloor a/g\rfloor\) full fixed groups plus at most one partial group.
For each output tilt, a log-domain matrix-polynomial dynamic program computes
the exact random packet-order transition for one row. It then raises that
transition through all \(B\) independently permuted rows and sums all active
outer-block counts.

Initial result at \(n=2^{20}\), relative distance 0.09:

```text
g=4, B=1024, sigma=15: provisional lambda +195.283900
```

The earlier S-Stripe point at the same parameters was lambda +110.033707.
The new dominant occupation is one active outer block: outer multiplicity
523.00 bits, inner tail -718.28 bits. That case is exact and does not use the
packing assumption. The two-active-block pointwise exponent is -437.35, about
242 bits below the dominant exponent -195.28.

The packing-domination lemma is now closed. One packing move changes
occupancies \((r,s)\) to \((r+1,s-1)\), where \(g>r\ge s>0\). It preserves
the probability of zero nonzero packets and lowers the probability of two.
Uniform row shuffling converts this count domination into an input-support
inclusion coupling. The averaged RandomStepConv output weight is monotone
under that inclusion.

An independent 80-decimal-digit interval calculation at fixed Chernoff
parameters first certified the sparse terms. The full checker now evaluates
all 2048 occupations with upward-rounded mantissa-plus-exponent arithmetic.
It proves

```text
distance threshold:       188743 of 2097152 output bits
lambda lower bound:       195.283900426038041150079401119
target lambda:            40
dominant occupation:      1 active outer block
```

Therefore the random setup has minimum distance at least 188744, which
exceeds 9% of the output length, except with probability at most
`2^-195.283900426038`.

Artifacts:

- `constructions/riffle_transpose_packetshuffle_randomstepconv/`;
- `constructions/riffle_transpose_bitshuffle_randomstepconv/`;
- `scripts/analyze_riffle_transpose_packetshuffle.py`;
- `scripts/certify_riffle_transpose_packetshuffle_sparse.py`;
- `scripts/certify_riffle_transpose_packetshuffle_full.py`;
- `constructions/riffle_transpose_packetshuffle_randomstepconv/proof/PACKING_DOMINATION.md`.
- `constructions/riffle_transpose_packetshuffle_randomstepconv/proof/DISTANCE_CERTIFICATE.md`.

Recommended next goal: audit the certified ensemble against the intended
encoder cost model, then search for the largest certified distance at the same
\(g=4,B=1024,\sigma=15\) parameters.

## 2026-08-26 certified g=4 outer-memory tradeoff

Keep \(n=2^{20}\), \(g=4\), and distance threshold 188743 fixed. A sparse
one-block scan located the power-of-two boundary cells. Complete
outward-rounded calculations then certified all active-block occupations.

```text
B=1024, sigma=13: lambda 83.087891407589
B=512,  sigma=15: lambda 57.641932041422
B=256,  sigma=18: lambda 45.400358851794
```

The adjacent lower-memory cells have one-block bounds below 40:

```text
B=1024, sigma=12: lambda -2.29
B=512,  sigma=14: lambda 34.97
B=256,  sigma=17: lambda 37.28
```

For \(B=128\), the current one-block bound approaches only about 21.38 bits
even as \(\sigma\) grows. This is a limitation of the present proof bound, not
a refutation of the construction.

Using dense binary bit products as a cost proxy gives:

```text
B=1024, sigma=13: outer 1,073,741,824; inner 151,519,232; total 1,225,261,056
B=512,  sigma=15: outer   536,870,912; inner 189,267,968; total   726,138,880
B=256,  sigma=18: outer   268,435,456; inner 253,755,392; total   522,190,848
```

Thus \(B=256,\sigma=18\) is cheapest under this proxy. The exact power-of-two
ladder has no certified point with both \(B<1024\) and \(\sigma<15\).
Intermediate padded outer lengths between 512 and 1024 are the natural place
to search for such a point.

Artifacts:

- `constructions/riffle_transpose_packetshuffle_randomstepconv/TRADEOFF_G4.md`;
- `scripts/analyze_riffle_transpose_packetshuffle_sparse_curve.py`;
- receipts `g4_b1024_sigma13_full_interval.json`,
  `g4_b512_sigma15_full_interval.json`, and
  `g4_b256_sigma18_full_interval.json`.

Recommended next goal: benchmark the three certified cells with identical
low-level kernels, or define the padded intermediate-\(B\) construction and
search for a point that reduces both outer length and inner state.

## 2026-08-26 packet g=8 and bit-transpose comparison

The packet-transpose construction does not scale cleanly from \(g=4\) to
\(g=8\). Fixed groups let a bulk set of active outer blocks remain packed in
the same packets in every row. At \(B=1024,\sigma=12\), the dominant event is
near 424 active blocks and lambda is about -49417. Complete certificates find
these first positive power-of-two cells:

```text
B=1024,  sigma=102: lambda  140.289349191542; dominant occupation 424
B=2048,  sigma=54:  lambda  601.228272652283; dominant occupation 1
B=4096,  sigma=30:  lambda 1011.906471547732; dominant occupation 1
B=8192,  sigma=17:  lambda 2238.062870150003; dominant occupation 1
B=16384, sigma=11:  lambda 3626.465383813988; dominant occupation 1
B=32768, sigma=8:   lambda 5191.959819043644; dominant occupation 1
```

The distinct **Riffle TransposeBitShuffle-RandomStepConv g** construction
shuffles all \(L\) individual bits independently in every transposed row and
only then forms \(g\)-bit packets. For \(a\) active blocks its exact row
transition is the degree-\(a\) coefficient of

```text
[sum_r C(g,r) x^r M_r]^(L/g) / C(L,a).
```

The new outward-rounded checker evaluates these coefficients using scaled
binary64 arithmetic. It computes high occupations through the reversed
polynomial. Exhaustive small-instance tests cover both directions.

Complete bit-transpose g=8 certificates give:

```text
B=1024, sigma=12: lambda 83.201070356493
B=512,  sigma=14: lambda 57.659313166668
B=256,  sigma=17: lambda 45.402590000096
```

Every cell is dominated by one active outer block. The corresponding g=4
boundary is \((13,15,18)\), with margins \((83.0879,57.6419,45.4004)\). Thus
g=8 saves one inner-state bit at every tested outer length. Under the dense
bit-product proxy, the cheapest point is bit-transpose
\(g=8,B=256,\sigma=17\), at 432,275,456 total products. This excludes the
cost of shuffling individual bits.

Artifacts:

- `constructions/riffle_transpose_packetshuffle_randomstepconv/TRADEOFF_G8.md`;
- `constructions/riffle_transpose_bitshuffle_randomstepconv/TRADEOFF_G4_G8.md`;
- `scripts/analyze_riffle_transpose_bitshuffle.py`;
- `scripts/certify_riffle_transpose_bitshuffle_full.py`;
- direct receipts `g8_b1024_sigma12_full_interval.json`,
  `g8_b512_sigma14_full_interval.json`, and
  `g8_b256_sigma17_full_interval.json` under the bit-shuffle construction;
- packet g=8 full receipts under the packet-shuffle construction.

Recommended next goal: search padded bit-transpose outer lengths between 256
and 512, then benchmark the g=4 and g=8 bit-permutation plus inner kernels.
The proof proxy favors g=8, but only a low-level benchmark can price the extra
bit-shuffle memory traffic.

## 2026-08-27 one-block tail sharpening

The complete bit-transpose certificates are dominated by one active outer
block. Their Chernoff bounds omit a lattice tail prefactor. The exact
one-block generating function is

```text
R(z) = (1/P) sum_j M0(z)^j M1(z) M0(z)^(P-1-j)
F(z) = e0^T R(z)^B 1.
```

The new script differentiates this transfer to locate the saddle. It also
inverts the optimally tilted generating function with a complex FFT. On three
scaled exact instances, the lattice estimate is within 0.031 bits. The FFT
reproduces one exact tail within 1.4e-14 bits. At full size, transforms of
lengths 2^20 and 2^21 agree within 6.4e-12 bits.

The sharper numerical frontier is:

```text
B=1024: g=4 sigma13 certified; g=8 sigma12 certified
B=512:  g=4 sigma14 lambda 40.2973 numerical; g=8 sigma13 lambda 40.3259 numerical
B=256:  g=4 sigma17 lambda 42.4010 numerical; g=8 sigma16 lambda 42.4043 numerical
```

The next lower states fail by wide margins. Full spectrum runs at the four
new cells remain dominated by one active block. The sums over all occupations
with at least two active blocks have exponents below -79 at B=256 and below
-129 at B=512.

The four new cells are not yet certificates because the FFT lacks an
outward-rounded error bound. The clean first certification target is
bit-transpose g=8, B=256, sigma16, which has 2.40 bits of room. B=512,
sigma13 has only 0.326 bits of room.

Artifacts:

- `constructions/riffle_transpose_bitshuffle_randomstepconv/ONEBLOCK_TAIL_SHARPENING.md`;
- `scripts/analyze_riffle_transpose_bitshuffle_oneblock_saddle.py`;
- one-block saddle and Fourier receipts under the bit-transpose construction;
- full-spectrum exploratory receipts for g4/B256/sigma17,
  g4/B512/sigma14, g8/B256/sigma16, and g8/B512/sigma13.

Recommended next goal: implement an outward-rounded upper bound for the
tilted Fourier tail at g=8, B=256, sigma16. Then combine it with a complete
outward-rounded sum over occupations a>=2.

## 2026-08-27 fixed-spectrum outer with bit transpose

The active outer-replacement candidate is **Riffle
SpectrumPerm-TransposeBitShuffle-RandomStepConv g**. Fix one binary
\([B,B/2]\) code \(C\) with weight spectrum \(A_w\). Reuse \(C\) in all outer
blocks. Setup samples an independent coordinate permutation for every outer
block. It then applies the existing row bit permutations and random inner
maps.

For a fixed message whose active outer words have weights
\((w_1,\ldots,w_a)\), each block permutation produces an independent uniform
row subset of the specified size. Thus the first-moment proof needs the outer
code only through \(A_w\). The block permutations must be independent; one
shared permutation would retain relative support structure.

Goal 01 handles one active outer block. If \(R_0(z)\) is an inactive-row
transfer and \(R_1(z)\) is an active-row transfer, the exact conditioned
matrix is

```text
[u^w] (R0(z) + u R1(z))^B / C(B,w).
```

At n=2^20, B=256, and distance threshold 188743, the ideal expected random
\([256,128]\) spectrum gives:

```text
g=4 sigma18: one-block lambda 46.6215; dominant weight 26
g=4 sigma17: one-block lambda 38.1852; dominant weight 35
g=8 sigma17: one-block lambda 46.6239; dominant weight 26
g=8 sigma16: one-block lambda 38.1886; dominant weight 35
```

The weight-26 expected multiplicity is below one. It represents rare bad
members of the random-code ensemble and disappears in a fixed code with
larger minimum distance.

A modeled complement-symmetric even \([256,128,38]\)-shaped spectrum gives:

```text
g=4 sigma18: one-block lambda 51.9477
g=4 sigma17: one-block lambda 40.5495
g=8 sigma17: one-block lambda 51.9530
g=8 sigma16: one-block lambda 40.5539
```

Weight 38 dominates. The higher-state cells tolerate about 12.25 bits of
extra weight-38 multiplicity if other modeled weights remain fixed. The
lower-state cells tolerate only about 0.68 bits. The spectrum is real-valued
and does not assert the existence of an explicit code.

Artifacts:

- `constructions/riffle_spectrumperm_transpose_bitshuffle_randomstepconv/`;
- `scripts/analyze_riffle_spectrumperm_bitshuffle_oneblock.py`;
- `GOAL_01_RESULTS.md` and Goal 01 receipts in the construction folder.

Recommended next goal: implement the exact two-active-block evaluator. Given
weights \((w_1,w_2)\), compress the two random row subsets by their
intersection size. Continue using the distance-38 model until an explicit
constituent spectrum is selected.

## 2026-08-27: first structured outer-and-inner design

The active family is **Riffle
BCHPerm-TransposeBitShuffle-StructuredStepConv(t,s)**. Its stages are:

1. a fixed binary \([256,128]\) outer code with a supplied spectrum;
2. an independent coordinate permutation on every outer word;
3. bit transpose into 256 regions;
4. an independent bit permutation inside every region; and
5. either FieldMulStepConv\((t,s)\) or ToeplitzStepConv\((t,s)\).

The inner step width \(t\) is independent of the state width \(s\). Both
structured inner families have the exact dense-random one-vector transition.
FieldMul uses \(t+s\) coefficient bits per step. Toeplitz uses
\(2(t+s)-1\).

For the modeled even, complement-symmetric \([256,128,38]\)-shaped spectrum,
the 9% one-active 40-bit boundary is:

```text
(t,s)=(4,17):   lambda 40.5495
(t,s)=(8,16):   lambda 40.5539
(t,s)=(16,15):  lambda 40.5634
(t,s)=(32,14):  lambda 40.5911
```

The preferred preliminary point is \((t,s)=(32,32)\). Its 9% one-active
margin is 80.8314 bits. Its one-active 40-bit distance frontier lies between
14.05% and 14.10%. These are numerical diagnostics, not complete distance
claims.

The transposed inner implementation processes \(2^{21}\) 128-bit elements
in reverse order. It uses an AVX2 128-by-128 transpose and two-lane
`VPCLMULQDQ`. On the Intel Core i7-13700H benchmark host, 15-trial medians at
\((32,32)\) are:

```text
FieldMul:  27.1398 ms
Toeplitz:  25.8285 ms
```

All benchmarks were serial. The implementation excludes the outer encoder
and both permutation layers. The C++ self-test matches dense references at
widths 31, 32, 46, 47, 48, and 64. A separate checker verifies every field
modulus and sampled one-vector ranks.

Artifacts:

- `constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/`;
- `scripts/sweep_bchperm_bitshuffle_stepconv.py`;
- `scripts/bench_transposed_structured_step.cpp`;
- `scripts/verify_structured_step_families.py`.

Recommended next proof goal: compute the multi-active outer-spectrum sum for
the bit-transpose geometry. The current parameter sweep covers exactly one
active outer block.

## 2026-08-27: physical-transpose implementation audit

The construction requires the transpose of the encoder as a linear map. It
does not require the inner kernel to physically transpose a 128-by-128 data
tile. The benchmark now contains two equivalent evaluators:

1. a bitsliced evaluator that physically transposes the tile and batches 128
   scalar carryless multiplications; and
2. a direct evaluator that keeps the 128-bit elements in place and applies
   the binary matrix with four-bit XOR-combination tables.

Both evaluators match dense reference maps and complete reverse chains at all
supported widths. The measured crossover is width-dependent:

```text
(t,s)       direct Field   bitslice Field   direct Toeplitz   bitslice Toeplitz
(16,16)       37.3897         56.4151          15.6509           48.8910 ms
(32,14)       30.4932         29.2352          17.9484           25.8625 ms
(32,15)       30.6888         29.4721          18.4905           25.9712 ms
(32,32)       42.7634         27.1398          34.3959           25.8285 ms
```

Therefore, a physical transpose is not the general answer. Direct Toeplitz
is best at the narrow widths near the 40-bit one-active boundary. Bitslicing
wins at width 64 because 128-way carryless-multiplication batching repays the
fixed transpose cost. The mathematical construction and distance model are
unchanged.

Artifacts:

- updated `scripts/bench_transposed_structured_step.cpp`, with
  `--implementation direct|bitsliced`;
- `constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/direct_kernel_comparison.json`;
- corrected implementation discussion in `CONSTRUCTION.md` and `RESULTS.md`.

Recommended next implementation goal: specialize direct Toeplitz for widths
46 and 47, then measure the complete transposed encoder, including the two
permutation layers. Recommended next proof goal remains the exact
multi-active outer-spectrum sum.

## 2026-08-27: block-resident Toeplitz algorithm comparison

The inner benchmark now includes three Toeplitz algorithms that preserve
128-bit block form throughout:

1. Four Russians with four-input XOR tables;
2. AVX2 XORs along the nonzero Toeplitz diagonals; and
3. compile-time Karatsuba polynomial convolution.

Every implementation matches the dense map at widths 31, 32, 46, 47, 48,
and 64. Complete reverse-chain tests pass at every supported parameter pair.
The benchmarks were run serially.

```text
(t,s)       width   Four Russians   diagonal AVX2   Karatsuba
(16,16)       32       14.8980          39.0523       36.7084 ms
(32,14)       46       18.9647          33.3514       55.0741 ms
```

Karatsuba cutoffs of 8, 16, and 32 gave 105.0004, 76.6637, and 55.0741 ms,
respectively, at width 46. Padding to 64 coefficients and materializing
temporary products dominates at this size. The diagonal method has
contiguous accesses but performs more block XORs. Four Russians remains the
fastest tested block-resident evaluator. Its repeated 18.9647 ms result is
consistent with the earlier 17.9484 ms measurement.

Artifacts:

- updated `scripts/bench_transposed_structured_step.cpp`, with
  `--implementation four-russians|diagonal|karatsuba|bitsliced`;
- `constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/block_toeplitz_algorithm_comparison.json`;
- updated construction and results notes.

Recommended next implementation goal: measure a no-map chain floor and
profile the Four-Russians kernel to separate table construction, indexed
loads, and state movement. Then optimize the dominant component or benchmark
the complete encoder if the inner is no longer expected to dominate.

## 2026-08-27: Four-Russians cost diagnosis

The width-46, ((t,s)=(32,14)) block-resident kernel was split into phases.
A recurrent movement floor touches every input block, carries the state, and
uses one block XOR per output instead of the dense Toeplitz map. Setup-row
expansion and sampled serialized timestamp profiling were also implemented.

```text
recurrent movement floor:                  5.7071 ms
compact Four Russians before local tuning: 19.2238 ms
compact Four Russians after local tuning:  18.5588 ms
pre-expanded-row Four Russians:            17.4525 ms
```

Pre-expanded rows occupy 24,117,248 bytes. They save 1.11 ms despite the
extra sequential metadata reads. The compact tuning streams the Toeplitz
tail bits and initializes outputs from their first table selection. A
partial-final-group specialization regressed to 22.6792 ms and was reverted.

The sampled compact-kernel profile reports:

```text
state assembly and commit:       53.2 ticks/step   6.2%
Toeplitz-row generation:        104.9 ticks/step  12.3%
table construction:             147.2 ticks/step  17.3%
indexed table application:      546.6 ticks/step  64.2%
```

One of every 256 steps was sampled. The minimum timestamp overhead was
subtracted. The profile is for attribution; the wall-clock benchmark remains
the authoritative timing. The result rules out benchmark overhead and state
movement as the main explanation. The dense block-XOR application dominates.

Artifacts:

- updated `scripts/bench_transposed_structured_step.cpp`, including
  `--implementation movement-floor|preexpanded-rows` and
  `--profile-phases`;
- `constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/four_russians_cost_profile.json`;
- updated `RESULTS.md`.

Recommended next implementation decision: either accept approximately
17.5 ms for this dense block-resident Toeplitz inner and benchmark the full
encoder, or change the inner family to one with a substantially smaller XOR
circuit. Micro-optimizing row generation cannot provide a large additional
gain because table application alone is about two thirds of the kernel.

## 2026-08-27: FieldCheckpointAccumulate proof exploration

The active low-cost inner candidate is now **Riffle
FieldCheckpointAccumulate t=32 s=64 K=32**. It retains the fixed-spectrum
outer, per-outer-block coordinate permutations, transpose, and independent
bit permutation within each of the 256 regions. The inner uses 64 parallel
accumulator lanes. After every 32 batches of 32 bits, it applies the transpose
of multiplication by an independent uniform nonzero element of
\(\mathrm{GF}(2^{64})\) to the state. The transposed evaluator therefore uses
one ordinary field multiplication per 1024 output bits.

The checkpoint makes the inner close exactly on two state classes: zero and
uniform nonzero. Exact polynomial formulas now compute the tilted epoch
transfer for every fixed impulse count from zero through 1024. Small cases
match exhaustive enumeration within floating-point roundoff. The target
\(s=64\), 1024-bit transfer is finite, nonnegative, and stochastic at
\(z=1\) to within \(1.2\times10^{-15}\).

For one active outer block and the modeled distance-38 spectrum, the margin
at 9% relative distance is 81.2667 bits. The 40-bit frontier is approximately
14.049%: margins at 14.049% and 14.050% are 40.0007 and 39.9944 bits. Similar
one-active margins across several state/checkpoint geometries indicate that
late placement of the outer support, rather than checkpoint extinction, is
the leading failure mode.

The multi-active proof remains open. The inner transfer for every fixed
region occupancy is now exact and two-dimensional. The remaining problem is
to average the ordered product of 256 region matrices over the correlated
region occupancies induced by several independently permuted outer words.

Artifacts:

- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/`;
- `scripts/analyze_riffle_fieldcheckpoint_accumulate_oneblock.py`.

Recommended next proof goal: build every exact region-occupancy matrix and
test a common positive-potential domination against exact one- and two-active
calculations. If it is tight, use it to scalarize the full outer-spectrum
sum.

## 2026-08-27: checkpoint-accumulator cancellation audit

The local weight-two concern is real. Two impulses in the same epoch-lane
return that lane to zero. Across a complete 8192-bit region, a uniform
weight-two support fails to activate with probability \(2^{-9.09}\). The
worst consecutive pair emits weight one and occurs with probability
\(2^{-12.09}\).

The outer transpose appears to suppress this mechanism globally. For two
independently permuted weight-38 outer words, the support intersection has
mean 5.64 and mode five. The modal configuration has 66 singleton regions;
every singleton activates the state. The first occupied region fails with
probability \(2^{-12.72}\). Nonactivation through region 210 has exponent
240.23 bits. After charging the modeled 73.55-bit multiplicity of all such
two-block messages, 166.68 bits remain before bounding the output tail.

This is an exact activation diagnostic, not a two-active distance
certificate. The next proof goal is the exact two-active Chernoff calculation
using \(R_0(z)\), \(R_1(z)\), and \(R_2(z)\).

Artifacts:

- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/CANCELLATION_ANALYSIS.md`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal03_cancellation_diagnostic.json`;
- `scripts/analyze_riffle_fieldcheckpoint_cancellations.py`.

## 2026-08-27: full-certificate route

The exact two-active calculation is complete under the modeled spectrum. It
sums all 4278 unordered outer-weight pairs and all support intersections. The
margin is 159.8280 bits, with \((38,38)\) dominant.

A spectrum-density envelope compresses arbitrary regular outer-weight tuples.
For every modeled support other than the all-one support, its expected mass
after coordinate permutation is at most twice the mass induced by a uniform
256-bit outer word. The factor-two loss applies per active block. The outer
sum therefore reduces to one regular active-block count \(a\).

The floating-point envelope sum for \(1\le a\le64\) has 55.9008 bits of
margin and is dominated by \(a=1\). At \(a=64\), the pointwise exponent is
about \(-3393\) bits. Exploratory values through \(a=128\) remain thousands
of bits below zero.

This is not yet a full certificate. The remaining ranges are regular
occupations above 64 and configurations containing one or more all-one outer
words. The final checker also requires outward-rounded scaled arithmetic.

Artifacts:

- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/proof/FULL_CERTIFICATE_PLAN.md`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/proof/SPECTRUM_ENVELOPE.md`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal05_two_active_fullspectrum.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal06_regular_envelope_a64.json`;
- `scripts/analyze_riffle_fieldcheckpoint_twoactive_fullspectrum.py`;
- `scripts/analyze_riffle_fieldcheckpoint_regular_envelope.py`.

Recommended next proof implementation: compute \(\overline R_a(z)\) with
mantissa-plus-exponent coefficients beyond occupation 128, then overlap that
range with an analytic bulk bound. Treat the all-one occupation by shifting
the region occupancy in the same checker.

## 2026-08-27: FieldCheckpointAccumulate middle-density result

The exponent-safe region checker is implemented in
`scripts/analyze_riffle_fieldcheckpoint_regular_bulk_logdp.py`. It computes
the exact hypergeometric recurrence across the eight epochs of a region in
the log semiring. It matches the earlier ordinary-arithmetic calculation on
the overlap range and matches a small direct reference within
\(10^{-15}\). The FFT bulk experiment is not used because its coefficient
floor is too large.

The full-occupation calculation changes the distance picture. At 9% relative
distance, 1024 regular active outer blocks have pointwise first-moment
exponent about \(+16039\) bits. The Chernoff optimum is inside the tested
grid. The spectrum-density envelope loses only about one bit per active
block, so exact spectrum accounting cannot plausibly recover the deficit.
This is evidence against the 9% first-moment proof route, not a refutation of
the code.

At 6% relative distance, every regular occupation closes. Occupations 1025
through 5466 have 834.4399 bits of aggregate margin. The low range inherits
the existing 55.9008-bit bound at 9%. Starting at occupation 5467, the
invertible-inner Hamming-ball bound has 107.1058 bits of aggregate margin.

The first exceptional all-one layer also closes. The identity

\[
 F_{a,1}=2F_{a+1,0}-F_{a,0}
\]

reproduces the independent low-occupation mixed enumerator and gives
1161.5591 bits of margin for 65 through 5487 regular blocks. Starting at 5488
regular blocks, the dense bound sums every possible all-one count and retains
91.3725 bits.

Artifacts:

- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal10_regular_distance_scan_a1024.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal11_regular_middle_logdp_delta06.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal13_one_allone_fullmiddle_delta06.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/receipts/goal14_dense_volume_delta06.json`;
- `scripts/analyze_riffle_fieldcheckpoint_dense_volume.py`.

Recommended next proof goal: compute the \(b=2,3,\ldots\) all-one layers
from \(F_{a,b+1}=2F_{a+1,b}-F_{a,b}\), validate them against the direct mixed
enumerator, and look for a monotone tail bound in \(b\). The only remaining
combinatorial range for the conditional 6% certificate is \(b\ge2\) with
fewer than 5488 regular blocks. Outward rounding follows after that range is
closed.

## 2026-08-27: frequent refresh recovers the 9% regular bound

The distinct candidate **Riffle FieldCheckpointAccumulate t=32 s=64 K=8**
refreshes the 64-bit state every 256 output bits. The parent K=32 candidate
refreshes every 1024 bits.

The intermediate K=16 experiment refreshes every 512 bits. It improves the
pointwise exponent at 1024 regular active blocks from a 16039-bit deficit to
6706 bits of room. The obstruction moves to larger occupations. At 2048
active blocks, K=16 still has a 16750-bit deficit.

K=8 removes the regular middle-density obstruction at 9% relative distance.
The modeled-spectrum density envelope gives the following aggregate margins.

| Regular active blocks | Aggregate margin |
|---:|---:|
| 1–511 | 55.9507 bits |
| 512–4096 | 19399.3836 bits |
| 4097–7241 | 28938.7956 bits |
| 7242–8192 | 168.1408 bits |

The first three rows use the exact log-domain region recurrence. The final
row uses invertibility and Hamming-ball volume. The weakest middle point is
near 3240 active blocks and retains approximately 19399 bits. The complete
regular sum is dominated by one active block.

The result is conditional on the modeled outer spectrum and uses nearest
binary64 arithmetic. The unique all-one support remains outside the density
envelope. The dense bound covers every all-one count once at least 7249
regular blocks are active. The other all-one configurations remain open.

The shorter epoch requires approximately four times as many field
multiplications as K=32. Accumulator work is unchanged.

Artifacts:

- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k8/`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k8/receipts/goal01_regular_middle_delta09.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k8/receipts/goal02_regular_low_delta09.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k8/receipts/goal03_regular_high_delta09.json`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k8/receipts/goal04_dense_volume_delta09.json`.

Recommended next proof goal: keep the 9% target and evaluate the all-one
layers for K=8. Start with \(b=1\), which follows from
\(F_{a,1}=2F_{a+1,0}-F_{a,0}\). Then extend to \(b\ge2\) and search for a
monotone tail after charging all-one placement multiplicity. Outward rounding
follows after the exceptional support closes.

## 2026-08-27: four visits per lane, not K alone

Two wider-state variants test whether K=8 is intrinsically necessary. Both
variants preserve four accumulator visits per state lane:

| Candidate | Mixing work per 1024 output bits | 9% middle margin |
|---|---:|---:|
| \(s=64,K=8\) | \(4M_{64}\) | 19399.3836 bits |
| \(s=128,K=16\) | \(2M_{128}\) | 21135.6506 bits |
| \(s=256,K=32\) | \(M_{256}\) | 21609.2409 bits |

Each calculation covers 512 through 4096 regular active outer blocks. The
same modeled spectrum, exact log-domain recurrence, Chernoff grid, and 9%
threshold are used. All three candidates have ample margin. The result
identifies visits per lane as the important distance parameter. It permits
the original 1024-bit checkpoint interval when the state has 256 lanes.

The current benchmark has no optimized transposed field kernels above width
64. Performance therefore remains undecided. A superlinear circuit favors
four small maps, while fewer transposes can favor a wider map.

Artifacts:

- `explorations/riffle_fieldcheckpoint_four_visits_landscape.md`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s128_k16/`;
- `constructions/riffle_fieldcheckpoint_accumulate_t32_s256_k32/`.

Recommended next implementation goal: add optimized width-128 and width-256
transposed field kernels and benchmark \(4M_{64}\), \(2M_{128}\), and
\(M_{256}\) in isolation. Do not run those benchmarks concurrently. The proof
work should remain on the 9% target.
