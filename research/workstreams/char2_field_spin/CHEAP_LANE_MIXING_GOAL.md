# Cheap lane mixing: goal and proof ledger

## Active objective

Explore the [implemented b(theta+a) refresh](CHEAP_LANE_MIXING_20261016.md)
while keeping the certified 3.8-ms K20, rate-1/2 design unchanged.
First seek a complete 35% certificate for the modified code; then climb
above 35%, with 40% as the next checkpoint and distance near 50% as the
longer-term objective. The working targets are at least 40 bits of
ideal-setup failure margin and roughly at most 5% matched timing overhead.

These are research targets, not established guarantees. The goal remains
active. There is no new whole-code certificate, production promotion,
paper edit, commit, or push.

## Progress, 2026-10-17

| Part | Current result | Remaining work |
|---|---|---|
| Implementation | Independent scalar checks; about 0--3% measured overhead | Rebenchmark only after a relevant implementation change |
| Local refresh lemmas | Independently audited; correct with fixed-message conditioning | Build a joint continuation bound |
| True all-quiet histories | New directed 40% bound, more than 76 bits | Closed at the next distance checkpoint; preserve receipt |
| One actual active source group | Directed 40% bound, more than 44 bits, all lane ranks | Closed for all 8192 group choices; extend actual-support accounting |
| Two through six actual active groups | Directed 35% bounds: more than 450, 250, 100, 50, and 45 bits | Preserve these stronger small-class bounds |
| Seven through 32 actual active groups | Directed 35% joint bound: at least 178 bits, all supports and trajectories | Extend the size-dependent proof witnesses beyond this band |
| 33 through 224 actual active groups | Directed 35% band bounds, each exceeding 64 bits, with no gaps | Preserve these stronger nested band receipts |
| Combined small-source cover | Every support of size1--751, all paths, at35%: more than42 bits jointly; nested1--507 exceeds43 bits | Extend without spending the remaining whole-code failure budget |
| Selected 512-group profile | Directed35% bound:107.5446806 bits for (16,480) in both matchings | Preserve the stronger restricted-class bound |
| Clustered 512-group profile | Exact regional coefficient closes (223,66) in both matchings at35% with82.3387645 bits; directed256-/384-bit scale-invariance replay | Reproduced inside the complete-size calculation |
| Every 512-group support | Directed35% bound:60.8663748 bits over all66048 ordered profiles and every state path | Extend shared-table coverage to neighboring source sizes |
| Every support of sizes508--515 | Directed35% joint bound:43.5679475 bits, all profiles and state paths | Preserve the stronger60.8664-bit bound at512; test broader ranges next |
| Every516-group support | Directed35% bound:1205.5528627 bits within the256-bit coefficient replay; earlier current-pair singleton357.8913170 bits | Preserve stronger singleton evidence within the new band |
| Every support of sizes508--548 | Directed35% joint bound: more than43 bits; new516--548 band alone842.2922173 bits | Extend beyond this band without extrapolating |
| Every1024-group support | Directed35% bound:430.0100235 bits, all263168 ordered profiles and every state path | Extend the improved pair operator to broad source-size intervals |
| Every2048-group support | Directed35% bound:86.7424999 bits, all1050624 ordered profiles and every state path | Preserve this stronger singleton bound alongside joint bands |
| Every support of sizes2047--2049 | Directed35% joint bound:74.2104338 bits, all1050625 ordered envelope pairs and every state path | Preserve this stronger nested-band receipt |
| Every support of sizes2040--2056 | Directed35% joint bound:60.2491326 bits, all1058840 ordered envelope pairs and every state path | Preserve stronger nested bounds; investigate the larger-support regime |
| Every3072-group support | Directed35% bound:3356.7491176 bits, all2362368 ordered profiles and every state path | Preserve this stronger singleton receipt within the new continuous band |
| Every support of sizes3008--3392 | Directed35% joint bound: more than196 bits, combining four disjoint bands | Reuse verified regional-bound banks; recalibrate beyond their useful range |
| Every3264-group support | Directed35% bound:3473.8373402 bits from the recalibrated residue bank | Preserve this stronger singleton control without adding it to the band union |
| Every4096-group support | Directed35% bound:4621.1552024 bits; cheaper ordinary marker envelope:256.2333998 bits, both complete | Preserve stronger receipt; extend cheap verification toward size bands |
| Every6144-group support | Directed35% ordinary-Cauchy bound:49748.3343231 bits, all1050624 ordered profiles and every state path | Preserve the stronger singleton receipt within the new joint band |
| Every support of sizes5632--6528 | Directed35% joint bound:4134.3211433 bits, all897 integer sizes and every state path; nested5888--6400 bound:25386.7993133 bits | Recalibrate the proof witnesses outside this band; no construction change |
| Every4352-group support | Directed35% bound:497.8311219 bits | Preserve the stronger singleton receipt within the new band |
| Every support of sizes4352--4383 | Directed35% joint bound:221.7124546 bits, using cell-specific run counts | Bridge the remaining gaps without extrapolating |
| Every support of sizes4544--8192 | Directed35% joint bound: more than780 bits, combining nine disjoint bands including the dense bound | Closed for all3649 sizes; bridge the lower gap without extrapolating |
| Dense actual supports | Directed 35% bound for 7552--8192 groups, more than 781 bits jointly | Retain output-image and feedback ranks below this band |
| Proper input-lane spans | Exact activation and rare-return bounds below | Count histories and messages; cover terminal tails |
| Full-rank inputs | One-group class closed; audited general rank interface | Cover nonquiet multi-group cases without losing source rank |
| Short quiet excursions | Directed 76-bit bound for all size-256 supports under the 44-packet-window restriction; consecutive-size 196--256 tail has 45 bits | Retain these restricted classes; close complementary trajectories |
| Quiet-heavy unrestricted paths | Four actual size-256 profiles have floating log2 union -510.227 at births<=512, quiet steps>=24576 | Cover complementary branches before expanding profiles or directed verification |
| All trajectories, every size-255 or size-256 support | Directed 35% bounds:105.2677757 and93.5574055 bits respectively | Cover the other source sizes; no profile restriction remains at255 or256 |
| All trajectories, every size-257 support | Directed 35% bound:94.287033 bits after audited source-factor improvement | Use joint missing-group comparison at larger sizes |
| Whole code | Not certified | At the40-bit target, nonquiet sizes752--1023,1025--2039,2057--3007,3393--4095,4097--4351,4384--4543 remain open |

Three independent read-only reviews checked local algebra, quiet histories,
and full-rank inputs. The local audit found no counterexample, but clarified
the degree-window initialization and the order of conditioning. Those
corrections are in the construction note. All bounds concern the actual
limited-support refresh coefficients, not uniform GF128 coefficients.

## Retaining D1 randomness gives a stronger local bound

Use E=GF256 and F=E(theta) as defined in the construction note. Fix a
proper E-linear subspace U of F and a message whose input stream lies in U.
Condition on the outer/route and all earlier refreshes, but not the current
independent nonzero diagonal D1. Suppose the current quotient state has
support h>0. Set V=Q+CX and W=H D1 V. V has the same quotient support h.

For a set J of desired zero quotient coordinates, choose an E-linear
functional ell on F/U that detects V. Choose ell before sampling D1.
Let I be the support of the scalar vector ell(V). The event
W[J] in U^J implies

    H[J,I] diag(ell(V)[I]) d[I] = 0,

where d is the diagonal of D1. Superregularity gives rank
min(|J|,|I|). If |I|<=|J|, the event is impossible because d[I] has no
zero coordinates. Otherwise, condition on nonpivot diagonal entries;
the remaining |J| entries have at most one admissible assignment. Thus

    Pr[W[J] in U^J] <= 255^(-|J|).

Write g=wt(W mod U). The MDS support bound gives g>=17-h. A union over
zero-coordinate sets therefore gives

    Pr[g(W)=g] <= binom(16,g) 255^(-(16-g)).

Now expose W, leaving the new independent a bytes fresh. Each of its g
outside coordinates can return to U for at most one byte value. Therefore

    Pr[Q_next in U^16]
      <= sum_{g=17-h}^{16} binom(16,g) 255^(-(16-g)) 256^(-g)
      <= (1/255 + 1/256)^16 - 255^(-16)
      =: epsilon_return.

This upper bound is approximately 1.987216721030805e-34, or
2^(-111.9548060081843). For h=1, the sharper bound is exactly 2^-128.
The current D2 entries can be exposed throughout: nonzero E scaling
does not change membership in U.

This bound averages D1 instead of first fixing its worst realization.
No additional mixing or implementation work is needed.

### A tail bound for quotient support

Let eta=1/255+1/256. Fix a set J of m coordinates required to land in U
after refresh. Partition J into R, the coordinates outside U before the
new multiplier, and J minus R. The latter impose |J minus R| quotient
zero constraints on W. The former require |R| exceptional a bytes.
Dropping the exact-status constraints only enlarges the event, giving

    Pr[Q_next[J] in U^J]
       <= sum_{R subset J} 255^(-|J minus R|)256^(-|R|)
       = eta^m.

Consequently, still conditional on a nonzero incoming quotient state,

    Pr[h_next <= j] <= min(1, binom(16,j) eta^(16-j)).

For j=0, use the slightly sharper epsilon_return above. In particular,
Pr[h_next<=12] <= 6.833361642673765e-6. If h_next>=13, all four
RS[16,4] state blocks are nonzero modulo U, so the next emitted packet
has at least 52 nonzero records. Summing the four support-tail bounds
gives an expected next-packet lower bound of 51.9999111663 out of 64.

This is an activated, fixed-message local result. It does not mean the
code has 81.25% distance. A codeword can activate late, remain quiet,
or have full lane span; the whole-code union is still absent.

For a fixed stream, a union over at most 32767 transitions before the
first return gives probability at most 6.511513029801638e-30 of such a
return. That union is over time, not over messages or subspaces.

The small [exact-rational checker](theory/cheap_lane_refresh_bounds.py)
verifies the closed forms, monotonicity, target-set binomial expansion,
and displayed numerical ceilings. It uses no sampled trajectories and
writes no files. It verifies arithmetic, not the entire coding theorem.

## True quiet histories transfer exactly

Fix the E-valued outer encoder O and paired route P. Let R=N/64=32768.
Let M impose C X[t]=0 for t=0,...,R-2 on X=P O m. The last state update
does not affect any emitted output, so it need not be constrained.
M is a matrix over E. Define the quiet message space over a field L by

    V_L = ker_L(M),  L in {E,F}.

Then V_F is the scalar extension of V_E. For every quiet message, all
incoming states are zero, and the emitted word is P O m for either the
old refresh or the new one. This holds for every realization of refresh
randomness.

To compare distances, expand m in an E-basis of F. Every component lies
in V_E. A nonzero component produces a nonzero outer word whose support
is contained in the F-word's support. Conversely, any E quiet word embeds
in F. Therefore the two quiet-subcode distances are exactly equal.
There is no factor-sixteen or subspace-counting loss.

If an outer/route setup has a quiet word of weight at most 734003, the
old scalar encoder is bad for every old refresh setup. Hence the old
[complete certificate](HD_SHIFTED_CLOSED_CERTIFICATE.md) implies

    Pr_outer,route[quiet distance <= 734003]
        <= 2^(-55.092738508223).

The same bound holds for the new construction's quiet bad event. This
reuses the preserved complete certificate; it is not a fresh coefficient
replay and does not certify nonquiet messages.

True zero-state prefixes also have E-defined constraint matrices, whose
ranks are unchanged by extension to F. A component can preserve the first
activation time and cannot increase prefix support. This does not bound
the later mixed suffix. For higher-dimensional U, quotient-quiet states
inside U^16 need not be zero; do not conflate those with V_L.

### An actual distance ceiling

Restrict messages to any 4096 outer groups. There are 4096*128=524288
source coordinates over E and at most 32767*16=524272 quiet equations.
Thus the restricted quiet subcode has E-dimension d>=16 and support
contained in at most n=4096*256=1048576 output coordinates.

For a uniform nonzero word of an E-linear d-dimensional code, each
nonzero coordinate functional is nonzero with probability

    (256^d - 256^(d-1)) / (256^d - 1).

The average word weight is therefore at most

    n * (255/256) / (1-256^(-d)) < 1044481.

Some nonzero quiet word has weight at most 1044480. For every setup,
the candidate consequently has relative distance at most

    1044480 / 2097152 = 255/512 = 49.8046875%.

This is an actual limitation of the unchanged GF256 quiet subcode, not
a failure of the new proof method. It leaves room for the 40% checkpoint
and for distance near, but not arbitrarily close to, 50%.

Full-rank inputs are not exempt from the quiet obstruction. Choose 16
E-independent quiet messages and combine them with coefficients
1,theta,...,theta^15. Their coordinate matrix has rank 16, so the resulting
message has full E-lane span. It remains quiet and is supported within
n=N/2 outputs. This establishes a full-rank witness at the weaker 50%
ceiling; the averaging witness at 49.8046875% need not have full rank.

## Next proof targets

1. Repair the exact lower-domain cover beyond751 toward1024.
   Keep256-bit coefficient rounding and retune the proof witnesses when
   a completed band loses margin. Preserve indexed counts and the output-
   only terminal, and check every size; witness proposals do not supply
   coverage before verification. The648--776 aggregate fails despite a
   fully completed passing prefix648--751. First diagnose its strong
   odd/even rounding pattern with the matched size776/run55 tail control.
   Preserve the45-bit new-band budget gate and fresh own-size numerators.
2. Bridge the remaining intervals using recalibrated proof witnesses and
   verified regional maxima. Complete source-size certificates now exist
   at512,1024,2048,and3072; the goal is coverage between them, not another
   selected-profile improvement. Retain every size's own source factor
   and both ordered terminal chains. The failed2944 and3393--3456 bounds
   identify limits of the current witnesses, not codeword counterexamples.
   For1025--2039, first generalize the existing17-size counted-band driver
   with explicit witnesses. Its live residue tables already support fresh
   actual-size targets; do not replace them by universal residue maxima.
3. Only after a complete35% cover, climb toward40% and higher. Local
   quotient-state bounds and restricted-support certificates do not
   imply the whole-code claim.

The outer's GF256 Toeplitz entries do not become uniform F entries.
Combinatorial routing counts and E-matrix rank identities transfer, but
scalar output densities cannot be replaced by the corresponding formula
with field size 2^128. Rank-sensitive outer bounds are a separate task.

For the quiet-only calculation, a smaller scalar operator is available.
For l potential positions in one emission block, its exact quiet kernel is

    k_l(z) = sum_{c=0}^l binom(l,c) z^(l-c) (1-z)^c
                         256^(-c-min(4,l-c)),  0<z<1.

This follows by expanding z^wt as a sum over additional zero-coordinate
constraints. It equals 256^(-l) times the weight enumerator of the
shortened feedback kernel. A packet contributes the product for its four
block occupancies. The final packet has no feedback constraint and uses
(z+(1-z)/256)^(total occupancy) instead. The paired-route support counting
can then remain entirely scalar. The second round below completed this
quiet-only checkpoint at 40%, with the free final packet explicit.

## Second round: the quiet class closes at 40%

The new calculation proves, for the unchanged outer and actual paired
route,

    Pr[there is a nonzero true-quiet word of weight <= 838860]
      < 8394752 / 2^100 < 2^(-76).

Here 838860=floor(0.40*2097152). Probability is over ideal independent
outer mixers, regional pair permutations, and fair pair swaps. The event
does not depend on refresh coefficients. The exact scalar-extension
argument above therefore transfers this bound to the GF128 quiet subcode.
This is not a certificate for nonquiet words or for a benchmark seed.

The [bounded helper](theory/cheap_lane_quiet_screen.py) and
[selected receipt](theory/cheap_lane_quiet40_cover.json) retain the proof
parameters. Fresh optimizer-free replays at 192 and 256 bits both passed.
Each checked 24 disjoint source-count bands, 104 clipped polygons, all
8192 integer source counts, and 413 directed vertices. Every vertex is
strictly below the required -100-bit cap. The number of feasible nonempty
histograms is binom(4098,2)-1=8394752; multiplying by this number gives
the stated rational union bound.

Source SHA256:
`426f13d3198f5913d9385839dd4333e8b7afdb8ca7811871fe3649c263455cf8`.
Receipt SHA256:
`d8c9b03a142b08f1d8bdac31e5e177552f5bfbd14357d573ae540231a284f360`.

From the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_quiet_screen.py --verify-receipt research/workstreams/char2_field_spin/theory/cheap_lane_quiet40_cover.json --precision 256
```

The verifier reads at most 256 KiB, checks the fixed scope and rational
geometry, and recomputes all vertex bounds. It does not invoke the
optimizer or trust recorded floating bounds or hashes. Hashes identify
the selected evidence; they do not replace the mathematical argument.

### Why the quiet calculation is small

The formula k_l above is the uniform-input expectation of z^wt(X) times
the zero-feedback indicator. Its rank min(4,l-c) follows from the RS
evaluation matrix, including the implemented coefficient basis. It is
not the invalid claim that every row/column submatrix is superregular.

A half packet has sixteen routed pairs feeding two such blocks. If it
contains s singleton pairs and f full pairs, fair swaps put u of the
singletons in the first block with probability binom(s,u)/2^s. Define

    P_(z,alpha)(x,y) = sum_(s+f<=16) multinomial(16-s-f,s,f) x^s y^f
       * sum_(u=0)^s binom(s,u)/2^s * (k_(f+u)(z) k_(f+s-u)(z))^alpha.

A region contains 4096 pairs, or 256 half packets. For its histogram
(e,s,f), its fractional quiet mean R is the coefficient of x^s y^f
in P^256 divided by multinomial(4096;e,s,f). This is exact for uniform
pair permutations and independent swaps.

Fix a number g of active source groups, and use the same z and alpha
for every support of that size. Let R_0(S),R_1(S) be the regional means
for the two adjacent-pair matchings. Cauchy--Schwarz gives

    sum_(|S|=g) R_0(S)^128 R_1(S)^128
      <= sum_(|S|=g) R_0(S)^256.

The two marginal sums are equal by cyclically relabeling the source
groups. This step does not assume independent matching geometries.
In the right-hand sum a histogram with s=g-2f and e=4096-g+f has
multinomial(4096;e,s,f)*2^s supports.

Use the proved scalar outer density beta=256^160/255^32. The final
packet has no feedback constraint; k_l(z) is at least 256^(-4) times
the corresponding unconstrained moment. Removing all four final-block
constraints therefore costs at most 256^(16 alpha). Charge this before
applying the matching reduction. The resulting quiet bound is

    (beta^g * 256^16 / (255*z^L))^alpha
       * sum_f multinomial(4096;e,s,f) * 2^s * R(g,f)^256.

For positive coefficient markers x=2^a and y=2^b, each summand has
base-two logarithm at most

    alpha*(g log2(beta)-log2(255)-L log2(z)+128)
      + s - 255 log2(multinomial(4096;e,s,f))
      + 65536 log2(P_(z,alpha)(2^a,2^b)) - 256*(s*a+f*b).

At fixed z,alpha,a,b, this is convex in (g,f): its only nonlinear
histogram term is a positive sum of log-Gamma functions. Exact rational
polygon vertices therefore bound their interiors. Markers can vary
between polygons; z and alpha stay common throughout each complete
g band. The verifier enforces that distinction and checks continuous
f-interval coverage on every integer g slice.

### What this round changed for nonquiet words

The rank interface now includes a global reduction to scalar constraint
counts. It uses the actual coefficient set and covers full E-lane-rank
inputs. A positive zero-subset expansion removes an unnecessary
output-support completion factor.

However, the first raw scalar-shadow operator admits artificial state
histories. Even one scalar source-group contribution has a relaxed
expression of at least about 2^1916.98, for every allowed fractional
exponent and Chernoff mark. An elementary exact inequality proves a
floor above 2^1015. This is a failure of the upper-bound method, not
a low-weight word. A large coefficient search for that operator would
not help, so none was launched.

A support-preserving refinement is available for fewer than 256 actual
active groups, and is exact for one group. It permits legitimate empty
input/state restrictions and adds at most eight bits of transfer loss
through 255 groups. For larger supports, a bounded-missing-groups scalar
shadow is another proved interface. The fourth round below closes a
dense band using input-column deletion; a tighter intermediate-support
count remains unfinished. No encoder change is justified by the current
obstruction.

## Third round: one actual active group closes at 40%

The [density-free rank argument](CHEAP_LANE_RANK_INTERFACE.md) now proves

    Pr[there is a nonzero one-group word of weight <= 838860] < 2^-44.

The directed base-two bound is -44.29481729861079297183875309738. It
includes all 8192 group choices and every GF128 message within a group,
including messages with full GF256 lane span. Independent 192- and
256-bit replays passed; the main agent also replayed the 256-bit check.
The simpler 35% checkpoint gives more than 180 bits. These are bounds
for the actual limited-support refresh distribution, not uniform GF128
coefficients or a transfer of the old encoder's certificate.

The important change is to the proof. A wholly quiet outer column forces
four independent RS equations. A partly quiet column supplies up to three
more. Preimposing those deterministic constraints avoids counting scalar
shadows that lose the actual source support. The rank potential can use
tau=65280 directly, without pretending there is a field of that size.

At 40%, the loose operator was dominated by late activation after 31
quiet columns and three quiet regions. Charging the three additional
source equations removes almost 48 bits of slack. Averaging the exact
activation and extinction restrictions supplies the remaining tightening.
No larger state, new outer, extra randomizer, or additional encoding work
was needed.

The source is [cheap_lane_sparse_gate.py](theory/cheap_lane_sparse_gate.py),
16,199 bytes, SHA256
`0de7630749a9f06b921e8faec8b894d12ce1a7efaf623bc6890b9d2a10d1ba2a`.
Its fixed-parameter replay is:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_sparse_gate.py --method density-free --verify --tau 65280 --z 0.9968 --distance 0.4 --activation-d1 --extinction-position --partial-quiet --precision 256 --target-bits 44
```

The floating proposal was slightly optimistic because small intermediate
terms underflowed. Only the directed result is used above. The failed
scalar-shadow mode remains explicitly diagnostic; it is not the proof
used by this checkpoint.

The new one-group event and the all-quiet event may overlap, but their
union can be bounded by adding their probabilities. This still leaves
nonquiet messages with two or more active groups uncovered at 40%.
There is no whole-code certificate yet.

### Two through four active groups also close at 35%

The next extension replaces whole-column counting by the best thirty-two
columns' quiet-region ranks. It averages the actual paired routing with
a small matrix polynomial. The proof and its conditioning are in the
[rank interface](CHEAP_LANE_RANK_INTERFACE.md). All messages over F are
covered within each stated source-count class; there is no lane-rank
restriction.

| Active groups r | Four fixed marks z_0,z_1,z_2,z_3 | Worst log2 union bound | Verified margin |
|---|---|---:|---:|
| 2 | .993, .994, .996, .997 | -463.0965409282736 | >450 bits |
| 3 | .993, .994, .996, .998 | -264.5425512849358 | >250 bits |
| 4 | .993, .995, .9965, .99835 | -110.5164679537663 | >100 bits |

All use delta=7/20 and tau=65280. Each bound already includes
binom(8192,r) source supports and the worst fixed pair-matching profile.
Independent 192- and 256-bit directed replays passed for every row.
The final mathematical and source audit found no blocking issue.

Source: [cheap_lane_few_group_gate.py](theory/cheap_lane_few_group_gate.py),
21,479 bytes, SHA256
`a5fb6a71a81d2652b9127c199d8423be845363fbc15de0b847771e31ac1426c7`.
For example:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_few_group_gate.py --groups 4 --tau 65280 --threshold-z .993 .995 .9965 .99835 --rank-mode top32 --activation-d1 --distance 7/20 --verify --precision 256 --target-bits 100
```

For r=2 or 3, substitute the corresponding marks and target margin from
the table. These runs are bounded, write no result files, and invoke no
optimizer. At 192 bits they took approximately 4.1, 4.9, and 9.2 seconds.

The remaining 35% gap is nonquiet messages with at least five active
groups. The current quiet-region rank formula cannot simply be extended:
five inputs in one RS feedback block can have a nonzero kernel. The next
operator should retain actual feedback rank and the deterministic input
constraints imposed at activation. The audited quantity is
rank[M,N]-rank(M), not an assumption that a quiet packet has zero input.
This identifies the proof refinement used in the fourth round below.

## Fourth round: five groups and a dense band close at 35%

The new rank-debit operator handles quiet packets whose allowed inputs
can cancel. It retains c=sum_b min(4,ell_b) feedback ranks, and at every
state transition retains the lower bound max(0,c-h-h'). The exact paired
route averages those ranks together with the local transition caps.
The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) gives the argument.

For exactly five actual active source groups, the resulting bound is

    Pr[there is a nonzero word of weight <= 734003] < 2^-50.

The directed log2 bound is -51.30995315679061752812820293990353. It
already includes binom(8192,5) source supports and the worst of all nine
fixed matching profiles. All twenty source-rank threshold classes are
summed. Thresholds 0--18 use a positive cutoff marker; threshold 19
retains the exact thirty-two-counter cutoff. The final cutoff was
necessary: the initial marker-only screen left a positive log2 bound.

Source: [cheap_lane_rank_debit_gate.py](theory/cheap_lane_rank_debit_gate.py),
19,665 bytes, SHA256
`04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053`.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_rank_debit_gate.py --verify-five --precision 256 --target-bits 50
```

An independent source audit passed. Fixed 192- and 256-bit replays agree;
the main agent independently replayed the 256-bit result. No floating
proposal value enters the verifier.

The other new certificate covers the union of all actual source-count
classes from 7688 through 8192, with no restriction on GF256 lane rank:

    Pr[there is such a nonzero word of weight <= 734003] < 2^-98.

Its directed log2 bound is -98.62666920007299253281040466936204. The
proof uses the actual 65280-element refresh distribution, not uniform
GF128 randomness. First it charges missing scalar-shadow groups by
matrix-column deletion. Then it pads the remaining inputs to full
packets and bounds the seventeen-state continuation by an explicit
positive vector. The factor for all source supports is included.
Only 32767 refresh transitions are charged; the final packet is free
of a next-state constraint.

Source: [cheap_lane_dense_band.py](theory/cheap_lane_dense_band.py),
6,266 bytes, SHA256
`4d220621ab64b431f2e29c804ccf6655a49f2747e78ed6401571bb3091df46e7`.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_dense_band.py --precision 256
```

The dense proof and source passed independent review. Both directed
precisions passed, followed by the main agent's independent 256-bit
replay. The helper contains the complete fixed expression and its
derivation; it writes no result data and performs no parameter search.

The remaining 35% gap is nonquiet messages with 6--7687 active groups.
A bounded six-group floating screen does not close that class: nearly
singleton matching profiles suggest roughly 41--42 bits, but pair-rich
profiles fail at the chosen marks. The worst tested log2 bound was about
+1703, dominated by threshold zero rather than the last cutoff. This is
an upper-bound diagnostic, not a low-weight codeword or an impossibility
result. A focused local-cap and parameter diagnosis is preferable to
blindly extending the sparse sweep.

The focused follow-up tested four z marks and four cutoff marks on
profile (2,3). Its best bound remained about +1121 bits. Enforcing the
exact threshold-zero minimizer condition made no material difference.
The dominant relaxed paths activate from clustered pairs and return to
zero while paying almost no deterministic input rank. This identifies
lost source constraints, not an observed codeword.

Two independently audited local refinements are now in the rank interface.
The first averages activation feasibility jointly with the rank it imposes,
removing the projective-input union factor when that rank is charged.
The second retains the weight of an affine feedback coset in a scalar
return bound. Neither has yet been integrated into a closing middle-range
operator. These are proof changes; they require no encoder work.

Bounded scalar-kernel diagnostics at r=256, 1024, and 4096 also remained
positive. Keeping a seventeen-state matrix improved on a common-vector
row bound, and the affine-coset refinement improved the r=4096 test by
about 60000 bits. None closed. These were floating repeated-full-pair
stress profiles, not a directed cover of actual alternating profiles.
The fixed marker grid does not prove an obstruction; under fractional
power alpha, the natural quiet-pair compensation is q^(2 alpha).
Do not launch a broad search or change the code based on this grid.

The next round below integrates joint activation ranks. Exact extinction
geometry remains a possible tightening for larger support classes.

## Fifth round: six groups close at 35%

The joint activation-rank operator now covers every nonzero message with
exactly six actual active source groups, over the full GF128 alphabet:

    Pr[there is such a word of weight <= 734003] < 2^-45.

The directed log2 bound is
-46.2061996903376335129611596212772137709738474612436454363698208581421173,
with interval radius below 5.02e-71 at 256-bit precision. It includes
quiet histories, all binom(8192,6) supports, all fifteen fixed matching
profiles, and all twenty-four source-rank thresholds. The worst matching
profile is (0,3). The probability is over the ideal independent setup
distribution implemented by the construction, not a uniform GF128 refresh.

The improvement charges each activation's deterministic source rank
together with its probability. A four-block rank floor excludes further
impossible activations. Live transitions use a conservative quotient-rank
bound. The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) derives these
positive rank tags and their conditioning. No encoder work was added.

The verifier fixes tau=65280 and delta=7/20. For thresholds v=0,...,22,
except v=12, it uses z=1-27(768-32v)/2000000 and cutoff marker 1/sqrt(2).
Threshold 12 uses z=3981/4000 and marker 1. Threshold 23 uses z=3997/4000
and the exact thirty-two-counter cutoff, including its terminal factor.
The two largest terms are thresholds 23 and 22, with log2 contributions
approximately -47.07386119 and -47.35191499, respectively.

Source: [cheap_lane_joint_activation_gate.py](theory/cheap_lane_joint_activation_gate.py),
24,041 bytes, SHA256
`da67674f450b6aef3a057a0e1431b6fe61314009fcf91ae196bc62d222982e6d`.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_joint_activation_gate.py --verify-six --precision 256 --target-bits 45
```

Independent 192- and 256-bit directed replays passed. The main agent
independently replayed the 256-bit result in about 37 seconds. A separate
mathematical and source audit found no blocking issue. The fixed verifier
does not use floating proposal values or run an optimizer.

An earlier half-total-rank screen failed even after excluding wholly
quiet histories. It lost the exact exclusion of nonzero messages after
thirty-two full-rank columns. Its positive log2 bounds were overcounts,
not codewords. Restoring that cutoff and retaining activation-rank tags
closed the class. A final-threshold minimizer-marker test at 2, 4, and 16
worsened the bound; none of those choices enters the certificate.

The remaining 35% gap is nonquiet messages with **7--7687 active groups**.
There is still no whole-code certificate. At 40%, only the all-quiet and
one-group classes are closed. The next sparse refinement should retain
the shared RS points of full pairs when computing extinction feasibility
and rank. A bounded local example improves that transition by about
20.094 bits; forgetting the shared points retains only about 2.006 bits.
This example is not a new class certificate. For intermediate supports,
the separate priority is to charge requested output zeros jointly with
feedback rank while retaining the source cutoff. These are proof
refinements of the same encoder, not reasons to add mixing work.

## Sixth round: extend the dense band without an encoder change

A direct output-image bound replaces the old input-padding estimate.
For an allowed input set I, incoming state mask H, and outgoing zero
set J, let r_f=min(|J|,|H|), g=z+(1-z)/256. The scalar auxiliary kernel
satisfies

    K_I(z) <= 256^-r_f g^|I|.

The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) derives this from the
feedback kernel's output-image dimension and a systematic coordinate
projection. The bound is uniform in route and D1. It does not compare
exact kernels under deletion; such a comparison with cost 1/g is false.
Applying the direct bound to each scalar shadow gives the source factor
A0^|S|, where A0=beta*g^256.

For all actual source counts 7584--8192 jointly, the directed calculation
now proves a 35% distance bound with **more than 768 bits** of ideal-setup
margin. The log2 bound is
-773.4740494720689951607533237444136437180207059124686994558351310175682373,
with radius below 4.95e-71 at 256-bit precision. The fixed mark is
z=3/1250. The proof includes all supports, all lane ranks, at most
thirty-two missing scalar groups, and exactly 32767 state transitions.
No route-profile search or coefficient data is needed.

Source: [cheap_lane_dense_gpadding.py](theory/cheap_lane_dense_gpadding.py),
SHA256
`7bb0b9dcdbc3c575c53ba182ff0296c393047816cbbf61911e66a979a7e8a423`.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_dense_gpadding.py --precision 256
```

Two independent mathematical reviews and the independent source review
passed. Both directed precisions passed; the main agent also replayed
the 256-bit check. Each fixed replay took less than 0.1 seconds. The
old 7688--8192 certificate and all other frozen helpers are unchanged.
The new band covers 104 additional actual source counts. Nonquiet
messages with **7--7583 active groups** remain uncovered at 35%.

The next bound should retain state-generated output rank at unoccupied
positions, not merely improve the padding constant. Two compact
source-cutoff shortcuts were tested first. The half-total-rank bound
with an exact count of deficient columns gives about +5392 bits on the
six-group diagnostic. Conditioning its averaging fraction on the number
of fully ranked columns improves that only to about +5382 bits. Nearly
full columns, not just full columns, must retain their contribution.
These are floating diagnostics with possible underflow, not certificates
or counterexamples to the code. No broad middle-range sweep was launched.

A fixed middle-rank threshold also failed the bounded six-group test.
At thresholds 12 and 18, common mark z=.995 gives a log2 bound about
+5396. The zero-high-column class separately gives about +5381, so
discarding fully ranked ghost histories alone cannot repair this shortcut.
The derivation and the ordering needed for weighted rank tags are in the
rank interface. The selected class certificates still use their audited
operators; none was replaced by these failed shortcuts.

### Output-zero constraints give a useful local improvement

The new [output-rank helper](theory/cheap_lane_output_rank_tags.py) retains
positive source-rank tags after accounting for requested output zeros.
It combines output and feedback ranks through a weighted cumulative
bound, not an independence assumption. It also excludes zero subsets
that force an actually active state block to zero. This exclusion does
not apply to the dense proof's auxiliary scalar states.

To interpolate with the old output-floor estimate, split a fixed mark
z=a*b and use a^floor before expanding b^weight. Gamma=1 recovers the
old tags and gamma=0 uses the new tags. An attempted local minimum of
the two bounds was rejected and removed: they bound different relaxed
objects. No certificate depends on that attempted combination. The
global mark must remain fixed across output-rank classes unless a new
argument permits different marks.

For incoming and outgoing state weights h=j=16, the floating half-rank
moment improves as follows. These are local proof gains, not encoding
speedups or whole-code margins.

| Four block occupancies | Gain at z=.5 | Gain at z=.8 |
|---|---:|---:|
| (1,1,0,0) | 13.999 bits | 4.507 bits |
| (2,2,0,0) | 15.999 bits | 5.150 bits |
| (2,2,1,1) | 17.999 bits | 5.794 bits |
| (8,8,8,8) | 43.999 bits | 14.164 bits |

The screens use gamma=0 and the valid feasibility filter. At the old
sparse-case mark z=.995, gains were only about 0.1--0.3 bits. Testing
marks appropriate to larger supports changed the next action: this
local operator warrants a bounded intermediate-range experiment.
Some other outgoing masks still favor the old gamma=1 endpoint.

Source: 16,158 bytes, SHA256
`258791f3bd16cfbf87f9320776d6a31c3455300a89a6967bdeadce5a1a057c3c`.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_output_rank_tags.py --self-test
python -B research/workstreams/char2_field_spin/theory/cheap_lane_output_rank_tags.py --screen --incoming 16 --z .5 .8 .9 --gamma 0 1 --exclude-impossible
```

Independent source and mathematical review passed. Exact tests covered
1572 synthetic matrix/zero-set cases, 17340 rational cumulative
inequalities, and 30156 continuation comparisons. The main agent
replayed those tests and the selected local screen. These tests establish
local algebra and implementation consistency, not a global certificate.

Next implement a bounded pilot for actual source counts 256--512 using
the actual paired-route distribution and a common global mark. Retain
enough source-rank information; the tested low-dimensional shortcuts
are not adequate substitutes for the best-thirty-two-column count.
The pilot must cover every matching profile and its source multiplicity,
the final packet boundary, and every integer count in the claimed band.
Do not infer whole-band coverage from a few successful endpoint tests.
Keep exact shared-point extinction as a separate sparse-side refinement.

## Seventh round: the coarse intermediate rank cutoff is insufficient

The [actual-route pilot](theory/cheap_lane_middle_pilot.py) now combines
the new quiet and full-state output ranks with paired routing. It retains
all occupancy tails, all matching profiles, a common output mark, and
separate positive quiet/nonquiet components. Its final packet has emission
equations only. It counts full and deficient columns, but replaces the
best-thirty-two-column rank by half the total rank. This last relaxation
is deliberate and does not replace any selected certificate.

The pilot fails by a large margin. At r=256, z=.98, tau=256, and marker
scale .1, the floating log2 bounds are +258902.36 for nonquiet histories
and +186962.79 for quiet histories. All 129 individual matching profiles
are included through a conservative regional envelope. The dominant
regional terms use the pair-rich profile f=127. Its full-state nonquiet
diagonal has log2 mass 921.23, and its quiet deficient diagonal has
log2 mass 640.26. This is a diagnostic of the bound's overcounting, not
evidence of a low-weight codeword.

The positive activity split matters: removing true quiet histories alone
does not repair this relaxation. A bounded comparison also tested the
larger potential parameter 65280, several marks and coefficient scales,
and r=512. None was close to useful. Tau is a proof parameter, not the
implemented coefficient distribution. Further marker tuning of this
particular shortcut is not justified.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_middle_pilot.py --r 256 --z .98 --tau 256 --scale .1 --log2-y 0,8,16
```

The helper is 13,798 bytes, SHA256
`1b73c7736e1b8be40c860928d2d71454555ed2b283c3357fcb4f57f78be425b6`.
Independent route, local-rank, activity-flag, and terminal-boundary review
found no blocking issue. The main replay reproduced the selected values
in about 2.5 seconds, using 47,776,324 bytes of source arrays. These are
floating diagnostics; there is no new directed certificate or covered band.

Two inexpensive local refinements were checked separately. Full-state
output ranks factor into four block polynomials. At activation, requested
output zeros first delete source coordinates; the joint activation/rank
cap then applies to the residual occupancies. This preserves their
dependence and the original route markers. Exact in-memory checks passed
514 local rank cases and small paired-route substitution identities.
The formulas and their scope are in the
[rank interface](CHEAP_LANE_RANK_INTERFACE.md#factorized-output-constraints-at-the-two-extreme-state-supports).

An online column-selection shortcut was also reviewed and rejected as an
exact cutoff replacement. It can select a deficient column before learning
that enough full columns will arrive later. A four-column counterexample
is recorded in the rank interface. No implementation or certificate uses
the rejected exact-cutoff claim.

The next comparison returns to the scalar-shadow kernel while preserving
the actual source support. A sharper hyperplane-union lemma now preserves
all groups even at r=256, where the previous elementary bound vanished.
Its potential-to-witness factor is (256/255)^2. Independent review and
exact small-field enumeration checked the argument. This is a new lemma,
not a transfer of the old scalar certificate.

Initial scalar-kernel diagnostics also failed. They used representative
profiles rather than an all-profile cover. Bounds based on feedback rank
and output-image dimension remained positive by roughly 180,000 bits in
the best small tested set. Reusing markers optimized for the quiet class
was much worse: high-occupancy live packets dominated the numerator, and
conditioning back to the required input count incurred a large penalty.
The quiet-only helper, by contrast, gives very negative floating bounds
at those source sizes. Its parameters therefore cannot be imported into
the mixed-state operator without checking that operator.

The next test should use the new exact-support lemma to remove impossible
actual activation and output-zero patterns before scalar relaxation.
This must not be justified by pretending auxiliary scalar states are
nonzero. If the resulting bound still loses the regional input budget,
retain that budget by a univariate coefficient extraction, leaving only
the full-pair count as a marker. The rank interface specifies the slot
polynomial. No new class is closed by these exploratory calculations;
the nonquiet 7--7583-group gap remains unchanged.

## Eighth round: filtered scalar bounds isolate another counting loss

The [filtered scalar pilot](theory/cheap_lane_filtered_scalar.py) implements
the exact-support comparison for r<=256. It removes impossible actual
activation and output-zero patterns before introducing auxiliary scalar
vectors. It combines output-image, filtered output-moment, and feedback
bounds on the same positive kernel. The auxiliary states may still have
zero coordinates; their nonzeroness is never assumed.

The actual output floor matters independently of those scalar moments.
The pilot now splits a global mark z=a*b, charges a^floor on the actual
path, and expands only b into zero subsets. At the floor-only endpoint,
this retains output spreading and feedback constraints together. The
proof and its order of operations passed independent review. An exact
small-field test checked 448 filtered RS output cases; the main agent
replayed it with rational arithmetic.

The refined diagnostic still fails. The first screen suggested that
late activation dominated. A positive partition by the first column
containing a nonzero incoming state showed that short active excursions
surrounded by quiet packets matter too. Activation in a column's final
update is counted as incoming activity in the following column.

There was also a concrete coefficient-bound loss. Markers chosen for
long active paths poorly bound quiet paths. Recomputing the same regional
coefficients with quiet-specific markers and taking entrywise minima
fixes much of that loss. This is valid because both evaluations bound
the same coefficient. It does not make the remaining nonquiet sum useful.

The final helper also supports an activity marker v, outside the local
fractional power. If A counts nonzero incoming packet states, a cutoff B
allows separate bounds for 1<=A<=B and A>=B+1. The former uses v<1 and
pays v^(-B); the latter uses v>1 and pays v^(-B-1). The few-active class
may use z=1 and ignore output weight. These formulas passed independent
review, including the final packet.

Selected floating diagnostics at r=256 use only the matching profiles
f=0,64,128, with gamma=1 and independently refined quiet markers:

| z | alpha | v | Raw marked nonquiet log2 sum | Use |
|---|---|---|---|---|
| .85 | .523 | 1 | +129622.685 | Ordinary nonquiet diagnostic |
| 1 | .523 | .01 | +27540.295 | Add 6.643856*B for the few-active class |
| .8 | .5 | 4 | +157065.593 | Subtract 2*(B+1) for the many-active class |

None of these activity witnesses closes its intended class. In the first
row the quiet component is -8707.347, but the nonquiet component remains
positive. Quiet and active contributions are tracked positively, not by
subtracting nearly equal totals. These are failed selected-profile
diagnostics, not an all-profile bound or evidence of a low-distance word.

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_filtered_scalar.py --r 256 --z .85 --alpha .523 --gamma 1 --profiles 0,64,128 --optimize-markers --refine-quiet-markers
```

For the second and third rows, additionally set `--activity-mark .01`
or `--activity-mark 4` and use the displayed z and alpha. The main agent
replayed all three final witnesses. The source is 20,454 bytes, SHA256
`4198a26f54c3074fc953bb3f752c236de1fb0ab0e2d9c9c4b74bb10bd9405d78`.
Independent source and mathematical audits found no blocking issue.

At the eighth-round checkpoint, the next diagnostic was to compute
the regional input-count coefficient exactly for the pair-only profile
f=r/2. The [rank interface](CHEAP_LANE_RANK_INTERFACE.md#next-diagnostic-exact-coefficients-for-pair-only-profiles)
specifies a degree-f positive matrix polynomial and its terminal boundary.
This bounded comparison will distinguish route coefficient loss from
remaining local-kernel loss. It must not use FFT cancellation or discard
small positive terms. Extend profile coverage only if that test helps.
The nonquiet 7--7583-group gap and all selected certificates are unchanged.

## Ninth round: exact routing coefficients and local support filters

The pair-only diagnostic is now implemented. It computes the exact
regional average of the current local upper bounds by positive logarithmic
coefficient convolution. It neither uses FFT subtraction nor drops small
quiet contributions. Independent review and direct exact slot enumeration
passed: 917 configurations and 186 logarithmic coefficient comparisons.

This test has a restricted scope. A nonempty proper support cannot be
fully paired under both alternating matchings. At r=256, the repeated
f=128 marginal is therefore an artificial diagnostic, not a source class
or an established worst case. Its regional coefficients are valid for
the single prescribed marginal. The nearest feasible joint profiles
(128,127) and (127,128) have not been evaluated by this helper.

At r=256, z=.85, alpha=.523, exact coefficients lower the artificial
nonquiet log2 expression from +129622.685 to +123049.174. Allowing a
local choice of gamma in {0,.5,1} lowers it to +118351.16. These gains
are insufficient: route coefficient slack is not the principal remaining
loss at these test points.

We then strengthened the joint output/feedback bound. A path with nonzero
outgoing state requires a nonzero vector in the output-constrained
pre-refresh image. Averaging this necessary condition over D1 gives
the audited factor min(1, [256/255]255^(-(16-j))) times the output-only
kernel for j>0. It applies to live-to-live transitions, not just activation.
It does not apply to extinction or the terminal. An exact GF4 test checks
252 rank configurations and 28 averages, including zero-dimensional images.

The final selected diagnostic, with both refinements, is

| Quantity | log2 value |
|---|---:|
| Regional quiet-to-quiet contribution | -957.749203 |
| Regional nonquiet-to-quiet contribution | -496.644906 |
| Artificial whole-chain nonquiet expression | +105541.908209 |
| Positive restart contribution within that expression | +101629.967278 |

The restart contribution selects a nonquiet zero-to-zero path in every
region. It survives the valid exclusion requiring at least 33 nonquiet
outer columns. Thus that exclusion alone cannot fix this particular bound.
These are diagnostics of an upper expression, not failure probabilities
or evidence of a low-distance codeword.

Reproduction of the selected diagnostic and exact tests:

```text
python -B research/workstreams/char2_field_spin/theory/cheap_lane_paired_coefficient.py --self-test
python -B research/workstreams/char2_field_spin/theory/cheap_lane_paired_coefficient.py --r 256 --z .85 --alpha .523 --gamma-mix --joint-live-filter
```

The source is 20,058 bytes, SHA256
`4397d6b39f56786243fd68347ab01cacbd19970068480061b804951976a39969`.
The main agent replayed both commands after independent source review.
The selected run took 6.7 seconds; generated output files: zero.

The next structural refinement is now proved locally. Required nonzero
state coordinates define at most h+j proper hyperplanes in the local
E-linear constraint kernel. Actual F-valued feasibility implies that
none of these forms vanishes identically. A hyperplane union bound then
compares the existing scalar mass to its exact-support restriction at
factor 256/(256-h-j), at most 8/7 before alpha. The
[local support comparison](CHEAP_LANE_RANK_INTERFACE.md#enforcing-auxiliary-support-with-an-explicit-local-loss)
also gives a five-form variant at factor at most 256/251, retaining only
one nonzero coordinate per active emission block and one outgoing coordinate.
Both variants passed independent mathematical review. They preserve the
implemented coefficient distribution and require no encoding change.
The standalone support-filter checker passed 883 exact small-field
comparisons and named extinction/terminal cases. The main agent replayed
it in under 0.01 seconds. Its source is 11,029 bytes, SHA256
`3d9bce383376146d6506294f471d4e13acb2cbf50d1a39cd1ca1cdf0bd6a8940`.
Run it with `python -B` and `--self-test`. The lemma does not permit
conditioning on exact support without paying its explicit local factor.

The inexpensive block-support version was tested in a final bounded
wrapper. It combines a nonzero RS block's outside-input output floor
with the uniform auxiliary input moment and D1 feedback bound. The
incoming-coordinate normalization cancels most of the hyperplane loss.
Nevertheless, the current joint-live/gamma caps already dominate it:
the largest local improvement is 0.0000435 bits, and the regional and
chain diagnostics above are unchanged to displayed precision. Independent
source review passed; the main replay took 8.7 seconds. Source: 8,005 bytes,
SHA256 `6c826ae6ff8e5cbfb256313b85f8155ab34b6f348a97604c4b0309a02e540b71`.
Reproduce with `python -B` on `theory/cheap_lane_block_support.py --run`.
Its 175 block-placement checks and inherited exact tests all pass.

Extinction has a more concrete geometric refinement. With one allowed
input per emission block, active extinction forces u=C_a*x with x!=0.
The actual incoming mask is fixed by point a; incompatible output-zero
requests can be discarded before relaxation. Exact AES-GF256 arithmetic
gives output weight 15 at points 0,2,3 and 16 at all other points, versus
the previous floor 12. An independent algebraic check confirmed the
table using the implemented normalized basis. A preliminary unscaled
GF16 calculation was rejected and is not used. Full pairs must average
their shared point jointly. This refinement is recorded but not yet
integrated into the regional calculation.

The next step from the ninth round was to integrate this exact singleton extinction/output cap and test its
effect before extending to two-input blocks. Do not launch a broad parameter
scan or extend profile axes before this local refinement helps. No new
distance class is certified in this round; the nonquiet 7--7583-group gap
is unchanged.

## Tenth round: distinguish extinction from sustained live paths

The singleton extinction refinement is now integrated. Full pairs average
their shared point before mask aggregation; a shared point does not require
both block amplitudes to be nonzero. The helper preserves quiet entries,
live destinations, and the output-only terminal. Independent source review
and the main replay passed.

At the same artificial repeated-marginal diagnostic r=256, z=.85,
alpha=.523, the selected results are:

| Quantity | Previous | Singleton refinement |
|---|---:|---:|
| Regional nonquiet zero-to-zero log2 term | -496.644906 | -498.522095 |
| Regional nonquiet 16-to-16 log2 term | -499.933122 | -501.809097 |
| Artificial chain log2 expression | +105541.908209 | +105301.664615 |
| Selected restart contribution | +101629.967278 | +101149.407034 |

The gain is real but too small. These numbers concern an upper expression
with incompatible repeated route marginals, not a feasible source class,
failure probability, or codeword. No new distance class is certified.

An oracle diagnostic deletes every regular live-to-zero term. This is
**not a probability upper bound**. Even this deletion leaves the artificial
expression at +25163.442971. Thus extinction-only improvements cannot close
this fixed operator at this fixed parameter point. This conclusion does
not cover other fractional-moment parameters or other proof operators.
The two-input extinction formula is preserved in the rank note, but its
implementation is not the next priority.

We also proved and tested a nonzero RS output moment for sustained active
states. It uses the full positive requested-zero expansion, rather than
only the block's minimum weight. A local support comparison removes the
auxiliary zero block at explicit factor q/(q-b-1), or q/(q-b) at the terminal.
The implementation keeps the old joint bound separate from these modified
terminal moments. Two mathematical reviews and an independent source
review passed.

The new moment gives no meaningful improvement at the selected point.
The chain expression remains +105301.664615; the largest local gain is
0.0000435 bits. The existing gamma-zero cap already captures almost all
of this output-only information. This is a tested dead end for that local
refinement, not evidence against the construction.

Reproduce the fixed diagnostics with `python -B`:

```text
research/workstreams/char2_field_spin/theory/cheap_lane_singleton_extinction.py --run
research/workstreams/char2_field_spin/theory/cheap_lane_nonzero_moment.py --run
```

The singleton helper is 11,518 bytes, SHA256
`992e6b7bd02e4785b9bfc4ae61437625afabb28a430e085fe32d60a3ab4a3453`.
The nonzero-moment helper is 11,913 bytes, SHA256
`075aa24ee05467fd04352d9c1cbe05cf0bd75e9d439b95850006c9db576d76b6`.
The latter includes 3,840 exact GF8 mask/puncture comparisons over 6,545
nonzero words, including 256 full-dimensional equalities. It also replays
the singleton and positive-convolution tests. Main replay took 9.0 seconds.

The nonquiet 7--7583-group gap at 35% remains open. The next proof work
must preserve information across packets or strengthen the joint
output/feedback calculation; further output-only floor variants are not
justified by these tests. The encoder, historical measurements, and
certified 3.8-ms baseline remain unchanged.

One endpoint check also replaced alpha=.523 by alpha=1, with r and z
unchanged. It worsened the artificial expression to +259466.690614:
the source prefactor rose to +435926.662268, outweighing stronger local
terms. This was one bounded run, not an optimization over alpha. It
does not exclude other fractional-moment parameters.

The next structural lemma is now independently reviewed: retain the
possible state subspaces on both sides of a live packet boundary. A
full-column-rank leading-theta matrix rules out some low-weight paths
for every byte a and every F-valued message. Unlike the previous
excursion filter, it applies without a return to zero. The
[two-live-packet interface](CHEAP_LANE_RANK_INTERFACE.md#two-live-packets-keep-the-state-subspaces-at-the-cut)
states the exact spaces and conditioning. Next, construct one positive
two-packet bound before any profile sweep. A dimension-only probability
estimate requires spaces independent of the current D1; diagonal
randomization alone cannot mix arbitrary coordinate subspaces.

## Eleventh round: group masks across packet boundaries

The previous round was progress: it identified which local refinements
were dominated and supplied a new two-packet feasibility lemma. The
next inspection found a limit to that lemma at z=.85. In a full-state,
output-only pair, its removable scalar mass is at most 2^-97.9884 of
the empty-zero contribution. This is not a useful leading refinement
for that channel. It does not rule out more selective joint kernels.

Instead, the new operator pools intermediate state masks before the
fractional power. Independent derivations and source review passed.
Only pair-boundary masks are split by subadditivity. Raw alpha=1
kernels supply the internal mask sum; endpoint and route counts stay
outside alpha. The code and setup distribution do not change.

The fixed r=256, z=.85, alpha=.523 diagnostic gives:

| Quantity | Same-raw per-packet bound | Grouped pairs |
|---|---:|---:|
| Regional quiet zero-to-zero log2 term | -957.749203 | -957.749203 |
| Regional nonquiet zero-to-zero log2 term | -498.484523 | -521.106226 |
| Regional nonquiet 16-to-16 log2 term | -501.771660 | -524.614229 |
| Artificial chain log2 expression | +105306.469756 | +99521.276831 |
| Selected restart contribution | +101159.025435 | +95367.869346 |

Pooling saves 5785.192926 bits in the artificial chain expression, but
does not close it. The comparison uses identical raw local caps on both
sides. It differs slightly from the previous +105301.664615 checkpoint
because the singleton point average is taken at alpha=1 before powering;
the dominated extra output-only caps are omitted. These remain artificial
repeated-marginal diagnostics, not feasible source-class bounds.

The helper `theory/cheap_lane_grouped_masks.py` is 12,391 bytes, SHA256
`4526522c47c9920482f8990bb104243eada12318774ee2d0dbb68fda3704ab36`.
Run it with `python -B` and `--run`. Exact tests expand actual mask
multiplicities, check positive activity classes and the output-only
terminal, and recover ordinary polynomial convolution at alpha=1.
Independent source review passed. The main replay took 17.2 seconds;
its polynomial arrays occupy 254,016 bytes and no files were generated.

A reviewed route-Jensen normalization allows larger disjoint tiles:
divide each raw degree-d coefficient by its binomial route count and
endpoint-mask count before powering, then restore both counts. This is
less sharp than explicit occupancy averaging, but avoids exponentially
many occupancy axes. The bounded comparison gives:

| Tile | Artificial chain log2 expression | Selected restart contribution |
|---|---:|---:|
| Two packets, route Jensen | +100097.612672 | +96517.697437 |
| Four packets, route Jensen | +98082.918905 | +95147.277285 |
| Sixteen packets, route Jensen | +103490.776101 | +103442.866144 |

Four packets improve the exact-pair score by 1438.357926 bits. Sixteen
packets regress by 5407.857196 bits against four: route averaging loses
more than the additional mask pooling gains. Stop this tile-size scan.
The eight-packet case was not run. The quiet terms also become looser;
the separate quiet-class certificate is unchanged and remains available.

The wrapper `theory/cheap_lane_mask_tiles.py` is 9,900 bytes, SHA256
`a96bcfc95c927501264f9787692e509b03f64a8f968570bc0ac631f581d8f34e`.
Independent proof and source review passed. Tests cover 260 alpha=1
comparisons, exact binomial route counts, positive activity classes,
the terminal, and 80 explicit-occupancy/Jensen inequalities. Reproduce
one selected case with `python -B` and `--tile 4`. The main independent
replay reproduced +98082.918905 in about eight seconds. Its tile arrays
occupy 334,368 bytes; estimated peak scratch is below 64 MiB.

The certified support coverage remains unchanged. The next inspection
retains the outer code's input-weight floor inside the scalar comparison.
The current density relaxation admits auxiliary inputs that need not
satisfy that floor; local mask tuning does not address this loss.

### An outer-weight marker needs sharper live kernels

Two independent reviews establish a weight floor for exact-T scalar
shadows: at least 33r nonzero input records. Insert the positive marker
v^(wt(X)-33r), for v>=1, before density domination. It contributes the
global factor v^(-33r alpha). The rank-interface note gives the argument
and exact quiet kernels, including a positive MDS weight polynomial.

A bounded first implementation multiplies each ordinary raw packet cap
by v^ell, where ell counts its allowed input coordinates. Quiet-to-quiet
and quiet-terminal caps instead use their exact marked moments. Both
the v=1 baseline and v=1.01 case use those same exact replacements.

| Input mark v | Raw artificial log2 expression | Threshold correction | Adjusted expression |
|---|---:|---:|---:|
| 1 | +98082.918905 | 0 | +98082.918905 |
| 1.01 | +98568.584139 | -63.426048 | +98505.158091 |

The marked ceiling worsens by 422.239186 bits: the live-packet inflation
exceeds the threshold rebate. Stop this coarse marker test without a
sweep. It does not refute using the true joint input/output moment.
The main independent replay reproduced both expressions in 17.6 seconds.

Source: `theory/cheap_lane_input_marker.py`, 6,612 bytes, SHA256
`9fd71bffd1beac7cb612a4456387af4c149cc30cf1326f4a79d92a3befceb71b`.
Run with `python -B` and no arguments. Independent source review passed.
Tests verify 153 nonnegative integer polynomial coefficients, 17 exact
kernel cardinalities, and the inherited mask/route identities. At v=1,
the exact quiet replacement agrees with the selected old cap to rounding.

This round improved the proof operator, not its support coverage. The
next step is a genuinely weighted live-packet kernel or a comparison
retaining correlations imposed by the outer constituent. Keep the new
pooled-mask operator as a comparison, but do not keep tuning its tile
size or coarse v^ell marker. The remaining nonquiet 7--7583-group classes
are still open at 35%; no whole-code claim transfers to the lane mixer.

## Twelfth round: isolate quiet excursions before source counting

The [excursion analysis](CHEAP_LANE_EXCURSIONS.md) identifies a concrete
mechanism behind the loose scalar majorant. For a selected regional
16-to-16 entry, restricting every intermediate state to full support gives
-1000.287280; allowing any nonempty support gives -943.647442; admitting
zero states gives -525.633675. These restricted operators are diagnostics,
not complementary probabilities. The full-state local moment is already
tight, so further tuning of that cap is not the current priority.

An audited surplus-pivot lemma exploits the full 65280-point coefficient
set when charged equations outnumber possible independent pivots. Its
coarse feedback-rank marker does not help this path family: t=.99 worsens
the artificial score from +98082.918905 to +98299.703806. The ledger stops
that scan. Enlarging the projective potential and then collapsing it to
one scalar moment also loses the intended gain; two outer shadows cannot
be treated as independent under the same Toeplitz map.

There is a real bad-route example. Route all inputs from 256 consecutive
groups into the first 44 packets of each region. The 32768-dimensional
source space loses at most 4096 dimensions by imposing zero state after
these windows. Thus some word has weight at most 34.375% of N, for every
other setup realization. Its actual source support can be 224--256 groups,
not necessarily all 256. The fixed route event has floating log2
probability -52149.139093. By contrast, an analogous artificial fixed-path
term retains a huge source-count prefactor and gives +73742.613560.
This is an overcounting mechanism, not evidence of frequent bad codes.

Two new exact bounds now separate useful event classes:

| Event class | Exact probability upper bound |
|---|---:|
| Some set of at least seven groups has every allowed position within a separately chosen 44-packet window in every region | 2^-592 |
| Some word has exact support equal to a consecutive 256-group set and all endpoint-live regular steps within such windows | 2^-579 |

The first uses a seven-subset reduction and a direct routing bound. The
second allows coordinates outside the windows: quiet steps force C X=0,
so their actual outer-code constraints are retained. Conditional on the
route, the necessary-condition probability is bounded by
min(1,eta_r beta^r 256^-D_Q), where D_Q is the sum of quiet feedback ranks.
A separate fractional exponent can then be chosen before route averaging.

The consecutive-support computation includes all 8192 source intervals,
all 85^256 window choices and both feasible matching profiles (128,0) and
(127,2), while omitting the unused terminal update. It uses positive integer
coefficients at y=3/4 and exact rational comparisons, not floating
logarithms. Independent proof/source review and the main replay passed.
The resulting 579 bits concern this restricted state/support event, not
the code's global failure margin. No complete active-group-size band is
newly closed: nonquiet 7--7583-group distance classes remain open.

The next step is to reuse the same quiet-constraint table across
consecutive source sizes, then admit several source runs with the actual
joint matching counts. After that, combine this event bound with the live
output analysis. Do not resume coarse tile or full-state-cap tuning, and
do not cap individual transfer entries at one.

The new helpers are reproduced with `python -B` and no arguments:

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_path_decomposition.py` | 6013 | `c8207a99d9a3482ea37f8ab5bde19fec908fe0387b2f8b336655b291f37d6f25` |
| `cheap_lane_surplus_pivots.py` | 5419 | `53c9877a3763b379e8f00b7de1eeeb7ac69f84af766bd2d9c0edd28c5afb293a` |
| `cheap_lane_window_diagnostic.py` | 8084 | `a76233f72dea9e0578ea53521f402df68a6de92ede2202ed76a225ef1ba3cd5f` |
| `cheap_lane_window_cover.py` | 6797 | `b4c7df086502b283cb875e84a0b258df104fee4a36ede4ff786eb18df8d4f47a` |
| `cheap_lane_quiet_window_cover.py` | 8811 | `e6c058696a5383b97df74c65a1f5a18866854fe010c202da075f950fa6c2f69d` |

The last helper's unreduced integer numerator and denominator have hashes
`e2fa7ad9fdc5e392ba8565d3314b4df61f8c3eb99fb96ab92e443a81371e7b34`
and `f8aa0328b1faaafbab18ddeddbb92e71cc3d67e6791b671eaa4286df1d9e4f18`,
respectively. Hash the minimal big-endian byte representation of each
positive integer. The integers are recomputed in memory, not saved.

## Thirteenth round: arbitrary supports in the quiet-window class

The [excursion analysis](CHEAP_LANE_EXCURSIONS.md) now extends the exact
consecutive-support calculation across every size 7--256. Sizes 196--256
jointly have a 45-bit bound for the same 44-packet-window trajectory
class. Including size 195 reduces this dyadic bound to 36 bits. The
helper checks every size and all even/odd matching profiles; it does not
infer monotonicity. Its 270 normalization and 278 cyclic-profile checks
passed, as did independent source review and the main exact replay.

A sharper quiet-only scalarization removes the live-potential comparison
factor. For exact support T with |T|<=256, existence of an F-valued
quiet-kernel witness implies at least 255 E-valued exact-T witnesses.
Thus the conditional event bound is min(1,beta^|T|*256^-D_Q/255).
For larger source sets, a scalar witness misses at most ceil(|T|/256)-1
groups. Union over these smaller supports directly; do not count their
possible supersets. This does not transfer output weight or live state
behavior from the F witness to the scalar witness.

The arbitrary-support screen includes all 16640 size-256 ordered
matching profiles. Three complete rank-mark bounds improve its log2
expression from +511.841 to +119.297. A new singleton-deletion count
reduces this only to +118.728: correcting the former worst profile moves
the bottleneck elsewhere. The count passed 2387 exact small-cycle
comparisons and four large boundary checks.

The useful next refinement extracts independent categorical draws from
the quiet-half polynomial. The conditional mixture is exact and retains
the shared coin for each half; treating all its slots as independent
Bernoulli selections would be incorrect. A compound-tail and multinomial
maximum-atom bound lowers the complete floating expression to
-76.03854856. All profiles are negative; the largest term is about
-85.0600 at (f0,f1)=(115,55). The construction, encoder, and setup
distribution are unchanged. Two independent proof reviews and a source
audit passed, with 459 coefficient, nine mixture, 36 compound-tail, and
32 mode checks. Directed 192- and 256-bit replays now prove the complete
restricted-class bound below 2^-76. The verifier checks 774 multinomial
modes exactly and evaluates all 387 regional witnesses and 16640 source
profiles with directed arithmetic. The main 256-bit replay took 3.85
seconds and reproduced the same rational-witness hash. Its log2 interval
has radius at most 5.22e-70. This is not an all-trajectory size-256 bound
or a whole-code certificate.

The next global test is now explicit: partition all trajectories by the
number of state activations and quiet steps. Combine a two-state
quiet-equation bound with the marked activation/return and persistent-live
output bounds. Preserve matching order when multiplying regional matrices;
they need not commute. Test a small predeclared threshold grid at r=256.
The success criterion is an all-trajectory 40-bit class bound, not another
restricted-window gain. If the best split still fails by thousands of
bits, stop coefficient polishing and identify which branch needs a new
structural argument. The nonquiet 7--7583-group whole-code gap remains.

The new helpers are reproduced with `python -B`:

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_quiet_window_band.py` | 7378 | `6cac5b9f735d8169f31bfc98ec82da485a8ed9f5be674d3fd990aa961cb40a5a` |
| `cheap_lane_quiet_profile_screen.py` | 11244 | `28083dbe863b0408c2f5b7c22cbc09302f444bdb4682d1d2bb9ce1af7072e151` |
| `cheap_lane_quiet_profile_marks.py` | 6676 | `5d1402f5869d7cd3ffcfd97e8414fd58f983a6721b45ec723aa934b869bbeb37` |
| `cheap_lane_quiet_profile_counts.py` | 8463 | `ec52c2706499e7722482283fc194534362f70ebc058a749e9f87c21e433178eb` |
| `cheap_lane_quiet_profile_smoothing.py` | 10480 | `473428234351fd0813a4d5cf3e9bcf4db96179ef7e1cd0b884c232cd0e8a2a73` |
| `cheap_lane_quiet_profile_verify.py` | 10219 | `9f8cc18e38634190c9cb5da20a7d1252cc13fadea751d698ab5521ea7895b791` |

The directed verifier accepts `--precision 192` or `--precision 256`.
Both replays select rational witnesses with SHA256
`7894707571b15d0a5a7e5c3e0a7f9c20d797291044315a0e46e88abf900de69f`.

The exact interval-tail numerator has SHA256
`0b29bfc30d054627d94c0a652bd7eee573e17d8bb37287059f4a63d081106228`;
its denominator is 2^579. The integer has 534 bits. The complete per-size
margin table has hash
`7557e6df49c8bf5ca458c474f2f64c8f64d78d1cbe3c73fc04cc6ac495ba5357`.
These values are recomputed in memory, not saved as a dataset.

## Fourteenth round: testing an all-trajectory split

The new quiet pilot sums all binary trajectories, rather than imposing
regional windows. It retains the true alternating matching order and
32767 regular updates. A uniform independent-slot extraction bounds the
coefficient mass for each path and hence the positive path sum. It does
not assume the guaranteed free slots from the earlier window calculation.
Independent source reviews and exact toy checks passed; the main replay
took 2.84 seconds and reproduced all nine threshold cells.

For the four predeclared size-256 profiles, the floating union is below
2^-3141 at K=256,Q=24576 and below 2^-510 at K=512,Q=24576. Here K bounds
the number of state activations and Q is a lower bound on quiet regular
steps; these K values are not message lengths. At K=256,Q=20480 the
expression remains positive, with log2 value +3085.348. These selected
floating results are not certificates. Do not expand to all profiles
until the complementary branches can use compatible thresholds.

A second helper checks fixed positive contributions to the complementary
scalar majorant. It uses exact oriented route coefficients for actual
consecutive-support profiles, not repeated artificial marginals. Two
explicit paths lie outside the already certified window class: separated
excursions and excursions straddling regional boundaries. At the original
three parameter anchors, they reject the selected majorant by thousands
of bits, regardless of admissible count marks. But both contributions
become negative at alpha=.05,z=.9 or .99. Thus those first failures were
parameter-dependent, not a structural impossibility or a code attack.

An analytical all-live contribution still rules out the two high-z
anchors. One orientation contributes at least
128*(256^32768*z^1363149)^alpha. At z>=.9,alpha=1/20 this exceeds
2^1754.625. The derivation uses the actual frozen caps and the event that
the first paired slot is full. It does not lower-bound actual failure.
The all-live path is actually impossible at the target output weight;
the obstruction applies only to the untruncated majorant. A cheap global
floor retains this missing restriction: if S sums ceil(h/4) over actual
incoming masks, including the terminal packet, then output weight is at
least 13S-256r. At r=256, a bad path must have S<=61503. A raw mark w^S
and correction w^(-alpha*61503) impose this constraint. The next test
sums intermediate supports and trajectories while retaining this mark
and the branch marks before raw pooling and fractional powers. More
quiet-only gains are not enough.

The fixed-path main replay took 6.27 seconds. It passed 36 coefficient
geometry checks, 256 direct slot assignments, 37 rational convolution
checks, six terminal checks, 18 marker inequalities, and five concrete
out-of-window path checks.

The subsequent 17-state pilot includes the global floor, all state paths,
the actual alternating profiles, and raw 16-packet pooling. It averages
route placements before the fractional power. None of its eight selected
anchor/branch witnesses closes. The best log2 values are +37827.018 for
quiet steps<=24575 and +39277.686 for births>=513, both at z=.85,alpha=.05.
The bounded optimizer reaches its evaluation or iteration limits; these
are selected witness values, not lower bounds on the best possible bound.
The main replay reproduced both in 13.66 seconds. Independent source
review and 96 exact route normalizers, 2500 output-floor comparisons,
12 route/mask entries, and 12 Jensen inequalities passed.

This identifies a specific next A/B test, not a construction failure.
Keep each tile's full-pair and singleton degrees until after the fractional
power. Evaluated-route Jensen can overweight rare low-constraint routes;
degree-wise Jensen preserves their multiplicity. Compare one fixed
anchor and branch against the current operator before another parameter
scan. If that gain is insufficient, a separate compact refinement tracks
(h,b), where b is the number of nonzero RS state blocks. There are only
35 feasible classes. For h=8, only six masks occupy two blocks, yet the
current 17-state bound applies that worst floor to all C(16,8)=12870 masks.
The exact class multiplicity is C(4,b)*[x^h]((1+x)^4-1)^b. This is a
candidate refinement, not an implemented or certified gain.

No new source-support class is certified in this round. The whole-code
nonquiet 7--7583-group gap remains; the active goal is not complete.

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_quiet_path_pilot.py` | 12666 | `cd73966925964f1da5a93deb69d84aa8cf578446b12a5869ad6a8e871f2eb750` |
| `cheap_lane_path_marks.py` | 17117 | `e5bd5beda090ca63e6648a07106181d733327a4d5614a9dfb8a340a98f07c483` |
| `cheap_lane_marked_state_pilot.py` | 15796 | `7ca2f572886d3b19c6a37ed2513cde9c29895e9c8424219f8f547317db0d09b8` |

Run the two pilots with no arguments and the fixed-path helper with
`--run`, using `python -B`. All are bounded, stdout-only diagnostics.

## Fifteenth round: distinguish counting losses at fixed parameters

The degree-conditioned comparison keeps the actual consecutive-support
profiles and the previous z,alpha,quiet,floor witnesses fixed. It retains
full-pair and singleton degrees through the raw 16-packet product, then
applies the fractional power separately to each coefficient. Terminal
handling and source accounting are unchanged.

The main replay reproduced +37827.01773005103 for evaluated-route Jensen
and +36376.29710010885 for degree-conditioned Jensen at identical route
markers. A short route-marker optimization gives +36327.86231064020;
evaluated-route Jensen at those same new markers gives +37787.75257500522.
The improvement is about 1451 bits. The earlier speculative much larger
route-averaging loss is not supported by this controlled test.

The comparison also safely deletes impossible tile degrees. It therefore
measures that deletion together with degree conditioning, not Jensen
alone. Independent review and the main 29.24-second replay passed the
exact toy coefficient, normalization, Jensen, and terminal checks.
No additional support class is certified.

The next two tests keep this anchor fixed. One pools longer raw segments
to reduce boundary-mask counting. The other distinguishes the number of
active RS state blocks, rather than applying the worst block placement
to every mask of the same weight. Restricted zero/full-state channels
must be reported separately. Neither test changes the implementation.

Both tests have now completed. Larger pooling gives full17-state log2
expressions +19350.27 at 64 packets and +16552.60 at 128 packets. The
all-packet zero/full16 control remains positive at +13187.07947 for the
128-packet pool. The changes also alter degree truncation and route
Jensen, so the roughly 19824-bit improvement is not pure mask-count loss.
The run passed independent source review, took 56.15 seconds, and peaked
at 98672640 bytes. Restricted channels are not full-code bounds.

The separate 35-state evaluated-route test tracks (h,b), where b counts
active RS state blocks. Its compatible-placement caps and exact mask
multiplicities passed independent mathematical and source review. Raw
matrix products factor through dimension 17; no 35-dimensional powers
are needed. The main replay improves +37827.01773005103 to
+37145.25255760601, only 681.77 bits. The zero/full16 control is unchanged
at +17650.47667666854. All 35 class counts, 3125 placement/rank checks,
factorization, terminal, and old-score controls passed in 3.22 seconds.

An exact regional-coefficient control then removes Cauchy evaluation
without changing the raw caps or selected parameters. The all-packet
zero/full16 expression falls from +13187.079472798836 to
+11617.188069039783. Independent review and the 2.78-second main replay
passed. This remaining positive control is unaffected by intermediate
mask subdivisions; it does not establish low distance for the code.

The omitted branch restrictions now provide a more specific next test.
For zero/full16 paths, the intended branch has at most 512 births and
between 7680 and 15375 active incoming packets. The current control drops
the birth cap and selects the neutral output-floor mark. Compute its
tilted counts by positive dual convolution before deciding whether to
restore these restrictions or change the Chernoff anchor. Do not infer
that the full parameter family fails from the current witness.

That derivative test is now complete and independently audited. Positive
dual convolution gives tilted means of 148.09 births, 25579.98 quiet
steps, and 7039.92 active incoming packets. The exact path-count identity
holds, including the terminal indicator. The main 2.61-second replay
passes 432 rational dual entries and four finite-chain derivatives.
The corrected slopes for birth, quiet, and floor marks are -26.25029,
+72.49375, and -2405.21101 bits per natural-log mark. Convexity rules out
improvement by moving either the birth or floor mark alone below one
at this fixed witness. The quiet penalty is the parameter currently
miscalibrated for the exact regional control.

Next: solve the single convex quiet-mark problem until E[q0]=24575,
then recheck the other slopes. Keep z=.85,alpha=.05 and the local caps
fixed for that isolation. If a large positive control remains, revisit
the zero/full-state argument rather than add more mask subdivisions.
This round made diagnostic progress but certified no additional support
class. The nonquiet 7--7583-group gap remains; the goal stays active.

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_degree_conditioned_tile.py` | 13806 | `264a216fd5beb0c2b2d18826e5fc3f7cd4638561ff2aa37d9862382e36bbd0d2` |
| `cheap_lane_pooling_scale_pilot.py` | 9097 | `2c306223f5b1f876670aede55ca80b905f16fb98364a840a645fea3b19cc32ad` |
| `cheap_lane_active_blocks.py` | 17086 | `b703bd658a852a62212e61947efff43be6bf005e08560c28c7147d84ee05e46c` |
| `cheap_lane_exact_region_control.py` | 5909 | `e2dbb5fdd206feb0453dee426c17f3eb95283433542e411733c65238f0175586` |
| `cheap_lane_branch_sensitivity.py` | 12111 | `84b6b757a3fc4e85f28c0c3e9f59f6063ae8daf70cc3704aca6feb9ebe27600b` |

## Sixteenth round: preserve quiet routing before pooling live paths

The quiet-mark calibration finishes the previous diagnostic. Its best
floating score is +11592.38616256, improving only 24.80 bits. At its
stationary point the birth and floor slopes remain negative. Joint
convexity rules out these three marks as the repair at this fixed anchor.

A controlled fixed-path comparison identifies a much larger loss.
The same budget-admissible boundary-straddling zero/full16 path gives
-4598.42402459 before regional route Jensen and +10679.55735626 after
it. There is no intermediate-mask sum in this comparison. The difference,
15277.98138 bits, is purely the order of route averaging and the fractional
power. Its counts satisfy all three branch budgets, so imposing those
budgets cannot remove this term from that selected pooled majorant.
This is a failure of the tested bound, not evidence of low code distance.

The replacement control retains route occupancies before alpha and sums
all zero/full16 paths by a two-state matrix polynomial. Exact regional
coefficients give +537.70899537 at neutral marks. The pure-quiet term is
-4480.15361210. Its tilted count means violate the birth and floor budgets,
so those constraints are now potentially useful; the old pool's quiet
mark is not. A single comparison at the earlier z=.9 anchor regresses
to +3546.74229505. No wider parameter search was performed.

Most importantly, weighted route Jensen retains the recovered gain while
allowing a form of pooling. Split a fixed path's raw product into B*V,
where B contains quiet/birth/return factors and V contains sustained-live
factors and the terminal. The valid bound

    E[(B V)^alpha] <= E[B^alpha]^(1-alpha) E[B^alpha V]^alpha

matches the route-first fixed-path result to 3e-11 bits. It recovers the
entire diagnosed loss in this test. The proof, source, exact-rational
inequality tests, route normalization, and terminal handling passed
independent review. This is floating diagnostic evidence, not a directed
distance certificate.

Next: fix the quiet/live schedule and birth/return endpoint masks, then
pool partial nonzero masks only within sustained-live segments. Keeping
those endpoints fixed makes B common to the pooled terms. If endpoints
are later summed, their factors must be retained as indices or moved
inside V; only quiet factors are common to an arbitrary binary schedule.
This tests whether the corrected route order also handles mask entropy.
It changes neither the encoder nor the actual limited refresh distribution.

No new source-support class closes. The nonquiet 7--7583-group gap remains,
and the goal stays active. The certified 3.8-ms baseline, candidate encoder,
and previous proof evidence are unchanged. The detailed comparison and
the next compositional interface are in `CHEAP_LANE_EXCURSIONS.md`.

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_quiet_stationarity.py` | 7259 | `57ab6c41e0b7e43b0f08db9033c5f876290ceb7ec09efe7060ffe25750d843b1` |
| `cheap_lane_fixed_path_jensen.py` | 5598 | `1214da897e45e40bd5196cff467454bb5ced0535a0c064febbd014c76415e12c` |
| `cheap_lane_route_first_control.py` | 6936 | `fbf6d6b6649e0962feb885db0cdaff5348f7300b11e8ec03116016da7c8053b8` |
| `cheap_lane_route_first_sensitivity.py` | 5373 | `bd8f8192ffbe674bd93192bbcb9d6435a8e3fc68eec55d276dc4f4266f0dc52b` |
| `cheap_lane_quiet_tilt_jensen.py` | 7858 | `26b019f4d967603f3bd7e8aa14229fe15d2fec4214993e507302f3c377c319fa` |
| `cheap_lane_route_first_anchor.py` | 4127 | `56c18ee6f4107888d5ea505064c19e7bda25512c422782b2ac2f869214074744` |

## Seventeenth round: an all-trajectory witness for consecutive supports

The fixed-schedule test confirms that sustained-live mask pooling is
cheap when quiet routing is retained: the selected log2 expression changes
from -4598.42402459 to -4544.22732062. This costs only 54.20 bits, with
birth/return endpoints and regional cut masks still fixed. It is a
diagnostic for one schedule, not a source-class bound.

A new fixed positive cone then sums every actual nonzero state mask and
every binary zero/live schedule before the fractional power. The two-state
operator retains exact regional route coefficients. Its first score is
+4312.54152383 for the actual consecutive 256-group source class at
z=.85,alpha=.05. This is not the old restricted zero/full16 calculation.

Three audited proof refinements improve this all-mask expression:

- Enforcing exact incoming support at return makes those auxiliary
  supports disjoint. The score falls to +3736.19672317.
- Retaining nonzero-syndrome RS output floors lowers it to +3262.21324769.
- Remembering the first live packet after a birth lowers it to
  +3133.11957238. Birth in b input blocks requires at least 17-4b active
  state coordinates; a following packet in b' blocks cannot return when
  b+b'<=4. The small gain shows that immediate excursions are not the
  full remaining loss.

The exponent inherited from the older route-pooled operator is poorly
calibrated. Positive dual arithmetic gives derivative -172170.31 at
alpha=.05 and +16628.65 at alpha=.1. At the latter point the two-state
score is only +137.11622265. Six bounded bisections select the evaluated
witness alpha=119/1280, with score +78.57605118 at z=.85. The same-alpha
newborn score is +71.54683757. These are floating calculations; the
root bracket is not a directed optimality claim.

Rebuilding every local cap at the nearby output marks .84,.85,.86 gives
two-state scores -218.66953880,+78.57605118,+613.63398816. At z=.84 the
newborn refinement gives -227.35829705. The simpler two-state witness
already has about 219 bits of numerical room. It covers all state masks,
all binary schedules, and all 8192 consecutive source sets of size 256,
including both actual alternating matching profiles and honest terminal
handling. It has no quiet-window restriction and uses the implemented
limited refresh distribution. The encoder is unchanged.

This is meaningful progress: the difficult consecutive-support test now
has a negative all-trajectory bound, not merely a negative fixed-path
term. It is still **FLOAT only**. No directed certificate, arbitrary
size-256 support cover, or whole-code distance guarantee follows yet.
The nonquiet 7--7583-group gap therefore remains open. The earlier closed
classes, certified 3.8-ms baseline, and remote preservation bundle are
unchanged. The goal remains active.

Next: verify z=21/25,alpha=119/1280 using directed arithmetic for the
simpler two-state floor cone. Recompute all local caps and exact route
coefficients, not just the final floating score. Then test neighboring
actual profiles such as (128,126), with their correct source counts.
Singleton routing may strengthen the conditional bound while the number
of source sets increases; no monotonicity argument is assumed. Do not
add another construction change or more short-pair machinery first.

The directed replay can omit the optimizer, derivatives, and newborn
states. At rational z=21/25, compute quiet moments and their subtractions
exactly. For a minimum, selecting any independently valid cap is safe;
for a maximum, the directed bound must dominate every candidate. Use
positive polynomial products truncated to full-pair degree128 and
singleton degree2. Keep the live source correction eta=(256/255)^2,
not the quiet-only /255 correction. The per-mask return normalization
remains q^-h q/(q-h), not (q-1)^-h.

The detailed actual-mask cone, exact-support return comparison, affine
coset bound, and first-step restriction are now in
`CHEAP_LANE_RANK_INTERFACE.md`. Independent proof/source reviews passed.
The main agent ran all new cone/return/calibration checks; the separately
audited live-segment run was performed by its implementing subagent.

| Helper under `theory/` | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_live_segment_tilt.py` | 12213 | `6edc953350a67b16dc8dba211494b3edb5a5614aeaa6968c4aa83d14d888d8ac` |
| `cheap_lane_live_cone.py` | 8990 | `f450d2bef7520363011e68e7600e9324b176a8fb20bc9350d059b00e2cdfc1d2` |
| `cheap_lane_disjoint_return.py` | 7808 | `a3eff28456bad943f64932c8f10219ad316d9baede9220c3fcc380d080361abf` |
| `cheap_lane_return_floor.py` | 11130 | `2143482d17afe3624bd0f445773c561a47c60cf7ec45ff8bcef57ff3417b3588` |
| `cheap_lane_newborn_cone.py` | 8762 | `0d4b1d7f403b9bde6a55aaf1d81c502995f69cfd4eb4ab341607c27e5637641a` |
| `cheap_lane_cone_alpha.py` | 10022 | `c0b5976135dc91f317209b9618abb5266ad6d626b74feb89bce7a7a929761700` |
| `cheap_lane_cone_stationarity.py` | 7844 | `ea5a4670f5135a690d8b8fd06a6bce2a4a4b94b8b6fda28b722a8762a6f27996` |
| `cheap_lane_cone_output_anchor.py` | 4838 | `4bac712a2933825de4ade73ad2db3701326894988749d6fc7a79fcfd27c73356` |

## Eighteenth round: directed all-trajectory checks and source runs

The fixed cone admits a much smaller independent implementation. Keep
exact quiet entries, the positive return-floor recurrence, and the
original output-only live/terminal bound. Drop the other feedback
refinements. At z=21/25 and alpha=119/1280 the floating score changes
from -218.6695388033404 to -218.66952702337585, a loss of only
0.00001178 bits. The simpler operator is independently valid; it need
not reproduce the full17-state historical tables.

The new verifier reconstructs every local cap, exact oriented route
count, truncated regional coefficient, finite alternating chain, and
source factor with Arb. Rational variable scales cancel exactly in the
coefficient normalizers. At192 and256 bits, the directed log2 interval
for all8192 consecutive size256 source supports contains
-218.66952702337101557165546909408039234. Both prove a200-bit bound
at35% distance. All fine state masks and binary schedules are covered.
The encoder and limited-support refresh distribution are unchanged.

The256-bit replay exposed a test-interface issue, not an arithmetic
failure: `contains(fmpq)` can ask whether the ball contains a separately
rounded enclosure of the rational. A narrow valid4/25 enclosure failed
that stronger test. A small pinned adapter checks the exact rational
against exact outward endpoints instead. It replaces only the toy test;
the arithmetic verifier remains frozen and unchanged. Independent review
checked the repair and the complete verifier.

A positive boundary-parity recurrence now counts each ordered matching
profile exactly. It preserves indexed rotations and passes5440 exhaustive
small-cycle masks,748 profile comparisons, and42 run-total checks.
For weight256, write full-pair counts as(128-a,128-b). The support has
c=a+b cyclic runs and singleton counts(2a,2b). A singleton truncation
at2c_max therefore covers every profile through c_max runs.

| Covered source supports, all trajectories | Directed log2 union |
|---|---:|
| Exactly one cyclic run | -218.669527023371 |
| At most two runs | -208.301408590818 |
| At most eight runs | -166.933170386442 |
| At most sixteen runs | -133.814902217474 |

The last result includes152 ordered matching profiles. Its192- and256-bit
replays agree and prove the explicit100-bit gate; the displayed directed
interval also lies below-133.8. The largest single contribution is the
sixteen-run profile(a,b)=(8,8), with log2 bound-135.759322655040.
Both orientations and all terminal states are included exactly once.
This first calculation closes a restricted source class. The full cover
below removes the run restriction. No profile monotonicity is assumed.

Carry-safe polynomial packing then makes all marginal profiles affordable.
The radix is2F+1 for retained full degreeF, preventing full-degree overflow
from being misread as singleton degree. Four downward-closed rectangles
recover the129 pairs(full,single)=(128-a,2a), for0<=a<=128. Each rectangle's
exact variable scales cancel before cross-rectangle matrix composition.
The packed verifier reproduces the previous16-run union and passes72
exact toy coefficient checks, including forbidden carries.

A four-vertex directed graph gives a faster exact source count. Vertices
record position parity and activity. A last-exit-tree count sums O(r)
nonnegative terms per profile, avoiding the expanding parity histogram.
The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) gives the formula.
Tests cover5440 masks,1370 all-weight profiles,86 run totals,279 tree
cofactors, and27 comparisons with the first counter.

The full384-bit replay sums all16640 ordered profiles and every state
trajectory for size256 supports. For each cyclic run count c, the profile
counts sum to (8192/c) C(255,c-1) C(7935,c-1). Their grand total equals
C(8192,256). The directed log2 bound is

    -93.5574055068572906550312346546009957059259358406723846743
    (interval radius below3.3e-110).

Thus this source-size class has more than93.55 bits of ideal-setup
failure margin at35% distance. The explicit40-bit gate passes; the100-bit
gate does not. No interval intersection with the nonnegative half-line
was needed. The consecutive-source control reproduces the earlier result.

| Exact cyclic run count | Directed log2 contribution, all its profiles |
|---|---:|
| 16 | -133.968327297632 |
| 32 | -102.309697456541 |
| 43, largest contribution | -96.981220741206 |
| 64 | -112.165984822375 |
| 128 | -297.310090310758 |
| 256 | -1168.270705415009 |

The worst run class is intermediate, not consecutive or fully fragmented.
An independent conservative count uses the total number of c-run sources
times the largest conditional chain for that c. Its combined log2 bound
is-93.5114948669477, agreeing closely here without assuming any general
bound on that counting loss. Both versions include all32767 used refreshes
and the output-only terminal; neither charges an unused final refresh.

This closes every size256 support, not the whole code. The encoder and
refresh distribution have not changed.

A second256-bit pass extracts the entire16641-profile triangle2f+s<=256
once. It then evaluates nine predeclared source weights, with exact
weight-dependent counts and eta_r=q^2/[(q-1)^2(q-r+1)]. It includes odd
weights and uses beta^r, not beta^256. All48152 nonempty ordered chains
pass the exact per-run/source-count tests. The size256 bound overlaps
the384-bit interval above; the size255 bound is

    -105.267775674949284958082024995097323392898376572993055829
    (interval radius below5.1e-71).

Both source sizes therefore close at35% with more than40 bits. The same
fixed z=.84,alpha=119/1280 witness gives positive log2 bounds for
r=7,8,16,32,64,128,192: approximately16576.75,16491.62,15805.83,
14419.75,11620.88,6000.65,875.10. These are uninformative upper bounds,
not counterexamples or arithmetic failures. The fixed output penalty
z^(-alpha L) is particularly wasteful at small r. No interpolation in r
is claimed. The current nonquiet gap is7--254 and257--7583 groups.

For r>256, a separate analytic review identifies the precise scalar
support obstruction and a rank-sensitive extension in the
[rank interface](CHEAP_LANE_RANK_INTERFACE.md). It preserves the existing
cone for initial kernels with at most228 distinct rank-one projection
kernels, uniformly through r7583. Exceptional terms remain open; their
dimension deficit is not an independent probability penalty. Two scalar
shadows give a valid initial comparison but require joint support tracking
locally. None of these new analytic interfaces is a numerical certificate
for an additional source size.

Selected reproduction commands are `python -B` followed by:

- `theory/cheap_lane_cone_interval_tests.py --precision 256`;
- `theory/cheap_lane_run_profile_counts.py --max-runs 6`;
- `theory/cheap_lane_cone_run_verify.py --max-runs 16 --precision 256`;
- `theory/cheap_lane_cone_packed.py --precision 192`;
- `theory/cheap_lane_fast_profile_counts.py`;
- `theory/cheap_lane_cone_profile_cover.py --precision 384`;
- `theory/cheap_lane_cone_weight_pilot.py --precision 256`.

Next retune z toward1 at smaller source weights, using smaller coefficient
rectangles and exact counts. A bounded anchor schedule can start from
z_r=1-r/1600, but each selected witness still requires verification.
Before implementing a two-shadow operator, also test the existing
bounded-missing-group comparison near257 with calibrated marks. At the
present alpha its one-missing-group penalty is about95.3 bits; the
unmodified size256 margin is93.6 bits, so a40-bit certificate does not
follow. Neither tuning suggestion is a construction change or a claim
about untested source sizes.

| Frozen round18 helper | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cone_simple.py` | 7559 | `79df00e5bbae3f2e9e533166091fe7e1103d6621eeba2a7c89dd9057c683f44d` |
| `cheap_lane_cone_verify.py` | 14581 | `7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01` |
| `cheap_lane_cone_interval_tests.py` | 2470 | `83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478` |
| `cheap_lane_run_profile_counts.py` | 8186 | `25e0d9221a62ffdb79d1d625e6e858ffbb8c0b671078d2bde5d2c6fb27e43ed3` |
| `cheap_lane_cone_run_verify.py` | 6717 | `afdfd902dac457c6f49e93dbe19d9841faed1db9a4ed0895f64a2ad31133e9cc` |
| `cheap_lane_cone_packed.py` | 11055 | `9379de72b355d04e646a52bd560f538d77cf7095e0e8785b7e5bf982a1292293` |
| `cheap_lane_fast_profile_counts.py` | 8185 | `4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8` |
| `cheap_lane_cone_profile_cover.py` | 10476 | `84ac66e9e3b49b854ce10adf2b2f244d74a530524252736d35d1e72cdfe6e7fb` |
| `cheap_lane_cone_weight_pilot.py` | 8680 | `dd982dac03e11a533df5c72930430f1bb019e8822b66b5a63d6513df04cb6e2b` |

## Nineteenth round: size-dependent witnesses and a contiguous small band

The encoder and setup distribution remain unchanged. The fixed round18
witness was badly calibrated for sparse sources. A new bounded verifier
accepts exact rational z and alpha, sizes its carry-safe coefficient
triangle to the actual source weight, and retains every matching profile.
It applies the corresponding source factor only once. Exact boundary
tests retain degree r and reject overflow, including full degrees above128.
Independent source review found no blocking issue.

At r=7, the initial z=1-r/1600 and alpha=119/1280 give log2 bound
+238.430. Moving z to1-r/1280 and then increasing alpha to1/2 gives
-179.25333948392939. At r=32, z=39/40 and alpha=3/16 give
-644.48151188837153. These are proof parameters, not construction
parameters or changes to the randomness distribution. Increasing alpha
is not uniformly beneficial: at r=128,z=.92, changing alpha from119/1280
to1/8 worsens the bound from-28.2471 to+166.2594. Different source profiles
can become dominant, so neither a universal alpha rule nor profile
monotonicity is inferred from these tests.

The fixed band verifier uses z_r=1-r/1280, alpha=1/2 for7<=r<=15,
and alpha=3/16 for16<=r<=32. All26 directed256-bit computations pass.
For each printed log2 interval with upper endpoint u, it replaces the
probability bound by2^(-floor(-u)), then adds those dyadic bounds exactly.
The resulting joint bound is below2^-178. This covers every indexed
source support and every state trajectory with7--32 active groups at35%
distance; it is not merely26 individually passing40-bit tests.

The conservative per-weight integer margins, in increasing source size,
are

    179,419,659,843,929,1013,1097,1181,1164,404,430,455,481,
    507,532,558,583,608,633,659,683,708,725,702,675,644.

The exact joint dyadic sum has denominator2^1181 and a1003-bit numerator.
Its rounding, serialization, deadlines, source pins, and completion gates
passed an independent audit. Every child checks its per-run counts and
their sum C(8192,r), and includes32767 used refreshes and the output-only
terminal. The complete band run took63.78seconds and generated no files.

An additional256-bit witness at r=128,z=91/100,alpha=119/1280 gives

    -196.40889813805544234468109983245069395728736438857706235
    (interval radius below6.8e-72).

All4224 ordered profile chains are included. The largest run class has
39runs and contributes log2 bound-199.685610046331. The run took24.27seconds;
all count tests pass and no nonnegative interval intersection is needed.
The neighboring z=.90 and z=.92 witnesses are substantially weaker,
confirming that output-mark calibration matters here.
A second,384-bit arithmetic replay reproduces the r=128 interval
with radius below6.0e-110 in30.18seconds. This replays the complete4224
profile sum, not merely the displayed logarithm.

The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) now also proves two
extensions. First, the one-missing-group comparison for257<=r<=511 can
feed the refined cone if padding happens before actual-support filters.
Its complete source factor is

    beta^r*(1+r*q^256/beta) / ((q-1)*(1-r/(2q))).

Second, paired regional matrix envelopes can cover cells of source sizes
and run counts without interpolation. Both derivations passed independent
review. They license numerical checks; neither alone closes another size.

The first one-missing-group check is now complete at r=257,
z=1343/1600,alpha=119/1280. The256-bit directed log2 union is
+0.89564762496111756 (radius below8.4e-72), over all16641 ordered
profiles. The largest run class has43runs. All exact counting and terminal
tests pass; the run took158.65seconds within its180-second budget.
This is a failed distance bound, not a counterexample. The extra padding
cost is now an explicit proof bottleneck immediately above256 groups.
Reducing that loss by about41 bits would meet the40-bit gate at this
single source size, but would not close any untested larger size.

### Structural improvement: average the missing outer maps last

The preceding failure led to a new independently audited source lemma.
Retain the condition that each missing group has a nonzero projection in
the initial kernel. Conditional on the other outer maps, the missing
maps remain independent. Their scalar outer counting bound charges for
the necessary kernel intersection and cancels the ambient rank increment
that caused the padding penalty. This compares the same positive
actual-T kernel before the fractional power; all local cone bounds and
filters remain unchanged.

For one missing group, even the conservative factor

    Omega_r = beta^r*(r+1) / [(q-1)*(1-r/(2q))]

suffices. Its ratio to the old factor at r=257, raised to alpha, reduces
the log2 bound by95.18268059806034. The adjusted directed result is

    -94.28703297309922250062530870719105672365788703283799616
    (interval radius below3.4e-71).

Thus every size257 source support and every state trajectory is now
certified at35% with more than94 bits. A small separately pinned receipt
checks the algebraic adjustment at256 and384 bits. It reuses the completed
16641-profile computation above; it does not claim a fresh coefficient
recomputation or a tighter enclosure than that original input interval.

The proof also handles a missing set jointly. If at most q groups are
missing, a kernel with nonzero projections onto each contains a scalar
vector nonzero on all of them. Average these outer maps together, not
as separately conditioned events. The complete source factor becomes

    beta^r / [(q-1)*c_h]
      * [1 + sum_(j=1)^h C(r,j)/(q-1)],
    c_h=1-r/[q*(h+1)]>0, h<r, h<=q.

For r<=8192, h=floor(r/q) is at most32 and satisfies these conditions
when used above256. Two independent proof reviews checked the joint
argument and its transcription in the rank interface. This removes the
old source-comparison obstruction, but the new larger-size regional
bounds have not yet been evaluated. No new distance claim for those
sizes follows from the lemma alone.

Selected reproduction commands are `python -B` followed by:

- `theory/cheap_lane_cone_small_band.py --precision 256`;
- `theory/cheap_lane_cone_weight_witness.py --weight 128 --z 91/100 --precision 256`;
- `theory/cheap_lane_cone_weight_witness.py --weight 257 --allow-missing-one --precision 256`;
- `theory/cheap_lane_missing_group_receipt.py --precision 384`.

At the end of round19 the nonquiet gap was33--127,129--254,and258--7583. Prior quiet,
one-group, two-through-six,255--256,and dense-band certificates remain
preserved. These gains constitute progress, not a whole-code certificate.
Next integrate the joint missing-group comparison into the size-cell
verifier and build contiguous coverage. Do not optimize already-large
margins at isolated source sizes or infer untested sizes from this result.

| Frozen round19 helper | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cone_weight_witness.py` | 12607 | `e32069c849efa8e3abe1d3eee0a252a386e8acb3c161725d9ec83314a7b1678b` |
| `cheap_lane_cone_small_band.py` | 9905 | `32c3feb9245bca524855ec82fa713d8b7d0b90840e5a040b4f06fee969a8c1d4` |
| `cheap_lane_missing_group_receipt.py` | 5131 | `13bd94571c1d54d6e695d410df96e530ab97852611edc43b94297278d34d2cf8` |

## Twentieth round: contiguous support bands and calibrated route markers

The encoder, outer distribution, and cheap lane-refresh distribution are
unchanged. A new verifier shares regional coefficient calculations across
a contiguous band of source sizes. It extracts every required marginal
coefficient, forms the audited paired envelopes, and sums exact indexed
run counts. It uses the joint missing-group source comparison above256.
Each numerical result below covers all source supports and all state paths
in its stated band, not only sampled profiles or endpoint sizes.

All witnesses use rho=sigma=1 and singleton size/run cells. Directed256-bit
arithmetic gives the following conservative integer margins. A margin b
means the joint band probability is at most2^-b. The bands64--80 and49--64
overlap at64; summing both remains a valid upper bound.

| Actual source sizes | z | alpha | Conservative joint margin |
|---|---|---|---:|
| 33--34 | 39/40 | 3/16 | 639 |
| 35--48 | 97/100 | 1/8 | 106 |
| 49--64 | 24/25 | 1/8 | 313 |
| 64--80 | 19/20 | 1/8 | 188 |
| 81--88 | 47/50 | 119/1280 | 172 |
| 89--104 | 163/175 | 119/1280 | 92 |
| 105--120 | 327/355 | 119/1280 | 132 |
| 121--136 | 41/45 | 119/1280 | 141 |
| 137--152 | 329/365 | 119/1280 | 132 |
| 153--168 | 33/37 | 119/1280 | 119 |
| 169--184 | 331/375 | 119/1280 | 104 |
| 185--200 | 83/95 | 119/1280 | 86 |

Thus the previous small-source certificates now extend continuously
through200 groups. Together these give more than43 bits for the union of
all1--200-group events at35%. There are19 component bounds: the one-group
bound is at most2^-44, the six-group bound at most2^-45, and the other17
are each at most2^-50. Their sum is at most113*2^-50<2^-43. The harmless
double counting at64 does not invalidate this union bound.
At192 groups, the separately pinned exact-weight helper
with z=7/8,alpha=119/1280 gives log2 bound-154.31349577608313.
The contiguous-band cover, not that isolated result, covers its neighbors.
Every cell run checks
that its support counts sum to C(8192,r) at each r. Exact small-cycle,
packing, terminal, and cell-envelope tests pass. No nonnegative interval
intersection was required in these selected runs.

Two unsuccessful band witnesses explain why the cover uses several marks.
At33--48, z=.97,alpha=1/8 gives log2 bound+126.763; retuning33--34 closes
those sizes while the unchanged witness closes35--48. At81--96, z=.94
gives only10.221 bits. Narrowing the first band and moving z for89--104
closes both. These are failures of fixed proof witnesses, not codewords
showing low distance. At r64, complete-chain controls show that paired
envelopes cost about2.30 bits at z=.95 and3.61 bits at z=.96.

### Extend the dense endpoint without changing the local kernel

The joint missing-group comparison also sharpens the frozen dense bound.
Use q=256, h=32, z=3/1250, g=z+(1-z)/q, A=beta*g^256<1, and

    J=1+sum_(j=1)^32 C(8192,j)/255,
    rho=(1+1/q)^16+65536*q^-8.

Since c_h>=1/33, the joint failure probability for all supports7552--8192
is bounded by

    2^8192 * (33/255) * A^7552 * J * q^8
      * rho^32767 * z^-734003.

The local output/feedback majorant and final-packet treatment are unchanged.
Directed256- and384-bit evaluations both give log2 bound
-781.46273555593582010034250920066589371886893339917364016452291442873664165.
This replaces the previous dense endpoint7584 in the current coverage.

### Diagnose larger-support coefficient bounds before expanding them

Four fixed positive coefficient bounds avoid large regional tables but
fail at their initial route markers. They cover one actual equal-matching
profile each, not every profile at the listed size:

| r | Each profile (full,singleton) | z | Initial log2 upper bound |
|---:|---|---|---:|
| 512 | (16,480) | 7/10 | +12102.3821 |
| 1024 | (64,896) | 1/2 | +18139.9099 |
| 2048 | (256,1536) | 1/4 | +33219.9140 |
| 4096 | (1024,2048) | 1/16 | +41470.4544 |

Here alpha=119/1280 and the initial markers are the unweighted categorical
centers. The [rank interface](CHEAP_LANE_RANK_INTERFACE.md) derives the
finite-chain coefficient bound and its exact marker derivatives. At r512,
the tilted means per region are6.21486249 full pairs and305.62658510
singletons, versus the required16 and480. The derivatives of log2 bound
with respect to ln(u),ln(v) are-3613.94416 and-64401.32120. These directed
results show that the initial markers are not stationary. They justify
calibration, not by themselves a larger-support certificate. The derivative helper also
reproduces the original bound with an independent regional assembly.

The bounded calibration then takes seven FLOAT proposal evaluations and
rounds the selected markers to u=6943/990029 and v=152028/880327. All other
parameters are unchanged. Rebuilding the full finite-chain bound from
exact route weights in directed arithmetic gives

    log2 U = -107.54468059066184461847619149009820523220523531929169.

Both256- and384-bit replays agree; the latter interval radius is below
4.3e-107. The selected profile therefore has more than107 bits at35%
distance. The earlier positive bound lost12209.9268 bits solely from
marker miscalibration. This is not a claim about every512-group profile.
The proposal optimizer supplies no proof inputs except the final explicit
rationals, whose complete bound is evaluated again with directed arithmetic.

Next cover every512-group profile using a small library of calibrated
marks. The rank interface proves the profile-cell bound, including the
ordered finite terminal chain and exact run-count accounting. Marker reuse
does not add a new Jensen step. Rare worst-case profiles may still make
run-count times maximum too loose; exact marginal-count capacities give
an audited inexpensive refinement. Do not widen to larger sizes before
testing this complete fixed-size cover.

A small constant-marker library can itself lose too much through marker
mismatch. A bounded alternative is to calibrate a few anchors, interpolate
only their correction to the natural markers, and round a witness for
each full count f. Reevaluate each resulting rational witness; do not
interpolate probability bounds. At r512 there are257 marginal full counts
and66048 degree-feasible ordered pairs. A cached local polynomial table
can avoid large regional arrays. Exact profile counts are preferable if
run-count maxima assign excessive weight to rare imbalanced profiles.
The cheapest proposed first pass needs only the successful anchor:

    u_f=(1562175/990029)*f/(3584+f),
    v_f=(1140210/880327)*(512-2f)/(3584+f), 0<=f<=256.

Replace an endpoint's zero marker by2^-30. This curve reproduces the
selected f16 markers and tracks the natural categorical proportions.
At the end of round20 this was an untested witness proposal, not an assumption about optimality
or probability interpolation. Recompute each rational witness and every
compatible ordered chain; refine only where the resulting cover fails.

Reproduce a selected band with `python -B theory/cheap_lane_cone_cells.py`
and `--rlo LO --rhi HI --z Z --alpha ALPHA` from the table. Other commands:

- `theory/cheap_lane_cone_weight_witness.py --weight 192 --z 7/8 --alpha 119/1280`;
- `theory/cheap_lane_dense_joint.py --precision 384`;
- `theory/cheap_lane_cone_cauchy_anchors.py --precision 256`;
- `theory/cheap_lane_cauchy_marker_gradient.py --precision 256`.
- `theory/cheap_lane_cauchy_marker_calibration.py --precision 384`.

| Frozen round20 helper | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cone_cells.py` | 15533 | `0ffd65ea35b6797715b4f55c7fde00b7985e780aa60ca94266259634de0fa11a` |
| `cheap_lane_dense_joint.py` | 3216 | `494a6a9cf8514cb11b8109e719e0d493a39a2dea3e37f86fe3026d420759d86f` |
| `cheap_lane_cone_cauchy_anchors.py` | 11616 | `319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9` |
| `cheap_lane_cauchy_marker_gradient.py` | 10452 | `2bdf98517687e3fa549f412bab29db45752e3b4e6b5e8984f50838c9916f5631` |
| `cheap_lane_cauchy_marker_calibration.py` | 14119 | `f0358e6296357f13e1d0b18b676e06b9672a1ac4ecbb4bbe65872a3391e8054e` |

## Twenty-first round: the complete512-group screen identifies clustered profiles

The previous round made verified progress. This round first extends the
exact contiguous source-size cover through216. The encoder and all frozen
helpers remain unchanged. With z=333/385 and alpha=119/1280, the shared-cell
verifier covers every201--216-group support and all state paths at35%:

    log2 U = -64.31454722444925336858501242026791187012406202139371545
    (directed256-bit interval radius below4.5e-71).

The worst cell is r216 with44 runs. All exact count and terminal tests
pass; the run took99.44 seconds. Adding this band to the previous19 bounds
gives at most114*2^-50<2^-43 for every source size1--216 jointly.

A second directed band, with z=267/311 and the same alpha, covers every
217--224-group support and every state path at35%:

    log2 U = -101.92925218481550128448163117337481492179777487928139595
    (directed256-bit interval radius below2.1e-71).

Its worst cell is r224 with43 runs. The run took108.11 seconds, checking
888 marginals and1764 source-size/run-count pairs. Both new bands pass
the exact count and terminal checks, with no nonnegative interval
intersections. The current continuous cover is1--224, jointly more than
43 bits: the previous bound113*2^-50 plus these two bands is at most
115*2^-50<2^-43. No source-size interpolation is used.

### Evaluate the entire512-group union before adding refinements

The new all-profile screen evaluates all257 marginal witnesses from the
recorded correction curve. It then sums all66048 ordered chains with exact
indexed source multiplicities. Only source-count symmetry is reused;
the two terminal orientations remain distinct. The summed counts equal
C(8192,512), and every cyclic-run subtotal matches its exact formula.
No whole-profile run-total upper bound is substituted for these counts.

The completed FLOAT screen gives log2 union+2916.2331536374513. This is
not a certificate or a low-weight codeword. The largest individual term
has f_0=f_1=223, s_0=s_1=66, hence66 runs, and log2 value+2910.52673407.
The67-run union is slightly larger than the66-run union. In contrast,
the480-run union is-108.62955661, near the previously calibrated profile.

| Selected run count | FLOAT log2 union |
|---:|---:|
| 1 | +2033.3865 |
| 32 | +2838.0140 |
| 67 | +2912.8921 |
| 128 | +2793.2745 |
| 256 | +2111.6350 |
| 384 | +1053.6759 |
| 480 | -108.6296 |
| 512 | -1356.1041 |

The dominant inputs in this calculation are clustered, not the fragmented
profile that supplied the first marker calibration. Next separate marker
mismatch at the dominant profile from coefficient-bound slack. The
[rank interface](CHEAP_LANE_RANK_INTERFACE.md) proves a positive regional
degree-residue refinement that can remove some coefficient slack without
changing the encoder or introducing another Jensen step.

### Calibrate the clustered profile, then measure coefficient loss

For the actual equal profile f_0=f_1=223,s_0=s_1=66, retain z=7/10,
alpha=119/1280, and its exact indexed source count. Five FLOAT proposals
select positive rational markers u=29411/311125,v=10441/460035. A fresh
directed256-bit replay gives log2 bound+2760.1551811470249958720, versus
+2910.5267340694476332030 initially. The gain is only150.3715529224 bits.
The selected tilted regional means are223.00001194 and66.00000105,
nearly stationary at this fixed z and alpha. This does not assert a global
optimum. The larger run-count control changes the log bound by only
2.79675 bits, so source-profile counting is not the dominant loss here.
The replay took4.15 seconds with a55,435,264-byte measured process peak.

At these same rational markers, the audited positive degree-residue
calculation gives:

| Regional degree modulus | Directed log2 bound | Gain over Cauchy |
|---:|---:|---:|
| 1 | +2760.1551811470 | 0 |
| 4 | +2248.1551811470 | 512 bits |
| 8 | +1992.1551811470 | 768 bits |
| 16 | +1736.1551811470 | 1024 bits |

These are selected-profile majorants, not certificates or counterexamples.
The residues are computed, not assumed uniform. Each completed region is
filtered before composing the256-region chain. All aliases are retained;
the terminal counts its input route but performs no final update. Exact
small polynomial tests and the modulus-one/sum-of-residues controls pass.
The entire directed run took1.80 seconds with a39,256,064-byte measured
process peak. Retaining both full and singleton degrees is the next
coefficient diagnostic; gains need not continue linearly once the modulus
exceeds the width of the tilted degree distribution.

### Two occupancy degrees remove most, but not all, of the diagnosed loss

The new bivariate positive-residue calculation retains full and singleton
counts separately. It passes independent source/proof review,60 exact
packed-product tests,561 exact local degree normalizers, and342 complete
regional collapse checks against the frozen one-axis implementation.
At the same f223 markers, z=7/10 and alpha=119/1280, directed256-bit
bounds are:

| Modulus on each occupancy count | Log2 bound | Gain over ordinary Cauchy |
|---:|---:|---:|
| 8 | +1224.1551813107 | 1535.9999998363 bits |
| 16 | +712.0656154537 | 2048.0895656933 bits |
| 32 | +311.4317631894 | 2448.7234179576 bits |

The default run takes3.85 seconds with a45,555,712-byte measured process
peak and no nonnegative interval intersections. It is still a positive
majorant for one selected profile, not a new certificate. The exact
source count, source comparison, and terminal remain unchanged.

A three-point fixed-marker sensitivity then tests only
z=224/325,7/10,91/128. Its modulus32 log2 values are respectively
+626.7776403203,+311.4317631894,+1031.9291582103. One safeguarded
quadratic proposal, rounded to z=17447/25000, gives+304.3534420771
in directed arithmetic. The approximate prediction was+291.6233;
that prediction is not used as evidence. Route markers are stationary
only for the earlier ordinary Cauchy objective, not automatically for
this residue-projected objective. No envelope-theorem claim is made.

At that selected z, the bounded alpha comparison gives:

| Alpha | Modulus32 directed log2 bound |
|---:|---:|
| 119/1300 | +307.3477355073 |
| 119/1280 | +304.3534420771 |
| 1547/16384 | +337.6048202059 |

These are fixed-marker sensitivities, not an optimized family of bounds.
The local caps are independent of alpha, so the underlying positive
finite-path expression is log-convex in alpha at fixed z,u,v and modulus.
This fact is stronger than the local quadratic approximation in z.
Independent review confirms a global obstruction to alpha-only tuning
of this fixed expression. Write a_-=119/1300,a_0=119/1280,
a_+=1547/16384 and f_i=log2 U(a_i). Their two secant slopes have
opposite signs. Outside[a_-,a_+], convexity gives f>=min(f_-,f_+).
Inside, extend the opposite secant line to obtain

    f(alpha) >= min(f_-, f_+, (129 f_0-64 f_+)/65,
                                 (129 f_0-65 f_-)/64).

Conservative outward decimal endpoints from the directed runs give
the last two lower bounds exactly13580681180591/50000000000 and
192839912199649/640000000000. Exact rational arithmetic verifies that
both exceed271. Thus no admissible alpha makes this fixed, unclipped
proof expression negative. This does not constrain new route markers,
finer residues, exact coefficients, or the actual distance of the code.
The selected calculations take at most4.12 seconds each and less than
46 MB measured process peak. All are serial, directed, and stdout-only.
Next isolate the remaining alias mass with one modulus64 or exact
selected-coefficient calculation before expanding the profile search.
Neither residue flattening nor a saddle heuristic bounds the possible
remaining coefficient gain.

The screen took48.12 seconds with a measured49,991,680-byte process peak.
It wrote no files. Reproduction uses `python -B` followed by
`theory/cheap_lane_cauchy_profile_screen.py`. The sparse band uses
`theory/cheap_lane_cone_cells.py --rlo 201 --rhi 216 --z 333/385 --alpha 119/1280`.
The second band uses
`theory/cheap_lane_cone_cells.py --rlo 217 --rhi 224 --z 267/311 --alpha 119/1280`.
The clustered-profile checks use
`theory/cheap_lane_cauchy_profile_calibration.py --full 223 --run-count-control`
and `theory/cheap_lane_cauchy_modular.py --full 223 --u 29411/311125 --v 10441/460035 --include-16`.
The bivariate default is `theory/cheap_lane_cauchy_bivariate_modular.py`.
Use its `--z` and `--alpha` options for exactly the rational values above;
the f223 route markers remain fixed by default.

| Frozen round21 helper | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_profile_screen.py` | 13590 | `6cd1f551da748c461c995f7927144345655cc741156d0550625388fcd6876861` |
| `cheap_lane_cauchy_profile_calibration.py` | 10523 | `3ef21698a0e245466c4d7f4b5d2cb83b940e0b8d5b35dc2bdece9bcc90315170` |
| `cheap_lane_cauchy_modular.py` | 11639 | `128ac26691d4e046267e277516c07594094ae7130e8777f3b1d35f99e4e1f014` |
| `cheap_lane_cauchy_bivariate_modular.py` | 15119 | `01356ae20b8b00a3354d76c4c7978bcb3d77d73293b110ea2a9cd4dc9b495d95` |

## Twenty-second round: exact extraction closes the clustered profile

The previous round made verified progress. This round tests whether the
remaining clustered-profile gap comes from coefficient aliases. It keeps
the actual refresh distribution, local caps, source comparison, and encoder
unchanged. The source set has512 active groups and profile(f,s)=(223,66)
in both alternating matchings. Its exact indexed multiplicity is charged
once. All state trajectories and the output-only terminal are included.

At z=17447/25000 and alpha=119/1280, with the previous positive route
markers, the directed256-bit calculations give:

| Regional coefficient method | Log2 bound |
|---|---:|
| Residues modulo32 on both counts, control replay | +304.3534420771222 |
| Residues modulo64 on both counts | +44.0843077780031 |
| Exact coefficient | -82.3387645228212 |

The exact calculation retains the downward-closed rectangle
0<=f<=223,0<=s<=66 after each positive product. Packed stride447 prevents
a discarded full degree from carrying into the singleton degree. Seven
doublings computeP^128 andP^127T. Extracting coefficient(223,66), then
normalizing and composing the256 regions, removes all degree aliases.
This is an exact coefficient of the existing local-cap upper bound,
not the code's exact failure probability.

Independent proof/source review passes. The wrapper checks89 exact
rectangle-product coefficients,60 cyclic-product coefficients,561 exact
local degree normalizers, the indexed source counts, and the terminal.
The residue64 run takes11.26 seconds with a75,309,056-byte process peak.
The exact256-bit run takes56.21 seconds with a198,328,320-byte peak and
no nonnegative interval intersections.

### Independent precision and scale replay

For exact extraction, every regional target coefficient acquires precisely
the factoru^f v^s. Its normalization cancels that factor. Therefore the
exact bound is independent of positive route scales, unlike a residue
bound that still includes aliases.

A separate384-bit replay changes the scales to u=3/32,v=3/128. It gives

    log2 U = -82.33876452282115505512336387669992409691654620662418
    (interval radius below7.4e-110).

The replay overlaps the256-bit result and takes61.66 seconds with a
242,810,880-byte measured process peak. Both explicitly prove U<2^-40.
The unchanged ordinary Cauchy evaluation moves from+2759.8441 to+2796.1268
under the scale change, confirming that the cancellation is specific to
the exact coefficient. The selected clustered class is now closed with
more than82 bits at35% distance. No complete size512 or whole-code
certificate follows from this single profile pair.

Reproduce from the repository root with `C:/Python314/python.exe -B`
followed by:

- `research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_fine.py --control32`;
- `research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_fine.py --mode exact`;
- `research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_fine.py --mode exact --precision 384 --u 3/32 --v 3/128`.

| Frozen round22 helper | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_fine.py` | 12339 | `e4fcb2bbb4a73183d31f52ab2b86d5ee9e4906045dd045158fd75966d9aba1a2` |

### Complete source-size screen

The next screen retains both occupancy residues for all257 marginal
profiles, then sums all66048 ordered profile pairs. It uses the same
z,alpha, local caps, indexed counts, and output-only terminal. The fixed
positive marker curve is the round21 curve; its endpoints use2^-30.
FFT evaluation is a diagnostic, not directed verification.

| Modulus in each occupancy count | FLOAT log2 union | Positive profile terms | Worst ordered profile |
|---|---:|---:|---|
|64|+398.9453696442|5222|(250,12),(250,12)|
|128|-59.9351744667|0|(235,42),(235,42)|

The residue128 worst profile contributes log2-66.0487306221. Every
source is counted: the exact counts sum to C(8192,512), and each exact
run total agrees independently. Six noncommuting terminal tests and480
exact toy coefficients pass. Each marginal's zero-frequency matrix
replays ordinary evaluation, and the inverse FFT passes per-entry
imaginary, negative-noise, and total-mass checks without clipping.

The residue64 selected-profile reference uses the curve markers
u=38707225/418782267,v=929060/41375369, not the earlier optimized scales.
A directed256-bit run gives log2+230.9064716006121906391795876;
the matching FLOAT result is+230.90647160060277. The directed reference
takes10.11 seconds with a74,059,776-byte process peak. Complete residue64
and128 screens take51.08 and57.82 seconds, with measured peaks50,569,216
and55,820,288 bytes. They produce no files. The128 screen gives a viable
whole-size target, but does not yet certify that source size.

Reproduce either screen with `C:/Python314/python.exe -B` followed by
`research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_residue_screen.py --z 17447/25000 --modulus 64`
or the same command with `--modulus 128`. The directed64 reference uses
`theory/cheap_lane_cauchy_fine.py --u 38707225/418782267 --v 929060/41375369`
from this workstream. The next step reuses a fixed marker bank, as justified
in the rank-interface note, before directed complete-size verification.

| Frozen round22 screen | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_residue_screen.py` |14221| `0b28a10417eec7a1f6669137f2f12fc805e1b694388777fd81040d1c816a1ed9` |

### Fixed marker bank

The audited cell wrapper assigns each marginal to its nearest fixed
marker center; ties use the smaller center. It extracts the target
residue from that center's table and uses the target's own denominator.
The exact target term cancels its marker powers. Only additional residue
terms depend on the chosen center; this is not interpolation.

Ten centers(16,48,80,112,144,176,208,224,240,252) give FLOAT log2 union
+16.3895550145. The only positive profile is(246,20) in both matchings,
assigned to center240. That approximation is too loose. Adding centers
232,236,244,248,254,256 yields

    centers=16,48,80,112,144,176,208,224,232,236,240,244,248,252,254,256
    complete r512 FLOAT log2 union = -59.24900641557496.

All66048 ordered terms are negative. The largest is the pair(234,44),
(234,44), with log2-64.7909120505. Both runs retain all257 marginals and
pass the exact count, terminal, scaling, and FFT controls. The selected
minimum coefficient/DC ratio is3.80e-9, above the enforced1e-10 floor;
no coefficient is clipped or omitted. The ten- and sixteen-center runs
take49.10 and50.40 seconds with peaks53,100,544 and53,170,176 bytes.

Run `theory/cheap_lane_cauchy_cells.py` with `C:/Python314/python.exe -B`
from this workstream for the failed ten-center control. Add `--centers`
and the sixteen-center list above for the successful FLOAT screen.
A separately audited directed replay is still required before treating
this numerical margin as a complete source-size certificate.

| Frozen round22 marker-cell screen | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_cells.py` |12389| `0ed129b48ee91cfc456c2e434ea5e279bd6ad777358f836d6152d2d2c9c813b1` |

The first directed bank implementation passed two independent source
reviews. Its first two cells completed in96.81 and162.15 seconds cumulatively,
with a263,974,912-byte measured process peak. The parent stopped the job
after confirming its process identity: sixteen such tables appeared unlikely
to finish within the declared900-second budget. No final sum or certificate
was produced, and no generated file needed deletion. This was a resource
decision, not an arithmetic failure or a failed distance bound.

The next plan uses positive truncation to avoid unnecessary coefficients.
For centers16,48,80, keep exact full degrees through32,64,96 respectively
and retain singleton degrees modulo128. Keep the four cyclic tables at
centers112,144,176,208. Replace all nine remaining tables by one exact
rectangle F<=256,S<=78, which contains every target for217<=f<=256.
The exact rectangle can use center232's scales because normalization cancels
them. Both changes remove positive aliases and cannot weaken the passing
sixteen-center bound. This comparison passed independent proof review;
the new directed computation remains necessary.

| Frozen round22 directed bank attempt | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_cells_directed.py` |12648| `3e4fe22d9caf9419f05d29ec001a7468b3a64a085844e34ac008f5537f439768` |

The hybrid verifier implements this eight-band plan without changing any
frozen helper. Two independent source reviews pass. Its new mixed-axis
tests compare unwrapped rational products against exact full-degree
truncation and cyclic singleton folding. The frozen cyclic and rectangle
tests remain in the replay. All257 marginals are assigned exactly once;
the ordered count sum and terminal checks remain unchanged. Reproduce
with `C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_hybrid_directed.py --precision 256`
from the repository root. Its fixed budget is900 seconds and512 MiB,
with a20-GiB free-space preflight and no generated files.

| Frozen round22 hybrid verifier | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_cauchy_hybrid_directed.py` |14695| `985c3fe04f055811976d8060a4b5a93063a6cabd81dd6f4c285769f60e096727` |

### Directed complete-size result

The256-bit hybrid replay completes with

    log2 U_512 = -60.86637479476222931534153913373382070564090522274946
    (interval radius below2.5e-71).

Thus, at35% distance, the ideal-setup failure contribution from all
messages supported on exactly512 actual source groups is less than2^-60.
The sum includes all66048 ordered matching profiles and all state paths.
It is not a whole-code certificate. The exact indexed source counts sum
to C(8192,512), and all512 run-count sums agree. The source factor uses
h=2,c_h=1/3,eta=43861/7225 exactly once before its fractional power.

The largest individual term has f_0=f_1=234 and log2 bound-66.8641086165.
The largest run-count subtotal has44 runs and log2 bound-64.3620171377.
The selected exact(223,66) pair gives-82.3387645228211550551233638767,
overlapping both earlier independent-precision and scale replays. There
are no nonnegative interval intersections. Checks include32 mixed-axis,
89 rectangle, and60 cyclic exact toy coefficients,8 terminal-order checks,
561 local route normalizers, and24 cyclic plus24 truncated mass controls.

The complete run takes411.89 seconds, with a273,104,896-byte process peak,
within its900-second/512-MiB budget. Its witness digest is
`6ab2ef3d2a972c4f0200a9ce99207e0fbd88822f6837d711a669431af78fbad3`.
No generated data is written. The improvement changes coefficient counting,
not the encoder or its setup distribution.

Next, test the conservative source-size/run-count envelope at512 before
extending a short neighboring band, for example508--516. The packet
polynomial can be reused at fixed z,alpha, but exact profile counts currently
matter: replacing them by run-total maxima may consume the margin. Refine
profile cells if that control fails. The remaining nonquiet source sizes
are225--254,258--511,513--7551; neither35% for the whole code nor a higher
whole-code distance is established.

## Twenty-third round: reuse coefficients without losing exact counts

The previous round made progress: its directed calculation closed every
size512 support, not just selected profiles. This round tests whether the
same coefficient computation can cover neighboring sizes. The encoder,
local caps, z=17447/25000, and alpha=119/1280 remain unchanged.

### The run-count shortcut loses too much

The new FLOAT screen reuses a fixed residue128 marker bank. It compares
two sums built from identical conditional chains: the exact profile-count
sum, and the sum of each run count times its largest conditional chain.
At size512 the results are

    log2 exact-profile sum = -59.24900641551881,
    log2 run-count envelope = -32.13045484730855.

Thus the counting relaxation costs about27.12 bits at this witness. The
dominant envelope uses highly unbalanced profiles, such as(256,216) for
40 runs, whereas the exact-count sum is led by profiles near(234,234).
The coarse run-count envelope is not suitable for the40-bit target.
These are FLOAT diagnostics, not new certificates. The directed hybrid
size512 certificate remains60.8663748 bits.

The explicit15-center bank ends at254. Omitting center256 avoids its
nearzero singleton marker when extracting positive-singleton neighboring
targets. No failed coefficient is skipped or silently retried. Each
selected target retains its own route denominator. The size512 screen
takes49.77 seconds and peaks at52,609,024 bytes; all exact support/run
counts and FFT controls pass. No files are generated.

A second serial screen extracts508,513,516 from the same tables. Its
complete exact-profile FLOAT sums are respectively-75.7534598358,
-54.0669635076, and-36.1703189190 in log2 scale. The size516 witness
therefore misses40 bits at this stage. Its run-count envelope is worse:
+21.6838178473, led by the rare(0,258) profile in the258-run class.
The exact sum for that run class is only2^-752.4023. This isolates the
counting loss rather than suggesting a low-weight codeword. All three
sizes pass their exact count checks, including the odd size513. The run
takes140.17 seconds, peaks at52,609,024 bytes, and generates no files.
These three diagnostics do not cover uncomputed intermediate sizes.

The rank interface now states an independently audited alternative:
group ordered triples(r,f_0,f_1), retaining their exact multiplicities.
Positive matrix envelopes can then bound cells without replacing a rare
profile's count by its entire run-class count. Singleton cells recover
the existing majorant. The first directed band instead keeps the full
exact profile sum; no unmeasured aggregation loss is accepted.

A simpler independently audited alternative is now stated in the rank
interface: sort conditional profile bounds, compute exact counts from
the largest terms downward, and retain the exact residual run-class count.
Its tail is bounded by that residual count times the next largest proved
chain endpoint. The post-source tail budget bounds the additive loss.
This preserves rare-profile penalties without necessarily evaluating every
BEST count. The separate FLOAT control below measures its potential
usefulness; it is not substituted into the directed band run.

### Shared-table directed band

For508<=r<=516, seven low-full-degree tables are unchanged. Enlarge only
the high exact rectangle from(F,S)=(256,78) to(258,82). It contains every
requested high target(f,r-2f) with f>=217. Extract all marginal matrices
while each table is alive, then discard the polynomial table. Recompute
each size's exact source factor, route denominator, and indexed profile
counts. Keep the source size fixed throughout every state trajectory.

The new driver checks every odd and even source size separately, including
the source-factor boundary at511/512. If512 is included, its complete
bound and selected f223 bound must overlap the frozen directed receipts.
The final claim requires the sum across every requested size, not merely
successful endpoints. All three new sources passed independent source review.

| Frozen round23 sources | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_neighbor_screen.py` |15706| `c7cce6b8763dd8145e525061c82a9c2ad51139f2d0d66d8c870de476bb374870` |
| `cheap_lane_cauchy_neighbors_directed.py` |14119| `ea39279ceae78310e832f56db5d62a73e7abe08551113bae40e97ac3374f2605` |
| `cheap_lane_profile_peeling.py` |10937| `f0de323a986f98c2d1928cfd2f6606b289f87fc39f5e3cbf993615e34f8441dd` |

### Directed band result and its limit

The complete256-bit run checks every source size508 through516. It uses
2,311 separate marginal normalizations, each size's exact source factor,
and all 593,424 ordered profile counts. Every per-run subtotal agrees with the
indexed cyclic-support count; every size's total is C(8192,r). The512
complete bound and its selected f223 control both reproduce the frozen
receipts. No nonnegative interval intersection is needed.

| Active source groups r | Directed log2 failure upper bound, rounded |
|---:|---:|
|508| -79.0518801540 |
|509| -75.0548010956 |
|510| -70.6806919865 |
|511| -65.9073910975 |
|512| -60.8663747948 |
|513| -55.4658079472 |
|514| -49.7064712431 |
|515| -43.5889654613 |
|516| -37.1138857521 |

The union over all nine sizes has log2 upper-bound interval centered at
-37.097523492662400346062365030503773683131883407931752,
with radius below1.6e-71. This is a valid35% bound but misses the40-bit
working target. A separate256-bit directed sum of the first eight
reported intervals gives

    log2 U_[508,515] = -43.56794746769952382212351260544255371741935104240845
    (interval radius below1.0e-70).

Thus every support of size508--515 is covered jointly with more than43
bits. Size516 still needs a stronger witness for the40-bit target.
The table is not a monotonicity claim and does not certify uncomputed
sizes. The whole-code guarantee remains open.

Reproduce the band with
`C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_cauchy_neighbors_directed.py --rlo 508 --rhi 516 --seconds 1200 --precision 256`.
The successful run takes834.60 seconds and peaks at280,961,024 bytes.
It writes no files. Its witness digest is
`a6d00c43445cc7436a87691e33bb11e8deb10780cdb8fd71476ff2c41677faae`.
The subtotal uses exp(log2_bound*log(2)) for each reported Arb interval,
then sums and takes log2; including the ninth interval reproduces the
complete band receipt.

The peeling diagnostic instruments the existing complete FLOAT control,
then simulates on-demand count queries from the cached chains and counts.
It reports hypothetical count savings, not a measured faster verifier.
Its stopping allowance is2^-80 per source size after the source factor;
allowances across different sizes must be added. Its source passed an
independent review. The complete original control report is unchanged.

At size512, the control performs all33,153 unique count queries. The
simulated peeling needs only1,263 unique queries (3.80961%) and peels
2,490 ordered profiles. Reverse orientations reuse their exact counts,
not their terminal chains. The FLOAT bound changes from
-59.24900641551881 to-59.24900640651474 in log2 scale, a loss of about
9.0e-9 bits. The summed post-source residual tail is2^-85.1759037445,
below the prescribed2^-80 allowance. Twenty exact prefix-bound tests
and fourteen unique-query checks pass, including zero counts, tied
bounds, reordered loose bounds, and asymmetric terminal chains.

The complete diagnostic takes49.08 seconds. It performs no extra count
queries or finite-chain calls beyond the control and writes no files.
This is evidence for a directed on-demand implementation, not a timing
measurement of that implementation or a new certificate. Count savings
alone do not reduce the cost of constructing coefficient tables.

The next strategic test should cover broad ranges rather than continue
only with neighboring-size certificates. The rank interface now gives
an independently audited positive-vector bound for whole source-size
intervals, using the current local caps and exact total occupancy256r.
It needs no profile counts. A separate calibrated reconnaissance at
1024,2048,4096 should distinguish marker mismatch, residue aliasing,
and count relaxation before treating a positive bound as a local-cone
obstruction. The old uncalibrated Cauchy failures do not establish one.

## Twenty-fourth round: retain routing in the middle range

The previous round made progress: it certified the508--515 source-size
band jointly and measured potential count-query savings. This round
tests broader ranges, without changing the encoder or its refresh law.

An independent proof audit identifies a structural limit of the proposed
route-independent positive-vector bound. Empty and full occupancies alone
force its per-size majorant above one for every512<=r<=7525, for all
admissible z,alpha,t,w. The rank interface gives the algebraic proof.
The obstruction applies even at size512, whose route-averaged directed
bound already has60.8664 bits. It is therefore a loss from discarding
routing information, not evidence of inadequate code distance.

A separate conditional construction explains why uniform-over-route
reasoning cannot cover the low-middle range. An allowed route can place
all occupied pairs of a contiguous2048-group set into44-packet windows.
Imposing zero state at each window end leaves a nonzero source subspace.
Some message then uses2016--2048 groups and has relative output weight
at most34.375%, for every realization of the remaining setup. This is a
rare-route obstruction, not a claim that random setup fails appreciably.

The stronger obstruction uses the exact full-occupancy live cap
U_full=P g^64, P>=1. Its entropy endpoint at381/1000 is positive by the
exact integer comparison255^381*1000^1000>2^4000*381^381*619^619.
The main agent checked this comparison in0.0016 seconds. A256-bit
directed evaluation gives0.004593889756726406 bits per coordinate at
that endpoint. The witness digest, hashing the exact numerator and
denominator separated by a slash, is
`3a566e5e5eea32655d8daf04f7147dcbd962fed187acbb551d36072afd3da023`.
The cutoff7525 follows by exact integer division, not a fitted crossing.

Before this stronger obstruction was derived, a FLOAT high-range test
used z=3/1250,alpha=1 and optimized the two vector parameters. The
4301--6143 band bound was+2752731.8543 in log2 scale. Its active update
constraints were empty occupancy, one allowed coordinate, and full
occupancy with live state. The two-band run took3.375 seconds, peaked
at47,480,832 bytes, and passed391 small exact checks. It generated no
files. Its earlier512--4300 analytic exclusion remains a valid subset
of the stronger result; its frozen source is not modified.

| Frozen round24 sources | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_route_free_screen.py` |17378| `d963313ebeb5b64c315b21d71b0555e4e790313c757278dfab2ece34b977780d` |
| `cheap_lane_broad_peeling_screen.py` |24855| `df955f7f8458b0e22314c82311b2b713c2142cd9c24b1d188441edf1a208c33f` |
| `cheap_lane_broad_profile_calibration.py` |23907| `e9f90bdda919247df80e71e5d2555f24ef7ca12323d95b1df622b22556327c56` |
| `cheap_lane_singleton_endpoint_directed.py` |11984| `88f5156ba031f348b23d2e824fd37886822ac1b63874634063272014d2d5577a` |
| `cheap_lane_return_decomposition.py` |16329| `93eaa89b8bfe36bb1b54e0c3e640b6c1db852d64227fa6242510046cc262387c` |
| `cheap_lane_exact_small_return.py` |15028| `62eff6f2f831647696f451fdc3e736e63de0484298781725b2ac3fe89f132a3d` |

The next numerical tests therefore retain exact paired-route averaging
for the middle range. Only7526--7551 remains potentially useful for
route-independent extension below the existing dense certificate.
That narrow extension is not the priority. Neither a failed numerical
witness nor the conditional route example is a whole-code failure claim.

### On-demand counting and the first broader diagnostic

The new FLOAT helper streams ordered profiles one run class at a time.
It evaluates their finite noncommuting matrix chains in batches, then
queries exact profile counts only as needed. It retains each run's exact
remaining count. If the count budget expires, the largest unpeeled chain
bounds the entire residual mass; no source profiles disappear.

At size512, the new implementation reproduces the preceding diagnostic:
log2 bound-59.24900640653069, with1,263 actual unique count queries and
post-source residual tail2^-85.1759. All66,048 ordered chains are included.
It takes3.3004 seconds, versus49.08 seconds for the earlier full-count
control. This is measured proof-screen cost, not an encoding improvement
or a directed certificate. The process peak is52,715,520 bytes.
Exact small tests cover prefix bounds, run-count recurrences, noncommuting
terminals, source counts, FFT coefficients, and scale cancellation.

The first size1024 screen uses z=(17447/25000)^2, alpha=119/1280,
the scaled inherited marker bank, and residues128 in both degrees.
These are uncalibrated witnesses for the new size. The run includes all
263,168 ordered profiles and finishes in8.6971 seconds, peaking at
53,669,888 bytes. Its10,000-query budget expires;826 run classes retain
residual tails above their2^-90 targets. The complete FLOAT envelope is
positive, with log2 bound+9980.7636.

That number combines distinct losses. Five boundary marginals, f=0..4,
have unresolved selected FFT coefficients. The explicitly requested
ordinary-Cauchy fallback retains them, and the all-singleton profile
(f_0,f_1)=(0,0) dominates the resulting envelope. Separately, the actual
profile(384,500) has a positive exact-count weighted bound+2841.9963;
neither of its marginals uses fallback. Its markers come from centers352
and496 rather than its own target degrees. This contribution is not
explained by exhausted counting alone.

### Separate boundary error, aliasing, and the remaining bound

The selected-profile control holds the encoder, local caps, z, alpha,
and exact indexed source count fixed. It compares inherited markers,
target-centered markers, and fitted markers for(384,500) and(384,384).
Unequal profiles use four marker derivatives of the actual ordered
finite chain; equal profiles share two markers. Both fits converge in
four proposals. Independent source review checks the derivatives,
normalizers, exact counts, and output-only terminal.

At modulus128, fitting lowers the first profile's log2 bound from
+2841.9963 to+1848.3323. At modulus256, all three marker choices give
approximately+1839.5129. The equal profile similarly gives+1609.92435.
A separate modulus512 replay gives+1839.512913925 and+1609.924346284
with fitted markers. Across all three choices, the512-residue results
agree within2e-9 bits. This numerical stability suggests that aliases
are no longer the main loss for these profiles. It is not an exact
coefficient identity or a directed certificate.

The calibration run takes4.0933 seconds and peaks at68,022,272 bytes;
the modulus512 replay takes7.0339 seconds and peaks at118,128,640 bytes.
Both pass548 moment, exact small-case, and frozen FFT checks. The latter
imports the unchanged helper and sets its MODULI tuple to(512,) before
calling main; no source file changes between runs.

The all-singleton boundary is different. Setting the full-pair variable
to zero exactly, the new positive univariate verifier extracts degree1024
from P^128 and P^127 T. It then evaluates the finite256-region chain.
The indexed support count is8192/1024*C(7167,1023), independently checked
against the general counter. Directed256-bit arithmetic gives

    log2 bound = -4318.5936051561484532764338628160283772690173944
                 (interval radius below6e-71).

This closes only that ordered all-singleton class, not all size1024
supports. It confirms that the earlier+9980.7636 envelope was dominated
by a loose numerical fallback. The run takes1.3405 seconds and peaks at
41,807,872 bytes. All183 tiny exact coefficient/count tests and33 local
normalizer checks pass; no nonnegative interval intersections are needed.

Four bounded parameter controls test the two non-boundary profiles:

| z | alpha |(384,500), fitted residue256 log2|(384,384), fitted residue256 log2|
|---|---|---:|---:|
|304397809/625000000|119/1280|+1839.5130|+1609.9243|
|3/5|119/1280|+11921.0285|+12325.7184|
|2/5|119/1280|+7858.6109|+7848.9144|
|304397809/625000000|3/40|+3494.8515|+3742.3675|
|304397809/625000000|9/80|+2333.1346|+1468.1104|

Each control takes less than4.1 seconds and peaks below68 MB. These are
FLOAT witnesses, not an optimization theorem. They do not exclude better
joint parameters. They do show that neither repaired boundary coefficients
nor these nearby parameter choices close the selected clustered classes.
### Return paths explain the selected gap

A controlled decomposition multiplies the already powered live-to-quiet
channel by nu, keeping every other entry, marker, count, and source factor
fixed. A path with k returns receives nu^k. There is no compensating
global factor, so nu<1 is a diagnostic submajorant, not a code bound.
The default nu=1 control reproduces both preceding residue512 results.

For nu=2^-16, ordinary Cauchy extraction gives log2 values-513.18884
and-924.94551 for(384,500) and(384,384). For nu=2^-32, they are-513.34490
and-925.11914. The zero-return component is coefficientwise no larger
than either positive marked sum. Thus these FLOAT upper-bound controls
identify paths containing returns as the source of the positive original
majorant. Pure quiet paths and paths with one activation and no return
do not explain that gap. The marked residue coefficients and separate
pure-quiet coefficient remain numerically unresolved at these fixed
markers; no values are inferred from them. No nu=0 limit is claimed.

The decomposition passes two independent source audits and530 exact-small
or frozen FFT checks. Its run takes7.6362 seconds and peaks at119,160,832
bytes. A parameter guard for extreme user-supplied witnesses would also
need to reject underflow after applying nu; the selected default runs
use representable coefficients. No extreme-witness claim is made here.

The next control genuinely improves the bound rather than suppressing
paths. In the actual return matrix S=I+AC, a single allowed input produces
15 or16 nonzero output coordinates. The old generic RS floor charged12.
For two allowed inputs, a projective-line count gives the exact nonzero
moment for each of120 point pairs. Take the largest moment across those
pairs before the unchanged positive four-block recurrence. This preserves
the actual route dependencies and changes only raw E, before alpha.
The rank interface records the proof, independently reviewed.

The exact field checks verify255 inverses, CA=0, S^2=I, rank16, all120
two-column C ranks, and1,285 independently enumerated projective
directions. Eighteen rational moment checks and480 frozen FFT tests pass.
There are four distinct two-input weight histograms; pair(3,14) maximizes
the new moment at the selected z. The row-major16-by16 matrix digest is
`fa88801ea2159d27b7975aa547f9b9cf6b1343fa887b18cc7d24a24153675b02`.

At the same source counts and own-profile markers, the FLOAT residue512
replay gives:

| Ordered full-pair profile | Original log2 bound | Refined log2 bound | Improvement |
|---|---:|---:|---:|
|(384,500)|+1839.512913925|+1105.679800520|733.833113405 bits|
|(384,384)|+1609.924346284|+812.933632156|796.990714128 bits|

The raw return cap improves at1,785 of4,845 occupancy shapes. Every other
channel is unchanged. The7.8341-second run peaks at120,074,240 bytes.
Independent source review passes. Neither selected profile closes, and
these FLOAT improvements do not extend the whole-code certificate.

The next recommended test pools the internal zero/live state of adjacent
packets before alpha. For a fixed route, use(M_0 M_1)[i,j]^alpha instead
of(M_0^alpha M_1^alpha)[i,j]. This keeps a short activation-and-return
alternative inside one power. The raw fixed cone composes, so this step
is valid and pointwise no worse before averaging routes. A degree-conditioned
two-packet implementation adds its own Jensen loss; compare it with the
old same-degree operator and retain a proved entrywise minimum. Preserve
the exact regional degrees, source counts, and output-only terminal.
This is proof work for the current code, not a proposal to add encoder work.

## Twenty-fifth round: tighten actual returns and pool binary paths

This round retains the construction, actual setup distribution, rate,
35% target, source factor, and finite-chain terminal. It tests two
independent losses diagnosed in round24. Numerical statements below are
FLOAT diagnostics, not additional directed certificates.

First pool the internal zero/live state of two packets before alpha.
Conditional Jensen at a fixed total route degree gives a valid cap, but
using it alone is much worse: selected log2 values exceed29000. Its
entrywise minimum with the old pair cap gains223/225 bits. The underlying
pooling is useful; averaging occupancies before alpha loses too much.

A separable concavity correction avoids that loss. For powered path
terms ab and cd, subtract(2-2^alpha) times a lower bound on min(ab,cd).
For each fixed rho, min(rho*a,c)min(b/rho,d) is such a lower bound.
Compute the minima per occupancy before averaging routes, then convolve.
The coefficientwise maximum over nine fixed rho values remains valid.
The rank interface gives the scalar proof and terminal conventions.

Separately, enumerate the exact nonzero output moments of every three-
input submatrix of the actual return matrix S=I+AC. Its proper row flats
give positive integer weight histograms. All560 triples pass rank,
closure, and zero-incidence checks. Direct enumeration of65793 projective
directions independently verifies one triple. There are15 distinct
histograms. Their full indexed-list digest is
`5b3eba025df0a794f016bc495a2660edae47fca2363026db1ca1a84a447bdd7a`.
The uniform maximum sharpens D3, retaining the earlier D1/D2 caps.
Only raw return E changes; it is then powered once.

All comparisons use r=1024, z=304397809/625000000, alpha=119/1280,
residue modulus512, each profile's own markers, and its exact source count.

| Bound | Profile(384,500), log2 | Profile(384,384), log2 |
|---|---:|---:|
| Round24 exact D1/D2 |+1105.679801|+812.933632|
| D1/D2 and degree-Jensen minimum |+883.142688|+587.638732|
| D1/D2 and separable correction |-80.755091|-423.995065|
| D1/D2 and correction/Jensen minimum |-101.005322|-443.496628|
| Exact D1/D2/D3, no pair improvement |+30.711613|-288.298115|
| D1/D2/D3 and degree-Jensen minimum |-143.663247|-464.029499|
| D1/D2/D3 and separable correction |-1053.388984|-1417.227693|
| D1/D2/D3 and correction/Jensen minimum |-1068.237077|-1431.422133|

Thus two genuine refinements remove the positive bounds on these
selected profiles without increasing encoding work. They do not establish
a complete source-size bound. Independent proof and full-source audits
pass for each helper. The composed replay also matches the independent
single-packet engine before applying its pair improvement.

| Frozen round25 sources | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_binary_pair_pooling.py` |18504|`6286943ea9ccfb6acc1f0470117e2510ad2625c80c0bd51cef57482f97183c09`|
| `cheap_lane_exact_triple_return.py` |15341|`d03f60009fcc9f0f757a19681e9647609aad7f501493b1d89980fca8bf1da6b5`|
| `cheap_lane_pair_concavity.py` |18344|`804eeb5b9e7df90cba0ffa57651b0e411784e5bb775e6e3d62b532e994728cf1`|
| `cheap_lane_pair_concavity_triple.py` |10210|`57c3abe184ce6ee643779c48106b79f1d39c85bfb4175e28404ec86b88bc7169`|
| `cheap_lane_pair_bank_screen.py` |17736|`2b5099dc1794a97a2c72c4d09584a62ee248b4982b9889bc791e01e5bd685742`|
| `cheap_lane_old_boundary_regions.py` |11470|`fdb3094bda1f319827c7f50d63cd68720c460cc66498f54575e81c5e2e518698`|
| `cheap_lane_pair_bank_boundary_screen.py` |10741|`15608a5d808df9a180ef7c74fb6c741c69d58c4d998d80b26d2fc60791098927`|

### The complete screen isolates seven numerical boundary entries

The pair-bank helper sums all263168 ordered r1024 profiles with the
frozen exact counting/tail algorithm. Every residual source count remains
in the bound. The first run adds endpoints0/512 to the fifteen-center
bank. Those extreme markers underflow positive local terms, forcing
ordinary fallbacks at f=0..16,511,512. The resulting log2+84089.383 bound
is dominated by those explicit fallbacks and is uninformative. It takes
57.137 seconds with a124854272-byte process peak. The query budget expires,
but all residual counts remain included.

Removing the two extreme centers restores the original fifteen-center
witness, without retuning z, alpha, or the encoder. The complete scan
takes37.410 seconds and peaks at124297216 bytes. It queries2320 unique
exact profile counts; all remaining source-count tails sum to at most
2^-91.230753 in FLOAT arithmetic. No resolved-marginal profile requires
an exact-count query. The log2+6565.993 total comes from the seven
unresolved boundary marginals f=0..6, which use explicit ordinary Cauchy
fallbacks. Thus the resolved interior profiles no longer explain the gap.
This is not a distance certificate or a codeword counterexample.

The boundary helper retains full degree exactly through6 and singleton
degree modulo512, using the older, looser cone in directed arithmetic.
Its6650-slot packed polynomials preserve both the actual terminal and
each target's own denominator. It passes215 exact toy checks,210 local
degree normalizers, six DC controls, and six comparisons with an exact
univariate singleton calculation. No nonnegative interval intersection
is needed. Runtime is7.419 seconds; process peak is81006592 bytes.
The exact f0 conditional log2 control is-176883.1960413302815, agreeing
with the earlier all-singleton calculation after source/count factors.

### Complete size1024 FLOAT cover

The mixed-region wrapper replaces exactly f0..6 with those directed
old-cone bounds. The other506 marginals use the refined/corrected pair
operator. Each replacement retains its own markers and denominator;
source counting, alpha, and the finite-chain terminal are unchanged.
Independent full-source reviews pass for the bank, boundary, and wrapper.

The resulting complete FLOAT log2 bound is **-688.340391794714** at35%
distance. It covers all263168 ordered profiles and every state path for
source size r=1024, at the fixed K=2^20. There are no unresolved marginal
fallbacks. Every run's residual count envelope already clears the tail
target, so no exact profile-count query is needed. The coarser sum of
each exact run mass times its largest conditional chain gives the same
bound. Its leading run is41, with profile(512,471) and log2 envelope
-691.917263211223. This is a worst run envelope, not the exact contribution
of that profile. The wrapper takes26.394 seconds and peaks at124280832
bytes. It generates no files.

This is progress from selected profiles to a complete source-size
calculation, but it is **not a directed certificate**. The local pair
coefficients, cyclic interior extraction, and final chain/count join
still use FLOAT arithmetic. No certified coverage interval is enlarged
by this diagnostic, and no whole-code failure bound is asserted.

Next verify this complete r1024 sum in directed arithmetic. The run-mass
envelope avoids costly exact profile counts. First determine whether a
smaller positive residue modulus retains enough slack, then reconstruct
the actual local correction with outward arithmetic. In particular,
subtract a lower correction from an upper old-pair bound; never treat
floating subtraction as a certificate. Extend to other source sizes
only after this checkpoint, without assuming monotonicity.

Independent review identifies a compact directed implementation. Pack
packet degrees at f+65s: all paired full degrees are at most64, so the
products are carry-safe. Twelve polynomial products build the old pair,
twelve build the raw Jensen pair, and54 build the nine rho corrections.
Use one coherent outward raw-cap table throughout. The directed final
sum can use exact run totals times the largest ordered-chain upper bound;
no individual profile-count queries are necessary for this witness.

There is also an independently reviewed local refinement for ell>=4.
Surjectivity of C_I makes its syndrome uniform in GF256^4. Discard the
inside-I weights and use the exact positive moment of the outside
punctured RS word. This gives a uniform cap without larger-flat enumeration;
intersect it with the old D_ell. It is recorded in the rank interface but
is not included in the preceding numerical table.

## Twenty-sixth round: directed reconstruction of the paired bound

Keep the construction, K=2^20, distance35%, r=1024, z=304397809/625000000,
alpha=119/1280, and the fifteen regional marker centers fixed. Only the
paired interior residue modulus changes in the first control; the old
directed f0..6 boundary still uses singleton modulus512.

| Paired modulus | Complete FLOAT log2 bound | Time | Peak process bytes |
|---|---:|---:|---:|
|256|-688.2131891740006|16.79 s|91111424|
|128|+1371.8674810925527|17.08 s|92401664|

Both controls retain all263168 ordered profiles and have no unresolved
fallback after boundary replacement. Modulus256 needs no exact profile
queries: exact run totals times the maximum conditional chain suffice.
Its leading run envelope is c=41 at(512,471), log2-691.790152049585.
Modulus128 loses too much to aliasing; even2933 exact count queries do not
repair the positive(384,384) contribution. Use256 for directed reconstruction.
Neither FLOAT control enlarges certified coverage.

The new local helper freezes one coherent raw-cap table to exact rational
upper endpoints before powering it. It forms occupancy-level minima before
route averaging, uses78 carry-safe polynomial products, subtracts only a
rational lower correction from an upper old-pair coefficient, and then
takes the minimum with an independent upper Jensen bound. It exports an
exact rational majorant, not an enclosure of the optimal pair coefficient.
Two independent full-source audits pass. The256-bit replay takes9.432 s,
peaks at67518464 bytes, and passes526 local exact tests and561 packet
normalizers, with zero nonnegative intersections. It improves12856 entries;
the Jensen minimum further improves3638. The raw-cap digest is
`5ceb40b034985ff6dac14dbf4efde0444544c30c2cf7934b42105c48035fe0ec`;
the exported pair-table digest is
`745ac297c928162a4098391ee688525b02ade32631d79f192169d4147a53f457`.

The separate directed join takes normalized regional matrices and terminal
vectors. It shares a matrix power between the two ordered orientations,
but retains the correct terminal for each. For each run count it bounds
all degree-feasible chains by their maximum and multiplies by the exact
indexed cyclic-support count. It never drops unqueried profiles. Its36
Fraction checks and840 exhaustive small masks pass; an additional Arb
join contains the independent exact Fraction envelope on an eight-profile
toy. The common source factor remains the caller's responsibility.

| New source | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_pair_residue_control.py` |9318|`b2beb358bdd3d69b14718fec30b9b8f343ede5ea5134f2855f7b8321aa4dd801`|
| `cheap_lane_pair_local_directed.py` |16511|`7284d1cd125d2cf1d1a3c1a38aa3dfd8b2ef7ceae4386ea2975406ea3e64f81f`|
| `cheap_lane_directed_run_mass.py` |12360|`54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb`|
| `cheap_lane_pair_regions_directed.py` |14717|`d3a6e858437e19bede76e215a9d9b31f26a49b9b601ecd60705f1a2b99f0f213`|
| `cheap_lane_pair_regions_resources.py` |4552|`96cece110ff133891793cab23822ac42785465b9960a89ce66914301cb1834a5`|
| `cheap_lane_integer_regions.py` |19792|`9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf`|
| `cheap_lane_pair_regions_integer.py` |11610|`a2a2a6d784b5ae23f26b93e2c4e3caec87f7667598f8bd76d70ab3b49614e9a0`|
| `cheap_lane_integer_precision_control.py` |7130|`4a60a21d3b0a5d1b1fa9b5016aa357e1d9f6493c6224227fe3583d5ad662337d`|

At this intermediate stage, regional reconstruction and the complete join
had not yet run. The subsequent complete result is recorded below.

The first directed center224 pilot hits the512-MiB memory guard during
its third doubling. It reports the first two completed doublings at16.61
and36.76 seconds, with a301187072-byte peak at the latter checkpoint.
It produces no regional bound or certificate. The256-by256 cyclic table
has130561 packed slots per entry and products up to261121 slots, four
times the previous128-by128 table. Before another pilot, recheck physical
memory and use a separately recorded bounded resource launcher; do not
change the frozen mathematics or lower the20-GiB disk guard.

The audited regional driver uses eight tables: three exact-full/cyclic-single
tables at centers32,96,160; four cyclic256-by256 tables at224,288,352,416;
and one exact rectangle F<=512,S<=158 for targets f433..512. Six pair
doublings yield H^64 and H^63 J. The old directed f0..6 exports are lifted
as exact binary rational upper entries, times a directed exponential of
their integer scale minus their exact binary lower log-denominator. No
rounded displayed interval enters the proof. Independent source audits
pass. A separate launcher changes only the explicit memory budget and
requires at least that much free physical RAM plus2 GiB at launch.

The1536-MiB center224 pilot passes its memory guard but reaches the300-s
deadline during the fourth doubling, again with no regional bound. The
third doubling completes at240.3425 s with an803266560-byte reported peak.
Do not run the larger Arb rectangle or full eight-table calculation under
the earlier runtime estimate. The next arithmetic engine represents each
coefficient by a nonnegative integer times one shared dyadic scale per
matrix. Exact integer products and positive folding are followed only
by upward power-of-two rounding. This preserves coefficientwise upper
bounds without an FFT rounding assumption and bounds integer operand size.
The frozen Arb driver remains available as a reference; the construction,
local witness, terminal, and source count are unchanged.

The integer helper passes1181 independent Fraction comparisons across
four axis layouts and six doublings, including ceiling carries, exact
zeros, and terminal order. It takes0.157 s and31.6 MB. Both the helper
and the new driver pass independent full-source audits. The new center224
pilot completes in30.133 s with a212779008-byte peak and all six cyclic
DC upper checks. This is a regional bound, not yet the complete source
class. Every polynomial coefficient is rounded upward at a shared128-bit
mantissa scale; tiny coefficients are not assumed to have128 relative
bits. The explicit majorant, rather than an estimated rounding error,
enters the eventual source-count sum.

### Complete directed size1024 certificate

The exact high-rectangle pilot also passes, in34.770 s with a214544384-byte
peak. The complete integer replay then covers all513 regional marginals,
all263168 ordered matching profiles, every binary/fine state trajectory,
and all indexed supports of size1024. Its directed log2 union bound is

    [-430.01002352988228499310715207931119791608889956951705361777417564922488871
       +/- 9.87e-72].

Thus this entire actual source-size class has more than430 bits of
ideal-setup failure margin at35% distance and K=2^20. The source factor
is applied once; its log2 value is168329.5139862743557384. The exact run
counts sum to binom(8192,1024). The join evaluates131840 unordered powers
but retains both terminal orientations, and makes no profile-count
queries. The largest run envelope is c=41 with log2-433.5872863370173.
The replay takes154.506 s, peaks at230801408 bytes, and has zero
nonnegative intersections. It generates no files.

This rigorously closes a source class that was previously FLOAT-only.
It does not establish a whole-code distance bound. In particular, no
interpolation between certified sizes is assumed. The directed result
is about258 bits weaker than the earlier FLOAT sum, despite matching
the source factor. A fixed128-bit common coefficient scale need not
resolve a very small individual target coefficient accurately.

The independently reviewed precision control reuses one coherent local
table, the center472 exact rectangle, and both terminal orientations.
It changes only the shared integer mantissa from128 to192 bits. The
normalized f471 entries change by less than3e-31 bits. At f512,s0, the
regular return entry falls by3.7576757 bits; its repeated use improves
the two ordered chain bounds by258.3120972 and258.3299769 bits.
Charging each selected chain the full run41 mass gives log2 bounds
-691.9272055 and-691.9172632, compared with-433.6151083 and-433.5872863
at128 bits. Thus tail quantization explains the leading-profile loss.
These selected-chain diagnostics are not a run cover or an improved
whole-size certificate. The retained430-bit certificate includes every
upward-rounded coefficient and does not depend on this explanation.
The control takes66.205 s, peaks at287518720 bytes, and generates no files.

Next generalize the regional extraction and join to actual source sizes,
then test complete-size witnesses at2048,4096,6144 before filling intervals.
Initial screening marks are respectively z=1/4,1/16,1/64 with alpha=119/1280;
they are candidate witnesses, not claimed certificates. Above4096, empty/full
complementation can reduce the extracted degree: for h=8192-r,
f_h=4096-r+f and the singleton degree is unchanged. Retain the source factor
for actual r, never the hole count h. Audit that transformation before use.
Distinguish projection loss, run-envelope loss, and positive local bounds
before changing the construction. If larger occupancies dominate, test the
already proved outside-coordinate return cap for ell>=4. Do not interpolate
between closed source sizes; an interval requires an explicit joint bound.

## Twenty-seventh round: broader sizes and the return constraint

The preceding round is progress: it produced the complete directed
size1024 certificate. This round tests actual sizes2048,4096,6144 without
changing the encoder. The frozen paired caps feed all ordered profiles
and exact indexed run totals. The6144 calculation uses the independently
audited hole-degree reindexing, while retaining the actual6144 source
factor. All outputs below are FLOAT diagnostics, not new certificates.

The initial witnesses are alpha=119/1280 and z=1/4,1/16,1/64 respectively.
Each uses the scaled fifteen-center marker bank and residue modulus256.
The complete screens cover1050624,4198400,1050624 ordered profiles,
respectively, and pass5454 exact complement checks on1344 small cyclic
words. No profile count is discarded. Unresolved residue coefficients
receive a reported ordinary Cauchy bound, not a clipped FFT value.

The new inside-coordinate lemma bounds each affine syndrome coset by
the homogeneous kernel moment. It gives the uniform return cap
D_ell<=q^4 Q_ell B_(16-ell) for ell>=4, before local powering. Only raw E
changes; the pair correction is rebuilt from matching raw/powered maps.
Two independent source reviews pass. The parent replay passes704 exact
GF4 coset checks,44 moment checks, and54 GF256 rational checks. All twelve
block caps D4..D15 improve at the three witnesses;4775 occupancy shapes
improve at z=1/4.

| Actual size | Original complete log2 envelope | Coset complete log2 envelope | Original / coset resolved-profile envelope |
|---|---:|---:|---:|
|2048|+28070.2431|+21026.6017|+13748.0132 / +6497.5339|
|4096|+101366.5214|+82048.8904|+32312.1142 / +11168.9897|
|6144|+200918.9397|+202972.6574|+196567.0305 / +193957.3019|

Each resolved/fallback category separately charges its selected profiles
the full run mass; these diagnostic categories do not partition source
counts. The6144 complete envelope increases because the numerical
fallback set changes, not because the mathematical return cap worsens.
There are55/365/122 original fallback marginals and50/324/98 after the
coset refinement. The six complete screens take9.96--15.37 seconds each
and peak below76 MB. Positive upper bounds are not low-distance examples.

### Separate marker, coefficient, and local-cap losses

Fix actualr2048, profile(768,768), z=1/4, alpha=119/1280, and the exact
run512 mass for every comparison. The positive ordinary objective is
e0 H^16383 J divided by its regional normalizer to power256. It retains
32767 updates and the output-only terminal. A bounded twelve-proposal
fit selects positive rational markers; independent regional replays use
both residue256 and512. These are selected-profile controls, not run covers.

| Return cap | Inherited marker, residue256 | Own-profile marker, residue512 | Fitted marker, residue512 |
|---|---:|---:|---:|
|D1/D2/D3 refinement|+13748.0132|+9224.3277|+9224.3277|
|Inside-coordinate coset|+6497.5339|+2192.7058|+2192.7058|

The fitted coset markers are183307/436492 and90021/391109. At those
markers, residues256 and512 differ by only0.000017 bits. Thus ordinary
marker fitting repairs conditioning and aliases, but does not remove the
remaining positive bound. The standalone controls take7.47 and7.79 seconds
and peak below125 MB. Their finite-chain objective and terminal receive
independent source review.

Five explicit nearby parameter checks, each with a fresh bounded marker
fit, retain the same actual profile and run count:

| z | alpha | Fitted residue512 log2 bound |
|---|---|---:|
|1/5|119/1280|+3359.2389|
|5/16|119/1280|+11517.3673|
|1/4|9/128|+4939.5753|
|1/4|15/128|+2741.8263|
|23/100|13/128|+1353.3423|

These checks show useful parameter slack but do not establish an optimal
witness or a structural impossibility. The rank-interface note proves
rank(I+C_T A_T)>=2 from the implemented tensor-row identities. Retaining
this constraint gives a positive zero-subset bound with residual dimension
at most min(4-v,k+2,n+k-v). Its last term also uses invertibility of S.
The exact tests pass256 GF8 subset-rank checks,4096 rank-one tensors,
1280 zero-subset counts,20 moment checks,4096 direct input vectors, and
39 GF256 rational identities. Independent mathematical and source audits pass.

At z=1/4,alpha=119/1280, the new cap lowers the fitted selected bound
from+2192.70577625 to+1385.58340040. At z=23/100,alpha=13/128 it lowers
the bound from+1353.34230765 to+419.24126755. Fitted markers for the latter
are320890/682431 and251362/1026933. Residues256 and512 differ by less
than0.001 bits there, so increasing the residue modulus is not the next
useful test for that profile. These runs take8.81 and9.58 seconds, with
process peaks below123 MB. They change only the nonquiet return cap,
then rebuild the pair correction from coherent raw/powered tables.

### Retain the first packet in the pair correction

For powered left-row entries a,c and right-column entries b,d, the identity
min(rho*a,c)min(b/rho,d)=min(ab,cd) holds at rho=c/a when a,c>0. The old
bound chooses one rho after route averaging. The new bound instead
partitions first-packet occupancies by their powered ratio. Each bucket
uses a fixed rho and the full right-packet polynomial. Summing disjoint
bucket products preserves every route pair. Independent proof and source
audits pass; the rank-interface note gives the coefficientwise argument.

One matched run keeps z=23/100,alpha=13/128, the same actual profile,
source/run mass, terminal, and the markers fitted on the old9 bound.
The dense grid uses powers2^-16 through2^16. There are no clipped positive
ratios at this witness. Bucket corrections are summed, global choices
are maximized, and valid lower corrections are intersected before
subtraction. The degree-Jensen ceiling remains unchanged.

| Pair correction | Residue256 log2 bound | Residue512 log2 bound |
|---|---:|---:|
|Old nine global choices|+419.2419040|+419.2412675|
|33 global choices|+48.1314123|+48.1308099|
|33 first-occupancy buckets|−337.4088362|−337.4095671|

The grid recovers371.1105 bits and retaining the first occupancy recovers
another385.5404 bits. This identifies both losses without a parameter or
marker refit. All397 exact fraction/bucket controls,45 positive-moment
controls,480 FFT toys, and residual-rank tests pass. The run takes14.756 s
and peaks at126640128 bytes. It generates no files.

This is a negative selected-profile FLOAT bound, not a directed certificate
or a cover of all profiles in the run. No new source-support class closes.

### Complete size2048 control and the next obstruction

The same audited operator now feeds a complete size2048 screen at
residue512. The original fifteen-center bank gains center768, using its
fixed fitted markers above. Every target retains its own normalization.
All1050624 ordered profiles and exact indexed run masses are retained.
The selected anchor reproduces−337.409567132, but the complete FLOAT
envelope is+18329.882342445. It is dominated by profile(1024,0), run1024,
using an ordinary fallback;364 of1025 regional marginals need fallback.
There are436921 wholly resolved ordered profiles and613703 profiles
with at least one fallback marginal. These counts partition profiles,
not the source-support mass assigned to the diagnostic categories.

The separate resolved-profile envelope is+546.177078204. Its largest
term is profile(898,898), run252, at+545.798585659. Thus removing fallback
alone does not establish closure. The fixed bank can still retain large
positive residue aliases even for a numerically resolved coefficient.
This complete control takes26.709 s, peaks at125886464 bytes, passes
the exact route/pair/rank/FFT tests, and receives independent source review.

Next calibrate the actual(898,898) profile at unchanged z/alpha and
compare residues256/512. Then improve the regional marker bank and use
exact or directed extraction for stubborn boundaries. Do not yet spend
a full directed run on the current positive complete envelope. For6144,
first fit hole markers at its current z/alpha. The quiet channel suggests
counterweights q^(-2alpha),q^(-alpha) in the hole marks; this is a proposal,
not a new cap. A quiet/full-live proxy suggests z=1/85 as a later controlled
witness, not as a proven optimum. None of these diagnostics justifies
changing the construction. The whole-code goal remains active.

| New source | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_pair_sizes_screen.py` |20058|`7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd`|
| `cheap_lane_coset_return_refinement.py` |14865|`8fd5db13b0aca7ebb4149772827985dfba322a5ca02e17ea6fd2931a8199bd06`|
| `cheap_lane_coset_sizes_screen.py` |2992|`d031ac8ce17321bc86763c0803274d743b33c4b299c6549af4c35517f8b12321`|
| `cheap_lane_pair_profile_calibration.py` |15175|`5b5e18d2b62b8dbbf5614b5b07c258b4aa1fe6e53a554513518906c931212518`|
| `cheap_lane_coset_profile_control.py` |1753|`8d21a248b69d425d324494c2754ea9c9d0ff001def3519f88828807dd5a4a4f6`|
| `cheap_lane_pair_parameter_control.py` |2593|`b51560e95111a704abc09f4e87e0a2fb75b9509e687065dcf9fb3d447b492d84`|
| `cheap_lane_residual_rank_return.py` |14786|`3efa28f205a01cafe47a8058ca4d3a69e065d5e22e75ea00a8fe79ed4a5dcf7c`|
| `cheap_lane_residual_rank_control.py` |3225|`a50ad947f722ed6650f0425668fe9ce56d410a17bc1e5d34ee4af00387838f8b`|
| `cheap_lane_rank_parameter_control.py` |2412|`50509f0ca5758e0a59c7cd01cca9d8eb1a6ebb9d547df589e9bbbf65180a09fc`|
| `cheap_lane_rho_bucket_control.py` |15899|`57391fb3ab3ada8078c17a98eb921dbb5777452422c429960e3fa47aa07624b8`|
| `cheap_lane_rho_bucket_sizes.py` |9248|`548cbb77b5d9b01e884552b5fda097ff99cf17d62d8f49d4e94d4a449195e89a`|

## Twenty-eighth round: calibrated coverage and finer pair corrections

This round keeps the encoder, setup distribution, source factors, and
exact indexed run masses fixed. All new numerical results below are
FLOAT diagnostics at35% distance, not directed certificates. Source size
r=2048 refers to active outer groups, not K. The complete calculation
retains all1050624 compatible ordered profiles and both terminal orientations.

### Calibration removes numerical fallbacks, not the remaining proof gap

At z=23/100 and alpha=13/128, the selected profile f0=f1=898 has252
singleton positions per matching and252 runs. Fitting its own positive
rational markers gives(523842/992311,3983/34478). Its ordinary bound is
+3759.137708, whereas residues256 and512 give+545.840044 and+545.798586.
The same residue512 result holds at the inherited and own-curve markers.
Thus the resolved positive bound is not an ordinary-marker artifact.

Two fixed output-mark controls at the same alpha give+644.960646 at
z=11/50 and+791.991999 at z=6/25. At z=23/100, alpha=3/32 gives
+230.280757 and alpha=7/64 gives+1240.370058. Every comparison rebuilds
local caps, pair corrections, source factors, and marker fits coherently.
These are explicit witnesses, not a proof of a global parameter optimum.

A complete31-center calibrated bank eliminates every ordinary fallback.
Its fit centers are

    0,3,11,25,44,69,98,131,169,211,256,304,354,406,458,512,
    566,618,670,720,768,813,855,893,926,955,980,999,1013,1021,1024.

Endpoint fitting uses positive half-slot targets only to choose markers.
All extracted profiles, multinomial denominators, and source counts remain
the actual integers. Each fit has at most twelve proposals. The complete
residue512 bounds are:

| alpha | Complete log2 bound | Leading profile | Leading run contribution |
|---|---:|---|---:|
|13/128|+686.441608327|(957,957),134 runs|+682.179074658|
|3/32|+243.594238447|(915,915),218 runs|+239.005502972|

Both runs resolve all1025 regional marginals with zero fallbacks. They
take34.61 and36.58 seconds respectively. The difference from round27's
+18329.88 total is principally removal of numerical fallback bounds;
the remaining positive values are losses in the resolved proof bound.

Exact selected-profile counting cannot fix the earlier f898 gap alone.
Set u=1024,v=3072,a=126. The positive t=1024 term in the exact BEST
count gives

    N_run/N_profile <= uv/[2uv-a(u+v)]
        * C(2u,2a)/C(u,a)^2 * C(2v,2a)/C(v,a)^2 < 253^2 < 2^16.

Vandermonde's identity and log-concavity bound each binomial ratio by
2a+1=253; the prefactor is less than one. Independent derivations agree.
This rules out recovering the observed546-bit loss solely by counting
that profile exactly. It is not evidence of a low-distance codeword.

### Finer buckets retain more of the same local inequality

The [pair correction](CHEAP_LANE_RANK_INTERFACE.md) now has a129-point
rational grid:2^e times{1,5/4,3/2,7/4}, e=-16,...,15, followed by2^16.
Each left occupancy belongs to exactly one nearest-log bucket. The right
occupancies remain unrestricted. Sum the disjoint corrections and intersect
the resulting upper cap with the frozen coarse cap. This intersection is
needed: changing nearest-log buckets alone need not improve every entry.
Inside the grid, the correction retains at least sqrt(4/5) of the ideal
minimum term, versus1/sqrt(2) for dyadic spacing. Source and proof reviews
pass, including six regular/terminal channels and exact rational toy tests.

At f0=f1=915,z=23/100,alpha=3/32, fit only the coarse table and replay
both caps at the identical markers(274146/534529,38287/395245):

| cap | Residue256 log2 | Residue512 log2 |
|---|---:|---:|
|Coarse33 buckets|+239.015660635|+239.005502972|
|Fine129 buckets, intersected with coarse|+46.750010308|+46.739104212|

This recovers192.266399 bits at residue512 without adding encoding work.
It uses43 occupied buckets, with no clipped ratios;2439 local entries
improve. The18.88-second selected control is not a complete run-class or
source-size bound.

The complete fine-bank replay then fits all31 centers using that finished
fine cap. It gives log2+51.969373808, retaining all1050624 profiles with
zero fallbacks. Its leading run has209 runs and profile(920,919), with
log2+47.402773965. The(768,768) and(898,898) anchor controls give
-471.044951257 and+33.335824647. Runtime is39.20 seconds, with a126509056-byte
process peak. The remaining loss is not explained by missing profiles or
fallback extraction. Independent review checks the coherent cap substitution
before every fit and marginal extraction, as well as the unchanged counts.

One further bounded refinement uses eight rational subdivisions of each
factor-two interval. It retains the completed four-subdivision cap as an
independent upper ceiling. All fits then use the final cap. The complete
bound improves to+13.637053781, with zero fallbacks. The same(920,919)
profile leads at+9.076049285. It uses77 occupied buckets, takes44.22 seconds,
and peaks at124846080 bytes. The smaller38.33-bit improvement motivates
stopping grid densification here. The helper supports sixteen subdivisions,
but that option has not been run.

A stronger independently audited lemma orders both left and right ratios.
On different ratio buckets, the minimum of the two path contributions is
known exactly. Positive prefix/suffix convolutions retain those exact terms;
only the same-bucket rectangle needs the separable approximation. The
[rank interface](CHEAP_LANE_RANK_INTERFACE.md) records this argument,
an exact-rational plan for directed ordering, and an unimplemented
ratio-dependent deficit refinement. Neither local lemma by itself closes
a source-size class.

The matched ordered-bucket control uses f0=f1=920 and208 runs, at the
same z=23/100,alpha=3/32. Fit the four-subdivision table only, giving
markers(488277/948091,59265/641956). The residue512 bound improves from
+47.397495286 to-34.531093633; residue256 gives-34.517924619. This is an
81.928589-bit gain at identical markers and counts. The new correction
uses43 left buckets and126 positive cross-convolutions, with no clipped
ratios. All8571 improved local entries retain the prior upper ceiling.
Runtime is25.60 seconds, with a128339968-byte peak. The exact204-check
toy suite covers six channels, boundary ties, zero paths, and clipping.

The initial prototype stopped in that toy suite before constructing any
real cap: the toy degree(1,1) needed total-degree extent2, not coordinate
extent1. Correcting the extent leaves production degree32 geometry
unchanged. The corrected source passes the toy suite. A complete calibrated
bank using the ordered cap then gives log2-29.970305594 over all1050624
profiles, with zero fallbacks. The leading profile is(920,920), contributing
-34.531093633 when charged its full208-run mass. The anchors(768,768) and
(898,898) give-559.917373333 and-49.541075071. Runtime is45.09 seconds,
with a126672896-byte peak. Thus complete size2048 FLOAT coverage is now
negative, but it is below the40-bit target and still needs directed
reconstruction. No new directed source-size certificate is claimed.

One preselected midpoint control, alpha=23/256 at the same z=23/100,
rebuilds every cap, source factor, and marker witness. Its complete
residue512 FLOAT bound is **-86.746269063**, with all1050624 profiles and
zero fallbacks. The leading profile(895,895),258 runs, contributes
-91.464053769. Runtime is42.71 seconds, with a126189568-byte peak. Stop
parameter search at this witness and pursue directed reconstruction.
The exact positive rational markers chosen by the31 fits are frozen here
for that next step; these are proof witnesses, not encoder parameters:

    0: (255/646883, 623632/483565)
    3: (2419/1024041, 293785/228638)
    11: (8764/1015915, 76375/60142)
    25: (16268/835553, 990191/796033)
    44: (19097/562552, 1188291/982718)
    69: (20840/396259, 488307/419335)
    98: (13149/178492, 2571/2308)
    131: (81591/841483, 237530/224483)
    169: (27031/219901, 995648/999729)
    211: (86273/572830, 942700/1014063)
    256: (14896/83141, 497696/577845)
    304: (145616/698575, 268421/339217)
    354: (145441/611765, 436473/605084)
    406: (98552/369131, 566040/868703)
    458: (292881/992593, 333124/569733)
    512: (212837/658769, 338729/653828)
    566: (246817/705174, 255105/561733)
    618: (122742/327335, 268248/679163)
    670: (352629/883663, 268101/793189)
    720: (106361/252384, 138583/485930)
    768: (461297/1043126, 96565/408826)
    813: (458463/994147, 185747/968869)
    855: (14701/30731, 58540/386667)
    893: (470924/954065, 75529/651815)
    926: (45263/89353, 48961/571044)
    955: (457383/883348, 47125/788193)
    980: (483733/917311, 29909/791045)
    999: (450692/843147, 674/31575)
    1013: (422884/783397, 3992/427051)
    1021: (119281/219747, 649/255264)
    1024: (96022/176593, 432/1020359)

The endpoint half-slots affect only marker selection. Every extraction
still uses its actual integer full/single counts and exact multinomial
normalizer. The35% directed coverage listed at the top is unchanged;
the86.7-bit result is a complete FLOAT bound for size2048 only, not for
all source sizes and not for40% distance.

A matched residue256 replay of the same complete witness gives
log2-86.742499919, only0.003769144 bits weaker than residue512. All31
rational markers are identical, and all1050624 profiles remain resolved
without fallbacks. Runtime is29.70 seconds, with a104185856-byte peak.
This smaller residue grid fits the existing positive-integer engine's
packed-coefficient limit, unlike the512-by512 grid. The next experiment
is therefore a directed local-table reconstruction followed by a bounded
residue256 replay at these fixed witnesses. Start with selected regional
entries and a precision control; do not launch a bulk coefficient tree.

For that reconstruction, freeze powered local caps upward to exact
rationals, rebuild both the base pair sum and ordered correction from
that same table, and use lower bounds for every subtracted correction. Preserve every
profile's own denominator, the actual r2048 source factor, indexed run
counts outside alpha,32767 updates, and the output-only terminal.
The coefficient-source ordering argument is proved in the rank interface.
The optional weighted-deficit source is unrun research and is not part
of this86.7-bit witness; do not transfer results from it. Main and independent
source review find no blocker, but its toy suite and numerical control
have not been executed. Keep it secondary to the directed replay.

Frozen round28 sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_rho_profile_calibration.py` |9539|`97e237880824fec84527a73d400a811006ed7a492aeda774004bee0e54ade633`|
| `cheap_lane_rho_calibrated_bank.py` |10913|`e9ff3adda8aa363383386b8e8665b40da5f1075a8e43ddf8ffc1854ddea60776`|
| `cheap_lane_calibrated_bank_parameters.py` |1987|`f221400c61509b24853e217f5601cbc2d44e39fb7bd2f98a9b1d23280c79ee89`|
| `cheap_lane_rho_fine_buckets.py` |13059|`ca8a73442ff87ca46ce64192dd999967716cf145c65a9712553c0e81a9a69701`|
| `cheap_lane_rho_fine_bank.py` |4973|`a82044f15e0cfdf20dd85bd92818f49d42c5a24be707ad55eedec077c8f83625`|
| `cheap_lane_rho_grid_bank.py` |8601|`508607362da04637ead7fc7cdaac1e445a91bb8ee1e398ee1f0426b88d6e12e2`|
| `cheap_lane_rho_ordered_buckets.py` |15523|`f816addd50969a83ad951567944d905938ece966c9260175e3b4a08b1be53255`|
| `cheap_lane_rho_ordered_bank.py` |4711|`b80c002b02e5b0075b74cb8b34a56c16001f2201561a973f420c09eb7182978a`|

## Twenty-ninth round: directed ordered-ratio reconstruction

Keep the encoder, distance35%, actual source size2048, z23/100,
alpha23/256, and all31 rational marker pairs from round28 fixed.
This round reconstructs that witness; it does not fit new parameters or
change the refresh distribution. The baseline and prior certificates
remain untouched.

### Coherent local arithmetic

The new local constructor rebuilds the actual ell1/2/3, coset, and
residual-rank return bounds before taking the single fractional power.
It freezes each powered cap upward to a multiple of2^-256. The base pair
sum and every ordered-ratio correction use this same exact table.
Integer ratio comparisons, exact route numerators over2^32, and stride65
integer polynomial products avoid ambiguous interval ordering. One LCM
clears the diagonal rho denominators. Only a lower enclosure of2-2^alpha
is subtracted. A separately valid raw-cap Jensen bound is intersected
afterward; its table need not share the powered rounding.

The monotone-majorant argument in the rank interface justifies upward
freezing. Every grid-point rho witness is dominated by the ordered
correction on this same table, so the old9/global33 products need not be
reconstructed. This does not justify combining tables from different
rounding schemes. All six MM/MT channels and the output-only terminal
are retained. Main and three independent source reviews find no blocker.
Exact rational tests cover zeros, clipped ratios, boundary ties, terminal
order, carry-free polynomial packing, and monotonicity under freezing.

At the fixed witness, the constructor uses4845 occupancy classes,
29068 positive powered atoms and2 exact zeros, and617 integer correction
products. It records the actual freezing ratio rather than assuming it
is one. The powered-table hash is
`d16642f585fbbb1c02260831645d417a1773ce441cf732c1a82c9c35f4bbf512`.
Working Arb precision is384 bits; regional arithmetic uses exact integer
products followed by upward192- or256-bit dyadic quantization.

### Selected precision controls

The center893 cell covers full counts875--909. For the leading ordered
profile(895,895), which has258 runs and258 singleton positions in each
matching, charging the entire exact258-run mass gives directed log2 bounds

    192 bits: -91.462075527167962548907665004801705148342182904217...
    256 bits: -91.462075527167962548907665004801705148342182904221...

The control takes72.693 seconds and peaks at333684736 bytes. It agrees
with the earlier FLOAT value but is a bound on this selected profile,
not yet on its entire run class or source size. A separate center1024
control covers full counts1023--1024. All six reported normalized entries
at1023 agree between192 and256 bits beyond the first45 decimal places
in their logarithms. That replay takes68.038 seconds and peaks at
256655360 bytes. No r1024 boundary table is substituted.

### Complete directed size2048 certificate

The complete31-cell replay passes. At35% distance its directed log2
union bound is

    -86.7424999191238645958893203926089944792182381689367683005040598962483828704386481099728502852213938077752075787
      +/- 3.37e-110.

Thus the ideal-setup failure probability for this entire source-size
class is below2^-86.7424. It covers all1025 marginals, all1050624 ordered
profiles, every indexed source support of size2048, and all state paths.
The run join uses525824 unordered matrix powers, retains both terminal
orders, and checks that its exact support counts sum to binom(8192,2048).
Each target removes its own multinomial/marker denominator. The actual
source factor uses h8,c_h=1/9 and
eta=447126155859705527567/425; it is applied once, with counts outside alpha.
There are32767 updates, no ordinary fallback, no historical boundary
substitution, and no nonnegative intersection in the regional join.

The largest run contribution is the258-run class, log2-91.462075527168.
The full result agrees with the earlier FLOAT value-86.74249991909315
to about3.1e-11 bits, but the directed replay is the certificate.
Runtime is848.767 seconds, including18.708 seconds for the final join;
peak process memory is282468352 bytes. The two pilot replays and complete
run produce the same exact rational local-coefficient hash:
`e14e75ef2f3cb5fa1a01803ae27eb75befcc38bc5e28cda1ecf535471385d094`.
The new local tests make796 exact comparisons, with526 inherited local
checks; regional arithmetic makes1181 exact comparisons. The run-count
and noncommuting-terminal tests also pass.

Reproduce serially from the repository root, with at least20GiB free:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordered_regions_directed.py --complete --seconds 1200

This closes one entire source size at fixed K=2^20, not all source sizes
and not40% distance. No encoding operation, randomness distribution, or
baseline proof changes. The whole-code goal remains active.

### Next coverage test

At fixed z, alpha, and markers, the local caps and marked H64/H63J
polynomials are independent of actual source size. Reuse each table
while live to extract actual sizes2047--2049,
using each size's own singleton degree, normalizer, and source factor.
The existing normalized size2048 entries cannot supply those neighboring
coefficients by rescaling. The joint-missing factor changes at2048, so
no continuity or profile monotonicity is assumed.

A newly audited band lemma, recorded in the rank interface, distributes
Omega_r^alpha across the256 regional factors and then takes positive
entrywise envelopes across r. Exact run masses regroup by d=f0+f1=r-c.
This can share a single expensive matrix-power join across the band.
Its numerical loss has not been tested. First require singleton-band
reproduction, exact small examples, and the full band count identity.
If the entrywise envelope is too loose, retain the common source size
within adjacent pairs or use narrower bands.

After this neighboring-band check, test a complete current-operator
size4096 FLOAT witness before changing the construction. Its old positive
screens used weaker local bounds and incomplete calibration. The proved
route-free relaxation obstruction is not an obstruction to this routed
operator; the49.8047% quiet-subcode ceiling does not obstruct35% either.

Frozen round29 sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_ordered_local_directed.py` |20640|`1ed1d59d9beed36b9a1bd1a2bc65af043efe342e95e2afd5047f90c696c04c47`|
| `cheap_lane_ordered_regions_directed.py` |10452|`1faaf32f2cf1ced0e78fb9952534e1c79b138cf9f530c7c82b82607e650dd78d`|

## Thirtieth round: shared neighboring-source tables

Keep the round29 local witness, encoder, and35% target fixed. The new
extractor builds each marked cyclic table once, then extracts actual
sizes2047,2048,2049 with separate singleton indices and exact normalizers.
The source-scaled band join uses the lemma in the rank interface. It
checks every original marginal before accepting complete coverage.

Two independent source reviews and the main review pass. Before running,
the driver is corrected to use the helper's ordered-envelope count key
and to select only a pilot degree with nonzero exact support count.
The2049 neighbor introduces the ordered endpoint(1024,1024); the2048
exclusion is not reused. The complete band therefore needs1050625
ordered envelope pairs. Its counts must sum exactly to
binom(8192,2047)+binom(8192,2048)+binom(8192,2049).

### Selected band control

At center893, all420 lower/upper endpoint comparisons reproduce the
frozen size2048 extraction exactly. The new helper tests pass304 exact
fraction/interface checks and588 exhaustive small-cycle supports;
the extractor passes80 exact checks. Inherited integer and single-run
tests also pass. No whole-size claim is inferred from this pilot.

For the selected pair(895,895), d1790 corresponds to257,258,259 runs
at the three actual source sizes. Charging each size's full corresponding
run mass gives directed log2 bounds

    r2047: -93.823754871433435813651791198984...
    r2048: -91.462075527167962548907665004801...
    r2049: -89.021853721230130683252453893381...

The sum has log2-88.73484856311348. The simple three-size envelope gives
-78.93788810010936, losing9.79696046300412 bits in this control. This
is neither the maximum over all pairs with d1790 nor the complete band.
Runtime is68.687 seconds, peak process memory270168064 bytes, and
generated-data growth zero.

### Complete neighboring-band certificate

The complete31-cell replay passes. Its directed log2 union bound is

    -74.2104337688321901564597776892545863474841165098489863450043038998339566592412144435054269283112147068226834897
      +/- 2.88e-110.

This is a joint bound for every actual support of sizes2047--2049 and
all state trajectories, at35% distance. It retains1050625 ordered
envelope pairs, using525825 unordered powers and both terminal orders.
All per-size marginals are present; the exact count identity passes.
There are no nonnegative intersections in the envelope or join.
The largest degree contribution is d1791, with log2-78.929527237217.
The size2048-only86.7425-bit certificate remains independently available.

Runtime is856.263 seconds, including17.781 seconds for the run-mass
join; peak process memory is279838720 bytes. The exact summed-support
count has SHA-256
`59bcbc971dfe4293b5455cf944081f03e79ea3e1ce04e933a8b1000fed437b21`.
This digest identifies a count, not a separate certificate file.

Reproduce serially from the repository root, with at least20GiB free:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordered_band_directed.py --complete --seconds 1200

This certifies a band of source supports at fixed K=2^20, not a range
of code lengths and not the whole code. No encoder operation changes.

### Count-aware follow-up

The simple envelope can pair one source size's large conditional matrix
with another size's large run count. The newly audited degree-specific
lemma absorbs each exact count as N_r(r-d)^(1/256) before taking maxima.
The old source factor is already distributed as Omega_r^(alpha/256).
The resulting chain is multiplied by the number of eligible actual
sizes at that degree; only z^(-alpha L) remains outside. A marginal is
eligible only when both it and its partner index exist for that size.

Two independent mathematical/source reviews and the main review pass.
The A/B keeps every construction and proof parameter fixed, extracts
actual sizes2040--2056 from one shared table, and compares the same
selected pair against the old bound and the sum of individual bounds.
All62 exact-fraction/interface tests pass, including noncommuting terminal
order, singleton reproduction, count-root scaling, and partner validity.
The three-size old-envelope and individual-sum controls reproduce the
previous pilot exactly. All runs retain the same local-coefficient hash.

For the selected pair(895,895), d1790, the directed log2 bounds are:

| Actual sizes | Sum of individual controls | Source-only envelope | Count-aware envelope |
|---|---:|---:|---:|
|2047--2049|-88.73484856311348|-78.93788810010936|-87.43689122050897|
|2040--2056|-69.15172667473704|+11.10821317536483|-65.24592039769524|

The count-aware loss against the individual sum falls to1.298 bits in
the three-size case and3.906 bits in the17-size case. The latter recovers
76.3541 bits relative to the source-only envelope. It directly addresses
the identified mismatch between conditional bounds and support counts.
This is a selected profile, not the maximum over its degree or complete
source-size coverage. The wider band remains uncertified.

The selected run takes39.871 seconds and peaks at259727360 bytes; it
generates no files. Reproduce the selected controls with

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_counted_band_directed.py

Next run the complete2040--2056 calculation, using the same frozen
markers and separate actual targets. Its new f1025--1028 marginals must
be extracted explicitly. The bounded command is implemented and audited
but has not been run:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_counted_band_directed.py --complete --seconds 1800

After that range check, return to the current-operator size4096 FLOAT
screen before changing the construction. The generic count/routing
helpers support4096, but the current ordered calibrated-bank driver is
hardcoded to2048 and needs new orchestration. At4096,f0,s4096, the local
degree(0,64) is forced in every paired factor; H_(0,64)^64 and
H_(0,64)^63 J_(0,64) give a direct boundary extraction with denominator1.
This avoids an extreme marker fit at that endpoint. This source-level
plan has not been run and adds no certificate claim.

Frozen sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_band_regions.py` |11803|`ab74014cbd1e2e6334ab15acaf2fcbc9dc6d9bc97b19e5e6c495f0d2eb8a5852`|
| `cheap_lane_band_run_mass.py` |15833|`127a6a82fcb785bcb339115f0e0ea8d951841340ce32540655b3e86b06a20503`|
| `cheap_lane_ordered_band_directed.py` |12054|`e775786311a3840d21a708beb529720bf4a0e32a7a064fc69f038641555c3ce1`|
| `cheap_lane_counted_band_run_mass.py` |15369|`2c481652ef8a37a14c5bd88ce6af31e0517a3f51706b7a07e3aa6447906ab74b`|
| `cheap_lane_counted_band_directed.py` |8792|`f150f51ae101ec3f795f688621fe85b135b63b837d26bd3051e8cedd8f39c2f2`|

## Thirty-first round: complete count-aware band and larger supports

The complete2040--2056 replay closes the wider band. It keeps the
round29 local witness z23/100,alpha23/256, all31 frozen marker cells,
and the actual cheap-lane encoder unchanged. Every source size uses
its own coefficient, multinomial normalizer, and source factor. In
particular, the new f1025--1028 marginals are extracted, not interpolated.

The directed joint log2 union bound at35% distance is

    -60.2491325684610942192784539108966504960891355530198614280445127945051851129526387713407082868007083640191454587
      +/- 5.47e-110.

This covers every actual support of sizes2040--2056 and every state
trajectory at fixed K=2^20. It is not a range of code lengths or a
whole-code certificate. The join retains1058840 ordered envelope pairs
and529934 unordered powers, including both terminal orders. All2056
degree classes and original own-size marginals are present. The exact
counts sum to sum_(r=2040)^2056 binom(8192,r); that integer has SHA-256
`5dfcfa88dc271b3d882ac30e202c17392b92e99b4252603d27e5b7642325c440`.

All62 exact-fraction/interface tests pass. The local coefficients retain
SHA-256 `e14e75ef2f3cb5fa1a01803ae27eb75befcc38bc5e28cda1ecf535471385d094`.
Local and input-region nonnegative intersections are zero. The nested
selected controls reproduce round30, including log2-65.2459204 for the
17-size profile(895,895). The complete join uses count roots1/256,
not alpha/256, and applies the output threshold factor only once.

Runtime is981.428 seconds, including83.158 seconds for the join; peak
process memory is293822464 bytes. No data file is generated. Reproduce
serially, with at least20GiB free, using

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_counted_band_directed.py --complete --seconds 1800

The stronger86.7425-bit size2048 certificate and74.2104-bit nested-band
certificate remain available. The new joint band supersedes neither
receipt. It extends coverage to2040--2046 and2050--2056 without changing
the construction. Larger source-size gaps still require proof work.

### Current-operator size4096 test

The new bounded driver uses the same composition, residual-rank return
bound, and ordered129 pair correction as the improved size2048 screen.
It fits31 deterministic regional markers for actual size4096, retaining
each target's own normalizer. The all-singleton endpoint is extracted
directly as H_(0,64)^64 and H_(0,64)^63 J_(0,64). Two independent source
reviews and the main review pass; the boundary argument is recorded in
the rank interface. No frozen helper or encoder is edited.

This FLOAT diagnostic must cover all4198400 ordered profiles and retain
any unresolved coefficient through an explicit ordinary Cauchy fallback.
It reports resolved and fallback leaders separately. Neither category
partitions source counts: each diagnostic charges the full applicable
run mass. Negative FLOAT output is not a directed certificate.

At z1/16,alpha119/1280, the complete screen gives:

| Cyclic modulus | Complete log2 bound, FLOAT | Fallback marginals | Runtime | Peak bytes |
|---|---:|---:|---:|---:|
|256|-4621.155202411889|0|35.201s|105381888|
|512|-4621.176571168393|0|49.513s|126947328|

Both runs cover all4198400 ordered profiles and pass the exact run-count
identity. They use identical rational markers at all31 centers. The
minimum selected coefficient divided by its channel's DC mass exceeds
3.5679e-7; neither run uses the ordinary fallback. The78 exact boundary
checks and inherited local, rank, scale, and finite-chain tests pass.

The leading run has c2048 and ordered profile(2048,0). Its full-run
charged bound changes from log2-4621.6241621 to-4621.6455966. These are
already far below the working threshold despite the conservative count.
The two indexed alternating supports, profile(0,0), have a separate
diagnostic log2 bound of-14200.2952230. That diagnostic is already covered
by the complete sum and is not added again.

The higher-resolution replay improves the complete FLOAT bound by only
0.02137 bits. This is encouraging evidence, not a rounding certificate.
The old positive size4096 diagnostic used weaker local bounds and
different marker choices; it does not identify a remaining obstruction
for the current operator. The next step is directed replay of the
leading endpoint pair, then the complete fixed31-cell witness. Exact
count peeling and construction changes are not needed for this witness
unless directed arithmetic exposes a new loss.

### Directed leading-profile control at size4096

The new directed wrapper freezes all31 successful rational markers and
rebuilds the local coefficients independently. Main-agent comparison
confirms every marker matches the FLOAT receipt exactly. The main and
independent source reviews pass. The default selected run uses384-bit
Arb enclosures,192-bit upward integer coefficients, center2048, and the
independent marker-free f0 boundary. It extracts f2046--2048 with their
own actual4096 normalizers.

The ordered profiles(2048,0) and(0,2048), each charged with the entire
exact c2048 run mass, give directed log2 bounds

    (2048,0): -4621.624162121373180412029395328262035671315968007200284382584836199745374381361724634560107960545073823114037200
                +/- 3.38e-109;
    (0,2048): -4622.3215878249840216052527231990438895525417065709645914507416797678824421922480907458482763307141662064312175329
                +/- 7.88e-110.

The first reproduces the modulus256 FLOAT leader to numerical precision.
The exact boundary profile(0,0), including its two indexed supports,
has directed log2-14200.29522300782. Its regional and direct finite-chain
enclosures overlap. This diagnostic is not added again in the complete
join. These are selected-class bounds, not a complete c2048 run bound or
size4096 certificate. No other profile is inferred from them.

All inherited exact local, boundary, extraction, integer, and finite-chain
tests pass. Nonnegative intersections are zero. The new exact rational
local-coefficient hash is
`afebea4130f5cfa5c6ae6cb5abfc465358d42484f11673bf9f3b68f4ec548cba`.
This differs from the band receipt because z and alpha differ. Runtime
is37.066 seconds, peak process memory205504512 bytes, and generated-data
growth zero. Reproduce the selected control with

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordered_large_directed.py

Next run its audited complete mode, still unrun at this checkpoint:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordered_large_directed.py --complete --seconds 1200

That mode requires all2049 own-size marginals and4198400 ordered chains
before reporting coverage. A successful size4096 replay would still
leave the source-size gaps in the top table. The goal remains active.

Frozen new source:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_ordered_large_screen.py` |16303|`1a0655fde36dbb4f090a143998833bf8ccaf8c25f029df42b1b8837a1db670ed`|
| `cheap_lane_ordered_large_directed.py` |13753|`9b2cef00e449dc10c0cd58eca062caf1986ca3f7eed5003087abe32e53e6b6d6`|

## Thirty-second round: larger-support certificates and cheaper projection

The full frozen size4096 replay closes that entire support-size class,
not only its selected leader. It retains all2049 matching marginals,
4198400 ordered profiles, exact indexed run counts, and all state paths.
At z1/16 and alpha119/1280, its directed bound at35% is

    log2 U = -4621.155202411949657928030014748782785578604616381996458677907590926691020392815363955179151253663604761694887129
                +/- 5.49e-109.

The leader remains run2048 with selected ordered profile(2048,0).
All tests pass, with zero nonnegative intersections and no fallbacks.
The run takes941.572 seconds, including71.878 seconds for the join,
and peaks at283136000 bytes. Its local coefficient hash remains
`afebea4130f5cfa5c6ae6cb5abfc465358d42484f11673bf9f3b68f4ec548cba`.
The exact sum of indexed run counts has SHA-256
`20d556b0911b5fb3862c3eff8132091202f00f6221812c740f949c2e94fd8faf`.
Reproduce with the round31 helper's complete command above.

### Dropping regional projection is not uniformly safe for tightness

The new ordinary-Cauchy control keeps the same local caps, markers,
own-target normalizers, exact f0 endpoint, and complete directed join.
It replaces projected regional coefficients by H(u,v)^64 and H(u,v)^63 J.
This remains a valid upper bound, but its complete log2 bound is
**+262.5035244260433**, dominated by run2047. It proves no useful
guarantee. The earlier selected(2048,0) control still has3518.2974 bits;
that selected success cannot stand in for the complete join. The control
takes94.230 seconds and peaks at101576704 bytes, with all exact tests
passing and no nonnegative intersections. No encoder defect follows
from a positive proof majorant.

The complete profile diagnosis reproduces FLOAT log2+262.503524426142,
with directed checks on its eight leading channels. The actual leader
is(2013,36), with matching singleton counts70 and4024. Its selected
directed bound reproduces the complete run leader to numerical precision.
Other leaders include(2013,14),(2035,36),(2013,69),and(2013,113).
These sit at coarse marker-cell boundaries; the initial suspicion that
f1 caused the loss was incorrect. Runtime is29.559 seconds, peak110915584
bytes, with zero nonnegative intersections. The FLOAT complete join is
not a separate directed certificate; its selected checks are explicitly
restricted to their own ordered channels.

The full/empty coefficient transformation in the rank interface has a
bounded reusable implementation and380 passing exact rational checks.
The tests take0.049 seconds, create no data, and check truncation, terminal
order, and own normalizers at4095--4097. This partial regional helper
adds no size coverage. Its optional edge refinement is not run as a fix
for the ordinary leader, because that leader lies outside its small strip.
The matched test selects entrywise among all31 existing ordinary
marker witnesses instead of assigning only the nearest center. It succeeds:

    log2 U = -256.23339980531610985762432108225060146528414920979472320521922337399384477779302288930614010835863402419091008
                +/- 2.73e-108.

This complete directed bound retains all4198400 ordered profiles, the
exact f0 projection, and the same source factor. No markers, local caps,
or encoder operations are added. It examines381114 candidate entries
and changes only203 of the12288 non-f0 entries relative to nearest-center
assignment. Each chosen entry is a valid majorant, so no marker union
factor is needed. The old leader(2013,36) improves from log2+262.5035
to-577.6637799306597; its reverse gives-578.3338002397580. The new worst
run is1971, with directed log2-256.8047864000295.

The94 exact selection/recurrence controls,409 ordinary controls,
279 exact-denominator comparisons, and inherited boundary/join tests
pass. Nonnegative intersections are zero. Runtime is85.550 seconds,
peak101429248 bytes. Thus an approximately11-times cheaper proof
verification still exceeds the40-bit target for this whole size class.
The earlier4621.1552-bit coefficient-based receipt remains the stronger
bound. Neither result certifies neighboring sizes or the whole code.
Reproduce the cheaper complete verification with

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordinary_marker_envelope.py --seconds 300

### Actual size6144: complement coordinates, unchanged source factor

The new complete FLOAT screen rebuilds the current local operator at
z1/64,alpha119/1280 and fits31 positive rational markers. It reverses
the complete local degree table into hole coordinates once, retaining
the original terminal chronology. Actual size6144 has2048 holes.
Its1050624 ordered profiles give FLOAT log2-53312.93610343699, with
zero fallback marginals. The leading hole profile is(996,0), corresponding
to actual full counts(3044,2048), at run1052. This is reconnaissance,
not a directed certificate. Runtime is30.536 seconds, peak106409984 bytes.
The FFT controls have maximum imaginary/DC1.120e-17, negative/DC4.233e-18,
and minimum selected/DC1.591e-5.

Unlike size4096, the ordinary-Cauchy replay has enough slack here.
It freezes those same31 rational markers and builds directed local caps,
then evaluates ordinary H64/H63J without any coefficient projection or
boundary substitution. Every own hole normalizer is checked against the
corresponding actual multinomial. The join uses coordinate size2048 for
the identical support counts; the source factor is Omega6144, with
h24 and c_h1/25, applied once. The complete directed result is

    log2 U = -49748.334323071804335024723216219253931990546889739730175829676106033058025197900778945264105952652192617784073555
                +/- 8.07e-109.

This covers all actual6144-group supports and all histories at35%.
It does not cover adjacent sizes. The worst run is1062. The409 ordinary
rational checks,5454 complement checks over1344 toy cyclic supports,
36 finite-chain checks, and840 join masks pass. Nonnegative intersections
are zero. The local rational-coefficient hash is
`75ac478004e1572d716752cf989cee3cf1182fdf98e17ad67a401ba5b6c7bab3`.
The complete replay takes32.328 seconds and peaks at102699008 bytes.

Reproduce the two new size6144 calculations serially:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordered_hole_screen.py --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordinary_hole_directed.py --seconds 300

These results add genuine support-size coverage while leaving the
construction, encoder work, and baseline certificate unchanged. The
whole-code goal remains active; gaps are listed in the top table.

### Next coverage step

The rank interface now gives an independently reviewed ordinary marker-cell
band lemma, using exact rational polytope vertices and concave
log-multinomial denominators. It applies in actual or hole coordinates;
the density tilt always uses the actual source size. This lemma is not
yet numerically replayed. The6144 result motivates an initial broad hole
band such as actual sizes5888--6400. That run must cover hole indices
through1152, not stop at the6144 bank's1024 endpoint. Extend the final
marker cell explicitly or add witnesses, retain odd source sizes, and
sum the true actual-size source factors. Do not infer the band from
the singleton's margin. Likewise, an entrywise minimum over marker
witnesses is not automatically reproduced by one marker per band cell.

Frozen new sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_ordered_ordinary_directed.py` |11825|`cb2cc15b9209a95cdfd57b51d98baed99aa205a82099595075c2e72563b826ed`|
| `cheap_lane_ordered_hole_screen.py` |16205|`cb9e6201b3fda3d7598b66339be99ba4584c8c87524cd7400d96089f87320854`|
| `cheap_lane_ordinary_hole_directed.py` |8360|`a430c6c72bdc4c5b19ccf3d99da49fbc0e1385ee6361bc020469f8bc1806b217`|
| `cheap_lane_ordinary_profile_screen.py` |15058|`88005a2b493d50d1bf6f037b803d1e964e52d9e3d45ef97abb8c30b40bebd416`|
| `cheap_lane_edge_regions_directed.py` |11642|`2fbe2e1220f083914aa4419939e84927f7c61573163f6d56250e250c7a5a5c42`|
| `cheap_lane_ordinary_marker_envelope.py` |13457|`62ffd59bda5eccb6070ee8bcb7cd10be9f669f40257e5cb56f286504659439fe`|

## Thirty-third round: ordinary bounds over consecutive source sizes

The first band implementation uses exact rational vertices of the
nine-plane feasible region. Directed log-Gamma values bound the
source-tilted denominators from below. Both ordered terminal chains
are retained. The source count is summed once over actual sizes,
including odd sizes; no endpoint extrapolation or profile factor is used.

The fixed full-count cells expose a witness-selection problem:

| Actual sizes | Directed log2 upper bound | Worst hole vertex(k,a,b) | Interpretation |
|---|---:|---|---|
|5888--6400|+379700.0714671446|(2304,1023,1023)|No useful bound; extending the last cell assigns a near-full marker to258 singletons|
|6144--6400|+58476.1619984445|(1792,895,896)|No useful bound; nearly full profiles move into a marker calibrated for many singletons|
|6144--6160|-42710.6165412007|(2032,1015,1016)|Complete joint bound at35%, all17 sizes and every state path|

The successful narrow replay encloses

    log2 U = -42710.616541200724108195272662412622703018446654038599346201689978821828796698398892388496709866886056558909476978
                +/- 6.03e-109.

It uses496 nonempty unordered cells,961 ordered cells, and3970 exact
rational vertices. Their SHA-256 is
`0a75c05ef3611065959d99f4b3799ee911d358256a222cb133055a9e25168726`.
The indexed-count checksum is
`f47c344fe66972e0d097589e038567519857f2639510a55d9291ddb3512723dc`.
The eta-weighted count has log2 approximately6659.6291812; the output
correction has log2 exactly409436.0484375. The beta source factor
appears only in the256 denominator tilts. The local coefficient hash
remains `75ac478004e1572d716752cf989cee3cf1182fdf98e17ad67a401ba5b6c7bab3`.

All2066 geometry checks,541 scaling/count/coverage checks,23 directed
Gamma and singleton-cell checks,409 ordinary controls, and inherited
finite-chain controls pass. Nonnegative intersections are zero.
The selected replay takes13.207 seconds and peaks at103710720 bytes.
The two failed wide replays take14.281 and13.343 seconds. An initial
wide run also completes, but its console output is truncated; the
14.281-second repeat supplies the complete retained in-memory receipt.
These positive majorants do not exhibit low-weight codewords.

Reproduce the successful fixed-cell calculation with

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_ordinary_hole_band.py --first 6144 --last 6160 --seconds 300

### Align the proof cells with relative occupancy

The next calculation indexes the existing marker cells by x=2048a/k,
not by a alone. Half-integer cuts cover the complete real interval
[0,1024], including normalized targets from odd source sizes. The
resulting rational polytope still has nine planes, so the same directed
log-concavity argument applies. No marker, local operator, setup
distribution, or encoder work changes. This repair addresses the
observed mismatch. Independent source reviews pass before replay.

The matched directed A/B at5888--6400 now closes the whole interval:

    log2 U = -25386.799313271629211921633896057139011512901226035035590301670054845447884000219080235594723458376451289002512913
                +/- 7.10e-109.

The marker values, local coefficient hash, indexed-count checksum,
eta-weighted count interval, and output-correction interval are identical
to the failed fixed-cell run. Both programs use the same ordinary
numerator formula and frozen evaluator. Only the cell-to-marker
assignment and denominator minima change. The old bad profile now lies
in center893's cell, not the near-full center1024 cell.

This replay checks all513 integer hole sizes and525697 own-size
matching marginals. The worst normalized cell has centers(44,44);
its minimizing vertex is(k,a,b)=(1792,791/16,791/16), or actual6400.
The normalized full counts are exactly113/2. The bound retains this
fractional vertex rather than rounding to an integer support. Runtime
is16.441 seconds, peak102744064 bytes, and nonnegative intersections
are zero. The10057 normalized-geometry tests include80 cells,5409
feasible sample points,5212 integer marginals, and26 linear objectives.
The inherited exact and directed controls also pass.

A separate larger replay closes5632--6528, all897 actual sizes:

    log2 U = -4134.32114329785057052817840783057367510830739207205763936603626435321027198327728119850509712813169106732891317
                +/- 3.14e-108.

It checks947905 own-size marginals. The worst cell is(11,11), at
hole vertex(1664,481/32,481/32), actual6528, and normalized counts37/2.
The eta-weighted count has log2 approximately7352.4146080; the output
correction remains409436.0484375. Runtime is17.009 seconds and peak
102436864 bytes. Both successful normalized runs use496 nonempty
unordered cells,961 ordered cells,3970 exact vertices, and zero
nonnegative intersections. Their compact reproducibility checks are:

| Actual band | Vertex SHA-256 | Indexed-count SHA-256 |
|---|---|---|
|5888--6400|`38b3582a7e33b8ec5dd69812bbe1a851a932a80f12db0dafe4b7ec49b2495f90`|`bc2868c7cde74519edb71f9567b255f45a265b48e69e920326f9e3b12486eccb`|
|5632--6528|`8f2a8bdcba590bdddd40aa8f02564ddcd7718f6bf240fd9134669f36cc9cb7d6`|`d10f533c3209c0a9a095db07ddf4a050e2f6ae9dc1d967f3f4f1f116c82dd903`|

The next fixed-bank extension,5120--6528, gives directed log2
**+100297.2268357221**, not a useful guarantee. Its worst cell is(3,3),
at hole vertex(3072,9/4,9/4), actual5120, with normalized counts3/2.
It checks all1409 sizes and1669313 own-size marginals, takes17.550
seconds, and peaks at102928384 bytes. All tests pass. The positive
majorant marks a limit of these fixed witnesses, not an encoder defect
or a reason to discard the certified5632--6528 band.

Reproduce the successful normalized calculations serially:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_normalized_hole_band.py --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_normalized_hole_band.py --first 5632 --last 6528 --seconds 300

The next experiment should fit a new rational witness bank at actual
size5120, holes3072, with z1/32 and unchanged alpha119/1280, then
directed-test4864--5376. This z is an explicit proposed witness, not
an interpolated guarantee. Reuse the local operator and discard
no indexed supports. The frozen hole screen already contains the
bounded marker-fitting routine; no full coefficient table is needed to
propose the witnesses. Near4096, the ordinary entrywise-minimum result
does not promise enough slack for this coarser band envelope. Preserve
the stronger coefficient certificate and test that regime separately.

The897-size joint certificate is progress. No construction, alphabet,
refresh distribution, or encoding work changes. The whole-code goal
remains active, with its exact remaining gaps in the top table.

Frozen new sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_marker_band_geometry.py` |8089|`53685cda5b56f70f70820bb195bbd23d8453635a539360257b6d5f5a51dfa01b`|
| `cheap_lane_ordinary_hole_band.py` |18280|`ca4147fc8eb19aed544aca76eabc1ace357974b2bcd026474cac1465afb753f4`|
| `cheap_lane_normalized_marker_geometry.py` |9597|`00a4fde7841efed7c9c58b199e958c9a799bb190e4724ee9332d197e8ad9dee4`|
| `cheap_lane_normalized_hole_band.py` |14007|`8686ec42dd814974efeb3ec9230d9f4da898b1df83712fec38ed1acce2b7983e`|

## Thirty-fourth round: recalibrated witnesses close the larger-support gap

This round changes proof witnesses only. The actual outer, route, and
b(theta+a) refresh remain unchanged. A bounded floating-point fit proposes
31 positive rational marker pairs at a chosen actual source size. The
directed verifier then evaluates those witnesses using the same local
coefficient table. Floating-point fit values do not enter the final bound.
The table's exact rational hash is identical before and after fitting;
its Arb entries are preserved by identity. Hole coordinates relabel the
complete degree triangle once, without changing actual-input constraints.

All runs use alpha=119/1280 and 384-bit directed arithmetic. The normalized
cells retain their reference cuts, exact fractional vertices, and both
ordered terminal chains. Every actual size, including odd sizes, receives
its own source factor and support count. The count stays outside alpha;
the beta tilt appears in exactly 256 denominator factors.

| Actual sizes | Fitting anchor | z | Directed log2 upper bound, approximately | Result |
|---|---:|---:|---:|---|
|4864--5376|5120|1/32|+1097.9447477|No useful bound|
|4864--5312|5120|1/32|-7769.1704778|Complete joint bound|
|5313--5631|5504|1/48|-39803.8836851|Complete joint bound|
|6529--7040|7168|1/128|-8927.1121643|Complete joint bound|
|7041--7551|7168|1/128|+12327.6857486|No useful bound|
|7041--7296|7168|1/128|-40160.4544123|Complete joint bound|
|7297--7551|7424|1/256|-61719.8532553|Complete joint bound|
|4352--4863|4608|3/64|+28379.0780235|No useful bound|

Every displayed interval has radius less than 6.4e-108. The five successful
new bands therefore have conservative integer margins 7769, 39803, 8927,
40160, and 61719 bits. The previously verified 5632--6528 band contributes
4134 bits, and the dense 7552--8192 band contributes 781 bits. Their ranges
are disjoint and cover every size from 4864 through 8192. Adding their probability
ceilings gives

    U_(4864..8192) < 2^-781 + 6*2^-4134 < 2^-780.

This union requires no independence between bands. Different deterministic
markers and z values bound events for the same setup distribution. Do not
add nested singleton receipts or the earlier overlapping bands again.
The whole code remains uncertified because smaller support-size gaps remain.

The failed wide bounds are positive majorants, not low-weight codewords.
For 4864--5376, the worst denominator vertex is at actual size 5376; narrowing
the band to 5312 and fitting a separate bridge closes the gap. Likewise,
the 7168 bank's wider upper trial is worst at 7551. Splitting that range
and fitting a 7424 bank closes it. A global envelope can combine one
size's large support count with another size's worst conditional bound;
its failure does not establish failure for that endpoint alone.

The lower 4352--4863 attempt is also worst at 4863, with normalized full
counts 37/2 in both matchings. It adds no coverage. Next try a narrower
same-bank interval before changing the operator. Treat 4097--4351
separately: the complete 4096 result uses a stronger coefficient bound,
and its ordinary shortcut selects entries across all 31 markers. Neither
receipt transfers to this coarser one-marker-per-cell band envelope.

The successful replays cover 697185, 434079, 360704, 131200, and 98111 own-size
matching marginals respectively. Each checks 496 nonempty unordered cells
and 961 ordered cells. Their exact vertex counts are 3970, 3970, 3972, 3968,
and 3968. All wrapper, fitter, geometry, count, ordinary, directed-Gamma,
and finite-chain controls pass. No fit fails; nonnegative intersections
are zero. A fit may exhaust its proposal budget without invalidating its
positive rational witness: only the independent directed replay certifies it.

Compact reproducibility pins for the five successful new bands:

| Actual band | Vertex SHA-256 | Indexed-count SHA-256 |
|---|---|---|
|4864--5312|`016cc66a4b1fa55caf04ecc08ae9579a6cba6f203d020c56480554a30e8591ea`|`9c17823f77e7dadb00b4365c523c8faaaa558920fed831349bf899fed9988b17`|
|5313--5631|`37de6ba65da24978cc4f9e3755c01cbb434905405d59b282ded21469cd8f1b67`|`402e3488ffc9104570925791aa84f3ab39e96dc2ab6eaef3623997ab0747ac5a`|
|6529--7040|`fda29704fecb8f0279a8d2004e621f2e7ab2018b94bebbe7199b8f6bcfb155fc`|`96bbe697263a22a698dc2b6b17fc8ef248a242436056545df2fde509cca8e34b`|
|7041--7296|`266b377ca2a4575375449da2d3262b901482b763ef1292203151e40a2967db89`|`e819c5c934484d066e973c7c53d6f49a0657b322f1caeeaaa51da461ee462d49`|
|7297--7551|`c8720d0aedf68dd88d7564e6f169783bf9d51465cb89aa9414b587d661748adb`|`239a0e97bfa9a87c62483f3c500cc06b1360b4930aeb6343f7af86b0bd2e4bea`|

| Fitting anchor,z | Positive rational bank SHA-256 | Local coefficient SHA-256 |
|---|---|---|
|5120,1/32|`5be227ae51d8a53f2237dc99ebef26c0ad223c98b17a9b97bc10a75f104ebccd`|`723382ab49dc80624ecdedd0d5e835eb8337171b2e28e8196c5593eb61135653`|
|5504,1/48|`e70306d2bd6001717d23cadcf801ac1abab2f5159b629c986125d666b431d0f1`|`dd9a895c7090bc507a8b795e565eb6c7383b429fe1fa83ed215b705e72fd74f0`|
|7168,1/128|`2ef0d19d17fcadf1fc31439c620f0c7ad1df6ad60c3594d992a9da3a065c67f2`|`147c70f7204f41a8091cd86890831d79aa37e04fac943e1c4f8fc2f281906d6b`|
|7424,1/256|`99410f3c26082fc4c97c8d71c79009c5e363919d41d8ee7ceccdb41abb3702d8`|`6d4bfdd0045d95c6fadf956b3696b9b7f290b06b55ca22a7f1032651e37bffe8`|

An independent source-only audit also rechecked why the scalar source
factor applies to actual F-valued messages. The proof does not count
those messages or assume a support-preserving scalar projection. Its
potential P_q(d)=(q^d-1)/(q-1) uses the surviving F-dimension d, with
q=256. Fresh restricted coefficients contract that potential per pivot.
Only the initial E-defined constraint kernel identifies the potential
with an auxiliary scalar count. The
[global rank reduction](CHEAP_LANE_RANK_INTERFACE.md#a-global-reduction-to-scalar-constraint-counts)
and the joint missing-group argument provide this bridge; the scalar
outer density already includes the source-message count.

The audit found no new mathematical gap in this dependency chain.
Fixing the outer for the refresh-rank argument does not mean conditioning
its later average on the resulting kernel. State masks remain branch
labels, not observations of future randomness. The claims still require
the ideal independent setup, the specified field embedding, RS
information-set properties, and the cited invertible-Toeplitz count.
Neither a benchmark seed nor numerical interval arithmetic proves those
structural assumptions. No old full-field certificate is transferred.

To replay any row, run the following command serially, replacing the four
capitalized arguments with that row's exact values. The program prints
all 31 rational witnesses and their hash alongside the directed result.

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_hole_band.py --anchor ANCHOR --z Z --first FIRST --last LAST --seconds 300

Frozen new sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_hole_marker_fit.py` |15012|`d926a1283b4fb016687e6d5fdc66ce59d031db79a465cc015c3a81888ca8ec0d`|
| `cheap_lane_calibrated_hole_band.py` |13838|`95397e9c6a3fd51272d79b36d911fa0ef144927a479f5be44a7a3226dd2d9354`|

## Thirty-fifth round: actual coordinates and cell-specific run counts

This round changes only proof witnesses and support accounting. The encoder,
refresh distribution, and local pair operator remain unchanged. All results
below use 35% distance, alpha=119/1280, and 384-bit directed arithmetic.
The reported log2 bounds are rounded summaries; positive values provide no
certificate and are not codeword counterexamples.

First, narrower hole-coordinate bands extend the larger-support cover:

| Actual sizes | Fit anchor | z | Approximate log2 union bound |
|---|---:|---|---:|
|4352--4608|4608|3/64|+15182.2685109|
|4608|4608|3/64|-4666.0505301|
|4544--4672|4608|3/64|-1709.3807264|
|4673--4863|4736|5/128|-7164.6913905|
|4352--4479|4352|7/128|+13960.7233435|
|4352|4352|7/128|-497.8311219|

Together with the seven disjoint bands from round34, the two successful
new bands cover every size4544--8192. Their joint bound exceeds780 bits:
the dense term is below2^-781, and each of the other eight terms is
below2^-1709. The stronger singleton receipts remain separately indexed.
The failed wider bands show that a large singleton margin does not justify
interpolating neighboring sizes.

The new actual-coordinate driver retains the unreversed H/J table and
uses the source exponent r, not8192-r. Exact tests check that distinction,
odd sizes, the4096 endpoint, rational vertices, and both terminal orders.
The following complete checks use fit anchor3072:

| Actual sizes | z | Approximate log2 union bound |
|---|---|---:|
|2944--3200|1/8|+14895.1340005|
|3072|1/8|+11520.1977878|
|3072|1/10|+7548.0695023|
|3072|1/12|+8250.0518649|

None adds coverage. The last control worsens the bound, so this round
does not continue a blind z scan. Its final output was inspected; only
a compact receipt, rather than the full console report, is retained.

The next refinement retains how many indexed supports can occupy each
normalized cell. The
[run-count lemma](CHEAP_LANE_RANK_INTERFACE.md#actual-coordinates-and-cell-specific-run-counts)
uses exact cyclic run counts as upper bounds on cell populations. These
are not exact joint-profile counts. Each actual size keeps its own source
factor; off-diagonal cells charge both terminal orders. The new calculation
minimizes two completed bounds, never factors from different bounds.

Three independent source reviews found no blocking issue. The numerical
replay passes4636 exact checks, including1338 nonconstant binary masks,
21 weighted-union controls, and noncommuting terminal chains. Matched
controls reproduce the prior binomial-times-maximum bounds:

| Actual sizes | Anchor,z | Old log2 bound | Run-aware log2 bound |
|---|---|---:|---:|
|4352--4479|4352,7/128|+13960.7233435|+9366.1281058|
|3072|3072,1/10|+7548.0695023|+2015.5551480|
|4352--4383|4352,7/128|+195.2138941|-221.7124546|

The last row closes all32 sizes jointly, including16 odd sizes. Its
directed interval is centered at
-221.71245462134066631303730962092989386820943858358790343762988441365111137257400251803206609203791727176111406
with radius5.07e-108. This is new coverage, not just a stronger singleton.
At3072 the leading weighted cell moves from centers(1013,1013) to(855,855).
The improved bound remains positive; no whole-code claim follows.

Next preserve the relation between run count and conditional cost inside
the leading cells. Split selected cells into disjoint run intervals and
add c_lo<=k-a-b<=c_hi to their vertex geometry. Exact run-prefix masses
can count each interval without changing local tables or the construction.
Compare completed refined cell contributions with their existing bounds;
unrefined cells must remain in the total. This is proposed work, not
additional coverage.

The following witness hashes permit matched replay. The same local table
is preserved before and after fitting, including Arb object identities.

| Anchor,z | Rational bank SHA-256 | Local coefficient SHA-256 |
|---|---|---|
|4608,3/64|`eb1bcb6fafcf84131e010b597b2a72755257e4a6a72b286f1a24f3920d0dc2c7`|`90536197e6c03ec489234fe006680c4d83cb7df10a97408fd167d84ac788d425`|
|4736,5/128|`d7b2276669935b0d52db1dd021ded97fed78bca9f503321c25117f4b2ea16060`|`22b36d94bd743e8476523f4bb0ad37c37bd73a194d6936e6f0bdd5c8bfd305c2`|
|4352,7/128|`8e40f93f99c8e68890c984cac7f460dc09a2ef8474fadb28752fd31892736280`|`fb42e12f060ef89385599d3b5a7df93eb1dfd32b6e203ac263699f788af4d8c7`|
|3072,1/8|`4c00e5a3c13664b4daf37d45b6ab212c84037ba383ae1a571f90ca0f573d890d`|`06fe434245e9023cab83ee056a09a99eb2de56bc290250263ed214643d86297b`|
|3072,1/10|`3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7`|`636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d`|

| New certified band | Vertex SHA-256 | Indexed-count SHA-256 |
|---|---|---|
|4544--4672|`7d73680e42535591f6c758e0cb596e582b1bd5ed66802bbc93350dba3f63620f`|`fc7f7b1b514ffe1b8aaf4dbce950b6aaf06eb8038e23392077cc6e0c2398a915`|
|4673--4863|`33ee185b08bb636ae4d212f99ffe6dd2b4c16a957e0eb41d6bbd58f56f5437cb`|`fd6fbdf27bac1889d9e3c3577596f6710be8a3976ab9c9e2c9a52af295fd975d`|
|4352--4383|`22659a6a3bbfc3574b0d89c3edab71fb30e08d63c97b103a635ea53d64994b9e`|`ef47e7a171f0e2a84094eb7ba2e337dba3f69cbeb047feba0c2ee6f5b7286ecc`|

For4352--4383, the run-interval hash is
`56b65e392cdd26c38c7d2967c4dca736ad4754dfb2d9dbf5ce3d3789540d7612`,
and the unweighted cell-mass hash is
`a8028e04d0e0035a1a69342323562ce21eca134c8c9439dc86005a517f665270`.
All runs check each actual size's exact binomial checksum.

Replay the first table with the frozen hole-band driver from round34.
Use `cheap_lane_calibrated_actual_band.py` for the second table and
`cheap_lane_calibrated_run_band.py` for the third. Pass the table's
`--anchor`, `--first`, `--last`, and `--z`, with `--seconds 300`.
Both new drivers default to alpha119/1280 and384-bit precision.
Run serially with `C:/Python314/python.exe -B`; output goes to stdout.

Frozen new sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_calibrated_actual_band.py` |22399|`ca72174f0fc19848e02de159093f476eb61369629d8dc6f1271332171e7ae635`|
| `cheap_lane_cell_run_counts.py` |20201|`2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028`|
| `cheap_lane_calibrated_run_band.py` |15812|`19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2`|

## Thirty-sixth round: run intervals isolate coefficient-bound slack

This round keeps the encoder, refresh distribution, and local operator
unchanged. It first refines the complete size3072 calculation at35%
distance, then compares two bounds on one integer profile. The final
complete regional replay closes every3072-group support with3356.7491
bits. Neighboring sizes and whole-code closure remain unclaimed.

First, changing z from1/10 to1/12 worsens the complete run-aware log2
bound from+2015.5551480 to+2829.4417754. The old global maximum becomes
+8250.0518649. This control takes16.535 seconds and ends the z scan.

The new run-interval calculation implements the
[restricted-cell lemma](CHEAP_LANE_RANK_INTERFACE.md#keeping-run-counts-and-conditional-bounds-together).
It partitions each selected cell's count into disjoint integer run
intervals and restricts the conditional bound to that same interval.
Both terminal orders remain inside every off-diagonal contribution.
Unselected cells stay in the complete sum; minima combine only completed
bounds. The source factor uses actual r, and the output factor appears once.

The exact geometry checks cover singleton sizes and runs, fractional
vertices, and independent clipping controls. All14942 checks pass:
180 cells,45 full-run-range controls,3723 sampled feasible points,
26 linear objectives,68 integer partitions, and2244 integer targets,
including1140 odd targets. The standalone test takes3.924 seconds.

All following complete bounds use r3072,z1/10,alpha119/1280:

| Calculation | Approximate log2 bound | Time, seconds |
|---|---:|---:|
| Previous global maximum |+7548.0695023|Control in each replay|
| Previous run-aware cells |+2015.5551480|Control in each replay|
| Top16 cells, run-interval width32 |+1852.1997253|20.800|
| All496 unordered cells, width8 |+1546.8463520|61.472|
| Same all-cell replay, requested precision256 |+1546.8463520|61.377|

The final row is a reproducibility check, not an independent precision
check. The frozen local builder enforces at least384-bit Arb arithmetic,
so both requests use the same working precision. The first all-cell
console receipt was truncated; the compact replay retains a complete
receipt. There are no omitted cells or unresolved numerical fallbacks.

The complete all-cell interval is centered at
1546.84635201107033851008216750448238055656545412257231986148764325985241264489083937497523514057520581459683102
with radius5.19e-108. It remains positive. All6508 positive cell intervals
have nonempty geometry; the replay checks26003 rational vertices and12610
terminal orders. Exact interval masses partition the previous cell counts.
The partition hash is
`f601dcaf5ceace1c74cdb2e1aa45e668e86e198dae0a57c116996c78de3d923c`;
the vertex hash is
`667639c216ab7c084d7a64ad39329020a797387fc7e66b0a44b8ff1b49b92180`.

The leading cell has marker centers(813,980). Its largest interval is
433--440 runs, with minimizing fractional vertex
(r,a,b)=(3072,4743/4,5805/4), hence435 runs. Rather than continue narrowing
intervals blindly, the next test evaluates nearby integer full counts
(1186,1452), with434 runs, under those same two markers.

### A matched coefficient calculation repairs the selected profile

The selected calculation compares ordinary Cauchy evaluation with the
positive residue256 coefficient projection. Both use the same immutable
rational local table, regenerated31-marker bank, own target denominators,
actual-r source factor, and full434-run-class population. The latter is
an upper bound on the selected profile's population, not its exact count.
Both terminal orders are retained separately. Positive aliases and upward
integer rounding remain in the residue calculation; it is not exact
coefficient extraction.

The mandatory expected-hash checks in these replays use local table
`636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d`
and rational bank
`3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7`.
They match the complete run-interval calculation. No local cap or marker
is changed between the two bounds.

The accepted rational witnesses are preserved below as center:u,v.
They are proof parameters, not encoder coefficients. These values allow
the witness bank to be reconstructed independently of the FLOAT proposal
procedure; the normalized cuts remain those of the frozen31-center bank.

```text
0: 195/234398, 1817920/462587
3: 3466/464713, 3582419/917468
11: 477/17644, 3986783/1040888
25: 56905/944681, 2900597/783141
44: 93887/909154, 238815/67462
69: 154644/987503, 2274981/681757
98: 214179/999802, 3161965/1014259
131: 97112/353349, 330662/114549
169: 4764/14059, 1264777/478590
211: 291065/721632, 1035541/432019
256: 464915/997473, 2196558/1017953
304: 312346/592925, 1060460/550467
354: 4420/7567, 1747447/1023105
406: 424973/665854, 320495/213413
458: 63020/91663, 1001722/762307
512: 372912/507877, 390973/344112
566: 532160/684837, 213034/218851
618: 692653/849852, 477283/575573
670: 356029/418766, 332786/478201
720: 249254/282715, 328656/569743
768: 671389/737860, 462433/983691
813: 796549/852035, 46051/122444
855: 375581/392468, 70944/241961
893: 954314/977717, 16136/72713
926: 182602/184069, 120241/739199
955: 263671/262188, 63337/562898
980: 124343/122246, 39689/561608
999: 45884/44731, 33080/833413
1013: 1037660/1005399, 13871/801065
1021: 631177/609429, 4175/888433
1024: 607405/585798, 501/961088
```

| Terminal order | Ordinary log2 bound | Residue256 log2 bound,192-bit integer rounding |
|---|---:|---:|
|813 then980|+1468.8901601|-3379.9728676|
|980 then813|+1468.9143611|-3379.9492332|

Thus coefficient projection alone recovers about4848.86 bits for this
profile. The directed192-bit-integer replay takes66.054 seconds and
peaks at290136064 bytes. Its26 adapter checks,80 extraction checks,
1181 integer-product checks,409 ordinary checks,36 finite-chain checks,
840 enumerated supports, and18 bank checks pass. Source review was
independent of implementation. All probabilities still concern the
implemented restricted coefficient distribution.

A256-bit integer-rounding control reproduces both displayed residue
bounds. It takes76.118 seconds and peaks at350720000 bytes. Its local
coefficients, rational markers, ordinary bounds, source factors, and
run charge are identical. Only the regional integer quantization changes;
both runs use384-bit Arb arithmetic. Their small nonoverlap beyond the
displayed digits is consistent with two different upward-rounded upper
bounds, not an interval enclosure for one identical exact quantity.

This diagnoses a substantial loss in the ordinary coefficient bound.
It neither certifies the full434-run class nor every size3072 support.
That result motivates the complete regional calculation below. It retains
own denominators, source counts, both terminal orders, and a complete sum.

### Complete size3072 certificate

The complete driver builds31 cyclic tables sequentially and extracts
every own-weight marginal f=0,...,1536. Each target has s=3072-2f and
its own multinomial/marker denominator. No size4096 endpoint shortcut
applies: here f=0 does not force every local pair to have degree(0,64).
The driver uses no boundary replacement, ordinary fallback, or new alpha.

The exact run-mass join covers2362368 ordered matching profiles, sharing
1181952 matrix powers while keeping both terminal orders. Its indexed
counts sum to binom(8192,3072). The actual3072 source/output factor is
applied once afterward, with h=12,c_h=1/13 and
eta=18837534351525998328207900364912499/65025. The result is

    log2 U = -3356.7491176288611133963989569096663623099291764354924315215656796089412766978847349530488304140983050546941498112
                 +/- 8.53e-110.

This bounds the ideal-setup failure probability for all messages with
exactly3072 active source groups, over all indexed supports and state
histories, at35% distance. K remains2^20 and the rate remains1/2. It is
not a whole-code certificate. The largest run-class contribution is at
c=527, with log2-3362.5471535; c=528 and526 follow closely.

All186 cyclic-DC controls pass. The selected434-run controls reproduce
their earlier192-bit-integer intervals exactly from the complete tables;
they are not added separately to the join. The inherited exact tests,
both endpoints, source checksum, and complete assignment checks pass.
Nonnegative intersections are zero. The local-table and marker-bank
hashes are the same enforced values recorded above. The assignment hash is
`41366759d2f56cfa87f83f61150de6fc310bf9f3d5a4402ed226b19f04a35088`;
the indexed-count checksum is
`7c976f57a072427f4b0dda9c8d3a6b636b9f0de5747bab312ade647737602ea2`.

The replay takes875.753 seconds, of which40.366 seconds are the complete
join, and peaks at294932480 bytes. It uses384-bit Arb arithmetic and
192-bit regional integer rounding, with positive residue aliases retained.
Source and orchestration reviews are independent; the parent performs
the sole numerical replay. No encoder work, generated files, or baseline
certificate changes accompany this new coverage.

Next use the shared regional tables to test3064--3080 jointly. Every
weight needs fresh extraction of(f,r-2f), its own denominator, and its
own source factor; h changes at3072. The existing count-aware band
envelope can share terminal powers without extrapolating from this margin.
Retain a singleton3072 reproduction control before claiming a new band.
An entrywise maximum over positive residue coefficients is another valid
way to reuse the denominator geometry over broader intervals. Its gain
is at most16 bits per regional factor,4096 bits across the full chain,
relative to the same exact ordinary numerators. This untested envelope
may lose the selected-profile gain; test it before relying on it.

Reproduce these calculations serially:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 3072 --first 3072 --last 3072 --z 1/12 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_slabs.py --top-cells 496 --slab-width 8 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_selected_region_ab.py --bits 192 --expect-table-sha 636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d --expect-bank-sha 3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_complete3072_directed.py --bits 192 --seconds 1500

Frozen new sources:

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_run_slab_geometry.py` |11429|`9428a131913d6d6ac590393c4b73b988889d5e327bce54362ae1f6fecdeeb7b9`|
| `cheap_lane_calibrated_run_slabs.py` |24828|`786dc6cc9528c54e1921784d343db6694c68a0a88157012b62cfbb4336a66770`|
| `cheap_lane_selected_region_ab.py` |16035|`cf4d025277a281aecfd00b105bd7f791ca00c9864838b4849f97d1da7a2de26b`|
| `cheap_lane_complete3072_directed.py` |14670|`e599c0cd7ed7fd4f594f7b045aeb3702e9ded8fe0d2128ff0fd7dd6f47417b53`|

The rank interface also records an independently audited
[exact matching-imbalance enumerator](CHEAP_LANE_RANK_INTERFACE.md#optional-exact-count-for-matching-imbalance).
It tracks ones-run and zero-gap parity with two2-by-2 generating matrices,
retaining the imbalance discarded by the run count alone. This optional
count refinement is not implemented or numerically tested. It adds no
certificate and does not delay the complete coefficient calculation.

## Thirty-seventh round: reusable regional coefficient bounds

The new [common-numerator lemma](CHEAP_LANE_RANK_INTERFACE.md#a-common-numerator-for-all-regional-targets)
retains the positive cyclic256 calculation but replaces a selected residue
by the maximum over all65536 residues, separately for each of the six
regular/terminal entries. This loses some target-specific precision but
gives a numerator valid at every target degree. The existing actual-size
denominator, source factors, run counts, and ordered terminal chains remain
necessary. Neither this lemma nor a one-marker diagnostic certifies a
support class.

Two directed one-marker pilots use exactly the round36 local table
`636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d`
and rational marker bank
`3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7`,
with z=1/10,alpha=119/1280,192-bit upward regional quantization and
384-bit Arb arithmetic. The following ranges summarize all six entries:

| Marker center | log2(ordinary numerator / universal maximum) | log2(maximum / selected residue) | Seconds | Peak process bytes |
|---|---:|---:|---:|---:|
|813|12.8876333--12.8897259|5.9918965--6.0052253 at full1186|39.0848248|277454848|
|980|12.0625887--12.0639670|6.9245306--6.9382438 at full1452|38.1756582|258318336|

The maxima carry no target normalization, source factor, support count,
or additional alpha power. Each scan checks nonnegative legal entries,
zero padding, its own matrix exponent, and the exact inequalities
max<=DC<=65536*max. The directed comparison uses the same rational local
table before regional rounding. The ordinary/DC rounding discrepancies
are tiny but are not silently equated to zero. Both runs pass52 new exact
controls,1181 frozen integer comparisons,80 extraction controls,409
ordinary controls, and18 marker-bank controls, with zero generated bytes
and zero nonnegative-interval intersections.

The source-audited complete3064--3080 adapter is also preserved, but has
not been executed. The universal-max pilots motivate testing a wider
interval first. No coverage is inferred from the pilots or from the
previous singleton3072 margin.

Reproduction commands, from the repository root, one process at a time:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_maximum.py --center 813 --selected-control --bits 192 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_maximum.py --center 980 --selected-control --bits 192 --seconds 300

Frozen new sources (numerical use is recorded separately):

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_complete3072_band_directed.py` |16010|`f8cfe1932b1f56e15a8d81bbbe0bdab9548c5cae3758a7b4e33d03974c3c7d98`|
| `cheap_lane_residue_maximum.py` |20910|`d03a5dff0c607cac71f7b6da45c355c07654ca91bea879d726e58ce05c098036`|
| `cheap_lane_residue_max_band.py` |19617|`b2761383e6a929978764baa1af9f777e22907c662d505f352015e16f3d8a8938`|
| `cheap_lane_residue_max_replay.py` |20567|`77bd0acc081e68691810391c2cd4e9a5dad1f731f14296ec78d4de4831018867`|

The full-bank adapter uses31 such tables, exports186 exact dyadics, and
passes the unnormalized matrices directly to the existing finite-chain
cell evaluator. It does not raise H^64 or H^63 J to another regional
power. Each matching contributes128 denominator factors. The existing
beta tilt, actual eta_r^alpha, exact indexed run-interval masses, and
single output correction are unchanged. Independent source review checks
this interface. The singleton3072 control uses the weaker universal
numerators; it is not required to reproduce the stronger own-residue
certificate from round36.

The compact witness binds its coordinate convention, moduli, rounding
precision, z,alpha, local-table hash, marker-bank hash, and the31 exact
marker/maxima rows in one canonical hash. Later interval checks may reuse
these values without recomputing regional polynomials, but still depend
on the verified production of this witness. A hash alone does not prove
that arbitrary supplied numbers are valid upper bounds. Actual-coordinate
maxima must not silently be reused as a reversed-hole-coordinate table.

### Complete band and cached replay

The full31-marker calculation closes every indexed support of sizes
3008--3136 at35% distance. Its directed log2 union-bound interval is

    [-197.1555097232365514445287068507728162100500834404336490200839335778121628861118117848435528154042915644553794997 +/- 7.99e-110].

This includes all129 integer sizes and198241 marginal degrees, including
64 odd sizes. The count-aware join retains496 unordered cells,961 terminal
orders,3970 rational vertices, and each actual size's source factor.
In particular, h changes from11 to12 at3072. The old global-maximum
bound remains positive, about+5334.07: the coefficient maxima alone do
not close the band without the count-aware join. The same-table singleton
control gives1263.7385513 bits at3072. The stronger3356.7491176-bit
own-residue receipt remains valid and is preserved separately.

The generation run takes831.374608 seconds and peaks at291004416 bytes.
It checks186 DC comparisons,12189696 legal coefficients,12142080
zero-padding entries,372 matrix products,2232 polynomial products, and
434 quantizations. There are no nonnegative-interval intersections.
The band coverage hash is
`35e00dcd2da0aca11755184bef4b3bb0dfedd6ce557d73feb1c6f972886ef908`;
the indexed-count hash is
`92c6aa82c8e90ef78d679a644f71b2fa690d3e1bc748559d4967ff15c521a98a`.

The selected witness is retained in
`receipts/cheap_lane_residue_max3072_20261017.json` (13754 bytes).
Its canonical witness hash is
`7eb3fd910202b94c8dc1cbee2567d8c39c97f07651dba24ba659eefaae131dfb`;
the file-byte hash, including its final newline, is
`7feac146d45e16b4c70595474b0f7ffe0ea53a779565d6ab07f1261250c036a8`.
The file contains31 marker rows and186 exact positive dyadic upper bounds,
not the regional coefficient arrays. An independent metadata review
checks these fields and the canonical hash. A replay rebuilds the local
table and marker provenance, then recomputes target denominators, actual
source factors, run counts, geometry, and the final directed sum.
It does not independently reprove the retained regional upper bounds.

The cached3008--3136 replay reproduces the entire printed interval exactly
in18.394293 seconds, with zero regional polynomial products. The replay's
96 receipt-validation controls and inherited checks pass. Follow-up results
use this same bank; displayed logarithms below are rounded diagnostics.

| Actual support sizes | Directed log2 union bound | Seconds | Result |
|---|---:|---:|---|
|3137--3264|+2648.0163013|20.145595|No useful bound|
|3137--3168|-957.9651767|17.171926|Complete32-size certificate|
|2944--3007|+1225.7942184|17.020316|No useful bound|
|3264|+2646.3762002|15.655939|Endpoint also unresolved|
|2944|+1146.3153583|15.968623|Endpoint also unresolved|

For3137--3168 the full directed interval is

    [-957.965176655377302226171902919256540214850685906739306831443830566808079009901103922374161984835990518739906538 +/- 1.36e-109].

This includes50464 marginal degrees and16 odd sizes. Its coverage hash is
`bf483be39b93509477dc414b841c3cec8df183ee909d050352b52cf67c0d2237`.
The two successful bands are disjoint. Their union bounds are each below
2^-197, so their sum is below2^-196. Thus **every support size3008--3168
is jointly covered at35% with more than196 bits of ideal-setup margin**.
This establishes a161-consecutive-size band, not a whole-code certificate.
It adds160 previously uncovered sizes;3072 was already covered.
The stronger singleton3072 receipt is not added again.

The failed lower extension is led by marker cell(813,813); the upper
extension is led by(69,1013). Both outer endpoints remain positive when
tested individually. Consequently, removing cross-size envelope loss alone
cannot close either full extension with these witnesses. This does not
identify a low-distance word or an intrinsic obstruction. Next recalibrate
the proof witnesses near2944 and3264, then verify new regional maxima;
do not change the encoder merely because this fixed bank loses precision.

Reproduce the full generation and two successful cached checks, serially:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_max_band.py --first 3008 --last 3136 --bits 192 --seconds 1500
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_max_replay.py --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3072_20261017.json --expect-witness-sha 7eb3fd910202b94c8dc1cbee2567d8c39c97f07651dba24ba659eefaae131dfb --first 3008 --last 3136 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_max_replay.py --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3072_20261017.json --expect-witness-sha 7eb3fd910202b94c8dc1cbee2567d8c39c97f07651dba24ba659eefaae131dfb --first 3137 --last 3168 --seconds 300

The other cached checks differ only in first/last arguments from the table.
The full generator prints the exact witness to stdout; it writes no file.
Retain only its compact verified witness if regeneration is required.

## Thirty-eighth round: recalibration and marker selection

The previous round's new continuous band is progress. Its failed2944 and
3264 endpoint controls motivate recalibration, not an encoder change.
The construction, ideal setup distribution, and distance target remain
fixed. z and alpha below are proof parameters only. First rebuild the
same local operator and fit positive rational markers at each actual
endpoint. The frozen count-aware ordinary verifier gives these complete
directed results, with alpha119/1280 and384-bit Arb arithmetic:

| Actual support size | Fit anchor | z | log2 union bound, rounded | Seconds |
|---:|---:|---:|---:|---:|
|2944|2944|1/10|+3216.0814969|16.076951|
|2944|2944|1/9|+3027.3858620|15.828425|
|3264|3264|1/10|+2632.3703623|17.354916|
|3264|3264|1/12|-197.9821379|16.785346|

The last row closes **every3264-group support at35% with197.9821 bits**,
without any residue tables. Its full interval is

    [-197.98213792111999477769472041477517453225931728098819730190091704627730964047532631349008408114679288812658168 +/- 1.24e-108].

This is an ideal-setup bound over every indexed support and state history
in that size class. It neither extends neighboring sizes by inference nor
certifies a particular benchmark seed. The wide3169--3328 trial at the
same3264,z1/12 witness remains positive:

    [1998.04955520785245668849555612057716223250239343149590549056257071818317351657072324307032235487703611081493745 +/- 1.76e-108].

It takes20.097490 seconds and covers260000 marginal degrees. The new
recalibrated adapter's ordinary-only mode reproduces both complete printed
intervals exactly in21.977509 seconds. It exports no residue witness in
that mode. Its23 parameter checks and inherited directed/exact checks
pass. These are actual numerical controls, not a prediction that maxima
will close the larger band.

The four exact witness identities are:

| Anchor,z | Rational local-table SHA-256 | Positive rational marker-bank SHA-256 |
|---|---|---|
|2944,1/10|`636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d`|`d5f5e3823f8af0f97f6d10b39ee85b3ea6e5214bced6ef134386b0be9832f5d9`|
|2944,1/9|`ffbc03f0c930c5526443493261d64ebe87b221389a7dbebd956823aec0651e86`|`caba90b137e3d8d2cac14af466318811862c94a1ecadb6defe7b1f94001b4f87`|
|3264,1/10|`636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d`|`b46b1893992cd83e2a6b835c6182ce0c8bf5cd617a8fb42283d7131a45a637c5`|
|3264,1/12|`3576ff5128c4c51036873fcefd4f84317991922004ffbfe2675c3f8fbe55a3bb`|`99c97d71039b001257046d2ab1d1371a01a1f83cd28480b780c965d73b6aa3f3`|

The new adapter passes requested z and alpha explicitly through the local
caps, target denominators, source weights, and final output factor. It
never edits frozen module constants. Independent source review confirms
this interface and the distinction between ordinary and residue modes.
The subsequent3264,z1/12 maximum-bank calculation independently verifies
all31 regional tables. It closes3169--3328 jointly with the interval

    [-1308.630594006832561232456324565213948653921636038688198054900568170591891361986839130116837957261159932001236070 +/- 3.57e-109].

The same-bank3264 control gives3473.8373402 bits, strengthening its ordinary
receipt. Generation takes837.569032 seconds and peaks at295571456 bytes.
All186 DC comparisons,12189696 legal-coefficient checks,12142080 padding
checks,372 matrix products,2232 polynomial products, and434 quantizations
pass. No nonnegative-interval intersections occur. Every actual size,
indexed support, and state history in the requested band is covered;
the source factor changes at3328 and is evaluated at each actual size.
Its coverage hash is
`d53016b65541183e1c0226d3601f1614dee3970578931840435609e3c526a255`;
the indexed-count hash is
`0eabad2b3bea4ab8882b3a34a3c6065f9724b3ea1963f64feed7b43035d59da7`.

The selected witness is
`receipts/cheap_lane_residue_max3264_20261017.json` (13785 bytes).
Its canonical witness hash is
`76f2a88516c7a63a8f7cdf1968a8896ac4a6220a3d513a970530c2da6176e6ee`;
its file-byte hash is
`5f3224eed72cbd8c306f66646a5116086298d4ed8bf487de34443e69252a2578`.
The dynamic cached replay rebuilds the requested-z/alpha local table and
the anchor-fitted bank, checks every marker and the exact dyadic layout,
and reproduces the full printed band interval exactly in19.571413 seconds.
Its98 metadata controls and inherited tests pass. It performs no regional
polynomial products; validity of the retained maxima depends on the
verified generator receipt, not merely the matching hashes.

Cached upper-band checks then give:

| Actual support band | Directed log2 union bound | Seconds | Result |
|---|---:|---:|---|
|3329--3456|+635.5441593|19.397656|No useful bound|
|3329--3392|-2383.1517556|17.620009|Complete64-size certificate|
|3393--3456|+635.5441593|17.736772|No useful bound|

The successful extension has the full directed interval

    [-2383.151755612877875367172869754481495275111373030451225092151655444906392836616727130004273337843935536900554974 +/- 4.28e-109].

It checks107584 integer marginals across64 sizes, including32 odd sizes.
Coverage SHA-256 is
`f326339502e7bf0544f902a7b41e620768c6943be1daa86da16705e31fd83705`;
the indexed-count hash is
`aa6d65f2c588a951f640d440724f3a551edc91a98cb79d6d198f62a600b6670d`.
These are complete band checks with the verified3264 bank, not sampled
profiles. The positive remainder is an unresolved upper bound.

Together with the two previous disjoint bands, this closes **all385
consecutive sizes3008--3392 at35%, jointly with more than196 bits**.
The3008--3136 bound is below2^-197. Each of the other three disjoint
band bounds is below2^-199, so their sum is below2^-196. This adds224
new sizes relative to round37. Do not add singleton controls again.
The remaining whole-code gaps are separate obligations.

Reproduction commands, serially from the repository root:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 2944 --first 2944 --last 2944 --z 1/10 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 2944 --first 2944 --last 2944 --z 1/9 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 3264 --first 3264 --last 3264 --z 1/10 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 3264 --first 3264 --last 3264 --z 1/12 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_calibrated_run_band.py --anchor 3264 --first 3169 --last 3328 --z 1/12 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_max_band.py --anchor 3264 --first 3169 --last 3328 --z 1/12 --ordinary-only --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_max_band.py --anchor 3264 --first 3169 --last 3328 --z 1/12 --bits 192 --seconds 1500
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_replay.py --anchor 3264 --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3264_20261017.json --expect-witness-sha 76f2a88516c7a63a8f7cdf1968a8896ac4a6220a3d513a970530c2da6176e6ee --first 3169 --last 3328 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_replay.py --anchor 3264 --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3264_20261017.json --expect-witness-sha 76f2a88516c7a63a8f7cdf1968a8896ac4a6220a3d513a970530c2da6176e6ee --first 3329 --last 3456 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_replay.py --anchor 3264 --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3264_20261017.json --expect-witness-sha 76f2a88516c7a63a8f7cdf1968a8896ac4a6220a3d513a970530c2da6176e6ee --first 3329 --last 3392 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_recalibrated_residue_replay.py --anchor 3264 --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3264_20261017.json --expect-witness-sha 76f2a88516c7a63a8f7cdf1968a8896ac4a6220a3d513a970530c2da6176e6ee --first 3393 --last 3456 --seconds 300

| New source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_recalibrated_residue_max_band.py` |20314|`ed6745fb570200c27300a0da5e4c179ae4c998483c026592163f1428272fd619`|
| `cheap_lane_recalibrated_residue_replay.py` |18270|`d65b6701c9e22a97a586a3a3702c2634651cd4e9023243b094014e7779f88934`|
| `cheap_lane_residue_marker_envelope.py` |26473|`e66c63641a9a95f3b8880106fe9f526bc96738cfeefc2d93f7706bef01faf5a9`|

The [marker-selection lemma](CHEAP_LANE_RANK_INTERFACE.md#selecting-regional-witnesses-without-changing-support-cells)
provides another proof-only refinement. Keep each support cell and its
counts, but select among complete conditional bounds from different saved
marker pairs. Each candidate retains the same joint(r,a,b) polytope and
its own ordered-chain numerator and denominator. Including the assigned
pair makes the selected upper endpoint no worse. This audited lemma
does not itself claim a numerical gain or additional coverage.

The complete2944 test now checks all31-by31 marker pairs for each of496
unordered support-cell pairs, retaining both off-diagonal terminal orders.
This tests476160 alternatives beyond the assigned candidates. It improves
61 minor ordered cells but leaves the leading813/813 cell unchanged.
The completed log2 interval also remains exactly unchanged:

    [1146.3153582822572802066342800929646571995810492644669138915671383731967556675075054021335888057223698763130899222 +/- 7.05e-110].

The204 marker-selection controls and inherited tests pass. The run takes
36.067462 seconds and peaks at103739392 bytes. It rebuilds the local
coefficients and frozen assigned-cell controls but generates no regional
polynomials. This rules out a gain from this particular completed-pair
selection at2944; it does not rule out new markers or tighter cell bounds.

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_residue_marker_envelope.py --witness research/workstreams/char2_field_spin/receipts/cheap_lane_residue_max3072_20261017.json --expect-witness-sha 7eb3fd910202b94c8dc1cbee2567d8c39c97f07651dba24ba659eefaae131dfb --first 2944 --last 2944 --radius 31 --seconds 600

### Extending the actual-coordinate domain

Source review identifies641 as a historical wrapper limit, inherited
from the approved hole-coordinate interval. It is not a hypothesis of
the local pair bound. The local table depends on z and alpha, not global
support size. Actual denominators, normalized geometry, unique integer
marginal assignment, and indexed run counts already accept sizes1--4096.
The endpoint fitting proxies are positive for r>=2. Every currently
missing smaller class starts at225, so those formulas remain well-defined.

A new lower-domain adapter must test source-factor jumps and empty assigned
cells rather than silently relaxing a frozen guard. In particular,
for r<256 the existing source interface has h=0 and

    eta_r = 256 / [255 (256-r)].

This positive factor can be below one and must not be clipped. For each
new actual size, retain its own denominator and eta_r^alpha, with exact
counts outside alpha. Replay coverage and binomial checks at225,254--258,
507--508,511--512,515--516,640--641,767--768,and1023 before claiming new
intervals. Below342, endpoint proposal proxies may cross a neighboring
center's proposal target. No proof requires those targets to be monotone;
marker validity and assigned-cell coverage remain separate checks.
This is an extension plan supported by source inspection, not additional
numerical coverage or a lower-domain implementation.

## Thirty-ninth round: fill the small-source gap and isolate extraction loss

This round changes proof witnesses and coefficient evaluation, not the
encoder or its setup distribution. The old exact-coefficient verifier
closes the remaining225--254 gap at35% distance. The entire1--257 range
first has a joint bound exceeding42 bits. The current-pair refinement
below raises that joint bound above43 bits. This is partial coverage, not a
whole-code certificate.

### Complete small-source calculations

All three runs use `cheap_lane_cone_cells.py`, alpha119/1280, singleton
weight/run cells, and a180-second cooperative deadline. The reported
512-MiB figure in this frozen helper is an estimate, not an enforced or
measured peak. All requested profiles and state paths are included.

| Actual weights | z | Directed log2 union upper bound | Seconds | Use in the cover |
|---|---|---|---|---|
|225--232|267/311|+25.14708438|118.1833503|Only separately completed225--228 contributions|
|229--240|17/20|-80.59030544|125.8551791|Complete joint band|
|241--254|21/25|-115.86754259|147.1668053|Complete joint band|

The first band does not pass. Its fully computed singleton bounds at
225,226,227,228 are respectively-89.09322242,-75.38319174,-60.80121373,
and-45.34854186 in log2 scale; each is strictly below-45. The second
band replaces the overlapping229--232 contributions rather than adding
another copy. Its full directed interval is

    [-80.5903054395885259308553114267844654408548928339230468210525143433777191 +/- 2.38e-71].

The third band's full interval is

    [-115.8675425922578709218171727712287373997211210111523706390106627936938750 +/- 1.33e-71].

For transparent joint accounting, the earlier1--224 bound is at most
115*2^-50. The four new singletons contribute less than128*2^-50;
the two new joint bands contribute less than2*2^-50; and the existing
255,256,257 certificates contribute less than3*2^-50. Thus the complete
1--257 contribution is less than248*2^-50<2^-42. Keep the stronger
nested1--224 and individual255--257 receipts.

The first and third runs check respectively920/1743 marginal profiles,
1828/3465 disjoint weight/run points, and105808/217056 degree-feasible
ordered pairs. Each checks the exact indexed binomial count at every
weight,32767 terminal updates, and the inherited rational, support,
cofactor, polynomial-product, and normalizer controls.

### Ordinary evaluation is valid below the old wrapper cutoff, but loose

The new `cheap_lane_lower_actual_band.py` accepts actual weights below
the historical641 guard without changing the frozen helpers. Its exact
controls cover source-factor jumps, odd sizes, empty cells, indexed
run counts, and each actual size's own denominator. All6825 lower-domain
checks,4636 inherited run-count checks,1338 exhaustive support masks,
21 weighted cases, and409 ordinary fraction controls pass.

A3072 regression reproduces186 ordinary regional entries and the frozen
assembler's bound. The old assembler is only a regression oracle; its
reported interval is not an input to the new proof calculation. The
complete lower-domain trials are positive:

| Actual weight | z | Directed log2 upper bound | Seconds | Measured peak bytes |
|---|---|---|---|---|
|3072 control|1/10|+2015.55514796|20.0358047|102039552|
|225|267/311|+2404.54392161|13.1651470|86446080|
|516|17447/25000|+2063.20202624|13.6825624|87617536|

These are failed upper bounds, not low-weight codewords. At225 the old
exact-coefficient method already passes at the same z and alpha. The
next calculation therefore keeps current ordered pair caps but replaces
ordinary evaluation with exact downward coefficient rectangles. It must
retain odd weights, both terminal orders, and each actual denominator.

### Current pair caps with exact downward coefficient extraction

`cheap_lane_lower_exact_band.py` constructs the current ordered pair
operator once and shares four exact downward rectangles across the requested
sizes. Six positive polynomial doublings give H^64 and H^63 J. A discarded
degree cannot contribute to a lower retained target. Each extracted target
(f,r-2f) receives its own positive marker denominator and multinomial route
normalizer. Integer ceilings enlarge retained coefficients; degree
extraction is exact, but the rounded coefficients are upper bounds.

For each run count, the complete join takes the maximum over all feasible
ordered matching profiles, retains both terminal orders, and multiplies by
the exact number of indexed supports. Each actual size receives the current
joint-missing source factor once. In particular, the helper does not silently
substitute the stronger older source comparison below256. Independent source
review checks this distinction, odd weights, separate R/T scales, and all
assigned degrees. Runtime and memory checks are cooperative between operations;
they do not interrupt an in-progress FLINT call.

The225--232 pilot at z267/311 and alpha119/1280 now has directed joint bound

    [-334.5700202771383752879875807871115668777653208904041824172028934120674200773237054929683950802151634584735733047 +/- 5.39e-110].

It uses192-bit upward coefficient rounding and384-bit Arb arithmetic,
takes24.5096071 seconds, and peaks at89145344 bytes. Its singleton225
bound is-442.94384819. The local coefficient hash equals that of the
ordinary225 trial: regional extraction and the complete joining method
change between these two calculations. The exact control suite passes1181
fraction comparisons,80 normalizer checks,36 chain/cover checks,840 exhaustive
masks, and310 wrapper checks. The first23.7709-second run finished, but its
console report was truncated; only the complete second report supplies this
record.

Using this stronger band for225--232 and the229--240 result restricted to
233--240 replaces the four weak singletons above. Hence the current1--257
joint bound is less than(115+1+1+1+3)*2^-50=121*2^-50<2^-43. The earlier
42-bit accounting remains a valid weaker intermediate result.

The same helper at actual516, z17447/25000 and alpha119/1280, closes
all67080 ordered profiles and state paths with directed bound

    [-357.8913170366745189490989953182467975937968120040772405450013283710895879090507255081510861950270246392773563354 +/- 2.20e-110].

This192-bit coefficient/384-bit Arb run takes69.4445356 seconds and peaks
at192856064 bytes. The local hash matches the ordinary516 diagnostic.
The worst run class has34 runs, with log2 contribution-361.22592088.
Independent comparison against the historical37.1139-bit calculation
confirms identical source factor, threshold, route normalizers, and update
chronology. The gain is not attributed solely to extraction: the historical
calculation used the older one-packet cone and mixed cyclic/exact tables;
the new one uses stronger ordered pair caps and exact extraction everywhere.
The new run-mass envelope is conservative relative to exact profile weighting.

The257--289 trial at z21/25 completes in39.4997128 seconds, peaking at
104833024 bytes. Its full log2 bound is+260.79962803, so the full band
does not close. Nevertheless every separately completed contribution at
258--279 is below2^-60; the last is

    [-60.9024780070933240384133443226235998140906250100279590092366609832813517031653968727334974556086079052677130223 +/- 4.99e-110].

Thus these22 new sizes have a joint bound less than22*2^-60<2^-55.
Together with the preceding121*2^-50 accounting, every size1--279 has
a joint bound exceeding43 bits. Size257 has a stronger singleton control
of504.3825 bits, but is not counted again. Size280 gives only32.2748 bits,
281 gives2.9041 bits, and282--289 have positive bounds with this witness.
None of those sizes is added to the40-bit cover.

### A matched precision control closes the516--548 band

The516--548 extension at the same z17447/25000 and alpha119/1280 initially
gives log2+408.31617174 with192-bit coefficient rounding. It takes119.2858732
seconds and peaks at208207872 bytes. The separately completed516 bound
improves to-1205.55286268, despite the unchanged local coefficients. Several
larger even sizes then deteriorate sharply:542 gives-36.6801,544 gives+110.0794,
546 gives+258.4382, and548 gives+408.3162, while intervening odd sizes still
pass by hundreds of bits. The full band is not certified by this run.

The changed516 result is a numerical-majorant effect, not a mathematical
benefit from adding source sizes. Exact target coefficients are independent
of larger truncation limits, and exact marker factors cancel. Changing the
rectangle schedule changes the scales used for upward coefficient rounding.
Rare boundary coefficients can receive large relative inflation even when
the overall polynomial rounding loss is small. Both upper bounds remain valid.

The matched replay changes only `--bits 192` to `--bits 256`. The local
coefficient hash, four rectangle/marker witnesses, and all per-size source
factors are identical. It closes the entire516--548 band with interval

    [-842.2922173165695605971451734130813853487072430166441795273046860715520815058967021104132841627026073754427935997 +/- 3.26e-110].

The run takes136.6345318 seconds and peaks at264601600 bytes. All exact
controls pass. Its per-size bounds vary smoothly from-1205.55286268 at516
to-842.29222889 at548. Thus the matched experiment identifies rounding
inflation in the failed192-bit bound; it does not identify a code weakness.
Combining this new band with the earlier43.5679-bit508--515 certificate
gives more than43 bits jointly for508--548. There is no duplicated516 term.

### Recalibrated small-source band

The next33-size band,280--312, uses z677/825, alpha119/1280, and256-bit
coefficient rounding. Its z is the proposal
(5120+r0)/(5120+5*r0) at midpoint r0=296, continuing the successful small-
size witnesses. The proposal supplies no guarantee until replayed.
The complete directed replay gives

    [-198.4058008946141340530101082355784795820909478436549454749047517150604638330040111941727060247243061516850748417 +/- 1.25e-110].

It takes52.0444366 seconds and peaks at130834432 bytes. The local rational
coefficient hash is `c55a18c90b7b8e84b417c9229cad4b331ba6e171cd9ee9a5e7fa6ac53295a594`.
All exact controls pass; both endpoint sizes and every integer between
them are evaluated. Adding this band to the preceding cover gives

    U(1--312) < 121*2^-50 + 22*2^-60 + 2^-198 < 2^-43.

The new consecutive ranges are therefore1--312 and508--548, each with
more than43 bits jointly at35%. The seven remaining gaps in the progress
table are still required for whole-code closure. This round adds118
actual source sizes at the40-bit target, without changing the encoder.

Reproduce the completed runs from the repository root, serially:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_cone_cells.py --rlo 225 --rhi 232 --z 267/311 --alpha 119/1280 --seconds 180
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_cone_cells.py --rlo 229 --rhi 240 --z 17/20 --alpha 119/1280 --seconds 180
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_cone_cells.py --rlo 241 --rhi 254 --z 21/25 --alpha 119/1280 --seconds 180
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_actual_band.py --checks-only
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_actual_band.py --anchor 3072 --first 3072 --last 3072 --z 1/10 --ordinary-only
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_actual_band.py --anchor 225 --first 225 --last 225 --z 267/311 --ordinary-only
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_actual_band.py --anchor 516 --first 516 --last 516 --z 17447/25000 --ordinary-only
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 225 --last 232 --z 267/311 --alpha 119/1280 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 516 --last 516 --z 17447/25000 --alpha 119/1280 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 257 --last 289 --z 21/25 --alpha 119/1280 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 516 --last 548 --z 17447/25000 --alpha 119/1280 --bits 192 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 516 --last 548 --z 17447/25000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 280 --last 312 --z 677/825 --alpha 119/1280 --bits 256 --seconds 300

| Source | Bytes | SHA-256 |
|---|---:|---|
| `cheap_lane_lower_actual_band.py` |26835|`8481e9853a66973a1fd80551a7e7eee8d26333be92fff18afdaaf44760ec0b8a`|
| `cheap_lane_lower_exact_band.py` |19366|`fec69ce07aa584a63987de825b5cbdde9a507bdd7fdee902317e8d8426e3d115`|

The proof note also records an optional minimization over the joint
missing-set comparison parameter h. This is an independently audited
lemma, not yet used by the numerical helpers. The current runs keep
h=floor(r/256); do not credit them with the proposed improvement.

## Fortieth round: continue the exact cover and retune proof witnesses

The preceding round added complete source-size coverage and identified
upward-rounding inflation in a failed band calculation. That is verified
progress. This round continues the same actual-support proof at 35%
distance. The encoder, setup distribution, and baseline certificate remain
unchanged. All calculations use the frozen `cheap_lane_lower_exact_band.py`,
alpha119/1280, 256-bit coefficient rounding, and 384-bit Arb arithmetic.

The first new band,313--345, uses z5449/6765. Its complete directed bound is

    [-127.1541694839246331332508205305985651069043745067677706857205818789855411281426139522494123482068290890901228502 +/- 5.40e-110].

It takes63.3205926 seconds and peaks at141275136 bytes. The local rational
coefficient hash is `2b02ae29769aaaa48e93b24aae340906da670614acb0d574236492dab5e30f70`.
All exact controls pass. Adding less than2^-127 to the previous small-source
bound preserves more than43 bits jointly for every size1--345.

The next band,346--378, initially uses z2741/3465, continuing the same
small-size witness formula. Its complete directed bound is

    [-26.8648068531682126588701171265872183878488850090317658004392255847371092485467568166499761667318332889856990294 +/- 4.52e-110].

This is a valid bound but does not reach the40-bit working target. The run
takes72.3623004 seconds and peaks at159989760 bytes. Its first-size bound
has678.9681 bits, whereas the last-size bound has26.8648 bits. The local
hash is `3864f91c4d5a1a301dd5834cd3a2500e5ce926117d05d6cdfe5deff583a64c78`.
Rather than narrow the band, the next trial moves z toward the calibrated
weight512 witness while preserving all33 sizes, alpha, and arithmetic.
Interpolating witness choices supplies a proposal only; bounds are never
interpolated between certified sizes.

The retuned346--378 replay uses z783/1000, a rounded interpolation between
the successful z677/825 witness at296 and z17447/25000 at512. It closes the
same33 source sizes with directed bound

    [-440.5422463884899421264949756907376074367417186475676321875397556486828366007079770178527599371270262713959348464 +/- 2.88e-110].

It takes71.1860138 seconds and peaks at159707136 bytes. The local coefficient
hash is `fd125da3b8f10ff56417ece14409da1f378bf43de8e9628ca7dcbd9105f0c672`.
No source size, trajectory, or profile is removed. Only the proof parameter
changes; this is not a construction modification.

The next band,379--411, uses the same affine witness rule at midpoint395:
z15134261/19800000. Its complete directed bound is

    [-573.0880283646579450432310811417280719152086390376764890613907916543599493274925755200286520959254354922981951504 +/- 5.27e-110].

It takes83.2205090 seconds and peaks at173228032 bytes. The local coefficient
hash is `335a66de6e77480151e0cc9e1d5acd6274f14905ea299e6ef6b43349b8744a2f`.
All exact controls pass in both runs. Continuous small-source coverage now
reaches411, still exceeding43 bits jointly.

The affine witness rule also closes412--444. At midpoint428 it gives
z11072261/14850000. The complete directed log2 bound is
-720.8975050723885366864583474506930521989149438734386819777212181387914386743867108558331428411275515117904367981,
with interval radius6.20e-110. The run takes96.9965110 seconds and peaks at
201740288 bytes. Its local coefficient hash is
`50e26da8184b805c49d39f92d56a8e3fb184b68f49dd570672233cf150ff9523`.
All exact controls pass; the small-source cover now reaches444.

At midpoint461 the proposal becomes z8635061/11880000. It closes445--477
with directed log2 interval

    [-875.6321467219112806292626258807863088869365904988858327836406533597907206738963793460797409438145213957434161092 +/- 6.38e-110].

This run takes105.1004185 seconds and peaks at216047616 bytes. Its local
coefficient hash is `2a876d514521840839892aa243a997457934135228e93ccbd48e1107dedd213d`.
The source-size cover now reaches477; all exact controls still pass.

The final bridge band,478--507, uses z126285947/178200000. Its directed
log2 interval is

    [-1047.3891473854679853115551269219773794240949728507123202898187765722250845164875317082377297142312166186443948901 +/- 2.66e-110].

The run takes117.2665093 seconds and peaks at243200000 bytes. Its local
coefficient hash is `51fab6527e0a00233a353a47a1b0bde024d177c446deaa3d587654315702ccf5`.
All exact controls pass. This closes all195 previously open sizes313--507,
joining the older508--548 band into continuous coverage1--548. The nested
1--507 union exceeds43 bits; the combined1--548 union exceeds42 bits.

The next band,549--581, reuses z17447/25000 from the older516--548 run.
Its complete directed interval is

    [-86.9978854592273890757484905226401781263491809452938574511883576369523694616848837323156784708669817751135974374 +/- 3.12e-110].

The run takes153.1823579 seconds and peaks at282513408 bytes. Its local
coefficient hash is `538eef13aa632580fcedb464d33d994dda6c292ff3cc7864a3868c2b3133e25e`.
All exact controls pass. The first-size bound is825.0107 bits, but the
last-size bound is86.9979 bits. This drop motivates retuning z before
the next582--614 band, not extrapolating the current margin. Round40 adds
228 previously open sizes in total and extends continuous coverage to581.

### Joint failure budget for the covered classes

All probabilities here concern the cheap-lane code's ideal setup.
Let Q be the bad-distance event witnessed by a nonzero true-quiet history.
Let S be the union of bad-distance events from the currently certified
actual-source classes. Q may overlap S; an ordinary union bound suffices.
The40% quiet receipt bounds its35% subevent by2^-76.
For the older508--515 band,43.5679 bits implies a bound below3*2^-45,
since1/sqrt(2)<3/4. Together with the recorded disjoint size bands,

    Pr[Q union S] < 218*2^-50 < 2^-42.

This includes1--581,1024,2040--2056,3008--3392,4096,4352--4383,
and4544--8192, with no duplicate charge for nested source-size receipts.
It does not include unresolved sizes and is not a whole-code certificate.
The exact accounting uses the prior121*2^-50+22*2^-60+2^-198 small-source
bound, the new2^-127,2^-440,2^-573,2^-720,2^-875,2^-1047,and2^-86 bands, and the stronger existing
larger-size receipts. The leading508--515 term adds96*2^-50; all remaining
terms, including quiet, fit below the next unit of2^-50.

As a prospective budget rule, eight additional disjoint bands below2^-45
would preserve a covered-class union below474*2^-50<2^-41. This is a
conditional allowance, not coverage. The26.8648-bit trial is not admitted.

Reproduce completed runs serially from the repository root:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 313 --last 345 --z 5449/6765 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 346 --last 378 --z 2741/3465 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 346 --last 378 --z 783/1000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 379 --last 411 --z 15134261/19800000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 412 --last 444 --z 11072261/14850000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 445 --last 477 --z 8635061/11880000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 478 --last 507 --z 126285947/178200000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 549 --last 581 --z 17447/25000 --alpha 119/1280 --bits 256 --seconds 300

## Forty-first round: reuse exact tables across wider size bands

The preceding round adds228 previously unresolved source sizes. This is
verified progress, not whole-code closure. The present round continues at
35% distance with the same encoder and actual limited-support refresh.

The frozen exact helper closes582--614 at z105994380787/160000000000,
alpha119/1280,256-bit coefficient rounding, and384-bit Arb arithmetic.
The complete directed log2 union-bound interval is

    [-1151.0231463741918299304978853756903876143767888939239926070726567971974884960387300788745499269040948190015219972 +/- 4.41e-110].

The run takes173.0248511 seconds and peaks at330518528 bytes. Its local
coefficient hash is `3a9591447495de3082f1d2f8885734dc4a6e039774919009e38b19389449d399`.
All exact controls pass. The proof witness is an affine interpolation
between successful witnesses at512 and1024; the bound itself is computed
for every included size. Adding less than2^-1151 preserves the recorded
218*2^-50 covered-class budget. Continuous coverage now reaches614.

The33-size wrapper limit is not a mathematical hypothesis. Each exact
rectangle and its six doublings depend on the largest requested size and
marker pair, not on the number of targets extracted. A wider wrapper can
reuse those live tables, but must retain each actual size's coefficient,
denominator, source factor, run counts, and both terminal orders. Wider
coverage is admitted only after the complete sum passes its numerical gate.

While the wider wrapper is reviewed, the frozen helper closes615--647
with z207640114271/320000000000 and unchanged alpha/arithmetic:

    [-1172.7311927011228502592276227928264955431517166323429043798020384581013406413658628659911836944782510255084600030 +/- 5.68e-110].

It takes191.0660790 seconds and peaks at346181632 bytes. The local
coefficient hash is `050d336c28a6b99cb453883101891a60762dac9a3b005549330771cd3138c736`.
All exact controls pass. This extends continuous coverage to647 without
changing the covered-class218*2^-50 budget.

### Wider wrapper and exact controls

The new helper imports the frozen lower-exact helper, rather than editing
its historical limits. It admits at most129 actual sizes through1024.
It reuses the same live rectangles but performs separate complete joins;
there is no cross-size matrix envelope. It retains256-bit coefficient
rounding, the280000-position packed limit, and512-MiB peak-memory guard.
The300/600/1200-second guards are cooperative. Failure before the final
sum supplies no whole-band claim.

Independent source review checks the wrapper, own-target normalization,
once-only source factor, complete explicit-weight joins, output-only
terminal, and storage limits. Root's tests-only run passes1181 integer,
80 normalization,36 chain,409 ordinary-marker,310 inherited-wrapper,
and2397 new wide-index/rational checks, plus840 exhaustive masks.
The new checks cover218995 target indices, odd/even sizes, source-factor
jumps, and maximum-width schedules. No production table is built by this
test. It takes0.9973341 seconds and peaks at34652160 bytes.

| Source | Bytes | SHA256 |
|---|---:|---|
| `cheap_lane_lower_wide_band.py` |17104|`f10d185e8befc4f8bc1a0a961f0d518b2bb94b2ab317cc362adb5fca9b36a03a`|

### Complete wide trial and passing prefix

The648--776 pilot uses z984830809/1600000000 and alpha119/1280. All129
per-size calculations complete, but the whole-band directed interval is

    [1978.4852167995842554826377300428653702386674145559601169032649894626535609227916437837001919518233742393826035649 +/- 6.49e-110].

It takes473.2364901 seconds and peaks at473026560 bytes. The local rational
coefficient hash is `67645fba568b51e440913c178621c8ed7d13a68daaa894114f082857639dcfc2`.
All exact controls pass, and the run completes every requested marginal,
ordered profile, and source count. The positive aggregate does not certify
the whole band. It also does not exhibit a low-weight codeword.

Each of the104 completed contributions648--751 has log2 upper endpoint
strictly below-76. This comparison uses exact rational parsing of every
reported Arb interval, not binary floating-point selection. The largest
upper endpoint occurs at750:

    [-76.8755834139133772737504034494633570179489630715197652116947371189065941595198203097632525764062361281062524622 +/- 4.19e-110].

Consequently the prefix union is below104*2^-76<2^-69. Only that complete
prefix is admitted; the remaining sizes are not inferred from it. Together
with the two earlier bands, round41 adds170 previously open sizes and
extends continuous coverage through751. Adding2^-1151,2^-1172,and2^-69
to the preceding accounting preserves218*2^-50<2^-42 for all covered
classes, including the true-quiet bad-distance event. Six gaps remain.

The final bounds have a strong odd/even pattern: at767,768,769 the log2
bounds are respectively+411.514,+1326.589,+571.895. At776 the leading
run count is55. This motivates a matched rounding control, not a code
change. For r776 and c55, the two matching full-pair counts satisfy
a+b=721 and a,b<=388, hence333<=a,b<=388. A smaller exact rectangle
with full limit388 and singleton limit110 retains every required target.
After own-marker normalization, it bounds the same coefficients as the
older full388/singleton188 rectangle. The proposed control compares all56
ordered profiles, both terminals, and the same indexed run mass/source
factor. A control result would concern run55 only, not the entire size776.

For the larger gaps, source review identifies the existing counted-band
pipeline as the smaller extension. `one_table_many` already extracts fresh
targets for up to17 actual sizes through4096. The current driver fixes its
witness and wrapper to a neighborhood of2048; a new orchestrator can reuse
the positive tables, own-size normalizers, and count-aware join. This is a
proposal, not new coverage. Merely changing coordinates to t=2f+s does not
solve storage: retaining all legal coefficients through t=2048 still needs
1050625 entries per channel. Modulo256, this coordinate change preserves
the existing alias classes. Subtracting an upper bound on unwanted aliases
would be unsound; any such correction needs a valid lower bound instead.

Reproduce the completed run from the repository root:

    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 582 --last 614 --z 105994380787/160000000000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_exact_band.py --first 615 --last 647 --z 207640114271/320000000000 --alpha 119/1280 --bits 256 --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_wide_band.py --tests-only --seconds 300
    C:/Python314/python.exe -B research/workstreams/char2_field_spin/theory/cheap_lane_lower_wide_band.py --first 648 --last 776 --z 984830809/1600000000 --alpha 119/1280 --bits 256 --seconds 600

## Resource and progress discipline

This first goal round made no encoder change and ran no benchmark or
large numerical search. It added only proof text and a small exact-rational
checker. The checker used bounded in-memory arithmetic, `python -B`, and
stdout only, after verifying more than 20 GiB of local free space.
The previous 3.65-MB remote preservation bundle remains unchanged.

The second round likewise changed no encoder and ran no benchmark. Its
bounded proposal search took about 13.5 seconds; each independent receipt
replay took about one second. It added two small helpers, a 45,959-byte
selected certificate receipt, and proof/index text. No bulk output, build,
transport archive, or cache was generated. The local volume remains above
50 GiB free. The preexisting remote bundle and certified baseline are
unchanged.

The third round added two bounded proof helpers and proof/index text.
It generated no result dataset, build, archive, randomness tape, or cache.
The one-group and few-group directed checks were independently replayed
at 256 bits, as were the frozen quiet receipt and existing local checks.
The candidate encoder's SHA256 remains
`28bee7d5ec10cd7491e9f9ec69070f498ef8e5e9b4952bd4ca1b7e757cbc4969`.
The certified baseline, measurements, and remote preservation bundle are
unchanged. No benchmark, commit, push, promotion, or paper edit occurred.
At that checkpoint, the cheap-lane proof notes, helpers, and retained
quiet receipt totaled less than 190 KiB; no new generated dataset was saved.
C: still had more than 50 GiB free. Both newly selected helper hashes were
checked after replay.

The fourth round adds two small proof helpers and updates existing proof
notes. No encoder, benchmark, build, generated dataset, archive, or cache
was added. The main five-group replay took about eleven seconds; the
dense-band replay uses one fixed scalar expression. The source and proof
review were independent of the numerical replay. C: remains above
50 GiB free, and all previously pinned evidence remains unchanged.
The selected cheap-lane proof files, helpers, and quiet receipt occupy
less than 230 KiB in total. The two added helpers total 25,931 bytes;
all numerical output was streamed to stdout.

The fifth round adds one 24,041-byte proof helper and compact proof/index
text. It generated no dataset, build, archive, randomness tape, or cache.
All numerical checks used bounded memory and stdout only. No encoder
change or benchmark occurred. The certified baseline and all earlier
hash-pinned proof evidence remain unchanged. The selected cheap-lane
notes, helpers, and quiet receipt remain below 270 KiB in total; C: has
more than 50 GiB free. No commit, push, promotion, or paper edit occurred.

The sixth round adds two source helpers totaling 21,837 bytes and compact
proof/index text. It creates no generated dataset, build, archive, tape,
or cache. All checks use bounded in-memory arithmetic and stdout.
The dense replay takes less than 0.1 seconds; the exact local tests and
selected screens take less than two seconds per invocation. The failed
cutoff diagnostics each take less than four seconds. No encoder change,
benchmark, commit, push, promotion, or paper edit occurred. The baseline,
candidate encoder, remote preservation bundle, and all earlier frozen
proof helpers remain unchanged.
The selected cheap-lane notes, sources, and quiet receipt total less than
300 KiB after this round. C: retains more than 50 GiB free. All selected
source hashes were checked after the main replays; generated-data growth
is zero.

The seventh round adds one 13,798-byte diagnostic source and compact
proof/index updates. Its arrays use less than 48 MB; selected replays
take about 2.5--4.7 seconds. The scalar follow-ups and exact local checks
use bounded in-memory arithmetic and stdout only. No encoder change,
benchmark, build, dataset, archive, cache, commit, push, promotion, or paper
edit occurred. All eleven selected helper/receipt hashes and the unchanged
candidate encoder hash passed the final check. The selected cheap-lane
sources, notes, and retained receipt occupy less than 325 KiB in total;
generated-data growth remains zero. C: has more than 50 GiB free.

The eighth round adds one 20,454-byte diagnostic helper and compact proof
and index text. Its source arrays occupy 47,815,084 bytes; the three main
replays each took about 0.9--1.1 seconds. Exact local tests used in-memory
rational arithmetic. No encoder, benchmark, build, dataset, archive,
randomness tape, cache, commit, push, promotion, or paper edit occurred.
Generated-data growth is zero. The certified baseline and earlier frozen
proof evidence remain unchanged.
All thirteen selected helper/receipt/encoder hash checks passed. The
cheap-lane sources, notes, and retained receipt occupy less than 352 KiB;
C: retains more than 50 GiB free.

The ninth round adds three source helpers totaling 39,092 bytes and
bounded proof/index updates. Positive coefficient arrays occupy about
2 MiB, with scratch bounded below 256 MiB. Individual selected runs take
less than ten seconds; exact support tests take less than 0.01 seconds.
No generated data, build, cache, archive, tape, or benchmark was added.
There was no encoder change, promotion, commit, push, or paper edit.
All earlier pinned evidence and the remote preservation bundle remain
unchanged. The local volume retains about 49.9 GiB free, above the 20-GiB
floor. Generated-data growth remains zero; selected sources, notes, and
the retained quiet receipt remain below 420 KiB.

The tenth round adds two source helpers totaling 23,431 bytes and compact
proof/index updates. Numerical invocations took less than ten seconds
each, ran serially, and generated no files or caches. There was no
encoder change, benchmark, build, archive, tape, promotion, commit, push,
or paper edit. All eighteen selected helper/receipt/encoder hash checks
passed. The selected cheap-lane notes, sources, and retained receipt
occupy less than 465 KiB in total; generated-data growth remains zero.
C: retains about 49.9 GiB free. The old certificate and remote preservation
bundle remain unchanged.

The eleventh round adds three source helpers totaling 28,903 bytes and
bounded proof/index updates. All numerical runs were serial, took less
than eighteen seconds each, and generated no files or caches. Positive
coefficient arithmetic used bounded in-memory scratch. The pair operator,
tile normalization, and input-weight lemma passed independent proof and
source review; the main agent replayed the selected pair, tile4, and
two-mark results. All twenty-one selected helper/receipt/encoder hash
checks passed. Selected cheap-lane sources, notes, and the retained quiet
receipt remain below 500 KiB; generated-data growth is zero. C: retains
about 49.9 GiB free. No encoder change, benchmark, build, archive,
promotion, commit, push, or paper edit occurred. The baseline certificate
and remote preservation bundle remain unchanged. The goal stays active;
the new numerical gains do not certify another source-support class.

The twelfth round adds five small helpers and one compact proof note.
All numerical jobs ran serially, used bounded in-memory arithmetic, and
generated no files or caches. The selected full-state decomposition and
surplus replays took less than sixteen seconds each; each new window
replay took less than one second. The new exact class calculations do
not require a benchmark, build, or encoder change. No archive, tape,
promotion, commit, push, or paper edit occurred. Earlier certificates,
the candidate encoder, and remote preservation bundle remain unchanged.
All 26 selected helper/receipt/encoder hash checks passed. The 30 selected
cheap-lane sources, notes, and retained receipt occupy less than 560 KiB;
C: retains about 49.9 GiB free. Generated-data growth is zero. The goal
remains active: the new restricted-class certificates do not supply the
missing whole-code cover.

The thirteenth round adds six compact helpers totaling 54,460 bytes and
updates the existing proof/index notes. The exact interval-tail replay
took less than nine seconds; the arbitrary-profile directed replay took
less than four seconds. All jobs ran serially with the 20-GiB preflight,
bounded in-memory arithmetic, bytecode disabled, and stdout-only results.
Generated-data growth is zero. All 32 selected helper/receipt/encoder
hash checks passed. The 36 selected cheap-lane sources, notes, and retained
receipt occupy less than 625 KiB; C: retains about 49.9 GiB free. No encoder
change, benchmark, build, archive, promotion, commit, push, or paper edit
occurred. The 3.8-ms baseline and candidate implementation are unchanged.
The goal remains active: the new 76-bit result closes a restricted path
class, not the whole-code distance claim.

The fourteenth round adds three bounded diagnostic helpers totaling
45,579 bytes and updates existing notes/index text. Main replays took
less than 14 seconds each; every numerical job ran serially, with the
20-GiB preflight and bytecode disabled. Generated-data growth is zero.
All 31 earlier proof/receipt hashes and the unchanged encoder hash passed;
the three new helper hashes were also verified. The 39 selected cheap-lane
sources, notes, and retained receipt occupy less than 700 KiB. C: retains
about 49.9 GiB free. No encoder change, benchmark, build, archive,
promotion, commit, push, or paper edit occurred. No numerical job remains
running. The complementary proof branches remain open; the next round
should isolate degree-wise route counting before further parameter scans.

The fifteenth round adds five diagnostic helpers totaling 58009 bytes
and bounded proof/index text. Every numerical job ran serially with the
20-GiB preflight and bytecode disabled. The longest selected run took
56.15 seconds; peak working set stayed below 95 MiB. All 39 selected
proof/receipt/helper hashes and the unchanged candidate encoder hash
passed. The selected cheap-lane sources, notes, and retained receipt
remain below 800 KiB; C: retains more than 49 GiB free. Generated-data
growth is zero. No encoder change, benchmark, build, archive, promotion,
commit, push, or paper edit occurred. No numerical job remains running.
The old certificate and remote preservation bundle are unchanged.

The sixteenth round adds six bounded diagnostic helpers totaling 37151
bytes and compact proof/index updates. All numerical jobs ran serially
with bytecode disabled and a 20-GiB free-space guard. The longest was the
16.68-second quiet root solve; the new route-first checks took less than
five seconds each. The main agent also replayed the weighted-Jensen test
in 1.57 seconds. All 39 historical evidence hashes, six new helper hashes,
and the unchanged candidate encoder hash passed. Generated-data growth
remains zero. No encoder change,
benchmark, build, archive, promotion, commit, push, or paper edit occurred.
The selected sources, notes, and retained receipt remain below 850 KiB;
C: retains more than 49 GiB free. No numerical job remains running.

The seventeenth round adds eight bounded proof helpers totaling 71607 bytes
and compact proof/index updates. Every numerical job ran serially, with
bytecode disabled, a 20-GiB free-space guard, and stdout-only output. The
longest run took 27.91 seconds; the live-segment test measured a 95,035,392-byte
peak, and the smaller cone jobs budget 128 MiB. Generated-data growth is
zero. Historical table pins, all eight new helper pins, and the unchanged
encoder hash pass. Selected sources, notes, and the retained receipt
remain below 925 KiB; C: retains more than 49 GiB free. No encoder change,
benchmark, build, archive, tape, promotion, commit, push, or paper edit
occurred. No numerical job remains running. The active goal is not complete.

The eighteenth round adds nine bounded helpers totaling77909 bytes and
compact proof/index updates. The longest run took162.13 seconds inside
its180-second cooperative budget. The full-cover and weight-pilot runs
had384-MiB and512-MiB memory budgets. Sampled process records showed
192745472- and211496960-byte peaks at those times, not measurements of
their final peaks. All numerical jobs ran serially with
bytecode disabled, a20-GiB free-space preflight, and stdout-only results.
Generated-data growth is zero. All42 hash/size pins tabulated in this
ledger pass, along with the unchanged encoder hash. The selected cheap-lane
sources and notes occupy965005 bytes before this final resource edit;
C: retains more than49 GiB free. No numerical job remains running.
No encoder change, benchmark, build, archive, tape,
promotion, commit, push, or paper edit occurred. The full-source-size
certificate is meaningful progress, but the whole-code goal remains active.

The nineteenth round adds three small sources totaling27643 bytes and
compact proof/index updates. The longest computation took158.65seconds
within its180-second budget; all jobs ran serially with512-MiB maximum
declared scratch, bytecode disabled, stdout-only results, and the20-GiB
free-space guard. The contiguous band took63.78seconds. The two receipt
adjustments each took less than0.1seconds. Generated-data growth is zero.
All45 tabulated source hash/size pins and the unchanged encoder hash
pass; C: retains more than49 GiB free. No numerical job remains running.
There was no encoder change, benchmark, build, archive, tape, promotion,
commit, push, or paper edit. The actual support coverage and joint
missing-group lemma are progress; the whole-code goal remains active.

The twentieth round adds five bounded sources totaling54936 bytes and
proof/index updates. All numerical jobs ran serially with bytecode disabled,
a20-GiB free-space guard, and stdout-only results. The longest new run
took81.11 seconds; coefficient jobs had a512-MiB estimated scratch budget.
The calibrated-profile replays took5.32 and3.39 seconds, with measured
process peaks below56 MiB. Generated-data growth is zero. There was no
encoder change, benchmark, build, archive, tape, promotion, commit, push,
or paper edit. All50 tabulated source hash/size pins and the unchanged
encoder hash pass. Cheap-lane sources and notes occupy about1.1 MB;
generated-data growth remains zero, and C: retains
more than49 GiB free. No numerical job remains running.
Continuous support coverage, the dense extension, and
the calibrated512-profile certificate are progress; the whole-code goal
remains active.

The twenty-first round adds four bounded sources totaling50871 bytes and
updates existing proof/index text. The longest numerical run took108.11
seconds; the new calibration and residue checks took at most4.16 seconds
each, apart from the48.12-second complete FLOAT profile screen. All jobs
ran serially with bytecode disabled, the20-GiB free-space guard, and
stdout-only results. New residue runs measured less than46 MB process
peak; the selected calibration measured55,435,264 bytes. The source-size
band helper retained its512-MiB declared scratch budget. No generated
dataset, build, archive, tape, or cache was added. All54 tabulated
hash/size pins and the unchanged encoder hash pass. There was no encoder
change, benchmark, promotion, commit, push, or paper edit. The volume
retains more than49 GiB free, and generated-data growth is zero. Selected
cheap-lane sources and notes occupy about1.2 MB in total.
No numerical job remains running. The continuous cover through224,
measured coefficient gains, and alpha-only obstruction are progress;
the whole-code goal remains active.

The twenty-second round adds five bounded helpers totaling 66,292 bytes
and compact proof/index updates. All numerical jobs ran serially with
bytecode disabled, a 20-GiB free-space preflight, and stdout-only output.
The successful directed whole-size check took 411.89 seconds and measured
a 273,104,896-byte process peak. An earlier 16-table directed attempt was
intentionally stopped after two completed cells (162.15 seconds cumulative)
because its projected runtime exceeded the budget; it produced no final
bound or certificate. Generated-data growth is zero. All 59 tabulated
source hash/size pins and the unchanged encoder hash pass. The 84 selected
cheap-lane sources and notes occupy about 1.28 MB; C: retains more than
49 GiB free. No numerical job remains running. There was no encoder change,
benchmark, build, archive, tape, promotion, commit, push, or paper edit.
The complete size-512 certificate is progress, not whole-code closure;
the goal remains active.

The twenty-third round adds three helpers totaling40,762 bytes and
updates existing proof/index notes. All numerical jobs run serially,
with bytecode disabled, a20-GiB free-space preflight, and stdout-only
results. The longest run takes834.60 seconds with a280,961,024-byte
process peak. The directed subtotal calculation takes less than one
second; the peeling diagnostic takes49.08 seconds. Generated-data growth
is zero. No encoder change, benchmark, build, archive, tape, promotion,
commit, push, or paper edit occurs. The baseline and frozen proof helpers
remain unchanged. The508--515 joint certificate and identified counting
loss are progress. The whole-code goal remains active; the next tests
should address broad source-size intervals.
All62 tabulated source hash/size pins and the unchanged encoder hash
pass. The87 selected sources and notes occupy about1.33 MB; C: retains
about49.96 GiB free. No numerical job remains running.

The twenty-fourth round adds six bounded helpers totaling109481 bytes and
updates proof/index text. Numerical jobs run serially, with bytecode
disabled,20-GiB free-space preflights, and stdout-only results. All new
selected runs take less than nine seconds. The largest measured process
peak is120,074,240 bytes. Generated-data growth remains zero: no dataset,
cache, build, archive, or tape is added. The stronger analytic obstruction,
exact boundary result, and genuine return-cap gain are progress; the
whole-code goal remains active. No encoder change, benchmark, promotion,
commit, push, or paper edit occurs. Historical helpers remain frozen.
All68 tabulated source hash/size pins and the unchanged encoder hash pass.
The93 selected sources and notes occupy about1.46 MB; C: retains about
49.97 GiB free. No numerical job remains running. The next binary-pair
pooling test is specified in the rank interface but has not been run.

The twenty-fifth round adds seven bounded helpers totaling102346 bytes
and compact proof/index updates. All numerical jobs run serially, with
bytecode disabled, a20-GiB free-space preflight, and stdout-only output.
The longest run takes57.14 seconds; the largest measured process peak
is125497344 bytes. Generated-data growth remains zero. The boundary
calculation uses directed arithmetic, but the complete-coverage result
is explicitly FLOAT-only. The genuine local improvements and complete
size1024 diagnostic are progress; the whole-code goal remains active.
There is no encoder change, benchmark, build, cache, tape, archive,
promotion, commit, push, or paper edit. Historical helpers remain frozen.
All75 tabulated source hash/size pins and the unchanged encoder hash pass.
The100 selected sources and notes occupy about1.58 MB in total; C: retains
about48.72 GiB free. No numerical job remains running. The CWC-guided
notes separate proved local lemmas, FLOAT coverage, and certificate gaps.

The twenty-sixth round adds eight bounded helpers totaling95990 bytes
and updates existing proof/index notes. All numerical jobs run serially,
with bytecode disabled, a20-GiB free-space preflight, and stdout-only output.
Two Arb pilots stop at their explicit memory or runtime limits and yield
no certificate. The replacement exact-integer full replay takes154.506 s
and peaks at230801408 bytes; its selected precision control takes66.205 s
and peaks at287518720 bytes. Generated-data growth remains zero.
All83 tabulated source hash/size pins and the unchanged encoder hash pass.
The109 selected sources, notes, and retained receipt occupy less than1.71 MB;
C: retains about48.69 GiB free. No numerical job remains running. There is
no encoder change, benchmark, build, cache, tape, archive, promotion,
commit, push, or paper edit. The entire1024-group class now has a directed
430-bit bound. This is progress, not whole-code closure; the goal remains
active and the next work is broad source-size coverage.

The twenty-seventh round adds eleven bounded sources totaling103006 bytes
and updates existing proof/index notes. Every numerical job runs serially,
with bytecode disabled, a20-GiB free-space preflight, and stdout-only output.
The new complete bucket screen takes26.709 s; the matched bucket control
peaks at126640128 bytes. All94 tabulated source hash/size pins and the
unchanged encoder hash pass. The120 selected sources, notes, and retained
receipt occupy less than1.84 MB; C: retains about48.24 GiB free.
Generated-data growth is zero. No numerical job remains running.
There is no encoder change, benchmark, build, cache, tape, archive,
promotion, commit, push, or paper edit. CWC-guided notes separate the new
local lemmas and negative selected FLOAT bound from unchanged certificate
coverage. Those proof improvements are progress, but the complete2048
screen remains positive. The goal stays active; the next experiment is
the calibrated(898,898) control, followed by regional-bank repairs.

The twenty-eighth round adds nine bounded sources totaling86858 bytes
and updates proof/index notes. Eight sources have numerical replay and
frozen pins; the weighted-deficit prototype remains unrun and optional.
Every numerical job runs serially, with bytecode disabled, a20-GiB
free-space preflight, and stdout-only output. The longest completed run
takes45.09 seconds; the maximum measured process peak is128692224 bytes.
One initial toy test stops before any real cap construction and is fixed
as recorded above. Generated-data growth remains zero. There is no encoder
change, benchmark, build, cache, tape, archive, promotion, commit, push,
or paper edit. The129 selected sources and notes occupy less than1.95 MB;
C: retains about48.22 GiB free. No numerical job remains running.

All102 tabulated source size/hash pins and the unchanged encoder hash pass.
CWC-guided proof notes distinguish audited local lemmas, complete FLOAT
coverage, and directed certificates. The new86.7-bit complete size2048
candidate is progress, not whole-code closure. The goal remains active.
Next: use the frozen z23/100,alpha23/256 markers and ordered pair operator
for bounded directed residue256 verification before expanding source sizes.

The twenty-ninth round adds two bounded helpers totaling31092 bytes and
updates existing proof/index notes. All three numerical jobs run serially,
with bytecode disabled, a20-GiB free-space preflight, and stdout-only
output. The selected controls take72.693 and68.038 seconds; the complete
directed replay takes848.767 seconds. The maximum process peak across
the three runs is333684736 bytes. Generated-data growth is zero.
All104 tabulated source size/hash pins and the unchanged encoder hash
pass. The131 selected sources and notes occupy less than2.00 MB;
C: retains about48.21GiB free. No numerical job remains running.
There is no encoder change, benchmark, build, cache, tape, archive,
promotion, commit, push, or paper edit. The CWC-guided ledger distinguishes
the new complete size2048 directed certificate from the untested band
envelope and the remaining whole-code gaps. This is progress:86.7425 bits
at35% for the entire size2048 class. The goal remains active; next test
shared-table coverage of actual sizes2047--2049.

The thirtieth round adds five bounded sources totaling63851 bytes and
updates existing proof/index notes. The three numerical runs are serial,
with bytecode disabled,20-GiB free-space preflights, and stdout-only
results. They take68.687,856.263,and39.871 seconds; the largest process
peak is279838720 bytes. Generated-data growth remains zero. All109
tabulated source hash/size pins and the unchanged encoder hash pass.
The136 selected sources and notes occupy less than2.06 MB; C: retains
about47.90GiB free. No numerical process remains running. No encoder
edit, benchmark, build, cache, archive, tape, promotion, commit, push,
or paper edit occurs. CWC-guided notes distinguish the complete
74.2104-bit neighboring-band certificate from the selected wider-band
improvement. This is progress, not whole-code closure; the goal stays
active. Next run complete count-aware coverage of2040--2056, then test
the current operator at4096 before changing the construction.

The thirty-first round adds two bounded sources totaling30056 bytes and
updates existing proof/index notes. Four numerical jobs run serially:
the complete17-size directed replay, two complete size4096 FLOAT screens,
and one directed selected-profile control. They take981.428,35.201,
49.513,and37.066 seconds. The largest process peak is293822464 bytes.
Every run uses bytecode-disabled, stdout-only execution and a20-GiB
free-space preflight. Generated-data growth is zero. All111 tabulated
source size/hash pins and the unchanged encoder hash pass. The138
selected sources and notes occupy less than2.10 MB; C: retains about
47.91GiB free. No numerical process remains running. There is no
encoder edit, benchmark, build, cache, archive, tape, promotion, commit,
push, or paper edit. CWC-guided notes distinguish the60.2491-bit joint
band certificate, complete FLOAT larger-support evidence, and selected
directed controls. This is progress, not whole-code closure. The active
goal continues with complete directed verification at actual size4096.

The thirty-second round adds six bounded sources totaling76547 bytes and
compact proof/index updates. Seven numerical jobs run serially, with
bytecode disabled,20-GiB free-space preflights, and stdout-only output:
the full4096 directed replay, ordinary4096 control,6144 FLOAT screen,
ordinary6144 directed replay, complete ordinary-profile diagnosis,
exact endpoint tests, and complete all-marker directed comparison.
They take941.572,94.230,30.536,32.328,29.559,0.049,and85.550 seconds.
The maximum process peak is283136000 bytes. Generated-data growth is
zero. All117 tabulated source size/hash pins and the unchanged encoder
hash pass. The144 selected sources and notes occupy less than2.21 MB;
C: retains about47.87 GiB free. No numerical process remains running.
No encoder edit, benchmark, build, cache, tape, archive, promotion,
commit, push, or paper edit occurs. CWC-guided notes separate the new
complete4096/6144 certificates, FLOAT diagnostics, the corrected cause
of the ordinary bound's loss, and the untested band extension.
The goal remains active. Next test a genuinely broad ordinary hole band,
then bridge the smaller support-size gaps; do not infer those intervals
from the large singleton margins.

The thirty-third round adds four bounded sources totaling49973 bytes and
updates compact proof/index text. Seven complete band replays and one
small directed-Gamma API control run serially. The band replays take
14.515,14.281,13.343,13.207,16.441,17.009,and17.550 seconds; the first
console receipt is truncated and is replaced by the full second receipt.
The largest process peak is103710720 bytes. All runs use bytecode-disabled,
stdout-only execution after20-GiB free-space preflights. Generated-data
growth is zero. All121 tabulated source size/hash pins and the unchanged
encoder hash pass. The148 selected files occupy less than2.26 MB;
C: retains about47.74GiB free. No numerical process remains running.
No encoder change, benchmark, build, cache, tape, archive, promotion,
commit, push, or paper edit occurs. CWC-guided notes separate the new
897-size joint certificate from failed fixed-bank extensions and the
remaining whole-code gaps. This is progress; the goal remains active.
Next fit and verify a new normalized witness bank around actual5120.

The thirty-fourth round adds two bounded sources totaling 28850 bytes and
updates existing proof/index notes. Eight directed replays run serially,
with bytecode disabled, 20-GiB free-space preflights, and stdout-only output.
They take 18.439, 18.321, 17.859, 16.845, 16.721, 17.351, 16.086, and 18.642 seconds.
The maximum measured process peak is 103305216 bytes. Generated-data growth
is zero: no dataset, cache, tape, build, or archive is added. No encoder
change, benchmark, promotion, commit, push, or paper edit occurs. The new
continuous 4864--8192 bound is progress, not whole-code closure. The goal
remains active; next narrow the 4608-bank test before bridging the 4096 edge.
All 123 tabulated source hash/size pins and the unchanged encoder hash
pass. The 150 selected sources and notes occupy less than 2.30 MB; C:
retains more than 50 GiB free. No numerical process remains running.

The thirty-fifth round adds three bounded sources totaling58412 bytes and
updates existing proof/index notes. Thirteen directed replays run serially,
with bytecode disabled,20-GiB free-space preflights, and stdout-only output.
They take17.913,17.578,17.885,18.131,17.677,18.178,16.463,16.049,
16.276,16.658,20.914,16.931,and18.332 seconds. The maximum process peak
is102555648 bytes. Generated-data growth is zero. No encoder change,
benchmark, build, cache, tape, archive, promotion, commit, push, or paper
edit occurs. All126 tabulated source hash/size pins and the unchanged
encoder hash pass. The153 selected sources and notes occupy less than
2.45MB; C: retains more than50GiB free. No numerical process remains
running. CWC-guided notes separate new band certificates from positive
upper bounds and the proposed next refinement. This is progress, not
whole-code closure. The goal remains active; next refine the leading
3072-group cells by run count while retaining complete coverage.

The thirty-sixth round adds four bounded sources totaling66962 bytes and
updates the existing proof/index notes. Eight numerical jobs run serially,
with bytecode disabled,20-GiB free-space preflights, and stdout-only
results. They take16.535,3.924,20.800,61.472,61.377,66.054,76.118,
and875.753 seconds. The largest measured process peak is350720000 bytes.
Generated-data growth is zero. All130 tabulated source hash/size pins
pass, and the encoder hash is unchanged. The selected156 theory sources
and top-level cheap-lane notes, including the workstream README, occupy
less than2.50MB. C: retains more than50GiB free; no numerical process
remains running. No encoder edit, benchmark, build, cache, tape, archive,
promotion, commit, push, or paper edit occurs.

The CWC-guided notes distinguish the initially selected-profile control
from the subsequent complete3356.7491-bit size3072 certificate. They also
record that the requested256-bit Arb control actually ran at the frozen
384-bit minimum. Exact rational witnesses and source pins are retained;
raw console dumps are not added to the repository. This is new coverage,
not whole-code closure. The goal remains active. Next verify a complete
neighboring source-size band using the same regional tables and each
size's own source factor and normalizers.

The thirty-seventh round adds four bounded sources totaling77104 bytes
and updates existing proof/index notes. Nine numerical jobs run serially,
with bytecode disabled,20-GiB free-space preflights, and stdout-only
output. They take39.085,38.176,831.375,18.394,20.146,17.172,17.020,
15.656,and15.969 seconds. The largest process peak is291004416 bytes.
The only new selected generated data is the13754-byte regional-bound
witness; no bulk dataset, coefficient array, cache, tape, build, or
archive is retained. All134 tabulated source hash/size pins pass.
The encoder hash remains unchanged. The161 selected theory sources and
top-level cheap-lane notes, including the workstream README, occupy less
than2.63MB; the selected witness is counted separately. C: retains more
than50GiB free, and no numerical process remains running.

The CWC-guided record separates the audited reusable-numerator lemma,
verified generation, conditional cached replay, new joint band coverage,
and unsuccessful endpoint bounds. Independent final review confirms the
joint-margin calculation, exact remaining gaps, reproduction commands,
receipt length and both hashes, and the cached replay's trust boundary.
The continuous3008--3168 certificate
is progress, not whole-code closure. The goal remains active. Next
recalibrate proof witnesses near2944 and3264 before regenerating bounded
regional tables. No encoder change, benchmark, promotion, commit, push,
or paper edit occurs.

The thirty-eighth round adds three bounded sources totaling65057 bytes
and updates existing proof/index notes. Twelve numerical jobs run serially,
with bytecode disabled,20-GiB free-space preflights, and stdout-only output.
They take16.077,15.828,17.355,16.785,20.097,21.978,837.569,19.571,
36.067,19.398,17.620,and17.737 seconds. The largest process peak is
295571456 bytes. Only the13785-byte selected witness is retained as new
generated data. The two selected regional witnesses total27539 bytes;
no bulk data, coefficient arrays, caches, tapes, builds, or archives are added.

All137 tabulated source hash/size pins pass, and the encoder hash is
unchanged. The163 selected files comprise cheap_lane Python sources,
top-level CHEAP_LANE notes, and the workstream README. They occupy less
than2.67MB; selected receipts are counted separately. C: retains more
than50GiB free, and no Python
process remains running. The preserved baseline manifest and its two
verification receipts also pass their existing hash checks; this is
provenance verification, not a fresh baseline coefficient replay.

Independent source audits cover the three helpers, marker-selection
lemma, new witness metadata, and the proposed lower-domain extension.
Final independent review confirms the coverage, remaining gaps, joint
margin, commands, all137 source pins, and the unchanged encoder hash.
Its correction to the marker-trial counting unit is incorporated above.
The CWC-guided record separates verified generation, conditional cached
replay, complete new support bands, and unsuccessful upper bounds. This
round adds224 sizes at35%, extending continuous coverage to3008--3392
jointly with more than196 bits. The whole-code goal remains active.
Next test the lower-domain adapter, starting at225 and the source-factor
jumps near256; do not extrapolate the current bands. There is no encoder
change, benchmark, promotion, commit, push, or paper edit.

The thirty-ninth round adds two bounded helpers totaling46201 bytes and
updates existing proof/index notes. Fourteen numerical runs execute serially,
with bytecode disabled,20-GiB free-space preflights, and stdout-only output.
Their times are118.183,125.855,0.377,20.036,13.165,13.683,147.167,
23.771,24.510,69.445,39.500,119.286,136.635,and52.044 seconds. The
23.771-second report is truncated and replaced by its complete replay.
The maximum measured process peak among the new helpers is264601600 bytes.
The three old cone runs provide an estimated512-MiB figure only; no measured
peak is inferred for them. The new helpers' time/memory guards are cooperative,
not external hard limits.

Generated-data growth is zero. The two existing selected regional witnesses
still total27539 bytes. No coefficient arrays, bulk receipts, caches, tapes,
builds, archives, or benchmark data are retained. All139 tabulated source
hash/size pins pass. The encoder, preserved baseline manifest, and its two
verification receipts retain their hashes. The165 selected source/note files
occupy less than2.75MB; selected receipts are counted separately. C: retains
more than50GiB free. No numerical process remains running.

Independent reviews cover both helpers, odd/even exact extraction and joining,
the lower-domain source comparison, the old-versus-new516 factors, and the
joint small-source accounting. CWC-guided notes distinguish failed aggregate
witnesses from separately complete passing subsets, and isolate the matched
rounding refinement from local-operator improvements. New complete coverage
reaches1--312 and508--548, adding118 sizes at35%/40-bit target. The whole-code
goal remains active. Next extend the256-bit exact cover at313--345 before
bridging the remaining gaps. No encoder edit, benchmark, promotion, commit,
push, or paper edit occurs.

The fortieth round adds no source helper or generated-data file. Eight
directed runs execute serially after20-GiB free-space checks, with bytecode
disabled and stdout-only results. They take63.321,72.362,71.186,83.221,
96.997,105.100,117.267,and153.182 seconds. The maximum measured process
peak is282513408 bytes. The existing512-MiB and300-second helper guards
are cooperative, not external interruption limits. Every run passes the
exact integer, normalization, chain, and odd/even wrapper controls.

All139 tabulated helper byte/hash pins, the encoder, the preserved baseline
manifest, and both baseline verification receipts match independently
checked hashes. Generated-data growth is zero; the two retained regional
witnesses still total27539 bytes. The165 selected source/note files remain
below2.80MB. C: retains more than50GiB free. No numerical process remains
running. No encoder edit, benchmark, build, cache, archive, tape, promotion,
commit, push, or paper edit occurs.

The CWC-guided notes distinguish retuned proof witnesses from construction
changes and specify the bad-distance events in the joint failure budget.
Round40 adds228 previously open sizes, connecting coverage through581.
The covered-class ideal-setup failure bound remains below2^-42; six gaps
still prevent a whole-code certificate. This is verified progress and
the goal remains active. Next test the retuned582--614 witness, then
continue toward1024 without extrapolating completed intervals.

Keep this ledger updated with closed claims, failed approaches, and exact
remaining cases. Benchmark only one process at a time. Preserve source
and compact receipts locally; large numerical artifacts belong in bounded
remote storage. Do not mark the goal complete on the basis of local
lemmas, optimistic envelopes, or incomplete support coverage.
