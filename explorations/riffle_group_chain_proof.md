# Systematic packet-permuted group-chain Riffle: proof target

## 2026-08-10 correction: eight-bit permutation atoms, unchanged local codes

The active construction does **not** shrink either local code.  Both the
outer and recursive inner subcodes remain the committed binary extended BCH
`[128,64,22]` code.  The new parameter is only the granularity of the global
permutation between them.

Split every physical 64-bit outer group into eight contiguous 8-bit packets.
There are

```
M = N/8 = 262144 packets,
B = N/64 = 32768 recursive inner inputs.
```

Apply one uniform permutation to the `M` packets.  Concatenate each eight
consecutive permuted packets to obtain the next 64-bit input of the unchanged
recursive inner transform.  The internal order of a packet is initially left
unchanged.  Thus the construction preserves correlations among at most eight
binary coordinates, rather than among all 64 coordinates as in the earlier
whole-group permutation.

The random lane bijection already present in each physical outer group
factorizes into a random assignment of lanes to its eight packets and
independent random orders inside those packets.  This is preprocessing
randomness in the binary code construction; it is not a runtime Hadamard or
`8 x 8` mixer.  Consequently, conditional on its weight, each individual
packet has a uniform internal support, while packet weights originating in
the same outer group remain correlated.  The proof must retain that packet
weight correlation.

### Lemma 1 (transitivity on packet-weight profiles)

Fix a binary outer word and a packet width `g` dividing 64.  Independently in
each physical outer group, sample a uniform bijection from its 64 logical
lanes to its 64 physical lane positions.  Split the physical positions into
`64/g` ordered packets, and independently apply a uniform permutation to all
packets.  Let

```
a_j = number of final packets having Hamming weight j,  0 <= j <= g.
```

Conditional on `a=(a_0,...,a_g)`, the final ordered binary word is uniform on
the complete packet-weight profile class, whose size is

```
Q_g(a) = M! / product_j a_j! * product_j C(g,j)^a_j.
```

Indeed, a uniform 64-lane bijection factorizes uniquely into (i) an ordered
assignment of the 64 logical lanes to `64/g` unordered `g`-subsets and (ii)
an independent uniform order of the `g` lanes inside every assigned packet.
The second factor makes the concrete support of every packet, conditional on
its weight `j`, independent and uniform among the `C(g,j)` supports.  The
independent uniform packet permutation makes the ordered packet-weight
sequence uniform among its `M!/product_j a_j!` arrangements.  Multiplying the
two fiber sizes gives `Q_g(a)`.  Correlations among the packet *weights*
created by a common outer group are not discarded; the statement is only
conditional on their final counts.

This lemma is the precise justification for every division by `Q_g(a)` in
the grouped-permutation certificate.  A packet permutation by itself does
not imply it.

**Frozen implementation invariant.**  The encoder setup must sample either
the full independent 64-lane bijections above or the equivalent factorized
data (packet assignment plus independent within-packet orders), independently
of the global packet permutation and recursive state permutations.  These
bijections are preprocessing state and may be fused into the outer scatter;
they are not repeated work per encoded vector.  The experimental
`RifflePacketInnerChain` implements only the second-stage packet permutation
and deliberately preserves packet interiors, so it satisfies this lemma only
when fed by an outer layout implementing this invariant.

For a fixed outer word let `p` be its number of nonzero packets and let `G`
be the number of nonzero 64-bit recursive inputs after the packet
permutation.  Conditional on `p`, the occupied packet slots are a uniform
`p`-subset of the `M=8B` slots.  Therefore

```
Pr[G=g | p]
  = C(B,g) [x^p] ((1+x)^8-1)^g / C(8B,p),

E[q^G | p]
  = [x^p] (1 + q((1+x)^8-1))^B / C(8B,p).
```

These identities are exact and preserve all capacity-eight collisions.  The
ordered-slot recurrence in `scripts/certify_packet8_permutation.py` checks
the same law with integer arithmetic.

Every nonzero outer word activates at least one punctured local BCH block.
That block has at least 21 surviving coordinates, and the striped outer
layout places its distinct coordinates in distinct physical 64-bit groups.
They consequently lie in at least 21 distinct 8-bit packets, even in the
presence of every other active block.  If all nonzero packets lie in the last
`T=floor(d/64)=2949` recursive inputs, the output can have weight at most
`64T <= d`.  The probability of this former universal obstruction is now at
most

```
C(8T,21) / C(8B,21) < 2^-40.
```

The exact certificate reports the sharper exponent.  Thus the old
four-whole-group late-suffix counterexample is eliminated without weakening
either BCH impulse.

The first activation from zero state is also exactly enumerable.  For every
nonzero 8-bit packet value and each of its eight possible slots inside the
64-bit accumulator input, `scripts/analyze_packet8_impulse.py` computes

```
V = Acc(U),       S = P V.
```

Every row has `wt(V)+wt(S)>=22`, as required by the unchanged inner BCH
distance.  The smallest next-state weight over all 2040 nonzero packet/slot
choices is 18.  The table exposes the residual packet correlation honestly:
odd packet weight leaves a long accumulator suffix and has mean emitted
weight 32.5, whereas even packet weight returns the accumulator to zero after
the packet and emits only 1--7 bits.  The even rows nevertheless leave mean
next-state weight about 29--32, so the impulse is carried into subsequent
steps rather than lost.  This is the exact phenomenon that the sparse
episode ledger must track; no internal packet mixer has been assumed.

This placement lemma is not yet the full distance proof.  For dense outer
words, counting only `G` discards the weights of the nonzero packets and can
again be too weak.  The next operator must condition on the packet-weight
profile

```
(a_0,...,a_8),  where a_j = number of packets of weight j,
```

and combine eight sampled packet weights into each 64-bit inner input.  The
unchanged 65-state inner transfer can then be indexed by the resulting input
weight `u=0,...,64`.  This is the corrected high-support route: a
packet-profile transfer bound, not a smaller inner code and not the obsolete
tile-component cutoff at whole-group support 106.

There is also a stronger global dense-profile bound that avoids local inner
mixing assumptions.  Let `A_v` be the number of packets equal to the concrete
binary value `v in F_2^8`.  Conditional on these counts, the uniform packet
permutation has orbit size

```
M! / product_v A_v!.
```

For every fixed choice of the state permutations, the recursive systematic
inner is a binary bijection: from the emitted blocks one reconstructs each
drive and hence each input block causally.  Therefore no more than the binary
Hamming-ball volume `Vol(N,d)` members of a packet orbit can produce output
weight at most `d`.

Conditional on the packet-weight counts `a_j`, internal packet patterns are
independent and uniform among the `c_j=C(8,j)` values of weight `j`.  The
multinomial identity

```
E[product_v A_v! | (a_j)]
  = product_j a_j! C(a_j+c_j-1,c_j-1) / c_j^a_j
```

then gives an exact inverse-orbit moment.  At the balanced profile
`a_j=(M/256)C(8,j)`, `scripts/certify_packet8_dense_orbit.py` checks with
integers that the Hamming-ball bound, unioned over all `2^K` messages and
without using the graph codimension, is below `2^-40`.  This supplies the
dense endpoint.  The remaining global task is to set an orbit-entropy
threshold and prove that every profile below it is covered by the
packet-weight transfer/outer-spectrum branch.

A first robust weight-only transfer probe is implemented in
`scripts/probe_packet8_weight_transfer.py`.  It does not assume that the
accumulator output is uniform on a Hamming slice.  Instead it intersects the
exact state/input XOR shells with the exact accumulator weight rows and the
systematic BCH split caps.  This is a valid local relaxation, but it is far
too pessimistic in the dense regime: at the balanced iid diagnostic it gives
only about 1489 bits of low-output suppression, while unioning over the
message space needs roughly one million.  Hardening that relaxation would
therefore be wasted work.  The proof split is now:

1. exact packet placement and one-packet impulse/episode bounds for sparse
   and intermediate profiles;
2. the inverse-orbit/Hamming-ball bound for profiles with sufficiently high
   packet entropy;
3. an explicit overlap inequality showing that the two profile regions cover
   every outer word.

The first low-entropy regression is encouraging.  If the only concrete
packet values are `0x00` and `0xff`, with the conservative nonzero density
`21/128`, the
packet permutation alone has low orbit entropy.  This models the dangerous
family in which all outer lanes carry the same minimum-weight BCH word.
`scripts/probe_packet8_constant_value_scalar.py` computes the accumulator
output MGF exactly for all 65 incoming-state weights and all 256 subsets of
the eight current packet slots.  A two-state OFF/LIVE transfer retains the
first activation and upper-bounds every exact turnoff back to OFF; it makes
no assumption about the BCH next-state distribution.  At output pole `0.8`
the current floating Cauchy optimization gives approximately

```
low-output probability                     2^-65369
after a 64-dimensional family and r=24      2^-65329.
```

The local DP counts are exact; only the pole optimization and two-by-two
matrix arithmetic still need outward hardening.  This endpoint shows that
low packet-value entropy is not automatically a true-distance obstruction:
the accumulator plus randomized state carries the impulse across zero
packets.  It also provides the second endpoint needed for the eventual
profile-overlap argument.

## 2026-08-10 packet-profile goal checkpoint

The three endpoint certificates above were rerun without changes.  The exact
late-suffix exponent remains `-72.965477613186`, the one-packet impulse table
still has total impulse at least 22 in every row, and the balanced inverse
orbit endpoint still has `130330.690597` bits of margin beyond `2^-40`.

The repeated-value diagnostic is now exhaustive.  The byte-trellis in
`scripts/scan_packet8_repeated_values.py` computes the same 65-by-9 local
coefficient table as the original exact histogram DP while summing all 256
packet masks by active count.  Eight cross-check rows agree to relative error
at most `1.13e-15`.  Scanning all 255 nonzero bytes at density `21/128` and the
full pole grid finds the current worst repeated byte

```
v=0x18, wt(v)=2,
conditional low-output exponent             -42156.902488,
after the provisional 64-24 family charge   -42116.902488.
```

The next six rows are the other adjacent two-bit patterns `0x30`, `0x0c`,
`0x06`, `0x60`, `0xc0`, and `0x03`.  This is diagnostic floating arithmetic,
not an outer-family certificate.  The generated full-scan CSV has SHA-256
`9030e7ce8a36e443e463fc9cb4a980c5b4575d25a0c807ebc3097d1a14ffea1d`.

`scripts/probe_packet8_small_alphabet.py` generalizes the same OFF/LIVE
operator to an arbitrary bounded concrete-byte profile.  One fugacity is used
per byte, the global normalization is the exact multinomial orbit, and the
turnoff row is the exact coefficient of

```
(1 + sum_{v != 0} x_v t^wt(v))^8 / C(64,s).
```

With a single value it exactly reproduces the older scalar probe.  The equal
`0x18/0x30` two-value profile has diagnostic exponent `-42046.724702` after a
provisional 128-dimensional family and the 24-bit graph charge.  The
three-nonzero-value linear alphabet generated by those bytes is safer, at
about `-45229.489303` after the same provisional charge.  No small-alphabet
scan currently points toward an internal-packet invariant or a true-distance
failure.

The more useful generalization is conditional only on the nine packet-weight
counts.  `scripts/probe_packet8_weight_profile_scalar.py` sums every concrete
byte of each weight class and divides by the exact number

```
Q(a) = M! / product_j a_j! * product_j C(8,j)^a_j
```

of ordered binary packet strings with profile `a`.  It reproduces the
repeated `0xff` row and gives a balanced-profile full-`2^K` union diagnostic of
about `-133050.56`.  Multivariate optimization has several OFF/LIVE max-state
basins, so the script uses both the all-ones and empirical-distribution starts;
its optimized rows remain diagnostic.

There are now two exact high-profile tests in
`scripts/certify_packet8_profile_high_branches.py`.  For every supplied exact
profile it checks with integers:

1. the inverse-orbit/Hamming-ball bound already used at the balanced endpoint;
2. a full-message bijection moment at the exact pole `z=1/10`.

The second test uses only the fact that every frozen recursive inner map is a
binary bijection, hence the sum of `z^wt(Y)` over all binary inputs is exactly
`(1+z)^N`.  At the balanced profile its exact integer comparison reports
`-133156.250134`, while the orbit branch reports `-130370.690597`.

The current smallest explicitly inspected uncovered interpolation profile is

```
(33741,6963,24371,48742,60928,48743,24371,6963,7322).
```

It lies 85 percent of the way from the repeated `{0x00,0xff}` weight profile
to the balanced binomial profile.  The exact orbit and `z=1/10` full-family
tests are still above the target by respectively `8647.24` and `5919.57`
bits.  This is not evidence of construction failure: both tests pessimistically
union all `2^K` outer messages, while realizing that profile under the
three-band BCH layout has not yet been charged.  Along the same interpolation,
the full-family moment turns safe near 85.6 percent and the orbit branch near
85.8 percent, so the two high-profile mechanisms overlap numerically.

The active proof obligation is therefore narrower and more concrete than the
previous generic “profile overlap” label: bound the expected three-band outer
enumerator of packet-weight profiles below the exact high-profile threshold.
The independent within-band coordinate permutations and lane permutations
must be retained.  The exact BCH band projection spectra and the existing
cell-cap/Holder machinery are the intended inputs.  Evidence currently points
toward proof success or, at worst, a missing outer list-recovery bound; it does
not point toward a need for intra-packet randomization.

### Second checkpoint: global three-band enumerator and endpoint phase

The three-band outer enumerator now has a direct global Finner bound that is
substantially cleaner than a new component ledger.  For profile variables
`t_j`, define the exact random-lane group moment

```
R_w(t) = [y^w] (sum_j C(8,j) t_j y^j)^8 / C(64,w).
```

Each data-message block occurs in exactly one tile factor in each of the three
bands.  Applying Finner with exponent three to the independent uniform block
messages, and using the fact that every fixed-band projection is surjective,
gives

```
S_3(t) = 2^-64 sum_w C(64,w) R_w(t)^3,
A_outer(t) <= 2^K S_3(t)^(B/3).
```

This inequality is valid for the actual three-band incidence and every fixed
choice of within-band coordinate permutations; it does not replace the layout
by an independent configuration model.  Binary64 coefficient optimization is
implemented in `scripts/probe_packet8_three_band_profile_enumerator.py`.
The exact graph/puncture replacement has not yet been inserted, so the current
rows remain diagnostic.

On the repeated-to-balanced interpolation, combining this outer coefficient
with the conditional inner bound gives

```
balanced fraction 0.70     +24800 bits approximately (uncovered)
balanced fraction 0.75      -8426.95 bits
balanced fraction 0.80     -52129.64 bits
balanced fraction 0.85      safely covered (more than 9000 bits)
```

Thus the numerical intermediate branch now overlaps the exact high-profile
branches with thousands of bits to spare.  The remaining phase is genuinely
endpoint-heavy rather than a generic intermediate profile.

The exact endpoint core is linear.  If every packet is `0x00` or `0xff`, give
each packet one binary variable.  Using only the first two bands leaves
`680T` variables for `T` tiles.  The combined 85-coordinate EBCH projection
has rank 64 and 21 parity checks.  Substituting the actual independently
sampled coordinate permutations and balanced lane-to-packet maps gives an
exact binary constraint matrix.  `scripts/probe_packet8_constant_packet_rank.py`
finds nullity exactly one--the unavoidable global all-ones solution--for every
tested instance:

```
T = 1,2,4,8,16:  three seeds at each size
T = 32,64:       one seed at each size.
```

These are exact ranks of reduced random instances, not a full-size probability
certificate.  They strongly contradict the possibility of a persistent
packet-internal invariant.

The local code input to a rank/list-recovery proof is now exact.  The first-two
band projection has parameters `[85,64]` and dual dimension 21.
`scripts/analyze_ebch85_split_spectrum.py` enumerates its complete dual and
applies the bivariate MacWilliams transform.  The exact 42/43 split spectrum
has minimum pair weight 5, 1824 nonzero cells, and SHA-256

```
790894d6f7907d1c85dcbdec3e50984a02c7be80274a77df3d780175da4b9535.
```

For comparison, the ordinary two-edge-type configuration ensemble built from
this exact local spectrum has a strictly negative nontrivial weight exponent:
`-10.375` bits per block at density one half and about `-3.59` bits per block
even at the inspected one-percent endpoints.  This is diagnostic because the
real lift conditions every physical column to be a balanced eight-by-eight
lane partition.

The structured calculation has now been corrected to use the actual variable
degree.  Every constant-packet bit occurs in exactly eight block-row factors,
so Finner uses exponent eight, not the earlier 64-way Holder relaxation.  For
positive factors satisfying

```
p(a,b) <= f_0(a) f_1(b),
```

one band tile contributes

```
2^(8n) E[f(Bin(n,1/2))^8]^8.
```

`scripts/probe_packet8_constant_structured_upper.py` applies this to the
exact 42/43 split spectrum.  A uniformly repaired feasible binary64 witness
bounds the complete first-two-band expected solution count by about
`2^4096.001`, rather than `2^152832`.  With one common fugacity and the exact
zero assignment removed, it also bounds all nonzero projected solutions of
support at most `h`.  The resulting pointwise outer-plus-65-state diagnostics
are unsafe at `h=200` (`-31.52` bits) and safe at `h=208` (`-43.93` bits), then
improve rapidly.  Thus the generic constant-packet branch has a concrete
handoff at 208 packets; interval and outward certification remain pending.

The apparent failures at densities `63/64`, `1/128`, and `127/128` were
parameter/counting slack.  Minimizing the actual 65-state operator over its
pole grid changes the `63/64` universal-outer row from `+379.90` to
`-1686.31`.  Nearer endpoints require the cumulative support-dependent outer
bound, not an equal-band-density coefficient.  In particular, equal band
densities are useful diagnostics but do not cover arbitrary support splits.

The finite sparse residue has now been removed.  The exact-length construction
has 128 distinct data rows punctured once in band zero.  If a final constant-
packet word has at most 39 active packets, its active logical vertices--data
blocks plus the graph block when nonzero--number at most 14.  Any two active
eight-lane packets then intersect in at least two vertices.  Pair capacity one
forces all active packets into one physical band, while the exact fixed-band
outside distances `(5,4,5)` rule out every nonzero data block even after its
one possible band-zero puncture.  Thus
`scripts/certify_packet8_constant_min_support.py` proves the exact final
support floor 40.

The punctures and graph are also incorporated into the structured count.
From the exact 42/43 spectrum `A[a,b]`, a uniformly punctured band-zero row has
the exact 41/43 pair count

```
P[a,b] = (42-a) A[a,b] + (a+1) A[a+1,b].
```

`scripts/analyze_ebch84_punctured_band01_split_spectrum.py` validates its
integer mass, complement symmetry, and digest
`1e5e7bf5679098703b0d90c41e7718be1f1ceb4c9224ae7dcd7f976ea55d51e9`.
`scripts/probe_packet8_constant_punctured_structured_upper.py` uses 16256
ordinary factors and exactly 128 punctured factors under the same degree-8
Finner cover.  For each fixed nonzero data message and fixed final packet
assignment, the random `24 x K` graph map is uniform and at most one of its
`2^24` inputs supplies all 128 required replacement bits.  This supplies a
valid pointwise `-24` graph match charge for the constant-packet event; it is
not the formerly invalid bare graph charge for a general low-output MGF.

The OFF-start two-state renewal bound is much stronger than the stationary
65-state witness in the sparse regime.  At supports 21, 24, 32, and 40 its
diagnostic conditional exponents are approximately `-30.64`, `-35.70`,
`-50.59`, and `-65.61`.  Together with the exact support floor 40 and the
exact-length outer/graph count, the low constant endpoint has wide numerical
coverage rather than a `21..207` hole.

Two exact continuation inputs have been added.  The band-two-zero subcode has
dimension 21, minimum weight 22, 357 nonzero 42/43 split cells, and exact
split digest

```
d9f9a66f0bf8f16a543b8cdeff8ee33f4febeaf6ba733591cdb50dc80379252b.
```

The band-(1,2) projection is `[86,64]`, has dual dimension 22 and minimum pair
weight 5, and its exact split digest is

```
859217d417f433e1361aba51a40e05110b55d25c2a2a6473015754a76ff79b17.
```

A band-two-zero structured witness gives roughly `-18` to `-21` bits before
the graph/puncture correction.  A geometric-mean three-band hybrid pays an
unacceptable extra free-variable baseline (`2^6144` total), while a valid
determined-band-two weight tilt only improves the old `h=21` relaxed row to
`+5.70`.  These are retained as proof-method diagnostics, but the exact
support-floor and renewal argument supersede them for the endpoint ledger.

Finally, the all-ones endpoint itself is not an inner-distance obstruction.
`scripts/probe_packet8_repeated_value_state_transfer.py` retains all 65 state
weights, the exact accumulator histograms, and rigorous systematic-split caps.
For the fixed all-ones input sequence it gives a diagnostic conditional
low-output exponent `-681253.48` at pole `0.1`.  The same state refinement only
improves the density-`21/128` repeated `0xff` row from about `-65369` to
`-71016`, confirming that the remaining issue is the outer count, not a hidden
low-output trajectory.

The near-all-ones discontinuity is now closed exactly rather than by floating
continuity.  `scripts/probe_packet8_all_active_minplus.py` builds the 65-state
support graph from exact accumulator supports and rigorous BCH split caps.  An
integer potential verifies for every packet mask

```
emitted >= 6 + phi(next) - phi(current) - 3 * zero_packets.
```

The initial state minimizes `phi`, so `q` zero-packet defects leave at least
`196608-3q` emitted bits.  All `q <= 2621` are therefore deterministically
above 188743; at the endpoint the certificate reports 188745.  The
support-dependent 65-state/outer row at `q=2622` is already safe by thousands
of bits, so the deterministic and probabilistic regions overlap.

### One-conditioned-row outer bound for `g=4`

The `g=4` proof can retain one complete BCH row before applying the linear
Brascamp--Lieb inequality.  This removes most of the gap left by the earlier
row-independent outer bound.

Fix packet fugacities `t_0,...,t_4>0`, with `t_0=1`.  Define

```text
P_t(x) = sum_(j=0)^4 binom(4,j) t_j x^j,
R_w(t) = [x^w] P_t(x)^16 / binom(64,w).
```

The value `R_w(t)` is the packet moment conditional on a 64-bit column of
weight `w`.  Let `p in [1/2,1)` and set

```text
(p_0,p_1,p_2) = (1-p,p,p).
```

These coefficients satisfy the three-band linear BL dimension condition.
More generally, coefficients `p_0,p_1,p_2 in (0,1]` are valid when

```text
p_0+p_1 >= 1,  p_0+p_2 >= 1,  p_1+p_2 >= 1.
```

The first conditioned-row optimizer used only the symmetric boundary
`(1-p,p,p)`.  That restriction is not part of the lemma.  The current
optimizer also searches `p_0+p_1=1` with `p_2>=p_1`.
For `b in {0,1,2}` and `s in {0,1}`, define

```text
H_(b,s)(t) =
  (2^-63 sum_(u=0)^63 binom(63,u) R_(u+s)(t)^(1/p_b))^p_b.
```

Now fix one of the 64 BCH rows in a tile.  Its bit at one physical coordinate
is `s`.  The other 63 rows contribute `H_(b,s)(t)` at a coordinate in band
`b`.  The 63 free messages also contribute the factor `2^(63*64)` per tile.

Let `C` be the systematic extended BCH `[128,64,22]` code.  For `c in C`, let
`w_b(c)` be its weight in band `b`.  Put

```text
r_b = H_(b,1)(t) / H_(b,0)(t).
```

After conditioning on the selected BCH row, the remaining sum over that row
contains the three-band enumerator

```text
A_012(r_0,r_1,r_2)
  = sum_(c in C) product_(b=0)^2 r_b^w_b(c).
```

The repository has exact split spectra for bands `(0,1)` and `(1,2)`.  For
every `theta in [0,1]`, Cauchy--Schwarz gives

```text
A_012(r_0,r_1,r_2)
 <= sqrt(
      A_01(r_0^2,r_1^(2 theta))
      A_12(r_1^(2(1-theta)),r_2^2)
    ).
```

This inequality follows by factoring each summand into two nonnegative
terms.  The first term contains bands zero and one.  The second contains bands
one and two.  Their band-one exponents sum to the original exponent.

The 128 graph holes require a separate conditioned tile.  Each hole replaces
one band-zero coordinate in a distinct data row and a distinct tile.  In such
a tile, choose that punctured data row as the conditioned row.  The remaining
84 coordinates use the exact puncture-averaged `(0,1)` split spectrum,

```text
P_01(a,b) = (42-a) A_01(a,b) + (a+1) A_01(a+1,b).
```

The enumerator is divided by 42 because the punctured coordinate is uniform.
At the hole coordinate, the graph bit selects `H_(0,0)` or `H_(0,1)`.
Conditioned on a fixed nonzero data message, the random 24-bit graph syndrome
is uniform.  Its encoded graph word therefore has the exact stored graph-code
weight spectrum.  Averaging that spectrum jointly over all 128 distinct holes
retains their correlation; it does not charge 128 independent worst cases.

At `t_j=1`, all `H_(b,s)` equal one.  Each normal or hole tile then has total
mass `2^4096`.  The complete outer expression has mass `2^(256*4096)=2^K`.
This identity is both a normalization check and a required verifier invariant.

The discovery implementation is
`scripts/probe_packet_group_g4_conditioned_row_outer.py`.  The independent
outward implementation is in
`scripts/certify_packet_group_triangle_ledger.py`.  The latter reconstructs
the exact spectra, converts every frozen binary64 parameter to its exact
dyadic value, and evaluates all logarithms with directed intervals.

At the former leading profile

```text
(429359,24531,38395,28039,3964),
```

the optimized conditioned-row outer is approximately `267000.634654` bits.
With the frozen 20,000-step inner witness, the outward combined interval is

```text
[-186.3499452387, -186.3498302196].
```

The uniform per-profile target is about `-111.4150650164`.  Thus this fixed
witness closes that profile by at least `74.9347652031` bits.  This is a
single-profile result.  An end-to-end `g=4` claim still requires a complete
profile cover and an outward cell-local union ledger.

The asymmetric search is much stronger in the adjacent sparse region.  At
profile

```text
(443439,24075,32591,18739,5444),
```

it selects coefficients close to

```text
(p_0,p_1,p_2) = (0.478066,0.521934,0.999).
```

The independent outward implementation gives the combined interval

```text
[-76983.9807140, -76983.9805989].
```

The lower bound on the margin against the uniform target exceeds 76,872
bits.  The binary64 discovery value is about 2,292 bits smaller because
negligible positive terms underflow at the extreme tilt.  Only the outward
interval is theorem-facing.

The exponent-three obstruction has now been removed by using the linear
structure of the band projections rather than only the factor-incidence
degree.  At the rounded 70-percent repeated-to-balanced profile

```
(66458,5735,20070,40141,50176,40141,20070,5734,13619),
```

the old exponent-3 coefficient was `859864.827247` bits and left
`+35380.312595` bits after the exact full-message bijection moment.  The three
fixed-band projection kernels in one 64-bit EBCH message space have exact
dimensions `22,21,21`; their 64 joined basis vectors have rank 64.  Hence

```
ker(P_0) direct-sum ker(P_1) direct-sum ker(P_2) = F_2^64.
```

The same decomposition holds blockwise in the global message space.  If
`L_(b,t)` is the linear map exposed to tile `t` in band `b`, then for every
message subspace `V`,

```
sum_(b,t) dim L_(b,t)(V)
  >= sum_b dim P_b(V)
   = 3 dim(V) - sum_b dim(V intersect ker(P_b))
  >= 2 dim(V).
```

This is exactly the dimension condition for finite-field linear
Brascamp--Lieb with coefficient `1/2`.  Because each individual fixed-band
projection is surjective, the tile norms again reduce to iid column moments,
but now

```
S_2(t) = 2^-64 sum_w C(64,w) R_w(t)^2,
A_outer(t) <= 2^K S_2(t)^(B/2).
```

The direct-sum identity is checked with exact binary linear algebra by
`scripts/certify_ebch_three_band_kernel_decomposition.py`, with basis digest
`f2c0a6e227151fac4ecdadd640c81272c2fd6888ffedbf9c3adc7317d3acc1be`.
The binary64 coefficient optimizer
`scripts/probe_packet8_three_band_linear_bl2.py` gives outer coefficient
`783074.549655` at the 70-percent profile and combined exponent
`-41409.964997`.  On the same interpolation it crosses zero between 63.0 and
63.5 percent, reports `-907.590907` at 63.5 percent, and reaches
`-10783.492481` at 65 percent.

Two discarded diagnostics clarify why this algebraic step was needed.
Optimizing unequal fractional Finner cover weights changes the old 70-percent
outer row by only about 0.35 bits.  Exact pair moments from the `[85,64]` and
`[86,64]` split spectra show no one-shared-block tile correlation at binary64
resolution for powers one through three: the meaningful dependence is
genuinely three-way, exactly what the kernel decomposition resolves.

The symmetric point is not the full linear Brascamp--Lieb polytope.  Let
`K_b=ker(P_b)` and `r_b=dim(V intersect K_b)`.  Since the three kernels are a
direct sum, `sum_b r_b <= dim(V)`.  Therefore band coefficients
`p_0,p_1,p_2` satisfy the linear BL dimension condition whenever

```
sum_b p_b - 1 >= max_b p_b,
```

equivalently whenever every pair of coefficients sums to at least one.  The
one-parameter boundary family

```
(p_0,p_1,p_2) = (1-p,p,p),   1/2 <= p < 1,
```

is particularly useful because band zero has only 42 coordinates.  The
resulting outer evaluation is

```
2^K product_b
  (2^-64 sum_w C(64,w) R_w(t)^(1/p_b))^(256 n_b p_b).
```

`scripts/probe_packet8_three_band_linear_bl2.py
--optimize-band-coefficients` optimizes this valid family jointly with the
profile fugacities.

There was a second numerical issue in the low-entropy overlap: the 65-state
probe inherited fugacities optimized for the different OFF/LIVE objective.
`scripts/probe_packet8_weight_profile_state_tuned.py` directly minimizes the
robust kernel's own Cauchy objective.  At the 40-percent interpolation row this
changes the conditional exponent by about 58,306 bits, from
`-412581.131150` to `-471045.463052`.  At 50 percent the improvement is about
95,615 bits.  Exact finite iteration of the 30-percent tuned kernel changes
the Perron row by only about eight bits, so the improvement is genuinely the
fugacity saddle rather than a transient artifact.

Combining fixed tuned inner witnesses with optimized valid linear BL
coefficients gives the following binary64 grid:

```
balanced fraction     outer linear BL      inner transfer      combined
0.10                    28591.888265       -146606.267293     -118014.379029
0.20                   221262.820468       -231578.802177      -10315.981710
0.30                   352367.269866       -354645.760973       -2278.491106
0.40                   463414.396837       -471045.463052       -7631.066215
0.50                   571886.406714       -634336.009369      -62449.602655
0.60                   678351.863905       -692302.522810      -13950.658905
```

The optimized 5-percent outer coefficient is already `-90126.453793` before
the inner charge.  Thus the inspected repeated-to-balanced interpolation is
numerically covered from the exact constant-packet endpoint through the
high-profile region.  The active proof gate is now to replace that path/grid
with a finite cover of the complete nine-dimensional packet-weight simplex,
then insert the exact graph/puncture correction and outward-harden the fixed
witnesses.  No construction failure or need for intra-packet randomization
has been observed.

The first complete-simplex landscape probe is
`scripts/probe_packet8_profile_simplex_landscape.py`.  For a fixed outer and
inner witness its profile exponent is

```
constant - sum_j a_j log2(t_j f_j) - log2 Q(a),
```

which is convex in `a`.  This makes fixed-witness polytope certification
possible by checking vertices, but it also makes missing simplex anchors
immediately visible.  The initial interpolation-only atlas produced large
positive values at pure packet classes; these were witness-placement holes.
`scripts/probe_packet8_pure_class_frontier.py` shows that every nonzero pure
class is actually safe by hundreds of thousands of bits, even with a positive
inactive-class floor of `0.1`.  Adding pair midpoints and adaptive witnesses
continues to lower the sampled worst rows.  The remaining sampled failures are
not bad words and are not a certificate.

A tempting one-witness reduction couples the fugacities by
`C(8,j)t_jf_j=1`.  Then all profile vertices have the same exponent and the
complete profile sum costs at most

```
sum_(r=1)^9 C(9,r) (M-r+1)/(M (r-1)!) < 114.002.
```

`scripts/probe_packet8_joint_universal_witness.py` verifies the reduction, but
the resulting minimax witness is numerically unusable: the tested symmetric
row is about `+960112` bits after the profile-sum factor.  The useful version
is hierarchical rather than universal: sharp witnesses cover a fixed support
face, while positive inactive fugacities thicken that face into a neighborhood.

The adaptive scan also found the exact profile

```
(150014,48151,48529,15450,0,0,0,0,0).
```

Full robust tuning gives inner exponent `-334107.673362`.  The linear-BL outer
coefficient `406756.390284` alone misses by `72648.716921` bits.  This does not
survive the second existing rigorous outer bound: the exact EBCH total-spectrum
envelope in `scripts/probe_packet8_outer_profile_enumerator.py` gives
`251821.017981`, so their combination is about `-82286.655382`.  Both outer
families are now included in the landscape and cutting-plane atlas.  This row
is a concrete warning that the final ledger must take their pointwise minimum;
it is not a construction obstruction.  The next proof step is a
support-stratified finite atlas and convex-cell verification, not another flat
random sample or a packet mixer.

The paired cutting-plane experiment makes the stratification requirement
concrete.  For

```
(205799,0,0,32424,0,0,0,23921,0),
```

a deliberately broad inactive-class floor of `0.2` leaves `+133722.39` bits,
whereas retuning the same active support with a negligible inactive floor
gives `-15288.89` bits.  The next five selected sharp companions also close,
with margins between roughly 21586 and 251096 bits.  Thus a flat collection of
broad witnesses is inefficient, while the support-local rows remain strongly
safe.  `scripts/probe_packet8_adaptive_simplex_atlas.py --sharp-companion`
records both types and caches their fixed linear charges.  The next verifier
must cover support faces first, then thicken them by bounded rare-class count
buckets and certify each resulting convex cell at its vertices.

The finite geometry is now explicit rather than sampled.  Fixing an ordering
of the nine class counts gives one simplex chamber; its vertices are uniform
profiles on nested nonempty subsets of that ordering.  Hence all `9!` chambers
use only `2^9-1=511` distinct uniform-subset vertices.  When class zero is the
first vertex, the zero outer word is excluded and the chamber is clipped at
minimum outer weight 21, producing only 255 additional distinct vertices.
`scripts/probe_packet8_profile_ordered_chambers.py` exhausts all 362880
chambers by integer witness bitmasks.  The binary64 atlas now covers every
clipped vertex and all `511/511` uniform vertices, up from `251/511` before the
targeted paired atlas.  The final worst uniform-vertex value was
`-252.198771` against target `-168.700990`.  These are exhaustive finite-set
diagnostics, not outward certificates.

A chamber need not have one witness valid at all its vertices.  Along every
edge a fixed witness has a convex exponent, hence one interval of validity;
different witnesses may hand off.  The sample-independent binary64 probe
`scripts/probe_packet8_profile_edge_interval_cover.py` computes and unions
those intervals.  The greedy edge atlas closes all `2287/2287` relevant Hasse
edges.  This was not the full chamber one-skeleton: every pair of strictly
nested prefix subsets is an edge of some ordered simplex.  There are 18405
such comparable edges after excluding the zero-word vertex.  The complete
audit initially covered 16656; the first 20 long-chord anchors raised this to
16917, and the second 20-anchor pass raised it to `17077/18405`.  All targeted
midpoints closed, with useful cross-coverage, but 1328 chords remain.  The
validated monotone residual cache is
`out/packet8_comparable_edge_residual.json`.

The two-dimensional warning is more important.  There are 198580 distinct
strict nested-prefix triangle faces.  The exhaustive binary64 landscape in
`scripts/probe_packet8_nested_face_landscape.py` found only 17526 triangles
with one branch safe at all three vertices, although `193121/198580`
barycenters were safe under the union of branches.  The first four greedy
face anchors closed their selected barycenters and 124 neighboring rows.  The
fifth selected profile did not close:

```
face     100 < 1e0 < 1ee
profile  (0,12483,12483,12483,0,34328,34328,34329,121710)
target   -168.700990
broad    +70977.391179
sharp    +23165.930932
```

The sharp row improves the previous value by about 71000 bits but remains
roughly 23335 bits above target.  This is not evidence that the construction
fails: the optimizer is coarse, the bounds are upper bounds, and no actual
low-weight family has been produced.  It is the first concrete evidence that
the generic weight-profile atlas may need a new low-value-alphabet lemma or a
different proof decomposition.  Atlas growth is paused rather than treating
more anchors as progress by default.

The next diagnostic direction is deliberately pre-proof.  Adapt the expected
input-output enumerator composition from `enumerator_paper/turbo.tex` and the
polynomial transition machinery from `enumerator_paper/AccPoly.tex` to small
legal Riffle instances.  For a concrete packet-value profile `A`, the random
packet interleaver denominator is the exact orbit size

```
M! / product_v A_v!,
```

not the scalar bit-weight binomial.  The resulting small-length spectral
shape and finite-size trend should indicate whether intact packets incur a
real distance loss before further investment in the full outward proof.  The
matrix-accumulator derivation in `enumerator_paper/mtx.tex` is explicitly
marked broken and is not an input to this plan.

This remains substantial evidence that the frozen construction may be good,
but it no longer supports a claim that the current atlas method is plainly
converging to the full `2^-40` theorem.

### Full-size hard-profile product-kernel diagnostic

The first post-atlas experiment now targets the non-closing face in an exact
sampled full-size outer layout.  The packet IOWE remains a useful independent
direction, but its exact orbit state is the 256-component concrete-byte
profile `A_v`, with composition denominator

```
M! / product_v A_v!.
```

The nine packet-weight counts alone do not specify that orbit.  Also, the
distinct-tile graph-hole construction has no literal small-size copy, so a
small end-to-end enumerator would be a surrogate.

The more direct attack begins from the legal outer word whose 16384 data EBCH
blocks are all ones.  In the committed cyclic input basis the corresponding
local message is `0x8e63cf44efd4fa21`, not the all-one 64-bit message.  Choose
121710 non-hole packets and require a message delta to vanish on all their
coordinates.  They remain `0xff`, and the constraints decompose into exact
local EBCH kernels.  The graph syndrome is retained globally and its encoded
128 bits are placed in the sampled holes, so every search move is still a
legal outer codeword.

`scripts/probe_packet8_hard_profile_product_kernel.py` samples all three-band
coordinate bijections, balanced lane-to-packet maps, distinct graph holes,
the graph-coordinate bijection, and the random 24-by-K graph map.  Exact GF(2)
elimination found:

```
seed       product kernel dimension   median local nullity   rigid blocks
20260810   88097                      5                      2724
20260811   87992                      5                      2811
```

This is stronger than the initial 75k-dimensional heuristic.  On seed
20260810, exact graph-coupled greedy basis descent reduced profile L1 error
from the almost-all-`0xff` baseline to 83528 and stopped at

```
(0,14,274,2076,2577,27649,40893,43808,144853).
```

The target was not attained: low weights one through three remain deficient,
2577 forbidden weight-four packets remain, and weight eight is high by 23143.
This is a single-move local minimum, not evidence that the profile is
unattainable.  The next search should enumerate whole local kernels where the
nullity is small and use paired/tabu moves for larger kernels.  Only an exact
profile hit followed by the frozen recursive inner can distinguish an atlas
relaxation from a true low-distance family.  If that search remains
inconclusive, conditional parallel-tempered thermodynamic integration is the
best way to measure the inner-relaxation free-energy slack; the current face
needs only about 0.712 bits per recursive group of improvement.

The whole-kernel and thermodynamic follow-ups have now been run.  Exhaustive
local replacement through nullity eight improved the first seed to profile

```
(0,28,381,2612,2165,31773,41262,41835,142088),
```

with L1 error 73966.  A hundred thousand improving direction-pair proposals
reached L1 73840, while an independent layout and random anchor set reached
75106.  Nullity-ten replacement converged to the same neighborhood.  This is
a robust local-search basin, especially in the shortage of weights one
through three, but it cannot rule out correlated choices of the 121710
weight-eight packet locations.

`scripts/probe_packet8_hard_profile_thermodynamic.py` attacks the other side
of the decomposition.  Conditional on a finite-size rounding of the hard
nine-class profile, it samples packet order, concrete byte patterns within
each class, and the independent state-coordinate permutations.  Every move
recomputes the exact accumulator and systematic EBCH parity recurrence from
the earliest affected group.  Parallel tempering estimates the normalized
partition function from

```
log Z(beta) = - integral_0^beta E_beta[W] d beta
```

and reports `beta*D + log Z(beta)` as a diagnostic Chernoff exponent.  At
beta 0.6 the per-inner-group values were -19.5102, -20.1885, and -20.7665 for
16, 32, and 64 groups.  An independent 32-group seed gave -20.1518.  Extending
its burn-in produced -19.9891 on the full measurement window and -19.9330 on
the stabilized last half, with ten complete temperature-ladder round trips.

For comparison, `probe_packet8_outer_profile_enumerator.py` bounds the hard
outer profile by 538470.5231 bits, or 16.4328 bits per full-size recursive
group.  Reaching the -168.701 per-profile target requires an inner conditional
exponent below -16.4380 bits per group.  The least favorable stabilized
finite-size estimate is therefore about 3.50 bits per group stronger than the
required handoff.  This reverses the interpretation of the face warning:
the available evidence now points strongly to slack in the robust 65-state
inner envelope.  It still does not certify the face, because Monte Carlo
mixing, numerical quadrature, and finite-size extrapolation have not been
made rigorous and no exact legal outer word with the target profile has been
found.

Everything below records the earlier whole-64-group branch.  Its local BCH
certificates and three-band outer combinatorics remain useful inputs, but its
final `support <= 106` ledger does not certify the corrected packet-permuted
construction and must not be quoted as its final margin.

Status: current construction candidate for the performance/proof co-design
branch.  This is not yet a theorem claim.  It is a low-data-movement repair of
the 11.17 ms contiguous-layout kernel; the measured stripe delta is about
1.25 ms, before the systematic-circuit saving and removal of the parity pass.

## Parameters and outer code

At the primary checkpoint, `K=2^20`, `N=2^21`, and `b=64`.  Let
`B=N/b=32768`.  Split the user message into 16384 words
`m_j in F_2^64`.  Sample a uniform `24 x K` binary matrix `R`, form the graph
check `r=Rm`, and embed it as the 64-bit graph word `(r,0^40)`.  Encode the
16384 data words and this one graph word with the committed extended BCH
`[128,64,22]` encoder.  There is no componentwise-XOR parity block.

Before puncturing these 16385 words have length `N+128`.  The physical layout
has 256 tiles, each with 64 data-block lanes and 128 physical groups.  Sample
a uniform assignment of the data words to the 16384 lanes.  Choose 128 of the
256 tiles uniformly without replacement.  In each chosen tile, choose one
data lane uniformly and puncture a uniformly chosen coordinate from its
42-coordinate band-zero cell.  Map its other 127 coordinates bijectively to
the other 127 groups in that tile.  Map every unpunctured data word
bijectively to its tile's 128 groups.  Finally, map the graph codeword's 128
coordinates bijectively into the 128 freed lane positions.  Thus every
physical group still has exactly 64 coordinates and the final length is
exactly `N`.

The distinct-tile rule is proof-relevant but preprocessing-only.  It replaces
the invalid bare `2^-24` graph charge (which allowed adversarial overlap
between graph holes and data support) by the exact Maclaurin bound

```
2^-24 sum_w A_w [q + (1-q) h/(256*42)]^w
```

at total group support `h`.  At `h=106` its nonzero correction is below
`2^-48` for `q=.181` and below `2^-85` for `q=.05`.  Puncturing reduces the
minimum support of an active data block from 22 to 21, so the low-support
component-size cap is 323 rather than 308.  These claims are checked exactly
by `scripts/certify_three_band_exact_length.py` and add no runtime work.

The coordinate bijection for each data word, the graph-to-hole bijection, and
the lane order inside every physical group are sampled independently.  A
fixed data word of outer weight `w` therefore hits `w` distinct groups, and
conditional on a group's weight its occupied lanes are a uniform subset.
The global chain order is an independent uniform permutation of the `B`
groups.  Feistel is deliberately outside proof steps 1--6.

## Systematic-output recursive inner

Apply a fixed basis change to the committed cyclic EBCH generator so that

```
E_sys(v) = (v, P v),  v in F_2^64.
```

The original left and right halves are both invertible, so this changes only
the encoder basis and not the EBCH code.  Before each step independently
permute the incoming state's 64 coordinates.  The outer lane randomization
already supplies the independent permutation of the current input group.
Conditional on its weight, their XOR is therefore a uniform support.  Apply
the fixed length-64 accumulator `Acc` to that XOR; this is one
permutation-plus-accumulate weight-mixing round, with the permutation already
supplied by the support randomization.  Starting at `S_-1=0`, define

```
W_i       = U_pi(i) + sigma_i(S_(i-1))
V_i       = Acc(W_i)
Y_pi(i)   = V_i
S_i       = P V_i
```

for `i=0,...,B-1`.  The terminal state is discarded.  The final codeword is
the physical-order concatenation of the `Y_j`.

This branch has no per-step Singer map, repeated PAP mixer, or random
64-of-128 split.
Its inner randomness is the group-order permutation, the outer-supplied lane
permutations, and the local state permutations `sigma_i`.

## Exact transpose implemented by the performance probe

For the original 11.17 ms baseline and dual input groups `Z_j`, the adjoint
was evaluated from right to left using the cyclic split:

```
state = 0
for i = B-1,...,0:
    result = L^T Z_pi(i) + R^T state
    Z_pi(i) = result
    state = result
```

The generated split-input circuit computes this `E^T` call in 913 field
XORs.  The systematic target circuit uses 803 XORs for `P^T`, 64 state/input
merges, and 63 accumulator XORs, for 930 XORs per group.  This is only 17 XORs
above the original hot step.  The circuit count is exact, but the systematic
state permutations and accumulator have not yet been integrated into the
end-to-end benchmark.  The striped outer probe adds 1.23--1.27 ms; removing
the parity pass can only be credited after a matching end-to-end measurement.

## Differences from the existing full-split proof

The theorem-facing full-split construction currently uses a uniform
coordinate interleaver, an independent transitive `A_i` at every inner step,
and an independent random split of each 128-coordinate BCH word.  None of
those three randomizations is available in the frozen fast target.

Consequently:

* a fixed outer word is randomized only through the order of its 64-bit group
  values, not through a uniform Hamming slice;
* the local transition depends on the actual 64-bit vectors, not only their
  weights;
* the ordinary BCH spectrum and the old hypergeometric split law are not, by
  themselves, a valid local kernel for this construction;
* the existing identity matrix ledger can supply proof techniques, but does
  not certify the group-chain construction without a new group-profile bound.

## Proof gates

The proof branch proceeds only if all of the following gates close:

1. no short-chain or autonomous-orbit low-weight obstruction is found;
2. the deterministic fixed-split transition admits a rigorous local bound;
3. that bound composes along one chain without assuming weight-slice
   exchangeability;
4. the randomized tiled layout supplies a sufficient group-occupancy bound
   for every outer word, with co-tiled-block correlations retained, and the
   independent group permutation supplies the chain-placement law;
5. the bound combines with the exact punctured EBCH ambient spectrum and the
   codimension-24 graph factor over every nonzero outer weight;
6. exact or outward arithmetic gives the requested finite-distance margin.

If gate 1 or 2 fails moderately, the first fallback is a fixed systematic
basis change of the same BCH encoder or a very small fused mixer.  A repair
that restores coordinate-scale data movement or a Singer-sized hot-path cost
is considered a performance-significant failure and must be reported before
the construction changes.

## Contiguous-layout obstruction

The frozen contiguous layout fails the first adversarial gate.  Put
`T=floor(d/64)=2949`.  If every nonzero input group of an outer word occurs in
the last `T` chain positions, the first `B-T` emitted groups are zero and the
remaining output has weight at most `64*T=188736<=d`, independently of all
local BCH behavior.

For each unpunctured 64-coordinate data-message block, the restriction of the
codimension-24 graph map to that block has a kernel of dimension at least 40.
Choose a nonzero message in this kernel and zero elsewhere.  The data EBCH
word and the parity EBCH word are identical.  Since both fixed BCH halves are
invertible, this outer word has exactly four nonzero contiguous-layout groups:
the two data halves and the two parity halves.

The two parity groups are common to all such one-block words.  Conditional on
both parity groups landing in the last `T` positions, it is overwhelmingly
likely that at least one of the 16128 disjoint unpunctured data-group pairs is
also wholly in that suffix.  Therefore the ensemble failure probability has
a rigorous lower bound close to

```
C(T,2) / C(B,2) ~= 2^-6.95,
```

far above the desired `2^-40`.  This is a construction obstruction, not slack
in the upper-bound ledger.

## Low-data-movement repair

The tiled layout above is the current repair.  It replaces the fatal
two-contiguous-group support of a local EBCH word by one distinct group per
nonzero coordinate.  It also removes the parity block entirely.  The graph
block does not create dense special groups: its coordinates occupy the 128
lanes freed by the punctures, one coordinate per group.

The transpose of the data layout is a cache-tiled 128-by-64 gather followed
by 64 ordinary outer BCH transposes.  It adds a streaming/tiled 32 MiB layout
cost, not random element-scale reads across the full buffer.  The first direct
probe measures 7.23 ms for the existing sequential outer and 8.49 ms for the
striped outer model.  In comparable full paths, the medians are 12.69 and
13.92 ms on the same noisy run.  Thus the measured stripe cost is 1.23--1.27
ms.  Applying that delta to the clean 11.17 ms baseline projects 12.4--12.5
ms before removing the parity pass or crediting the systematic-circuit
saving.

For spreading, coordinates from the same outer block can never collide in a
physical group.  Two coordinates from distinct data blocks collide with
pair probability

```
(63/16383) * (1/128) = 63/2097024,
```

because their lanes must first land in the same tile and then in the same
tile column.  A graph/data coordinate pair collides with probability at most
`1/B`; two graph coordinates never collide.  Pair probabilities alone are
not a collision-tail proof: several coordinates from the same pair of
co-tiled blocks are correlated.  Step 5 must retain that tile-level
dependence, rather than applying an invalid independent-birthday bound.

## Systematic-output local refinement

A fixed basis change in the same BCH code makes the emitted half systematic:

```
E_sys(v) = (v, P v).
```

Independently permute the 64 incoming state coordinates before every XOR, and
independently permute the 64 lanes placed in every striped input group.  These
permutations are local to hot state or one cache tile; they do not reintroduce
global coordinate-random memory access.  Conditional on the state and input
weights, their supports are independent uniform subsets.  Conditional on XOR
weight `t`, the drive is uniform on the weight-`t` slice.

If the pre-accumulator drive has weight `s`, the exact accumulator law sends
it to weight `t`; the emitted weight is exactly `t`.  If the next state has
weight `q`, the systematic split enumerator `C[t,q]` obeys the rigorous caps

```
C[t,q] <= min(A_EBCH[t+q], C(64,q)),
sum_q C[t,q] = C(64,t).
```

The second cap is the right-half column marginal: `P` is bijective, so exactly
`C(64,q)` inputs produce a state of weight `q`.  Allowing every unknown row to
choose its worst capped distribution independently gives a robust Bellman
operator.  The exact systematic split
rows through `t=10` and their complements are committed in
`scripts/ebch128_systematic_split_slices.csv`.  At the original pole
`2333/2373`, exact rational arithmetic reports

```
log2(rho)       = -0.624039853324
log2(prefactor) = -0.024525793210
Chernoff crossing at d=floor(.09 N) = 7419 groups.
```

A finite exact-rational burn-in is essential for the displayed prefactor.  A
generic Perron-vector domination loses more than ten artificial bits, whereas
the reached-vector argument follows the exact robust operator for 12 steps
and then proves `O(v)<=rho*v` for that exact rational vector.  The complete
episode envelope holds from length one onward.

The exact state-column extrema through weight 10 and their complements give a
second local fact: excluding `(t,q)=(64,64)`, every exact-support cancellation
transition has probability at most

```
1/C(64,11) = 2^-39.435727630104.
```

For the identity-mixer pole above, the reached-vector prefactor divided by
`rho` is below two, so the old single-pole group ledger could charge a
conservative `2^-38` for each nontrivial episode restart.  In the multipole
PA1 ledger the exact `2^-39` turnoff charge and the pole-specific
`prefactor/rho` restart factor are retained separately.
The all-ones transition is handled separately: an input group of weight 64
requires at least 63 coordinate collisions and belongs in the outer spreading
tail.

An independent gap-composition audit found that the previous `2^-61.77`
identity-mixer estimate had implicitly anchored the first active group and
was invalid.  Averaging its actual random position makes the one-block bound
only about `2^-36.88`.  A pole scan confirms that the spectrum-only systematic
relaxation permits an artificial 25.5-bit-per-live-group trajectory; more
poles or one additional exact boundary row do not repair it.

PA1 removes that proof defect.  The exact accumulator weight chain lifts the
robust live rate to about 32 bits per group without needing the unknown middle
split enumerator.  A direct uniform-gap episode sum, retaining the
pole-specific restart factors, gives

```
one-data-block, graph-zero family       2^-44.922833091179
all collision-free outer weights 21..106 2^-44.922832756341
```

The second line differs by less than `4e-7` bits: the low ledger is almost
entirely the one-block family, whose coordinates cannot collide under the
striped layout.  These are floating diagnostics.  The remaining step-5 task
is to bound the rare co-tiled multi-block collision tail; step 6 then replaces
the floating PA1 pole rows and placement sums by rational or outward
certificates.

## Correlated tile-subspace obstruction to the occupancy ledger

The pair-collision tail is not the hard case.  Retaining the exact marked and
unmarked tile geometry gives, through group support 106,

```
two data blocks, distinct tiles       2^-66.9653
two data blocks, same tile            2^-53.1409
```

Both are safely below the one-block `2^-44.92285` total.  However, an entire
tile contains 64 data blocks, hence a 4096-dimensional message subspace, while
all of its outer coordinates occupy only 128 physical groups.  A
within-tile-hypergeometric/exponential-partition calculation exposes the
result immediately: the occupancy-only partial terms through support 106 are
about `2^-41.83` for three active blocks and `2^-10.94` for four.  Larger
tile-local subspaces get still worse.

This does not prove that the fast tiled construction has bad distance.  The
loose terms count many values sharing the same small group support as if their
late-placement events were independent.  Closing the fast construction now
requires a group-weight- or subspace-level inner argument, rather than a
first-moment sum depending only on the number of occupied groups.

Several construction-side spreading repairs were performance-probed before
changing the theorem target.  A global affine-line incidence over `F_128`
gives exactly 64 coordinates per group and makes every pair of data blocks
share at most one group.  Its straightforward outer transpose is much too
expensive: 29.52 ms as a pull and 27.12 ms as a push, versus 7.18 ms for the
fast tiled outer.  Four `F_64` supertiles with pair intersection at most two
still cost 18.68--18.83 ms and are also rejected.

A two-band stripe was generated and fused, but its proof obstruction is real
at the first-moment level.  The exact joint half-weight calculation for the
valid capacity-two star pattern (all blocks collide in one band and form
pairs/singletons in the other) reaches about `2^-29.44` at seven active blocks,
`2^-11.35` at eight, and becomes positive at nine.  Capacity one improves this
by only four to six bits.  Thus repairing only the product-marginal relaxation
does not certify two bands.

The current low-overhead repair is instead a fixed three-band split of sizes
`42/43/43`.  Every band has 256 tiles of 64 blocks.  The scratch implementation
maps logical block `(t,l)` to tile labels

```
t,  t + 9 l,  t + 20 l       (mod 256),   0 <= l < 64.
```

Every pair of band labels identifies at most one block.  The map is also
Pasch/intercalate-free: a four-block configuration with all six pairs
co-tiled would force a lane difference of 128, outside the lane set.  Exhaustive
checking in `scripts/certify_three_band_tile_map.py` reports zero such
configurations.  The same certificate checks exact balance (64 blocks per
tile in every band), pair capacity one, and counts the next-densest
four-block family.  There are exactly 281856 K4-minus-one-edge sets, split as
79616, 122880, and 79360 according to which band contains the single
collision.

Separate exact transpose circuits for the three coordinate bands use
`212 + 463 + 238 = 913` circuit XORs; the two later bands also require 128
output merges.  The fused implementation was checked bit-for-bit against the
full transpose.  Three sequential paired runs at `K=2^20` measured outer-kernel
deltas of `+0.277`, `+0.481`, and `+0.319` ms over the one-band stripe.  This is
the current performance/proof target.  It is much cheaper than the earlier
materialized two-band fallback, whose outer cost was 3.75 ms above the stripe.

The fixed bands deliberately do **not** use the unjustified hypergeometric law
of a fresh random `42/43/43` coordinate split.  Exact binary checks give full
inside ranks `42/43/43`, full outside rank 64, and outside distances `5/4/5`.
The complete outside spectra are obtained by enumerating duals of dimensions
22, 21, and 21 and applying MacWilliams.  Exact meet-in-the-middle rows through
outside weight eight and exact inside weights zero/one (plus complements) are
committed in:

```
scripts/ebch128_fixed_band_projection_spectra.csv
scripts/ebch128_fixed_band_outside_slices.csv
scripts/ebch128_fixed_band_inside_boundary.csv
```

Using those rows, ordinary-spectrum transportation caps, and an auxiliary
Chernoff tilt for group support at most 106 closes every one-tile star with at
least two blocks.  An earlier diagnostic incorrectly used 64 columns for all
three bands.  The dimension-correct exact rational certificate, including the
puncture and graph corrections, gives

```
sum of all stars with at least two blocks       <= 2^-56.13995.
```

The common-pole (`q=.05`) tail beginning at size six is below `2^-53.8163`.
Both statements are checked by `scripts/certify_three_band_star.py`; the exact
one-block PA1 row remains the global low-weight leader.
Specialized Finner bounds also put a two-band three-block path at about
`2^-66.27` and a three-color triangle at about `2^-52.42`.

The missing full joint split enumerator can be avoided.  For a tile cluster
of size `d`, the coordinate OR kernel has the nonnegative rank-one
decomposition

```
q^OR(x_1,...,x_d) = q + (1-q) product_i 1[x_i=0].
```

Tensoring this identity over the coordinates and averaging fixed-weight
slices leaves a completely positive symmetric tensor.  Generalized Holder
therefore gives

```
K_d(a_1,...,a_d) <= product_i K_d(a_i,...,a_i)^(1/d).
```

Every collision pattern consequently separates into one-codeword factors.
For each unknown fixed-band split cell `N[a,b,c]`, every committed exact
diagonal, marginal, complement-projection, low-outside, or inside-boundary
row containing the cell is an individual upper cap.  On each exact
total-weight diagonal, greedily filling the largest factor cells up to the
smallest such cap is a rigorous relaxation of the missing joint enumerator.
This is implemented in `scripts/probe_bch_three_band_cell_cap_motifs.py`.

The final small-support certificate enumerates every compatible triple of
set partitions for sizes two through six and evaluates the cell-cap relaxation
with exact rational/outward arithmetic.  Its post-graph bounds are

```
s=2                                         2^-56.2689
s=3                                         2^-62.1571
s=4                                         2^-63.8136
s=5                                         2^-60.3606
s=6                                         2^-53.7541
sum s=2..6                                  2^-53.5042
```

These are deliberately broad all-pattern sums: they include collision-free,
star, and disconnected patterns, so reusing them as component envelopes is a
safe overcount.  `scripts/enumerate_three_band_small_profiles.cpp` performs
the packed profile enumeration and `scripts/certify_three_band_small_motifs.py`
certifies the arithmetic.

The one-band star operator remains separate and is certified as described
above.  For the connected tail, pair capacity implies the
per-block degree constraint `d0+d1+d2 <= s+2`; degree-`d` vertex counts are
multiples of `d`; and a `d`-cluster and `e`-cluster in distinct bands have at
most one common block.  The integer relaxation in
`scripts/probe_three_band_tail_degree_mip.py` uses only these necessary
conditions.  It selects an `(s-1)`-star with one external block and two
cross-band spokes.  Applying the low-support Chernoff tilt and then
multiplying by `Bell(s)^3`--which counts every triple of set partitions,
including incompatible ones--gives the conditional size-7..64 frontier

```
sum over connected sizes 7..64             2^-55.2991
```

with the frontier terms `2^-55.3138` at size seven and `2^-61.9411` at size
eight.  The computation is in
`scripts/probe_three_band_decorated_star_tail.py`.  Despite the old name,
maximum tile degree 64 does not bound a connected component: a component can
walk through several tiles.  Under group-support cutoff 106, punctured BCH
distance 21
and physical group capacity 64 give the deterministic restriction

```
21 s <= 64 * 106, hence s <= 323.
```

Thus sizes 65..323 still need a coarse connected-component envelope.
Disconnected patterns can be composed by the exponential formula only after
that envelope is available; the singleton component activity is about
`2^-20.91` before the one global 24-bit graph constraint.

The seven- and eight-block frontier MIPs have now been exported with an exact
arithmetic objective.  The diagonal cluster moments are exact rationals; each
`d`-th root is rounded upward to 160 dyadic bits by an integer comparison; and
the resulting factor logarithms are rounded upward on a `2^20` integer grid.
Exact SCIP 10.0.2, with presolving, separation, and conflict analysis disabled,
proves the same optimum for each of the three choices of required collision
bands.  Completing the transcript with VIPRCOMP and checking it with VIPR
verifies respectively

```
s=7 objective <= -244394546 / 2^20
s=8 objective <= -326550734 / 2^20
```

The conflict-analysis restriction is certificate hygiene: SCIP's exact primal
and dual bounds were unchanged, while its default conflict transcript contained
two floating aggregation rows that VIPRCOMP could not dominate by less than
`8e-15`.  With conflicts disabled, VIPRCOMP performs zero repairs and VIPR
checks the original unpresolved integer models.

This closes the arithmetic and finite combinatorial claims for sizes seven and
eight without changing the encoder.  It does **not** yet prove the asserted
decorated-star extremality for every size through 64.  The abstract
pair-capacity relaxation admits balanced linear-hypergraph profiles at larger
sizes, so the remaining tail proof should use an exchange inequality or a
coarser large-component envelope, not merely extrapolate the small MIP scan.
The actual affine tile map does contain a 64-block `4 x 4 x 4` regular
substructure (take 16 equally spaced first coordinates and four equally spaced
lanes); at tilt `0.05` its individual cell-cap exponent is below `-1892`, so it
is not a numerical obstruction, but it is a useful required test case for the
large-component lemma.

There is a more useful positive structural lemma.  Regard data blocks as the
edges of the three-partite linear hypergraph whose vertices are the tile
labels.  For an `s`-edge subhypergraph let `n_b` be the number of used tiles in
band `b`.  Pair capacity one gives `s <= n_b n_c` for every pair of bands.
Hence some band has at least `sqrt(s)` used tiles and therefore contains a tile
of degree at most `floor(sqrt(s))`.  Repeating after deletion gives a tile
peeling order of width at most 17 throughout the only relevant range
`s <= 323`.  This statement includes degree-one leaves: by itself it does
**not** promise a collision tile of degree at most 17.  A one-tile star is the
basic counterexample.  The final composition argument therefore permits a
peeled collision tile of any physical degree 3 through 64 rather than relying
on the width statement.

For `s <= 64`, pair capacity gives width at most eight.  Equality at 64 would
force exactly eight tiles of degree eight in every band and all 64 pairs
between the first two tile sets.  Normalize labels by the inverse `57` of
`9 mod 256`; an edge becomes

```
(a, b, c) = (a, a+l, 116 b - 115 a),  0 <= l < 64.
```

The complete first/second-band pair set puts the eight first-band labels in
one cyclic arc of length 64.  Let `U=-115 A` and `V=116 B`.  Eight third-band
labels would imply `|U+V|=|U|=8`.  Multiplication by 116 has fibers of size
four, so `|V|>=2`; two translates of `U` must coincide.  Thus `U`, and hence
`A`, has a nonzero period.  Every nontrivial subgroup of `Z/256` contains 128,
forcing an antipodal pair in `A`, impossible inside that arc.  Therefore the
equality case cannot occur and every set of at most 64 blocks has peeling
width at most seven.  The finite modular checks and the 4-regular regression
witness are in `scripts/certify_three_band_peeling.py`.

This suggests replacing global decorated-star extremality by a local transfer
certificate.  If a peeled tile contains `r` active blocks, pair
capacity makes the `r` attachment tiles distinct in each of the other two
bands.  Conditional on an attachment tile already having union weight `u`, a
new independently permuted fixed-weight support of weight `w` contributes the
exact incremental moment

```
H_n(u,w) = sum_i C(u,i) C(n-u,w-i) / C(n,w) * q^(w-i).
```

This is monotone increasing in `u`: couple the old union of size `u` to one of
size `u+1`; the number of new coordinates can only decrease, and `q^x`
increases as `x` decreases for `0<q<1`.  Generalized Holder splits the one
coupled `r`-way kernel in the peeled band into its rooted diagonal factors.
The fixed-band BCH cell caps then give a one-codeword transfer factor depending
only on the two attachment weights.  Across all `r` attachments the sum of the
`2r` old union weights is at most the single global support cutoff 106.

`scripts/probe_three_band_peeling_transfer.py` rounds every attachment weight
upward to a four-coordinate grid.  It charges each rounded value by the least
integer weight that rounds to it, so every real boundary vector remains in the
DP while the shared budget stays exactly 106.  Its current floating diagnostic
gives

```
r=3       transfer <= 2^-16.4706
r=4       transfer <= 2^-42.3673
r=5       transfer <= 2^-59.4732
r=6       transfer <= 2^-74.2380
r=7       transfer <= 2^-81.2787
r=8..17   transfer <= 2^-90.2505 or smaller
```

An exact-weight boundary grid improves the delicate degree-three row to
`2^-30.5461`.  Direct degree-four through degree-64 probes all contract as
well; selected large-degree rows are `2^-141.86`, `2^-211.41`, and
`2^-428.05` at degrees 18, 32, and 64.  Thus no width-17 assumption is needed.

These are local transfer factors at the Chernoff tilt `q=0.05`.  The conversion
from that tilt to the target support event is charged **once globally**, not
once per peeled cluster.  The figures are not yet theorem-facing: their BCH
cell-cap sums, rooted moments, logarithms, and DP maxima still need outward or
exact certification.  They nevertheless show ample slack for every peeled
degree from three through seventeen without changing the construction.

The residual case needs a separate statement.  Once all tile clusters of
degree at least three have been eliminated, every collision tile has degree
two.  Contracting each such tile to an edge produces a properly three-edge-
colored simple graph on the active blocks, with maximum degree three.  It is
therefore not merely a union of paths and cycles: cubic branching at a block
is possible.

Fortunately this entire core has a coarse exact bound.  Drop both connectivity
and the prohibition against pairing the same two blocks in two colors.  The
three color classes are then independent matchings.  At `q=1/20`, exact
rational rooted moments and the cell-cap greedy relaxation verify, for every
subset `S` of paired bands,

```
F_S <= 2^(|S|-57).
```

After extracting `2^-57` per block, every matching edge has activity four.
Thus the weighted matching count obeys the integer recurrence

```
J_0=J_1=1,      J_s=J_(s-1)+4(s-1)J_(s-2).
```

Using the connected placement bound `2^14 64^(s-1)`, division by `s!`, and
the one codimension-24 graph constraint gives

```
sum_(s=7)^323 2^14 64^(s-1) J_s^3 / (s! 2^(24+57s))
    < 2^-346.1471                         at q=1/20.
```

The single support-cutoff conversion `(181/50)^106` then gives a target-pole
bound below `2^-149`.  `scripts/certify_three_band_pair_core.py` performs all
factor comparisons, matching recurrences, summation, and the final comparison
with exact integers/rationals; only its printed logarithms are approximate.
This closes the isolated pair-collision core with enormous slack.

### Composition of high clusters

There is now a finite composition route that avoids both the false width
inference and decorated-star extremality at arbitrary size.  Repeatedly remove
all active blocks in any tile of current degree at least three.  The removed
batches are disjoint, each physical tile is a root at most once, and the
residual incidence hypergraph has only singleton and pair tiles.  Reverse the
removal order and choose a canonical spanning attachment for every nonroot
batch.  Given a parent block, a child `r`-batch has at most

```
3 * 63 * 2 * C(63,r-1) = 378 C(63,r-1)
```

encodings: choose the shared band, the neighboring anchor block, one of the
other two root bands, and the remaining root-tile lanes.  Including at most
323 choices of parent block, the exact-weight degree-three transfer is the
worst child row and still contributes below `2^-2.77`; every larger degree has
more slack.  This is a deliberately coarse connected-growth count.

The saturated-boundary transfer must not be used for the first batch, because
that would discard the rarity of the boundary itself.  The static Holder
factors give a cleaner root treatment.  If no tile contains blocks from two
different high batches, every high-batch block has its root degree `r` and
degree at most two in both external bands.  Directly assigning the pessimistic
type `(r,2,2)` to every block gives, after all parent and lane choices, a child
activity below `2^-113`; degree three is the worst row over `r=3..64`.
Residual pair-core blocks have type at most `(2,2,2)` and cost below `2^-54`
before their placement choice.

These static factors close every such connected component except the already
enumerated sizes at most six.  More explicitly, with one high batch and zero
residual blocks the component is an exact one-tile star.  With one high batch
and one or two residual blocks, only degrees three/four can evade the coarse
bound, giving at most six blocks.  With two high batches, the second batch
already supplies more than 113 additional bits.  Extra batches and residual
blocks only decrease the bound.

If two high batches do meet through a later root tile, the final shared tile
has degree `r_2+1`, rather than degree two.  The correct pessimistic joint
motif assigns `(r_1,r_2+1,2)` to the earlier anchor, `(r_1,2,2)` to its other
blocks, and `(2,r_2+1,2)` to the later batch.  Over degrees at least three the
worst case is `r_1=r_2=3`; its BCH factor is about `2^-256.999`.  After root,
lane, attachment, cutoff, and graph charges the six-block base is
`2^-39.838`.  That size is already in the exact small ledger.  Any component
of size at least seven has at least one further batch or residual block.  A
high batch contributes at most `2^-2.77` under the exact-weight conditional
transfer, while a residual block fits one of the already pessimistically
allocated external pair slots.  Summing an arbitrary nonempty continuation
geometrically leaves the complete high-high continuation below about
`2^-42.38`.

The safe size-7..13 MIPs below are consequently redundant cross-checks rather
than a premise of the large-component composition.  They are retained because
they exercise every small mixed degree profile at the same pole.

The safe degree-profile MIP at the common tilt `q=0.05`, including the full
`Bell(s)^3` overcount and cutoff conversion, gives

```
s=7   -52.7532       s=8   -61.9411
s=9   -67.9225       s=10  -72.1933
s=11  -75.5456       s=12  -78.3644
s=13  -80.8434
```

and selects the decorated-star degree profile in every row.  Therefore the
composition reduces the formerly missing sizes 65..323 to a finite list of
local inequalities; it no longer assumes decorated-star extremality at large
size.

The main local rows are now hardened.  Exact rational moments and outward
160-bit roots prove every fully decorated child activity below `2^-112`, every
fully decorated root below `2^-125`, every degree-at-least-five root below
`2^-165`, the corrected adjacent degree-three base below
`(3/4) 2^-39`, and the one-residual family below `2^-50`.
The exact-weight degree-three boundary DP uses exact moments plus explicitly
guarded outward binary64 accumulation and proves
`tau_3 <= (11/16) 2^-30`.  These checks are in
`scripts/certify_three_band_composition_local.py` and
`scripts/certify_three_band_degree3_transfer_outward.py`.

Substituting those rational envelopes into the positive case ledger gives

```
adjacent-high nonempty continuation       <= 2^-41.0917
sum of all high-cluster cases             <= 2^-41.0886.
```

The arithmetic is exact in
`scripts/certify_three_band_composition_ledger.py`.  The remaining conditional
rows are now outward-verified in
`scripts/certify_three_band_degree4_64_transfer_outward.py`: degree four is the
worst of that range at `2^-10.261` after attachment choices, versus the
certified degree-three envelope `2^-2.780`, and every larger degree is safer.
Thus the high-cluster composition has no remaining local premise.  The
size-7..13 MIPs can still be transcript-certified as an independent regression
layer; they are not needed by the composition after these local inequalities
close.

The intended proof order is consequently:

1. freeze the exact three-band construction and its implementation cost;
2. certify the EBCH fixed-band spectra, split rows, PA1 kernel, and turnoff;
3. certify pair capacity and the affine tile-map peeling facts;
4. certify local `r>=3` transfer factors with a shared boundary budget;
5. certify the subcubic degree-two core and all leaf attachments;
6. for total group support at most 106, sum connected components through size
   323, compose disconnected components, and assemble the outward failure
   ledger; then prove a separate global outer/group-profile bound for support
   above 106;
7. analyze Feistel separately only after the algebraic proof is closed.

The first half of step 6 is now exact.  The target-pole singleton activity is
`2^-20.89369` before the graph factor; the broad size-two-through-six component
envelope is `2^-28.81636` before the graph factor.  The common-pole tail
combines the high-cluster, star, and pair-core envelopes `2^-41`, `2^-53`, and
`2^-149`.  Charging the cutoff conversion and the corrected graph factor only
once globally, and composing arbitrary disconnected components by a positive
geometric-series upper bound, gives

```
support-at-most-106 total failure                <= 2^-40.90535
margin below 2^-40                                  0.90535 bits.
```

The exact arithmetic is in `scripts/certify_three_band_final_ledger.py`.  This
is a low-support certificate, not yet the final theorem ledger.

### Remaining high-support gate

Words with group support above 106 are not covered by the component ledger.
The exact scalar PA1 envelope proves `P_g <= .181^g` through the low-support
range (with its tightest ratio at `g=21`, about `2^-0.01196` below the target)
and continues to hold with slack for a substantial intermediate range, but it
eventually loses to `.181^g` for dense group inputs.  More local motif
enumeration cannot repair that: a large one-tile star can have a positive
target-pole moment, so the support cutoff is essential.

At that stage, the remaining proof needed a global consistency bound coupling the
randomized tiled outer layout to outer-code membership.  The promising route
is a group-profile analogue of the existing identity-matrix ledger: condition
on per-group or per-block weights, use the random lane and coordinate
permutations to charge EBCH membership densities, and combine that outer
profile with a matrix/transfer PA1 bound. The packet-profile certificate below
supplies that bound. No runtime construction change was needed.

## Complete packet-profile certificate for `g=4`

The packet-profile route closes the high-support gate without changing the
construction. Fix the frozen `g=4` ensemble described above. For a sampled
setup, define

```text
Z_d = #{nonzero messages m : wt(C(m)) <= d},
d = 188743.
```

The remaining probability is over the independent lane bijections, packet
permutation, and recursive state permutations sampled during setup. The
packet-profile orbit lemma gives the exact normalization

```text
Q_4(a) = M! / product_j a_j! * 4^(a1+a3) * 6^a2,
M = 524288.
```

The certificate partitions every feasible profile by its exact positive
support. There are 30 feasible supports. Each support has an exact integer
hull and an exact rational simplicial mesh. The verifier checks every hull,
facet incidence relation, orientation, boundary, and determinant volume.

Each unchanged root cell uses one fixed rational mixture of affine witness
bounds. Convexity extends the vertex inequalities across that cell. A root
cell whose fixed mixture is too weak receives an exact rational binary space
partition. Each BSP leaf uses one fixed witness at all of its vertices.
The verifier reconstructs every leaf from its half-space constraints and
checks that the BSP partitions the complete parent cell.

Sparse witnesses set fugacities outside their exact support to zero. This is
valid because the local generating functions have nonnegative coefficients.
Deleting monomials for absent packet classes cannot increase the moment.
The verifier rejects a sparse witness whenever its zero-fugacity class is
active on the cell under evaluation.

The outer witnesses combine the asymmetric one-conditioned-row bound with
the exact graph and puncture average. The inner witnesses use the shared-drive
65-state transfer bound. Discovery uses binary64 optimization, but discovery
values are not theorem inputs. The independent verifier treats the selected
parameters as exact dyadics, encloses local nonnegative arithmetic upward,
and evaluates the final Collatz inequalities with directed arithmetic.

For each root cell `c`, let `n_c` be the verifier's integer-profile count
upper bound. Let `U_c` be the largest outward endpoint among its fixed-mixture
vertices or all vertices of its BSP leaves. Closed boundaries may occur in
more than one cell, so the proof uses the safe overcount

```text
E[Z_d] <= sum_c n_c 2^U_c.
```

The independent replay audits 40902 root cells. It replaces 899 roots by
2132 exact BSP leaves and evaluates 2144 fixed witnesses. Directed
log-sum-exp gives

```text
log2 E[Z_d]
  <= -61.788151555343981861003442630812211839178495852356.
```

Markov's inequality now yields

```text
Pr[d_min(C) <= 188743]
  <= Pr[Z_d >= 1]
  <= E[Z_d]
  <= 2^-61.78815155534398186100.
```

Thus the frozen ensemble has `61.7881515553` bits of first-moment security.
The margin beyond the required 40-bit threshold is `21.7881515553` bits.
Consequently, at least one setup gives a binary
`[2^21,2^20,d_min>=188744]` code.

The canonical report is
`out/g4_end_to_end_bsp_outward_certificate.json`, with SHA-256
`9a416987feda6c0edf15ea2af91b30eb84d3be8b97a002831cf85f03e7086364`.
The report authenticates every source ledger, witness report, and evaluated
inequality. `PROOF_STATUS.md` gives the complete replay command and the hashes
of its four structural inputs.
