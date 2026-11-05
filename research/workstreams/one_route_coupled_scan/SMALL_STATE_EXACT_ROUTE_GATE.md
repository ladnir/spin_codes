# Exact paired-route gates for the small-state GF32 proposal

The whole-state counting bridge now gives directed finite bounds for
messages on seven through 128 actual outer groups. Together with the
retained direct gates, this covers every message using at most 128
groups at 35% symbol distance. Separate gates below cover all messages
on 2048 through 2073 groups and on 4096 through 4108 groups. The intervening and remaining
source classes are open. This is a partial proof, not a whole-code
certificate or an encoding-time result.

The construction is unchanged: K=2^20, N=2^21, four GF256-valued RS[64,32]
rows per outer group, invertible column mixers, the alternating paired
route, and a 16-symbol state emitting 64 symbols per batch. Both DHD
diagonals have independent uniform nonzero GF(2^32) entries. Setup is
sampled once and fixed for every message; the terminal state is discarded.
See [the construction and sparse gates](SMALL_STATE_PROOF_GATE.md).

## What changed in the proof

The earlier floating screen replaced a regional coefficient by a tilted
IID sum divided by the probability of the required histogram. That is a
valid majorant, but it retains contributions from the wrong pair counts.
The exact coefficient can be much smaller; removing just the histogram
divisor does not measure this improvement.

Keep the [whole-state bridge](SMALL_STATE_RANK_BRIDGE.md), including its
fourth root before the positive output-zero expansion and its separation
of actual GF32 sources from auxiliary GF256 sources. Only the route
average changes. No additional randomness or encoder operation is used.

For g actual groups and f full pairs in one matching, a region contains
s=g-2f singleton pairs and f full pairs among 4096 pair slots. A physical
pair straddles two different RS[16,4] blocks. For each batch, sum the
two-state local matrix over its exact striped occupancies, retaining
singleton and full-pair polynomial variables. Raise this matrix polynomial
to the 128th power in chronological order and extract coefficient (s,f).
Its exact normalizer is

    binom(4096,s) binom(4096-s,f) 2^s.

Write R_f(z,alpha) for the normalized regional matrix. The local matrix
is raised entrywise to alpha before this route average. A fixed source
set has fixed matching counts (f_even,f_odd); they are not resampled from
column to column. There are 256 alternating regions, so the positive
free-terminal envelope is

    [1,0] (R_f_even R_f_odd)^128 [1,1]^T.

The all-ones terminal uses the proved inequality f_I<=M_I*ones. It does
not require a terminal reset or charge an unused final refresh.

All the present gates use the missing-group cutoff h=0. With
beta=256^160/255^32 and c_g=1-g/256, the original source-count bound is

    binom(8192,g)
      * [2 beta^g z^(-734003)/(255 c_g)]^alpha
      * max_(f_even,f_odd) [1,0](R_f_even R_f_odd)^128[1,1]^T.

The factor 2 bounds the full-state refresh density over the entire word.
The binomial counts all source sets, so profiles are maximized, not added
as another union. Each count satisfies 0<=f<=floor(g/2). The union of the
two matchings is a cycle on 8192 groups. A proper nonempty source subset
has at most g-1 induced edges, hence f_even+f_odd<g. This removes only the
impossible simultaneous full-pair endpoint for even g in this range.

For g=65,...,128, partition the actual source sets by their number c of
occupied cyclic runs. Their exact multiplicity is

    N_g(c)=(8192/c) binom(g-1,c-1) binom(8192-g-1,c-1).

Each such set has f_even+f_odd=g-c. Replace the full binomial times the
worst profile by the sum over c of N_g(c) times the worst ordered profile
on that diagonal. The counts partition binom(8192,g) exactly and stay
outside alpha. This retains the actual-empty pruning and the free
terminal. [Bridge Section 11](SMALL_STATE_RANK_BRIDGE.md) proves the census.
No encoder operation or setup distribution changes.

## Fixed witnesses and directed checks

The newer [shared floating engine](theory/probe_small_state_shared_route.py)
computes all regional coefficients up to a chosen source count at once.
For n pair slots, put

    C_n(a,b)=binom(n,a) binom(n-a,b) 2^a.

Store each matrix coefficient divided by C_n. A squaring step combines
left degree (a,b) and right degree (c,d) with positive weight
C_n(a,b) C_n(c,d)/C_(2n)(a+c,b+d). The weights sum to one for each target
degree by the multinomial Vandermonde identity. The initial n=32 table
checks its occupancy multiplicities against C_32 exactly; seven squarings
reach n=4096. Truncating to a+2b<=G preserves every requested coefficient
because all polynomial degrees are nonnegative. Initial terms also obey
a+b<=32. None of these changes replaces the exact route by an IID model.

Nine comparisons against the original engine, at g7,17,32 and three pair
counts each, differed by at most 7.2e-15 in regional log2 entries. A G64
table and full profile screen took about 1.0--1.3 seconds and 51 MB peak RSS.
Positive underflow is rejected; floating normalizer residuals are only
diagnostics. The separate
[shared directed verifier](theory/verify_small_state_shared_route.py)
instead retains unnormalized Arb coefficients and divides by the exact
C_4096 once. Arbitrary exponents avoid the floating underflow problem.
Its g7,23,32 regression bounds agree with the original verifier; its new
g33,...,64 bounds have passed the directed replay below.

Let z7=9634151/9689600 and z16=2238127/2263291. For g=7,...,13 use

    z_g = z7 + (z16-z7)(g-7)/9,       alpha=1.

The remaining witnesses and proved thresholds are:

| Actual source groups | Output marker | Alpha | Verified failure bound |
|---|---|---:|---:|
| Each g=7,...,13 | Rational interpolation above | 1 | <2^-900 per g |
| 14 | 532495673/536870912 | 1 | <2^-900 |
| 15 | 136298763405/137438953472 | 1 | <2^-100 |
| 16 | 308022019811679423/311064346442596352 | 4/5 | <2^-100 |
| 17 | 99/100 | 1/2 | <2^-300 |
| 18 | 99/100 | 4/9 | <2^-300 |
| 19 | 12647/12800 | 7/20 | <2^-150 |
| Each g=20,...,22 | 101023/102400 | 3/10 | <2^-150 per g |
| Each g=23,...,32 | 49/50 | 1/5 | <2^-150 per g |
| Each g=33,34 | 49/50 | 1/5 | <2^-60 per g |
| Each g=35,...,48 | 97/100 | 3/20 | <2^-60 per g |
| Each g=49,...,53 | 193/200 | 1/8 | <2^-60 per g |
| Each g=54,...,60 | 24/25 | 1/8 | <2^-60 per g |
| Each g=61,...,64 | 191/200 | 1/8 | <2^-60 per g |
| Each g=65,...,77, run census | 191/200 | 1/8 | <2^-150 per g |
| Each g=78,...,100, run census | 47/50 | 1/10 | <2^-150 per g |
| Each g=101,...,128, run census | 23/25 | 1/10 | <2^-150 per g |

The [fixed-witness verifier](theory/verify_small_state_exact_route.py)
computes local majorants as exact rational numbers, then encloses every
matrix coefficient, fractional power, normalization, source factor, and
final bound in Arb. For alpha<1, the binomial source-set count remains
outside the fractional power, as in the formula above.
Acceptance compares each full profile bound directly against the exact
dyadic threshold. Displayed logarithms are not acceptance tests.

The root ran every listed witness serially at both 128- and 256-bit
precision. All profiles passed: respectively 16,24,25,35,36,48,49,63,64,80
profiles for g=7,...,16, or 440 profiles at each precision. The additional
g=17,...,32 checks cover 2784 profiles per precision. The shared g33,...,64
replays cover another 20672, for 23896 profiles per precision in total.
The g7 worst logarithm is enclosed near
-979.07857607390354 at both precisions. The two standalone g7 checks took
0.49 and 0.39 seconds with less than 60 MB peak RSS. The later checks
reused imported local-map tables and took 0.02--0.08 seconds each, below
58 MB RSS. The two g16 checks took 0.50 and 0.12 seconds, below 59 MB RSS;
the (8,7) profile's logarithm is enclosed near -135.28109951975655.
These are verifier times, not encoding times.

The additional fixed witnesses give worst-profile logarithms near -897.76
at g17, -333.87 at g18, and -576.35 at g32. The least negative score among
g19,...,32 is near -192.28 at g23. All acceptance tests still use the
exact thresholds in the table, not these displayed approximations.
The larger replays took at most 5.03 seconds per source count and remained
below 59 MB peak RSS. An initial 4-second cap interrupted g31; that run
was not counted. The completed g31/g32 replays used 10-second caps.

The five shared witness windows for g33,...,64 contain respectively
612,6524,3433,6007,4096 profiles. Every profile passes at both precisions.
Each window took at most 1.91 seconds and remained below 60 MB peak RSS.
The least negative new bound is near -67.78356667491629 at g64; acceptance
uses the exact 2^-60 threshold. These 32 source-count bounds sum to less
than 2^-55. The initial shared-verifier smoke test stopped before any
bound because dependency hashes were requested before lazy imports;
hash collection was moved after provider evaluation, and all reported
checks use the repaired source.

The run-census replays cover another 159104 ordered profiles per precision:
17105, 47359, and 94640 in the three new windows. Total coverage g7,...,128
is therefore 183000 profiles per precision. For each run class the verifier
converts every Arb upper endpoint to an exact rational and takes their
exact maximum. It then sums the exact integer counts times these maxima.
Acceptance compares this full rational sum against 2^-150; it does not
require every profile, charged the full binomial, to pass separately.
The count partition is checked exactly for every g.

Every new source size passes at both 128 and 256 bits. Selected logarithms
of the certified upper bounds are about -262.80359 at g77, -201.48917 at
g100, -203.77792 at g101, and -299.71506 at g128. The 64 new bounds sum
to less than 2^-144. Each window took at most 16.39 seconds and less than
60 MB peak RSS. These remain verifier costs, not encoding measurements.
An old-mode g7 regression retains the prior -979.078576 logarithm;
its optional census bound improves to about -987.813450. Exact checks
on cycles of lengths 4,6,8,10,12 verified 5446 source sets and 90 run
classes. The proof and implementation received separate source reviews.

The verifier records hashes of the actual local-provider modules. The
striped-geometry dependency is pinned to SHA-256
`04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053`.
Independent reviews checked the bridge, prefactor, route normalization,
free terminal, and directed verifier. Root reviewed the generalized
source-count interface and replayed its witnesses.

The existing per-class bounds for g=1,...,6 are, respectively,
2^-411, 2^-450, 2^-250, 2^-100, 2^-50, and 2^-45. Adding the new gates
still gives a combined bound below 2^-44 for all covered source counts.
This does not bound messages using more than 128 groups.

## Other checks retained from this round

The [direct rank-tagged gate](SMALL_STATE_SEVEN_GROUP_GATE.md) is a
separate valid interface. Its seven-group floating bound was +1473.40
in base-two logarithms at the inherited markers; bounded coordinate
tuning improved it to +815.65. Raising tau to 2^32-1 without retuning
made it worse. These failed upper bounds are not low-distance words.
The exact-route whole-state gate makes further seven-group tuning of
that larger operator unnecessary.

Exact small-field tests also passed 1875 source-rank fixtures, 339
full-rank cutoff cases, 625 top-half threshold identities, 36 diagonal
rank caps over 3510 assignments, and 12 projected rank caps over 27000
assignments. Eighteen MDS minors were checked in the projected cases.
These test the local algebra, not the full finite code.

The rank notes retain a refinement that sums exact four-block state-support
multiplicities. A separate Boolean full-column source cutoff is now
implemented as a floating control in
[probe_small_state_column_cutoff.py](theory/probe_small_state_column_cutoff.py).
It retains quiet input/output ranks before density domination and rejects
histories with 32 fully constrained columns. It does not charge those
ranks again. Selected concentrated profiles at g17,24,32 showed no useful
cutoff gain; fractional clean/dirty splitting instead cost about 0.20 and
0.07 bits in the latter two tests. This control is not used for any
directed claim above. At alpha=1 its regional clean-plus-dirty identity
agreed with the old operator within 3.3e-16 relative error.

## Reproduce serially

From this workstream, with NumPy and python-flint installed:

```text
python -B theory/verify_small_state_exact_route.py --groups 7 --precision 128
python -B theory/verify_small_state_exact_route.py --groups 7 --precision 256
```

Repeat for g=8,...,13 using the corresponding `--groups` value. For g=14,
add `--z 532495673/536870912`. For g=15, add
`--z 136298763405/137438953472 --target-bits 100`.
For g=16, add `--z 308022019811679423/311064346442596352 --alpha 4/5 --target-bits 100`.
For g=17,...,32 use the corresponding rational `--z`, `--alpha`, and
`--target-bits` from the table. Above g16 the verifier requires an explicit
marker; it does not extrapolate the earlier interpolation rule.
For g33,...,64 use the shared verifier, once per witness window. For
example, the final window is:

```text
python -B theory/verify_small_state_shared_route.py --groups 61 62 63 64 --z 191/200 --alpha 1/8 --target-bits 60 --precision 128
python -B theory/verify_small_state_shared_route.py --groups 61 62 63 64 --z 191/200 --alpha 1/8 --target-bits 60 --precision 256
```

Run both precisions separately for each witness. Use an external
55-second timeout and 1-GiB address-space cap. Each invocation also
enforces a 50-second limit and a 20-GiB free-space floor, and writes no
data files or bytecode caches.

For the three g65,...,128 windows, add `--run-census --target-bits 150`
and pass every integer in the chosen window to `--groups`. For example:

```text
python -B theory/verify_small_state_shared_route.py --groups 65 66 67 68 69 70 71 72 73 74 75 76 77 --z 191/200 --alpha 1/8 --run-census --target-bits 150 --precision 128
```

Repeat at precision 256, then repeat for windows 78..100 and 101..128
with their table markers. Run invocations serially. The verifier source
used for these replays has SHA-256
`25e77dce6b7b93fba74de62f70144a87ddefc4d1611b949c2ee760d160f322be`.

The bounded floating control is
`theory/probe_small_state_exact_route.py`; it is not a certificate.
It can additionally account for source sets by their number of occupied
runs on the cycle. For c runs their exact count is
`(8192/c) binom(g-1,c-1) binom(8192-g-1,c-1)`.
The original helper's displayed run-count bounds are floating diagnostics.
The newer shared verifier implements the directed census used above.

## A directed dense interval without changing the encoder

At the same 35% distance, all messages using between 4096 and 4108 actual
outer groups have combined setup-failure probability below 2^-100.
This is a new, disjoint covered interval, not coverage of the gap after 128.
Use the joint source correction and the raw empty-input matrices from
[bridge Sections 8, 9, and 12](SMALL_STATE_RANK_BRIDGE.md).

The fixed witnesses are z=370391/5732639, alpha=3/20, proposal tilt 3/2,
and pseudocount 1/65536. At the lower source size 4096, partition every
full-pair count 0..2048 into 121 intervals, using a maximum of 128 cells.
For each cell, the verifier computes exact integer striped-occupancy
weights and checks their complete mass. It averages the local
entrywise-alpha matrices in Arb and takes the 128-batch matrix power.
It then divides by a positive certified lower bound on both endpoint
histogram probabilities.

Take the entrywise maximum of certified regional upper endpoints over
all 121 cells. The resulting matrix U bounds every region and therefore
every ordered matching profile through U^256. Raw-operator monotonicity
extends this operator bound to the requested source interval. The source
sum still uses each actual size 4096,...,4108, including its own binomial
and density factors. No profile union or lower-endpoint source shortcut
is used.

The [directed verifier](theory/verify_small_state_dense_anchor.py) passes
at both 128 and 256 bits. Its log2 upper bound is approximately
-1396.44545192524, comfortably below the requested -100 threshold.
Both replays cover all 121 cells and all 13 source sizes. They take at
most 1.81 seconds and 61 MB; these are verification costs, not encoder
timings. Independent static review found no issue. The source SHA-256 is
`b0f0f3e01d6fe7cadb9a2b482791d4b7f9f35adf2c47b74f03d6be2618c968ec`.

Run separately under the same external resource limits as above:

```text
python -B theory/verify_small_state_dense_anchor.py --groups 4096 --end-groups 4108 --z 370391/5732639 --alpha 3/20 --precision 128 --target-bits 100
python -B theory/verify_small_state_dense_anchor.py --groups 4096 --end-groups 4108 --z 370391/5732639 --alpha 3/20 --precision 256 --target-bits 100
```

An independent occupancy calculation agrees within 1.2e-16 for three
proposal distributions across all 4845 sorted occupancies. An immediate
timeout correctly produces no passing claim. A complete failing control
at g=2048, z=1/5, alpha=1/10 returns a positive log2 bound near 25532.07 and
`passed=false`. The helper therefore distinguishes a complete failed
bound from incomplete computation and from successful verification.

## A directed intermediate interval with two-batch grouping

All messages using 2048..2051 actual outer groups have combined
setup-failure probability below 2^-100 at the same 35% symbol distance.
The encoder, setup distribution, and free terminal state are unchanged.
This interval does not fill the gaps between existing checkpoints.

The [grouped verifier](theory/verify_small_state_grouped_route.py) uses
bridge Sections 13 and 14. It groups two ordered batches before taking
the fractional power. Exact rational ratio bins yield positive upper
bounds on the sum over the intermediate state. This avoids enumerating
all 4845 squared local-profile pairs while retaining every occupancy.

Use z=23/100, alpha=19/200, proposal tilt 33/25, pseudocount 1/65536,
and log2 ratio-bin width 1/4. At anchor g2048, the complete partition has
157 cells: at most 128 bulk cells and 32 singleton cells at the high-f end.
Each cell uses the unconditioned two-batch expectation W, followed by
W^64 and one certified regional histogram divisor. Conditioning does
not preserve batch independence and does not raise that divisor to alpha.

The census evaluates all 24648 compatible ordered cell pairs and all
2048 run classes. Its exact integer counts partition all binom(8192,2048)
source sets. Averaging this raw anchor and using inclusion monotonicity
extends the bound to 2051. The source sum retains each actual size's
density and binomial factors outside the route average.

Both 128- and 256-bit directed replays pass. The log2 bound is approximately
-125.5887508163504 for the interval, and -423.3700050391628 at the anchor.
The maximum replay cost is 4.43 seconds and 79 MB, not an encoding time.
Independent static reviews found no correctness issue. The source hash is
`bf9f75a0f8f3984cd53edafbc84095b5af8d27ef6e27e8b551f40909ef5bc83e`.

Run these commands separately, with the same 55-second external timeout,
1-GiB address-space limit, and 20-GiB free-space floor:

```text
python -B theory/verify_small_state_grouped_route.py --groups 2048 --end-groups 2051 --z 23/100 --alpha 19/200 --tilt 33/25 --bin-width 1/4 --tail-counts 32 --precision 128 --target-bits 100
python -B theory/verify_small_state_grouped_route.py --groups 2048 --end-groups 2051 --z 23/100 --alpha 19/200 --tilt 33/25 --bin-width 1/4 --tail-counts 32 --precision 256 --target-bits 100
```

A complete negative control at z=23/100, alpha=1/10, tilt 4/3,
width 1, and no singleton tail cells returns a positive log2 bound
near 1237.8846 and `passed=false`. An immediate time cap returns neither
a passing claim nor a bound. A noncommuting g8 matrix fixture agrees
with direct ordered-profile enumeration within 2.6e-38 relative error.
The corresponding FLOAT implementation also passes a direct comparison
over all 23474025 ordered local-profile pairs; see the diagnostic note.

## Direct actual-size conditioning extends the intermediate interval

All messages on 2048..2073 actual groups now have combined setup-failure
probability below 2^-100 at 35% symbol distance. This includes and extends
the preceding four-size checkpoint; it is not an additional disjoint term.

The [direct interval verifier](theory/verify_small_state_grouped_interval.py)
reuses the same rational witnesses. It caches unconditioned two-batch
numerators separately for the even and odd size families. At each actual
size, it recomputes both endpoint histogram probabilities for every cell,
all compatible ordered cell pairs, and the exact source-run census.
The source prefactor also uses that actual size. This avoids paying the
source growth at a larger size while retaining only the smaller anchor's
route bound. Bridge Section 15 gives the positive-matrix calculation.

Both 128- and 256-bit replays pass, with log2 union bound approximately
-276.66625679167865. Each replay covers 26 sizes and 53,573 run classes;
each parity family starts with 157 complete cells. The maximum verification
cost is 10.06 seconds and 86 MB, not an encoder timing. The verifier hash is
`79644db7e509a017340918e847b0685f51fc8910c72edb84bcbbf7fd385711ff`.
Old verifiers and their checkpoints are unchanged.

Run separately under the same resource limits:

```text
python -B theory/verify_small_state_grouped_interval.py --groups 2048 --end-groups 2073 --precision 128 --target-bits 100
python -B theory/verify_small_state_grouped_interval.py --groups 2048 --end-groups 2073 --precision 256 --target-bits 100
```

Independent static review checked the conditioning, endpoint ordering,
complete partitions, and source count. Extending the same fixed proposal
families through 2074 produces a positive log2 bound near 43.0572 and
`passed=false`, despite complete coverage. An immediate timeout returns
no bound and no passing claim. This is still partial coverage, not an
all-source certificate.

## Next proof task

Fill source counts 129..2047, 2074..4095, and 4109..8192 before changing
the encoder. The new grouped gate retains intermediate-state sums instead
of paying for separate branches at every batch boundary.
The dense gate shows that the inherited fractional marker was part of
the apparent obstruction, not a necessary construction cost. Accurate
source counts closed the earlier sparse range. The bounded continuation driver
is [tune_small_state_route_range.py](theory/tune_small_state_route_range.py);
its floating output only proposes witnesses for the directed verifier.
Do not interpret an inherited marker's failure as a structural obstruction.

The shared positive triangular coefficient calculation removes the
per-profile duplication and already skips impossible single-batch
monomials. Its current API is deliberately capped at128; enlarging that
domain requires new bounded checks, not an inferred certificate. At
g>=256 replace the h=0 source correction, whose retained
fraction is no longer positive; the joint-weighting correction is already
proved. At g=8192 include the fully occupied cycle profile: the proper-set
edge exclusion used here does not apply there. The full union over all
8192 counts and an implementation timing remain outstanding.
