# Width-eight packet co-design

**Parked by decision on 2026-10-03.** Keep the certified four-bit design as
the main path. Width eight produced a complete K16 proof and modest small-size
gains, but no compelling improvement across the range. At K20 it is 3.7%
slower than the matched four-bit control, and K18/K20 proof coverage remains
incomplete. Preserve all sources, proof procedures, and performance records;
do not promote this construction or launch further proof sweeps by default.
The next exploration is a construction designed around carry-less polynomial
multiplication, rather than further tuning of this byte-packet branch.

Latest scaling study: [iteration7](iteration7/README.md) measures the same
byte construction at **0.373 ms for K18** and **3.389 ms for K20** after
size-specific page/store choices. The matched four-bit K20 control remains
faster at **3.268 ms**. K18's seven screened occupancies pass 40 bits at
10% distance; K20's q16 bound is still vacuous. Neither larger size has a
whole-code certificate. If this branch is resumed, K18 is the next proof
target; keep the older K20 baseline. [All performance observations](iteration7/PERFORMANCE.md) are retained.

Current optimized checkpoint: [iteration6](iteration6/README.md) measures
**93.895 us** at K16 with the same **68.893749-bit /10%-distance certificate**.
It simplifies the existing randomizer circuit without changing the code.
A new structured randomizer is also proved, but does not beat that kernel.
The [fifth iteration](iteration5/README.md) remains the frozen proof/control:
it closes the wider 256-to-512-bit outer with byte packets and a 24-bit state. Its complete
outward certificate gives **68.893749 bits at 10% distance** for K=2^16.
The faithful precomputed transpose measures **98.1625 us** on 128-bit elements,
against 99.591 us for the certified four-bit control in the same campaign.
An independent global replay checks every occupancy and the exact union.
This is an isolated research checkpoint; production code and the certified
four-bit baselines remain unchanged. The earlier
[floating proposal](iteration4/README.md) is retained as its search history.

## Earlier checkpoints

Status: implementation study and second proof checkpoint completed, 2026-10-02.
The byte-native K16 kernel is faster, but its tested whole-code bound is
vacuous at 10% distance. It is not a replacement for the
[certified baselines](../../CERTIFIED_BASELINES.md). K18 small-state tuning
remains parked. No production kernel or frozen checkpoint was changed.

Width eight means permuting **eight binary code coordinates together**.
It does not change the 128-bit payload type and is not a permutation of
only eight positions. A candidate can retain the present outer code while
grouping its randomized output into bytes rather than four-bit packets.

## What the first study established

The isolated [implementation](../../../spin/experiments/packet8_codesign/README.md)
keeps byte-packed values through the inner, routing, and outer field operations.
It uses the retained small RS16 outer, a genuine 32-region byte route, and
the two-byte state specified below. Its transpose uses binary adjoints of
field multiplication, not ordinary field multiplication in their place.

| K=65,536, 128-bit elements | Precomputed transpose | Proof status at 10% |
|---|---:|---|
| Certified four-bit control | 99.341 us | Complete 62.04-bit certificate |
| Byte-native, shared outer | 91.120 us | q1 passes; middle bound fails |

These are medians of eight run medians per variant from one serial ABBA
holdout on Peach. Four new seeds, 2,001 calls per run, normal pages, and
no stage clocks were used. Setup and allocation are excluded. The measured
time reduction is 8.3%; it is not a certified-code speedup.

Stage diagnostics put the byte-native inner plus route near 44.5 us and
its outer near 46 us. Write-prefetch reduces route-only time from roughly
43 to 40 us but does not improve the full call. A half-payload outer
reduces temporary storage and register spills, but its extra masked stores
raise whole-call time to roughly 102 us. Neither variant is selected.

All six native test configurations pass, including scalar/SIMD equality,
forward/transpose adjoints, alignment, guards, and in-place operation.
The 27 portable algebra and proof tests also pass. The frozen K16 checkpoint's
220 source pins and archived receipt were checked again after this study.

## Why this is not ready for a certificate

[PROOF.md](PROOF.md) specifies the maps, exact finite formulas, and numerical
scope. [AUDIT.md](AUDIT.md) gives independent algebra and bound checks.
The baseline has a roughly 55.44-bit floating occupancy-one bound. However,
a finer tilt search over occupancies 80 through 160 still has a weakest
margin of **-4023.60 bits**, at occupancy 126. Negative values mean a
vacuous upper bound, not an observed low-distance word.

Two actual map changes were screened, not implemented:

- Scaling only the expansion raises its minimum weight from 8 to 17,
  but improves the tested middle bound by only a few bits.
- Reducing the physical step from 64 to 32 coordinates retains the byte
  route and 16-bit state. Its tested middle bound remains near -4099 bits.

The dominant bound contributions include excursions that leave zero and
return after occupied steps. The always-zero trajectory is not the main
obstruction. An independent pointwise comparison also limits slack in the
entire local transition bound to about 50 bits at tilt 0.4. Thus recovering
the missing thousands of bits requires more than refining that local bound.
This comparison says nothing about slack in the outer counting envelope.

An explicitly counterfactual diagnostic reduces only return probabilities,
while retaining all actual 16-bit birth and emission distributions. Even
a 32-bit return denominator leaves the sampled middle bound negative.
It is neither a realizable larger-state construction nor a certificate;
it shows why increasing a denominator alone is not a sufficient design argument.

## Next decision

The [third iteration](iteration3/README.md) tests conditioned
trajectory diagnostics, a three-packet routing condition, and a fractional
moment over routing. The completed checkpoints below remain unchanged.

Keep the certified four-bit results as the defaults. Do not launch a full
outward replay of these negative byte-native bounds.

The follow-up tested a real 24-bit state, an outer with 64 byte regions,
and two input-refresh variants. None closes the tested middle bound.
The [larger-outer study](wider_outer/README.md) gains about 326 bits at a
matched active-group fraction, but still misses by roughly 3,697 bits.
There is no evidence here to justify their added implementation cost.

The stronger lead is now a proof change. A trajectory diagnostic found that
the failing comparison puts about 3,568 active bytes into only 549 physical
steps. Most remaining steps have both zero input and zero state. This is
rare routing concentration, shared by many message assignments, rather
than uniformly poor state mixing.

[The new probability split](PROOF.md#separating-routing-failures-from-message-counting)
bounds poor routing over subsets of groups before counting message values.
For example, the route-count helper proves an exact 60-bit bound for the
event that some 119-group subset hits fewer than 1,483 physical steps.
A second exact condition also limits excess packets beyond two per step.
The [new checkpoint](route_conditioning/README.md) records their use in the
message bound:

| q=119 proposal | Baseline 16-bit state | Actual 24-bit state |
|---|---:|---:|
| Condition on occupied-step count | -2,255.30 bits | -2,225.73 bits |
| Condition on first-two-packets count | -617.49 bits | -258.39 bits |

These are selected floating message bounds, not certificates or optimized
comparisons across all map parameters. The 24-bit row uses a different,
unimplemented inner. Combining both route conditions on the baseline did
not improve the tested bound; the best remained the second condition alone.

The unchanged baseline therefore gains thousands of bits in its proof
screen, but its middle range remains open. Next, refine the routing
statistics before spending more implementation work. For the 24-bit
fallback, counting the first three packets per step is a concrete next
gate because every three-byte feedback restriction has full rank.
Neither direction justifies promotion or a full outward replay yet.

The follow-up's 58 tests pass. Twenty new receipts have 184 matching source
pins, and the frozen K16 checkpoint still has all 220 pins and its receipt
intact. No encoder benchmark or production edit was needed for this proof
exploration; the prior 91.120 us byte-kernel timing remains the only matched
measurement for that candidate.

The follow-up studies are preserved separately:

- [Actual 24-bit maps and trajectory diagnostics](larger_state/PROOF.md).
- [Larger outer](wider_outer/README.md).
- [Input-refresh orders](input_refresh/README.md).
- [Exact route witnesses and conditioned message screens](route_conditioning/README.md).

The roughly 8 us implementation saving is a measured budget for any added
work, not a promise that a stronger construction will retain it. K18 and K20
should wait until there is a credible K16 proof/performance combination.

## Study order and preserved controls

The first study followed this bounded co-design sequence:

1. Specify the grouping, local maps, and shuffle distribution. Check input
   invertibility, single-packet feedback ranks, and obvious cancellation
   patterns before timing a candidate.
2. Implement an isolated precomputed transpose kernel. Start at K=2^16
   with 128-bit elements and the existing small RS outer. Keep the current
   t64 step size initially so the first comparison isolates packet grouping.
3. Check the new code against its own scalar oracle and forward/transpose
   adjoint. Compare runtime with the frozen four-bit baseline in serial,
   reversed-order, multi-seed A/B runs. Different constructions need not
   produce equal output checksums.
4. Measure whole encoding plus separate inner/routing and outer diagnostics.
   Count packing conversions, route-index accesses, and scratch traffic.
   Eight-bit packets reduce the number of routed packets but not the number
   of payload bytes; they do not promise a twofold speedup.
5. If the end-to-end gain is material, run sparse and middle-occupancy proof
   screens for the exact measured maps. Let that result guide feedback or
   local mixing changes. Only then attempt the complete outward union.

K18 and K20 can receive bounded scaling probes after the K16 kernel is
credible. These have different accepted outer/inner baselines, so neither
timings nor certificates should be extrapolated from K16.

## Prior evidence to reuse, not overwrite

The retained [packet8 experiment](../k16_codesign_100us/packet8/README.md)
already considered the small outer with 32 eight-bit-packet regions and a
t64/s16 inner. Its source is part of the frozen K16 checkpoint; new work
belongs here rather than in that directory.

Naively grouping adjacent coordinates gives feedback rank seven on a byte.
Four coordinate swaps repair every single-packet restriction to rank eight.
That version proves a **52.58-bit occupancy-one bound** at K16 and 10%
distance, but its complete bound remains open. At occupancy 128, the tested
bound is vacuous by thousands of bits. An independent outward replay confirms
that selected bound. This is a proof gap, not a low-distance counterexample.

Two-packet restrictions still have ranks 11 through 13 in the corrected map.
A globally shuffled coordinate layout improves those ranks substantially
but barely improves the tested bound. Therefore full single-packet rank
alone, or a more expensive coordinate shuffle alone, is not enough evidence
to invest in a full proof.

No encoder benchmark accompanied that old packet8 screen. The missing first
piece is its actual speed potential: byte-native packing and routing may
save work, but any required local mixer or conversion must be timed too.
Compare a rank-corrected control with a deliberately byte-friendly design;
do not treat a routing-only optimistic timing as a complete encoder result.

The [byte-native inner proposal](../k16_outer_routing/QUOTIENT_RANDOMIZATION.md#one-bounded-byte-native-inner-screen)
is another starting point. For eight distinct field points `alpha_i` in
GF256, it uses two byte-valued state coordinates and maps

```
A(a,b)_i = a + alpha_i*b
C(x) = (sum_i x_i, sum_i alpha_i*x_i).
```

The prior screen uses points 0 through 7, original-order emission, and a
fresh transitive map on the nonzero sixteen-bit states at each step.
Every one-byte restriction has rank eight, and every two-byte restriction
has rank sixteen. The retained RS16 outer gives a 55.44-bit floating q1
screen, but neither q2 nor the tail nor implementation timing is available.
The [instruction-level design note](../../../spin/experiments/k16_outer_routing/search/DESIGN_SPACE.md)
explains how the inner and outer could share the byte-packed representation.

By contrast, adjacent pairing of nibble packets is not automatically a
good proof design even with the newer Moore feedback. Its two natural
bytes contain points `{a,a+1,b,b+1}`. Equal nonzero nibble labels cancel all
four checks `(1,h,h^2,h^4)`. Any three points have rank twelve, so the
two-byte restriction has exactly rank twelve. Its zero-feedback probability
is `15/255^2 = 1/4335` for uniform active byte labels. This is a local
collision calculation, not a whole-code failure bound.

For the first performance control, a new byte route can use one address
and two adjacent 64-byte stores per packet. Generate a genuine 32-region
route; do not join independently routed old nibble destinations after setup.
The larger opportunity is keeping that packed format through the inner
and outer. Eight 128-bit elements still occupy 128 bytes, so neither
compulsory payload traffic nor the number of vector stores is halved.

The initial proof checks should retain the eight-bit label distribution
and the dependence introduced by grouping. Existing four-bit local census
tables and complete certificates do not transfer by relabeling metadata.
