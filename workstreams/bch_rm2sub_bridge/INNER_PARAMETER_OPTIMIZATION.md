# BCH-256 inner optimization: first bounded pass

Follow-up: [CERTIFICATE_SEARCH.md](CERTIFICATE_SEARCH.md) describes the bounded,
resumable witness-search controller. Its target is a 40-bit total union margin;
the 80-bit retained anchor bounds below are not per-occupancy requirements.

The curve-estimator work is deferred. The immediate objective is a cheaper
inner map with at least 40 bits of distance/setup margin, beginning at K=2^20.
The selected t64_s20 map remains the full certified baseline at 50.487298 bits.
The setup and code are fixed as in [CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md).
All new calculations use H=floor(2K/10)=209715 and fresh independent multipliers.

The initial candidates are the previously snapshotted t64_s16 and t128_s19 maps.
Their Q1 diagnostics are approximately 50.049 and 50.448 bits, respectively.
Matching t and s to another workstream's nested maps does not identify the same instance.
No runtime improvement has been measured; fewer state updates or state bits
are engineering motivations, not demonstrated encoder speedups.

| Selected map | Minimum nonzero A weight | Minimum nonzero kernel weight |
|---|---:|---:|
| t64_s20 | 16 | 8 |
| t64_s16 | 16 | 6 |
| t128_s19 | 48 | 6 |

The cheaper candidates change the kernel as well as the nominal state size.
The endpoint comparison does not isolate the effect of either change alone.

## Why the screening method matters

`inner_candidate_screen.py` uses the incoming activation-aware four-state transfer.
It applies our square-root region-density comparison and a concave envelope
of row counting costs. It samples occupancies, tilts, and mean active densities.
It is a binary64 screen, with no coverage between samples.

This cheap comparison gives a weak bound at Q=8192 even for the proved t64_s20
baseline. Therefore its weak candidate results do not identify an actual
obstruction. For t128_s19 it also gives a weak bound at Q=2620.
Those results motivate a tighter calculation, not rejection of the map.

`inner_candidate_boxes.py` then imports the finite-theory engine's fixed-type
calculation. The initial 63-node baseline run gives a weak bound too. Its
four nonzero weight categories cover all dense types, but the coarse box
budget and band grouping are insufficient for a useful bound here.
This wrapper checks the integer type count and pairwise disjointness exactly.
Its numerical bounds remain binary64 diagnostics. No expansion of this weak
box search is justified merely by the existence of a complete cover.

The productive fallback is `inner_candidate_fixed_weight.py`. It reuses our
frozen three-state, activation-aware transfer and exact fixed-weight region
polynomials. The wrapper supplies the candidate map and message length explicitly.
It chooses Bernoulli witnesses in binary64, then evaluates the adaptive
recurrence outward with 256-bit Arb and directed binary64 arithmetic.
Replay uses 512-bit Arb and the saved witness, without optimization.

This is the established inequality from [LARGER_STATE_GAP_CLOSURE.md](LARGER_STATE_GAP_CLOSURE.md),
not a transfer of numerical bounds from the baseline map. Map spectra are
enumerated and kernel spectra checked by the existing exact loaders.
The existing certified BCH shell caps are unchanged.

## Initial anchor results

An anchor bounds the expected bad-message count at exactly one occupancy Q.
It does not cover neighboring occupancies. Passing receipts retain only
80 bits of margin even when the computed bound is much smaller.

| Map | Q | Tilt log(lambda) | Outward anchor result |
|---|---:|---:|---|
| t64_s16 | 512 | -2.6 | Upper bound at most 2^-80 |
| t64_s16 | 8192 | 0.6 | Weak bound; retained exponent +57180 |
| t64_s16 | 8192 | 0.4 | Weak bound; retained exponent +55209 |
| t128_s19 | 2620 | -0.5 | Upper bound at most 2^-80 |
| t128_s19 | 8192 | 0.6 | Weak bound; retained exponent +270009 |
| t128_s19 | 8192 | 0.0 | Weak bound; retained exponent +240276 |

The t128_s19 Q2620 result is particularly useful: the exact fixed-weight
calculation resolves a point where the cheaper sampled bound was weak.
The old t128_s15 first-moment obstruction at that occupancy does not transfer
to the selected larger-state map.

Two endpoint tilts were tried for each candidate. The adjustments improved
the bounds but did not approach a useful endpoint result. Neither these upper
bounds nor the cheap screen prove a first-moment obstruction, actual setup
failure, or inability to certify another witness. No full optimized
configuration has been established by this pass.

Both useful anchors passed 512-bit replay without reoptimizing their witnesses.
Twelve tests passed, covering the new screens and receipt authentication,
exact integer type coverage, polynomial region arithmetic, the constant-row
band, and dyadic rounding. The same arithmetic formulas are used by producer
and higher-precision replay; the existing exact toy tests check those formulas
separately. This is not an independent proof of every mathematical reduction.

Receipts are local under `generated/inner_<map>_q<Q>_v1.json`, with independent
precision replay receipts named `*_replay.json`. The sampled screens and
coarse type-box run are retained separately; none contributes to a certificate.
All producers and outputs are versioned separately from frozen proof files.

## Next decisions

The endpoint, not Q1, is the immediate optimization bottleneck. First identify
the band and state trajectory driving the weak endpoint bounds, then apply
stronger fixed-type or joint-witness methods before expanding the occupancy sweep.
The small-BCH workstream has such refiners; its completed results do not
automatically supply witnesses for our BCH-256 caps and selected maps.

One concrete diagnostic is available without a spectrum model. At Q=8192,
the input weight is at least 38*8192=311296, above the output cutoff 209715.
A trajectory that stays in state zero throughout cannot be bad at this endpoint.
This deterministic observation does not eliminate trajectories that activate
and later return to zero, and it does not by itself improve the saved bound.
Check whether the relaxed counting measure spends material mass on such
excluded trajectories before investing in a new reduction.

After the endpoint is controlled, resolve the remaining occupancies, including
Q2 and the sparse-to-dense transition. A finite list of passing anchors is
not sufficient: every integer Q must enter the final union.

Keep t64_s16 as an alternative, but first improve its endpoint witness or
obtain stronger fixed-type bounds. If that remains costly, an intermediate
state size is a separate candidate requiring a pinned map and fresh checks.
Only after full closure should implementation benchmarks rank the choices.
Benchmarks and numerical exploration are run serially.
