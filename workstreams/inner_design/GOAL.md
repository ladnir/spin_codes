# Two-part inner-design goal

## Part 1: equivalent implementations

Preserve the exact quarter-rate encoder for every existing setup and input.
Search state bases that jointly reduce expansion and feedback costs; conjugate
every field transition. Check identities symbolically and against the original
dense oracle. Benchmark full in-place transposed encoding on Peach, serially.
The supported baseline remains unchanged unless a gain survives repeated,
alternating, uninstrumented comparisons.

The practical stopping gate is a completed comparison of identity/grouped,
feedback-reduced, emission-reduced, and jointly searched bases, followed by
independent search restarts. Stop this tier when additional candidates cease
to improve full-encoder time materially. This is not a proof of global optimality.

## Part 2: proof-guided inner design

After that gate, separate the proof's essential requirements from conveniences
of the current construction. Audit expansion spectra, kernel events, state
mixing, cancellation/termination probabilities, region dependence, dense
occupation coverage, and deterministic low-weight obstructions.

Select cheaper candidates using both measured implementation costs and proof
quality. Derive sufficient conditions, test failure modes, and try to close a
full certificate for the strongest viable candidate. Record proved statements,
numerical evidence, counterexamples, and unresolved obligations separately.
No candidate replaces the baseline on gate count or heuristic margin alone.

## Current status

Part 1 practical gate closed: 60 annealing restarts / 3.9 million moves,
nine compiled bases checked against the original dense oracle, and four
finalists compared to the unchanged baseline in three serial 101-trial runs.
The best finalist was indistinguishable from baseline (17.587 vs 17.596 ms).
No candidate adopted; this is not an exhaustive optimality claim.

Part 2 has a successful candidate: one random rank-one update, combined with
the exactly audited balanced-image map at t=128,s=19. The image weights lie in
[48,80], and the kernel minimum weight is five. The proof tracks weight-biased
states, lazy cancellation, zero-syndrome inputs, dense occupations, and
unflushed termination explicitly. The original-map mixer remains a separate
experiment without a full certificate.

The new map's 256-bit outward certificate covers every nonzero occupancy at
K=2^20. Full margins are 41.048169 bits at distance 16.5% and 30.033492 bits at
distance 19%. The 512-bit replay passes every saved bound. Neither the original
map nor its production implementation was overwritten.

Both isolated final-map kernels pass the dense oracle at K=2^16, 2^18, and
2^20. Three alternating serial 101-trial comparisons give 17.135 ms (sparse)
and 17.145 ms (fixed masked), versus 17.847 ms for the unchanged baseline.
Both save about 4% online time and 2.125 MiB of retained setup; the difference
between the two candidates is not significant in these measurements.

The two-part investigation is complete. The consolidated deliverable is
[BALANCED_RESULT.md](BALANCED_RESULT.md), including proof scope, failure modes,
timings, and reproduction commands. `verify_balanced_artifact.py` binds the
certificate to the exact map, tested source, outer span, and timing receipts.
Integration into the production default and paper is the recommended next step.
