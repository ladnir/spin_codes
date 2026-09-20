# Larger SPIN redesign exploration

This completed, bounded workstream evaluates three directions: factored macroblocks, accumulator
cascades, and changes to routing/locality. It keeps the quarter-rate BCH
[128,32,32] outer and K=2^20. No supported source or existing certificate is
modified. See `DESIGN_SPACE.md` for exact interfaces and proof obligations.

## Findings

The bounded search has not found a replacement that beats the asymmetric
t=128,s=19 control, which remains around 16.0–16.2 ms on this host.

| Direction | Current result | Decision |
|---|---|---|
| Factored t=256, two mixer rounds | Correct implementation; about 16.74 ms | No performance reason to extend to a full certificate |
| Factored t=512, two mixer rounds | Correct implementation; about 17.13 ms | Same |
| Exact-map direct gather | Same output, about 32.7 ms; software prefetch did not help | Retain the existing bucketed route |
| Two accumulator passes | 59.02 ms direct; 45.49 ms after bucketed-routing optimization | Not competitive in either tested implementation |
| Three accumulator passes | About 82 ms with precomputed direct routing | Not competitive in this implementation |
| Purely chunk-local permutations plus short state | Exact 6.01%-weight counterexample | Cannot meet either distance target |

Every reported timing now has three separate 101-call runs. The bucketed
two-pass accumulator uses 64 MiB of retained permutation schedules and
128 MiB of workspace. All candidate and diagnostic correctness checks pass;
the three-repeat output hashes agree, and both routing implementations retain
the control encoder's exact output hash. The existing balanced certificate's
integrity check also passes unchanged.

These are bounded experiments, not impossibility results for all macroblock
circuits, accumulator implementations, or locality-aware codes. In particular,
no exhaustive joint search over A, B, mixer, and route has been attempted.

## What the diagnostics say

Removing the inner but retaining unrolled routing took approximately 15.7 ms,
versus roughly 16.0 ms for the full encoder. A scalar-loop version was slightly
slower than the full encoder, demonstrating scheduling effects. Contiguous
outer-only work took about 5.6 ms.

These invalid-distance diagnostics are not runtime lower bounds or an additive
profile. They suggest that eliminating inner arithmetic alone is not a likely
route to a 20–30% win with the current data movement. Direct gathering moves
less intermediate data but roughly doubles runtime, so bytes moved alone are
also insufficient to predict performance.

## Proof status

`MACRO_Q1.json` reconstructs exact expansion and feedback spectra, checks
rank/distinct syndromes and low-weight kernel counts, and evaluates the
single-active-row bound. Its tilt is selected separately for each outer weight,
so its control numbers need not equal the earlier fixed-tilt Q=1 diagnostic.
All listed Q=1 results are binary64 calculations, not full certificates.

At t=256 and two rounds, the Q=1 screen gives about 41.0645 bits at 16.5%
distance and 30.0488 bits at 19%. At t=512, it gives about 41.0053 and
30.0043 bits. Other occupancies remain unchecked for these new maps. The
performance results do not justify spending the current bounded budget on
that remaining certificate work.

For accumulator cascades, an exact input/output weight transition and its
first-moment composition are documented. Exhaustive small-instance tests check
the formula, but no large finite-length margin is claimed.

For cache-local permutations, the obstruction is exact. A 32-dimensional outer
row has a nonzero message in the kernel of its 19-bit outgoing-state map.
When all of that row's coordinates stay in a chunk, that message cannot affect
later chunks. `LOCALITY_COUNTEREXAMPLE.json` records a concrete message whose
full output has weight 251,944 out of 4,194,304, with outgoing state exactly
zero. This refutes the restricted-locality proposal, not the existing SPIN
permutation ensemble. Larger boundary communication or coordinate exchange
between chunks can avoid this particular obstruction, at additional cost.

The design notes separate identities and the locality obstruction from
numerical screens, measured performance, and unproved distance claims.

## Files and reproduction

- `screen_macro.py`: five map candidates, exact spectra/low-kernel audits,
  and a bounded Q=1 tilt grid with one/two rounds and an exact-refresh control.
- `generate_macro.py`, `MacroInner.h`: two isolated C++ macroblock kernels.
- `generate_routing.py`: exact-map gather implementations with zero/64-position
  prefetch distance.
- `generate_headroom.py`: no-inner and outer-only diagnostics; not valid SPIN
  distance candidates. Independent references describe their actual maps.
- `accumulator_benchmark.cpp`: two/three-pass BCH-accumulator implementations,
  precomputed permutation schedules, optional bucketed routing, and correctness
  checks outside the timed region.
- `test_designs.py`: factored-map adjoints, mixer order, accumulator adjoints,
  exhaustive accumulator transition counts, and spectra beyond 128 positions.
- `locality_witness.py`: exact local-message cancellation witness.
- `measurements/`: raw process medians/quantiles, correctness logs, binary hashes.
  `redesign-measured-sources.sha256` binds the remote compiled sources to the
  local source files, including the unchanged benchmark and outer circuit.
- `summarize.py`: source-binding checks and measurement summary; not a distance
  certificate verifier. `--partial` reports missing measurement receipts.

From the repository root, with `../requirements.txt` installed:

```text
python workstreams/inner_design/redesign/screen_macro.py
python workstreams/inner_design/redesign/locality_witness.py
python -m unittest discover -s workstreams/inner_design/redesign -p "test_*.py" -v
python workstreams/inner_design/redesign/generate_macro.py
python workstreams/inner_design/redesign/generate_routing.py
python workstreams/inner_design/redesign/generate_headroom.py
```

Use the existing `workstreams/inner_design/CMakeLists.txt` to build the generated
`basis_redesign_*_bench` and corresponding `_test` targets, plus the existing
`basis_asymmetric_greedy3_2_sparse_bench` control. For accumulator prototypes:

```text
cmake -S workstreams/inner_design/redesign -B redesign-build -DCMAKE_BUILD_TYPE=Release
cmake --build redesign-build -j3
```

Run `run_headroom.sh`, `run_macro.sh`, `run_routing.sh`, `run_accumulator.sh`,
and `run_accumulator_bucketed.sh` sequentially, each with the repository/build
root as its argument. Never run them concurrently. They use CPU 15 and the
shared lock. Builds target znver4; adapt these settings for another machine.
Copy each run's `redesign-*` timing, correctness, and binary-hash files into
`measurements/`. Finish with `python workstreams/inner_design/redesign/summarize.py`.

Measurements use three warmups and 101 calls per process, in place, with setup
and copying excluded. Repeated calls intentionally reuse the modified input.
Controls are interleaved with candidates. The accumulator workload includes
the same BCH outer and N-to-K geometry, not merely an inner-only microbenchmark.

## Recommendation

Keep the asymmetric 128-position winner and finish its outward certificate
before promoting it. A future large-redesign effort should focus on measured
data-movement improvements that preserve sufficient cross-chunk communication.
Neither another inner-only arithmetic sweep nor short-state-only cache
localization currently has a convincing path to the desired large speedup.
