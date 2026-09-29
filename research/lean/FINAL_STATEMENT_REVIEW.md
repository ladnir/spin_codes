# Final native statement review

Reviewed 2026-09-28 after both the actual native theorem/pin controller and the unconditional corollary controller published PASS. No statement mismatch, residual hypothesis, nonstandard axiom, or current hash mismatch was found. This was an independent read-only review of actual compiler output, source definitions, and verification records. It launched no compiler and does not claim a fresh replay of the whole project.

## What the compiler printed

`Spin.Structured.ConcreteNativeTheoremPin.actual_minimum_distance`, printed in `scripts/map_data/native_theorem_check_ConcreteNativeTheoremPin.log`, states that

\[
\Pr_{(\mathrm{out},\mathrm{inner})\sim P_m}
 [\operatorname{minimumDistance}(m,\mathrm{out},\mathrm{inner})
   \le \lfloor 0.11\,\operatorname{Nsched}(m)\rfloor]
\longrightarrow 0,
\]

where the printed law is exactly

```lean
(nativeSeedLaw m).prod
  (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))
```

The theorem has no hypothesis binders. In particular, it does not assume a spectrum envelope, selection probability, route identity, kernel norm, or dense certificate. `native_minimum_distance_failure_tendsto` supplies the proof; `native_probBad_tendsto` connects the same result to the established concrete first-moment family.

The second full pin, `Spin.Structured.ConcreteNativeTheoremPin.actual_rate`, quantifies only over `m : ℕ`, `out : NativeSeed m`, and `inner : InnerSeed m`. For every such realization it states

```lean
(Module.finrank (ZMod 2) (realizedCode m out inner) : ℝ) /
  (Nsched m : ℝ) = 1 / 2
```

This is exact rate for every native index and every seed, rather than an asymptotic rate or a property restricted to selected seeds.

## Meaning of the objects and event

The source definitions ground the printed probability in the actual construction. `ConcreteOuter.seedLaw` samples one pair of independent uniform permutations for the Golay–BAA constituent; `rowSupports` reuses this same seed for every row. `ConcreteRoute.seedLaw` supplies independent uniform permutations for each row and each transposed region. `ConcreteRoutedEncoder.experimentLaw` combines that route law with the product of the actual transvection laws, one per encoder round. The outer and inner seeds in the final product are independent. The final probability is unconditioned; the good-outer event is an internal proof device, not a restriction on the theorem's sample space.

`ConcreteNativeCodeword.codeword` is the emitted native binary word. `realizedCode` is the range of its proved linear encoder. Injectivity establishes dimension `Nsched m / 2`. The native length is positive and even, which justifies the rate division and its exact value. `minimumDistance` is the minimum nonzero weight in this realized code; `minimumDistance_le_iff_pair` proves its equivalence with the usual minimum Hamming distance over distinct emitted codewords. Thus the theorem is not merely about the weight of one message or a proxy statistic.

The threshold is the natural-number floor of the exact real rational `0.11 * Nsched m`. `threshold_floor` connects this to the native integer threshold, and `relative_distance_gt_iff_threshold` proves the exact equivalence with strict relative distance above `11/100`. No rounding discrepancy changes the failure event. The schedule is the actual native schedule indexed by `m`; the statements do not quantify over every requested block length.

## Unconditional consequences

The actual output in `encoder_consequence_check_ConcreteNativeDistanceConsequencesFinalPin.log` prints two further theorems with no failure-limit hypothesis:

- `ConcreteNativeFamily.relative_distance_success_tendsto`: under the same actual product law, the probability of `11/100 < minimumDistance / Nsched` tends to one.
- `ConcreteNativeFamily.eventually_exists_rate_half_distance_gt_eleven_percent`: for every sufficiently large native index, there exist outer and inner seeds with strict relative distance above `11/100`, exact rate `1/2`, and positive mass under that actual product law.

The helper theorems do explicitly accept a failure-limit hypothesis; the final wrappers discharge it using `native_minimum_distance_failure_tendsto`. Complementation is over the same probability law, and the existence proof extracts a positive-mass seed from a positive success probability. There is no deterministic search or encoding-complexity claim in these declarations.

## Evidence checked

The base report `native_theorem_verification.json` has SHA256
`6ff120a0ce7c66a136ea2a7a7235d32f468c5fe02c6ac10095c8c8c9e2d314ac`.
Both newly compiled modules exited successfully. All four recursive axiom outputs—the two base theorems and the two full pins—contain only `propext`, `Classical.choice`, and `Quot.sound`.

The review independently rehashed all 5215 current source/object pairs in that report, the fixed-range and mixed-dense input reports, the separate dependency-provenance report, the original `SpinCodes/Pin.lean`, and the recorded paper files. All matched. The 68 dependency timestamp inversions are backed by explicit paired compilation evidence, including five fresh isolated outputs that byte-match the live objects; they are not accepted solely from a cache inventory.

The corollary report `encoder_native_distance_consequences_verification.json` has SHA256
`a06a417300815c64e53e9201c00d3adcb1f5eac30685058f53a3abac0492a51a`.
All four owned modules were freshly checked. The review independently verified their current source/object hashes, their link to the above base report, and all seven recursive axiom outputs, which are standard-only.

The 5215-entry base inventory records the dependency snapshot used for assembly. It is not evidence that every imported source was freshly compiled in that run. The separate coverage review distinguishes prior individual checks from source-only and dependency-snapshot evidence.

Final completion update: `native_final_invariants_verification.json` passed with exit code zero and `ALL CHECKS PASS`. After that result, the independent read-only consolidated checker completed with **PASS: 17 gates, 5220 current recorded source/object pairs, zero issues, and no missing or invalid gates**. Its record is `scripts/map_data/final_closure_evidence.json`. The final stability pass revalidated the consumed reports, source/object pairs, original pin and papers, and required verification artifacts. All mandatory completion gates are now satisfied. A fresh whole-project replay remains a separate reproducibility exercise, not part of this closure claim.
