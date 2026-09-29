import SpinCodes.Structured.ConcreteNativeTheorem
import SpinCodes.Structured.ConcreteNativeDistanceConsequencesFinal

/-!
# Native Structured SPIN: distance and rate

This is the public import for the completed native-length result. The declarations
keep their original names in `Spin.Structured.ConcreteNativeFamily`; this module
introduces no replacement definitions or additional hypotheses.

## Construction

At index `m`, `Spin.Structured.Nsched m` is the block length. One outer seed is
sampled from `Spin.Structured.ConcreteOuter.nativeSeedLaw m` and shared by all
rows. Independently, the actual routing permutations and IMT transvections are
sampled from `Spin.Structured.ConcreteRoutedEncoder.experimentLaw`. The image of
the resulting binary linear encoder is `realizedCode m out inner`, and
`minimumDistance m out inner` is its minimum nonzero Hamming weight. The theorem
`minimumDistance_le_iff_pair` identifies this with pairwise Hamming distance.

## Main declarations

All four names below belong to `Spin.Structured.ConcreteNativeFamily`.

* `native_minimum_distance_failure_tendsto`: under the original product seed
  law, the probability of distance at most `floor (0.11 * Nsched m)` tends to zero.
* `realizedCode_rate`: every realization has exact rate `1 / 2`.
* `relative_distance_success_tendsto`: the probability of relative distance
  strictly greater than `11 / 100` tends to one.
* `eventually_exists_rate_half_distance_gt_eleven_percent`: every sufficiently
  large native index admits such a rate-half code with positive setup probability.

The selection event is internal to the proof. None of these final declarations
assumes successful selection, a continuum estimate, or a numerical certificate.
`ConcreteNativeTheorem` supplies the fixed-occupation and dense-rate arguments
to the final assembly; the sparse and selection bounds are already discharged
in that assembly.

## Reading and checking

Start with `PROOF_GUIDE.md` for the dependency structure and `NATIVE_RESULT.md`
for the exact claim. `SpinCodes.NativePin` prints the four full statements and
their recursive axiom dependencies. `FINAL_REPRODUCTION.md` distinguishes checks
using existing dependency objects from a complete project-source replay.

This import concerns native-length distance and rate. It makes no claim about
encoder operation counts or a complete arbitrary requested-length wrapper.
-/
