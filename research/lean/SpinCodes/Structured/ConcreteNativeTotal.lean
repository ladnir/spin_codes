import SpinCodes.Structured.ConcreteNativeSelection
import SpinCodes.FinitePositiveSelection

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Filter ConcreteOuter

/-- The raw paper event when non-null; the whole space otherwise. This only makes
conditioning total. It changes neither the encoder nor its bad-distance event. -/
def selectedGood (m : ℕ) : NativeSeed m → Prop :=
  (nativeSeedLaw m).positiveSelection (nativeGood m)

theorem selectedGood_pos (m : ℕ) : 0 < (nativeSeedLaw m).prob (selectedGood m) :=
  FinPMF.positiveSelection_pos _ _

def concreteFamily : Family := family selectedGood selectedGood_pos

/-- The probability being bounded is independent of the auxiliary selection event. -/
theorem concrete_probBad_eq (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m : ℕ) :
    concreteFamily.probBad m = (family G hG).probBad m := rfl

/-- An explicit event for the actual sampled permutations and transvections. -/
theorem concrete_probBad_explicit (m : ℕ) : concreteFamily.probBad m =
    ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => ∃ x : Fin (Nsched m / 2) → ZMod 2, x ≠ 0 ∧
        ConcreteRoutedEncoder.weight (ConcreteRoute.streamWiring (native_round_divisibility m))
          (rows m ω.1 x) ω.2 ≤ threshold m) := by
  change (nativeSetup m).joint.prob (fun ω => 1 ≤ (nativeSetup m).Z (threshold m) ω) = _
  apply FinPMF.prob_congr
  intro ω
  rw [one_le_Z_iff]
  simp only [setup_inner_weight]

theorem selectedGood_eq (m : ℕ) (h : 0 < (nativeSeedLaw m).prob (nativeGood m)) :
    selectedGood m = nativeGood m := FinPMF.positiveSelection_eq _ _ h

theorem selectedGood_eventually
    (h : Tendsto (fun m => (nativeSeedLaw m).prob (fun ω => ¬ nativeGood m ω)) atTop (nhds 0)) :
    ∀ᶠ m in atTop, selectedGood m = nativeGood m :=
  FinPMF.positiveSelection_eq_eventually nativeSeedLaw nativeGood h

theorem concrete_probNotGood_tendsto
    (h : Tendsto (fun m => (nativeSeedLaw m).prob (fun ω => ¬ nativeGood m ω)) atTop (nhds 0)) :
    Tendsto concreteFamily.probNotGood atTop (nhds 0) :=
  FinPMF.positiveSelection_compl_tendsto nativeSeedLaw nativeGood h

end Spin.Structured.ConcreteNativeFamily
