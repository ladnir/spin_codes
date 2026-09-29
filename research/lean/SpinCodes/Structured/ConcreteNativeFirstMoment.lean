import SpinCodes.Structured.ConcreteNativeFamily
import SpinCodes.FiniteConditionedMoment

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteRoute ConcreteRoutedFirstMoment

theorem prob_decidable {Ω : Type*} [Fintype Ω] (P : FinPMF Ω) (s : Ω → Prop)
    (d₁ d₂ : DecidablePred s) : @FinPMF.prob Ω _ P s d₁ = @FinPMF.prob Ω _ P s d₂ :=
  congrArg (fun D => @FinPMF.prob Ω _ P s D) (Subsingleton.elim d₁ d₂)

theorem qd_le_matrixBound (m : ℕ) (out : NativeSeed m)
    (message : Fin (Nsched m / 2) → ZMod 2)
    (hq0 : 0 < density (rows m out message)) (hq1 : density (rows m out message) < 1)
    (d : ℕ) {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (nativeSetup m).qd d ((nativeSetup m).outer out message) ≤
      matrixBound (rounds m) d z (rows m out message) := by
  simp only [nativeSetup, Setup.qd, ConcreteBinaryEncoding.wordToRows_rowsToWord, rounds]
  convert low_weight_probability (Lsched_pos m) (bsched_pos m)
    (streamWiring (native_round_divisibility m)) (rows m out message) hq0 hq1 d hz0 hz1 using 1
  exact prob_decidable _ _ _ _

/-- The concrete selected-family first moment; only selected nonzero words need density bounds. -/
theorem EZ_le_matrix_sum (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m Q : ℕ)
    (hdensity : ∀ out, G m out → ∀ message : Fin (Nsched m / 2) → ZMod 2,
      message ≠ 0 → 0 < density (rows m out message) ∧ density (rows m out message) < 1)
    {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (family G hG).EZ m Q ≤
      ((nativeSeedLaw m).condition (G m) (hG m)).expect
        (fun out => ∑ message ∈ (nonzeroMsgs (Fin (Nsched m / 2) → ZMod 2)).filter
          (fun message => occ m out message = Q),
          matrixBound (rounds m) (threshold m) z (rows m out message)) := by
  apply Setup.expect_ZQ_condition_le (nativeSetup m) (threshold m) (occ m) Q
    (G m) (hG m) (fun out message => matrixBound (rounds m) (threshold m) z (rows m out message))
  intro out hg message hm _
  have hd := hdensity out hg message (Finset.mem_filter.mp hm).2
  exact qd_le_matrixBound m out message hd.1 hd.2 (threshold m) hz0 hz1

end Spin.Structured.ConcreteNativeFamily
