import SpinCodes.Structured.ConcreteFourierSpectrum

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteMaps ConcreteEncoder ConcreteScalar Spin.Imt
open scoped symmDiff

def targetCap (β z : ℝ) (d : ℕ) : ℝ :=
  scalarF β z d * (spectrumCap d (ratio0 β z) (ratio1 β z)/2 + 1/(2*524287))

theorem actLaw_indicator_le {q : State} (hq : q ≠ ∅) (r s : State) :
    actLaw q (r ∆ s) ≤ (if s = r ∆ q then (1 : ℝ) else 0)/2 + 1/(2*524287) := by
  have he : r ∆ s = q ↔ s = r ∆ q := by
    constructor <;> intro h
    · have hh := congrArg (fun a => r ∆ a) h
      simpa only [symmDiff_symmDiff_cancel_left] using hh
    · rw [h, symmDiff_symmDiff_cancel_left]
  rw [actLaw, Spin.prob_act hq, nonzeroStates_card_actual]
  simp only [he, Nat.cast_ofNat]
  split_ifs <;> norm_num

theorem bernoulliRow_expect (β z : ℝ) (hb0 : 0 ≤ β) (hb1 : β ≤ 1) (q r : State) :
    bernoulliRow β z q r =
      (poissonBinom (fun _ : Fin 128 => β) (fun _ => hb0) (fun _ => hb1)).expect
        (fun X => emitted z q X * actLaw q (r ∆ Cset X)) := by
  simp only [bernoulliRow, FinPMF.expect, poissonBinom_const_apply, mul_assoc]

theorem expect_two_parts {Ω : Type*} [Fintype Ω] (P : FinPMF Ω)
    (f g : Ω → ℝ) (a b : ℝ) :
    P.expect (fun X => f X*(g X/a+1/b)) =
      P.expect (fun X => f X*g X)/a + P.expect f/b := by
  simp only [FinPMF.expect, mul_add, ← mul_div_assoc, mul_one, sum_add_distrib, sum_div, mul_assoc]

/-- The lazy branch uses the tilted syndrome cap; the refresh branch costs 1/(2M). -/
theorem bernoulliRow_target_le {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z)
    {q : State} (hq : q ≠ ∅) (r : State) :
    bernoulliRow β z q r ≤ targetCap β z (Aset q).card := by
  let P := poissonBinom (fun _ : Fin 128 => β) (fun _ => hb0.le) (fun _ => hb1.le)
  have he := emission_expect β z hb0.le hb1.le q
  have hw := weighted_syndrome_spectrum_le hb0 hb1 hz (Aset q) (r ∆ q)
  rw [bernoulliRow_expect β z hb0.le hb1.le]
  change P.expect (fun X => emitted z q X * actLaw q (r ∆ Cset X)) ≤ _
  calc
    _ ≤ P.expect (fun X => emitted z q X *
        ((if Cset X = r ∆ q then (1 : ℝ) else 0)/2 + 1/(2*524287))) := by
      apply FinPMF.expect_mono
      intro X
      exact mul_le_mul_of_nonneg_left (actLaw_indicator_le hq r (Cset X))
        (emitted_nonneg hz.le q X)
    _ = P.expect (fun X => emitted z q X * (if Cset X = r ∆ q then 1 else 0))/2 +
        P.expect (emitted z q)/(2*524287) := expect_two_parts _ _ _ _ _
    _ ≤ scalarF β z (Aset q).card*spectrumCap (Aset q).card (ratio0 β z) (ratio1 β z)/2 +
        scalarF β z (Aset q).card/(2*524287) := by
      apply add_le_add
      · apply div_le_div_of_nonneg_right _ (by norm_num)
        exact hw
      · exact (div_le_div_of_nonneg_right he.le (by norm_num))
    _ = targetCap β z (Aset q).card := by unfold targetCap; ring

end Spin.Structured.ConcreteFourier
