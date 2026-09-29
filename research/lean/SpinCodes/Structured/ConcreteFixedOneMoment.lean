import SpinCodes.Structured.ConcreteFixedOneBound

/-! Unconditional actual one-mark region moment bound. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder ConcreteMarked Filter

theorem actual_one_mark_moment :
    ∀ᶠ R : Nat in atTop, ∀ b : Nat, ∀ marks : Fin b → (Fin 1 ↪ Fin (128 * R)),
      (Spin.piPMF (fun j => regionLaw (univ.map (marks j)))).expect
        (fun regions => inputMoment (Real.exp (-(2 / (128 * R)))) (regionStream regions).get ∅) ≤
          1000 * (1011 / 2000 : ℝ) ^ b := by
  have h := actual_fair_regionStream_margin (by norm_num : (0:ℝ) ≤ 2)
    (by norm_num : (0:ℝ) < 1/1000) (by norm_num : (0:ℝ) ≤ 101/100)
    (by norm_num : (0:ℝ) < 1/1000) continuum_one_weighted_bound
  filter_upwards [h] with R hR
  intro b marks
  convert hR b marks using 1 <;> norm_num <;> ring

end Spin.Structured.Placement
