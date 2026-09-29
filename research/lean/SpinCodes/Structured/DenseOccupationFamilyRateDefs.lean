import SpinCodes.Structured.DenseOccupationSelectedParameters
import SpinCodes.Structured.DenseOccupationOuterSupport
import SpinCodes.Structured.DenseNativeRemainder

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open ConcreteRoute ConcreteRoutedEncoder

/-- The common actual-experiment conclusion supplied by each numerical witness family. -/
def PointRate (α x m c η C : ℝ) : Prop :=
  ∀ {L b R d : ℕ}, 0 < L → 0 < b →
    ∀ (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b)),
    0 < totalWeight rows → totalWeight rows < L*b → density rows = α*x →
    128*R = L*b → (d:ℝ) ≤ (11/100)*((L:ℝ)*b) →
    ∀ K : ℝ, 0 ≤ K → K ≤ Real.exp (((L:ℝ)*b)*α*(m*x+c)) →
    K*(experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      (((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*C*Real.exp (-η*((L:ℝ)*b))

/-- A common rate on the full dense rectangle, with certified outer supports. -/
def DenseRates (η C : ℝ) : Prop :=
  ∀ α x : ℝ, α ∈ Set.Icc (1/10000) 1 → x ∈ Set.Icc (13/125) (112/125) →
    ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧ PointRate α x p.1 p.2 η C

theorem PointRate.mono_constant {α x m c η C C' : ℝ}
    (h : PointRate α x m c η C) (hC : C ≤ C') : PointRate α x m c η C' := by
  intro L b R d hL hb e rows hw0 hw1 hden hround hd K hK0 hK
  apply (h hL hb e rows hw0 hw1 hden hround hd K hK0 hK).trans
  exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hC (by positivity))
    (Real.exp_pos _).le

#print axioms PointRate.mono_constant
end Spin.Structured.DenseOccupationFixed
