import SpinCodes.Structured.ConcreteFourierTilt
import SpinCodes.Structured.ConcreteFourierOverlap

noncomputable section
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteMaps
open scoped symmDiff

theorem sum_by_weight_real {s t : ℕ} (F : Finset (Fin s) → Finset (Fin t)) (f : ℕ → ℝ) :
    ∑ q : Finset (Fin s), f (F q).card =
      ∑ w : Fin (t+1), (weightCounts F w : ℝ)*f w := by
  let g : Finset (Fin s) → Fin (t+1) := fun q =>
    ⟨(F q).card, Nat.lt_succ_of_le (by simpa using Finset.card_le_univ (F q))⟩
  have h := Finset.sum_fiberwise' (univ : Finset (Finset (Fin s))) g (fun w => f w.val)
  change (∑ w : Fin (t+1), ∑ q ∈ univ.filter (fun q => g q = w), f w.val) =
    ∑ q : Finset (Fin s), f (F q).card at h
  rw [← h]
  apply sum_congr rfl
  intro w _
  have he : univ.filter (fun q => g q = w) = univ.filter (fun q => (F q).card = w.val) := by
    ext q
    simp only [mem_filter, mem_univ, true_and, Fin.ext_iff, g]
  rw [he, sum_const]
  simp only [nsmul_eq_mul, weightCounts]

def spectrumCap (d : ℕ) (a b : ℝ) : ℝ :=
  (∑ w : Fin 129, (FiberNumerics.Data.spectrum.getD w 0 : ℝ)*overlapCap 128 d w a b) / 524288

/-- The Fourier cap is computed from the independently checked actual transpose spectrum. -/
theorem fourierCap_le_spectrum {p : Fin 128 → ℝ} (Y : Finset (Fin 128)) {a b : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hp : ∀ i, |1-2*p i| = if i ∈ Y then b else a) :
    fourierCap p ≤ spectrumCap Y.card a b := by
  unfold fourierCap spectrumCap
  apply div_le_div_of_nonneg_right _ (by norm_num)
  calc
    _ ≤ ∑ q : Finset (Fin 19), overlapCap 128 Y.card (CtransposeSet q).card a b := by
      apply sum_le_sum
      intro q _
      simp_rw [hp]
      exact overlap_product_le Y (CtransposeSet q) ha hb
    _ = _ := by
      rw [sum_by_weight_real CtransposeSet (fun w => overlapCap 128 Y.card w a b)]
      simp only [transpose_spectrum]

def ratio0 (β z : ℝ) : ℝ := |1-2*(β*z/(1-β+β*z))|
def ratio1 (β z : ℝ) : ℝ := |1-2*(β/(β+(1-β)*z))|

theorem emissionTilt_ratios (β z : ℝ) (Y : Finset (Fin 128)) (i : Fin 128) :
    |1-2*emissionTilt β z Y i| = if i ∈ Y then ratio1 β z else ratio0 β z := by
  by_cases hi : i ∈ Y
  · simp only [hi, ite_true, emissionTilt, tiltProb, tiltMass, emit0, emit1, mul_one, ratio1]
    rw [add_comm ((1-β)*z) β]
  · simp only [hi, ite_false, emissionTilt, tiltProb, tiltMass, emit0, emit1, mul_one, ratio0]

theorem emission_fourierCap_le (β z : ℝ) (Y : Finset (Fin 128)) :
    fourierCap (emissionTilt β z Y) ≤ spectrumCap Y.card (ratio0 β z) (ratio1 β z) :=
  fourierCap_le_spectrum Y (abs_nonneg _) (abs_nonneg _) (emissionTilt_ratios β z Y)

/-- The entering word influences this syndrome bound only through its weight. -/
theorem weighted_syndrome_spectrum_le {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z)
    (Y : Finset (Fin 128)) (s : Finset (Fin 19)) :
    (poissonBinom (fun _ : Fin 128 => β) (fun _ => hb0.le) (fun _ => hb1.le)).expect
      (fun X => z^(X ∆ Y).card * if Cset X = s then 1 else 0) ≤
      ConcreteScalar.scalarF β z Y.card * spectrumCap Y.card (ratio0 β z) (ratio1 β z) := by
  refine (weighted_syndrome_le hb0 hb1 hz Y s).trans ?_
  apply mul_le_mul_of_nonneg_left (emission_fourierCap_le β z Y)
  have : 0 ≤ 1-β := by linarith
  unfold ConcreteScalar.scalarF
  positivity

end Spin.Structured.ConcreteFourier
