import SpinCodes.Structured.ConcreteOuterMajorantSpectrum
import SpinCodes.Structured.ConcreteOuterCounting

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter

/-- The first certified affine support bounds the majorant globally. -/
theorem refined_le_first (x : ℝ) :
    Spin.Majorant.refined.toFun x ≤ (833/500)*x-43311/250000 := by
  have h := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member0 x
  norm_num at h ⊢
  linarith

/-- The finite shared-seed selection cost per occupied row. -/
def outerCost (b : ℕ) : ℝ := (b:ℝ)^2 * Real.exp (spectrumLogError b)

theorem outerCost_pos {b : ℕ} (hb : 0 < b) : 0 < outerCost b := by
  unfold outerCost
  positivity

theorem selected_spectrum_affine {k w : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (hw : w ∈ W) (hr : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125)) :
    (spectrum seed w:ℝ) ≤ outerCost (k*24) *
      Real.exp (((k*24:ℕ):ℝ)*((833/500)*((w:ℝ)/(k*24:ℕ))-43311/250000)) := by
  have he := expected_spectrum_refined hk hr
  have ha := refined_le_first ((w:ℝ)/(k*24:ℕ))
  have hb : 0 ≤ ((k*24:ℕ):ℝ) := by positivity
  have hh := mul_le_mul_of_nonneg_left ha hb
  have hs := (hg.2 w hw).trans (mul_le_mul_of_nonneg_left he (sq_nonneg _))
  apply hs.trans
  unfold outerCost
  rw [mul_assoc, ← Real.exp_add]
  apply mul_le_mul_of_nonneg_left _ (sq_nonneg _)
  apply Real.exp_le_exp.mpr
  linarith

/-- Exact affine averaging for a fixed row-weight profile of the selected shared seed. -/
theorem selected_profile_affine {ι : Type*} {k : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (S : Finset ι) (weights : ι → ℕ) {x : ℝ}
    (hw : ∀ i ∈ S, weights i ∈ W)
    (hr : ∀ i ∈ S, (weights i:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125))
    (hmean : ∑ i ∈ S, (weights i:ℝ)/(k*24:ℕ) = (S.card:ℝ)*x) :
    (∏ i ∈ S, (spectrum seed (weights i):ℝ)) ≤ outerCost (k*24)^S.card *
      Real.exp (((k*24:ℕ):ℝ)*(S.card:ℝ)*((833/500)*x-43311/250000)) := by
  calc
    (∏ i ∈ S, (spectrum seed (weights i):ℝ)) ≤
        ∏ i ∈ S, (outerCost (k*24) *
          Real.exp (((k*24:ℕ):ℝ)*((833/500)*((weights i:ℝ)/(k*24:ℕ))-43311/250000))) :=
      Finset.prod_le_prod₀ (fun _ _ => Nat.cast_nonneg _) (fun i hi =>
        selected_spectrum_affine hk seed W hg (hw i hi) (hr i hi))
    _ = _ := by
      rw [Finset.prod_mul_distrib, Finset.prod_const, ← Real.exp_sum]
      congr 2
      rw [← Finset.mul_sum, Finset.sum_sub_distrib, ← Finset.mul_sum, hmean]
      simp only [Finset.sum_const, nsmul_eq_mul]
      ring

#print axioms refined_le_first
#print axioms selected_spectrum_affine
#print axioms selected_profile_affine
end Spin.Structured.DenseOccupationFixed
