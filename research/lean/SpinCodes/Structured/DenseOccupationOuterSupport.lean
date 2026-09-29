import SpinCodes.Structured.DenseOccupationOuterAffine

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter

/-- Every certified support gives an actual selected-seed spectral bound. -/
theorem selected_spectrum_support {k w : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (hw : w ∈ W) (hr : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125))
    {p : ℚ × ℚ} (hp : p ∈ Spin.Majorant.refined.supports) :
    (spectrum seed w:ℝ) ≤ outerCost (k*24) *
      Real.exp (((k*24:ℕ):ℝ)*((p.1:ℝ)*((w:ℝ)/(k*24:ℕ))+(p.2:ℝ))) := by
  have he := expected_spectrum_refined hk hr
  have ha := Spin.Majorant.refined.toFun_le_support hp ((w:ℝ)/(k*24:ℕ))
  have hh := mul_le_mul_of_nonneg_left ha (show 0 ≤ ((k*24:ℕ):ℝ) by positivity)
  have hs := (hg.2 w hw).trans (mul_le_mul_of_nonneg_left he (sq_nonneg _))
  apply hs.trans
  unfold outerCost
  rw [mul_assoc, ← Real.exp_add]
  apply mul_le_mul_of_nonneg_left _ (sq_nonneg _)
  apply Real.exp_le_exp.mpr
  linarith

/-- Any certified rational affine support controls the exact product of realized spectra. -/
theorem selected_profile_support {ι : Type*} {k : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (S : Finset ι) (weights : ι → ℕ) {x : ℝ}
    (hw : ∀ i ∈ S, weights i ∈ W)
    (hr : ∀ i ∈ S, (weights i:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125))
    (hmean : ∑ i ∈ S, (weights i:ℝ)/(k*24:ℕ) = (S.card:ℝ)*x)
    {p : ℚ × ℚ} (hp : p ∈ Spin.Majorant.refined.supports) :
    (∏ i ∈ S, (spectrum seed (weights i):ℝ)) ≤ outerCost (k*24)^S.card *
      Real.exp (((k*24:ℕ):ℝ)*(S.card:ℝ)*((p.1:ℝ)*x+(p.2:ℝ))) := by
  calc
    (∏ i ∈ S, (spectrum seed (weights i):ℝ)) ≤
        ∏ i ∈ S, (outerCost (k*24) *
          Real.exp (((k*24:ℕ):ℝ)*((p.1:ℝ)*((weights i:ℝ)/(k*24:ℕ))+(p.2:ℝ)))) :=
      Finset.prod_le_prod₀ (fun _ _ => Nat.cast_nonneg _) (fun i hi =>
        selected_spectrum_support hk seed W hg (hw i hi) (hr i hi) hp)
    _ = _ := by
      rw [Finset.prod_mul_distrib, Finset.prod_const, ← Real.exp_sum]
      congr 2
      rw [← Finset.mul_sum, Finset.sum_add_distrib, ← Finset.mul_sum, hmean]
      simp only [Finset.sum_const, nsmul_eq_mul]
      ring

#print axioms selected_spectrum_support
#print axioms selected_profile_support
end Spin.Structured.DenseOccupationFixed
