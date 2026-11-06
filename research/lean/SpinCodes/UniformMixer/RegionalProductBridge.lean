import SpinCodes.UniformMixer.RegionalProductError

/-!
# Identifying the coarse regional product with rounded support gaps

The rounded product is exactly the coarse two-state product indexed by
the marked support. This identity precedes the analytic comparison with
the continuum product; it does not introduce an approximation.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.FiniteKernel (RowError)

theorem coarseRegionProduct_supportMarks (γ : ℝ) (R : ℕ)
    (xs : List (Fin (64 * R))) (prev : ℕ) :
    coarseRegionProduct
        (fun g => Real.exp (-γ * ((g : ℝ) / R)))
        (supportMarks R xs prev) (supportLastGap R xs prev) =
      roundedSupportProduct γ R xs prev := by
  induction xs generalizing prev with
  | nil => rfl
  | cons x tail ih =>
      simp only [supportMarks, supportLastGap, coarseRegionProduct,
        roundedSupportProduct, timeRegionGap]
      rw [ih]

theorem production_coarse_gap_factor (θ : ℝ) (R g : ℕ) :
    Real.exp (-((θ / (64 * (R : ℝ))) * liveWeightAverage coordMap)) ^ g =
      Real.exp (-(regionGamma θ) * ((g : ℝ) / R)) := by
  rw [← Real.exp_nat_mul]
  congr 1
  unfold regionGamma
  ring

theorem production_coarseRegionProduct_supportMarks (θ : ℝ) (R : ℕ)
    (xs : List (Fin (64 * R))) (prev : ℕ) :
    coarseRegionProduct
        (fun g => Real.exp
          (-((θ / (64 * (R : ℝ))) * liveWeightAverage coordMap)) ^ g)
        (supportMarks R xs prev) (supportLastGap R xs prev) =
      roundedSupportProduct (regionGamma θ) R xs prev := by
  have hfun :
      (fun g : ℕ => Real.exp
        (-((θ / (64 * (R : ℝ))) * liveWeightAverage coordMap)) ^ g) =
      (fun g : ℕ => Real.exp (-(regionGamma θ) * ((g : ℝ) / R))) := by
    funext g
    exact production_coarse_gap_factor θ R g
  rw [hfun]
  exact coarseRegionProduct_supportMarks _ _ _ _

theorem good_support_coarse_ideal_rowError {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    {γ : ℝ} (hγ : 0 ≤ γ) :
    RowError
      (coarseRegionProduct
        (fun g => Real.exp (-γ * ((g : ℝ) / R)))
        (supportMarks R (T.sort (· ≤ ·)) 0)
        (supportLastGap R (T.sort (· ≤ ·)) 0))
      (idealSupportProduct γ R (T.sort (· ≤ ·)) 0)
      (γ * (3 * (T.card + 1 : ℕ) / (R : ℝ))) := by
  rw [coarseRegionProduct_supportMarks]
  have h := supportProducts_rowError
    (goodRegionalSites_validSupport hR T hgood) hR hγ 0
    (Or.inl ⟨rfl, rfl⟩) (by omega : 0 < 64 * R)
  intro b
  exact (h b).trans
    (mul_le_mul_of_nonneg_left (good_support_gap_error hR T hgood) hγ)

end Spin.UniformMixer
