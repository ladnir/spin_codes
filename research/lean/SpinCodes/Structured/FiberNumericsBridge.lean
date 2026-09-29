import SpinCodes.Structured.FiberNumerics

/-! Apply checked integer rows to the concrete feedback map. The sole
spectrum premise is kept explicit until the exhaustive map replay supplies
it. The numerical kernel and cap claims are conclusions, not assumptions. -/

noncomputable section
namespace Spin.Structured.FiberNumerics
open Finset

theorem kernel_count_of_certificate (spectrum : List ℕ) {j k : ℕ} {values : List ℤ}
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = spectrum.getD w 0)
    (hv : (List.range 129).map (kraw j) = values)
    (hk : signedSum spectrum values = 524288 * (k : ℤ)) :
    (ConcreteMaps.kernelLayer j).card = k := by
  have h := ConcreteMaps.kernel_fourier_weights j
  simp only [hs] at h
  rw [← signedSum_eq spectrum hv, hk] at h
  have he : ((ConcreteMaps.kernelLayer j).card : ℤ) = k := by omega
  exact_mod_cast he

theorem pair_count_of_certificate (spectrum : List ℕ) {j V : ℕ} {values : List ℤ}
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = spectrum.getD w 0)
    (hv : (List.range 129).map (kraw j) = values)
    (hV : squareSum spectrum values = 524288 * (V : ℤ)) :
    (ConcreteMaps.equalSyndromePairs j).card = V := by
  have h := ConcreteMaps.feedback_parseval_weights j
  simp only [hs] at h
  rw [← squareSum_eq spectrum hv, hV] at h
  have he : ((ConcreteMaps.equalSyndromePairs j).card : ℤ) = V := by omega
  exact_mod_cast he

theorem fiber_cap_of_fourier_certificate (spectrum : List ℕ) {j cap : ℕ} {values : List ℤ}
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = spectrum.getD w 0)
    (hv : (List.range 129).map (kraw j) = values)
    (hc : absSum spectrum values < 524288 * ((cap : ℤ) + 1))
    (q : Finset (Fin 19)) :
    (ConcreteMaps.syndromeFiber j q).card ≤ cap := by
  apply ConcreteMaps.fiber_cap_fourier j cap q
  simp only [hs]
  rw [← absSum_eq spectrum hv]
  exact hc

end Spin.Structured.FiberNumerics
