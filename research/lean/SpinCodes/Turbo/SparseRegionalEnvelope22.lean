import SpinCodes.Turbo.SparsePhysicalKernel22
import SpinCodes.Turbo.PhysicalProduct22

/-! Reuse `physical_list_product_mass_le` once for the whole ordered
physical stream. Average over the joint region-count profiles only after
that comparison; this avoids any intermediate state sum or reset. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem sparsePhysicalNested_mass_le (R : ℕ) (z : ℝ)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hs : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j))
    (counts : Fin 64 → Fin R → Fin 17) :
    (∑ q : State22, (List.ofFn (fun region : Fin 64 =>
      (List.ofFn (fun step : Fin R => Matrix.of (turboPhysicalKernel22 (counts region step) z))).prod)).prod 0 q) ≤
    ∑ k : Fin 11, (List.ofFn (fun region : Fin 64 =>
      (List.ofFn (fun step : Fin R => B (counts region step))).prod)).prod 0 k := by
  let steps : List PhysicalPair22 := (List.ofFn (fun region : Fin 64 =>
    List.ofFn (fun step : Fin R =>
      (Matrix.of (turboPhysicalKernel22 (counts region step) z), B (counts region step))))).flatten
  have hsteps : ∀ p ∈ steps, physicalRowSound22 p.1 p.2 := by
    intro p hp
    obtain ⟨row, hrow, hp⟩ := List.mem_flatten.mp hp
    obtain ⟨region, hregion⟩ := List.mem_ofFn.mp hrow
    rw [← hregion] at hp
    obtain ⟨step, hstep⟩ := List.mem_ofFn.mp hp
    rw [← hstep]
    exact hs _
  have h := physical_list_product_mass_le steps hsteps
  simpa only [steps, List.map_flatten, List.map_ofFn, List.map_map,
    Function.comp_def, List.prod_flatten] using h

theorem sparseExactRegional_product_mass_le (R : ℕ) (z : ℝ)
    (T : Fin 64 → Finset (Fin (R * 16)))
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hs : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j)) :
    (∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R (T region) z)).prod 0 q) ≤
      ∑ k : Fin 11, (List.ofFn (fun region => sparsePhysicalRegionMatrix22 B R (T region))).prod 0 k := by
  let P := Spin.piPMF (fun region : Fin 64 => regionalStepProfileLaw22 R (T region))
  have hfull : P.expect (fun counts => ∑ q : State22,
      (List.ofFn (fun region : Fin 64 => (List.ofFn (fun step : Fin R =>
        Matrix.of (turboPhysicalKernel22 (counts region step) z))).prod)).prod 0 q) =
      ∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R (T region) z)).prod 0 q := by
    rw [FinPMF.expect_sum]
    change (∑ q : State22, SparseMatrixAverage22.average P
      (fun counts => (List.ofFn (fun region : Fin 64 => (List.ofFn (fun step : Fin R =>
        Matrix.of (turboPhysicalKernel22 (counts region step) z))).prod)).prod) 0 q) = _
    unfold P
    rw [SparseMatrixAverage22.average_pi_product
      (fun region : Fin 64 => regionalStepProfileLaw22 R (T region))
      (fun (_ : Fin 64) (counts : Fin R → Fin 17) =>
        (List.ofFn (fun step : Fin R => Matrix.of (turboPhysicalKernel22 (counts step) z))).prod)]
    simp only [sparseExactRegionKernel22_eq_physical_profile, sparsePhysicalRegionMatrix22]
  have henv : P.expect (fun counts => ∑ k : Fin 11,
      (List.ofFn (fun region : Fin 64 => (List.ofFn (fun step : Fin R => B (counts region step))).prod)).prod 0 k) =
      ∑ k : Fin 11, (List.ofFn (fun region => sparsePhysicalRegionMatrix22 B R (T region))).prod 0 k := by
    rw [FinPMF.expect_sum]
    change (∑ k : Fin 11, SparseMatrixAverage22.average P
      (fun counts => (List.ofFn (fun region : Fin 64 =>
        (List.ofFn (fun step : Fin R => B (counts region step))).prod)).prod) 0 k) = _
    unfold P
    rw [SparseMatrixAverage22.average_pi_product
      (fun region : Fin 64 => regionalStepProfileLaw22 R (T region))
      (fun (_ : Fin 64) (counts : Fin R → Fin 17) =>
        (List.ofFn (fun step : Fin R => B (counts step))).prod)]
    rfl
  rw [← hfull, ← henv]
  exact P.expect_mono (fun counts => sparsePhysicalNested_mass_le R z B hs counts)

#print axioms sparseExactRegional_product_mass_le
end Spin.Turbo


