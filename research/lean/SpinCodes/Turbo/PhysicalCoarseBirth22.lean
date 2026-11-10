import SpinCodes.Turbo.PhysicalGenericLive22
import SpinCodes.Turbo.PhysicalThinningRows22

/-!
# Coarse complete birth rows without syndrome-class census

The old `turboSelectedKernel22_zero_row_sound` splits births by expansion
weight. This alternative keeps zero mass and puts all nonzero birth mass
in the arbitrary-live coordinate. It changes only a proof witness, not the
encoder. Zero and total-mass upper bounds are explicit premises, and may
overlap because they bound different pieces of one actual measure.
-/
noncomputable section
namespace Spin.Turbo
open Spin.UniformMixer
attribute [local instance] Classical.propDecidable

def physicalCoarseBirthRow22 (Z T : ℝ) (i : Fin 11) : ℝ :=
  if i = 0 then Z else if i = 1 then T else 0

theorem physicalCoarseBirth22_parts (ν : State22 → ℝ) (Z T : ℝ)
    (hn : ∀ r, 0 ≤ ν r) (hZ : ν 0 ≤ Z) (hT : (∑ r, ν r) ≤ T) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, ν r = ∑ i, parts i r) ∧
      ∀ i, physicalPiece22 i (physicalCoarseBirthRow22 Z T i) (parts i) := by
  let p0 : State22 → ℝ := fun r => if r = 0 then ν r else 0
  let p1 : State22 → ℝ := fun r => if r = 0 then 0 else ν r
  have h0 : physicalPiece22 0 Z p0 := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro r; dsimp [p0]; split_ifs; exact hn r; exact le_rfl
    · intro r hr
      have hz : r = 0 := by
        by_contra h; exact hr (by simp only [p0, h, ite_false])
      simpa only [physicalClassAllows22, Fin.val_zero, ite_true] using hz
    · simpa only [p0, Finset.sum_ite_eq', Finset.mem_univ, ite_true] using hZ
    · intro hi; omega
  have h1 : physicalPiece22 1 T p1 := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro r; dsimp [p1]; split_ifs; exact le_rfl; exact hn r
    · intro r hr
      have hz : r ≠ 0 := by
        intro h; subst r; exact hr (by simp only [p1, ite_true])
      simpa only [physicalClassAllows22, Fin.val_one, Nat.one_ne_zero,
        ite_false, true_or, ite_true] using hz
    · apply (Finset.sum_le_sum (fun r _ => ?_)).trans hT
      dsimp [p1]
      split_ifs
      · exact hn r
      · exact le_rfl
    · intro hi; omega
  let parts : Fin 11 → State22 → ℝ := fun i r =>
    if i = 0 then p0 r else if i = 1 then p1 r else 0
  refine ⟨parts, ?_, ?_⟩
  · intro r
    by_cases hr : r = 0
    · subst r
      simp only [parts, p0, p1, ite_true, ite_self,
        Finset.sum_ite_eq', Finset.mem_univ]
    · simp only [parts, p0, p1, hr, ite_false]
      symm
      apply Finset.sum_eq_single 1
      · intro i _ hi
        by_cases hi0 : i = 0 <;> simp only [hi0, hi, ite_true, ite_false]
      · intro h; exact (h (Finset.mem_univ 1)).elim
  · intro i
    by_cases hi0 : i = 0
    · subst i
      simpa only [parts, physicalCoarseBirthRow22, ite_true] using h0
    · by_cases hi1 : i = 1
      · subst i
        simpa only [parts, physicalCoarseBirthRow22, hi0, ite_false, ite_true] using h1
      · simpa only [parts, physicalCoarseBirthRow22, hi0, hi1, ite_false] using
          physicalPiece22_zero i 0 (le_refl 0)

theorem physicalInputKernel22_coarse_birth_row (P : FinPMF PacketStep16)
    (z : ℝ) (hz : 0 ≤ z) (Z T : ℝ)
    (hZ : physicalInputKernel22 P z 0 0 ≤ Z)
    (hT : physicalInputMoment22 P z 0 ≤ T)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 0 b μ) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (physicalInputKernel22 P z) μ r = ∑ i, parts i r) ∧
      ∀ i, physicalPiece22 i (b * physicalCoarseBirthRow22 Z T i) (parts i) := by
  obtain ⟨parts, hsum, hp⟩ := physicalCoarseBirth22_parts
    (physicalInputKernel22 P z 0) Z T
    (physicalInputKernel22_nonneg P z hz 0) hZ
    (by rw [physicalInputKernel22_row_sum]; exact hT)
  refine ⟨(fun i r => μ 0 * parts i r), ?_, ?_⟩
  · intro r
    unfold physicalKernelAction22
    rw [Finset.sum_eq_single 0]
    · rw [hsum, Finset.mul_sum]
    · intro q _ hq
      rw [physicalZeroPiece22_supported hμ q hq, zero_mul]
    · intro h; exact (h (Finset.mem_univ 0)).elim
  · intro i
    have hc := physicalPiece22_scale (hp i) (hμ.1 0)
    apply physicalPiece22_mono hc
    have hi : 0 ≤ physicalCoarseBirthRow22 Z T i :=
      (Finset.sum_nonneg (fun r _ => (hp i).1 r)).trans (hp i).2.2.1
    exact mul_le_mul_of_nonneg_right (physicalZeroPiece22_mass hμ) hi

/-- Applies directly to each exactly-j-nonzero kernel, before thinning. -/
theorem turboPhysicalKernel22_coarse_birth_row (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (Z T : ℝ)
    (hZ : turboPhysicalKernel22 j z 0 0 ≤ Z)
    (hT : physicalInputMoment22 (packetStepLaw22 j) z 0 ≤ T)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 0 b μ) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ i, parts i r) ∧
      ∀ i, physicalPiece22 i (b * physicalCoarseBirthRow22 Z T i) (parts i) := by
  have hK : Matrix.of (turboPhysicalKernel22 j z) =
      physicalInputKernel22 (packetStepLaw22 j) z := rfl
  have hZ' : physicalInputKernel22 (packetStepLaw22 j) z 0 0 ≤ Z := by
    rw [← hK]
    exact hZ
  rw [hK]
  exact physicalInputKernel22_coarse_birth_row
    (packetStepLaw22 j) z hz Z T hZ' hT b μ hμ

#print axioms physicalCoarseBirth22_parts
#print axioms turboPhysicalKernel22_coarse_birth_row
end Spin.Turbo
