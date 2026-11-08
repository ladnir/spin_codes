import SpinCodes.Turbo.PhysicalRowDominance22

/-!
# Complete-row thinning

The retained numerical producer first bounds each exact nonzero-packet
occupancy, then applies binomial thinning. This file transfers soundness
of complete rows through that nonnegative mixture. The corresponding
full-state probability identity is checked in `PhysicalThinning22`.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

theorem physicalPiece22_scale {i : Fin 11} {b w : ℝ}
    {μ : State22 → ℝ} (h : physicalPiece22 i b μ)
    (hw : 0 ≤ w) :
    physicalPiece22 i (w * b) (fun q => w * μ q) := by
  constructor
  · intro q
    exact mul_nonneg hw (h.1 q)
  constructor
  · intro q hq
    have hμ : μ q ≠ 0 := by
      intro hz
      exact hq (by simp [hz])
    exact h.2.1 q hμ
  constructor
  · calc
      (∑ q : State22, w * μ q) =
          w * ∑ q : State22, μ q := by rw [Finset.mul_sum]
      _ ≤ w * b := mul_le_mul_of_nonneg_left h.2.2.1 hw
  · intro hi q
    calc
      w * μ q ≤ w * (b / ((2 : ℝ) ^ 22 - 1)) :=
        mul_le_mul_of_nonneg_left (h.2.2.2 hi q) hw
      _ = (w * b) / ((2 : ℝ) ^ 22 - 1) := by ring

def physicalWeightedKernel22 {ι : Type*} [Fintype ι]
    (w : ι → ℝ)
    (K : ι → Matrix State22 State22 ℝ) :
    Matrix State22 State22 ℝ :=
  Matrix.of fun q r => ∑ x : ι, w x * K x q r

def physicalWeightedRow22 {ι : Type*} [Fintype ι]
    (w : ι → ℝ)
    (B : ι → Matrix (Fin 11) (Fin 11) ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  Matrix.of fun i k => ∑ x : ι, w x * B x i k

theorem physicalRowSound22_weighted
    {ι : Type*} [Fintype ι]
    (w : ι → ℝ)
    (K : ι → Matrix State22 State22 ℝ)
    (B : ι → Matrix (Fin 11) (Fin 11) ℝ)
    (hw : ∀ x, 0 ≤ w x)
    (hB : ∀ x, physicalRowSound22 (K x) (B x)) :
    physicalRowSound22
      (physicalWeightedKernel22 w K)
      (physicalWeightedRow22 w B) := by
  classical
  intro i b μ hμ
  choose out hout hshape using
    (fun x : ι => hB x i b μ hμ)
  let parts : Fin 11 → State22 → ℝ :=
    fun k r => ∑ x : ι, w x * out x k r
  refine ⟨parts, ?_, ?_⟩
  · intro r
    unfold physicalKernelAction22 physicalWeightedKernel22
    change (∑ q : State22,
      μ q * ∑ x : ι, w x * K x q r) = _
    simp only [Finset.mul_sum]
    calc
      (∑ q : State22, ∑ x : ι,
          μ q * (w x * K x q r)) =
        ∑ x : ι, w x * physicalKernelAction22 (K x) μ r := by
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro x _
        unfold physicalKernelAction22
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro q _
        ring
      _ = ∑ x : ι, w x * ∑ k : Fin 11, out x k r := by
        simp only [hout]
      _ = ∑ k : Fin 11, parts k r := by
        simp only [parts, Finset.mul_sum]
        rw [Finset.sum_comm]
  · intro k
    have hpart : physicalPiece22 k
        (∑ x : ι, w x * (b * B x i k))
        (parts k) := by
      exact physicalPiece22_sum k
        (fun x => w x * (b * B x i k))
        (fun x r => w x * out x k r)
        (fun x => physicalPiece22_scale (hshape x k) (hw x))
    have hbound :
        (∑ x : ι, w x * (b * B x i k)) =
          b * physicalWeightedRow22 w B i k := by
      unfold physicalWeightedRow22
      change (∑ x : ι, w x * (b * B x i k)) =
        b * ∑ x : ι, w x * B x i k
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro x _
      ring
    rw [← hbound]
    exact hpart

theorem thinningWeight22_nonneg (j k : Fin 17) :
    0 ≤ thinningWeight22 j k := by
  unfold thinningWeight22
  positivity

def turboThinnedRow22 (j : Fin 17) (_z : ℝ)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  physicalWeightedRow22 (thinningWeight22 j) B

theorem turboThinnedRow22_sound (j : Fin 17)
    (z : ℝ)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ k, physicalRowSound22
      (Matrix.of (turboPhysicalKernel22 k z)) (B k)) :
    physicalRowSound22
      (Matrix.of (turboSelectedKernel22 j z))
      (turboThinnedRow22 j z B) := by
  have h := physicalRowSound22_weighted
    (thinningWeight22 j)
    (fun k => Matrix.of (turboPhysicalKernel22 k z))
    B (thinningWeight22_nonneg j) hB
  have hK :
      physicalWeightedKernel22 (thinningWeight22 j)
        (fun k => Matrix.of (turboPhysicalKernel22 k z)) =
      Matrix.of (turboSelectedKernel22 j z) := by
    ext q r
    change (∑ k : Fin 17,
      thinningWeight22 j k * turboPhysicalKernel22 k z q r) =
      turboSelectedKernel22 j z q r
    exact (turboSelectedKernel22_mixture j z q r).symm
  rw [← hK]
  exact h

#print axioms physicalPiece22_scale
#print axioms physicalRowSound22_weighted
#print axioms turboThinnedRow22_sound

end Spin.Turbo
