import SpinCodes.Turbo.PhysicalBirth22
import SpinCodes.UniformMixer.ScalableEnvelope

/-!
# Live-state rows of the Turbo physical kernel

For nonzero entering state, the independent uniform invertible update
maps that state uniformly over all nonzero values. The fixed feedback
then shifts the uniform law. The weighted output and feedback are
correlated, so the return-to-zero moment retains the feedback event
inside the packet-input expectation.

This follows the checked `UniformMixer/ScalableWeightedRows` argument,
with the actual zero-inclusive four-bit `selectedStepLaw22` replacing
the Bernoulli-bit input law.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer
attribute [local instance] Classical.propDecidable

def turboLiveMoment22 (j : Fin 17) (z : ℝ)
    (q : State22) : ℝ :=
  (selectedStepLaw22 j).expect (fun x =>
    z ^ stepWeight (packetStepInput22 x) q)

def turboLiveReturnMoment22 (j : Fin 17) (z : ℝ)
    (q : State22) : ℝ :=
  (selectedStepLaw22 j).expect (fun x =>
    if stepFeedback (packetStepInput22 x) = 0 then 0 else
      z ^ stepWeight (packetStepInput22 x) q)

theorem turboLiveMoment22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (q : State22) :
    0 ≤ turboLiveMoment22 j z q := by
  unfold turboLiveMoment22
  apply FinPMF.expect_nonneg
  intro x
  exact pow_nonneg hz _

theorem turboLiveReturnMoment22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (q : State22) :
    0 ≤ turboLiveReturnMoment22 j z q := by
  unfold turboLiveReturnMoment22
  apply FinPMF.expect_nonneg
  intro x
  by_cases h : stepFeedback (packetStepInput22 x) = 0
  · simp only [h, ite_true]
    exact le_refl 0
  · simp only [h, ite_false]
    exact pow_nonneg hz _

theorem turboSelectedKernel22_live_entry (j : Fin 17)
    (z : ℝ) (q r : State22) (hq : q ≠ 0) :
    turboSelectedKernel22 j z q r =
      (selectedStepLaw22 j).expect (fun x =>
        if r = stepFeedback (packetStepInput22 x) then 0 else
          z ^ stepWeight (packetStepInput22 x) q /
            (liveCount State22 : ℝ)) := by
  unfold turboSelectedKernel22
  apply congrArg
  funext x
  change z ^ stepWeight (packetStepInput22 x) q *
      update (stepFeedback (packetStepInput22 x)) q r = _
  rw [update_live _ q r hq]
  split_ifs <;> ring

theorem turboSelectedKernel22_live_zero (j : Fin 17)
    (z : ℝ) (q : State22) (hq : q ≠ 0) :
    turboSelectedKernel22 j z q 0 =
      turboLiveReturnMoment22 j z q /
        (liveCount State22 : ℝ) := by
  rw [turboSelectedKernel22_live_entry j z q 0 hq]
  unfold turboLiveReturnMoment22 FinPMF.expect
  simp only [div_eq_mul_inv, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro x _
  by_cases hx : stepFeedback (packetStepInput22 x) = 0
  · simp only [hx, ite_true]
    ring
  · simp only [hx, eq_comm, ite_false]
    ring

theorem turboSelectedKernel22_live_atom_le (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (q r : State22) (hq : q ≠ 0) :
    turboSelectedKernel22 j z q r ≤
      turboLiveMoment22 j z q /
        (liveCount State22 : ℝ) := by
  rw [turboSelectedKernel22_live_entry j z q r hq]
  calc
    _ ≤ (selectedStepLaw22 j).expect (fun x =>
        z ^ stepWeight (packetStepInput22 x) q /
          (liveCount State22 : ℝ)) := by
      apply FinPMF.expect_mono
      intro x
      split_ifs
      · exact div_nonneg (pow_nonneg hz _) (Nat.cast_nonneg _)
      · exact le_refl _
    _ = _ := by
      unfold turboLiveMoment22 FinPMF.expect
      simp only [div_eq_mul_inv, Finset.sum_mul, mul_assoc]

theorem turboSelectedKernel22_row_sum (j : Fin 17)
    (z : ℝ) (q : State22) :
    (∑ r : State22, turboSelectedKernel22 j z q r) =
      turboLiveMoment22 j z q := by
  calc
    _ = (selectedStepLaw22 j).expect (fun x =>
          ∑ r : State22,
            stepKernel z (packetStepInput22 x) q r) := by
      rw [FinPMF.expect_sum]
      rfl
    _ = _ := by
      unfold turboLiveMoment22
      apply congrArg
      funext x
      simp only [stepKernel, ← Finset.mul_sum, update_sum, mul_one]

theorem physicalLivePiece22_mass {i : Fin 11} {b : ℝ}
    {μ : State22 → ℝ} (hi : i ≠ 0)
    (hμ : physicalPiece22 i b μ) :
    (∑ q ∈ liveStateSet, μ q) ≤ b := by
  have hsplit := Finset.add_sum_erase
    (Finset.univ : Finset State22) μ
    (Finset.mem_univ (0 : State22))
  change μ 0 + (∑ q ∈ liveStateSet, μ q) =
    ∑ q : State22, μ q at hsplit
  rw [physicalPiece22_live_zero hi hμ, zero_add] at hsplit
  rw [hsplit]
  exact hμ.2.2.1

theorem physicalLivePiece22_max_bound {i : Fin 11}
    {b : ℝ} {μ : State22 → ℝ}
    (hi : i ≠ 0) (hμ : physicalPiece22 i b μ)
    (f : State22 → ℝ) (M : ℝ) (hM : 0 ≤ M)
    (hf : ∀ q, μ q ≠ 0 → f q ≤ M) :
    (∑ q : State22, μ q * f q) ≤ b * M := by
  have hzero : μ 0 = 0 := physicalPiece22_live_zero hi hμ
  have hsplit := Finset.add_sum_erase
    (Finset.univ : Finset State22)
    (fun q => μ q * f q) (Finset.mem_univ (0 : State22))
  change μ 0 * f 0 + (∑ q ∈ liveStateSet, μ q * f q) =
    ∑ q : State22, μ q * f q at hsplit
  rw [hzero, zero_mul, zero_add] at hsplit
  rw [← hsplit]
  calc
    (∑ q ∈ liveStateSet, μ q * f q) ≤
        ∑ q ∈ liveStateSet, μ q * M := by
      apply Finset.sum_le_sum
      intro q _
      by_cases hq : μ q = 0
      · simp [hq]
      · exact mul_le_mul_of_nonneg_left (hf q hq) (hμ.1 q)
    _ = (∑ q ∈ liveStateSet, μ q) * M := by
      rw [Finset.sum_mul]
    _ ≤ b * M := mul_le_mul_of_nonneg_right
      (physicalLivePiece22_mass hi hμ) hM

theorem physicalLiveDenom22 :
    (liveCount State22 : ℝ) = (2 : ℝ) ^ 22 - 1 := by
  rw [state22_liveCount]
  norm_num

theorem physicalUniformPiece22_mean_bound
    {b : ℝ} {μ : State22 → ℝ}
    (hμ : physicalPiece22 2 b μ)
    (f : State22 → ℝ) (hf : ∀ q, 0 ≤ f q)
    (M : ℝ)
    (hmean : (∑ q ∈ liveStateSet, f q) /
      (liveCount State22 : ℝ) ≤ M) :
    (∑ q : State22, μ q * f q) ≤ b * M := by
  have hb : 0 ≤ b :=
    (Finset.sum_nonneg (fun q _ => hμ.1 q)).trans hμ.2.2.1
  have hzero : μ 0 = 0 :=
    physicalPiece22_live_zero (by decide) hμ
  have hsplit := Finset.add_sum_erase
    (Finset.univ : Finset State22)
    (fun q => μ q * f q) (Finset.mem_univ (0 : State22))
  change μ 0 * f 0 + (∑ q ∈ liveStateSet, μ q * f q) =
    ∑ q : State22, μ q * f q at hsplit
  rw [hzero, zero_mul, zero_add] at hsplit
  rw [← hsplit]
  have hcap (q : State22) :
      μ q ≤ b / (liveCount State22 : ℝ) := by
    rw [physicalLiveDenom22]
    exact hμ.2.2.2 rfl q
  calc
    (∑ q ∈ liveStateSet, μ q * f q) ≤
        ∑ q ∈ liveStateSet,
          (b / (liveCount State22 : ℝ)) * f q := by
      apply Finset.sum_le_sum
      intro q _
      exact mul_le_mul_of_nonneg_right (hcap q) (hf q)
    _ = b * ((∑ q ∈ liveStateSet, f q) /
          (liveCount State22 : ℝ)) := by
      rw [← Finset.mul_sum]
      ring
    _ ≤ b * M := mul_le_mul_of_nonneg_left hmean hb

def turboMaxLiveClass22 (i : Fin 11)
    (f : State22 → ℝ) : ℝ :=
  liveStateSet.sup' liveStateSet_nonempty
    (fun q => if physicalClassAllows22 i q then f q else 0)

theorem turboMaxLiveClass22_nonneg (i : Fin 11)
    (f : State22 → ℝ) (hf : ∀ q, 0 ≤ f q) :
    0 ≤ turboMaxLiveClass22 i f := by
  obtain ⟨q, hq⟩ := liveStateSet_nonempty
  have hterm :
      0 ≤ if physicalClassAllows22 i q then f q else 0 := by
    split_ifs
    · exact hf q
    · exact le_refl 0
  exact hterm.trans (Finset.le_sup'
    (fun q => if physicalClassAllows22 i q then f q else 0) hq)

theorem turboMaxLiveClass22_bound (i : Fin 11)
    (f : State22 → ℝ) (q : State22)
    (hq : q ≠ 0) (ha : physicalClassAllows22 i q) :
    f q ≤ turboMaxLiveClass22 i f := by
  have hmem : q ∈ liveStateSet :=
    Finset.mem_erase.mpr ⟨hq, Finset.mem_univ _⟩
  simpa only [turboMaxLiveClass22, ha, ite_true] using
    (Finset.le_sup'
      (fun q => if physicalClassAllows22 i q then f q else 0)
      hmem)

def turboMeanLive22 (f : State22 → ℝ) : ℝ :=
  (∑ q ∈ liveStateSet, f q) / (liveCount State22 : ℝ)

theorem turboMeanLive22_nonneg (f : State22 → ℝ)
    (hf : ∀ q, 0 ≤ f q) : 0 ≤ turboMeanLive22 f := by
  unfold turboMeanLive22
  exact div_nonneg
    (Finset.sum_nonneg (fun q _ => hf q))
    state22_liveCount_real_pos.le

def turboLiveMomentRow22 (j : Fin 17) (z : ℝ)
    (i : Fin 11) : ℝ :=
  if i = 2 then turboMeanLive22 (turboLiveMoment22 j z)
  else turboMaxLiveClass22 i (turboLiveMoment22 j z)

def turboLiveReturnRow22 (j : Fin 17) (z : ℝ)
    (i : Fin 11) : ℝ :=
  if i = 2 then turboMeanLive22 (turboLiveReturnMoment22 j z)
  else turboMaxLiveClass22 i (turboLiveReturnMoment22 j z)

theorem turboLiveMomentRow22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) :
    0 ≤ turboLiveMomentRow22 j z i := by
  unfold turboLiveMomentRow22
  split_ifs
  · exact turboMeanLive22_nonneg _
      (turboLiveMoment22_nonneg j z hz)
  · exact turboMaxLiveClass22_nonneg i _
      (turboLiveMoment22_nonneg j z hz)

theorem turboLiveReturnRow22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) :
    0 ≤ turboLiveReturnRow22 j z i := by
  unfold turboLiveReturnRow22
  split_ifs
  · exact turboMeanLive22_nonneg _
      (turboLiveReturnMoment22_nonneg j z hz)
  · exact turboMaxLiveClass22_nonneg i _
      (turboLiveReturnMoment22_nonneg j z hz)

theorem turboLiveMomentRow22_weighted_bound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    (∑ q : State22, μ q * turboLiveMoment22 j z q) ≤
      b * turboLiveMomentRow22 j z i := by
  by_cases hi2 : i = 2
  · subst i
    change (∑ q : State22, μ q * turboLiveMoment22 j z q) ≤
      b * turboMeanLive22 (turboLiveMoment22 j z)
    exact physicalUniformPiece22_mean_bound hμ _
      (turboLiveMoment22_nonneg j z hz) _ (le_refl _)
  · simp only [turboLiveMomentRow22, hi2, ite_false]
    apply physicalLivePiece22_max_bound hi hμ _ _
      (turboMaxLiveClass22_nonneg i _
        (turboLiveMoment22_nonneg j z hz))
    intro q hq
    have hq0 : q ≠ 0 := by
      intro hzero
      subst q
      exact hq (physicalPiece22_live_zero hi hμ)
    exact turboMaxLiveClass22_bound i _ q hq0
      (hμ.2.1 q hq)

theorem turboLiveReturnRow22_weighted_bound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    (∑ q : State22, μ q * turboLiveReturnMoment22 j z q) ≤
      b * turboLiveReturnRow22 j z i := by
  by_cases hi2 : i = 2
  · subst i
    change (∑ q : State22, μ q * turboLiveReturnMoment22 j z q) ≤
      b * turboMeanLive22 (turboLiveReturnMoment22 j z)
    exact physicalUniformPiece22_mean_bound hμ _
      (turboLiveReturnMoment22_nonneg j z hz) _ (le_refl _)
  · simp only [turboLiveReturnRow22, hi2, ite_false]
    apply physicalLivePiece22_max_bound hi hμ _ _
      (turboMaxLiveClass22_nonneg i _
        (turboLiveReturnMoment22_nonneg j z hz))
    intro q hq
    have hq0 : q ≠ 0 := by
      intro hzero
      subst q
      exact hq (physicalPiece22_live_zero hi hμ)
    exact turboMaxLiveClass22_bound i _ q hq0
      (hμ.2.1 q hq)

#print axioms turboSelectedKernel22_live_entry
#print axioms turboSelectedKernel22_live_zero
#print axioms turboSelectedKernel22_live_atom_le
#print axioms turboSelectedKernel22_row_sum
#print axioms physicalLivePiece22_max_bound
#print axioms physicalUniformPiece22_mean_bound
#print axioms turboLiveMomentRow22_weighted_bound
#print axioms turboLiveReturnRow22_weighted_bound

end Spin.Turbo
