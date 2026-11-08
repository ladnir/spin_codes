import SpinCodes.Turbo.PhysicalLive22

/-!
# Complete live-state row for the Turbo physical kernel

The row keeps the exact feedback-conditioned return to zero. All remaining
output mass is represented by the uniformly dominated live coordinate.
This is a symbolic soundness result; it does not authenticate the numeric
producer's coefficients.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer
attribute [local instance] Classical.propDecidable

def turboLiveOutput22 (j : Fin 17) (z : ℝ)
    (μ : State22 → ℝ) : State22 → ℝ :=
  physicalKernelAction22 (Matrix.of (turboSelectedKernel22 j z)) μ

theorem turboSelectedKernel22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (q r : State22) :
    0 ≤ turboSelectedKernel22 j z q r := by
  unfold turboSelectedKernel22
  apply FinPMF.expect_nonneg
  intro x
  exact mul_nonneg (pow_nonneg hz _)
    (update_nonneg (stepFeedback (packetStepInput22 x)) q r)

theorem turboLiveOutput22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (μ : State22 → ℝ)
    (hμ : ∀ q, 0 ≤ μ q) (r : State22) :
    0 ≤ turboLiveOutput22 j z μ r := by
  unfold turboLiveOutput22 physicalKernelAction22
  change 0 ≤ ∑ q : State22, μ q * turboSelectedKernel22 j z q r
  apply Finset.sum_nonneg
  intro q _
  exact mul_nonneg (hμ q) (turboSelectedKernel22_nonneg j z hz q r)

theorem turboLiveOutput22_zero (j : Fin 17)
    (z : ℝ) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    turboLiveOutput22 j z μ 0 =
      (∑ q : State22, μ q * turboLiveReturnMoment22 j z q) /
        (liveCount State22 : ℝ) := by
  unfold turboLiveOutput22 physicalKernelAction22
  change (∑ q : State22, μ q * turboSelectedKernel22 j z q 0) = _
  calc
    _ = ∑ q : State22,
          (μ q * turboLiveReturnMoment22 j z q) /
            (liveCount State22 : ℝ) := by
      apply Finset.sum_congr rfl
      intro q _
      by_cases hq : q = 0
      · subst q
        rw [physicalPiece22_live_zero hi hμ]
        ring
      · rw [turboSelectedKernel22_live_zero j z q hq]
        ring
    _ = _ := by rw [Finset.sum_div]

theorem turboLiveOutput22_mass (j : Fin 17)
    (z : ℝ) (μ : State22 → ℝ) :
    (∑ r : State22, turboLiveOutput22 j z μ r) =
      ∑ q : State22, μ q * turboLiveMoment22 j z q := by
  unfold turboLiveOutput22 physicalKernelAction22
  rw [Finset.sum_comm]
  change (∑ q : State22,
    ∑ r : State22, μ q * turboSelectedKernel22 j z q r) = _
  simp only [← Finset.mul_sum, turboSelectedKernel22_row_sum]

theorem turboLiveOutput22_atom_le (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) (r : State22) :
    turboLiveOutput22 j z μ r ≤
      (b * turboLiveMomentRow22 j z i) /
        (liveCount State22 : ℝ) := by
  have hden : 0 ≤ (liveCount State22 : ℝ) :=
    state22_liveCount_real_pos.le
  calc
    turboLiveOutput22 j z μ r ≤
        ∑ q : State22,
          μ q * (turboLiveMoment22 j z q /
            (liveCount State22 : ℝ)) := by
      unfold turboLiveOutput22 physicalKernelAction22
      change (∑ q : State22,
        μ q * turboSelectedKernel22 j z q r) ≤ _
      apply Finset.sum_le_sum
      intro q _
      by_cases hq : q = 0
      · subst q
        rw [physicalPiece22_live_zero hi hμ]
        simpa only [zero_mul] using (le_refl (0 : ℝ))
      · exact mul_le_mul_of_nonneg_left
          (turboSelectedKernel22_live_atom_le j z hz q r hq)
          (hμ.1 q)
    _ = (∑ q : State22, μ q * turboLiveMoment22 j z q) /
          (liveCount State22 : ℝ) := by
      rw [Finset.sum_div]
      apply Finset.sum_congr rfl
      intro q _
      ring
    _ ≤ (b * turboLiveMomentRow22 j z i) /
          (liveCount State22 : ℝ) := by
      exact div_le_div_of_nonneg_right
        (turboLiveMomentRow22_weighted_bound j z hz i hi b μ hμ)
        hden

theorem turboLiveOutput22_zero_le (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    turboLiveOutput22 j z μ 0 ≤
      (b * turboLiveReturnRow22 j z i) /
        (liveCount State22 : ℝ) := by
  rw [turboLiveOutput22_zero j z i hi b μ hμ]
  exact div_le_div_of_nonneg_right
    (turboLiveReturnRow22_weighted_bound j z hz i hi b μ hμ)
    state22_liveCount_real_pos.le

theorem turboLiveOutput22_mass_le (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    (∑ r : State22, turboLiveOutput22 j z μ r) ≤
      b * turboLiveMomentRow22 j z i := by
  rw [turboLiveOutput22_mass]
  exact turboLiveMomentRow22_weighted_bound j z hz i hi b μ hμ

def turboLiveZeroPart22 (j : Fin 17) (z : ℝ)
    (μ : State22 → ℝ) (r : State22) : ℝ :=
  if r = 0 then turboLiveOutput22 j z μ r else 0

def turboLiveRestPart22 (j : Fin 17) (z : ℝ)
    (μ : State22 → ℝ) (r : State22) : ℝ :=
  if r = 0 then 0 else turboLiveOutput22 j z μ r

theorem turboLiveZeroPart22_piece (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    physicalPiece22 0
      (b * (turboLiveReturnRow22 j z i /
        (liveCount State22 : ℝ)))
      (turboLiveZeroPart22 j z μ) := by
  constructor
  · intro r
    unfold turboLiveZeroPart22
    split_ifs
    · exact turboLiveOutput22_nonneg j z hz μ hμ.1 r
    · exact le_refl 0
  constructor
  · intro r hr
    have hr0 : r = 0 := by
      by_contra hne
      exact hr (by simp [turboLiveZeroPart22, hne])
    simpa [physicalClassAllows22] using hr0
  constructor
  · have hsum :
        (∑ r : State22, turboLiveZeroPart22 j z μ r) =
          turboLiveOutput22 j z μ 0 := by
      simp [turboLiveZeroPart22]
    rw [hsum]
    convert turboLiveOutput22_zero_le j z hz i hi b μ hμ using 1 <;> ring
  · intro h
    omega

theorem turboLiveRestPart22_piece (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    physicalPiece22 2 (b * turboLiveMomentRow22 j z i)
      (turboLiveRestPart22 j z μ) := by
  constructor
  · intro r
    unfold turboLiveRestPart22
    split_ifs
    · exact le_refl 0
    · exact turboLiveOutput22_nonneg j z hz μ hμ.1 r
  constructor
  · intro r hr
    have hr0 : r ≠ 0 := by
      intro he
      exact hr (by simp [turboLiveRestPart22, he])
    simpa [physicalClassAllows22] using hr0
  constructor
  · calc
      (∑ r : State22, turboLiveRestPart22 j z μ r) ≤
          ∑ r : State22, turboLiveOutput22 j z μ r := by
        apply Finset.sum_le_sum
        intro r _
        by_cases hr : r = 0
        · subst r
          simp only [turboLiveRestPart22, ite_true]
          exact turboLiveOutput22_nonneg j z hz μ hμ.1 0
        · simp only [turboLiveRestPart22, hr, ite_false]
          exact le_refl _
      _ ≤ _ := turboLiveOutput22_mass_le j z hz i hi b μ hμ
  · intro _ r
    have hcap := turboLiveOutput22_atom_le j z hz i hi b μ hμ r
    by_cases hr : r = 0
    · subst r
      simp only [turboLiveRestPart22, ite_true]
      exact div_nonneg
        (mul_nonneg
          ((Finset.sum_nonneg (fun q _ => hμ.1 q)).trans hμ.2.2.1)
          (turboLiveMomentRow22_nonneg j z hz i))
        (by norm_num)
    · simpa only [turboLiveRestPart22, hr, ite_false,
        physicalLiveDenom22] using hcap

def turboLiveRow22 (j : Fin 17) (z : ℝ)
    (i k : Fin 11) : ℝ :=
  if k = 0 then
    turboLiveReturnRow22 j z i / (liveCount State22 : ℝ)
  else if k = 2 then turboLiveMomentRow22 j z i else 0

theorem turboSelectedKernel22_live_row_sound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 i b μ) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, turboLiveOutput22 j z μ r =
        ∑ k : Fin 11, parts k r) ∧
      ∀ k, physicalPiece22 k
        (b * turboLiveRow22 j z i k) (parts k) := by
  let parts : Fin 11 → State22 → ℝ := fun k r =>
    if k = 0 then turboLiveZeroPart22 j z μ r
    else if k = 2 then turboLiveRestPart22 j z μ r else 0
  refine ⟨parts, ?_, ?_⟩
  · intro r
    simp only [parts]
    by_cases hr : r = 0
    · subst r
      simp [turboLiveZeroPart22, turboLiveRestPart22]
    · simp only [turboLiveZeroPart22, turboLiveRestPart22,
        hr, ite_false]
      symm
      apply Finset.sum_eq_single 2
      · intro k _ hk2
        by_cases hk0 : k = 0 <;> simp [hk0, hk2]
      · intro h
        exact (h (Finset.mem_univ 2)).elim
  · intro k
    by_cases hk0 : k = 0
    · subst k
      simpa only [parts, turboLiveRow22, ite_true]
        using turboLiveZeroPart22_piece j z hz i hi b μ hμ
    · by_cases hk2 : k = 2
      · subst k
        simpa only [parts, turboLiveRow22, hk0, ite_false,
          ite_true] using
          turboLiveRestPart22_piece j z hz i hi b μ hμ
      · have hzrow : b * turboLiveRow22 j z i k = 0 := by
          simp only [turboLiveRow22, hk0, hk2, ite_false, mul_zero]
        simpa only [parts, hk0, hk2, ite_false, hzrow] using
          physicalPiece22_zero k 0 (le_refl 0)

#print axioms turboLiveOutput22_nonneg
#print axioms turboLiveOutput22_zero
#print axioms turboLiveOutput22_mass
#print axioms turboLiveOutput22_atom_le
#print axioms turboLiveOutput22_zero_le
#print axioms turboLiveOutput22_mass_le
#print axioms turboLiveZeroPart22_piece
#print axioms turboLiveRestPart22_piece
#print axioms turboSelectedKernel22_live_row_sound

end Spin.Turbo
