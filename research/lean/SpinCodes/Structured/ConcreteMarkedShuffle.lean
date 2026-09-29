import SpinCodes.Structured.ConcreteMarkedConditioning

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset Routing

theorem fairOnMarks_eq_mask {n : ℕ} (S : Finset (Fin n)) :
    fairOnMarks S = (iidBits n (1/2) (by norm_num) (by norm_num)).map
      (fun bits => S ∩ bits) := by
  rw [← support_coin, FinPMF.map_comp]
  unfold fairOnMarks
  congr 1
  funext bits
  ext i
  simp [Function.comp_def]

theorem iidBits_shuffle {n : ℕ} (σ : Equiv.Perm (Fin n))
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (iidBits n p hp0 hp1).map (shuffleSupport σ) = iidBits n p hp0 hp1 := by
  ext S
  change ((iidBits n p hp0 hp1).map (Equiv.finsetCongr σ.symm)).p S = _
  rw [FinPMF.map_equiv_apply]
  unfold iidBits
  rw [poissonBinom_const_apply, poissonBinom_const_apply]
  simp [Equiv.finsetCongr]

theorem shuffleSupport_inter {n : ℕ} (σ : Equiv.Perm (Fin n)) (S T : Finset (Fin n)) :
    shuffleSupport σ (S ∩ T) = shuffleSupport σ S ∩ shuffleSupport σ T := by
  ext i
  simp

theorem fairOnMarks_shuffle {n : ℕ} (σ : Equiv.Perm (Fin n)) (S : Finset (Fin n)) :
    (fairOnMarks S).map (shuffleSupport σ) = fairOnMarks (shuffleSupport σ S) := by
  rw [fairOnMarks_eq_mask, FinPMF.map_comp, fairOnMarks_eq_mask]
  conv_rhs => rw [← iidBits_shuffle σ (1/2) (by norm_num) (by norm_num)]
  rw [FinPMF.map_comp]
  congr 1
  funext bits
  exact shuffleSupport_inter σ S bits

theorem shuffledLaw_expect {n : ℕ} (P : FinPMF (Finset (Fin n)))
    (F : Finset (Fin n) → ℝ) :
    (shuffledLaw P).expect F = (FinPMF.uniform (Equiv.Perm (Fin n))).expect
      (fun σ => (P.map (shuffleSupport σ)).expect F) := by
  rw [shuffledLaw, FinPMF.expect_bind]
  simp_rw [← uniform_permutation_shuffleLaw, FinPMF.expect_map]
  simp only [FinPMF.expect, Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro σ _
  apply Finset.sum_congr rfl
  intro S _
  ring

/-- Shuffling fixed marked positions with fair bits equals first shuffling marks,
then supplying independent fair bits at the resulting positions. -/
theorem shuffled_fairOnMarks {n : ℕ} (S : Finset (Fin n)) :
    shuffledLaw (fairOnMarks S) = regionLaw S := by
  apply FinPMF.eq_of_expect_eq
  intro F
  rw [shuffledLaw_expect, regionLaw_eq_permutation, FinPMF.expect_bind]
  simp only [fairOnMarks_shuffle]

end Spin.Structured.ConcreteMarked
