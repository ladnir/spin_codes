import SpinCodes.UniformMixer.Map22PacketProfileSmallCounts

noncomputable section
namespace Spin.UniformMixer

theorem profile_sum_by_weight {A M : Type*} [Fintype A] [AddCommMonoid M]
    (n : Nat) (weight : A → Nat) (bound : ∀ a, weight a ≤ n) (f : Nat → M) :
    (∑ a, f (weight a)) = ∑ w : Fin (n+1),
      (∑ a, if weight a = w.val then 1 else 0) • f w.val := by
  classical
  calc
    _ = ∑ a, ∑ w : Fin (n+1), if weight a = w.val then f w.val else 0 := by
      apply Finset.sum_congr rfl
      intro a _
      symm
      rw [Finset.sum_eq_single (⟨weight a, by have := bound a; omega⟩ : Fin (n+1))]
      · simp
      · intro w _ h
        have hh : weight a ≠ w.val := by
          intro he
          apply h
          exact Fin.ext he.symm
        simp [hh]
      · simp
    _ = _ := by
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro w _
      rw [← Finset.sum_nsmul_assoc]
      apply Finset.sum_congr rfl
      intro a _
      by_cases h : weight a = w.val <;> simp [h]

theorem profileQuad4Weight_bound (q : ProfileQuad4) : profileQuad4Weight q ≤ 16 := by
  unfold profileQuad4Weight
  calc
    _ ≤ ∑ _v : Fin 16, 1 := Finset.sum_le_sum (fun v _ => by split_ifs <;> omega)
    _ = 16 := by simp

theorem profileQuad3Weight_bound (q : ProfileQuad3) : profileQuad3Weight q ≤ 8 := by
  unfold profileQuad3Weight
  calc
    _ ≤ ∑ _v : Fin 8, 1 := Finset.sum_le_sum (fun v _ => by split_ifs <;> omega)
    _ = 8 := by simp

theorem profileQuad2Weight_bound (q : ProfileQuad2) : profileQuad2Weight q ≤ 4 := by
  unfold profileQuad2Weight
  calc
    _ ≤ ∑ _v : Fin 4, 1 := Finset.sum_le_sum (fun v _ => by split_ifs <;> omega)
    _ = 4 := by simp

theorem profileQuad4_weight_sum {M : Type*} [AddCommMonoid M] (f : Nat → M) :
    (∑ q : ProfileQuad4, f (profileQuad4Weight q)) =
      ∑ w : Fin 17, profileQuad4WeightCounts.getD w.val 0 • f w.val := by
  rw [profile_sum_by_weight 16 profileQuad4Weight profileQuad4Weight_bound]
  simp_rw [profileQuad4_weight_count]

theorem profileQuad3_weight_sum {M : Type*} [AddCommMonoid M] (f : Nat → M) :
    (∑ q : ProfileQuad3, f (profileQuad3Weight q)) =
      ∑ w : Fin 9, profileQuad3WeightCounts.getD w.val 0 • f w.val := by
  rw [profile_sum_by_weight 8 profileQuad3Weight profileQuad3Weight_bound]
  simp_rw [profileQuad3_weight_count]

theorem profileQuad2_weight_sum {M : Type*} [AddCommMonoid M] (f : Nat → M) :
    (∑ q : ProfileQuad2, f (profileQuad2Weight q)) =
      ∑ w : Fin 5, profileQuad2WeightCounts.getD w.val 0 • f w.val := by
  rw [profile_sum_by_weight 4 profileQuad2Weight profileQuad2Weight_bound]
  simp_rw [profileQuad2_weight_count]

end Spin.UniformMixer
