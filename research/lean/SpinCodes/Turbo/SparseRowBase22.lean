import SpinCodes.Turbo.SparseRowBaseDefs22
import SpinCodes.Turbo.SparseRegionalPower22
import SpinCodes.Turbo.PhysicalProfileBirthFamily22
import SpinCodes.Turbo.PhysicalProfileCompleteRows22

/-! Physical rows j ≤ 2 of a retained local family are exact integer
fractions over the endpoint denominator. Hence an integer base stage at
scale den·mult can be checked exactly against choose(16,d)·B_d. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22

namespace SparseRowBase22

theorem profileFamilyRows22_small_eq (small : Fin 3 → Fin 11 → ℝ)
    (first returns : List (List Int)) (zeroHold : List Int) (den : Int)
    (births : Fin 3 → List Int)
    (hsmall : ∀ j k, small j k = ((births j).getD k.val 0 : ℝ) / den)
    (j : ℕ) (hj : j < 3) (i k : Fin 11) :
    profileFamilyRows22 (profileFamilyBirth22 small first zeroHold den) first returns den
        ⟨j, by omega⟩ i k =
      (profileRowInt22 (births ⟨j, hj⟩) (first.getD j []) (returns.getD j [])
        i.val k.val : ℝ) / den := by
  unfold profileFamilyRows22 profileOutwardRows22 profileFamilyBirth22 profileEndpoint22
    profileRowInt22
  simp only [Matrix.of_apply]
  rw [dif_pos (show (⟨j, by omega⟩ : Fin 17).val < 3 from hj)]
  by_cases hi : i = 0
  · subst hi
    simp only [if_true, Fin.val_zero]
    rw [hsmall]
  · have hi' : i.val ≠ 0 := fun h => hi (Fin.ext h)
    rw [if_neg hi, if_neg hi']
    by_cases hk : k = 0
    · subst hk
      simp only [if_true, Fin.val_zero]
    · have hk' : k.val ≠ 0 := fun h => hk (Fin.ext h)
      rw [if_neg hk, if_neg hk']
      by_cases hk2 : k = 2
      · subst hk2
        simp only [if_true]
        rfl
      · have hk2' : k.val ≠ 2 := fun h => hk2 (Fin.ext h)
        rw [if_neg hk2, if_neg hk2']
        simp

theorem binom16_eq (d : ℕ) (hd : d ≤ 2) : (sparseBinom16Small22 d : ℝ) = Nat.choose 16 d := by
  interval_cases d <;> norm_num [sparseBinom16Small22, Nat.choose]

theorem base_check_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (row : ℕ → ℕ → ℕ → ℤ) (den mult : ℤ) (hden : 0 < den) (hmult : 0 < mult)
    (hB : ∀ d (hd : d ≤ 2), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k = (row d i.val k.val : ℝ) / den)
    (base : Array ℤ) (hcheck : sparseRegionalBaseCheck22 row mult base = true) :
    ∀ d (hd : d ≤ 2), ∀ i k : Fin 11, (Nat.choose 16 d : ℝ) * B ⟨d, by omega⟩ i k ≤
      (sparseIntPolyMatrix22 base d i k : ℝ) / ((den * mult : ℤ) : ℝ) := by
  intro d hd i k
  have hc := of_decide_eq_true (all11 _ (all11 _ (all3 _ hcheck d hd) i) k)
  have hcR : (sparseBinom16Small22 d : ℝ) * (row d i.val k.val : ℝ) * (mult : ℝ) ≤
      (sparseIntPolyMatrix22 base d i k : ℝ) := by exact_mod_cast hc
  rw [binom16_eq d hd] at hcR
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have hmR : (0 : ℝ) < mult := by exact_mod_cast hmult
  rw [hB d hd i k, Int.cast_mul, mul_div_assoc', div_le_div_iff₀ hdR (mul_pos hdR hmR)]
  nlinarith [mul_le_mul_of_nonneg_right hcR hdR.le]

end SparseRowBase22

open SparseRowBase22 in
#print axioms profileFamilyRows22_small_eq
open SparseRowBase22 in
#print axioms base_check_sound
end Spin.Turbo
