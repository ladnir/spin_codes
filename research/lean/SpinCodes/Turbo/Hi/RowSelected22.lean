import SpinCodes.Turbo.Hi.PhysicalProfileBirthFamily22
import SpinCodes.Turbo.Hi.PhysicalProfileCompleteRows22
import SpinCodes.Turbo.PhysicalThinningRows22

/-! Integer form of every physical row (including the coarse birth rows,
j ≥ 3) and of every selected (thinned) row of a high-precision local family:
selected_d = Σ_x C(d,x)·15^x/16^d · physical_x, an exact integer over
16^16 · den. -/
noncomputable section
namespace Spin.Turbo.Hi
open Finset

def profileRowIntAll22 (births : Nat → List Int) (first returns : List (List Int))
    (zeroHold : List Int) (j i k : Nat) : Int :=
  if i = 0 then
    (if j < 3 then (births j).getD k 0
     else if k = 0 then zeroHold.getD j 0 else if k = 1 then (first.getD j []).getD 0 0 else 0)
  else if k = 0 then (returns.getD j []).getD i 0
  else if k = 2 then (first.getD j []).getD i 0 else 0

theorem rows_eq (small : Fin 3 → Fin 11 → ℝ) (first returns : List (List Int))
    (zeroHold : List Int) (den : Int) (births : Nat → List Int)
    (hsmall : ∀ (j : Fin 3) (k : Fin 11), small j k = ((births j.val).getD k.val 0 : ℝ) / den)
    (j : Fin 17) (i k : Fin 11) :
    profileFamilyRows22 (profileFamilyBirth22 small first zeroHold den) first returns den j i k =
      (profileRowIntAll22 births first returns zeroHold j.val i.val k.val : ℝ) / den := by
  unfold profileFamilyRows22 profileOutwardRows22 profileEndpoint22 profileRowIntAll22
  simp only [Matrix.of_apply]
  by_cases hi : i = 0
  · subst hi
    simp only [if_true, Fin.val_zero]
    unfold profileFamilyBirth22
    by_cases hj : j.val < 3
    · rw [dif_pos hj, if_pos hj, hsmall]
    · rw [dif_neg hj, if_neg hj]
      unfold physicalCoarseBirthRow22
      by_cases hk0 : k = 0
      · subst hk0
        simp
      · have hk0' : k.val ≠ 0 := fun h => hk0 (Fin.ext h)
        rw [if_neg hk0, if_neg hk0']
        by_cases hk1 : k = 1
        · subst hk1
          simp
        · have hk1' : k.val ≠ 1 := fun h => hk1 (Fin.ext h)
          rw [if_neg hk1, if_neg hk1']
          simp
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

def pascal17 : List (List Int) := [[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  [1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  [1,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
  [1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0],
  [1,5,10,10,5,1,0,0,0,0,0,0,0,0,0,0,0],
  [1,6,15,20,15,6,1,0,0,0,0,0,0,0,0,0,0],
  [1,7,21,35,35,21,7,1,0,0,0,0,0,0,0,0,0],
  [1,8,28,56,70,56,28,8,1,0,0,0,0,0,0,0,0],
  [1,9,36,84,126,126,84,36,9,1,0,0,0,0,0,0,0],
  [1,10,45,120,210,252,210,120,45,10,1,0,0,0,0,0,0],
  [1,11,55,165,330,462,462,330,165,55,11,1,0,0,0,0,0],
  [1,12,66,220,495,792,924,792,495,220,66,12,1,0,0,0,0],
  [1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1,0,0,0],
  [1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1,0,0],
  [1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365,455,105,15,1,0],
  [1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1]]

theorem pascal_eq : ∀ d, d < 17 → ∀ x, x < 17 →
    (pascal17.getD d []).getD x 0 = ((Nat.choose d x : ℕ) : Int) := by
  decide +kernel

def rowSel (physInt : Nat → Nat → Nat → Int) (d i k : Nat) : Int :=
  ∑ x : Fin 17, (pascal17.getD d []).getD x.val 0 * 15 ^ x.val * 16 ^ (16 - d) *
    physInt x.val i k

theorem selected_eq (z : ℝ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (physInt : Nat → Nat → Nat → Int) (den : Int) (hden : 0 < den)
    (hB : ∀ (x : Fin 17) (i k : Fin 11), B x i k = (physInt x.val i.val k.val : ℝ) / den)
    (d : ℕ) (hd : d < 17) (i k : Fin 11) :
    turboThinnedRow22 ⟨d, hd⟩ z B i k =
      (rowSel physInt d i.val k.val : ℝ) / ((16 : ℝ) ^ 16 * den) := by
  unfold turboThinnedRow22 physicalWeightedRow22 thinningWeight22 rowSel
  simp only [Matrix.of_apply, hB]
  rw [Int.cast_sum, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro x _
  rw [pascal_eq d hd x.val x.isLt]
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have h16 : (16 : ℝ) ^ d * (16 : ℝ) ^ (16 - d) = (16 : ℝ) ^ 16 := by
    rw [← pow_add, Nat.add_sub_cancel' (by omega)]
  have h16d : (0 : ℝ) < (16 : ℝ) ^ d := by positivity
  push_cast
  rw [← h16]
  field_simp

end Spin.Turbo.Hi

#print axioms Spin.Turbo.Hi.rows_eq
#print axioms Spin.Turbo.Hi.selected_eq
