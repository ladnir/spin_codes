import SpinCodes.Structured.DenseOccupationFixedDefs
import SpinCodes.Structured.DenseOccupationSound
import SpinCodes.Numeric.Fixed

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric Spin.Imt

lemma imin_cast (a b : Int) : ((imin a b:Int):ℝ) = min (a:ℝ) (b:ℝ) := by
  unfold imin
  split_ifs with h
  · simp [min_eq_left (show (a:ℝ)≤b by exact_mod_cast h)]
  · simp [min_eq_right (show (b:ℝ)≤a by exact_mod_cast (le_of_not_ge h))]

lemma imax_cast (a b : Int) : ((imax a b:Int):ℝ) = max (a:ℝ) (b:ℝ) := by
  unfold imax
  split_ifs with h
  · simp [max_eq_right (show (a:ℝ)≤b by exact_mod_cast h)]
  · simp [max_eq_left (show (b:ℝ)≤a by exact_mod_cast (le_of_not_ge h))]

lemma min_mem {a b : Fix} {x y : ℝ} (hx : Fix.Mem a x) (hy : Fix.Mem b y) :
    Fix.Mem (iminFix a b) (min x y) := by
  dsimp [Fix.Mem, iminFix] at *
  rw [imin_cast, imin_cast, ← min_div_div_right scaleR_pos.le, ← min_div_div_right scaleR_pos.le]
  exact ⟨min_le_min hx.1 hy.1, min_le_min hx.2 hy.2⟩

lemma max_mem {a b : Fix} {x y : ℝ} (hx : Fix.Mem a x) (hy : Fix.Mem b y) :
    Fix.Mem (imaxFix a b) (max x y) := by
  dsimp [Fix.Mem, imaxFix] at *
  rw [imax_cast, imax_cast, ← max_div_div_right scaleR_pos.le, ← max_div_div_right scaleR_pos.le]
  exact ⟨max_le_max hx.1 hy.1, max_le_max hx.2 hy.2⟩

lemma zero_mem : Fix.Mem Fix.zero 0 := by norm_num [Fix.Mem, Fix.zero]

lemma sum_mem {α : Type*} (xs : List α) (f : α → Fix) (g : α → ℝ)
    (h : ∀ a ∈ xs, Fix.Mem (f a) (g a)) : Fix.Mem (isum (xs.map f)) ((xs.map g).sum) := by
  induction xs with
  | nil => exact zero_mem
  | cons a xs ih =>
    exact Fix.add_mem (h a (by simp)) (ih (fun a ha => h a (by simp [ha])))

lemma fold_max_mem {α : Type*} (xs : List α) (f : α → Fix) (g : α → ℝ)
    {b : Fix} {v : ℝ} (hb : Fix.Mem b v) (h : ∀ a ∈ xs, Fix.Mem (f a) (g a)) :
    Fix.Mem ((xs.map f).foldr imaxFix b) ((xs.map g).foldr max v) := by
  induction xs with
  | nil => exact hb
  | cons a xs ih => exact max_mem (h a (by simp)) (ih (fun a ha => h a (by simp [ha])))

lemma fold_min_mem {α : Type*} (xs : List α) (f : α → Fix) (g : α → ℝ)
    {b : Fix} {v : ℝ} (hb : Fix.Mem b v) (h : ∀ a ∈ xs, Fix.Mem (f a) (g a)) :
    Fix.Mem ((xs.map f).foldr iminFix b) ((xs.map g).foldr min v) := by
  induction xs with
  | nil => exact hb
  | cons a xs ih => exact min_mem (h a (by simp)) (ih (fun a ha => h a (by simp [ha])))

lemma powers_mem {z : Fix} {x : ℝ} (hz : Fix.Mem z x) {xs : List Fix} {N : Nat}
    (h : checkPowers z xs N = true) (k : Nat) : Fix.Mem (powers z xs N k) (x^k) := by
  simp only [checkPowers, Bool.and_eq_true, beq_iff_eq, List.all_eq_true] at h
  have ht : ∀ k, k≤N → Fix.Mem (xs.getD k Fix.zero) (x^k) := by
    intro k
    induction k with
    | zero => intro _; rw [h.1]; simpa using Fix.ofInt_mem 1
    | succ k ih =>
      intro hk
      have hh := h.2 k (List.mem_range.mpr (by omega))

      rw [hh, pow_succ]
      exact Fix.mul_mem (ih (by omega)) hz
  unfold powers
  split_ifs with hk
  · exact ht k hk
  · exact Fix.pow_mem hz k

end Spin.Structured.DenseOccupationFixed

