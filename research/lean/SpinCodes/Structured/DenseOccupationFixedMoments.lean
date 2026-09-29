import SpinCodes.Structured.DenseOccupationFixedBasic

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric Spin.Imt SparsePolynomial

lemma sum_map_div {α : Type*} (xs : List α) (f : α → ℝ) (d : ℝ) :
    (xs.map (fun a => f a/d)).sum = (xs.map f).sum/d := by
  induction xs with
  | nil => simp
  | cons a xs ih => simp [ih, add_div]

lemma moment_mem {j : Nat} (hj : j≤128) (w : Nat) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) :
    Fix.Mem (moment w j zp) (Occupation.hyperMoment w j z) := by
  have hc : 0 < (polyChoose 128 j:Int) := by
    exact_mod_cast (show 0 < polyChoose 128 j from by rw [polyChoose_eq]; exact Nat.choose_pos hj)
  unfold moment Occupation.hyperMoment
  rw [← sum_map_div]
  apply sum_mem
  intro h _
  dsimp
  by_cases h1 : h≤j
  · by_cases h2 : polyChoose w h * polyChoose (128-w) (j-h) = 0
    · simp only [h1, h2, ne_eq, not_true_eq_false, and_false, ↓reduceIte]
      have hh : (w.choose h:ℝ)*((128-w).choose (j-h):ℝ) = 0 := by
        exact_mod_cast (show w.choose h * (128-w).choose (j-h)=0 by simpa [polyChoose_eq] using h2)
      simp only [hh, zero_mul, zero_div]
      exact zero_mem
    · simp only [h1, h2, ne_eq, not_false_eq_true, and_self, ↓reduceIte]
      have hh := Fix.mul_mem (Fix.ofFrac_mem (p := (polyChoose w h * polyChoose (128-w) (j-h):Nat)) hc)
        (hz (w+j-2*h))
      convert hh using 1 <;> push_cast [polyChoose_eq] <;> ring
  · simp [h1, zero_mem]

lemma pattern_mem {j : Nat} (hj : j≤128) (pat : List (Nat × Nat))
    (zp : Nat → Fix) (z : ℝ) (hz : ∀ n, Fix.Mem (zp n) (z^n)) :
    Fix.Mem (pattern j pat zp) (Occupation.Sparse.pattern j pat z) := by
  have hc : 0 < (polyChoose 128 j:Int) := by
    exact_mod_cast (show 0 < polyChoose 128 j from by rw [polyChoose_eq]; exact Nat.choose_pos hj)
  unfold pattern Occupation.Sparse.pattern
  rw [← sum_map_div]
  apply sum_mem
  rintro ⟨out,c⟩ _
  have hh := Fix.mul_mem (Fix.ofFrac_mem (p := (c:Int)) hc) (hz out)
  convert hh using 1 <;> push_cast [polyChoose_eq] <;> ring

lemma live_mem {j : Nat} (hj : j≤128) (d : WeightData) :
    Fix.Mem (live j d) (Occupation.Sparse.live j d) := by
  have hc : 0 < (polyChoose 128 j:Int) := by
    exact_mod_cast (show 0 < polyChoose 128 j from by rw [polyChoose_eq]; exact Nat.choose_pos hj)
  have h := Fix.ofFrac_mem (p := (polyChoose 128 j:Int)-d.kernel) hc
  simpa [live, Occupation.Sparse.live, polyChoose_eq] using h

lemma maximumMoment_mem {j : Nat} (hj : j≤128) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) :
    Fix.Mem (maximumMoment j zp) (Occupation.Sparse.maximumMoment j z) :=
  fold_max_mem _ _ _ zero_mem (fun w _ => moment_mem hj w zp z hz)

end Spin.Structured.DenseOccupationFixed
