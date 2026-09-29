import SpinCodes.Structured.SparseSelection
import SpinCodes.Structured.SparseValidationDefs

/-! The checked branch indices select upper bounds for the cancellation
minima in the fixed-weight occupation matrix. -/

noncomputable section

namespace Spin.Structured.SparsePolynomial

theorem getD_den_pos (ps : List RatPoly) (fallback : RatPoly) (i : ℕ)
    (hf : 0 < fallback.den) (hp : ∀ p ∈ ps, 0 < p.den) :
    0 < (ps.getD i fallback).den := by
  induction ps generalizing i with
  | nil => simpa using hf
  | cons p ps ih =>
    cases i with
    | zero => exact hp p (by simp)
    | succ i => exact ih i (fun q hq => hp q (by simp [hq]))

theorem pattern_den_pos {j : ℕ} (hj : j ≤ 128) (p : List (ℕ × ℕ)) :
    0 < (pattern j p).den := by
  apply scale_den_pos _ (choose128_pos hj)
  apply RatPoly.sum_den_pos
  intro q hq
  obtain ⟨⟨out, c⟩, _, rfl⟩ := List.mem_map.mp hq
  exact scale_den_pos _ (by decide) (SparsePowers.zPow_den_pos out)

theorem arbitrary_den_pos {j : ℕ} (hj : j ≤ 128) (d : WeightData) :
    0 < (arbitrary j d).den := by
  apply getD_den_pos _ _ _ (by decide)
  intro p hp
  obtain ⟨w, _, rfl⟩ := List.mem_map.mp hp
  exact moment_den_pos w hj

theorem counts_getD_pos (i : ℕ) : 0 < counts.getD i 1 := by
  have h : ∀ c ∈ counts, 0 < c := by simp [counts]
  generalize counts = cs at *
  induction cs generalizing i with
  | nil => simp
  | cons c cs ih =>
    cases i with
    | zero => exact h c (by simp)
    | succ i => exact ih i (fun a ha => h a (by simp [ha]))

theorem boundsD_den_pos {j : ℕ} (hj : j ≤ 128) (d : WeightData) :
    ∀ p ∈ boundsD j d, 0 < p.den := by
  intro p hp
  simp only [boundsD, List.mem_append, List.mem_cons, List.not_mem_nil, or_false] at hp
  rcases hp with (rfl | rfl | rfl) | hp
  · exact arbitrary_den_pos hj d
  · exact live_den_pos hj d
  · exact scale_den_pos _ (choose128_pos hj) (SparsePowers.zPow_den_pos _)
  · cases h : d.lowPatterns with
    | none => simp [h] at hp
    | some ps =>
      simp only [h, List.mem_singleton] at hp
      subst p
      apply getD_den_pos _ _ _ (by decide)
      intro q hq
      obtain ⟨p, _, rfl⟩ := List.mem_map.mp hq
      exact pattern_den_pos hj p

theorem boundsS_den_pos {j : ℕ} (hj : j ≤ 128) (d : WeightData) (i : ℕ) :
    ∀ p ∈ boundsS j d i, 0 < p.den := by
  intro p hp
  simp only [boundsS, List.mem_append, List.mem_cons, List.not_mem_nil, or_false] at hp
  rcases hp with (rfl | rfl) | hp
  · exact moment_den_pos _ hj
  · exact scale_den_pos _ (Nat.mul_pos (counts_getD_pos i) (choose128_pos hj))
      (SparsePowers.zPow_den_pos _)
  · cases h : d.lowPatterns with
    | none => simp [h] at hp
    | some ps =>
      simp only [h, List.mem_singleton] at hp
      subst p
      exact scale_den_pos _ (counts_getD_pos i) (pattern_den_pos hj _)

theorem selected_min_le (xs : List ℝ) (ps : List RatPoly) (base : ℝ)
    (x : ℝ) (i : ℕ) (h : List.Forall₂ (fun a p => a ≤ p.eval x) xs ps)
    (hi : i < xs.length) : xs.foldr min base ≤ (ps.getD i zero).eval x := by
  induction h generalizing i with
  | nil => simp at hi
  | @cons a p xs ps ha ht ih =>
    cases i with
    | zero => exact (min_le_left _ _).trans ha
    | succ i =>
      exact (min_le_right _ _).trans (ih i (by simpa using hi))

theorem boundsD_compare (j : ℕ) (d : WeightData) (x : ℝ)
    (ha : Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j d).eval x)
    (hlow : lowMaximum j d (1 - x / 6250) ≤ (selectedLowPolynomial j d).eval x) :
    List.Forall₂ (fun a p => a ≤ p.eval x)
      (Spin.Imt.Occupation.Sparse.cancellationBoundsD j d (1 - x / 6250))
      (boundsD j d) := by
  have hl : Spin.Imt.Occupation.Sparse.live j d ≤ (live j d).eval x :=
    (eval_live j d x).ge
  have hc : (d.cap : ℝ) / Nat.choose 128 j * (1 - x / 6250) ^ minDistance j ≤
      (scale d.cap (polyChoose 128 j) (SparsePowers.zPow (minDistance j))).eval x := by
    simp [eval_scale, SparsePowers.zPow_eval, polyChoose_eq]
  cases h : d.lowPatterns with
  | none => simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsD, boundsD, h]
      using List.Forall₂.cons ha (List.Forall₂.cons hl (List.Forall₂.cons hc .nil))
  | some ps =>
    have hh := hlow
    simp only [lowMaximum, selectedLowPolynomial, h] at hh
    simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsD, boundsD, h]
      using List.Forall₂.cons ha (List.Forall₂.cons hl
        (List.Forall₂.cons hc (List.Forall₂.cons hh .nil)))

theorem boundsS_compare (j : ℕ) (d : WeightData) (i : Fin 5) (x : ℝ)
    (hk : d.kernel ≤ polyChoose 128 j) :
    List.Forall₂ (fun a p => a ≤ p.eval x)
      (Spin.Imt.Occupation.Sparse.cancellationBoundsS j d (1 - x / 6250) i)
      (boundsS j d i) := by
  have hm := (eval_moment (levels.getD i 0) j x).ge
  have hk' : d.kernel ≤ Nat.choose 128 j := by simpa [polyChoose_eq] using hk
  have he : ((min (Nat.choose 128 j - d.kernel) (counts.getD i 1 * d.cap) : ℕ) : ℝ) =
      ((min ((polyChoose 128 j : ℤ) - d.kernel)
        ((counts.getD i 1 : ℤ) * d.cap) : ℤ) : ℝ) := by
    simp [polyChoose_eq, Nat.cast_min, Nat.cast_sub hk', Int.cast_min]
  have hc :
      (min (Nat.choose 128 j - d.kernel) (counts.getD i 1 * d.cap) : ℕ) /
          ((counts.getD i 1 : ℝ) * Nat.choose 128 j) *
          (1 - x / 6250) ^ distance (levels.getD i 0) j ≤
      (scale (min ((polyChoose 128 j : ℤ) - d.kernel)
          ((counts.getD i 1 : ℤ) * d.cap)) (counts.getD i 1 * polyChoose 128 j)
          (SparsePowers.zPow (distance (levels.getD i 0) j))).eval x := by
    rw [eval_scale, SparsePowers.zPow_eval, he]
    simp [Nat.cast_mul, polyChoose_eq]
  have hl : Spin.Imt.Occupation.Sparse.pattern j (d.lowShells.getD i [])
      (1 - x / 6250) / counts.getD i 1 ≤
      (scale 1 (counts.getD i 1) (pattern j (d.lowShells.getD i []))).eval x := by
    rw [eval_scale, eval_pattern]
    simp only [Int.cast_one]
    exact le_of_eq (by ring)
  cases h : d.lowPatterns with
  | none => simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsS, boundsS, h]
      using List.Forall₂.cons hm (List.Forall₂.cons hc .nil)
  | some ps => simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsS, boundsS, h]
      using List.Forall₂.cons hm (List.Forall₂.cons hc (List.Forall₂.cons hl .nil))

theorem checkData_ranges {j : ℕ} {d : WeightData} (h : checkData j d = true) :
    d.kernel ≤ polyChoose 128 j ∧ d.choices.freshD < 2 ∧
    d.choices.cancelD < 3 + (if d.lowPatterns.isSome then 1 else 0) ∧
    ∀ i : Fin 5, d.choices.cancelS.getD i 0 <
        2 + (if d.lowPatterns.isSome then 1 else 0) ∧ d.choices.freshS.getD i 0 < 2 := by
  simp only [checkData, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq,
    List.all_eq_true] at h
  exact ⟨h.1.1.1.1.1.1.1, h.1.1.1.1.1.2, h.1.1.1.1.2,
    fun i => h.1.2 i (List.mem_range.mpr i.isLt)⟩

theorem cancelD_le_selected (j : ℕ) (d : WeightData) (x : ℝ)
    (hv : checkData j d = true)
    (ha : Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j d).eval x)
    (hlow : lowMaximum j d (1 - x / 6250) ≤ (selectedLowPolynomial j d).eval x) :
    (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).cancelD ≤
      ((boundsD j d).getD d.choices.cancelD zero).eval x := by
  apply selected_min_le _ _ _ _ _ (boundsD_compare j d x ha hlow)
  have hh := (checkData_ranges hv).2.2.1
  cases h : d.lowPatterns <;>
    simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsD, h] using hh

theorem cancelS_le_selected (j : ℕ) (d : WeightData) (i : Fin 5) (x : ℝ)
    (hv : checkData j d = true) :
    (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).cancelS i ≤
      ((boundsS j d i).getD (d.choices.cancelS.getD i 0) zero).eval x := by
  apply selected_min_le _ _ _ _ _ (boundsS_compare j d i x (checkData_ranges hv).1)
  have hh := ((checkData_ranges hv).2.2.2 i).1
  cases h : d.lowPatterns <;>
    simpa [Spin.Imt.Occupation.Sparse.cancellationBoundsS, h] using hh

theorem fresh_le_selected (a b : ℝ) (p q : RatPoly) (x : ℝ) (i : ℕ)
    (hi : i < 2) (ha : a ≤ p.eval x) (hb : b ≤ q.eval x) :
    min a b ≤ ([p, q].getD i zero).eval x := by
  interval_cases i
  · exact (min_le_left a b).trans ha
  · exact (min_le_right a b).trans hb

theorem freshD_le_selected (j : ℕ) (d : WeightData) (x : ℝ)
    (hv : checkData j d = true)
    (ha : Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j d).eval x) :
    min (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).momentD
      (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).live ≤
      ([arbitrary j d, live j d].getD d.choices.freshD zero).eval x :=
  fresh_le_selected _ _ _ _ _ _ (checkData_ranges hv).2.1 ha (eval_live j d x).ge

theorem freshS_le_selected (j : ℕ) (d : WeightData) (i : Fin 5) (x : ℝ)
    (hv : checkData j d = true) :
    min ((Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).momentS i)
      (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250)).live ≤
      ([moment (levels.getD i 0) j, live j d].getD (d.choices.freshS.getD i 0) zero).eval x :=
  fresh_le_selected _ _ _ _ _ _ ((checkData_ranges hv).2.2.2 i).2
    (eval_moment _ j x).ge (eval_live j d x).ge

end Spin.Structured.SparsePolynomial
