import SpinCodes.Structured.SparseModel
import SpinCodes.Structured.PolySignIdentity

/-! Analytic facts used to turn the checked finite branch choices into
upper bounds for the numerical occupation matrix. -/

noncomputable section

namespace Spin.Imt.Occupation

theorem hyperMoment_nonneg (w j : ℕ) {z : ℝ} (hz : 0 ≤ z) :
    0 ≤ hyperMoment w j z := by
  apply div_nonneg _ (Nat.cast_nonneg _)
  apply List.sum_nonneg
  intro a ha
  obtain ⟨h, _, rfl⟩ := List.mem_map.mp ha
  split <;> positivity

namespace Sparse

theorem pattern_nonneg (j : ℕ) (p : List (ℕ × ℕ)) {z : ℝ} (hz : 0 ≤ z) :
    0 ≤ pattern j p z := by
  apply div_nonneg _ (Nat.cast_nonneg _)
  apply List.sum_nonneg
  intro a ha
  obtain ⟨⟨out, c⟩, _, rfl⟩ := List.mem_map.mp ha
  positivity

end Sparse
end Spin.Imt.Occupation

namespace Spin.Structured.SparsePolynomial

def lowMaximum (j : ℕ) (d : WeightData) (z : ℝ) : ℝ :=
  match d.lowPatterns with
  | none => 0
  | some ps => (ps.map fun p => Spin.Imt.Occupation.Sparse.pattern j p z).foldr max 0

def selectedLowPolynomial (j : ℕ) (d : WeightData) : RatPoly :=
  match d.lowPatterns with
  | none => zero
  | some ps => (ps.map (pattern j)).getD d.choices.lowMax zero

theorem eval_pattern (j : ℕ) (p : List (ℕ × ℕ)) (x : ℝ) :
    (pattern j p).eval x = Spin.Imt.Occupation.Sparse.pattern j p (1 - x / 6250) := by
  have hpos : ∀ q ∈ p.map (fun (out, c) => scale c 1 (SparsePowers.zPow out)),
      0 < q.den := by
    intro q hq
    obtain ⟨⟨out, c⟩, _, rfl⟩ := List.mem_map.mp hq
    exact scale_den_pos _ (by decide) (SparsePowers.zPow_den_pos out)
  rw [pattern, eval_scale, RatPoly.eval_sum _ hpos]
  simp only [List.map_map, Function.comp_def, eval_scale, SparsePowers.zPow_eval,
    Int.cast_one, one_div, Nat.cast_one, div_one, Int.cast_natCast, polyChoose_eq,
    Spin.Imt.Occupation.Sparse.pattern]
  ring

theorem moment_eval_nonneg (w j : ℕ) {x : ℝ} (hx : x ≤ 1) :
    0 ≤ (moment w j).eval x := by
  rw [eval_moment]
  exact Spin.Imt.Occupation.hyperMoment_nonneg w j (by linarith)

/-- A selected list entry bounds the minimum whenever its index is valid. -/
theorem foldr_min_le_getD (xs : List ℝ) (base fallback : ℝ) (i : ℕ)
    (hi : i < xs.length) : xs.foldr min base ≤ xs.getD i fallback := by
  induction xs generalizing i with
  | nil => simp at hi
  | cons a xs ih =>
    cases i with
    | zero => simpa using (min_le_left a (xs.foldr min base))
    | succ i =>
      have ht : i < xs.length := by simpa using hi
      exact (min_le_right a (xs.foldr min base)).trans (by simpa using ih i ht)

theorem foldr_max_le (xs : List ℝ) {base upper : ℝ} (hb : base ≤ upper)
    (h : ∀ a ∈ xs, a ≤ upper) : xs.foldr max base ≤ upper := by
  induction xs with
  | nil => exact hb
  | cons a xs ih =>
    exact max_le (h a (List.mem_cons_self ..))
      (ih fun b hb => h b (List.mem_cons_of_mem _ hb))

theorem maximumMoment_le_polynomial (j : ℕ) (x : ℝ) (hx : x ≤ 1) (u : RatPoly)
    (h48 : (moment 48 j).eval x ≤ u.eval x)
    (h56 : (moment 56 j).eval x ≤ u.eval x)
    (h64 : (moment 64 j).eval x ≤ u.eval x)
    (h72 : (moment 72 j).eval x ≤ u.eval x)
    (h80 : (moment 80 j).eval x ≤ u.eval x) :
    Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤ u.eval x := by
  have h0 : 0 ≤ u.eval x := (moment_eval_nonneg 48 j hx).trans h48
  have hh := max_le h48 (max_le h56 (max_le h64 (max_le h72 (max_le h80 h0))))
  simpa only [Spin.Imt.Occupation.Sparse.maximumMoment, levels,
    List.map_cons, List.map_nil, List.foldr_cons, List.foldr_nil, eval_moment] using hh

theorem maximumPattern_le_polynomial (j : ℕ) (ps : List (List (ℕ × ℕ)))
    (x : ℝ) (u : RatPoly) (h0 : 0 ≤ u.eval x)
    (h : ∀ p ∈ ps, (pattern j p).eval x ≤ u.eval x) :
    (ps.map fun p => Spin.Imt.Occupation.Sparse.pattern j p (1 - x / 6250)).foldr max 0
      ≤ u.eval x := by
  apply foldr_max_le _ h0
  intro a ha
  obtain ⟨p, hp, rfl⟩ := List.mem_map.mp ha
  simpa only [eval_pattern] using h p hp

end Spin.Structured.SparsePolynomial
