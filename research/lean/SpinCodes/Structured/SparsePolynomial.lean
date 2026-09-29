import SpinCodes.Structured.SparsePolynomialDefs
import SpinCodes.Structured.SparsePowersSound
import SpinCodes.Structured.SparseWitness

/-! Evaluation of the sparse polynomial program at the paper's parameters.
The binomial coefficients are computed by factorial division and proved equal
to `Nat.choose`. The affine witness is the same column as `SparseWitness`.
-/

noncomputable section

namespace Spin.Imt.Occupation

/-- Hypergeometric emitted moment at image weight `w` and input weight `j`.
The guard excludes negative choices of the second intersection count. -/
def hyperMoment (w j : ℕ) (z : ℝ) : ℝ :=
  ((List.range 129).map fun h => if h ≤ j then
    (w.choose h : ℝ) * ((128 - w).choose (j - h) : ℝ) * z ^ (w + j - 2 * h)
    else 0).sum / (128 : ℕ).choose j

end Spin.Imt.Occupation

namespace Spin.Structured.SparsePolynomial

theorem scale_den_pos (a : ℤ) {d : ℕ} {p : RatPoly} (hd : 0 < d)
    (hp : 0 < p.den) : 0 < (scale a d p).den := Nat.mul_pos hd hp

theorem eval_scale (a : ℤ) (d : ℕ) (p : RatPoly) (x : ℝ) :
    (scale a d p).eval x = ((a : ℝ) / d) * p.eval x := by
  rw [scale, RatPoly.eval_mul, RatPoly.eval_constant]

theorem choose128_pos {j : ℕ} (hj : j ≤ 128) : 0 < polyChoose 128 j := by
  rw [polyChoose_eq]
  exact Nat.choose_pos hj

theorem moment_terms_pos (w j h : ℕ) :
    0 < (let c := polyChoose w h * polyChoose (128 - w) (j - h)
      if h ≤ j ∧ c ≠ 0 then scale c 1 (SparsePowers.zPow (w + j - 2 * h)) else zero).den := by
  dsimp
  split
  · exact scale_den_pos _ (by decide) (SparsePowers.zPow_den_pos _)
  · decide

theorem moment_den_pos (w : ℕ) {j : ℕ} (hj : j ≤ 128) : 0 < (moment w j).den := by
  apply scale_den_pos _ (choose128_pos hj)
  apply RatPoly.sum_den_pos
  intro p hp
  obtain ⟨h, _, rfl⟩ := List.mem_map.mp hp
  exact moment_terms_pos w j h

theorem eval_moment (w j : ℕ) (x : ℝ) :
    (moment w j).eval x = Spin.Imt.Occupation.hyperMoment w j (1 - x / 6250) := by
  have hpos : ∀ p ∈ (List.range 129).map (fun h =>
      let c := polyChoose w h * polyChoose (128 - w) (j - h)
      if h ≤ j ∧ c ≠ 0 then scale c 1 (SparsePowers.zPow (w + j - 2 * h)) else zero),
      0 < p.den := by
    intro p hp
    obtain ⟨h, _, rfl⟩ := List.mem_map.mp hp
    exact moment_terms_pos w j h
  rw [moment, eval_scale, RatPoly.eval_sum _ hpos]
  simp only [List.map_map, Function.comp_def, Spin.Imt.Occupation.hyperMoment,
    Int.cast_one, one_div, polyChoose_eq]
  have he : ∀ h : ℕ,
      (if h ≤ j ∧ w.choose h * (128 - w).choose (j - h) ≠ 0 then
        scale ((w.choose h * (128 - w).choose (j - h) : ℕ) : ℤ) 1 (SparsePowers.zPow (w + j - 2 * h))
       else zero).eval x =
      if h ≤ j then (w.choose h : ℝ) * ((128 - w).choose (j - h) : ℝ) *
        (1 - x / 6250) ^ (w + j - 2 * h) else 0 := by
    intro h
    by_cases hh : h ≤ j
    · by_cases hc : w.choose h * (128 - w).choose (j - h) = 0
      · have hc' : (w.choose h : ℝ) * ((128 - w).choose (j - h) : ℝ) = 0 := by
          exact_mod_cast hc
        simp [hh, hc, hc', zero, RatPoly.eval_constant]
      · simp [hh, hc, eval_scale, SparsePowers.zPow_eval, Int.cast_mul,
          Int.cast_natCast, Nat.cast_mul]
    · simp [hh, zero, RatPoly.eval_constant]
  simp_rw [he]
  ring

theorem probability_den_pos (j : ℕ) : 0 < (probability j).den := by
  apply scale_den_pos _ (by decide)
  exact RatPoly.mul_den_pos _ _ (SparsePowers.betaPow_den_pos _)
    (SparsePowers.complementPow_den_pos _)

theorem eval_probability (j : ℕ) (x : ℝ) :
    (probability j).eval x = Spin.Imt.Occupation.probability 128 j (x / 12500) := by
  simp only [probability, eval_scale, RatPoly.eval_mul, SparsePowers.betaPow_eval,
    SparsePowers.complementPow_eval, polyChoose_eq, Int.cast_natCast,
    Nat.cast_one, div_one, Spin.Imt.Occupation.probability]
  ring

theorem eval_witness_Z (x : ℝ) :
    (witness 0).eval x = (Spin.Imt.Occupation.Sparse.witness (x / 10000)).Z := by
  norm_num [witness, one, RatPoly.eval_constant, Spin.Imt.Occupation.Sparse.witness]

theorem eval_witness_D (x : ℝ) :
    (witness 1).eval x = (Spin.Imt.Occupation.Sparse.witness (x / 10000)).D := by
  norm_num [witness, RatPoly.eval, listEval, Spin.Imt.Occupation.Sparse.witness]
  ring

theorem eval_witness_S (i : Fin 5) (x : ℝ) :
    (witness (i + 2)).eval x = (Spin.Imt.Occupation.Sparse.witness (x / 10000)).S i := by
  fin_cases i <;> norm_num [witness, RatPoly.eval, listEval,
    Spin.Imt.Occupation.Sparse.witness, Spin.Imt.Occupation.Sparse.correction] <;> ring

end Spin.Structured.SparsePolynomial



