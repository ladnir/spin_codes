import SpinCodes.Structured.SparseCancellation

/-! The six live fixed-weight rows are dominated by the checked polynomial
actions. The zero row is the exact identity in `SparseModel`. -/

noncomputable section
namespace Spin.Structured.SparsePolynomial

theorem witness_den_pos (i : ℕ) : 0 < (witness i).den := by
  unfold witness
  split <;> decide

theorem eval_envelope (a c f u v : RatPoly) (x : ℝ)
    (ha : 0 < a.den) (hc : 0 < c.den) (hf : 0 < f.den)
    (hu : 0 < u.den) (hv : 0 < v.den) :
    (RatPoly.add (RatPoly.add (scale 1 2 c) (scale 1 1048574 f))
      (scale 1 2 (RatPoly.mul a (RatPoly.add u v)))).eval x =
      c.eval x / 2 + f.eval x / 1048574 + a.eval x * (u.eval x + v.eval x) / 2 := by
  rw [RatPoly.eval_add _ _
    (RatPoly.add_den_pos _ _ (scale_den_pos _ (by decide) hc)
      (scale_den_pos _ (by decide) hf))
    (scale_den_pos _ (by decide) (RatPoly.mul_den_pos _ _ ha
      (RatPoly.add_den_pos _ _ hu hv))),
    RatPoly.eval_add _ _ (scale_den_pos _ (by decide) hc) (scale_den_pos _ (by decide) hf),
    eval_scale, eval_scale, eval_scale, RatPoly.eval_mul, RatPoly.eval_add _ _ hu hv]
  norm_num only [Int.cast_one]
  ring

theorem pair_den_pos (p q : RatPoly) (hp : 0 < p.den) (hq : 0 < q.den) (i : ℕ) :
    0 < ([p, q].getD i zero).den := by
  apply getD_den_pos _ _ _ (by decide)
  intro a ha
  simp only [List.mem_cons, List.not_mem_nil, or_false] at ha
  rcases ha with rfl | rfl <;> assumption

theorem eval_action_D {j : ℕ} (hj : j ≤ 128) (d : WeightData) (x : ℝ) :
    (action j d 1).eval x =
      ((boundsD j d).getD d.choices.cancelD zero).eval x / 2 +
      ([arbitrary j d, live j d].getD d.choices.freshD zero).eval x / 1048574 +
      (arbitrary j d).eval x *
        ((Spin.Imt.Occupation.Sparse.witness (x / 10000)).D + 1 / 1024) / 2 := by
  rw [action, eval_envelope _ _ _ _ _ _ (arbitrary_den_pos hj d)
    (getD_den_pos _ _ _ (by decide) (boundsD_den_pos hj d))
    (pair_den_pos _ _ (arbitrary_den_pos hj d) (live_den_pos hj d) _)
    (witness_den_pos 1) (by decide), eval_witness_D, RatPoly.eval_constant]
  norm_num

theorem live_eq_zero_iff {j : ℕ} (hj : j ≤ 128) (d : WeightData) :
    Spin.Imt.Occupation.Sparse.live j d = 0 ↔ polyChoose 128 j = d.kernel := by
  have hn : (Nat.choose 128 j : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.ne_of_gt (Nat.choose_pos hj))
  simp only [Spin.Imt.Occupation.Sparse.live, div_eq_zero_iff, hn, or_false,
    sub_eq_zero, Nat.cast_inj, polyChoose_eq]

theorem eval_action_S {j : ℕ} (hj : j ≤ 128) (d : WeightData) (i : Fin 5) (x : ℝ) :
    (action j d (i + 2)).eval x =
      ((boundsS j d i).getD (d.choices.cancelS.getD i 0) zero).eval x / 2 +
      ([moment (levels.getD i 0) j, live j d].getD (d.choices.freshS.getD i 0) zero).eval x /
        1048574 +
      Spin.Imt.Occupation.hyperMoment (levels.getD i 0) j (1 - x / 6250) *
        (1 / 1024 + if Spin.Imt.Occupation.Sparse.live j d = 0
          then (Spin.Imt.Occupation.Sparse.witness (x / 10000)).S i
          else (Spin.Imt.Occupation.Sparse.witness (x / 10000)).D) / 2 := by
  have hlazy : 0 < (if polyChoose 128 j = d.kernel then witness (i + 2) else witness 1).den := by
    split <;> apply witness_den_pos
  rw [action, eval_envelope _ _ _ _ _ _ (moment_den_pos _ hj)
    (getD_den_pos _ _ _ (by decide) (boundsS_den_pos hj d i))
    (pair_den_pos _ _ (moment_den_pos _ hj) (live_den_pos hj d) _)
    (by decide) hlazy, eval_moment, RatPoly.eval_constant]
  simp only [live_eq_zero_iff hj]
  split <;> simp_all only [if_pos, if_neg, eval_witness_S, eval_witness_D,
    Int.cast_one, Nat.cast_ofNat]

theorem fixed_col_D_le_action {j : ℕ} (hj : j ≤ 128) (d : WeightData) (x : ℝ)
    (hx : 0 ≤ x) (hv : checkData j d = true)
    (ha : Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j d).eval x)
    (hlow : lowMaximum j d (1 - x / 6250) ≤ (selectedLowPolynomial j d).eval x) :
    ((Spin.Imt.Occupation.fixed Spin.Imt.Occupation.Sparse.count 524287
      (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))).applyCol
        (Spin.Imt.Occupation.Sparse.witness (x / 10000))).D ≤ (action j d 1).eval x := by
  have h := Spin.Imt.Occupation.fixed_col_D_le Spin.Imt.Occupation.Sparse.count 524287
    (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))
    (Spin.Imt.Occupation.Sparse.witness (x / 10000)) (by norm_num) (by norm_num [Spin.Imt.Occupation.Sparse.witness])
    (Spin.Imt.Occupation.Sparse.live_pair_nonneg (by positivity)) ha
    (cancelD_le_selected j d x hv ha hlow) (freshD_le_selected j d x hv ha)
  rw [eval_action_D hj, Spin.Imt.Occupation.Sparse.liveAverage_witness] at *
  simpa only [Spin.Imt.Occupation.Sparse.witness, mul_one, show (2 : ℝ) * 524287 = 1048574 by norm_num] using h

theorem fixed_col_S_le_action {j : ℕ} (hj : j ≤ 128) (d : WeightData) (i : Fin 5)
    (x : ℝ) (hv : checkData j d = true) :
    ((Spin.Imt.Occupation.fixed Spin.Imt.Occupation.Sparse.count 524287
      (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))).applyCol
        (Spin.Imt.Occupation.Sparse.witness (x / 10000))).S i ≤ (action j d (i + 2)).eval x := by
  have h := Spin.Imt.Occupation.fixed_col_S_le Spin.Imt.Occupation.Sparse.count 524287
    (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))
    (Spin.Imt.Occupation.Sparse.witness (x / 10000)) i (by norm_num)
    (by norm_num [Spin.Imt.Occupation.Sparse.witness])
    (cancelS_le_selected j d i x hv) (freshS_le_selected j d i x hv)
  rw [eval_action_S hj]
  rw [Spin.Imt.Occupation.Sparse.liveAverage_witness] at h
  by_cases hl : Spin.Imt.Occupation.Sparse.live j d = 0 <;>
    simpa only [Spin.Imt.Occupation.Sparse.row, Spin.Imt.Occupation.Sparse.witness,
      hl, if_true, if_false, mul_one, show (2 : ℝ) * 524287 = 1048574 by norm_num] using h

end Spin.Structured.SparsePolynomial
