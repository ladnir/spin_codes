import SpinCodes.Structured.SparsePolynomial
import SpinCodes.Structured.SparseModelData

/-! The numerical occupation matrix with the frozen kernel and fiber data.
The definitions use the actual maxima and minima from the finite-transfer
formulas. Connecting the frozen counts to the concrete binary maps remains
a finite-counting obligation; no such identification is assumed here.
-/

noncomputable section

namespace Spin.Imt.Occupation.Sparse

open Spin.Structured

def live (j : ℕ) (d : SparsePolynomial.WeightData) : ℝ :=
  ((Nat.choose 128 j : ℝ) - d.kernel) / Nat.choose 128 j

def pattern (j : ℕ) (p : List (ℕ × ℕ)) (z : ℝ) : ℝ :=
  (p.map fun (out, c) => (c : ℝ) * z ^ out).sum / Nat.choose 128 j

def maximumMoment (j : ℕ) (z : ℝ) : ℝ :=
  (SparsePolynomial.levels.map fun w => hyperMoment w j z).foldr max 0

def cancellationBoundsD (j : ℕ) (d : SparsePolynomial.WeightData) (z : ℝ) : List ℝ :=
  [maximumMoment j z, live j d,
   (d.cap : ℝ) / Nat.choose 128 j * z ^ SparsePolynomial.minDistance j] ++
  match d.lowPatterns with
  | none => []
  | some ps => [(ps.map fun p => pattern j p z).foldr max 0]

def cancellationBoundsS (j : ℕ) (d : SparsePolynomial.WeightData) (z : ℝ)
    (i : Fin 5) : List ℝ :=
  let w := SparsePolynomial.levels.getD i 0
  let c := SparsePolynomial.counts.getD i 1
  [hyperMoment w j z,
   (min (Nat.choose 128 j - d.kernel) (c * d.cap) : ℕ) /
     ((c : ℝ) * Nat.choose 128 j) * z ^ SparsePolynomial.distance w j] ++
  match d.lowPatterns with
  | none => []
  | some _ => [pattern j (d.lowShells.getD i []) z / c]

def row (j : ℕ) (d : SparsePolynomial.WeightData) (z : ℝ) : RowData 5 where
  live := live j d
  zeroMoment := z ^ j
  momentD := maximumMoment j z
  momentS := fun i => hyperMoment (SparsePolynomial.levels.getD i 0) j z
  cancelD := (cancellationBoundsD j d z).foldr min (maximumMoment j z)
  cancelS := fun i => (cancellationBoundsS j d z i).foldr min
    (hyperMoment (SparsePolynomial.levels.getD i 0) j z)

end Spin.Imt.Occupation.Sparse

namespace Spin.Structured.SparsePolynomial

theorem eval_live (j : ℕ) (d : WeightData) (x : ℝ) :
    (live j d).eval x = Spin.Imt.Occupation.Sparse.live j d := by
  simp only [live, RatPoly.eval_constant, polyChoose_eq, Int.cast_sub, Int.cast_natCast,
    Spin.Imt.Occupation.Sparse.live]

theorem live_den_pos {j : ℕ} (hj : j ≤ 128) (d : WeightData) :
    0 < (live j d).den := choose128_pos hj

/-- The zero-coordinate polynomial is exactly the zero row of the occupation
matrix against the sparse witness, including both input-weight endpoints. -/
theorem eval_action_zero {j : ℕ} (hj : j ≤ 128) (d : WeightData) (x : ℝ) :
    (action j d 0).eval x =
      ((Spin.Imt.Occupation.fixed Spin.Imt.Occupation.Sparse.count 524287
        (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))).applyCol
          (Spin.Imt.Occupation.Sparse.witness (x / 10000))).Z := by
  have hchoose := choose128_pos hj
  have hlive := live_den_pos hj d
  have hz := SparsePowers.zPow_den_pos j
  have hw : 0 < (witness 1).den := by decide
  rw [action, RatPoly.eval_add _ _ (scale_den_pos _ hchoose hz)
    (RatPoly.mul_den_pos _ _ (RatPoly.mul_den_pos _ _ hlive hz) hw)]
  rw [eval_scale, RatPoly.eval_mul, RatPoly.eval_mul, eval_live,
    SparsePowers.zPow_eval, eval_witness_D, Spin.Imt.Occupation.fixed_col_Z]
  simp only [Spin.Imt.Occupation.Sparse.row, Spin.Imt.Occupation.Sparse.witness,
    Spin.Imt.Occupation.Sparse.live, polyChoose_eq, Int.cast_natCast, mul_one]
  have hn : (Nat.choose 128 j : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.ne_of_gt (Nat.choose_pos hj))
  field_simp
  ring

theorem eval_zero_contribution {j : ℕ} (hj : j ≤ 128) (d : WeightData) (x : ℝ) :
    (contribution j d 0).eval x =
      Spin.Imt.Occupation.probability 128 j (x / 12500) *
      ((Spin.Imt.Occupation.fixed Spin.Imt.Occupation.Sparse.count 524287
        (Spin.Imt.Occupation.Sparse.row j d (1 - x / 6250))).applyCol
          (Spin.Imt.Occupation.Sparse.witness (x / 10000))).Z := by
  rw [contribution, RatPoly.eval_mul, eval_probability, eval_action_zero hj]

end Spin.Structured.SparsePolynomial
