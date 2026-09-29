import SpinCodes.Structured.DenseOccupationDefs
import SpinCodes.Structured.PositiveFiniteTilted

noncomputable section
namespace Spin.Structured.DenseOccupation
open Spin.Imt SparsePolynomial

abbrev QCoords.real (v : QCoords) : Coords 5 :=
  ⟨v.Z, v.D, fun i => v.S i⟩

lemma cast_list_sum (xs : List ℚ) : ((xs.sum:ℚ):ℝ) = (xs.map (fun q : ℚ => (q:ℝ))).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih => simp [ih]

lemma cast_fold_max (xs : List ℚ) (b : ℚ) :
    ((xs.foldr max b:ℚ):ℝ) = (xs.map (fun q : ℚ => (q:ℝ))).foldr max (b:ℝ) := by
  induction xs with
  | nil => rfl
  | cons a xs ih => simp [ih, Rat.cast_max]

lemma cast_fold_min (xs : List ℚ) (b : ℚ) :
    ((xs.foldr min b:ℚ):ℝ) = (xs.map (fun q : ℚ => (q:ℝ))).foldr min (b:ℝ) := by
  induction xs with
  | nil => rfl
  | cons a xs ih => simp [ih, Rat.cast_min]

@[simp] theorem cast_live (j : ℕ) (d : WeightData) :
    ((live j d:ℚ):ℝ) = Occupation.Sparse.live j d := by
  simp [live, Occupation.Sparse.live]

@[simp] theorem cast_pattern (j : ℕ) (pat : List (ℕ × ℕ)) (z : ℚ) :
    ((pattern j pat z:ℚ):ℝ) = Occupation.Sparse.pattern j pat z := by
  simp [pattern, Occupation.Sparse.pattern, cast_list_sum, List.map_map, Function.comp_def]

@[simp] theorem cast_hyperMoment (w j : ℕ) (z : ℚ) :
    ((hyperMoment w j z:ℚ):ℝ) = Occupation.hyperMoment w j z := by
  simp [hyperMoment, Occupation.hyperMoment, cast_list_sum, List.map_map, Function.comp_def, apply_ite]

@[simp] theorem cast_maximumMoment (j : ℕ) (z : ℚ) :
    ((maximumMoment j z:ℚ):ℝ) = Occupation.Sparse.maximumMoment j z := by
  simp [maximumMoment, Occupation.Sparse.maximumMoment, cast_fold_max, List.map_map,
    Function.comp_def]

@[simp] theorem cast_boundsD (j : ℕ) (d : WeightData) (z : ℚ) :
    (cancellationBoundsD j d z).map (fun q : ℚ => (q:ℝ)) =
      Occupation.Sparse.cancellationBoundsD j d z := by
  unfold cancellationBoundsD Occupation.Sparse.cancellationBoundsD
  cases d.lowPatterns <;> simp [cast_fold_max, List.map_map, Function.comp_def]

@[simp] theorem cast_boundsS (j : ℕ) (d : WeightData) (z : ℚ) (i : Fin 5) :
    (cancellationBoundsS j d z i).map (fun q : ℚ => (q:ℝ)) =
      Occupation.Sparse.cancellationBoundsS j d z i := by
  unfold cancellationBoundsS Occupation.Sparse.cancellationBoundsS
  cases d.lowPatterns <;> simp

@[simp] theorem cast_liveAverage (v : QCoords) :
    ((liveAverage v:ℚ):ℝ) = Occupation.liveAverage Occupation.Sparse.count 524287 v.real := by
  simp only [liveAverage, Occupation.liveAverage, Rat.cast_div, Rat.cast_sum, Rat.cast_mul, Rat.cast_natCast, Rat.cast_ofNat]
  congr 1

@[simp] theorem cast_fixedColumn (j : ℕ) (d : WeightData) (z : ℚ) (v : QCoords) :
    (fixedColumn j d z v).real =
      (Occupation.fixed Occupation.Sparse.count 524287 (Occupation.Sparse.row j d z)).applyCol v.real := by
  apply Coords.ext
  · simp [QCoords.real, fixedColumn, Occupation.fixed_col_Z, Occupation.Sparse.row]
  · simp [QCoords.real, fixedColumn, Occupation.fixed_col_D, Occupation.Sparse.row,
      cast_fold_min, ← cast_liveAverage]
  · funext i
    have he : Occupation.Sparse.live j d = 0 ↔ live j d = 0 := by
      rw [← cast_live]; exact Rat.cast_eq_zero
    simp [QCoords.real, fixedColumn, Occupation.fixed_col_S, Occupation.Sparse.row,
      cast_fold_min, ← cast_liveAverage, apply_ite, he]

@[simp] theorem cast_probability (j : ℕ) (q : ℚ) :
    ((probability j q:ℚ):ℝ) = Occupation.probability 128 j q := by
  simp [probability, Occupation.probability]

@[simp] theorem cast_column (q z : ℚ) (v : QCoords) :
    (column q z v).real = (Occupation.Sparse.numericalMatrix q z).applyCol v.real := by
  rw [Occupation.Sparse.numericalMatrix, Occupation.matrix_col]
  apply Coords.ext
  · simp only [QCoords.real, column, Occupation.Sparse.count, Coords.sum, Coords.smul,
      Rat.cast_sum, Rat.cast_mul, cast_probability]
    apply Finset.sum_congr rfl
    intro j _
    exact congrArg (fun c : Coords 5 => Occupation.probability 128 j q*c.Z)
      (cast_fixedColumn j (Data.weight j) z v)
  · simp only [QCoords.real, column, Coords.sum, Coords.smul,
      Rat.cast_sum, Rat.cast_mul, cast_probability]
    apply Finset.sum_congr rfl
    intro j _
    exact congrArg (fun c : Coords 5 => Occupation.probability 128 j q*c.D)
      (cast_fixedColumn j (Data.weight j) z v)
  · funext i
    simp only [QCoords.real, column, Coords.sum, Coords.smul,
      Rat.cast_sum, Rat.cast_mul, cast_probability]
    apply Finset.sum_congr rfl
    intro j _
    exact congrArg (fun c : Coords 5 => Occupation.probability 128 j q*c.S i)
      (cast_fixedColumn j (Data.weight j) z v)

end Spin.Structured.DenseOccupation
