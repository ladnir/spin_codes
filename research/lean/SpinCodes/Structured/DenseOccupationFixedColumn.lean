import SpinCodes.Structured.DenseOccupationFixedMoments
import SpinCodes.Structured.DenseOccupationFixedColumnDefs

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric Spin.Imt SparsePolynomial

lemma fold_min_mem₂ {xs : List Fix} {ys : List ℝ} {b : Fix} {v : ℝ}
    (h : List.Forall₂ Fix.Mem xs ys) (hb : Fix.Mem b v) :
    Fix.Mem (xs.foldr iminFix b) (ys.foldr min v) := by
  induction h with
  | nil => exact hb
  | cons h _ ih => exact min_mem h ih

lemma count_pos (i : Fin 5) : 0 < counts.getD i 1 := by
  fin_cases i <;> decide

lemma boundsD_mem {j : Nat} (hj : j≤128) (d : WeightData) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) :
    List.Forall₂ Fix.Mem (boundsD j d zp) (Occupation.Sparse.cancellationBoundsD j d z) := by
  have hc : 0 < (polyChoose 128 j:Int) := by
    exact_mod_cast (show 0 < polyChoose 128 j from by rw [polyChoose_eq]; exact Nat.choose_pos hj)
  have hcap := Fix.mul_mem (Fix.ofFrac_mem (p := (d.cap:Int)) hc) (hz (minDistance j))
  simp only [Int.cast_natCast, polyChoose_eq] at hcap
  have hbase := List.Forall₂.cons (maximumMoment_mem hj zp z hz)
    (List.Forall₂.cons (live_mem hj d) (List.Forall₂.cons hcap List.Forall₂.nil))
  unfold boundsD Occupation.Sparse.cancellationBoundsD
  cases d.lowPatterns with
  | none => simpa [polyChoose_eq] using hbase
  | some ps =>
    simpa [polyChoose_eq, Int.cast_natCast] using List.rel_append hbase (List.Forall₂.cons
      (fold_max_mem ps _ _ zero_mem (fun p _ => pattern_mem hj p zp z hz)) List.Forall₂.nil)

lemma boundsS_mem {j : Nat} (hj : j≤128) (d : WeightData) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) (i : Fin 5) :
    List.Forall₂ Fix.Mem (boundsS j d zp i) (Occupation.Sparse.cancellationBoundsS j d z i) := by
  have hc : 0 < ((counts.getD i 1 * polyChoose 128 j:Nat):Int) := by
    exact_mod_cast Nat.mul_pos (count_pos i) (show 0 < polyChoose 128 j from by
      rw [polyChoose_eq]; exact Nat.choose_pos hj)
  have hcap := Fix.mul_mem (Fix.ofFrac_mem
    (p := ((min (polyChoose 128 j-d.kernel) (counts.getD i 1*d.cap):Nat):Int)) hc)
      (hz (distance (levels.getD i 0) j))
  simp only [Int.cast_natCast, Nat.cast_mul, polyChoose_eq] at hcap
  have hbase := List.Forall₂.cons (moment_mem hj (levels.getD i 0) zp z hz)
    (List.Forall₂.cons hcap List.Forall₂.nil)
  unfold boundsS Occupation.Sparse.cancellationBoundsS
  dsimp only
  cases d.lowPatterns with
  | none => simpa [polyChoose_eq] using hbase
  | some ps =>
    have hpat : Fix.Mem (Fix.divInt (pattern j (d.lowShells.getD i []) zp) (counts.getD i 1))
        (Occupation.Sparse.pattern j (d.lowShells.getD i []) z / (counts.getD i 1:ℝ)) := by
      simpa only [Int.cast_natCast] using
        (Fix.divInt_mem (k := (counts.getD i 1:Int))
          (pattern_mem hj (d.lowShells.getD i []) zp z hz) (by exact_mod_cast count_pos i))
    simpa [polyChoose_eq, Int.cast_natCast] using List.rel_append hbase
      (List.Forall₂.cons hpat List.Forall₂.nil)

def FCoords.Mem (v : FCoords) (w : Coords 5) : Prop :=
  Fix.Mem v.Z w.Z ∧ Fix.Mem v.D w.D ∧ ∀ i, Fix.Mem (v.S i) (w.S i)

lemma liveAverage_mem {v : FCoords} {w : Coords 5} (hv : v.Mem w) :
    Fix.Mem (liveAverage v) (Occupation.liveAverage Occupation.Sparse.count 524287 w) := by
  unfold liveAverage Occupation.liveAverage
  apply Fix.divInt_mem _ (by norm_num)
  have h := sum_mem (List.finRange 5)
    (fun i => Fix.mul (Fix.ofInt (counts.getD i 1)) (v.S i))
    (fun i => (counts.getD i 1:ℝ)*w.S i)
    (fun i _ => Fix.mul_mem (Fix.ofInt_mem _) (hv.2.2 i))
  have hc : ∀ i : Fin 5, (counts.getD i 1:ℝ) = Occupation.Sparse.count i := by
    intro i; fin_cases i <;> norm_num [counts, Occupation.Sparse.count]
  simpa only [← List.ofFn_eq_map, List.sum_ofFn, hc] using h

lemma fixedColumn_mem {j : Nat} (hj : j≤128) (d : WeightData) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) {v : FCoords} {w : Coords 5} (hv : v.Mem w) :
    (fixedColumn j d zp v).Mem
      ((Occupation.fixed Occupation.Sparse.count 524287 (Occupation.Sparse.row j d z)).applyCol w) := by
  have hl := live_mem hj d
  have hm := maximumMoment_mem hj zp z hz
  have ha := liveAverage_mem hv
  have hc : (polyChoose 128 j=d.kernel) ↔ Occupation.Sparse.live j d = 0 := by
    unfold Occupation.Sparse.live
    rw [div_eq_zero_iff]
    have hp : (Nat.choose 128 j:ℝ) ≠ 0 := by exact_mod_cast (Nat.choose_pos hj).ne'
    simp only [hp, or_false, sub_eq_zero]
    rw [polyChoose_eq]
    exact_mod_cast Iff.rfl
  refine ⟨?_,?_,fun i => ?_⟩
  · rw [Occupation.fixed_col_Z]
    exact Fix.add_mem (Fix.mul_mem (Fix.mul_mem (Fix.sub_mem (by simpa using Fix.ofInt_mem 1) hl) (hz j)) hv.1)
      (Fix.mul_mem (Fix.mul_mem hl (hz j)) hv.2.1)
  · rw [Occupation.fixed_col_D]
    norm_num only [Occupation.Sparse.row, show (2:ℝ)*524287=1048574 by norm_num]
    exact Fix.add_mem
      (Fix.mul_mem (Fix.add_mem (Fix.divInt_mem (fold_min_mem₂ (boundsD_mem hj d zp z hz) hm) (by norm_num))
        (Fix.divInt_mem (k := 1048574) (min_mem hm hl) (by norm_num))) hv.1)
      (Fix.divInt_mem (Fix.mul_mem hm (Fix.add_mem hv.2.1 ha)) (by norm_num))
  · rw [Occupation.fixed_col_S]
    norm_num only [Occupation.Sparse.row, show (2:ℝ)*524287=1048574 by norm_num]
    have hmi := moment_mem hj (levels.getD i 0) zp z hz
    have hbranch : Fix.Mem (if polyChoose 128 j=d.kernel then v.S i else v.D)
        (if Occupation.Sparse.live j d=0 then w.S i else w.D) := by
      by_cases he : polyChoose 128 j=d.kernel
      · simpa [he, hc.mp he] using hv.2.2 i
      · simpa [he, not_iff_not.mpr hc |>.mp he] using hv.2.1
    exact Fix.add_mem
      (Fix.mul_mem (Fix.add_mem (Fix.divInt_mem (fold_min_mem₂ (boundsS_mem hj d zp z hz i) hmi) (by norm_num))
        (Fix.divInt_mem (k := 1048574) (min_mem hmi hl) (by norm_num))) hv.1)
      (Fix.divInt_mem (Fix.mul_mem hmi (Fix.add_mem ha hbranch)) (by norm_num))

end Spin.Structured.DenseOccupationFixed
