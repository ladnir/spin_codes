import SpinCodes.Turbo.SparseListDefs22
import SpinCodes.Turbo.SparseFixedSound22

/-! Soundness of the list-structured checks. Semantic access uses `getD`
only in statements; the kernel never evaluates it at endpoint size. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22

def sparseLMat22 (M : List (List Int)) (i j : ℕ) : Int := (M.getD i []).getD j 0

def sparseLPoly22 (P : List (List (List Int))) (d i j : ℕ) : Int :=
  sparseLMat22 (P.getD d []) i j

namespace SparseListSound22

theorem gd {α : Type*} (l : List α) (d : α) {i : ℕ} (h : i < l.length) :
    l.getD i d = l[i] := (List.getElem_eq_getD d).symm

def MatShape (M : List (List Int)) : Prop := M.length = 11 ∧ ∀ r ∈ M, r.length = 11

theorem matShape_of_check (M : List (List Int)) (h : sparseMatShapeL22 M = true) :
    MatShape M := by
  simp only [sparseMatShapeL22, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true,
    sparseLen11L22] at h
  exact ⟨h.1, h.2⟩

theorem polyShape_of_check (P : List (List (List Int))) (h : sparsePolyShapeL22 P = true) :
    P.length = 3 ∧ ∀ d, d < 3 → MatShape (P.getD d []) := by
  simp only [sparsePolyShapeL22, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  refine ⟨h.1, fun d hd => ?_⟩
  apply matShape_of_check
  rw [gd _ _ (by omega)]
  exact h.2 _ (List.getElem_mem _)

theorem zipWith_getD {f : Int → Int → Int} {a b : List Int} {n j : ℕ}
    (ha : a.length = n) (hb : b.length = n) (hj : j < n) :
    (List.zipWith f a b).getD j 0 = f (a.getD j 0) (b.getD j 0) := by
  have hz : (List.zipWith f a b).length = n := by simp [ha, hb]
  rw [gd _ _ (by omega), gd _ _ (by omega),
    gd _ _ (by omega), List.getElem_zipWith]

theorem map_getD {g : Int → Int} {a : List Int} {j : ℕ} (hj : j < a.length) :
    (a.map g).getD j 0 = g (a.getD j 0) := by
  rw [gd _ _ (by simpa using hj), gd _ _ hj,
    List.getElem_map]

theorem fold_spec : ∀ (ps : List (Int × List Int)) (acc : List Int),
    acc.length = 11 → (∀ p ∈ ps, p.2.length = 11) →
    (ps.foldl (fun acc p => List.zipWith (· + ·) acc (p.2.map (p.1 * ·))) acc).length = 11 ∧
      ∀ j, j < 11 →
        (ps.foldl (fun acc p => List.zipWith (· + ·) acc (p.2.map (p.1 * ·))) acc).getD j 0 =
          acc.getD j 0 + (ps.map fun p => p.1 * p.2.getD j 0).sum
  | [], acc, hacc, _ => ⟨hacc, fun j _ => by simp⟩
  | p :: ps, acc, hacc, hps => by
    have hp : p.2.length = 11 := hps p List.mem_cons_self
    have hlen : (List.zipWith (· + ·) acc (p.2.map (p.1 * ·))).length = 11 := by
      simp [hacc, hp]
    obtain ⟨h1, h2⟩ := fold_spec ps _ hlen (fun q hq => hps q (List.mem_cons_of_mem _ hq))
    refine ⟨h1, fun j hj => ?_⟩
    rw [List.foldl_cons, h2 j hj, zipWith_getD hacc (by simp [hp]) hj,
      map_getD (by omega), List.map_cons, List.sum_cons]
    ring

theorem zip_sum_gen (g : Int → List Int → Int) :
    ∀ (n : ℕ) (x : List Int) (M : List (List Int)), x.length = n → M.length = n →
      ((x.zip M).map (fun p => g p.1 p.2)).sum =
        ∑ i : Fin n, g (x.getD i.val 0) (M.getD i.val [])
  | 0, x, M, hx, hM => by
      rw [List.eq_nil_of_length_eq_zero hx, List.eq_nil_of_length_eq_zero hM]
      simp
  | n + 1, a :: x, r :: M, hx, hM => by
      simp only [List.zip_cons_cons, List.map_cons, List.sum_cons, Fin.sum_univ_succ,
        Fin.val_zero, Fin.val_succ, List.getD_cons_zero, List.getD_cons_succ]
      rw [zip_sum_gen g n x M (by simpa using hx) (by simpa using hM)]
  | n + 1, [], _, hx, _ => by simp at hx
  | n + 1, _ :: _, [], _, hM => by simp at hM

theorem vecMat_spec (x : List Int) (M : List (List Int)) (hx : x.length = 11)
    (hM : MatShape M) :
    (sparseVecMatL22 x M).length = 11 ∧ ∀ j, j < 11 →
      (sparseVecMatL22 x M).getD j 0 =
        ∑ i : Fin 11, x.getD i.val 0 * sparseLMat22 M i.val j := by
  have hps : ∀ p ∈ x.zip M, p.2.length = 11 := fun p hp => hM.2 _ (List.of_mem_zip hp).2
  obtain ⟨h1, h2⟩ := fold_spec (x.zip M) sparseZero11L22 rfl hps
  refine ⟨h1, fun j hj => ?_⟩
  unfold sparseVecMatL22
  rw [h2 j hj]
  have h0 : sparseZero11L22.getD j 0 = 0 := by
    unfold sparseZero11L22; interval_cases j <;> rfl
  rw [h0, zero_add]
  exact zip_sum_gen (fun a r => a * r.getD j 0) 11 x M hx hM.1

theorem all_zipWith {α β : Type*} {g : α → β → Bool} {a : List α} {b : List β}
    (h : (List.zipWith g a b).all id = true) (i : ℕ) (ha : i < a.length) (hb : i < b.length) :
    g a[i] b[i] = true := by
  have hi : i < (List.zipWith g a b).length := by simp [ha, hb]
  have hm := List.all_eq_true.mp h _ (List.getElem_mem hi)
  simpa [List.getElem_zipWith] using hm

theorem leScaled_sound (s : Int) (lhs rhs : List Int) (hl : lhs.length = 11)
    (hr : rhs.length = 11) (h : sparseLeScaledL22 s lhs rhs = true) (j : ℕ) (hj : j < 11) :
    lhs.getD j 0 ≤ s * rhs.getD j 0 := by
  have := all_zipWith h j (by omega) (by omega)
  rw [gd _ _ (by omega), gd _ _ (by omega)]
  simpa using this

theorem matMul_spec (A B : List (List Int)) (hA : MatShape A) (hB : MatShape B) :
    MatShape (sparseMatMulL22 A B) ∧ ∀ i j, i < 11 → j < 11 →
      sparseLMat22 (sparseMatMulL22 A B) i j =
        ∑ k : Fin 11, sparseLMat22 A i k.val * sparseLMat22 B k.val j := by
  refine ⟨⟨by simp [sparseMatMulL22, hA.1], fun r hr => ?_⟩, fun i j hi hj => ?_⟩
  · obtain ⟨row, hrow, rfl⟩ := List.mem_map.mp hr
    exact (vecMat_spec row B (hA.2 row hrow) hB).1
  · have hiA : i < A.length := by rw [hA.1]; exact hi
    have hrow := (vecMat_spec _ B (hA.2 _ (List.getElem_mem hiA)) hB).2 j hj
    show ((A.map fun row => sparseVecMatL22 row B).getD i []).getD j 0 = _
    rw [gd (A.map fun row => sparseVecMatL22 row B) [] (by simpa using hiA),
      List.getElem_map, hrow]
    unfold sparseLMat22
    rw [gd A [] hiA]

theorem matAdd_spec (A B : List (List Int)) (hA : MatShape A) (hB : MatShape B) :
    MatShape (sparseMatAddL22 A B) ∧ ∀ i j, i < 11 → j < 11 →
      sparseLMat22 (sparseMatAddL22 A B) i j = sparseLMat22 A i j + sparseLMat22 B i j := by
  have hlen : (sparseMatAddL22 A B).length = 11 := by simp [sparseMatAddL22, hA.1, hB.1]
  refine ⟨⟨hlen, fun r hr => ?_⟩, fun i j hi hj => ?_⟩
  · obtain ⟨k, hk, rfl⟩ := List.getElem_of_mem hr
    have hkA : k < A.length := by rw [hA.1]; omega
    have hkB : k < B.length := by rw [hB.1]; omega
    simp only [sparseMatAddL22, List.getElem_zipWith, sparseVecAddL22, List.length_zipWith,
      hA.2 _ (List.getElem_mem hkA), hB.2 _ (List.getElem_mem hkB), min_self]
  · have hiA : i < A.length := by rw [hA.1]; exact hi
    have hiB : i < B.length := by rw [hB.1]; exact hi
    show ((List.zipWith sparseVecAddL22 A B).getD i []).getD j 0 =
      (A.getD i []).getD j 0 + (B.getD i []).getD j 0
    rw [gd (List.zipWith sparseVecAddL22 A B) [] (by simp [hA.1, hB.1, hi]),
      List.getElem_zipWith, gd A [] hiA, gd B [] hiB]
    exact zipWith_getD (hA.2 _ (List.getElem_mem hiA)) (hB.2 _ (List.getElem_mem hiB)) hj

theorem matLeScaled_sound (s : Int) (L R : List (List Int)) (hL : MatShape L)
    (hR : MatShape R) (h : sparseMatLeScaledL22 s L R = true) (i j : ℕ) (hi : i < 11)
    (hj : j < 11) : sparseLMat22 L i j ≤ s * sparseLMat22 R i j := by
  have hiL : i < L.length := by rw [hL.1]; exact hi
  have hiR : i < R.length := by rw [hR.1]; exact hi
  have hrow := all_zipWith h i hiL hiR
  unfold sparseLMat22
  rw [gd _ _ hiL, gd _ _ hiR]
  exact leScaled_sound s _ _ (hL.2 _ (List.getElem_mem hiL)) (hR.2 _ (List.getElem_mem hiR))
    hrow j hj

/-! ## Polynomial product -/

theorem polyMul_sound (scale : ℤ) (hs : 0 < scale) (a b u : List (List (List Int)))
    (ha : ∀ d, d < 3 → MatShape (a.getD d [])) (hb : ∀ d, d < 3 → MatShape (b.getD d []))
    (hcheck : sparsePolyMulCheckL22 scale a b u = true)
    (F G : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hF0 : ∀ d i j, 0 ≤ F d i j) (hG0 : ∀ d i j, 0 ≤ G d i j)
    (hF : ∀ d, d ≤ 2 → ∀ i j : Fin 11, F d i j ≤ (sparseLPoly22 a d i j : ℝ) / scale)
    (hG : ∀ d, d ≤ 2 → ∀ i j : Fin 11, G d i j ≤ (sparseLPoly22 b d i j : ℝ) / scale) :
    (∀ d, d < 3 → MatShape (u.getD d [])) ∧
      ∀ d, d ≤ 2 → ∀ i j : Fin 11,
        (∑ p ∈ antidiagonal d, F p.1 * G p.2) i j ≤ (sparseLPoly22 u d i j : ℝ) / scale := by
  simp only [sparsePolyMulCheckL22, Bool.and_eq_true] at hcheck
  obtain ⟨⟨⟨hu, h0⟩, h1⟩, h2⟩ := hcheck
  have huS := (polyShape_of_check u hu).2
  refine ⟨huS, ?_⟩
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hpairInt : ∀ r t, r ≤ 2 → t ≤ 2 → ∀ i j : Fin 11, (F r * G t) i j ≤
      (sparseLMat22 (sparseMatMulL22 (a.getD r []) (b.getD t [])) i j : ℝ) /
        ((scale : ℝ) * scale) := by
    intro r t hr ht i j
    rw [(matMul_spec _ _ (ha r (by omega)) (hb t (by omega))).2 i j i.isLt j.isLt,
      Matrix.mul_apply, Int.cast_sum]
    simp only [Int.cast_mul]
    exact dot_le _ hsR _ _ _ _ (fun k => hF0 r i k) (fun k => hG0 t k j)
      (fun k => hF r hr i k) (fun k => hG t ht k j)
  have mS : ∀ r t, r ≤ 2 → t ≤ 2 → MatShape (sparseMatMulL22 (a.getD r []) (b.getD t [])) :=
    fun r t hr ht => (matMul_spec _ _ (ha r (by omega)) (hb t (by omega))).1
  intro d hd i j
  interval_cases d
  · have hc := matLeScaled_sound scale _ _ (mS 0 0 (by norm_num) (by norm_num)) (huS 0 (by norm_num))
      h0 i j i.isLt j.isLt
    rw [Finset.Nat.antidiagonal_zero, Finset.sum_singleton]
    exact scale_close scale _ _ hs _ (hpairInt 0 0 (by norm_num) (by norm_num) i j) hc
  · have hS := matAdd_spec _ _ (mS 0 1 (by norm_num) (by norm_num)) (mS 1 0 (by norm_num) (by norm_num))
    have hc := matLeScaled_sound scale _ _ hS.1 (huS 1 (by norm_num)) h1 i j i.isLt j.isLt
    rw [antidiag1, Matrix.add_apply]
    refine scale_close scale _ _ hs _ ?_ hc
    show _ ≤ (sparseLMat22 (sparseMatAddL22 _ _) i j : ℝ) / _
    rw [hS.2 i j i.isLt j.isLt, Int.cast_add, add_div]
    exact add_le_add (hpairInt 0 1 (by norm_num) (by norm_num) i j)
      (hpairInt 1 0 (by norm_num) (by norm_num) i j)
  · have hS1 := matAdd_spec _ _ (mS 0 2 (by norm_num) (by norm_num)) (mS 1 1 (by norm_num) (by norm_num))
    have hS := matAdd_spec _ _ hS1.1 (mS 2 0 (by norm_num) (by norm_num))
    have hc := matLeScaled_sound scale _ _ hS.1 (huS 2 (by norm_num)) h2 i j i.isLt j.isLt
    rw [antidiag2, Matrix.add_apply, Matrix.add_apply]
    refine scale_close scale _ _ hs _ ?_ hc
    show _ ≤ (sparseLMat22 (sparseMatAddL22 (sparseMatAddL22 _ _) _) i j : ℝ) / _
    rw [hS.2 i j i.isLt j.isLt, hS1.2 i j i.isLt j.isLt, Int.cast_add, Int.cast_add,
      add_div, add_div]
    exact add_le_add (add_le_add (hpairInt 0 2 (by norm_num) (by norm_num) i j)
      (hpairInt 1 1 (by norm_num) (by norm_num) i j)) (hpairInt 2 0 (by norm_num) (by norm_num) i j)

/-! ## Normalization -/

theorem normDen_eq (R d : ℕ) (hd : d ≤ 2) :
    sparseNormalizeCheckL22.sparseRegionDenominator22' R d = (Nat.choose (R * 16) d : ℤ) := by
  rw [← regionDenominator_eq R d hd]
  interval_cases d <;> rfl

theorem normalize_sound (R : ℕ) (hR : 0 < R) (scale : ℤ) (hs : 0 < scale)
    (u result : List (List (List Int))) (hu : ∀ d, d < 3 → MatShape (u.getD d []))
    (hcheck : sparseNormalizeCheckL22 R u result = true)
    (X : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hX : ∀ d, d ≤ 2 → ∀ i j : Fin 11, X d i j ≤ (sparseLPoly22 u d i j : ℝ) / scale) :
    (∀ d, d < 3 → MatShape (result.getD d [])) ∧
      ∀ d, d ≤ 2 → ∀ i j : Fin 11,
        ((1 / (Nat.choose (R * 16) d : ℝ)) • X d) i j ≤ (sparseLPoly22 result d i j : ℝ) / scale := by
  simp only [sparseNormalizeCheckL22, Bool.and_eq_true] at hcheck
  obtain ⟨⟨⟨hr, h0⟩, h1⟩, h2⟩ := hcheck
  have hrS := (polyShape_of_check result hr).2
  refine ⟨hrS, fun d hd i j => ?_⟩
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hc : sparseLMat22 (u.getD d []) i j ≤
      sparseNormalizeCheckL22.sparseRegionDenominator22' R d * sparseLMat22 (result.getD d []) i j := by
    interval_cases d
    · exact matLeScaled_sound _ _ _ (hu 0 (by norm_num)) (hrS 0 (by norm_num)) h0 i j i.isLt j.isLt
    · exact matLeScaled_sound _ _ _ (hu 1 (by norm_num)) (hrS 1 (by norm_num)) h1 i j i.isLt j.isLt
    · exact matLeScaled_sound _ _ _ (hu 2 (by norm_num)) (hrS 2 (by norm_num)) h2 i j i.isLt j.isLt
  rw [normDen_eq R d hd] at hc
  have hpos : (0 : ℝ) < (Nat.choose (R * 16) d : ℝ) := by
    exact_mod_cast Nat.choose_pos (by omega)
  have hcR : (sparseLPoly22 u d i j : ℝ) ≤
      (Nat.choose (R * 16) d : ℝ) * (sparseLPoly22 result d i j : ℝ) := by
    unfold sparseLPoly22; exact_mod_cast hc
  rw [Matrix.smul_apply, smul_eq_mul]
  calc 1 / (Nat.choose (R * 16) d : ℝ) * X d i j
      ≤ 1 / (Nat.choose (R * 16) d : ℝ) * ((sparseLPoly22 u d i j : ℝ) / scale) :=
        mul_le_mul_of_nonneg_left (hX d hd i j) (by positivity)
    _ ≤ _ := by
        rw [div_mul_div_comm, one_mul, div_le_div_iff₀ (mul_pos hpos hsR) hsR]
        nlinarith [mul_le_mul_of_nonneg_right hcR hsR.le]

end SparseListSound22

open SparseListSound22 in
#print axioms polyMul_sound
open SparseListSound22 in
#print axioms normalize_sound
end Spin.Turbo
