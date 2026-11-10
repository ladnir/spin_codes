import SpinCodes.Turbo.SparseListSound22
import SpinCodes.Turbo.SparseRegionalPower22
import SpinCodes.Turbo.SparseRowBase22

/-! List-form regional certificate and q1 row DP soundness. Same structure as
the array versions (`SparseRegionalPower22`, `SparseFixedSound22`): one check
per stored stage plus an induction theorem. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 SparseRegionalPower22

namespace SparseListQ1Regional22

abbrev PolyBound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (scale : ℤ) (n : ℕ)
    (P : List (List (List Int))) : Prop :=
  (∀ d, d < 3 → MatShape (P.getD d [])) ∧
    ∀ d, d ≤ 2 → ∀ i j : Fin 11,
      regionalCountDP22 B n d i j ≤ (sparseLPoly22 P d i j : ℝ) / scale

theorem add_bound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (scale : ℤ) (hs : 0 < scale) (a b u : List (List (List Int))) (m n : ℕ)
    (hcheck : sparsePolyMulCheckL22 scale a b u = true)
    (ha : PolyBound B scale m a) (hb : PolyBound B scale n b) :
    PolyBound B scale (m + n) u := by
  have h := SparseListSound22.polyMul_sound scale hs a b u ha.1 hb.1 hcheck
    (regionalCountDP22 B m) (regionalCountDP22 B n)
    (count_nonneg B hB m) (count_nonneg B hB n) ha.2 hb.2
  refine ⟨h.1, fun d hd i j => ?_⟩
  unfold regionalCountDP22
  rw [orderedMatrixDP_add]
  exact h.2 d hd i j

theorem chain_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (base : List (List (List Int))) (hbase : PolyBound B scale 1 base) :
    ∀ (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
      (prev : List (List (List Int))) (e : ℕ),
      PolyBound B scale e prev →
      sparseRegionalPowerCheckL22 scale base prev trace = true →
      PolyBound B scale (trace.foldl sparseRegionalStepL22 e)
        (sparseRegionalLastL22 prev trace)
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((sparsePolyMulCheckL22 scale prev prev sq &&
        (if bit then sparsePolyMulCheckL22 scale sq base next else true)) &&
        sparseRegionalPowerCheckL22 scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := add_bound B hB scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : sparseRegionalPowerCheckL22 scale base sq tail = true := hrest
      have he : sparseRegionalStepL22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hsq' : PolyBound B scale (sparseRegionalStepL22 e (false, sq, next)) sq := by
        rw [he]; exact hsqB
      exact chain_sound B hB scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : sparsePolyMulCheckL22 scale sq base next = true := hbit
      have hrest' : sparseRegionalPowerCheckL22 scale base next tail = true := hrest
      have hnext := add_bound B hB scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStepL22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hnext' : PolyBound B scale (sparseRegionalStepL22 e (true, sq, next)) next := by
        rw [he]; exact hnext
      exact chain_sound B hB scale hs base hbase tail next _ hnext' hrest'

theorem base_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (row : ℕ → ℕ → ℕ → ℤ) (den mult : ℤ) (hden : 0 < den) (hmult : 0 < mult)
    (hrow : ∀ d (hd : d ≤ 2), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k = (row d i.val k.val : ℝ) / den)
    (base : List (List (List Int)))
    (hcheck : sparseRegionalBaseCheckL22 row mult base = true) :
    PolyBound B (den * mult) 1 base := by
  simp only [sparseRegionalBaseCheckL22, Bool.and_eq_true, List.all_eq_true,
    List.mem_range, decide_eq_true_eq] at hcheck
  obtain ⟨hshape, hent⟩ := hcheck
  refine ⟨(polyShape_of_check base hshape).2, fun d hd i k => ?_⟩
  have hc := hent d (by omega) i.val i.isLt k.val k.isLt
  have hcR : (sparseBinom16Small22 d : ℝ) *
      (row d i.val k.val : ℝ) * (mult : ℝ) ≤ (sparseLPoly22 base d i k : ℝ) := by
    unfold sparseLPoly22 sparseLMat22; exact_mod_cast hc
  rw [SparseRowBase22.binom16_eq d hd] at hcR
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have hmR : (0 : ℝ) < mult := by exact_mod_cast hmult
  rw [count_one B d (by omega), Matrix.smul_apply, smul_eq_mul, hrow d hd i k,
    Int.cast_mul, mul_div_assoc', div_le_div_iff₀ hdR (mul_pos hdR hmR)]
  nlinarith [mul_le_mul_of_nonneg_right hcR hdR.le]

/-- A passing list certificate bounds every normalized occupancy kernel
entry used by the sparse q1/q2 kernels. -/
theorem certificate_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (R : ℕ) (s den mult : ℤ) (hsd : s = den * mult)
    (hden : 0 < den) (row : ℕ → ℕ → ℕ → ℤ)
    (hrow : ∀ d (hd : d ≤ 2), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k = (row d i.val k.val : ℝ) / den)
    (base : List (List (List Int)))
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (result : List (List (List Int)))
    (hcheck : sparseRegionalCertificateCheckL22 R s row mult base trace result = true) :
    (∀ d, d < 3 → MatShape (result.getD d [])) ∧
      ∀ d, d ≤ 2 → ∀ i j : Fin 11,
        sparseRegionalDP22 B R d i j ≤ (sparseLPoly22 result d i j : ℝ) / s := by
  simp only [sparseRegionalCertificateCheckL22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨⟨⟨hs, hm, hR, hexp⟩, hbaseC⟩, hchain⟩, hnorm⟩ := hcheck
  have hbase : PolyBound B s 1 base := by
    rw [hsd]; exact base_sound B row den mult hden hm hrow base hbaseC
  have hc := chain_sound B hB s hs base hbase trace base 1 hbase hchain
  rw [hexp] at hc
  exact SparseListSound22.normalize_sound R hR s hs _ result hc.1 hnorm
    (regionalCountDP22 B R) hc.2

/-! ## q1 rows -/

def RowsShape (rows : List (List Int)) : Prop := rows.length = 65 ∧ ∀ r ∈ rows, r.length = 11

theorem zero11_getD (i : ℕ) : sparseZero11L22.getD i 0 = 0 := by
  unfold sparseZero11L22
  rcases i with _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | i <;> rfl

theorem q1Rows_sound (scale : ℤ) (hs : 0 < scale) (r0 r1 : List (List Int))
    (hr0S : MatShape r0) (hr1S : MatShape r1)
    (W : Bool → Matrix (Fin 11) (Fin 11) ℝ) (hW0 : ∀ b i j, 0 ≤ W b i j)
    (hr0 : ∀ i j : Fin 11, W false i j ≤ (sparseLMat22 r0 i j : ℝ) / scale)
    (hr1 : ∀ i j : Fin 11, W true i j ≤ (sparseLMat22 r1 i j : ℝ) / scale)
    (rows : ℕ → List (List Int)) (N : ℕ)
    (hinit : rows 0 = sparseQ1InitialRowsL22 scale)
    (hstep : ∀ n, n < N → sparseQ1StepCheckL22 scale r0 r1 (rows n) (rows (n + 1)) = true) :
    ∀ n, n ≤ N → RowsShape (rows n) ∧ ∀ v, v < 65 → ∀ j : Fin 11,
      orderedMatrixDP sparseBoolDegree22 W n v 0 j ≤
        (sparseLMat22 (rows n) v j : ℝ) / scale := by
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hdp0 : ∀ n q i j, 0 ≤ orderedMatrixDP sparseBoolDegree22 W n q i j := by
    intro n q i j
    rw [← orderedMatrixCoefficient_eq_dp]
    unfold orderedMatrixCoefficient
    rw [Matrix.sum_apply]
    apply Finset.sum_nonneg
    intro xs _
    split_ifs
    · exact listNonneg _ (fun M hM => by
        obtain ⟨step, rfl⟩ := List.mem_ofFn.mp hM
        exact hW0 _) i j
    · exact le_rfl
  intro n
  induction n with
  | zero =>
    intro _
    rw [hinit]
    refine ⟨⟨by simp [sparseQ1InitialRowsL22], fun r hr => ?_⟩, fun v hv j => ?_⟩
    · simp only [sparseQ1InitialRowsL22, List.mem_cons, List.mem_replicate] at hr
      rcases hr with rfl | ⟨_, rfl⟩
      · simp
      · rfl
    · simp only [orderedMatrixDP]
      rcases v with _ | v
      · rw [if_pos rfl, Matrix.one_apply]
        have hval : sparseLMat22 (sparseQ1InitialRowsL22 scale) 0 j =
            if (0 : Fin 11) = j then scale else 0 := by
          unfold sparseLMat22 sparseQ1InitialRowsL22
          fin_cases j <;> rfl
        rw [hval]
        split_ifs
        · rw [div_self hsR.ne']
        · simp
      · rw [if_neg (Nat.succ_ne_zero v), Matrix.zero_apply]
        have hval : sparseLMat22 (sparseQ1InitialRowsL22 scale) (v + 1) j = 0 := by
          unfold sparseLMat22 sparseQ1InitialRowsL22
          rw [List.getD_cons_succ, gd (List.replicate 64 sparseZero11L22) []
            (by rw [List.length_replicate]; omega), List.getElem_replicate]
          exact zero11_getD _
        rw [hval]
        simp
  | succ n ih =>
    intro hn
    obtain ⟨⟨hpl, hpr⟩, ihv⟩ := ih (by omega)
    have hc := hstep n (by omega)
    simp only [sparseQ1StepCheckL22, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true,
      sparseLen11L22] at hc
    obtain ⟨⟨hnl, hnr⟩, hall⟩ := hc
    refine ⟨⟨hnl, hnr⟩, fun v hv j => ?_⟩
    set prev := rows n
    set next := rows (n + 1)
    have hvP : v < prev.length := by omega
    have hvN : v < next.length := by omega
    have hzl : (next.zip (prev.zip (sparseZero11L22 :: prev))).length = 65 := by
      simp [hnl, hpl]
    have hvZ : v < (next.zip (prev.zip (sparseZero11L22 :: prev))).length := by omega
    have ht := hall _ (List.getElem_mem hvZ)
    simp only [List.getElem_zip] at ht
    have hprevRow : prev[v].length = 11 := hpr _ (List.getElem_mem hvP)
    have hnextRow : next[v].length = 11 := hnr _ (List.getElem_mem hvN)
    have hshift : (sparseZero11L22 :: prev)[v]'(by simp; omega) =
        if h : v = 0 then sparseZero11L22 else prev[v - 1]'(by omega) := by
      rcases v with _ | v
      · rfl
      · simp
    have hshiftLen : ((sparseZero11L22 :: prev)[v]'(by simp; omega)).length = 11 := by
      rw [hshift]; split_ifs
      · rfl
      · exact hpr _ (List.getElem_mem _)
    have hX := vecMat_spec _ r0 hprevRow hr0S
    have hY := vecMat_spec _ r1 hshiftLen hr1S
    have hle := leScaled_sound scale _ _ (by simp [sparseVecAddL22, hX.1, hY.1]) hnextRow ht
      j.val j.isLt
    rw [show sparseVecAddL22 = List.zipWith (· + ·) from rfl,
      zipWith_getD hX.1 hY.1 j.isLt, hX.2 j.val j.isLt, hY.2 j.val j.isLt] at hle
    have hnextv : sparseLMat22 next v j = next[v].getD j.val 0 := by
      unfold sparseLMat22; rw [gd next [] hvN]
    rw [hnextv]
    apply scale_close scale _ _ hs _ _ hle
    rw [Int.cast_add, add_div, Int.cast_sum, Int.cast_sum]
    simp only [Int.cast_mul]
    rcases v with _ | v
    · rw [dp_succ_zero, Matrix.mul_apply]
      have hz : ∀ i : Fin 11,
          (((sparseZero11L22 :: prev)[0]'(by simp)).getD i.val 0 : ℝ) = 0 := by
        intro i; show ((sparseZero11L22.getD i.val 0 : ℤ) : ℝ) = 0
        rw [zero11_getD]; simp
      simp only [hz, zero_mul, Finset.sum_const_zero, zero_div, add_zero]
      exact dot_le _ hsR _ _ _ _ (fun k => hdp0 n 0 0 k) (fun k => hW0 false k j)
        (fun k => by
          have := ihv 0 (by norm_num) k
          unfold sparseLMat22 at this; rwa [gd prev [] hvP] at this)
        (fun k => hr0 k j)
    · rw [dp_succ_succ, Matrix.add_apply, Matrix.mul_apply, Matrix.mul_apply]
      have hsv : (sparseZero11L22 :: prev)[v + 1]'(by simp; omega) = prev[v]'(by omega) := by
        simp
      simp only [hsv]
      exact add_le_add
        (dot_le _ hsR _ _ _ _ (fun k => hdp0 n (v + 1) 0 k) (fun k => hW0 false k j)
          (fun k => by
            have := ihv (v + 1) hv k
            unfold sparseLMat22 at this; rwa [gd prev [] hvP] at this)
          (fun k => hr0 k j))
        (dot_le _ hsR _ _ _ _ (fun k => hdp0 n v 0 k) (fun k => hW0 true k j)
          (fun k => by
            have := ihv v (by omega) k
            unfold sparseLMat22 at this; rwa [gd prev [] (by omega)] at this)
          (fun k => hr1 k j))

end SparseListQ1Regional22

open SparseListQ1Regional22 in
#print axioms certificate_sound
open SparseListQ1Regional22 in
#print axioms q1Rows_sound
end Spin.Turbo
