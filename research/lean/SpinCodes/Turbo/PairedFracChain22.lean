import SpinCodes.Turbo.SparseListPolyD22
import SpinCodes.Turbo.SparseRpowCert22
import SpinCodes.Turbo.MatrixPolynomial22

/-!
# Degree-D chain for the paired fractional profile DP (α = 4/5)

At α ≠ 1 the paired two-step polynomial is not a square of the single-step
polynomial, so the chain runs directly on
`orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B α)`.

* The squaring/product checks are the GREEN `sparsePolyMulCheckDL22` /
  `sparseRegionalPowerCheckDL22` (they are independent of the step family);
  `chainW_sound` re-proves the trace induction for a general nonnegative
  step family `W : Fin 17 × Fin 17 → Matrix`.
* Base stage: `orderedMatrixDP deg W 1 d = Σ_{deg ab = d} W ab` (`dp_one`).
* Pair certificates: `X_ab = R_a · R_b` (integer matrix product, scale `den²`)
  dominates `B a * B b`; `Y_ab` certifies `(X_ab/den²)^(4/5) ≤ Y_ab / s`
  (`sparseRpow45CheckL22`). The base polynomial satisfies
  `Σ_{a+b=d} C(16,a) C(16,b) Y_ab ≤ P_d` entrywise.

All kernel checks traverse lists; the only indexed accesses are into
17-entry lists (binomial literal, `R`, `Ys`, rows of `Ys`) and 11-row matrices.
-/
namespace Spin.Turbo

/-! ## Core-only definitions -/

/-- `C(16,a) · C(16,d-a)` when `a ≤ d` (zero once `d - a > 16`), else `0`. -/
def pairedCoefL22 (d a : Nat) : Int :=
  if a ≤ d then sparseBinom16L22.getD a 0 * sparseBinom16L22.getD (d - a) 0 else 0

/-- Entry `(i,k)` of `Σ_{a ≤ 16} C(16,a) C(16,d-a) · Y_{a,d-a}`. -/
def pairedBaseSumL22 (Ys : List (List (List (List Int)))) (d i k : Nat) : Int :=
  ((List.range 17).map fun a =>
    pairedCoefL22 d a * ((((Ys.getD a []).getD (d - a) []).getD i []).getD k 0)).sum

def pairedBaseGoL22 (Ys : List (List (List (List Int)))) :
    Nat → List (List (List Int)) → Bool
  | _, [] => true
  | d, m :: ms =>
    (sparseMatShapeL22 m &&
      (List.range 11).all fun i => (List.range 11).all fun k =>
        decide (pairedBaseSumL22 Ys d i k ≤ (m.getD i []).getD k 0)) &&
      pairedBaseGoL22 Ys (d + 1) ms

/-- Base stage: `D+1` shaped matrices `P_d` with
`Σ_{a+b=d} C(16,a) C(16,b) Y_ab ≤ P_d` entrywise (for `d > 32` this is `0 ≤ P_d`). -/
def pairedBaseCheckDL22 (D : Nat) (Ys : List (List (List (List Int))))
    (base : List (List (List Int))) : Bool :=
  decide (base.length = D + 1) && pairedBaseGoL22 Ys 0 base

/-- Pair certificates: `0 < den`, each `R_a` shaped, and for all `a, b ≤ 16`,
`Y_ab` certifies `((R_a R_b)/den²)^(4/5) ≤ Y_ab / sy` entrywise (with `Y_ab ≥ 0`). -/
def pairedRpowCheckL22 (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) : Bool :=
  decide (0 < den) && (List.range 17).all fun a =>
    sparseMatShapeL22 (R.getD a []) && (List.range 17).all fun b =>
      sparseRpow45CheckL22 (den * den) sy
        (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b [])

/-- Row matrices from a row function (bridge from the `baseD_sound` convention). -/
def pairedRowMatL22 (row : Nat → Nat → Nat → Int) (a : Nat) : List (List Int) :=
  (List.range 11).map fun i => (List.range 11).map fun k => row a i k

def pairedRowMatsL22 (row : Nat → Nat → Nat → Int) : List (List (List Int)) :=
  (List.range 17).map (pairedRowMatL22 row)

end Spin.Turbo

/-! ## Soundness -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 SparseListPolyD22

namespace PairedFracChain22

theorem dp_nonneg (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hW : ∀ ab i j, 0 ≤ W ab i j) (n d : ℕ) (i j : Fin 11) :
    0 ≤ orderedMatrixDP pairedProfileDegree22 W n d i j := by
  rw [← orderedMatrixCoefficient_eq_dp]
  unfold orderedMatrixCoefficient
  rw [Matrix.sum_apply]
  apply Finset.sum_nonneg
  intro xs _
  split_ifs
  · exact listNonneg _ (fun M hM => by
      obtain ⟨step, rfl⟩ := List.mem_ofFn.mp hM
      exact hW _) i j
  · exact le_rfl

abbrev PolyBoundW (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (scale : ℤ) (n : ℕ)
    (P : List (List (List Int))) (D : ℕ) : Prop :=
  (∀ d, d ≤ D → MatShape (P.getD d [])) ∧
    ∀ d, d ≤ D → ∀ i j : Fin 11,
      orderedMatrixDP pairedProfileDegree22 W n d i j ≤ (sparseLPoly22 P d i j : ℝ) / scale

theorem add_boundW (D : ℕ) (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hW : ∀ ab i j, 0 ≤ W ab i j) (scale : ℤ) (hs : 0 < scale)
    (a b u : List (List (List Int))) (m n : ℕ)
    (hcheck : sparsePolyMulCheckDL22 D scale a b u = true)
    (ha : PolyBoundW W scale m a D) (hb : PolyBoundW W scale n b D) :
    PolyBoundW W scale (m + n) u D := by
  have h := polyMulD_sound D scale hs a b u ha.1 hb.1 hcheck
    (orderedMatrixDP pairedProfileDegree22 W m) (orderedMatrixDP pairedProfileDegree22 W n)
    (fun d => dp_nonneg W hW m d) (fun d => dp_nonneg W hW n d) ha.2 hb.2
  refine ⟨h.1, fun d hd i j => ?_⟩
  rw [orderedMatrixDP_add]
  exact h.2 d hd i j

theorem chainW_sound (D : ℕ) (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hW : ∀ ab i j, 0 ≤ W ab i j) (scale : ℤ) (hs : 0 < scale)
    (base : List (List (List Int))) (hbase : PolyBoundW W scale 1 base D) :
    ∀ (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
      (prev : List (List (List Int))) (e : ℕ),
      PolyBoundW W scale e prev D →
      sparseRegionalPowerCheckDL22 D scale base prev trace = true →
      PolyBoundW W scale (trace.foldl sparseRegionalStepL22 e)
        (sparseRegionalLastL22 prev trace) D
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((sparsePolyMulCheckDL22 D scale prev prev sq &&
        (if bit then sparsePolyMulCheckDL22 D scale sq base next else true)) &&
        sparseRegionalPowerCheckDL22 D scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := add_boundW D W hW scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : sparseRegionalPowerCheckDL22 D scale base sq tail = true := hrest
      have he : sparseRegionalStepL22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hsq' : PolyBoundW W scale (sparseRegionalStepL22 e (false, sq, next)) sq D := by
        rw [he]; exact hsqB
      exact chainW_sound D W hW scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : sparsePolyMulCheckDL22 D scale sq base next = true := hbit
      have hrest' : sparseRegionalPowerCheckDL22 D scale base next tail = true := hrest
      have hnext := add_boundW D W hW scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStepL22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hnext' : PolyBoundW W scale (sparseRegionalStepL22 e (true, sq, next)) next D := by
        rw [he]; exact hnext
      exact chainW_sound D W hW scale hs base hbase tail next _ hnext' hrest'

/-! ## Base stage -/

theorem dp_one (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (d : ℕ) :
    orderedMatrixDP pairedProfileDegree22 W 1 d =
      ∑ ab, if pairedProfileDegree22 ab = d then W ab else 0 := by
  simp only [orderedMatrixDP]
  refine Finset.sum_congr rfl fun ab _ => ?_
  by_cases h1 : pairedProfileDegree22 ab ≤ d
  · by_cases h2 : pairedProfileDegree22 ab = d
    · rw [if_pos h1, if_pos h2, h2, Nat.sub_self, if_pos rfl, mul_one]
    · rw [if_pos h1, if_neg h2, if_neg (by omega), mul_zero]
  · rw [if_neg h1, if_neg (by omega)]

theorem list_sum_range (f : ℕ → ℤ) :
    ∀ n, ((List.range n).map f).sum = ∑ i ∈ range n, f i
  | 0 => by simp
  | n + 1 => by
    rw [List.range_succ, List.map_append, List.sum_append, list_sum_range f n,
      Finset.sum_range_succ]
    simp

theorem baseGo_sound (Ys : List (List (List (List Int)))) :
    ∀ (ms : List (List (List Int))) (d : ℕ), pairedBaseGoL22 Ys d ms = true →
      ∀ e, e < ms.length → MatShape (ms.getD e []) ∧ ∀ i k, i < 11 → k < 11 →
        pairedBaseSumL22 Ys (d + e) i k ≤ sparseLMat22 (ms.getD e []) i k
  | [], _, _ => fun e he => by simp at he
  | m :: ms, d, hc => by
    have hc' : ((sparseMatShapeL22 m &&
        (List.range 11).all fun i => (List.range 11).all fun k =>
          decide (pairedBaseSumL22 Ys d i k ≤ (m.getD i []).getD k 0)) &&
        pairedBaseGoL22 Ys (d + 1) ms) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hm, hent⟩, hrest⟩ := hc'
    simp only [List.all_eq_true, List.mem_range, decide_eq_true_eq] at hent
    intro e he
    cases e with
    | zero =>
      refine ⟨by simpa using matShape_of_check m hm, fun i k hi hk => ?_⟩
      exact hent i hi k hk
    | succ e =>
      have h := baseGo_sound Ys ms (d + 1) hrest e (by simp at he; omega)
      rw [List.getD_cons_succ]
      refine ⟨h.1, fun i k hi hk => ?_⟩
      have h2 := h.2 i k hi hk
      rwa [show d + 1 + e = d + (e + 1) by omega] at h2

/-- Base soundness for any step family dominated by the `Y` certificates. -/
theorem base_sound (D : ℕ) (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (s : ℤ) (hs : 0 < s) (Ys : List (List (List (List Int))))
    (hW : ∀ (a b : Fin 17) (i k : Fin 11), W (a, b) i k ≤
      ((Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b.val : ℝ)) *
        (sparseLMat22 ((Ys.getD a.val []).getD b.val []) i.val k.val : ℝ) / s)
    (base : List (List (List Int))) (hc : pairedBaseCheckDL22 D Ys base = true) :
    PolyBoundW W s 1 base D := by
  simp only [pairedBaseCheckDL22, Bool.and_eq_true, decide_eq_true_eq] at hc
  obtain ⟨hlen, hgo⟩ := hc
  have hG := baseGo_sound Ys base 0 hgo
  refine ⟨fun d hd => (hG d (by omega)).1, fun d hd i k => ?_⟩
  have hck := (hG d (by omega)).2 i k i.isLt k.isLt
  rw [Nat.zero_add] at hck
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  rw [dp_one, Matrix.sum_apply]
  have hsum : (∑ ab : Fin 17 × Fin 17,
      (if pairedProfileDegree22 ab = d then W ab else 0) i k) =
      ∑ a : Fin 17, ∑ b : Fin 17, (if a.val + b.val = d then W (a, b) i k else 0) := by
    rw [Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    by_cases h : a.val + b.val = d
    · rw [if_pos (show pairedProfileDegree22 (a, b) = d from h), if_pos h]
    · rw [if_neg (show ¬ pairedProfileDegree22 (a, b) = d from h), if_neg h]; rfl
  rw [hsum]
  have hinner : ∀ a : Fin 17,
      (∑ b : Fin 17, (if a.val + b.val = d then W (a, b) i k else 0)) ≤
        (pairedCoefL22 d a.val : ℝ) *
          (sparseLMat22 ((Ys.getD a.val []).getD (d - a.val) []) i.val k.val : ℝ) / s := by
    intro a
    by_cases hab : a.val ≤ d ∧ d - a.val ≤ 16
    · let b0 : Fin 17 := ⟨d - a.val, by omega⟩
      rw [Finset.sum_eq_single b0]
      · rw [if_pos (show a.val + b0.val = d by simp only [b0]; omega)]
        have hcf : (pairedCoefL22 d a.val : ℝ) =
            (Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b0.val : ℝ) := by
          unfold pairedCoefL22
          rw [if_pos hab.1, binom16L_eq _ (by omega), binom16L_eq _ hab.2]
          push_cast
          rfl
        rw [hcf]
        exact hW a b0 i k
      · intro b _ hb
        rw [if_neg]
        intro h
        apply hb
        ext
        simp only [b0]
        omega
      · intro h
        exact absurd (Finset.mem_univ _) h
    · rw [Finset.sum_eq_zero (fun b _ => if_neg (by
        intro h
        apply hab
        have := b.isLt
        omega))]
      have hcf : pairedCoefL22 d a.val = 0 := by
        unfold pairedCoefL22
        split_ifs with h1
        · rw [binom16L_high (d - a.val) (by omega), mul_zero]
        · rfl
      rw [hcf]
      simp
  refine (Finset.sum_le_sum fun a _ => hinner a).trans ?_
  have hlist : (pairedBaseSumL22 Ys d i.val k.val : ℝ) =
      ∑ a : Fin 17, (pairedCoefL22 d a.val : ℝ) *
        (sparseLMat22 ((Ys.getD a.val []).getD (d - a.val) []) i.val k.val : ℝ) := by
    unfold pairedBaseSumL22
    rw [list_sum_range, Int.cast_sum, ← Fin.sum_univ_eq_sum_range]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Int.cast_mul]
    rfl
  rw [← Finset.sum_div, ← hlist]
  apply div_le_div_of_nonneg_right _ hsR.le
  unfold sparseLPoly22
  exact_mod_cast hck

/-! ## Pair certificates -/

theorem pairedProfile_nonneg (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (α : ℝ) (ab : Fin 17 × Fin 17) (i j : Fin 11) :
    0 ≤ pairedProfileMatrix22 B α ab i j := by
  unfold pairedProfileMatrix22
  rw [Matrix.smul_apply, smul_eq_mul]
  apply mul_nonneg (by positivity)
  apply Real.rpow_nonneg
  rw [Matrix.mul_apply]
  exact Finset.sum_nonneg fun k _ => mul_nonneg (hB _ _ _) (hB _ _ _)

theorem rpowCheck_entry (den sy : ℤ) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (h : pairedRpowCheckL22 den sy R Ys = true)
    (a b : ℕ) (ha : a < 17) (hb : b < 17) :
    0 < den ∧ MatShape (R.getD a []) ∧ MatShape (R.getD b []) ∧
      sparseRpow45CheckL22 (den * den) sy
        (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b []) = true := by
  simp only [pairedRpowCheckL22, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true,
    List.mem_range] at h
  obtain ⟨hden, hall⟩ := h
  exact ⟨hden, matShape_of_check _ (hall a ha).1, matShape_of_check _ (hall b hb).1,
    (hall a ha).2 b hb⟩

/-- Each paired step matrix at α = 4/5 is dominated by its `Y` certificate. -/
theorem pairedProfile_le (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (den s : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowCheckL22 den s R Ys = true)
    (a b : Fin 17) (i k : Fin 11) :
    pairedProfileMatrix22 B ((4 : ℝ) / 5) (a, b) i k ≤
      ((Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b.val : ℝ)) *
        (sparseLMat22 ((Ys.getD a.val []).getD b.val []) i.val k.val : ℝ) / s := by
  obtain ⟨hden, hRa, hRb, hck⟩ := rpowCheck_entry den s R Ys hrp a.val b.val a.isLt b.isLt
  have hdenR : (0 : ℝ) < den := by exact_mod_cast hden
  have hM0 : ∀ i j, 0 ≤ (B a * B b) i j := fun i j => by
    rw [Matrix.mul_apply]
    exact Finset.sum_nonneg fun k _ => mul_nonneg (hB _ _ _) (hB _ _ _)
  have hMX : ∀ i j : Fin 11, (B a * B b) i j ≤
      (sparseLMat22 (sparseMatMulL22 (R.getD a.val []) (R.getD b.val [])) i.val j.val : ℝ) /
        ((den * den : ℤ) : ℝ) := by
    intro i j
    rw [Matrix.mul_apply, (matMul_spec _ _ hRa hRb).2 i j i.isLt j.isLt, Int.cast_sum,
      Int.cast_mul]
    simp only [Int.cast_mul]
    exact dot_le (den : ℝ) hdenR _ _ _ _ (fun k => hB a i k) (fun k => hB b k j)
      (fun k => hR a i k) (fun k => hR b k j)
  have hY := SparseRpowCert22.sparseRpow45CheckL22_entrywiseRpow hck (B a * B b) hM0 hMX i k
  unfold pairedProfileMatrix22
  rw [Matrix.smul_apply, smul_eq_mul, mul_div_assoc]
  exact mul_le_mul_of_nonneg_left hY (by positivity)

/-- Full chain: pair certificates + base check + squaring trace give a degree-`D`
bound on the paired fractional DP after `trace.foldl sparseRegionalStepL22 1` pairs. -/
theorem pairedChain_sound (D : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (s : ℤ) (hs : 0 < s) (den : ℤ)
    (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowCheckL22 den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (hchain : sparseRegionalPowerCheckDL22 D s base base trace = true) :
    PolyBoundW (pairedProfileMatrix22 B ((4 : ℝ) / 5)) s
      (trace.foldl sparseRegionalStepL22 1) (sparseRegionalLastL22 base trace) D := by
  have hW0 := pairedProfile_nonneg B hB ((4 : ℝ) / 5)
  have hb := base_sound D (pairedProfileMatrix22 B ((4 : ℝ) / 5)) s hs Ys
    (pairedProfile_le B hB den s R hR Ys hrp) base hbase
  exact chainW_sound D _ hW0 s hs base hb trace base 1 hb hchain

/-- Bridge from the row-function convention of `baseD_sound`. -/
theorem rowMats_get (row : ℕ → ℕ → ℕ → ℤ) (a i k : ℕ) (ha : a < 17) (hi : i < 11)
    (hk : k < 11) : sparseLMat22 ((pairedRowMatsL22 row).getD a []) i k = row a i k := by
  unfold sparseLMat22 pairedRowMatsL22
  rw [gd ((List.range 17).map (pairedRowMatL22 row)) [] (by simp [ha]), List.getElem_map,
    List.getElem_range]
  unfold pairedRowMatL22
  rw [gd ((List.range 11).map fun i => (List.range 11).map fun k => row a i k) []
      (by simp [hi]), List.getElem_map, List.getElem_range,
    gd ((List.range 11).map fun k => row a i k) 0 (by simp [hk]), List.getElem_map,
    List.getElem_range]

theorem rowMats_hyp (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (row : ℕ → ℕ → ℕ → ℤ)
    (den : ℤ) (hrow : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / den) :
    ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 ((pairedRowMatsL22 row).getD a.val []) i.val k.val : ℝ) / den := by
  intro a i k
  rw [rowMats_get row a.val i.val k.val a.isLt i.isLt k.isLt]
  exact hrow a.val (by omega) i k

end PairedFracChain22

open PairedFracChain22 in
#print axioms chainW_sound
open PairedFracChain22 in
#print axioms base_sound
open PairedFracChain22 in
#print axioms pairedProfile_le
open PairedFracChain22 in
#print axioms pairedChain_sound
open PairedFracChain22 in
#print axioms rowMats_hyp
end Spin.Turbo

/-! ## Kernel smoke tests (all-zero toy data) -/
namespace Spin.Turbo.PairedFracChain22Smoke

def zMat : List (List Int) := List.replicate 11 (List.replicate 11 0)
def zR : List (List (List Int)) := List.replicate 17 zMat
def zYs : List (List (List (List Int))) := List.replicate 17 (List.replicate 17 zMat)

example : Spin.Turbo.pairedRpowCheckL22 1 1 zR zYs = true := by decide +kernel
example : Spin.Turbo.pairedBaseCheckDL22 2 zYs (List.replicate 3 zMat) = true := by
  decide +kernel
example : Spin.Turbo.sparseRegionalPowerCheckDL22 1 1 (List.replicate 2 zMat)
    (List.replicate 2 zMat) [(false, List.replicate 2 zMat, List.replicate 2 zMat)] = true := by
  decide +kernel

end Spin.Turbo.PairedFracChain22Smoke
