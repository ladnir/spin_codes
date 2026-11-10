import SpinCodes.Turbo.SparseListQ1Regional22

/-!
# Degree-D list-structured polynomial-matrix checks

Generalizes the degree-≤2 checks of `SparseListDefs22` to an arbitrary
truncation degree `D`. A polynomial is a list of `D+1` matrices (11×11, lists
of rows). All checks traverse lists in lockstep:

* the product check walks `u` and `b` together, keeping the reversed prefix
  `b_d, …, b_0` of `b`; the degree-`d` convolution is then a lockstep walk of
  `a` against that reversed prefix (no indexed access into long lists);
* the base check walks the base stage once; the only indexed access is into
  the 17-entry binomial literal and into 11-row matrices.

Soundness follows `SparseListSound22` / `SparseListQ1Regional22`: one check
per stored stage plus an induction over the MSB-first squaring trace.
-/
namespace Spin.Turbo

/-! ## Core-only definitions -/

def sparseZeroMatL22 : List (List Int) := List.replicate 11 sparseZero11L22

/-- `acc + Σ_t xs_t · ys_t`, walking both lists in lockstep (stops at the
shorter one). -/
def sparsePolyConvL22 : List (List (List Int)) → List (List (List Int)) →
    List (List Int) → List (List Int)
  | x :: xs, y :: ys, acc => sparsePolyConvL22 xs ys (sparseMatAddL22 acc (sparseMatMulL22 x y))
  | _, _, acc => acc

/-- Walk `u` (output degree `d`) and `b` together. `rev = [b_{d-1}, …, b_0]`. -/
def sparsePolyMulGoDL22 (s : Int) (a : List (List (List Int))) :
    List (List (List Int)) → List (List (List Int)) → List (List (List Int)) → Bool
  | _, _, [] => true
  | _, [], _ :: _ => false
  | rev, bm :: bs, m :: us =>
    (sparseMatShapeL22 m &&
      sparseMatLeScaledL22 s (sparsePolyConvL22 a (bm :: rev) sparseZeroMatL22) m) &&
      sparsePolyMulGoDL22 s a (bm :: rev) bs us

/-- `u` has `D+1` shaped matrices and, for every `d ≤ D`,
`Σ_{r ≤ d} a_r · b_{d-r} ≤ s · u_d` entrywise. `a`, `b` are assumed shaped
for degrees `≤ D` (checked where they are produced). -/
def sparsePolyMulCheckDL22 (D : Nat) (s : Int) (a b u : List (List (List Int))) : Bool :=
  decide (u.length = D + 1) && sparsePolyMulGoDL22 s a [] b u

/-- Degree-D squaring trace; exponent bookkeeping and last stage reuse
`sparseRegionalStepL22` (`2e + bit`) and `sparseRegionalLastL22`. -/
def sparseRegionalPowerCheckDL22 (D : Nat) (s : Int) (base prev : List (List (List Int))) :
    List (Bool × List (List (List Int)) × List (List (List Int))) → Bool
  | [] => true
  | (bit, sq, next) :: tail =>
    (sparsePolyMulCheckDL22 D s prev prev sq &&
      (if bit then sparsePolyMulCheckDL22 D s sq base next else true)) &&
      sparseRegionalPowerCheckDL22 D s base (if bit then next else sq) tail

/-- `choose 16 d` for `d ≤ 16` as a literal (no `Nat.choose` in kernel code). -/
def sparseBinom16L22 : List Int :=
  [1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1]

def sparseRegionalBaseGoDL22 (row : Nat → Nat → Nat → Int) (mult : Int) :
    Nat → List (List (List Int)) → Bool
  | _, [] => true
  | d, m :: ms =>
    (sparseMatShapeL22 m &&
      (List.range 11).all fun i => (List.range 11).all fun k =>
        decide (sparseBinom16L22.getD d 0 * row d i k * mult ≤ (m.getD i []).getD k 0)) &&
      sparseRegionalBaseGoDL22 row mult (d + 1) ms

/-- Base stage: `D+1` shaped matrices with `choose(16,d)·row d i k·mult ≤ base_d[i][k]`
for `d ≤ min(D,16)` and `0 ≤ base_d[i][k]` for `16 < d ≤ D`. -/
def sparseRegionalBaseCheckDL22 (D : Nat) (row : Nat → Nat → Nat → Int) (mult : Int)
    (base : List (List (List Int))) : Bool :=
  decide (base.length = D + 1) && sparseRegionalBaseGoDL22 row mult 0 base

end Spin.Turbo

/-! ## Soundness -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 SparseRegionalPower22

namespace SparseListPolyD22

theorem getD_ge {α : Type*} (x : α) : ∀ (l : List α) (n : ℕ), l.length ≤ n → l.getD n x = x
  | [], _, _ => by simp
  | _ :: _, 0, h => by simp at h
  | _ :: l, n + 1, h => by
    rw [List.getD_cons_succ]
    exact getD_ge x l n (by simp at h; omega)

theorem len_of_shape (l : List (List (List Int))) (d : ℕ) (h : MatShape (l.getD d [])) :
    d < l.length := by
  by_contra hc
  rw [getD_ge [] l d (by omega)] at h
  have := h.1
  simp at this

theorem drop_cons {α : Type*} (x : α) : ∀ (l : List α) (d : ℕ) (y : α) (ys : List α),
    l.drop d = y :: ys → l.getD d x = y ∧ l.drop (d + 1) = ys
  | [], _, _, _, h => by simp at h
  | z :: l, 0, y, ys, h => by
    simp only [List.drop_zero, List.cons.injEq] at h
    obtain ⟨rfl, rfl⟩ := h
    simp
  | _ :: l, d + 1, y, ys, h => by
    simp only [List.drop_succ_cons, List.getD_cons_succ] at h ⊢
    exact drop_cons x l d y ys h

theorem zeroMat_shape : MatShape sparseZeroMatL22 := by
  refine ⟨by simp [sparseZeroMatL22], fun r hr => ?_⟩
  simp only [sparseZeroMatL22, List.mem_replicate] at hr
  rw [hr.2]; rfl

theorem zeroMat_get (i j : ℕ) (hi : i < 11) : sparseLMat22 sparseZeroMatL22 i j = 0 := by
  unfold sparseLMat22 sparseZeroMatL22
  rw [gd (List.replicate 11 sparseZero11L22) [] (by simpa using hi), List.getElem_replicate]
  exact SparseListQ1Regional22.zero11_getD j

theorem conv_spec : ∀ (n : ℕ) (a rev : List (List (List Int))) (acc : List (List Int)),
    n ≤ a.length → rev.length = n →
    (∀ t, t < n → MatShape (a.getD t []) ∧ MatShape (rev.getD t [])) → MatShape acc →
    MatShape (sparsePolyConvL22 a rev acc) ∧ ∀ i j, i < 11 → j < 11 →
      sparseLMat22 (sparsePolyConvL22 a rev acc) i j = sparseLMat22 acc i j +
        ∑ t ∈ range n, sparseLMat22 (sparseMatMulL22 (a.getD t []) (rev.getD t [])) i j
  | 0, a, rev, acc, _, hr, _, hacc => by
    rw [List.eq_nil_of_length_eq_zero hr]
    have h : sparsePolyConvL22 a [] acc = acc := by cases a <;> rfl
    rw [h]
    exact ⟨hacc, fun i j _ _ => by simp⟩
  | n + 1, x :: xs, y :: ys, acc, ha, hr, hsh, hacc => by
    have h0 := hsh 0 (by omega)
    simp only [List.getD_cons_zero] at h0
    have hxy := matMul_spec x y h0.1 h0.2
    have hadd := matAdd_spec acc _ hacc hxy.1
    have ih := conv_spec n xs ys _ (by simp at ha; omega) (by simpa using hr)
      (fun t ht => by
        have := hsh (t + 1) (by omega)
        simpa only [List.getD_cons_succ] using this) hadd.1
    have hdef : sparsePolyConvL22 (x :: xs) (y :: ys) acc =
        sparsePolyConvL22 xs ys (sparseMatAddL22 acc (sparseMatMulL22 x y)) := by
      rw [sparsePolyConvL22]
    rw [hdef]
    refine ⟨ih.1, fun i j hi hj => ?_⟩
    rw [ih.2 i j hi hj, hadd.2 i j hi hj, Finset.sum_range_succ']
    simp only [List.getD_cons_succ, List.getD_cons_zero]
    ring
  | _ + 1, [], _, _, ha, _, _, _ => by simp at ha
  | _ + 1, _ :: _, [], _, _, hr, _, _ => by simp at hr

theorem go_sound (D : ℕ) (s : ℤ) (a b : List (List (List Int)))
    (ha : ∀ t, t ≤ D → MatShape (a.getD t [])) (hb : ∀ t, t ≤ D → MatShape (b.getD t []))
    (haL : D < a.length) :
    ∀ (u : List (List (List Int))) (d : ℕ) (rev bRest : List (List (List Int))),
      d + u.length ≤ D + 1 → rev.length = d →
      (∀ t, t < d → rev.getD t [] = b.getD (d - 1 - t) []) → bRest = b.drop d →
      sparsePolyMulGoDL22 s a rev bRest u = true →
      ∀ e, e < u.length → MatShape (u.getD e []) ∧ ∀ i j, i < 11 → j < 11 →
        (∑ t ∈ range (d + e + 1),
          sparseLMat22 (sparseMatMulL22 (a.getD t []) (b.getD (d + e - t) [])) i j) ≤
          s * sparseLMat22 (u.getD e []) i j
  | [], _, _, _, _, _, _, _, _ => fun e he => by simp at he
  | m :: us, d, rev, bRest, hlen, hrl, hrv, hbr, hc => by
    cases bRest with
    | nil => simp [sparsePolyMulGoDL22] at hc
    | cons bm bs =>
      have hc' : ((sparseMatShapeL22 m &&
          sparseMatLeScaledL22 s (sparsePolyConvL22 a (bm :: rev) sparseZeroMatL22) m) &&
          sparsePolyMulGoDL22 s a (bm :: rev) bs us) = true := hc
      rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
      obtain ⟨⟨hm, hle⟩, hrest⟩ := hc'
      obtain ⟨hbm, hbs⟩ := drop_cons [] b d bm bs hbr.symm
      have hlen' : d + 1 + us.length ≤ D + 1 := by simp at hlen; omega
      have hrv' : ∀ t, t < d + 1 → (bm :: rev).getD t [] = b.getD (d - t) [] := by
        intro t ht
        cases t with
        | zero => rw [List.getD_cons_zero, Nat.sub_zero, hbm]
        | succ t =>
          rw [List.getD_cons_succ, hrv t (by omega)]
          congr 1
          omega
      have IH := go_sound D s a b ha hb haL us (d + 1) (bm :: rev) bs hlen'
        (by simp [hrl]) hrv' hbs.symm hrest
      intro e he
      cases e with
      | zero =>
        have hmS := matShape_of_check m hm
        refine ⟨by simpa using hmS, fun i j hi hj => ?_⟩
        have hS := conv_spec (d + 1) a (bm :: rev) sparseZeroMatL22 (by omega)
          (by simp [hrl])
          (fun t ht => ⟨ha t (by omega), by rw [hrv' t ht]; exact hb _ (by omega)⟩)
          zeroMat_shape
        have hcmp := matLeScaled_sound s _ _ hS.1 hmS hle i j hi hj
        rw [hS.2 i j hi hj, zeroMat_get i j hi, zero_add] at hcmp
        show (∑ t ∈ range (d + 1),
          sparseLMat22 (sparseMatMulL22 (a.getD t []) (b.getD (d - t) [])) i j) ≤
          s * sparseLMat22 m i j
        calc (∑ t ∈ range (d + 1),
              sparseLMat22 (sparseMatMulL22 (a.getD t []) (b.getD (d - t) [])) i j)
            = ∑ t ∈ range (d + 1),
              sparseLMat22 (sparseMatMulL22 (a.getD t []) ((bm :: rev).getD t [])) i j :=
              Finset.sum_congr rfl fun t ht => by rw [hrv' t (Finset.mem_range.mp ht)]
          _ ≤ _ := hcmp
      | succ e =>
        have h := IH e (by simp at he; omega)
        rw [List.getD_cons_succ]
        refine ⟨h.1, fun i j hi hj => ?_⟩
        have h2 := h.2 i j hi hj
        rwa [show d + 1 + e = d + (e + 1) by omega] at h2

/-- Degree-D analog of `SparseListSound22.polyMul_sound`. -/
theorem polyMulD_sound (D : ℕ) (scale : ℤ) (hs : 0 < scale) (a b u : List (List (List Int)))
    (ha : ∀ d, d ≤ D → MatShape (a.getD d [])) (hb : ∀ d, d ≤ D → MatShape (b.getD d []))
    (hcheck : sparsePolyMulCheckDL22 D scale a b u = true)
    (F G : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hF0 : ∀ d i j, 0 ≤ F d i j) (hG0 : ∀ d i j, 0 ≤ G d i j)
    (hF : ∀ d, d ≤ D → ∀ i j : Fin 11, F d i j ≤ (sparseLPoly22 a d i j : ℝ) / scale)
    (hG : ∀ d, d ≤ D → ∀ i j : Fin 11, G d i j ≤ (sparseLPoly22 b d i j : ℝ) / scale) :
    (∀ d, d ≤ D → MatShape (u.getD d [])) ∧
      ∀ d, d ≤ D → ∀ i j : Fin 11,
        (∑ p ∈ antidiagonal d, F p.1 * G p.2) i j ≤ (sparseLPoly22 u d i j : ℝ) / scale := by
  simp only [sparsePolyMulCheckDL22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨hlen, hgo⟩ := hcheck
  have haL : D < a.length := len_of_shape a D (ha D le_rfl)
  have hGo := go_sound D scale a b ha hb haL u 0 [] b (by omega) rfl
    (fun t ht => absurd ht (Nat.not_lt_zero _)) (List.drop_zero).symm hgo
  refine ⟨fun d hd => (hGo d (by omega)).1, fun d hd i j => ?_⟩
  have hI := (hGo d (by omega)).2 i j i.isLt j.isLt
  rw [Nat.zero_add] at hI
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hpair : ∀ r t, r ≤ D → t ≤ D → (F r * G t) i j ≤
      (sparseLMat22 (sparseMatMulL22 (a.getD r []) (b.getD t [])) i j : ℝ) /
        ((scale : ℝ) * scale) := by
    intro r t hr ht
    rw [(matMul_spec _ _ (ha r hr) (hb t ht)).2 i j i.isLt j.isLt,
      Matrix.mul_apply, Int.cast_sum]
    simp only [Int.cast_mul]
    exact dot_le _ hsR _ _ _ _ (fun k => hF0 r i k) (fun k => hG0 t k j)
      (fun k => hF r hr i k) (fun k => hG t ht k j)
  rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Matrix.sum_apply]
  refine scale_close scale _ _ hs _ ?_ hI
  rw [Int.cast_sum, Finset.sum_div]
  apply Finset.sum_le_sum
  intro t ht
  have ht' := Finset.mem_range.mp ht
  exact hpair t (d - t) (by omega) (by omega)

/-! ## Squaring trace -/

abbrev PolyBoundD (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (scale : ℤ) (n : ℕ)
    (P : List (List (List Int))) (D : ℕ) : Prop :=
  (∀ d, d ≤ D → MatShape (P.getD d [])) ∧
    ∀ d, d ≤ D → ∀ i j : Fin 11,
      regionalCountDP22 B n d i j ≤ (sparseLPoly22 P d i j : ℝ) / scale

theorem add_boundD (D : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (a b u : List (List (List Int))) (m n : ℕ)
    (hcheck : sparsePolyMulCheckDL22 D scale a b u = true)
    (ha : PolyBoundD B scale m a D) (hb : PolyBoundD B scale n b D) :
    PolyBoundD B scale (m + n) u D := by
  have h := polyMulD_sound D scale hs a b u ha.1 hb.1 hcheck
    (regionalCountDP22 B m) (regionalCountDP22 B n)
    (count_nonneg B hB m) (count_nonneg B hB n) ha.2 hb.2
  refine ⟨h.1, fun d hd i j => ?_⟩
  unfold regionalCountDP22
  rw [orderedMatrixDP_add]
  exact h.2 d hd i j

theorem chainD_sound (D : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (base : List (List (List Int))) (hbase : PolyBoundD B scale 1 base D) :
    ∀ (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
      (prev : List (List (List Int))) (e : ℕ),
      PolyBoundD B scale e prev D →
      sparseRegionalPowerCheckDL22 D scale base prev trace = true →
      PolyBoundD B scale (trace.foldl sparseRegionalStepL22 e)
        (sparseRegionalLastL22 prev trace) D
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((sparsePolyMulCheckDL22 D scale prev prev sq &&
        (if bit then sparsePolyMulCheckDL22 D scale sq base next else true)) &&
        sparseRegionalPowerCheckDL22 D scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := add_boundD D B hB scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : sparseRegionalPowerCheckDL22 D scale base sq tail = true := hrest
      have he : sparseRegionalStepL22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hsq' : PolyBoundD B scale (sparseRegionalStepL22 e (false, sq, next)) sq D := by
        rw [he]; exact hsqB
      exact chainD_sound D B hB scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : sparsePolyMulCheckDL22 D scale sq base next = true := hbit
      have hrest' : sparseRegionalPowerCheckDL22 D scale base next tail = true := hrest
      have hnext := add_boundD D B hB scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStepL22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hnext' : PolyBoundD B scale (sparseRegionalStepL22 e (true, sq, next)) next D := by
        rw [he]; exact hnext
      exact chainD_sound D B hB scale hs base hbase tail next _ hnext' hrest'

/-! ## Base stage -/

theorem binom16L_eq (d : ℕ) (hd : d ≤ 16) :
    sparseBinom16L22.getD d 0 = (Nat.choose 16 d : ℤ) := by
  rw [Nat.choose_eq_factorial_div_factorial hd]
  interval_cases d <;> rfl

theorem binom16L_high (d : ℕ) (hd : 16 < d) : sparseBinom16L22.getD d 0 = 0 :=
  getD_ge 0 sparseBinom16L22 d (by simp [sparseBinom16L22]; omega)

theorem count_one_high (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (d : ℕ) (hd : 16 < d) :
    regionalCountDP22 B 1 d = 0 := by
  unfold regionalCountDP22
  simp only [orderedMatrixDP]
  apply Finset.sum_eq_zero
  intro j _
  have hj : d - j.val ≠ 0 := by have := j.isLt; omega
  simp [hj]

theorem baseGo_sound (row : ℕ → ℕ → ℕ → ℤ) (mult : ℤ) :
    ∀ (ms : List (List (List Int))) (d : ℕ), sparseRegionalBaseGoDL22 row mult d ms = true →
      ∀ e, e < ms.length → MatShape (ms.getD e []) ∧ ∀ i k, i < 11 → k < 11 →
        sparseBinom16L22.getD (d + e) 0 * row (d + e) i k * mult ≤
          sparseLMat22 (ms.getD e []) i k
  | [], _, _ => fun e he => by simp at he
  | m :: ms, d, hc => by
    have hc' : ((sparseMatShapeL22 m &&
        (List.range 11).all fun i => (List.range 11).all fun k =>
          decide (sparseBinom16L22.getD d 0 * row d i k * mult ≤ (m.getD i []).getD k 0)) &&
        sparseRegionalBaseGoDL22 row mult (d + 1) ms) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hm, hent⟩, hrest⟩ := hc'
    simp only [List.all_eq_true, List.mem_range, decide_eq_true_eq] at hent
    intro e he
    cases e with
    | zero =>
      refine ⟨by simpa using matShape_of_check m hm, fun i k hi hk => ?_⟩
      exact hent i hi k hk
    | succ e =>
      have h := baseGo_sound row mult ms (d + 1) hrest e (by simp at he; omega)
      rw [List.getD_cons_succ]
      refine ⟨h.1, fun i k hi hk => ?_⟩
      have h2 := h.2 i k hi hk
      rwa [show d + 1 + e = d + (e + 1) by omega] at h2

theorem baseD_sound (D : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j)
    (row : ℕ → ℕ → ℕ → ℤ) (den mult : ℤ) (hden : 0 < den) (hmult : 0 < mult)
    (hrow : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / den)
    (base : List (List (List Int)))
    (hcheck : sparseRegionalBaseCheckDL22 D row mult base = true) :
    PolyBoundD B (den * mult) 1 base D := by
  have _hB := hB
  simp only [sparseRegionalBaseCheckDL22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨hlen, hgo⟩ := hcheck
  have hG := baseGo_sound row mult base 0 hgo
  refine ⟨fun d hd => (hG d (by omega)).1, fun d hd i k => ?_⟩
  have hc := (hG d (by omega)).2 i k i.isLt k.isLt
  rw [Nat.zero_add] at hc
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have hmR : (0 : ℝ) < mult := by exact_mod_cast hmult
  have hdm : (0 : ℝ) < ((den * mult : ℤ) : ℝ) := by
    rw [Int.cast_mul]; exact mul_pos hdR hmR
  by_cases h16 : d ≤ 16
  · rw [binom16L_eq d h16] at hc
    have hcR : (Nat.choose 16 d : ℝ) * (row d i.val k.val : ℝ) * (mult : ℝ) ≤
        (sparseLPoly22 base d i k : ℝ) := by
      unfold sparseLPoly22; exact_mod_cast hc
    rw [count_one B d (by omega), Matrix.smul_apply, smul_eq_mul, le_div_iff₀ hdm,
      Int.cast_mul]
    have hC : (0 : ℝ) ≤ (Nat.choose 16 d : ℝ) := Nat.cast_nonneg _
    have h1 := mul_le_mul_of_nonneg_left (hrow d h16 i k) hC
    have h2 : (Nat.choose 16 d : ℝ) * ((row d i.val k.val : ℝ) / den) * ((den : ℝ) * mult) =
        (Nat.choose 16 d : ℝ) * (row d i.val k.val : ℝ) * mult := by
      field_simp
    calc (Nat.choose 16 d : ℝ) * B ⟨d, by omega⟩ i k * ((den : ℝ) * mult)
        ≤ (Nat.choose 16 d : ℝ) * ((row d i.val k.val : ℝ) / den) * ((den : ℝ) * mult) :=
          mul_le_mul_of_nonneg_right h1 (mul_pos hdR hmR).le
      _ = _ := h2
      _ ≤ _ := hcR
  · rw [binom16L_high d (by omega), zero_mul, zero_mul] at hc
    rw [count_one_high B d (by omega), Matrix.zero_apply]
    apply div_nonneg _ hdm.le
    unfold sparseLPoly22
    exact_mod_cast hc

end SparseListPolyD22

open SparseListPolyD22 in
#print axioms polyMulD_sound
open SparseListPolyD22 in
#print axioms chainD_sound
open SparseListPolyD22 in
#print axioms baseD_sound
end Spin.Turbo
