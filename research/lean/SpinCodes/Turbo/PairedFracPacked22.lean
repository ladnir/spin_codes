import SpinCodes.Turbo.PairedFracChain22
import SpinCodes.Turbo.SparseListPacked22

/-! Packed-check variants of the paired fractional chain (`PairedFracChain22`).
Proofs are the green originals with `sparsePolyMulCheckDL22` replaced by the
packed `pkPolyMulCheckL22` (soundness `SparseListPacked22.polyMulPk_sound`). -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 SparseListPolyD22

namespace PairedFracChain22

theorem add_boundWPk (D v w : ℕ) (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hW : ∀ ab i j, 0 ≤ W ab i j) (scale : ℤ) (hs : 0 < scale)
    (a b u : List (List (List Int))) (m n : ℕ)
    (hcheck : pkPolyMulCheckL22 D v w scale a b u = true)
    (ha : PolyBoundW W scale m a D) (hb : PolyBoundW W scale n b D) :
    PolyBoundW W scale (m + n) u D := by
  have h := SparseListPacked22.polyMulPk_sound D v w scale hs a b u ha.1 hb.1 hcheck
    (orderedMatrixDP pairedProfileDegree22 W m) (orderedMatrixDP pairedProfileDegree22 W n)
    (fun d => dp_nonneg W hW m d) (fun d => dp_nonneg W hW n d) ha.2 hb.2
  refine ⟨h.1, fun d hd i j => ?_⟩
  rw [orderedMatrixDP_add]
  exact h.2 d hd i j

theorem chainWPk_sound (D v w : ℕ) (W : Fin 17 × Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hW : ∀ ab i j, 0 ≤ W ab i j) (scale : ℤ) (hs : 0 < scale)
    (base : List (List (List Int))) (hbase : PolyBoundW W scale 1 base D) :
    ∀ (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
      (prev : List (List (List Int))) (e : ℕ),
      PolyBoundW W scale e prev D →
      pkRegionalPowerCheckL22 D v w scale base prev trace = true →
      PolyBoundW W scale (trace.foldl sparseRegionalStepL22 e)
        (sparseRegionalLastL22 prev trace) D
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((pkPolyMulCheckL22 D v w scale prev prev sq &&
        (if bit then pkPolyMulCheckL22 D v w scale sq base next else true)) &&
        pkRegionalPowerCheckL22 D v w scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := add_boundWPk D v w W hW scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : pkRegionalPowerCheckL22 D v w scale base sq tail = true := hrest
      have he : sparseRegionalStepL22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hsq' : PolyBoundW W scale (sparseRegionalStepL22 e (false, sq, next)) sq D := by
        rw [he]; exact hsqB
      exact chainWPk_sound D v w W hW scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : pkPolyMulCheckL22 D v w scale sq base next = true := hbit
      have hrest' : pkRegionalPowerCheckL22 D v w scale base next tail = true := hrest
      have hnext := add_boundWPk D v w W hW scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStepL22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hnext' : PolyBoundW W scale (sparseRegionalStepL22 e (true, sq, next)) next D := by
        rw [he]; exact hnext
      exact chainWPk_sound D v w W hW scale hs base hbase tail next _ hnext' hrest'

/-! ## Base stage -/

theorem pairedChainPk_sound (D v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (s : ℤ) (hs : 0 < s) (den : ℤ)
    (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowCheckL22 den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (hchain : pkRegionalPowerCheckL22 D v w s base base trace = true) :
    PolyBoundW (pairedProfileMatrix22 B ((4 : ℝ) / 5)) s
      (trace.foldl sparseRegionalStepL22 1) (sparseRegionalLastL22 base trace) D := by
  have hW0 := pairedProfile_nonneg B hB ((4 : ℝ) / 5)
  have hb := base_sound D (pairedProfileMatrix22 B ((4 : ℝ) / 5)) s hs Ys
    (pairedProfile_le B hB den s R hR Ys hrp) base hbase
  exact chainWPk_sound D v w _ hW0 s hs base hb trace base 1 hb hchain


theorem rpowCheck_of (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (h0 : 0 < den)
    (h : ∀ a, a < 17 → (sparseMatShapeL22 (R.getD a []) && (List.range 17).all fun b =>
      sparseRpow45CheckL22 (den * den) sy
        (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b [])) = true) :
    pairedRpowCheckL22 den sy R Ys = true := by
  unfold pairedRpowCheckL22
  rw [decide_eq_true h0, Bool.true_and, List.all_eq_true]
  intro a ha
  exact h a (List.mem_range.mp ha)

theorem baseGo_append (Ys : List (List (List (List Int)))) :
    ∀ (l1 l2 : List (List (List Int))) (d : Nat),
      pairedBaseGoL22 Ys d (l1 ++ l2) =
        (pairedBaseGoL22 Ys d l1 && pairedBaseGoL22 Ys (d + l1.length) l2)
  | [], l2, d => by simp [pairedBaseGoL22]
  | m :: ms, l2, d => by
    rw [List.cons_append, pairedBaseGoL22, pairedBaseGoL22, baseGo_append Ys ms l2 (d + 1),
      List.length_cons, show d + (ms.length + 1) = d + 1 + ms.length by omega]
    simp only [Bool.and_assoc]

end PairedFracChain22

open PairedFracChain22 in
#print axioms pairedChainPk_sound
end Spin.Turbo
