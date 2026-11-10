import SpinCodes.Turbo.PairedFracABDefs22
import SpinCodes.Turbo.PairedFracPacked22
import SpinCodes.Turbo.SparseRpowABCert22

/-! Paired fractional chain at a general exponent α = a/b (packed squaring
checks). Proofs are the green α = 4/5 originals (`PairedFracChain22`,
`PairedFracPacked22`) with the exponent and certificate generalized. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 SparseListPolyD22

namespace PairedFracChain22

theorem rpowABCheck_entry (ea eb : ℕ) (den sy : ℤ) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (h : pairedRpowABCheckL22 ea eb den sy R Ys = true)
    (a b : ℕ) (ha : a < 17) (hb : b < 17) :
    0 < den ∧ MatShape (R.getD a []) ∧ MatShape (R.getD b []) ∧
      sparseRpowABCheckL22 ea eb (den * den) sy
        (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b []) = true := by
  simp only [pairedRpowABCheckL22, pairedRpowABRowL22, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true,
    List.mem_range] at h
  obtain ⟨hden, hall⟩ := h
  exact ⟨hden, matShape_of_check _ (hall a ha).1, matShape_of_check _ (hall b hb).1,
    (hall a ha).2 b hb⟩

theorem pairedProfileAB_le (ea eb : ℕ) (heb : eb ≠ 0) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (den s : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den s R Ys = true)
    (a b : Fin 17) (i k : Fin 11) :
    pairedProfileMatrix22 B ((ea : ℝ) / eb) (a, b) i k ≤
      ((Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b.val : ℝ)) *
        (sparseLMat22 ((Ys.getD a.val []).getD b.val []) i.val k.val : ℝ) / s := by
  obtain ⟨hden, hRa, hRb, hck⟩ := rpowABCheck_entry ea eb den s R Ys hrp a.val b.val a.isLt b.isLt
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
  have hY := SparseRpowABCert22.sparseRpowABCheckL22_entrywiseRpow ea eb heb hck (B a * B b) hM0 hMX i k
  unfold pairedProfileMatrix22
  rw [Matrix.smul_apply, smul_eq_mul, mul_div_assoc]
  exact mul_le_mul_of_nonneg_left hY (by positivity)

theorem pairedChainABPk_sound (ea eb : ℕ) (heb : eb ≠ 0) (D v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (s : ℤ) (hs : 0 < s) (den : ℤ)
    (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (hchain : pkRegionalPowerCheckL22 D v w s base base trace = true) :
    PolyBoundW (pairedProfileMatrix22 B ((ea : ℝ) / eb)) s
      (trace.foldl sparseRegionalStepL22 1) (sparseRegionalLastL22 base trace) D := by
  have hW0 := pairedProfile_nonneg B hB ((ea : ℝ) / eb)
  have hb := base_sound D (pairedProfileMatrix22 B ((ea : ℝ) / eb)) s hs Ys
    (pairedProfileAB_le ea eb heb B hB den s R hR Ys hrp) base hbase
  exact chainWPk_sound D v w _ hW0 s hs base hb trace base 1 hb hchain


end PairedFracChain22

open PairedFracChain22 in
#print axioms pairedChainABPk_sound
end Spin.Turbo
