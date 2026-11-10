import SpinCodes.Turbo.TailTermFrac22
import SpinCodes.Turbo.PairedFracPacked22

/-! Packed-chain variant of `TailTermFrac22.tail_term_frac_le` (identical proof;
only the squaring-chain hypothesis is the packed check). -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22
namespace TailTermFrac22

theorem tail_term_frac_pk_le (v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (D q : ℕ) (hq : q ≤ D) (s : ℤ) (S : ℕ) (hS : s = 2 ^ S)
    -- pair certificates and degree-D chain over 64 pairs
    (den : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowCheckL22 den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (htr : trace.foldl sparseRegionalStepL22 1 = 64)
    (hchain : pkRegionalPowerCheckL22 D v w s base base trace = true)
    -- regional normalization and domination
    (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose 2048 q : ℝ)) (hCpos : 0 < Cq)
    (U0 : List (List Int)) (hU0S : sparseMatShapeL22 U0 = true)
    (hdom : sparseMatLeScaledL22 Cq ((sparseRegionalLastL22 base trace).getD q []) U0 = true)
    -- 64 regions
    (stages : List (Nat × Nat × List (List Int))) (hlen : stages.length = 6)
    (hsq : sparseSquareChainCheckL22 s U0 stages = true)
    -- prefactor
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ 80)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum 80 52428 ztrace = true)
    (vtrace : List (Bool × Nat × Int))
    (hvc : sparseScaledPowerCheck22 65535 0 (8 * q) vtrace = true) (yn G : ℤ)
    (hpre : pow45Pow2Check22
      (256 * (q : ℤ) + sparseScaledPowerExp22 80 ztrace + sparseScaledPowerExp22 0 vtrace) G
      (sparseScaledPowerMant22 ztrace * sparseScaledPowerMant22 vtrace) yn = true)
    -- final comparison
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailFracFinalCheck22 Cq yn
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j)) a G
      (SparseListSquare22.chainExp 0 stages) Kt S = true) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ 52428) ^ ((4 : ℝ) / 5) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((4 : ℝ) / 5) 64 q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
  have hs : (0 : ℤ) < s := by rw [hS]; positivity
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hCR : (0 : ℝ) < Cq := by exact_mod_cast hCpos
  -- degree-D chain on the paired DP
  have hlast := PairedFracChain22.pairedChainPk_sound D v w B hB s hs den R hR Ys hrp base hbase
    trace hchain
  rw [htr] at hlast
  set last := sparseRegionalLastL22 base trace with hlastdef
  have hW0 := PairedFracChain22.pairedProfile_nonneg B hB ((4 : ℝ) / 5)
  -- paired coefficient
  have hA : pairedRegionDP22 B ((4 : ℝ) / 5) 64 q = (1 / (Cq : ℝ)) •
      orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B ((4 : ℝ) / 5)) 64 q := by
    unfold pairedRegionDP22
    rw [show (32 * 64 : ℕ) = 2048 from rfl, hC]
  set A := pairedRegionDP22 B ((4 : ℝ) / 5) 64 q with hAdef
  have hA0 : ∀ i j, 0 ≤ A i j := by
    intro i j
    rw [hA, Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (PairedFracChain22.dp_nonneg _ hW0 64 q i j)
  have hU0 := matShape_of_check U0 hU0S
  have hlastq := hlast.1 q hq
  have hAU : ∀ i j : Fin 11, A i j ≤ (2 : ℝ) ^ (0 : ℤ) * (sparseLMat22 U0 i j : ℝ) / s := by
    intro i j
    have hX := hlast.2 q hq i j
    have hd := matLeScaled_sound Cq _ _ hlastq hU0 hdom i j i.isLt j.isLt
    have hdR : (sparseLPoly22 last q i j : ℝ) ≤ (Cq : ℝ) * sparseLMat22 U0 i j := by
      unfold sparseLPoly22; exact_mod_cast hd
    rw [hA, Matrix.smul_apply, smul_eq_mul, zpow_zero, one_mul]
    calc 1 / (Cq : ℝ) * orderedMatrixDP pairedProfileDegree22
          (pairedProfileMatrix22 B ((4 : ℝ) / 5)) 64 q i j
        ≤ 1 / (Cq : ℝ) * ((sparseLPoly22 last q i j : ℝ) / s) :=
          mul_le_mul_of_nonneg_left hX (by positivity)
      _ ≤ 1 / (Cq : ℝ) * (((Cq : ℝ) * sparseLMat22 U0 i j) / s) :=
          mul_le_mul_of_nonneg_left (div_le_div_of_nonneg_right hdR hsR.le) (by positivity)
      _ = (sparseLMat22 U0 i j : ℝ) / s := by field_simp
  -- 64 regions
  have hpow := SparseListSquare22.chain_sound s hs stages U0 0 A hA0 hU0 hAU hsq
  rw [hlen] at hpow
  set E := SparseListSquare22.chainExp 0 stages
  set Ul := sparseSquareChainLastL22 U0 stages
  set sumU := sparseSum11Int22 (fun j => sparseLMat22 Ul 0 j)
  have hmass : (∑ k : Fin 11, (A ^ 64) 0 k) ≤ (2 : ℝ) ^ E * (sumU : ℝ) / s := by
    have hsum : (sumU : ℝ) = ∑ k : Fin 11, (sparseLMat22 Ul 0 k : ℝ) := sum11_cast _
    rw [hsum, Finset.mul_sum, Finset.sum_div]
    exact Finset.sum_le_sum fun k _ => by simpa [mul_div_assoc] using hpow 0 k
  have hmass0 : 0 ≤ ∑ k : Fin 11, (A ^ 64) 0 k :=
    Finset.sum_nonneg fun k _ => markerPowNonneg22 A hA0 64 0 k
  -- prefactor
  obtain ⟨hX0, hXp⟩ := prefactor_le q z zNum hz ztrace hzc vtrace hvc yn G hpre
  set X := turboOldUniformFactor22 ^ q / z ^ 52428
  have hY0 : 0 ≤ (yn : ℝ) * (2 : ℝ) ^ G := (Real.rpow_nonneg hX0 _).trans hXp
  -- assemble
  have hfin := final_real Cq yn sumU a G E Kt S hfinal
  rw [← hC]
  have h2G : (2 : ℝ) ^ G ≠ 0 := zpow_ne_zero _ two_ne_zero
  have h2E : (2 : ℝ) ^ E ≠ 0 := zpow_ne_zero _ two_ne_zero
  calc (Cq : ℝ) * (X ^ ((4 : ℝ) / 5) * ∑ k : Fin 11, (A ^ 64) 0 k)
      ≤ (Cq : ℝ) * (((yn : ℝ) * (2 : ℝ) ^ G) * ((2 : ℝ) ^ E * (sumU : ℝ) / s)) :=
        mul_le_mul_of_nonneg_left (mul_le_mul hXp hmass hmass0 hY0) hCR.le
    _ = ((Cq : ℝ) * yn * sumU * (2 : ℝ) ^ (G + E + Kt - S)) / 2 ^ Kt := by
        rw [hS]
        push_cast
        have e1 : (2 : ℝ) ^ (G + E + Kt - S) =
            (2 : ℝ) ^ G * (2 : ℝ) ^ E * (2 : ℝ) ^ Kt / (2 : ℝ) ^ S := by
          rw [zpow_sub₀ (by norm_num), zpow_add₀ (by norm_num), zpow_add₀ (by norm_num),
            zpow_natCast, zpow_natCast]
        rw [e1]
        field_simp
    _ ≤ (a : ℝ) / 2 ^ Kt := div_le_div_of_nonneg_right hfin (by positivity)

end TailTermFrac22
#print axioms TailTermFrac22.tail_term_frac_pk_le
end Spin.Turbo
