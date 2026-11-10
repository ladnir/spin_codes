import SpinCodes.Turbo.SparseListPolyD22
import SpinCodes.Turbo.SparseListSquare22
import SpinCodes.Turbo.SparseScaledPower22
import SpinCodes.Turbo.PairedSquare22
import SpinCodes.Turbo.RS8UniformReference22

/-! One exact-coefficient tail occupation term (α = 1, K18: 64 pairs, 64
regions, d = 52428), bounded by integer certificates:
degree-D chain bound on the single-step power A^128 → paired coefficient
(`pairedRegionDP22_one_eq`) → domination A ≤ U0/s → scaled squaring A^64 →
scaled z^d lower bound → one integer comparison with all powers of two
folded into a single exponent. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22

/-- Final integer comparison; T collects every power of two. -/
def tailFinalCheck22 (C sumU a : Int) (q : Nat) (Ez E : Int) (Kt S : Nat) (m : Int) : Bool :=
  if 0 ≤ (256 * (q : Int) + Ez + E + Kt - S) then
    decide (C * sumU * 2 ^ (256 * (q : Int) + Ez + E + Kt - S).toNat ≤ a * 65535 ^ (8 * q) * m)
  else
    decide (C * sumU ≤ a * 65535 ^ (8 * q) * m * 2 ^ (-(256 * (q : Int) + Ez + E + Kt - S)).toNat)

namespace TailTerm22

theorem final_real (C sumU a : Int) (q : Nat) (Ez E : Int) (Kt S : Nat) (m : Int)
    (h : tailFinalCheck22 C sumU a q Ez E Kt S m = true) :
    (C : ℝ) * sumU * (2 : ℝ) ^ (256 * (q : Int) + Ez + E + Kt - S) ≤
      (a : ℝ) * (65535 : ℝ) ^ (8 * q) * m := by
  unfold tailFinalCheck22 at h
  set T : Int := 256 * (q : Int) + Ez + E + Kt - S with hT
  split_ifs at h with hT0
  · have h' := of_decide_eq_true h
    have hR : ((C * sumU * 2 ^ T.toNat : ℤ) : ℝ) ≤ ((a * 65535 ^ (8 * q) * m : ℤ) : ℝ) := by
      exact_mod_cast h'
    push_cast at hR
    have hz : (2 : ℝ) ^ T = (2 : ℝ) ^ T.toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg hT0]
    rw [hz]; exact hR
  · have h' := of_decide_eq_true h
    have hR : ((C * sumU : ℤ) : ℝ) ≤ ((a * 65535 ^ (8 * q) * m * 2 ^ (-T).toNat : ℤ) : ℝ) := by
      exact_mod_cast h'
    push_cast at hR
    have hpos : (0 : ℝ) < (2 : ℝ) ^ (-T).toNat := by positivity
    have hz : (2 : ℝ) ^ T = 1 / (2 : ℝ) ^ (-T).toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg (by omega), _root_.zpow_neg, one_div, inv_inv]
    rw [hz, mul_one_div, div_le_iff₀ hpos]
    exact hR

theorem tail_term_le (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (D q : ℕ) (hq : q ≤ D) (s : ℤ) (S : ℕ) (hS : s = 2 ^ S)
    (last : List (List (List Int))) (hlast : SparseListPolyD22.PolyBoundD B s 128 last D)
    (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose 2048 q : ℝ)) (hCpos : 0 < Cq)
    (U0 : List (List Int)) (hU0S : sparseMatShapeL22 U0 = true) (down0 up0 : ℕ)
    (hdom : sparseMatLeScaledL22 (Cq * 2 ^ up0) (sparseMatScaleL22 (2 ^ down0) (last.getD q [])) U0 = true)
    (stages : List (Nat × Nat × List (List Int))) (hlen : stages.length = 6)
    (hchain : sparseSquareChainCheckL22 s U0 stages = true)
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ 80)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum 80 52428 ztrace = true)
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailFinalCheck22 Cq
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j)) a q
      (sparseScaledPowerExp22 80 ztrace) (SparseListSquare22.chainExp ((up0 : ℤ) - down0) stages) Kt S
      (sparseScaledPowerMant22 ztrace) = true) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
  have hs : (0 : ℤ) < s := by rw [hS]; positivity
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hCR : (0 : ℝ) < Cq := by exact_mod_cast hCpos
  -- paired coefficient
  have hA : pairedRegionDP22 B 1 64 q = (1 / (Cq : ℝ)) • regionalCountDP22 B 128 q := by
    have h := PairedSquare22.pairedRegionDP22_one_eq B 64 q
    rw [show (32 * 64 : ℕ) = 2048 from rfl, show (2 * 64 : ℕ) = 128 from rfl] at h
    rw [h, hC]
  set A := pairedRegionDP22 B 1 64 q with hAdef
  have hA0 : ∀ i j, 0 ≤ A i j := by
    intro i j
    rw [hA, Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (SparseRegionalPower22.count_nonneg B hB 128 q i j)
  have hU0 := matShape_of_check U0 hU0S
  have hlastq := hlast.1 q hq
  have hAU : ∀ i j : Fin 11,
      A i j ≤ (2 : ℝ) ^ ((up0 : ℤ) - down0) * (sparseLMat22 U0 i j : ℝ) / s := by
    intro i j
    have hX := hlast.2 q hq i j
    have hsc := SparseListSquare22.matScale_spec (2 ^ down0) _ hlastq
    have hd := matLeScaled_sound _ _ _ hsc.1 hU0 hdom i j i.isLt j.isLt
    rw [hsc.2 i j i.isLt j.isLt] at hd
    have hdR : (2 : ℝ) ^ down0 * (sparseLPoly22 last q i j : ℝ) ≤
        (Cq : ℝ) * 2 ^ up0 * sparseLMat22 U0 i j := by
      unfold sparseLPoly22; exact_mod_cast hd
    have h2d : (0 : ℝ) < (2 : ℝ) ^ down0 := by positivity
    have hz : (2 : ℝ) ^ ((up0 : ℤ) - down0) = (2 : ℝ) ^ up0 / (2 : ℝ) ^ down0 := by
      rw [zpow_sub₀ (by norm_num), zpow_natCast, zpow_natCast]
    rw [hA, Matrix.smul_apply, smul_eq_mul, hz]
    calc 1 / (Cq : ℝ) * regionalCountDP22 B 128 q i j
        ≤ 1 / (Cq : ℝ) * ((sparseLPoly22 last q i j : ℝ) / s) :=
          mul_le_mul_of_nonneg_left hX (by positivity)
      _ ≤ 1 / (Cq : ℝ) * (((Cq : ℝ) * 2 ^ up0 * sparseLMat22 U0 i j / 2 ^ down0) / s) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply div_le_div_of_nonneg_right _ hsR.le
          rw [le_div_iff₀ h2d]; linarith
      _ = (2 : ℝ) ^ up0 / (2 : ℝ) ^ down0 * (sparseLMat22 U0 i j : ℝ) / s := by field_simp
  -- 64 regions
  have hpow := SparseListSquare22.chain_sound s hs stages U0 ((up0 : ℤ) - down0) A hA0 hU0 hAU hchain
  rw [hlen] at hpow
  set E := SparseListSquare22.chainExp ((up0 : ℤ) - down0) stages
  set Ul := sparseSquareChainLastL22 U0 stages
  set sumU := sparseSum11Int22 (fun j => sparseLMat22 Ul 0 j)
  have hmass : (∑ k : Fin 11, (A ^ 64) 0 k) ≤ (2 : ℝ) ^ E * (sumU : ℝ) / s := by
    have hsum : (sumU : ℝ) = ∑ k : Fin 11, (sparseLMat22 Ul 0 k : ℝ) := sum11_cast _
    rw [hsum, Finset.mul_sum, Finset.sum_div]
    exact Finset.sum_le_sum fun k _ => by simpa [mul_div_assoc] using hpow 0 k
  have hmass0 : 0 ≤ ∑ k : Fin 11, (A ^ 64) 0 k :=
    Finset.sum_nonneg fun k _ => markerPowNonneg22 A hA0 64 0 k
  -- z^d
  obtain ⟨hm, hzlow⟩ := SparseScaledPower22.scaledPower_sound zNum 80 52428 ztrace hzc
  set m := sparseScaledPowerMant22 ztrace
  set Ez := sparseScaledPowerExp22 80 ztrace
  have hmR : (0 : ℝ) < m := by exact_mod_cast hm
  have hzd : (m : ℝ) * (2 : ℝ) ^ (-Ez) ≤ z ^ 52428 := by rw [hz]; exact hzlow
  have hzdpos : 0 < z ^ 52428 := lt_of_lt_of_le (by positivity) hzd
  -- assemble
  have hfin := final_real Cq sumU a q Ez E Kt S m hfinal
  have hbeta : turboOldUniformFactor22 ^ q = (2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) := by
    unfold turboOldUniformFactor22; rw [div_pow, ← pow_mul, ← pow_mul]
  rw [Real.rpow_one, ← hC]
  have h65 : (0 : ℝ) < (65535 : ℝ) ^ (8 * q) := by positivity
  have hnum0 : (0 : ℝ) ≤ (2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) :=
    div_nonneg (pow_nonneg (by norm_num) _) h65.le
  have hbeta0 : 0 ≤ turboOldUniformFactor22 ^ q / z ^ 52428 := by
    rw [hbeta]; exact div_nonneg hnum0 hzdpos.le
  have hmz : (0 : ℝ) < (m : ℝ) * (2 : ℝ) ^ (-Ez) := mul_pos hmR (zpow_pos (by norm_num) _)
  calc (Cq : ℝ) * (turboOldUniformFactor22 ^ q / z ^ 52428 * ∑ k : Fin 11, (A ^ 64) 0 k)
      ≤ (Cq : ℝ) * ((2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) / ((m : ℝ) * (2 : ℝ) ^ (-Ez)) *
          ((2 : ℝ) ^ E * (sumU : ℝ) / s)) := by
        apply mul_le_mul_of_nonneg_left _ hCR.le
        apply mul_le_mul _ hmass hmass0 (div_nonneg hnum0 hmz.le)
        rw [hbeta]
        exact div_le_div_of_nonneg_left hnum0 hmz hzd
    _ = ((Cq : ℝ) * sumU * (2 : ℝ) ^ (256 * (q : Int) + Ez + E + Kt - S)) /
          ((65535 : ℝ) ^ (8 * q) * m * 2 ^ Kt) := by
        rw [hS]
        push_cast
        have e1 : (2 : ℝ) ^ (256 * (q : Int) + Ez + E + Kt - S) =
            (2 : ℝ) ^ (256 * q) * (2 : ℝ) ^ Ez * (2 : ℝ) ^ E * (2 : ℝ) ^ Kt / (2 : ℝ) ^ S := by
          rw [show (256 * (q : ℤ)) = ((256 * q : ℕ) : ℤ) by push_cast; ring,
            zpow_sub₀ two_ne_zero, zpow_add₀ two_ne_zero, zpow_add₀ two_ne_zero,
            zpow_add₀ two_ne_zero]
          simp only [zpow_natCast]
        rw [e1, _root_.zpow_neg]
        have h2E : (0 : ℝ) < (2 : ℝ) ^ Ez := zpow_pos (by norm_num) _
        field_simp
    _ ≤ (a : ℝ) / 2 ^ Kt := by
        rw [div_le_div_iff₀ (by positivity) (by positivity)]
        have := mul_le_mul_of_nonneg_right hfin (show (0 : ℝ) ≤ 2 ^ Kt by positivity)
        linarith

end TailTerm22

#print axioms TailTerm22.tail_term_le
end Spin.Turbo
