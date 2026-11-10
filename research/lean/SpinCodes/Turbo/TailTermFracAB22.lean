import SpinCodes.Turbo.TailTermFrac22
import SpinCodes.Turbo.PairedFracABChain22

/-! One exact-coefficient tail occupation term at a general fractional
exponent α = a/b and general geometry (`n` two-step pairs per region, 64
regions, `32 n` groups, cutoff `d`, `z = zNum / 2^zb`), bounded by integer
certificates with the packed squaring chain. Proofs are the green
`TailTermFrac22` (α = 4/5, K18) with exponent and geometry generalized. -/
namespace Spin.Turbo

/-- `0 ≤ yn` and `2^(a T - b G) ≤ yn^b · M^a` (sign-split on the exponent). -/
def powABPow2Check22 (ea eb : Nat) (T G M yn : Int) : Bool :=
  decide (0 ≤ yn) &&
    (if 0 ≤ (ea : Int) * T - (eb : Int) * G then
      decide (2 ^ ((ea : Int) * T - (eb : Int) * G).toNat ≤ yn ^ eb * M ^ ea)
    else
      decide (1 ≤ yn ^ eb * M ^ ea * 2 ^ (-((ea : Int) * T - (eb : Int) * G)).toNat))

end Spin.Turbo

noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22

namespace TailTermFrac22

theorem powAB_pow2_le (ea eb : ℕ) (heb : eb ≠ 0) (T G M yn : ℤ) (hM : 0 < M)
    (h : powABPow2Check22 ea eb T G M yn = true) :
    ((2 : ℝ) ^ T / M) ^ ((ea : ℝ) / eb) ≤ (yn : ℝ) * (2 : ℝ) ^ G := by
  unfold powABPow2Check22 at h
  rw [Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨hyn, h⟩ := h
  have hMR : (0 : ℝ) < M := by exact_mod_cast hM
  have hY0 : 0 ≤ (2 : ℝ) ^ T / M := div_nonneg (zpow_nonneg (by norm_num) _) hMR.le
  have hR0 : 0 ≤ (yn : ℝ) * (2 : ℝ) ^ G :=
    mul_nonneg (by exact_mod_cast hyn) (zpow_nonneg (by norm_num) _)
  apply SparseRpowABCert22.rpowAB_le_of_pow ea eb heb hY0 hR0
  set e : ℤ := (ea : ℤ) * T - (eb : ℤ) * G with he_def
  have key : (2 : ℝ) ^ e ≤ (yn : ℝ) ^ eb * (M : ℝ) ^ ea := by
    split_ifs at h with he
    · have h' := of_decide_eq_true h
      have hR : ((2 ^ e.toNat : ℤ) : ℝ) ≤ ((yn ^ eb * M ^ ea : ℤ) : ℝ) := by
        exact_mod_cast h'
      push_cast at hR
      have hz : (2 : ℝ) ^ e = (2 : ℝ) ^ e.toNat := by
        rw [← zpow_natCast, Int.toNat_of_nonneg he]
      rw [hz]; exact hR
    · have h' := of_decide_eq_true h
      have hR : ((1 : ℤ) : ℝ) ≤ ((yn ^ eb * M ^ ea * 2 ^ (-e).toNat : ℤ) : ℝ) := by
        exact_mod_cast h'
      push_cast at hR
      have hpos : (0 : ℝ) < (2 : ℝ) ^ (-e).toNat := by positivity
      have hz : (2 : ℝ) ^ e = 1 / (2 : ℝ) ^ (-e).toNat := by
        rw [← zpow_natCast, Int.toNat_of_nonneg (by omega), _root_.zpow_neg, one_div, inv_inv]
      rw [hz, div_le_iff₀ hpos]; exact hR
  have e1 : ((2 : ℝ) ^ T) ^ ea = (2 : ℝ) ^ e * (2 : ℝ) ^ ((eb : ℤ) * G) := by
    rw [← zpow_add₀ two_ne_zero, he_def, sub_add_cancel, ← zpow_natCast, ← _root_.zpow_mul,
      mul_comm]
  have e2 : ((yn : ℝ) * (2 : ℝ) ^ G) ^ eb = (yn : ℝ) ^ eb * (2 : ℝ) ^ ((eb : ℤ) * G) := by
    rw [mul_pow, ← zpow_natCast ((2 : ℝ) ^ G) eb, ← _root_.zpow_mul, mul_comm G]
  rw [div_pow, e1, e2, div_le_iff₀ (pow_pos hMR ea)]
  have h5 := mul_le_mul_of_nonneg_right key
    (zpow_nonneg (by norm_num : (0 : ℝ) ≤ 2) ((eb : ℤ) * G))
  calc (2 : ℝ) ^ e * (2 : ℝ) ^ ((eb : ℤ) * G)
      ≤ (yn : ℝ) ^ eb * (M : ℝ) ^ ea * (2 : ℝ) ^ ((eb : ℤ) * G) := h5
    _ = (yn : ℝ) ^ eb * (2 : ℝ) ^ ((eb : ℤ) * G) * (M : ℝ) ^ ea := by ring

theorem prefactorAB_le (ea eb : ℕ) (heb : eb ≠ 0) (zb d : ℕ) (q : ℕ) (z : ℝ) (zNum : ℤ)
    (hz : z = (zNum : ℝ) / 2 ^ zb)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum zb d ztrace = true)
    (vtrace : List (Bool × Nat × Int))
    (hvc : sparseScaledPowerCheck22 65535 0 (8 * q) vtrace = true) (yn G : ℤ)
    (hpre : powABPow2Check22 ea eb
      (256 * (q : ℤ) + sparseScaledPowerExp22 zb ztrace + sparseScaledPowerExp22 0 vtrace) G
      (sparseScaledPowerMant22 ztrace * sparseScaledPowerMant22 vtrace) yn = true) :
    0 ≤ turboOldUniformFactor22 ^ q / z ^ d ∧
      (turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) ≤ (yn : ℝ) * (2 : ℝ) ^ G := by
  obtain ⟨hm, hzlow⟩ := SparseScaledPower22.scaledPower_sound zNum zb d ztrace hzc
  obtain ⟨hmv, hvlow⟩ := SparseScaledPower22.scaledPower_sound 65535 0 (8 * q) vtrace hvc
  set m := sparseScaledPowerMant22 ztrace
  set Ez := sparseScaledPowerExp22 zb ztrace
  set mv := sparseScaledPowerMant22 vtrace
  set Ev := sparseScaledPowerExp22 0 vtrace
  have hmR : (0 : ℝ) < m := by exact_mod_cast hm
  have hmvR : (0 : ℝ) < mv := by exact_mod_cast hmv
  have hzd : (m : ℝ) * (2 : ℝ) ^ (-Ez) ≤ z ^ d := by rw [hz]; exact hzlow
  have hv : (mv : ℝ) * (2 : ℝ) ^ (-Ev) ≤ (65535 : ℝ) ^ (8 * q) := by
    rw [show ((65535 : ℤ) : ℝ) / (2 : ℝ) ^ (0 : ℕ) = (65535 : ℝ) by norm_num] at hvlow
    exact hvlow
  have hzlo : 0 < (m : ℝ) * (2 : ℝ) ^ (-Ez) := mul_pos hmR (zpow_pos (by norm_num) _)
  have hvlo : 0 < (mv : ℝ) * (2 : ℝ) ^ (-Ev) := mul_pos hmvR (zpow_pos (by norm_num) _)
  have hbeta : turboOldUniformFactor22 ^ q = (2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) := by
    unfold turboOldUniformFactor22; rw [div_pow, ← pow_mul, ← pow_mul]
  have hX0 : 0 ≤ turboOldUniformFactor22 ^ q / z ^ d := by
    rw [hbeta]
    exact div_nonneg (div_nonneg (pow_nonneg (by norm_num) _) (pow_nonneg (by norm_num) _))
      (hzlo.le.trans hzd)
  have hX : turboOldUniformFactor22 ^ q / z ^ d ≤
      (2 : ℝ) ^ (256 * (q : ℤ) + Ez + Ev) / ((m : ℝ) * mv) := by
    rw [hbeta, div_div]
    calc (2 : ℝ) ^ (256 * q) / ((65535 : ℝ) ^ (8 * q) * z ^ d)
        ≤ (2 : ℝ) ^ (256 * q) / (((mv : ℝ) * (2 : ℝ) ^ (-Ev)) * ((m : ℝ) * (2 : ℝ) ^ (-Ez))) :=
          div_le_div_of_nonneg_left (by positivity) (mul_pos hvlo hzlo)
            (mul_le_mul hv hzd hzlo.le (hvlo.le.trans hv))
      _ = (2 : ℝ) ^ (256 * (q : ℤ) + Ez + Ev) / ((m : ℝ) * mv) := by
          have e1 : (2 : ℝ) ^ (256 * (q : ℤ) + Ez + Ev) =
              (2 : ℝ) ^ (256 * q) * (2 : ℝ) ^ Ez * (2 : ℝ) ^ Ev := by
            rw [zpow_add₀ two_ne_zero, zpow_add₀ two_ne_zero, ← zpow_natCast,
              show ((256 * q : ℕ) : ℤ) = 256 * (q : ℤ) by push_cast; ring]
          rw [e1, _root_.zpow_neg, _root_.zpow_neg]
          have h2a : (2 : ℝ) ^ Ez ≠ 0 := zpow_ne_zero _ two_ne_zero
          have h2b : (2 : ℝ) ^ Ev ≠ 0 := zpow_ne_zero _ two_ne_zero
          field_simp
  have hpow := powAB_pow2_le ea eb heb (256 * (q : ℤ) + Ez + Ev) G (m * mv) yn (mul_pos hm hmv) hpre
  rw [Int.cast_mul] at hpow
  exact ⟨hX0, (SparseRpowABCert22.rpowAB_mono ea eb hX0 hX).trans hpow⟩

theorem tail_term_fracAB_le (ea eb : ℕ) (heb : eb ≠ 0) (n d zb v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (D q : ℕ) (hq : q ≤ D) (s : ℤ) (S : ℕ) (hS : s = 2 ^ S)
    -- pair certificates and degree-D chain over 64 pairs
    (den : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (htr : trace.foldl sparseRegionalStepL22 1 = n)
    (hchain : pkRegionalPowerCheckL22 D v w s base base trace = true)
    -- regional normalization and domination
    (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose (32 * n) q : ℝ)) (hCpos : 0 < Cq)
    (U0 : List (List Int)) (hU0S : sparseMatShapeL22 U0 = true)
    (hdom : sparseMatLeScaledL22 Cq ((sparseRegionalLastL22 base trace).getD q []) U0 = true)
    -- 64 regions
    (stages : List (Nat × Nat × List (List Int))) (hlen : stages.length = 6)
    (hsq : sparseSquareChainCheckL22 s U0 stages = true)
    -- prefactor
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ zb)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum zb d ztrace = true)
    (vtrace : List (Bool × Nat × Int))
    (hvc : sparseScaledPowerCheck22 65535 0 (8 * q) vtrace = true) (yn G : ℤ)
    (hpre : powABPow2Check22 ea eb
      (256 * (q : ℤ) + sparseScaledPowerExp22 zb ztrace + sparseScaledPowerExp22 0 vtrace) G
      (sparseScaledPowerMant22 ztrace * sparseScaledPowerMant22 vtrace) yn = true)
    -- final comparison
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailFracFinalCheck22 Cq yn
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j)) a G
      (SparseListSquare22.chainExp 0 stages) Kt S = true) :
    (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
  have hs : (0 : ℤ) < s := by rw [hS]; positivity
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hCR : (0 : ℝ) < Cq := by exact_mod_cast hCpos
  -- degree-D chain on the paired DP
  have hlast := PairedFracChain22.pairedChainABPk_sound ea eb heb D v w B hB s hs den R hR Ys hrp base hbase
    trace hchain
  rw [htr] at hlast
  set last := sparseRegionalLastL22 base trace with hlastdef
  have hW0 := PairedFracChain22.pairedProfile_nonneg B hB ((ea : ℝ) / eb)
  -- paired coefficient
  have hA : pairedRegionDP22 B ((ea : ℝ) / eb) n q = (1 / (Cq : ℝ)) •
      orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B ((ea : ℝ) / eb)) n q := by
    unfold pairedRegionDP22
    rw [hC]
  set A := pairedRegionDP22 B ((ea : ℝ) / eb) n q with hAdef
  have hA0 : ∀ i j, 0 ≤ A i j := by
    intro i j
    rw [hA, Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (PairedFracChain22.dp_nonneg _ hW0 n q i j)
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
          (pairedProfileMatrix22 B ((ea : ℝ) / eb)) n q i j
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
  obtain ⟨hX0, hXp⟩ := prefactorAB_le ea eb heb zb d q z zNum hz ztrace hzc vtrace hvc yn G hpre
  set X := turboOldUniformFactor22 ^ q / z ^ d
  have hY0 : 0 ≤ (yn : ℝ) * (2 : ℝ) ^ G := (Real.rpow_nonneg hX0 _).trans hXp
  -- assemble
  have hfin := final_real Cq yn sumU a G E Kt S hfinal
  rw [← hC]
  have h2G : (2 : ℝ) ^ G ≠ 0 := zpow_ne_zero _ two_ne_zero
  have h2E : (2 : ℝ) ^ E ≠ 0 := zpow_ne_zero _ two_ne_zero
  calc (Cq : ℝ) * (X ^ ((ea : ℝ) / eb) * ∑ k : Fin 11, (A ^ 64) 0 k)
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


theorem tail_term_fracAB_shift_le (ea eb : ℕ) (heb : eb ≠ 0) (n d zb v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (D q : ℕ) (hq : q ≤ D) (s : ℤ) (S : ℕ) (hS : s = 2 ^ S)
    -- pair certificates and degree-D chain over 64 pairs
    (den : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den s R Ys = true)
    (base : List (List (List Int))) (hbase : pairedBaseCheckDL22 D Ys base = true)
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (htr : trace.foldl sparseRegionalStepL22 1 = n)
    (hchain : pkRegionalPowerCheckL22 D v w s base base trace = true)
    -- regional normalization and domination
    (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose (32 * n) q : ℝ)) (hCpos : 0 < Cq)
    (U0 : List (List Int)) (hU0S : sparseMatShapeL22 U0 = true)
    (down0 up0 : ℕ)
    (hdom : sparseMatLeScaledL22 (Cq * 2 ^ up0)
      (sparseMatScaleL22 (2 ^ down0) ((sparseRegionalLastL22 base trace).getD q [])) U0 = true)
    -- 64 regions
    (stages : List (Nat × Nat × List (List Int))) (hlen : stages.length = 6)
    (hsq : sparseSquareChainCheckL22 s U0 stages = true)
    -- prefactor
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ zb)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum zb d ztrace = true)
    (vtrace : List (Bool × Nat × Int))
    (hvc : sparseScaledPowerCheck22 65535 0 (8 * q) vtrace = true) (yn G : ℤ)
    (hpre : powABPow2Check22 ea eb
      (256 * (q : ℤ) + sparseScaledPowerExp22 zb ztrace + sparseScaledPowerExp22 0 vtrace) G
      (sparseScaledPowerMant22 ztrace * sparseScaledPowerMant22 vtrace) yn = true)
    -- final comparison
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailFracFinalCheck22 Cq yn
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j)) a G
      (SparseListSquare22.chainExp ((up0 : ℤ) - down0) stages) Kt S = true) :
    (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
  have hs : (0 : ℤ) < s := by rw [hS]; positivity
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hCR : (0 : ℝ) < Cq := by exact_mod_cast hCpos
  -- degree-D chain on the paired DP
  have hlast := PairedFracChain22.pairedChainABPk_sound ea eb heb D v w B hB s hs den R hR Ys hrp base hbase
    trace hchain
  rw [htr] at hlast
  set last := sparseRegionalLastL22 base trace with hlastdef
  have hW0 := PairedFracChain22.pairedProfile_nonneg B hB ((ea : ℝ) / eb)
  -- paired coefficient
  have hA : pairedRegionDP22 B ((ea : ℝ) / eb) n q = (1 / (Cq : ℝ)) •
      orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B ((ea : ℝ) / eb)) n q := by
    unfold pairedRegionDP22
    rw [hC]
  set A := pairedRegionDP22 B ((ea : ℝ) / eb) n q with hAdef
  have hA0 : ∀ i j, 0 ≤ A i j := by
    intro i j
    rw [hA, Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (PairedFracChain22.dp_nonneg _ hW0 n q i j)
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
    calc 1 / (Cq : ℝ) * orderedMatrixDP pairedProfileDegree22
          (pairedProfileMatrix22 B ((ea : ℝ) / eb)) n q i j
        ≤ 1 / (Cq : ℝ) * ((sparseLPoly22 last q i j : ℝ) / s) :=
          mul_le_mul_of_nonneg_left hX (by positivity)
      _ ≤ 1 / (Cq : ℝ) * (((Cq : ℝ) * 2 ^ up0 * sparseLMat22 U0 i j / 2 ^ down0) / s) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply div_le_div_of_nonneg_right _ hsR.le
          rw [le_div_iff₀ h2d]; linarith
      _ = (2 : ℝ) ^ up0 / (2 : ℝ) ^ down0 * (sparseLMat22 U0 i j : ℝ) / s := by field_simp
  -- 64 regions
  have hpow := SparseListSquare22.chain_sound s hs stages U0 ((up0 : ℤ) - down0) A hA0 hU0 hAU hsq
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
  -- prefactor
  obtain ⟨hX0, hXp⟩ := prefactorAB_le ea eb heb zb d q z zNum hz ztrace hzc vtrace hvc yn G hpre
  set X := turboOldUniformFactor22 ^ q / z ^ d
  have hY0 : 0 ≤ (yn : ℝ) * (2 : ℝ) ^ G := (Real.rpow_nonneg hX0 _).trans hXp
  -- assemble
  have hfin := final_real Cq yn sumU a G E Kt S hfinal
  rw [← hC]
  have h2G : (2 : ℝ) ^ G ≠ 0 := zpow_ne_zero _ two_ne_zero
  have h2E : (2 : ℝ) ^ E ≠ 0 := zpow_ne_zero _ two_ne_zero
  calc (Cq : ℝ) * (X ^ ((ea : ℝ) / eb) * ∑ k : Fin 11, (A ^ 64) 0 k)
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

#print axioms TailTermFrac22.tail_term_fracAB_le
#print axioms TailTermFrac22.tail_term_fracAB_shift_le
end Spin.Turbo
