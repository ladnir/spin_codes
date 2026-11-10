import SpinCodes.Turbo.TailTermFracAB22
import SpinCodes.Turbo.TailMarker22

/-! One tail occupation term at fractional exponent α = a/b, bounded through
a positive marker (the "dense" components; generalizes `TailMarker22`, which
is α = 1 with the single-step marker squared).

Chain of bounds (geometry: `n` two-step pairs per region, 64 regions):
* `pairedRegionDP22_marker_pow64_le`: `(DP_q)^64 ≤ (1/C(32n,q)/t^q)^64 · M_t^(64n)`
  with the paired marker `M_t = Σ_{a,b} t^(a+b) W(a,b)`,
  `W = pairedProfileMatrix22 B α`;
* pair certificates (`pairedProfileAB_le`): `W(a,b) ≤ C(16,a) C(16,b) Y_ab / sy`;
* base check `tailMarkerPairBaseCheck22`: `M_t ≤ 2^(up-down) U0 / s`;
* `log2(64 n)` scaled squarings (`SparseListSquare22.chain_sound`);
* prefactor `(β^q / z^d)^α ≤ yn · 2^G` (`TailTermFrac22.prefactorAB_le`);
* one final integer comparison `tailMarkerFracFinalCheck22`. -/

namespace Spin.Turbo

/-! ## Core-only definitions -/

/-- `C(16,a) C(16,b) tn^(a+b) td^(32-(a+b)) · Y_ab[i][k]`. -/
def tailMarkerPairTerm22 (Ys : List (List (List (List Int)))) (tn td : Int)
    (i k a b : Nat) : Int :=
  sparseBinom16L22.getD a 0 * sparseBinom16L22.getD b 0 * tn ^ (a + b) *
    td ^ (32 - (a + b)) * ((((Ys.getD a []).getD b []).getD i []).getD k 0)

def tailMarkerPairRow22 (Ys : List (List (List (List Int)))) (tn td : Int)
    (i k a : Nat) : Nat → Int
  | 0 => 0
  | b + 1 => tailMarkerPairRow22 Ys tn td i k a b + tailMarkerPairTerm22 Ys tn td i k a b

def tailMarkerPairSum22 (Ys : List (List (List (List Int)))) (tn td : Int)
    (i k : Nat) : Nat → Int
  | 0 => 0
  | a + 1 => tailMarkerPairSum22 Ys tn td i k a + tailMarkerPairRow22 Ys tn td i k a 17

/-- `U0` shaped and `2^down · s · Σ_{a,b} term ≤ sy · td^32 · 2^up · U0[i][k]`. -/
def tailMarkerPairBaseCheck22 (Ys : List (List (List (List Int)))) (tn td sy s : Int)
    (down up : Nat) (U0 : List (List Int)) : Bool :=
  sparseMatShapeL22 U0 &&
    (List.range 11).all fun i => (List.range 11).all fun k =>
      decide (2 ^ down * s * tailMarkerPairSum22 Ys tn td i k 17 ≤
        sy * td ^ 32 * 2 ^ up * (U0.getD i []).getD k 0)

/-- `yn · sumU · td0^(64q) · 2^T ≤ a · C^63 · tn0^(64q)` (sign-split on `T`). -/
def tailMarkerFracFinalCheck22 (C yn sumU a tn0 td0 : Int) (q : Nat) (T : Int) : Bool :=
  if 0 ≤ T then
    decide (yn * sumU * td0 ^ (64 * q) * 2 ^ T.toNat ≤ a * C ^ 63 * tn0 ^ (64 * q))
  else
    decide (yn * sumU * td0 ^ (64 * q) ≤ a * C ^ 63 * tn0 ^ (64 * q) * 2 ^ (-T).toNat)

end Spin.Turbo

/-! ## Soundness -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22

namespace TailMarkerFracAB22

theorem pairRow_eq (Ys : List (List (List (List Int)))) (tn td : ℤ) (i k a : ℕ) :
    ∀ m, tailMarkerPairRow22 Ys tn td i k a m = ∑ b ∈ range m, tailMarkerPairTerm22 Ys tn td i k a b
  | 0 => rfl
  | m + 1 => by rw [tailMarkerPairRow22, pairRow_eq Ys tn td i k a m, Finset.sum_range_succ]

theorem pairSum_eq (Ys : List (List (List (List Int)))) (tn td : ℤ) (i k : ℕ) :
    ∀ m, tailMarkerPairSum22 Ys tn td i k m =
      ∑ a ∈ range m, ∑ b ∈ range 17, tailMarkerPairTerm22 Ys tn td i k a b
  | 0 => rfl
  | m + 1 => by
    rw [tailMarkerPairSum22, pairSum_eq Ys tn td i k m, Finset.sum_range_succ, pairRow_eq]

theorem pairTerm_bound (ea eb : ℕ) (heb : eb ≠ 0) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (den sy : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den sy R Ys = true)
    (hsy : 0 < sy) (tn td : ℤ) (htn : 0 < tn) (htd : 0 < td) (a b : Fin 17) (i k : Fin 11) :
    ((tn : ℝ) / td) ^ (a.val + b.val) * pairedProfileMatrix22 B ((ea : ℝ) / eb) (a, b) i k ≤
      (tailMarkerPairTerm22 Ys tn td i.val k.val a.val b.val : ℝ) / ((sy : ℝ) * (td : ℝ) ^ 32) := by
  have ha : a.val ≤ 16 := by have := a.isLt; omega
  have hb : b.val ≤ 16 := by have := b.isLt; omega
  have hW := PairedFracChain22.pairedProfileAB_le ea eb heb B hB den sy R hR Ys hrp a b i k
  have hsyR : (0 : ℝ) < sy := by exact_mod_cast hsy
  have htdR : (0 : ℝ) < td := by exact_mod_cast htd
  have htnR : (0 : ℝ) < tn := by exact_mod_cast htn
  unfold tailMarkerPairTerm22
  rw [SparseListPolyD22.binom16L_eq a.val ha, SparseListPolyD22.binom16L_eq b.val hb]
  push_cast
  have h32 : (td : ℝ) ^ 32 = (td : ℝ) ^ (a.val + b.val) * (td : ℝ) ^ (32 - (a.val + b.val)) := by
    rw [← pow_add, Nat.add_sub_cancel' (by omega)]
  have hpa : (0 : ℝ) < (td : ℝ) ^ (a.val + b.val) := pow_pos htdR _
  have hpb : (0 : ℝ) < (td : ℝ) ^ (32 - (a.val + b.val)) := pow_pos htdR _
  have e : (Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b.val : ℝ) * (tn : ℝ) ^ (a.val + b.val) *
      (td : ℝ) ^ (32 - (a.val + b.val)) *
      (((((Ys.getD a.val []).getD b.val []).getD i.val []).getD k.val 0 : ℤ) : ℝ) /
      ((sy : ℝ) * (td : ℝ) ^ 32) =
      ((tn : ℝ) / td) ^ (a.val + b.val) * (((Nat.choose 16 a.val : ℝ) * (Nat.choose 16 b.val : ℝ)) *
        (sparseLMat22 ((Ys.getD a.val []).getD b.val []) i.val k.val : ℝ) / sy) := by
    unfold sparseLMat22
    rw [h32, div_pow]
    field_simp
  rw [e]
  exact mul_le_mul_of_nonneg_left hW (pow_nonneg (div_pos htnR htdR).le _)

/-- Base check soundness: `M_t ≤ 2^(up-down) U0 / s` for `t = tn/td`. -/
theorem base_sound (ea eb : ℕ) (heb : eb ≠ 0) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (den sy : ℤ) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den sy R Ys = true)
    (hsy : 0 < sy) (tn td : ℤ) (htn : 0 < tn) (htd : 0 < td) (s : ℤ) (hs : 0 < s)
    (down up : ℕ) (U0 : List (List Int))
    (hcheck : tailMarkerPairBaseCheck22 Ys tn td sy s down up U0 = true) :
    MatShape U0 ∧ ∀ i k : Fin 11,
      orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B ((ea : ℝ) / eb))
        ((tn : ℝ) / td) i k ≤ (2 : ℝ) ^ ((up : ℤ) - (down : ℤ)) * (sparseLMat22 U0 i k : ℝ) / s := by
  simp only [tailMarkerPairBaseCheck22, Bool.and_eq_true, List.all_eq_true, List.mem_range,
    decide_eq_true_eq] at hcheck
  obtain ⟨hsh, hent⟩ := hcheck
  refine ⟨matShape_of_check U0 hsh, fun i k => ?_⟩
  have hc := hent i.val i.isLt k.val k.isLt
  rw [pairSum_eq] at hc
  have hsyR : (0 : ℝ) < sy := by exact_mod_cast hsy
  have htdR : (0 : ℝ) < td := by exact_mod_cast htd
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hD : (0 : ℝ) < (sy : ℝ) * (td : ℝ) ^ 32 := mul_pos hsyR (pow_pos htdR 32)
  set L : ℤ := ∑ a ∈ range 17, ∑ b ∈ range 17, tailMarkerPairTerm22 Ys tn td i.val k.val a b
    with hL
  have hcR : (2 : ℝ) ^ down * (s : ℝ) * (L : ℝ) ≤
      (sy : ℝ) * (td : ℝ) ^ 32 * (2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ) := by
    unfold sparseLMat22; exact_mod_cast hc
  have hsum : orderedMatrixMarker22 pairedProfileDegree22
      (pairedProfileMatrix22 B ((ea : ℝ) / eb)) ((tn : ℝ) / td) i k ≤
      (L : ℝ) / ((sy : ℝ) * (td : ℝ) ^ 32) := by
    unfold orderedMatrixMarker22
    rw [Matrix.sum_apply, hL, Int.cast_sum, Finset.sum_div, Fintype.sum_prod_type]
    rw [← Fin.sum_univ_eq_sum_range (fun a => ((∑ b ∈ range 17,
      tailMarkerPairTerm22 Ys tn td i.val k.val a b : ℤ) : ℝ) / ((sy : ℝ) * (td : ℝ) ^ 32)) 17]
    refine Finset.sum_le_sum fun a _ => ?_
    rw [Int.cast_sum, Finset.sum_div,
      ← Fin.sum_univ_eq_sum_range (fun b => ((tailMarkerPairTerm22 Ys tn td i.val k.val a.val b : ℤ) : ℝ) /
        ((sy : ℝ) * (td : ℝ) ^ 32)) 17]
    refine Finset.sum_le_sum fun b _ => ?_
    simp only [Matrix.smul_apply, smul_eq_mul, pairedProfileDegree22]
    exact pairTerm_bound ea eb heb B hB den sy R hR Ys hrp hsy tn td htn htd a b i k
  refine hsum.trans ?_
  have e : (2 : ℝ) ^ ((up : ℤ) - (down : ℤ)) * (sparseLMat22 U0 i k : ℝ) / s =
      ((2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ)) / ((2 : ℝ) ^ down * s) := by
    rw [zpow_sub₀ (by norm_num), zpow_natCast, zpow_natCast]
    field_simp
  rw [e, div_le_div_iff₀ hD (mul_pos (pow_pos (by norm_num) _) hsR)]
  calc (L : ℝ) * ((2 : ℝ) ^ down * s) = (2 : ℝ) ^ down * (s : ℝ) * (L : ℝ) := by ring
    _ ≤ (sy : ℝ) * (td : ℝ) ^ 32 * (2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ) := hcR
    _ = _ := by ring

theorem final_real (C yn sumU a tn0 td0 : ℤ) (q : ℕ) (T : ℤ)
    (h : tailMarkerFracFinalCheck22 C yn sumU a tn0 td0 q T = true) :
    (yn : ℝ) * sumU * (td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ T ≤
      (a : ℝ) * (C : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) := by
  unfold tailMarkerFracFinalCheck22 at h
  split_ifs at h with hT0
  · have h' := of_decide_eq_true h
    have hR : ((yn * sumU * td0 ^ (64 * q) * 2 ^ T.toNat : ℤ) : ℝ) ≤
        ((a * C ^ 63 * tn0 ^ (64 * q) : ℤ) : ℝ) := by exact_mod_cast h'
    push_cast at hR
    have hz : (2 : ℝ) ^ T = (2 : ℝ) ^ T.toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg hT0]
    rw [hz]; exact hR
  · have h' := of_decide_eq_true h
    have hR : ((yn * sumU * td0 ^ (64 * q) : ℤ) : ℝ) ≤
        ((a * C ^ 63 * tn0 ^ (64 * q) * 2 ^ (-T).toNat : ℤ) : ℝ) := by exact_mod_cast h'
    push_cast at hR
    have hpos : (0 : ℝ) < (2 : ℝ) ^ (-T).toNat := by positivity
    have hz : (2 : ℝ) ^ T = 1 / (2 : ℝ) ^ (-T).toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg (by omega), _root_.zpow_neg, one_div, inv_inv]
    rw [hz, mul_one_div, div_le_iff₀ hpos]
    exact hR

/-- One tail term at α = a/b bounded through the marker
`t = tn0·2^kp / (td0·2^kn)`. Every certificate condition is a Boolean check. -/
theorem tail_marker_fracAB_le (ea eb : ℕ) (heb : eb ≠ 0) (n d zb : ℕ)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (q : ℕ) (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose (32 * n) q : ℝ)) (hCpos : 0 < Cq)
    (den sy : ℤ) (hsy : 0 < sy) (R : List (List (List Int)))
    (hR : ∀ (a : Fin 17) (i k : Fin 11),
      B a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / den)
    (Ys : List (List (List (List Int)))) (hrp : pairedRpowABCheckL22 ea eb den sy R Ys = true)
    (tn0 td0 : ℤ) (htn : 0 < tn0) (htd : 0 < td0) (kp kn : ℕ)
    (s : ℤ) (S : ℕ) (hS : s = 2 ^ S) (down up : ℕ) (U0 : List (List Int))
    (hbase : tailMarkerPairBaseCheck22 Ys (tn0 * 2 ^ kp) (td0 * 2 ^ kn) sy s down up U0 = true)
    (stages : List (Nat × Nat × List (List Int))) (hlen : 2 ^ stages.length = n * 64)
    (hchain : sparseSquareChainCheckL22 s U0 stages = true)
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ zb)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum zb d ztrace = true)
    (vtrace : List (Bool × Nat × Int))
    (hvc : sparseScaledPowerCheck22 65535 0 (8 * q) vtrace = true) (yn G : ℤ)
    (hpre : powABPow2Check22 ea eb
      (256 * (q : ℤ) + sparseScaledPowerExp22 zb ztrace + sparseScaledPowerExp22 0 vtrace) G
      (sparseScaledPowerMant22 ztrace * sparseScaledPowerMant22 vtrace) yn = true)
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailMarkerFracFinalCheck22 Cq yn
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j))
      a tn0 td0 q
      (G + SparseListSquare22.chainExp ((up : ℤ) - (down : ℤ)) stages + (Kt : ℤ) +
        ((64 * q * kn : ℕ) : ℤ) - (S : ℤ) - ((64 * q * kp : ℕ) : ℤ)) = true) :
    (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
  have hs : (0 : ℤ) < s := by rw [hS]; positivity
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hCR : (0 : ℝ) < Cq := by exact_mod_cast hCpos
  have htnI : (0 : ℤ) < tn0 * 2 ^ kp := mul_pos htn (by positivity)
  have htdI : (0 : ℤ) < td0 * 2 ^ kn := mul_pos htd (by positivity)
  have htn0R : (0 : ℝ) < tn0 := by exact_mod_cast htn
  have htd0R : (0 : ℝ) < td0 := by exact_mod_cast htd
  obtain ⟨t, htdef⟩ : ∃ t : ℝ, t = ((tn0 * 2 ^ kp : ℤ) : ℝ) / ((td0 * 2 ^ kn : ℤ) : ℝ) :=
    ⟨_, rfl⟩
  have ht : 0 < t := by
    rw [htdef]; exact div_pos (by exact_mod_cast htnI) (by exact_mod_cast htdI)
  set W := pairedProfileMatrix22 B ((ea : ℝ) / eb) with hWdef
  set M := orderedMatrixMarker22 pairedProfileDegree22 W t with hMdef
  -- paired marker base
  obtain ⟨hU0, hMU⟩ := base_sound ea eb heb B hB den sy R hR Ys hrp hsy _ _ htnI htdI s hs
    down up U0 hbase
  rw [← htdef] at hMU
  have hW0 := PairedFracChain22.pairedProfile_nonneg B hB ((ea : ℝ) / eb)
  have hM0 : ∀ i j, 0 ≤ M i j := by
    intro i j
    rw [hMdef]
    unfold orderedMatrixMarker22
    rw [Matrix.sum_apply]
    exact Finset.sum_nonneg fun ab _ => by
      rw [Matrix.smul_apply, smul_eq_mul]; exact mul_nonneg (pow_nonneg ht.le _) (hW0 ab _ _)
  -- squarings
  have hpow := SparseListSquare22.chain_sound s hs stages U0 ((up : ℤ) - (down : ℤ))
    M hM0 hU0 hMU hchain
  rw [hlen] at hpow
  set E := SparseListSquare22.chainExp ((up : ℤ) - (down : ℤ)) stages
  set Ul := sparseSquareChainLastL22 U0 stages
  set sumU := sparseSum11Int22 (fun j => sparseLMat22 Ul 0 j)
  have hmass : (∑ k : Fin 11, (M ^ (n * 64)) 0 k) ≤ (2 : ℝ) ^ E * (sumU : ℝ) / s := by
    have hsum : (sumU : ℝ) = ∑ k : Fin 11, (sparseLMat22 Ul 0 k : ℝ) := sum11_cast _
    rw [hsum, Finset.mul_sum, Finset.sum_div]
    exact Finset.sum_le_sum fun k _ => by simpa [mul_div_assoc] using hpow 0 k
  -- paired coefficient through the marker
  set c : ℝ := 1 / (Cq : ℝ) / t ^ q with hc
  have hc0 : 0 ≤ c := div_nonneg (one_div_nonneg.mpr hCR.le) (pow_nonneg ht.le _)
  have hDP : ∀ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k ≤
      c ^ 64 * (M ^ (n * 64)) 0 k := by
    intro k
    have h := pairedRegionDP22_marker_pow64_le B ((ea : ℝ) / eb) hB n q t ht 0 k
    rw [← hC] at h
    exact h
  have hDPsum : (∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤
      c ^ 64 * ((2 : ℝ) ^ E * (sumU : ℝ) / s) := by
    calc (∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k)
        ≤ ∑ k : Fin 11, c ^ 64 * (M ^ (n * 64)) 0 k := Finset.sum_le_sum fun k _ => hDP k
      _ = c ^ 64 * ∑ k : Fin 11, (M ^ (n * 64)) 0 k := (Finset.mul_sum _ _ _).symm
      _ ≤ _ := mul_le_mul_of_nonneg_left hmass (pow_nonneg hc0 _)
  have hD0 : ∀ i j, 0 ≤ pairedRegionDP22 B ((ea : ℝ) / eb) n q i j := by
    intro i j
    unfold pairedRegionDP22
    rw [Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (PairedFracChain22.dp_nonneg _ hW0 n q i j)
  have hDP0 : 0 ≤ ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k :=
    Finset.sum_nonneg fun k _ => markerPowNonneg22 _ hD0 64 0 k
  -- prefactor
  obtain ⟨hX0, hXp⟩ := TailTermFrac22.prefactorAB_le ea eb heb zb d q z zNum hz ztrace hzc
    vtrace hvc yn G hpre
  set X := turboOldUniformFactor22 ^ q / z ^ d
  have hY0 : 0 ≤ (yn : ℝ) * (2 : ℝ) ^ G := (Real.rpow_nonneg hX0 _).trans hXp
  -- assemble
  set T : ℤ := G + E + (Kt : ℤ) + ((64 * q * kn : ℕ) : ℤ) - (S : ℤ) - ((64 * q * kp : ℕ) : ℤ)
    with hTdef
  have hfin := final_real Cq yn sumU a tn0 td0 q T hfinal
  rw [← hC]
  have hK : (0 : ℝ) < (2 : ℝ) ^ Kt := pow_pos (by norm_num) _
  have hden : (0 : ℝ) < (Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt :=
    mul_pos (mul_pos (pow_pos hCR _) (pow_pos htn0R _)) hK
  calc (Cq : ℝ) * (X ^ ((ea : ℝ) / eb) *
        ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k)
      ≤ (Cq : ℝ) * (((yn : ℝ) * (2 : ℝ) ^ G) * (c ^ 64 * ((2 : ℝ) ^ E * (sumU : ℝ) / s))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul hXp hDPsum hDP0 hY0) hCR.le
    _ = ((yn : ℝ) * sumU * (td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ T) /
          ((Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt) := by
        have ht64 : (t ^ q) ^ 64 = ((tn0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ (64 * q * kp)) /
            ((td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ (64 * q * kn)) := by
          rw [htdef, ← pow_mul, show q * 64 = 64 * q by ring]
          push_cast
          rw [div_pow, mul_pow, mul_pow, ← pow_mul, ← pow_mul,
            show kp * (64 * q) = 64 * q * kp by ring, show kn * (64 * q) = 64 * q * kn by ring]
        have e2 : (2 : ℝ) ^ T = (2 : ℝ) ^ G * (2 : ℝ) ^ E * (2 : ℝ) ^ Kt *
            (2 : ℝ) ^ (64 * q * kn) / ((2 : ℝ) ^ S * (2 : ℝ) ^ (64 * q * kp)) := by
          rw [hTdef, zpow_sub₀ (by norm_num), zpow_sub₀ (by norm_num), zpow_add₀ (by norm_num),
            zpow_add₀ (by norm_num), zpow_add₀ (by norm_num), zpow_natCast, zpow_natCast,
            zpow_natCast, zpow_natCast]
          ring
        rw [hc, div_pow, ht64, e2, hS]
        push_cast
        have hCne : (Cq : ℝ) ≠ 0 := hCR.ne'
        have h2G : (2 : ℝ) ^ G ≠ 0 := zpow_ne_zero _ two_ne_zero
        have h2E : (2 : ℝ) ^ E ≠ 0 := zpow_ne_zero _ two_ne_zero
        field_simp
    _ ≤ ((a : ℝ) * (Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q)) /
          ((Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt) :=
        div_le_div_of_nonneg_right hfin hden.le
    _ = (a : ℝ) / 2 ^ Kt := by
        have h1 : (Cq : ℝ) ^ 63 ≠ 0 := pow_ne_zero _ hCR.ne'
        have h2 : (tn0 : ℝ) ^ (64 * q) ≠ 0 := pow_ne_zero _ htn0R.ne'
        field_simp

end TailMarkerFracAB22

#print axioms TailMarkerFracAB22.tail_marker_fracAB_le
end Spin.Turbo
