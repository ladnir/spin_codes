import SpinCodes.Turbo.SparseListPolyD22
import SpinCodes.Turbo.SparseListSquare22
import SpinCodes.Turbo.SparseScaledPower22
import SpinCodes.Turbo.PairedSquare22
import SpinCodes.Turbo.RS8UniformReference22
import SpinCodes.Turbo.MatrixMarker22

/-! One K18 tail occupation term (α = 1, 64 pairs, 64 regions, d = 52428) for
large `q`, bounded through a positive marker instead of an exact coefficient.

Chain of bounds:
* `pairedRegionDP22_marker_pow64_le`: the 64-region product of the paired
  coefficient is at most `(1/C(2048,q)/t^q)^64 · (marker_t)^(64·64)`;
* at α = 1 the paired marker is `A_t * A_t` with
  `A_t = Σ_a t^a C(16,a) • B a` (`marker_one_eq`), so the marker power is
  `A_t ^ 8192 = A_t ^ (2^13)`;
* `A_t ≤ 2^(up-down) U0 / s` from integer upper bounds on `B` (`base_sound`,
  check `tailMarkerBaseCheck22`);
* 13 scaled squarings (`SparseListSquare22.chain_sound`);
* `z^d ≥ m 2^(-Ez)` (`SparseScaledPower22.scaledPower_sound`);
* one final integer comparison (`tailMarkerFinalCheck22`) with every power of
  two folded into one exponent `T` (`tailMarkerExp22`).

Marker parameter. `t = (tn0 · 2^kp) / (td0 · 2^kn)` with positive integers
`tn0, td0` and natural shifts `kp, kn`. The base check works with
`tn = tn0·2^kp`, `td = td0·2^kn` raised to powers `≤ 16` only. The final check
raises `tn0` and `td0` (NOT the shifts, which are folded into `T`) to the power
`64·q ≤ 131072`, and `C(2048,q)` to the power 63. Sizes: `tn0^(64q)` has about
`64·q·log2 tn0` bits, so keep `tn0, td0 ≤ 2^12` (≤ ~1.6·10^6 bits at q = 2048);
the large dynamic range of `t` (q near 2048 needs `t ≫ 2^13`) goes into
`kp`/`kn`, which cost nothing. All integer powers are `Int.pow`, which reduces
to the GMP-accelerated `Nat.pow` in the kernel. -/

namespace Spin.Turbo

/-! ## Core-only definitions -/

/-- One weighted base term `C(16,a) · tn^a · td^(16-a) · row a i k`. -/
def tailMarkerTerm22 (row : Nat → Nat → Nat → Int) (tn td : Int) (i k a : Nat) : Int :=
  sparseBinom16L22.getD a 0 * tn ^ a * td ^ (16 - a) * row a i k

/-- `Σ_{a < n}` of `tailMarkerTerm22`. -/
def tailMarkerSumAux22 (row : Nat → Nat → Nat → Int) (tn td : Int) (i k : Nat) : Nat → Int
  | 0 => 0
  | n + 1 => tailMarkerSumAux22 row tn td i k n + tailMarkerTerm22 row tn td i k n

/-- `U0` is shaped and `2^down · s · Σ_a C(16,a) tn^a td^(16-a) row a i k ≤
den · td^16 · 2^up · U0[i][k]`, i.e. `A_t ≤ 2^(up-down) U0 / s` for
`t = tn/td`, `B a ≤ row a / den`. -/
def tailMarkerBaseCheck22 (row : Nat → Nat → Nat → Int) (tn td den s : Int) (down up : Nat)
    (U0 : List (List Int)) : Bool :=
  sparseMatShapeL22 U0 &&
    (List.range 11).all fun i => (List.range 11).all fun k =>
      decide (2 ^ down * s * tailMarkerSumAux22 row tn td i k 17 ≤
        den * td ^ 16 * 2 ^ up * (U0.getD i []).getD k 0)

/-- Total binary exponent: `256q` (β^q), `Kt` (target), `64q·kn` and `-64q·kp`
(marker shifts), `Ez` (z^d), `E` (squaring chain), `-S` (chain scale). -/
def tailMarkerExp22 (q : Nat) (Ez E : Int) (Kt S kp kn : Nat) : Int :=
  ((256 * q + Kt + 64 * q * kn : Nat) : Int) + Ez + E - ((S + 64 * q * kp : Nat) : Int)

/-- `sumU · td0^(64q) · 2^T ≤ a · 65535^(8q) · m · C^63 · tn0^(64q)`. -/
def tailMarkerFinalCheck22 (C sumU a tn0 td0 : Int) (q : Nat) (T m : Int) : Bool :=
  if 0 ≤ T then
    decide (sumU * td0 ^ (64 * q) * 2 ^ T.toNat ≤
      a * 65535 ^ (8 * q) * m * C ^ 63 * tn0 ^ (64 * q))
  else
    decide (sumU * td0 ^ (64 * q) ≤
      a * 65535 ^ (8 * q) * m * C ^ 63 * tn0 ^ (64 * q) * 2 ^ (-T).toNat)

end Spin.Turbo

/-! ## Soundness -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22

namespace TailMarker22

/-- The single-step marker `A_t = Σ_a t^a C(16,a) • B a`. -/
def singleMarker22 (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (t : ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  ∑ a : Fin 17, t ^ a.val • ((Nat.choose 16 a.val : ℝ) • B a)

theorem singleMarker22_nonneg (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (t : ℝ) (ht : 0 ≤ t) (i j : Fin 11) :
    0 ≤ singleMarker22 B t i j := by
  unfold singleMarker22
  rw [Matrix.sum_apply]
  refine Finset.sum_nonneg fun a _ => ?_
  simp only [Matrix.smul_apply, smul_eq_mul]
  exact mul_nonneg (pow_nonneg ht _) (mul_nonneg (Nat.cast_nonneg _) (hB a i j))

/-- At α = 1 the paired marker is the square of the single-step marker. -/
theorem marker_one_eq (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (t : ℝ) :
    orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B 1) t =
      singleMarker22 B t * singleMarker22 B t := by
  unfold orderedMatrixMarker22 singleMarker22
  rw [Finset.sum_mul_sum, ← Finset.univ_product_univ, Finset.sum_product]
  apply Finset.sum_congr rfl
  intro a _
  apply Finset.sum_congr rfl
  intro b _
  rw [PairedSquare22.pairedProfileMatrix22_one, smul_mul_smul_comm]
  simp only [pairedProfileDegree22, pow_add]
  rw [smul_mul_smul_comm, smul_mul_smul_comm]

theorem marker_pow_eq (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (t : ℝ) :
    orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B 1) t ^ (64 * 64) =
      singleMarker22 B t ^ (2 ^ 13) := by
  rw [marker_one_eq, ← sq, ← pow_mul]
  norm_num

theorem sumAux_eq (row : ℕ → ℕ → ℕ → ℤ) (tn td : ℤ) (i k : ℕ) :
    ∀ n, tailMarkerSumAux22 row tn td i k n =
      ∑ a ∈ range n, tailMarkerTerm22 row tn td i k a
  | 0 => rfl
  | n + 1 => by rw [tailMarkerSumAux22, sumAux_eq row tn td i k n, Finset.sum_range_succ]

theorem term_bound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (row : ℕ → ℕ → ℕ → ℤ) (den : ℤ) (hden : 0 < den)
    (hrow : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / den)
    (tn td : ℤ) (htn : 0 < tn) (htd : 0 < td) (a : Fin 17) (i k : Fin 11) :
    ((tn : ℝ) / td) ^ a.val * ((Nat.choose 16 a.val : ℝ) * B a i k) ≤
      (tailMarkerTerm22 row tn td i.val k.val a.val : ℝ) / ((den : ℝ) * (td : ℝ) ^ 16) := by
  have ha : a.val ≤ 16 := by have := a.isLt; omega
  have hb : B a i k ≤ (row a.val i.val k.val : ℝ) / den := hrow a.val ha i k
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have htdR : (0 : ℝ) < td := by exact_mod_cast htd
  have htnR : (0 : ℝ) < tn := by exact_mod_cast htn
  unfold tailMarkerTerm22
  rw [SparseListPolyD22.binom16L_eq a.val ha]
  push_cast
  have h16 : (td : ℝ) ^ 16 = (td : ℝ) ^ a.val * (td : ℝ) ^ (16 - a.val) := by
    rw [← pow_add, Nat.add_sub_cancel' ha]
  have hpa : (0 : ℝ) < (td : ℝ) ^ a.val := pow_pos htdR _
  have hpb : (0 : ℝ) < (td : ℝ) ^ (16 - a.val) := pow_pos htdR _
  have e : (Nat.choose 16 a.val : ℝ) * (tn : ℝ) ^ a.val * (td : ℝ) ^ (16 - a.val) *
      (row a.val i.val k.val : ℝ) / ((den : ℝ) * (td : ℝ) ^ 16) =
      ((tn : ℝ) / td) ^ a.val * ((Nat.choose 16 a.val : ℝ) * ((row a.val i.val k.val : ℝ) / den)) := by
    rw [h16, div_pow]
    field_simp
  rw [e]
  apply mul_le_mul_of_nonneg_left _ (pow_nonneg (div_pos htnR htdR).le _)
  exact mul_le_mul_of_nonneg_left hb (Nat.cast_nonneg _)

/-- Soundness of the base check: `A_t ≤ 2^(up-down) U0 / s` for `t = tn/td`. -/
theorem base_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (row : ℕ → ℕ → ℕ → ℤ) (den : ℤ) (hden : 0 < den)
    (hrow : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / den)
    (tn td : ℤ) (htn : 0 < tn) (htd : 0 < td) (s : ℤ) (hs : 0 < s) (down up : ℕ)
    (U0 : List (List Int)) (hcheck : tailMarkerBaseCheck22 row tn td den s down up U0 = true) :
    MatShape U0 ∧ ∀ i k : Fin 11, singleMarker22 B ((tn : ℝ) / td) i k ≤
      (2 : ℝ) ^ ((up : ℤ) - (down : ℤ)) * (sparseLMat22 U0 i k : ℝ) / s := by
  simp only [tailMarkerBaseCheck22, Bool.and_eq_true, List.all_eq_true, List.mem_range,
    decide_eq_true_eq] at hcheck
  obtain ⟨hsh, hent⟩ := hcheck
  refine ⟨matShape_of_check U0 hsh, fun i k => ?_⟩
  have hc := hent i.val i.isLt k.val k.isLt
  rw [sumAux_eq, Finset.sum_range] at hc
  have hdR : (0 : ℝ) < den := by exact_mod_cast hden
  have htdR : (0 : ℝ) < td := by exact_mod_cast htd
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hD : (0 : ℝ) < (den : ℝ) * (td : ℝ) ^ 16 := mul_pos hdR (pow_pos htdR 16)
  set L : ℤ := ∑ a : Fin 17, tailMarkerTerm22 row tn td i.val k.val a.val with hL
  have hcR : (2 : ℝ) ^ down * (s : ℝ) * (L : ℝ) ≤
      (den : ℝ) * (td : ℝ) ^ 16 * (2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ) := by
    unfold sparseLMat22; exact_mod_cast hc
  have hsum : singleMarker22 B ((tn : ℝ) / td) i k ≤ (L : ℝ) / ((den : ℝ) * (td : ℝ) ^ 16) := by
    unfold singleMarker22
    rw [Matrix.sum_apply, hL, Int.cast_sum, Finset.sum_div]
    refine Finset.sum_le_sum fun a _ => ?_
    simp only [Matrix.smul_apply, smul_eq_mul]
    exact term_bound B row den hden hrow tn td htn htd a i k
  refine hsum.trans ?_
  have e : (2 : ℝ) ^ ((up : ℤ) - (down : ℤ)) * (sparseLMat22 U0 i k : ℝ) / s =
      ((2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ)) / ((2 : ℝ) ^ down * s) := by
    rw [zpow_sub₀ (by norm_num), zpow_natCast, zpow_natCast]
    field_simp
  rw [e, div_le_div_iff₀ hD (mul_pos (pow_pos (by norm_num) _) hsR)]
  calc (L : ℝ) * ((2 : ℝ) ^ down * s) = (2 : ℝ) ^ down * (s : ℝ) * (L : ℝ) := by ring
    _ ≤ (den : ℝ) * (td : ℝ) ^ 16 * (2 : ℝ) ^ up * (sparseLMat22 U0 i k : ℝ) := hcR
    _ = _ := by ring

theorem final_real (C sumU a tn0 td0 : ℤ) (q : ℕ) (T m : ℤ)
    (h : tailMarkerFinalCheck22 C sumU a tn0 td0 q T m = true) :
    (sumU : ℝ) * (td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ T ≤
      (a : ℝ) * (65535 : ℝ) ^ (8 * q) * m * (C : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) := by
  unfold tailMarkerFinalCheck22 at h
  split_ifs at h with hT0
  · have h' := of_decide_eq_true h
    have hR : ((sumU * td0 ^ (64 * q) * 2 ^ T.toNat : ℤ) : ℝ) ≤
        ((a * 65535 ^ (8 * q) * m * C ^ 63 * tn0 ^ (64 * q) : ℤ) : ℝ) := by
      exact_mod_cast h'
    push_cast at hR
    have hz : (2 : ℝ) ^ T = (2 : ℝ) ^ T.toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg hT0]
    rw [hz]; exact hR
  · have h' := of_decide_eq_true h
    have hR : ((sumU * td0 ^ (64 * q) : ℤ) : ℝ) ≤
        ((a * 65535 ^ (8 * q) * m * C ^ 63 * tn0 ^ (64 * q) * 2 ^ (-T).toNat : ℤ) : ℝ) := by
      exact_mod_cast h'
    push_cast at hR
    have hpos : (0 : ℝ) < (2 : ℝ) ^ (-T).toNat := by positivity
    have hz : (2 : ℝ) ^ T = 1 / (2 : ℝ) ^ (-T).toNat := by
      rw [← zpow_natCast, Int.toNat_of_nonneg (by omega), _root_.zpow_neg, one_div, inv_inv]
    rw [hz, mul_one_div, div_le_iff₀ hpos]
    exact hR

theorem exp_real (q : ℕ) (Ez E : ℤ) (Kt S kp kn : ℕ) :
    (2 : ℝ) ^ tailMarkerExp22 q Ez E Kt S kp kn =
      (2 : ℝ) ^ (256 * q) * (2 : ℝ) ^ Kt * (2 : ℝ) ^ (64 * q * kn) * (2 : ℝ) ^ Ez *
        (2 : ℝ) ^ E / ((2 : ℝ) ^ S * (2 : ℝ) ^ (64 * q * kp)) := by
  unfold tailMarkerExp22
  rw [zpow_sub₀ (by norm_num), zpow_add₀ (by norm_num), zpow_add₀ (by norm_num),
    zpow_natCast, zpow_natCast, pow_add, pow_add, pow_add]

/-- One K18 tail term bounded through the marker `t = tn0·2^kp / (td0·2^kn)`.
All certificate conditions are Boolean checks for `decide +kernel`. -/
theorem tail_marker_le (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (q : ℕ) (Cq : ℤ) (hC : (Cq : ℝ) = (Nat.choose 2048 q : ℝ)) (hCpos : 0 < Cq)
    (row : ℕ → ℕ → ℕ → ℤ) (den : ℤ) (hden : 0 < den)
    (hrow : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11,
      B ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / den)
    (tn0 td0 : ℤ) (htn : 0 < tn0) (htd : 0 < td0) (kp kn : ℕ)
    (s : ℤ) (S : ℕ) (hS : s = 2 ^ S) (down up : ℕ) (U0 : List (List Int))
    (hbase : tailMarkerBaseCheck22 row (tn0 * 2 ^ kp) (td0 * 2 ^ kn) den s down up U0 = true)
    (stages : List (Nat × Nat × List (List Int))) (hlen : stages.length = 13)
    (hchain : sparseSquareChainCheckL22 s U0 stages = true)
    (z : ℝ) (zNum : ℤ) (hz : z = (zNum : ℝ) / 2 ^ 80)
    (ztrace : List (Bool × Nat × Int)) (hzc : sparseScaledPowerCheck22 zNum 80 52428 ztrace = true)
    (a : ℤ) (Kt : ℕ)
    (hfinal : tailMarkerFinalCheck22 Cq
      (sparseSum11Int22 (fun j => sparseLMat22 (sparseSquareChainLastL22 U0 stages) 0 j))
      a tn0 td0 q
      (tailMarkerExp22 q (sparseScaledPowerExp22 80 ztrace)
        (SparseListSquare22.chainExp ((up : ℤ) - (down : ℤ)) stages) Kt S kp kn)
      (sparseScaledPowerMant22 ztrace) = true) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k) ≤ (a : ℝ) / 2 ^ Kt := by
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
  -- single-step marker
  obtain ⟨hU0, hAU⟩ := base_sound B row den hden hrow _ _ htnI htdI s hs down up U0 hbase
  rw [← htdef] at hAU
  have hA0 := singleMarker22_nonneg B hB t ht.le
  -- 13 squarings
  have hpow := SparseListSquare22.chain_sound s hs stages U0 ((up : ℤ) - (down : ℤ))
    (singleMarker22 B t) hA0 hU0 hAU hchain
  rw [hlen] at hpow
  set E := SparseListSquare22.chainExp ((up : ℤ) - (down : ℤ)) stages
  set Ul := sparseSquareChainLastL22 U0 stages
  set sumU := sparseSum11Int22 (fun j => sparseLMat22 Ul 0 j)
  have hmass : (∑ k : Fin 11, (singleMarker22 B t ^ (2 ^ 13)) 0 k) ≤
      (2 : ℝ) ^ E * (sumU : ℝ) / s := by
    have hsum : (sumU : ℝ) = ∑ k : Fin 11, (sparseLMat22 Ul 0 k : ℝ) := sum11_cast _
    rw [hsum, Finset.mul_sum, Finset.sum_div]
    exact Finset.sum_le_sum fun k _ => by simpa [mul_div_assoc] using hpow 0 k
  -- paired coefficient through the marker
  set c : ℝ := 1 / (Cq : ℝ) / t ^ q with hc
  have hc0 : 0 ≤ c := div_nonneg (one_div_nonneg.mpr hCR.le) (pow_nonneg ht.le _)
  have hDP : ∀ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k ≤
      c ^ 64 * (singleMarker22 B t ^ (2 ^ 13)) 0 k := by
    intro k
    have h := pairedRegionDP22_marker_pow64_le B 1 hB 64 q t ht 0 k
    rw [show (32 * 64 : ℕ) = 2048 from rfl, marker_pow_eq, ← hC] at h
    exact h
  have hDPsum : (∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k) ≤
      c ^ 64 * ((2 : ℝ) ^ E * (sumU : ℝ) / s) := by
    calc (∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k)
        ≤ ∑ k : Fin 11, c ^ 64 * (singleMarker22 B t ^ (2 ^ 13)) 0 k :=
          Finset.sum_le_sum fun k _ => hDP k
      _ = c ^ 64 * ∑ k : Fin 11, (singleMarker22 B t ^ (2 ^ 13)) 0 k := (Finset.mul_sum _ _ _).symm
      _ ≤ _ := mul_le_mul_of_nonneg_left hmass (pow_nonneg hc0 _)
  have hD0 : ∀ i j, 0 ≤ pairedRegionDP22 B 1 64 q i j := by
    intro i j
    rw [PairedSquare22.pairedRegionDP22_one_eq, Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (by positivity) (SparseRegionalPower22.count_nonneg B hB _ q i j)
  have hDP0 : 0 ≤ ∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k :=
    Finset.sum_nonneg fun k _ => markerPowNonneg22 _ hD0 64 0 k
  -- z^d
  obtain ⟨hm, hzlow⟩ := SparseScaledPower22.scaledPower_sound zNum 80 52428 ztrace hzc
  set m := sparseScaledPowerMant22 ztrace
  set Ez := sparseScaledPowerExp22 80 ztrace
  have hmR : (0 : ℝ) < m := by exact_mod_cast hm
  have hzd : (m : ℝ) * (2 : ℝ) ^ (-Ez) ≤ z ^ 52428 := by rw [hz]; exact hzlow
  -- assemble
  set T := tailMarkerExp22 q Ez E Kt S kp kn with hTdef
  have hfin := final_real Cq sumU a tn0 td0 q T m hfinal
  have hbeta : turboOldUniformFactor22 ^ q = (2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) := by
    unfold turboOldUniformFactor22; rw [div_pow, ← pow_mul, ← pow_mul]
  rw [Real.rpow_one, ← hC]
  have h65 : (0 : ℝ) < (65535 : ℝ) ^ (8 * q) := pow_pos (by norm_num) _
  have h2Ez : (0 : ℝ) < (2 : ℝ) ^ Ez := zpow_pos (by norm_num) _
  have h2E : (0 : ℝ) < (2 : ℝ) ^ E := zpow_pos (by norm_num) _
  have hK : (0 : ℝ) < (2 : ℝ) ^ Kt := pow_pos (by norm_num) _
  have hY : (0 : ℝ) < (65535 : ℝ) ^ (8 * q) * m * (Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt :=
    mul_pos (mul_pos (mul_pos (mul_pos h65 hmR) (pow_pos hCR _)) (pow_pos htn0R _)) hK
  have hmz : (0 : ℝ) < (m : ℝ) * (2 : ℝ) ^ (-Ez) := mul_pos hmR (zpow_pos (by norm_num) _)
  calc (Cq : ℝ) * (turboOldUniformFactor22 ^ q / z ^ 52428 *
        ∑ k : Fin 11, ((pairedRegionDP22 B 1 64 q) ^ 64) 0 k)
      ≤ (Cq : ℝ) * ((2 : ℝ) ^ (256 * q) / (65535 : ℝ) ^ (8 * q) / ((m : ℝ) * (2 : ℝ) ^ (-Ez)) *
          (c ^ 64 * ((2 : ℝ) ^ E * (sumU : ℝ) / s))) := by
        apply mul_le_mul_of_nonneg_left _ hCR.le
        apply mul_le_mul _ hDPsum hDP0 (div_nonneg (div_nonneg (by positivity) h65.le) hmz.le)
        rw [hbeta]
        exact div_le_div_of_nonneg_left (div_nonneg (by positivity) h65.le) hmz hzd
    _ = ((sumU : ℝ) * (td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ T) /
          ((65535 : ℝ) ^ (8 * q) * m * (Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt) := by
        have ht64 : (t ^ q) ^ 64 = ((tn0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ (64 * q * kp)) /
            ((td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ (64 * q * kn)) := by
          rw [htdef, ← pow_mul, show q * 64 = 64 * q by ring]
          push_cast
          rw [div_pow, mul_pow, mul_pow, ← pow_mul, ← pow_mul,
            show kp * (64 * q) = 64 * q * kp by ring, show kn * (64 * q) = 64 * q * kn by ring]
        rw [hc, div_pow, ht64, hTdef, exp_real, hS, _root_.zpow_neg]
        push_cast
        have hCne : (Cq : ℝ) ≠ 0 := hCR.ne'
        have hmne : (m : ℝ) ≠ 0 := hmR.ne'
        have htnne : (tn0 : ℝ) ≠ 0 := htn0R.ne'
        have htdne : (td0 : ℝ) ≠ 0 := htd0R.ne'
        have hEzne : (2 : ℝ) ^ Ez ≠ 0 := h2Ez.ne'
        have hEne : (2 : ℝ) ^ E ≠ 0 := h2E.ne'
        field_simp
    _ ≤ (a : ℝ) / 2 ^ Kt := by
        rw [div_le_div_iff₀ hY hK]
        calc (sumU : ℝ) * (td0 : ℝ) ^ (64 * q) * (2 : ℝ) ^ T * 2 ^ Kt
            ≤ (a : ℝ) * (65535 : ℝ) ^ (8 * q) * m * (Cq : ℝ) ^ 63 * (tn0 : ℝ) ^ (64 * q) * 2 ^ Kt :=
              mul_le_mul_of_nonneg_right hfin hK.le
          _ = _ := by ring

end TailMarker22

/-! ## Kernel smoke tests on toy data (each new Boolean check reduces). -/

/-- Toy base: `row ≡ 1`, `tn = 3·2^2`, `td = 5·2^1`, `den = 7`, `s = 2^4`. Then
`L = Σ_a C(16,a) 12^a 10^(16-a) = 22^16`, and the check
`2^0 · 16 · 22^16 ≤ 7 · 10^16 · 2^4 · U0` holds iff `U0 ≥ 22^16 / (7·10^16) ≈ 43567`. -/
example : tailMarkerBaseCheck22 (fun _ _ _ => 1) (3 * 2 ^ 2) (5 * 2 ^ 1) 7 (2 ^ 4) 0 4
    (List.replicate 11 (List.replicate 11 45000)) = true := by decide +kernel

/-- The base check rejects an undersized `U0`. -/
example : tailMarkerBaseCheck22 (fun _ _ _ => 1) (3 * 2 ^ 2) (5 * 2 ^ 1) 7 (2 ^ 4) 0 4
    (List.replicate 11 (List.replicate 11 40000)) = false := by decide +kernel

example : tailMarkerExp22 2048 (-5) 7 256 704 11 0 = 256 * 2048 + 256 - 5 + 7 - 704 - 64 * 2048 * 11 := by
  decide +kernel

/-- Final check, negative-exponent branch. -/
example : tailMarkerFinalCheck22 1 1 1 1 1 1 (-5) 1 = true := by decide +kernel

/-- Final check, nonnegative-exponent branch. -/
example : tailMarkerFinalCheck22 1 1 1 1 1 1 5 1 = true := by decide +kernel

/-- Endpoint-size stress test: `q = 2048`, `C ≈ 2^2045` (so `C^63` ≈ 1.3·10^5 bits),
`tn0 = 4095`, `td0 = 4093` (each raised to `64q = 131072`, ≈ 1.57·10^6 bits). -/
example : tailMarkerFinalCheck22 (3 ^ 1290) 1 1 4095 4093 2048 (-10) 1 = true := by
  decide +kernel

#print axioms TailMarker22.marker_one_eq
#print axioms TailMarker22.base_sound
#print axioms TailMarker22.tail_marker_le
end Spin.Turbo
