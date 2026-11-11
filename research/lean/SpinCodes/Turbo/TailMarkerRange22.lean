import SpinCodes.Turbo.TailMarkerFracAB22

/-! Range form of the fractional marker tail bound.

For a fixed marker `t > 0`, the marker bound on the `q`-th occupation term is
`F(q) = C(L,q) · (β^q/z^d)^α · (1/C(L,q)/t^q)^64 · W`, where `L = 32 n` and
`W = Σ_k (M_t^(64 n))_{0k}` does not depend on `q`. Binomial coefficients are
log-concave, so `F` is log-convex in `q`, and on any interval `[a, b] ⊆ [0, L]`
it is maximized at an endpoint. Two endpoint certificates therefore bound
every term of the range; this keeps occupations up to `L = 2^19` tractable. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22 TailMarkerFracAB22

namespace TailMarkerRange22

/-! ## A discrete log-convex sequence is maximized at an endpoint -/

theorem step_up (f : ℕ → ℝ) (x : ℕ) (hpos : 0 < f (x + 1)) (hf0 : 0 ≤ f x)
    (hcv : f (x + 1) * f (x + 1) ≤ f x * f (x + 2)) (h : f x ≤ f (x + 1)) :
    f (x + 1) ≤ f (x + 2) := by
  by_contra hc
  push_neg at hc
  have : f x * f (x + 2) < f (x + 1) * f (x + 1) := by
    calc f x * f (x + 2) ≤ f (x + 1) * f (x + 2) :=
          mul_le_mul_of_nonneg_right h (le_of_lt (lt_of_lt_of_le (by
            by_contra h2; push_neg at h2
            exact absurd hcv (by nlinarith)) (le_refl _)))
      _ < f (x + 1) * f (x + 1) := mul_lt_mul_of_pos_left hc hpos
  linarith

theorem le_max_of_logconvex (f : ℕ → ℝ) (a b : ℕ)
    (hpos : ∀ q, a ≤ q → q ≤ b → 0 < f q)
    (hcv : ∀ q, a ≤ q → q + 2 ≤ b → f (q + 1) * f (q + 1) ≤ f q * f (q + 2)) :
    ∀ q, a ≤ q → q ≤ b → f q ≤ max (f a) (f b) := by
  -- increasing from x on: f x ≤ f (x+1) propagates to the right end
  have up : ∀ x, a ≤ x → x + 1 ≤ b → f x ≤ f (x + 1) →
      ∀ m, x + 1 + m ≤ b → f x ≤ f (x + 1 + m) ∧ f (x + m) ≤ f (x + 1 + m) := by
    intro x hx hxb h m
    induction m with
    | zero => intro _; exact ⟨by simpa using h, by simpa using h⟩
    | succ m ih =>
      intro hm
      obtain ⟨h1, h2⟩ := ih (by omega)
      have key : f (x + m + 1) ≤ f (x + m + 2) := by
        have hc := hcv (x + m) (by omega) (by omega)
        have hp := hpos (x + m + 1) (by omega) (by omega)
        have hp0 := (hpos (x + m) (by omega) (by omega)).le
        have h2' : f (x + m) ≤ f (x + m + 1) := by simpa [Nat.add_right_comm] using h2
        have hc' : f (x + m + 1) * f (x + m + 1) ≤ f (x + m) * f (x + m + 2) := by
          simpa [Nat.add_assoc] using hc
        by_contra hk
        push_neg at hk
        have := mul_le_mul_of_nonneg_right h2' (hpos (x + m + 2) (by omega) (by omega)).le
        nlinarith [mul_lt_mul_of_pos_left hk hp]
      refine ⟨?_, ?_⟩
      · calc f x ≤ f (x + 1 + m) := h1
          _ ≤ f (x + 1 + (m + 1)) := by
            have := key; rw [show x + 1 + m = x + m + 1 by omega,
              show x + 1 + (m + 1) = x + m + 2 by omega]; exact this
      · rw [show x + (m + 1) = x + m + 1 by omega, show x + 1 + (m + 1) = x + m + 2 by omega]
        exact key
  -- decreasing up to x: f (x+1) ≤ f x propagates to the left end
  have down : ∀ x, a ≤ x → x + 1 ≤ b → f (x + 1) ≤ f x →
      ∀ m, a + m ≤ x → f x ≤ f (x - m) ∧ f (x - m + 1) ≤ f (x - m) := by
    intro x hx hxb h m
    induction m with
    | zero => intro _; exact ⟨by simp, by simpa using h⟩
    | succ m ih =>
      intro hm
      obtain ⟨h1, h2⟩ := ih (by omega)
      have key : f (x - m) ≤ f (x - (m + 1)) := by
        set y := x - (m + 1) with hy
        have hy1 : x - m = y + 1 := by omega
        have hy2 : x - m + 1 = y + 2 := by omega
        rw [hy1]
        rw [hy2, hy1] at h2
        have hc := hcv y (by omega) (by omega)
        have hp := hpos (y + 1) (by omega) (by omega)
        have hp0 := (hpos y (by omega) (by omega)).le
        by_contra hk
        push_neg at hk
        have := mul_le_mul_of_nonneg_left h2 hp0
        nlinarith [mul_lt_mul_of_pos_left hk hp]
      refine ⟨h1.trans key, ?_⟩
      rw [show x - (m + 1) + 1 = x - m by omega]
      exact key
  intro q hqa hqb
  rcases Nat.eq_or_lt_of_le hqb with hb | hb
  · subst hb; exact le_max_right _ _
  by_cases hinc : f q ≤ f (q + 1)
  · have := (up q hqa hb hinc (b - q - 1) (by omega)).1
    rw [show q + 1 + (b - q - 1) = b by omega] at this
    exact this.trans (le_max_right _ _)
  · push_neg at hinc
    have := (down q hqa hb hinc.le (q - a) (by omega)).1
    rw [show q - (q - a) = a by omega] at this
    exact this.trans (le_max_left _ _)

/-! ## Binomial log-concavity -/

theorem choose_logconcave (L q : ℕ) (hq : q + 2 ≤ L) :
    (Nat.choose L q : ℝ) * (Nat.choose L (q + 2) : ℝ) ≤
      (Nat.choose L (q + 1) : ℝ) * (Nat.choose L (q + 1) : ℝ) := by
  have e1 : (Nat.choose L (q + 1) : ℝ) * (q + 1) = (Nat.choose L q : ℝ) * ((L : ℝ) - q) := by
    have h := Nat.choose_succ_right_eq L q
    have h' : ((Nat.choose L (q + 1) * (q + 1) : ℕ) : ℝ) = ((Nat.choose L q * (L - q) : ℕ) : ℝ) := by
      rw [h]
    push_cast [Nat.cast_sub (show q ≤ L by omega)] at h'
    linarith
  have e2 : (Nat.choose L (q + 2) : ℝ) * (q + 2) =
      (Nat.choose L (q + 1) : ℝ) * ((L : ℝ) - (q + 1)) := by
    have h := Nat.choose_succ_right_eq L (q + 1)
    have h' : ((Nat.choose L (q + 1 + 1) * (q + 1 + 1) : ℕ) : ℝ) =
        ((Nat.choose L (q + 1) * (L - (q + 1)) : ℕ) : ℝ) := by rw [h]
    push_cast [Nat.cast_sub (show q + 1 ≤ L by omega)] at h'
    linarith
  have hL : (q : ℝ) + 2 ≤ L := by exact_mod_cast hq
  have c0 : (0 : ℝ) ≤ Nat.choose L q := Nat.cast_nonneg _
  have c1 : (0 : ℝ) ≤ Nat.choose L (q + 1) := Nat.cast_nonneg _
  have c2 : (0 : ℝ) ≤ Nat.choose L (q + 2) := Nat.cast_nonneg _
  have hpos : (0 : ℝ) < ((L : ℝ) - q) * ((q : ℝ) + 2) := by
    apply mul_pos <;> linarith
  -- C0 C2 (L-q)(q+2) = C1 (q+1) · C1 (L-q-1) ≤ C1^2 (L-q)(q+2)
  have key : (Nat.choose L q : ℝ) * (Nat.choose L (q + 2) : ℝ) * (((L : ℝ) - q) * ((q : ℝ) + 2)) ≤
      (Nat.choose L (q + 1) : ℝ) * (Nat.choose L (q + 1) : ℝ) * (((L : ℝ) - q) * ((q : ℝ) + 2)) := by
    have lhs : (Nat.choose L q : ℝ) * (Nat.choose L (q + 2) : ℝ) * (((L : ℝ) - q) * ((q : ℝ) + 2)) =
        ((Nat.choose L q : ℝ) * ((L : ℝ) - q)) * ((Nat.choose L (q + 2) : ℝ) * ((q : ℝ) + 2)) := by ring
    rw [lhs, ← e1, e2]
    have : ((q : ℝ) + 1) * ((L : ℝ) - (q + 1)) ≤ ((L : ℝ) - q) * ((q : ℝ) + 2) := by nlinarith
    nlinarith [mul_nonneg c1 c1]
  exact le_of_mul_le_mul_right key hpos

/-! ## The marker bound as a function of `q` -/

/-- `W = Σ_k (M_t^(64 n))_{0k}` (independent of `q`). -/
def markerMass (ea eb n : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (t : ℝ) : ℝ :=
  ∑ k : Fin 11, ((orderedMatrixMarker22 pairedProfileDegree22
    (pairedProfileMatrix22 B ((ea : ℝ) / eb)) t) ^ (n * 64)) 0 k

/-- `F(q) = C(L,q) · (β^q/z^d)^α · (1/C(L,q)/t^q)^64 · W`. -/
def markerF (ea eb n d : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (z t : ℝ) (q : ℕ) : ℝ :=
  (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
    ((1 / (Nat.choose (32 * n) q : ℝ) / t ^ q) ^ 64 * markerMass ea eb n B t))

theorem markerM_nonneg (ea eb : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (t : ℝ) (ht : 0 ≤ t) (i j : Fin 11) :
    0 ≤ orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B ((ea : ℝ) / eb)) t i j := by
  unfold orderedMatrixMarker22
  rw [Matrix.sum_apply]
  exact Finset.sum_nonneg fun ab _ => by
    rw [Matrix.smul_apply, smul_eq_mul]
    exact mul_nonneg (pow_nonneg ht _) (PairedFracChain22.pairedProfile_nonneg B hB _ ab _ _)

theorem markerMass_nonneg (ea eb n : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (t : ℝ) (ht : 0 ≤ t) : 0 ≤ markerMass ea eb n B t :=
  Finset.sum_nonneg fun k _ => markerPowNonneg22 _ (markerM_nonneg ea eb B hB t ht) _ 0 k

theorem beta_pos : 0 < turboOldUniformFactor22 := by
  unfold turboOldUniformFactor22; positivity

/-- The actual term is at most `F(q)`. -/
theorem term_le_markerF (ea eb n d : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (z t : ℝ) (hz : 0 < z) (ht : 0 < t) (q : ℕ) :
    (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤
      markerF ea eb n d B z t q := by
  unfold markerF markerMass
  apply mul_le_mul_of_nonneg_left _ (Nat.cast_nonneg _)
  apply mul_le_mul_of_nonneg_left _ (Real.rpow_nonneg (div_nonneg
    (pow_nonneg beta_pos.le _) (pow_nonneg hz.le _)) _)
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun k _ => pairedRegionDP22_marker_pow64_le B _ hB n q t ht 0 k

/-- Closed form for `0 < C(L,q)`. -/
theorem markerF_eq (ea eb n d : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (z t : ℝ)
    (q : ℕ) (hC : (0 : ℝ) < Nat.choose (32 * n) q) (ht : 0 < t) :
    markerF ea eb n d B z t q = (turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      markerMass ea eb n B t / ((Nat.choose (32 * n) q : ℝ) ^ 63 * (t ^ q) ^ 64) := by
  unfold markerF
  have h1 : (Nat.choose (32 * n) q : ℝ) ≠ 0 := hC.ne'
  have h2 : t ^ q ≠ 0 := pow_ne_zero _ ht.ne'
  rw [div_pow, div_pow, one_pow]
  field_simp

theorem markerF_logconvex (ea eb n d : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (z t : ℝ) (hz : 0 < z) (ht : 0 < t) (q : ℕ)
    (hq : q + 2 ≤ 32 * n) :
    markerF ea eb n d B z t (q + 1) * markerF ea eb n d B z t (q + 1) ≤
      markerF ea eb n d B z t q * markerF ea eb n d B z t (q + 2) := by
  have hCpos : ∀ r, r ≤ 32 * n → (0 : ℝ) < Nat.choose (32 * n) r := fun r hr => by
    exact_mod_cast Nat.choose_pos hr
  rw [markerF_eq _ _ _ _ _ _ _ _ (hCpos q (by omega)) ht,
    markerF_eq _ _ _ _ _ _ _ _ (hCpos (q + 1) (by omega)) ht,
    markerF_eq _ _ _ _ _ _ _ _ (hCpos (q + 2) (by omega)) ht]
  set W := markerMass ea eb n B t
  have hW : 0 ≤ W := markerMass_nonneg ea eb n B hB t ht.le
  set α : ℝ := (ea : ℝ) / eb
  set X : ℕ → ℝ := fun r => (turboOldUniformFactor22 ^ r / z ^ d) ^ α with hX
  have hXr : ∀ r, 0 ≤ turboOldUniformFactor22 ^ r / z ^ d := fun r =>
    div_nonneg (pow_nonneg beta_pos.le _) (pow_nonneg hz.le _)
  have hXmul : X (q + 1) * X (q + 1) = X q * X (q + 2) := by
    simp only [hX]
    rw [← Real.mul_rpow (hXr _) (hXr _), ← Real.mul_rpow (hXr _) (hXr _)]
    congr 1
    have hz' : z ^ d ≠ 0 := pow_ne_zero _ hz.ne'
    field_simp
    ring
  set C0 := (Nat.choose (32 * n) q : ℝ)
  set C1 := (Nat.choose (32 * n) (q + 1) : ℝ)
  set C2 := (Nat.choose (32 * n) (q + 2) : ℝ)
  have c0 : 0 < C0 := hCpos q (by omega)
  have c1 : 0 < C1 := hCpos (q + 1) (by omega)
  have c2 : 0 < C2 := hCpos (q + 2) (by omega)
  have hlc : C0 * C2 ≤ C1 * C1 := choose_logconcave (32 * n) q hq
  have ht' : ∀ r, 0 < (t ^ r) ^ 64 := fun r => pow_pos (pow_pos ht _) _
  have htmul : (t ^ (q + 1)) ^ 64 * (t ^ (q + 1)) ^ 64 = (t ^ q) ^ 64 * (t ^ (q + 2)) ^ 64 := by
    rw [← pow_mul, ← pow_mul, ← pow_mul, ← pow_add, ← pow_add]; congr 1; ring
  have hC63 : (C0 * C2) ^ 63 ≤ (C1 * C1) ^ 63 :=
    pow_le_pow_left₀ (mul_pos c0 c2).le hlc 63
  have hX0 : ∀ r, 0 ≤ X r := fun r => Real.rpow_nonneg (hXr r) _
  rw [div_mul_div_comm, div_mul_div_comm]
  apply div_le_div₀ (mul_nonneg (mul_nonneg (hX0 q) hW) (mul_nonneg (hX0 _) hW))
  · show X (q + 1) * W * (X (q + 1) * W) ≤ X q * W * (X (q + 2) * W)
    have : X (q + 1) * W * (X (q + 1) * W) = X q * W * (X (q + 2) * W) := by
      calc X (q + 1) * W * (X (q + 1) * W) = (X (q + 1) * X (q + 1)) * (W * W) := by ring
        _ = (X q * X (q + 2)) * (W * W) := by rw [hXmul]
        _ = X q * W * (X (q + 2) * W) := by ring
    exact this.le
  · exact mul_pos (mul_pos (pow_pos c0 _) (ht' _)) (mul_pos (pow_pos c2 _) (ht' _))
  · calc C0 ^ 63 * (t ^ q) ^ 64 * (C2 ^ 63 * (t ^ (q + 2)) ^ 64)
        = (C0 * C2) ^ 63 * ((t ^ q) ^ 64 * (t ^ (q + 2)) ^ 64) := by ring
      _ ≤ (C1 * C1) ^ 63 * ((t ^ q) ^ 64 * (t ^ (q + 2)) ^ 64) :=
          mul_le_mul_of_nonneg_right hC63 (mul_pos (ht' _) (ht' _)).le
      _ = C1 ^ 63 * (t ^ (q + 1)) ^ 64 * (C1 ^ 63 * (t ^ (q + 1)) ^ 64) := by
          rw [← htmul]; ring

/-- Every term of a range is at most the larger endpoint marker bound. -/
theorem term_le_endpoints (ea eb n d : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (z t : ℝ) (hz : 0 < z) (ht : 0 < t) (a b : ℕ)
    (hb : b ≤ 32 * n) (q : ℕ) (hqa : a ≤ q) (hqb : q ≤ b) :
    (Nat.choose (32 * n) q : ℝ) * ((turboOldUniformFactor22 ^ q / z ^ d) ^ ((ea : ℝ) / eb) *
      ∑ k : Fin 11, ((pairedRegionDP22 B ((ea : ℝ) / eb) n q) ^ 64) 0 k) ≤
      max (markerF ea eb n d B z t a) (markerF ea eb n d B z t b) := by
  refine (term_le_markerF ea eb n d B hB z t hz ht q).trans ?_
  set W := markerMass ea eb n B t
  rcases (markerMass_nonneg ea eb n B hB t ht.le).eq_or_lt with hW0 | hWpos
  · -- W = 0: every F vanishes
    have hz0 : ∀ r, markerF ea eb n d B z t r = 0 := fun r => by
      unfold markerF; rw [← hW0]; simp
    rw [hz0 q, hz0 a, hz0 b, max_self]
  have hFpos : ∀ r, r ≤ 32 * n → 0 < markerF ea eb n d B z t r := fun r hr => by
    have hC : (0 : ℝ) < Nat.choose (32 * n) r := by exact_mod_cast Nat.choose_pos hr
    rw [markerF_eq _ _ _ _ _ _ _ _ hC ht]
    apply div_pos (mul_pos (Real.rpow_pos_of_pos (div_pos (pow_pos beta_pos _)
      (pow_pos hz _)) _) hWpos)
    exact mul_pos (pow_pos hC _) (pow_pos (pow_pos ht _) _)
  exact le_max_of_logconvex (markerF ea eb n d B z t) a b
    (fun r _ hr => hFpos r (by omega))
    (fun r _ hr => markerF_logconvex ea eb n d B hB z t hz ht r (by omega)) q hqa hqb

/-- Endpoint certificate: the integer checks of `tail_marker_fracAB_le`
bound `F(q)` itself. -/
theorem markerF_le (ea eb : ℕ) (heb : eb ≠ 0) (n d zb : ℕ)
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
    markerF ea eb n d B z (((tn0 * 2 ^ kp : ℤ) : ℝ) / ((td0 * 2 ^ kn : ℤ) : ℝ)) q ≤ (a : ℝ) / 2 ^ Kt := by
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
  rw [← htdef]
  unfold markerF markerMass
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
  have hM0' : 0 ≤ ∑ k : Fin 11, (M ^ (n * 64)) 0 k :=
    Finset.sum_nonneg fun k _ => markerPowNonneg22 _ hM0 _ 0 k
  calc (Cq : ℝ) * (X ^ ((ea : ℝ) / eb) * (c ^ 64 * ∑ k : Fin 11, (M ^ (n * 64)) 0 k))
      ≤ (Cq : ℝ) * (((yn : ℝ) * (2 : ℝ) ^ G) * (c ^ 64 * ((2 : ℝ) ^ E * (sumU : ℝ) / s))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul hXp (mul_le_mul_of_nonneg_left hmass
          (pow_nonneg hc0 _)) (mul_nonneg (pow_nonneg hc0 _) hM0') hY0) hCR.le
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

end TailMarkerRange22

#print axioms TailMarkerRange22.term_le_endpoints
#print axioms TailMarkerRange22.markerF_le
end Spin.Turbo
