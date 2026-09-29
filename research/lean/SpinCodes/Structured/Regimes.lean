/-
The four-regime summation of `thm:structured-spin-scalable`
(paper: "Completing the Argument").

The structured proof bounds `𝔼[Z_{⌊0.11N⌋,Q} | 𝒢_b]` by four different
arguments depending on the occupation `Q`, and then claims their sum over
`1 ≤ Q ≤ L` is `o(1)`.  That last step is where an explicit cutoff matters:
the paper stresses that the bound is "a summable union, not merely a
statement about each occupation sequence".  This file discharges exactly
that step, from the three regime envelopes as hypotheses.

Regimes, in the paper's own order:

* `1 ≤ Q < 4096`  — fixed occupation, continuum kernels.  Each `Q` gets its
  own bound with a `Q`-dependent remainder, so all we may assume is that
  each of the finitely many terms vanishes.
* `4096 ≤ Q ≤ cut` — uniform sparse occupation, `𝔼 ≤ e^{-0.006 Q b}`.
  Repackaged as `ρ_m ^ Q` with `ρ_m = e^{-0.006 b_m} → 0`.
* `cut < Q ≤ L`   — positive occupation, `≤ L e^{-ηN + o(N)}`.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset Filter

/-- Geometric tail above a cutoff: `Σ_{Q ≥ a} ρ^Q ≤ ρ^a / (1 - ρ)`. -/
lemma sum_geom_tail_le {ρ : ℝ} (hρ0 : 0 ≤ ρ) (hρ1 : ρ < 1) (a K : ℕ) :
    ∑ Q ∈ Ico a K, ρ ^ Q ≤ ρ ^ a / (1 - ρ) := by
  have h1ρ : 0 < 1 - ρ := by linarith
  rw [Finset.sum_Ico_eq_sum_range]
  have hfac : ∀ i, ρ ^ (a + i) = ρ ^ a * ρ ^ i := fun i => pow_add ρ a i
  rw [Finset.sum_congr rfl (fun i _ => hfac i), ← Finset.mul_sum]
  refine mul_le_mul_of_nonneg_left ?_ (by positivity)
  have hmul : (∑ i ∈ range (K - a), ρ ^ i) * (1 - ρ) = 1 - ρ ^ (K - a) := by
    linear_combination -geom_sum_mul ρ (K - a)
  rw [← one_div, le_div_iff₀ h1ρ, hmul]
  have : 0 ≤ ρ ^ (K - a) := by positivity
  linarith

/-- **Fixed occupation.**  Finitely many `Q`, each with its own `Q`-dependent
remainder; their sum vanishes. -/
lemma fixed_regime (EZ : ℕ → ℕ → ℝ)
    (h : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => EZ m Q) atTop (nhds 0)) :
    Tendsto (fun m => ∑ Q ∈ Ico 1 4096, EZ m Q) atTop (nhds 0) := by
  have := tendsto_finsetSum (Ico 1 4096) (fun Q hQ => h Q hQ)
  simpa using this

/-- **Uniform sparse occupation.**  A geometric bound with a base tending to
zero, uniform over `4096 ≤ Q ≤ cut m`, sums to `o(1)` no matter how `cut`
grows.  This is the content of the paper's "explicit cutoff avoids any
assumption about how fast `Q` grows". -/
lemma sparse_regime (EZ : ℕ → ℕ → ℝ) (ρ : ℕ → ℝ) (cut : ℕ → ℕ)
    (hnn : ∀ m Q, 0 ≤ EZ m Q)
    (hρ0 : ∀ m, 0 ≤ ρ m) (hρ1 : ∀ m, ρ m < 1)
    (hρ : Tendsto ρ atTop (nhds 0))
    (hbd : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1), EZ m Q ≤ ρ m ^ Q) :
    Tendsto (fun m => ∑ Q ∈ Ico 4096 (cut m + 1), EZ m Q) atTop (nhds 0) := by
  have hlow : ∀ m, 0 ≤ ∑ Q ∈ Ico 4096 (cut m + 1), EZ m Q :=
    fun m => Finset.sum_nonneg fun Q _ => hnn m Q
  have hup : ∀ m, (∑ Q ∈ Ico 4096 (cut m + 1), EZ m Q) ≤ ρ m ^ 4096 / (1 - ρ m) :=
    fun m => le_trans (Finset.sum_le_sum (hbd m))
      (sum_geom_tail_le (hρ0 m) (hρ1 m) 4096 (cut m + 1))
  have hmaj : Tendsto (fun m => ρ m ^ 4096 / (1 - ρ m)) atTop (nhds 0) := by
    have hnum : Tendsto (fun m => ρ m ^ 4096) atTop (nhds 0) := by
      simpa using hρ.pow 4096
    have hden : Tendsto (fun m => 1 - ρ m) atTop (nhds 1) := by
      simpa using (tendsto_const_nhds (x := (1:ℝ)) (f := atTop (α := ℕ))).sub hρ
    have h := hnum.div hden one_ne_zero
    rw [zero_div] at h
    exact h.congr fun m => rfl
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hmaj hlow hup

/-- **Positive occupation.**  `L e^{-ηN + o(N)} → 0`, using only `L ≤ N`.
The schedule gives `L = N / b` with `b = Θ(log N)`, so `L ≤ N` is immediate;
no finer growth information is needed for the sum to vanish. -/
lemma dense_regime (S : ℕ → ℝ) (Lm Nm : ℕ → ℕ) (η : ℝ) (ε : ℕ → ℝ)
    (hη : 0 < η) (hnn : ∀ m, 0 ≤ S m)
    (hL : ∀ m, (Lm m : ℝ) ≤ (Nm m : ℝ))
    (hbd : ∀ m, S m ≤ (Lm m : ℝ) * Real.exp (-η * Nm m + ε m * Nm m))
    (hε : Tendsto ε atTop (nhds 0))
    (hN : Tendsto (fun m => (Nm m : ℝ)) atTop atTop) :
    Tendsto S atTop (nhds 0) := by
  -- eventually `ε ≤ η/2`, so `S ≤ N · exp(-(η/2) N)`
  have hev : ∀ᶠ y : ℝ in nhds (0 : ℝ), y ≤ η / 2 := by
    filter_upwards [Iio_mem_nhds (show (0:ℝ) < η / 2 by linarith)] with y hy
    exact le_of_lt hy
  have hhalf : ∀ᶠ m in atTop, ε m ≤ η / 2 := hε.eventually hev
  have hNpos : ∀ᶠ m in atTop, (0:ℝ) ≤ (Nm m : ℝ) :=
    Eventually.of_forall fun m => by positivity
  have hup : ∀ᶠ m in atTop, S m ≤ (Nm m : ℝ) * Real.exp (-((η / 2) * (Nm m : ℝ))) := by
    filter_upwards [hhalf, hNpos] with m hm hNm
    refine le_trans (hbd m) ?_
    refine mul_le_mul (hL m) (Real.exp_le_exp.mpr ?_) (le_of_lt (Real.exp_pos _)) (by positivity)
    nlinarith [hNm]
  -- `x · e^{-cx} → 0`, by rescaling `x ↦ (η/2) x` in `u · e^{-u} → 0`
  have hmaj : Tendsto (fun m => (Nm m : ℝ) * Real.exp (-((η / 2) * (Nm m : ℝ))))
      atTop (nhds 0) := by
    have hbase : Tendsto (fun u : ℝ => u * Real.exp (-u)) atTop (nhds 0) := by
      simpa using Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 1
    have hscale : Tendsto (fun m => (η / 2) * (Nm m : ℝ)) atTop atTop :=
      Tendsto.const_mul_atTop (by linarith) hN
    have hcomp := hbase.comp hscale
    have h2 := hcomp.const_mul (2 / η)
    rw [mul_zero] at h2
    refine h2.congr fun m => ?_
    have hηne : η ≠ 0 := by linarith
    simp only [Function.comp_apply]
    field_simp
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hmaj
    (Eventually.of_forall hnn) hup

/-- **The union over all occupations.**

Given the three regime envelopes and a cutoff that eventually sits between
`4095` and `L`, the full sum `Σ_{Q=1}^{L} 𝔼[Z_{d,Q} | 𝒢_b]` is `o(1)`. -/
theorem total_regime_sum (EZ : ℕ → ℕ → ℝ) (Lm : ℕ → ℕ) (cut : ℕ → ℕ)
    (hcut : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ Lm m)
    (hfixed : Tendsto (fun m => ∑ Q ∈ Ico 1 4096, EZ m Q) atTop (nhds 0))
    (hsparse : Tendsto (fun m => ∑ Q ∈ Ico 4096 (cut m + 1), EZ m Q) atTop (nhds 0))
    (hdense : Tendsto (fun m => ∑ Q ∈ Ico (cut m + 1) (Lm m + 1), EZ m Q)
      atTop (nhds 0)) :
    Tendsto (fun m => ∑ Q ∈ Ico 1 (Lm m + 1), EZ m Q) atTop (nhds 0) := by
  have hsum := (hfixed.add hsparse).add hdense
  rw [add_zero, add_zero] at hsum
  refine hsum.congr' ?_
  filter_upwards [hcut] with m ⟨h1, h2⟩
  have e1 : (4096 : ℕ) ≤ cut m + 1 := by omega
  have e2 : cut m + 1 ≤ Lm m + 1 := by omega
  have s1 := Finset.sum_Ico_consecutive (EZ m) (by omega : (1:ℕ) ≤ 4096) e1
  have s2 := Finset.sum_Ico_consecutive (EZ m) (by omega : (1:ℕ) ≤ cut m + 1) e2
  rw [← s2, ← s1]


/-! ## Splitting a sum at a threshold

`eq:structured-ba-sparse-sum` bounds `S_b = Σ_{1≤j≤ηb} [C_*(j/b)^3]^j` by
splitting at `j = √b`: below the threshold the summand is at most
`(C_* b^{-3/2})^j`, above it at most `(C_* η^3)^j ≤ 0.657^j`.  Both halves are
then geometric.  This is the abstract form of that argument; it is also the
shape of the regime sums above, so it is stated once here. -/

/-- A sum whose summand is dominated by one geometric series below a threshold
and another above it is bounded by the two geometric tails. -/
theorem sum_split_geom {r q : ℝ} (hr0 : 0 ≤ r) (hr1 : r < 1) (hq0 : 0 ≤ q) (hq1 : q < 1)
    (f : ℕ → ℝ) (J M : ℕ) (hJ1 : 1 ≤ J) (hJM : J ≤ M + 1)
    (hlow : ∀ j ∈ Ico 1 J, f j ≤ r ^ j)
    (hhigh : ∀ j ∈ Ico J (M + 1), f j ≤ q ^ j) :
    ∑ j ∈ Ico 1 (M + 1), f j ≤ r / (1 - r) + q ^ J / (1 - q) := by
  have hsplit := Finset.sum_Ico_consecutive f hJ1 hJM
  rw [← hsplit]
  refine add_le_add ?_ ?_
  · calc ∑ j ∈ Ico 1 J, f j ≤ ∑ j ∈ Ico 1 J, r ^ j := Finset.sum_le_sum hlow
      _ ≤ r ^ 1 / (1 - r) := sum_geom_tail_le hr0 hr1 1 J
      _ = r / (1 - r) := by rw [pow_one]
  · calc ∑ j ∈ Ico J (M + 1), f j ≤ ∑ j ∈ Ico J (M + 1), q ^ j :=
        Finset.sum_le_sum hhigh
      _ ≤ q ^ J / (1 - q) := sum_geom_tail_le hq0 hq1 J (M + 1)

/-- `√ → ∞`.  Mathlib has no such lemma, so it is proved here and reused. -/
theorem tendsto_sqrt_atTop : Tendsto Real.sqrt atTop atTop := by
  refine tendsto_atTop_atTop.mpr fun M => ⟨max M 0 ^ 2, fun x hx => ?_⟩
  have h0 : (0 : ℝ) ≤ max M 0 := le_max_right _ _
  have hle : max M 0 ≤ Real.sqrt x := by
    rw [show max M 0 = Real.sqrt ((max M 0) ^ 2) from (Real.sqrt_sq h0).symm]
    exact Real.sqrt_le_sqrt hx
  exact le_trans (le_max_left _ _) hle

/-- A geometric factor with a `√x` exponent beats every polynomial.  This is
why the high half of the split is `O(b^{-3/2})` even though its base is a
fixed constant below one. -/
theorem pow_mul_rpow_sqrt_tendsto {q : ℝ} (hq0 : 0 < q) (hq1 : q < 1) (k : ℕ) :
    Tendsto (fun x : ℝ => x ^ k * q ^ Real.sqrt x) atTop (nhds 0) := by
  have hc : 0 < -Real.log q := by
    have := Real.log_neg hq0 hq1
    linarith
  set c := -Real.log q with hcdef
  -- rewrite `q ^ √x` as `exp (-c √x)`
  have hrw : ∀ x : ℝ, 0 ≤ x → q ^ Real.sqrt x = Real.exp (-(c * Real.sqrt x)) := by
    intro x hx
    rw [Real.rpow_def_of_pos hq0, hcdef]
    ring_nf
  -- `u ^ (2k) * exp (-c u) → 0`, then substitute `u = √x`
  have hbase : Tendsto (fun u : ℝ => u ^ (2 * k) * Real.exp (-u)) atTop (nhds 0) :=
    Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero (2 * k)
  have hscale : Tendsto (fun u : ℝ => c * u) atTop atTop :=
    Tendsto.const_mul_atTop hc tendsto_id
  have hcomp := hbase.comp hscale
  have hconst := hcomp.const_mul (1 / c ^ (2 * k))
  rw [mul_zero] at hconst
  have hmid : Tendsto (fun u : ℝ => u ^ (2 * k) * Real.exp (-(c * u))) atTop (nhds 0) := by
    refine hconst.congr fun u => ?_
    simp only [Function.comp_apply]
    rw [mul_pow]
    field_simp
  have hfin := hmid.comp tendsto_sqrt_atTop
  refine hfin.congr' ?_
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with x hx
  simp only [Function.comp_apply]
  rw [hrw x hx]
  congr 1
  rw [pow_mul, Real.sq_sqrt hx]


/-! ## `S_b = O(b^{-3/2})`

Both halves of the split vanish faster than `1/b`, which is what the
downstream union needs (`b · S_b → 0`).  The low half is `C/(b√b)` by
construction; the high half is `q^{√b}` with `q = C_*η³ < 0.657` fixed, and
that still beats `1/b` by `pow_mul_rpow_sqrt_tendsto`. -/

theorem tendsto_natSqrt_atTop : Tendsto (fun b : ℕ => Real.sqrt (b : ℝ)) atTop atTop :=
  tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop

/-- The low half: `b · C/(b√b) = C/√b → 0`. -/
theorem mul_div_sqrt_tendsto (C : ℝ) :
    Tendsto (fun b : ℕ => (b : ℝ) * (C / ((b : ℝ) * Real.sqrt (b : ℝ)))) atTop (nhds 0) := by
  have hinv : Tendsto (fun b : ℕ => C / Real.sqrt (b : ℝ)) atTop (nhds 0) := by
    have h := tendsto_natSqrt_atTop.inv_tendsto_atTop.const_mul C
    rw [mul_zero] at h
    exact h.congr fun b => (div_eq_mul_inv C _).symm
  refine hinv.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with b hb
  have hbpos : (0 : ℝ) < (b : ℝ) := by exact_mod_cast hb
  have hspos : (0 : ℝ) < Real.sqrt (b : ℝ) := Real.sqrt_pos.mpr hbpos
  field_simp

/-- The high half: `b · q^{√b} → 0` for a fixed `0 < q < 1`. -/
theorem mul_qpow_sqrt_tendsto {q : ℝ} (hq0 : 0 < q) (hq1 : q < 1) :
    Tendsto (fun b : ℕ => (b : ℝ) * q ^ Real.sqrt (b : ℝ)) atTop (nhds 0) := by
  have h := (pow_mul_rpow_sqrt_tendsto hq0 hq1 1).comp
    (tendsto_natCast_atTop_atTop (R := ℝ))
  refine h.congr fun b => ?_
  simp [Function.comp_apply]

/-- **`b · S_b → 0`.**

The hypothesis is exactly the shape `sum_split_geom` produces, with the low
half's `r/(1-r)` absorbed into `2r` (valid once `r ≤ 1/2`). -/
theorem sparse_sum_mul_tendsto {C q : ℝ} (hq0 : 0 < q) (hq1 : q < 1)
    (S : ℕ → ℝ) (hSnn : ∀ b, 0 ≤ S b)
    (hbound : ∀ᶠ b : ℕ in atTop,
      S b ≤ 2 * (C / ((b : ℝ) * Real.sqrt (b : ℝ))) + q ^ Real.sqrt (b : ℝ) / (1 - q)) :
    Tendsto (fun b : ℕ => (b : ℝ) * S b) atTop (nhds 0) := by
  have h1z : 0 < 1 - q := by linarith
  have hlow : Tendsto (fun b : ℕ => (b : ℝ) * (2 * (C / ((b : ℝ) * Real.sqrt (b : ℝ)))))
      atTop (nhds 0) := by
    have h := (mul_div_sqrt_tendsto C).const_mul 2
    rw [mul_zero] at h
    exact h.congr fun b => by ring
  have hhigh : Tendsto (fun b : ℕ => (b : ℝ) * (q ^ Real.sqrt (b : ℝ) / (1 - q)))
      atTop (nhds 0) := by
    have h := (mul_qpow_sqrt_tendsto hq0 hq1).const_mul (1 / (1 - q))
    rw [mul_zero] at h
    exact h.congr fun b => by field_simp
  have hsum := hlow.add hhigh
  rw [add_zero] at hsum
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hsum ?_ ?_
  · filter_upwards [eventually_ge_atTop 0] with b _
    have : (0 : ℝ) ≤ (b : ℝ) := by positivity
    exact mul_nonneg this (hSnn b)
  · filter_upwards [hbound, eventually_ge_atTop 0] with b hb _
    have hbn : (0 : ℝ) ≤ (b : ℝ) := by positivity
    calc (b : ℝ) * S b
        ≤ (b : ℝ) * (2 * (C / ((b : ℝ) * Real.sqrt (b : ℝ)))
            + q ^ Real.sqrt (b : ℝ) / (1 - q)) := by
          exact mul_le_mul_of_nonneg_left hb hbn
      _ = (b : ℝ) * (2 * (C / ((b : ℝ) * Real.sqrt (b : ℝ))))
            + (b : ℝ) * (q ^ Real.sqrt (b : ℝ) / (1 - q)) := by ring

end Spin.Structured
