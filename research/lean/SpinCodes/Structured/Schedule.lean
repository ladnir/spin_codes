/-
The native-length schedule (`eq:structured-native-schedule`).

    L_m := 128 m,   b_m := least positive multiple of 24 with b ≥ (39/4) log₂(L_m b),
    N_m := L_m b_m, K_m := N_m / 2.

Two things need proof before the theorem can even be non-vacuous: that such a
`b_m` exists at all, and that `b_m → ∞`.  Both are here, and they discharge
five `ScalableCertificate` fields (`N_eq`, `b_pos`, `b_tendsto`, `N_tendsto`,
`L_le_N`).

Existence avoids any numeric bound on `log 2`.  Comparing `(39/4) log₂` against
a linear function would need `(39/4)/log 2 < 24`; instead we use
`log x ≤ 2√x`, which makes the constraint `A√k ≤ 24k` and settles it for
`k ≥ (A/24)²` whatever `A` is.

Index shift: the paper takes `m` positive, so we use `L_m = 128(m+1)` and every
index is usable.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Filter

-- `SchedOK` is a real inequality, so `Nat.find` needs a classical decision.
-- `bsched` is therefore `noncomputable`; nothing downstream evaluates it.
attribute [local instance] Classical.propDecidable

/-- `log x ≤ 2√x` for positive `x`.  Weaker than `log x ≤ x - 1`, but it is
sublinear with an absolute constant, which is what the schedule needs. -/
lemma log_le_two_mul_sqrt {x : ℝ} (hx : 0 < x) : Real.log x ≤ 2 * Real.sqrt x := by
  have hs : 0 < Real.sqrt x := Real.sqrt_pos.mpr hx
  have h1 : Real.log (Real.sqrt x) ≤ Real.sqrt x - 1 := Real.log_le_sub_one_of_pos hs
  have h2 : Real.log (Real.sqrt x) = Real.log x / 2 := Real.log_sqrt hx.le
  rw [h2] at h1
  linarith

/-- Outer length `L_m = 128(m+1)`. -/
def Lsched (m : ℕ) : ℕ := 128 * (m + 1)

lemma Lsched_pos (m : ℕ) : 0 < Lsched m := by unfold Lsched; positivity

/-- The schedule constraint `b ≥ (39/4) log₂(L_m b)`. -/
def SchedOK (m b : ℕ) : Prop :=
  (39 / 4 : ℝ) * Real.logb 2 ((Lsched m : ℝ) * (b : ℝ)) ≤ (b : ℝ)

/-- Some positive multiple of 24 satisfies the schedule constraint. -/
lemma exists_sched (m : ℕ) : ∃ n : ℕ, SchedOK m (24 * (n + 1)) := by
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  set C : ℝ := (Lsched m : ℝ) * 24 with hC
  have hLpos : (0:ℝ) < (Lsched m : ℝ) := by exact_mod_cast Lsched_pos m
  have hCpos : 0 < C := by rw [hC]; linarith
  set A : ℝ := 39 / 2 * Real.sqrt C / Real.log 2 with hA
  have hA0 : 0 ≤ A := by
    rw [hA]; positivity
  obtain ⟨k0, hk0⟩ := exists_nat_ge ((A / 24) ^ 2)
  refine ⟨k0, ?_⟩
  set k : ℝ := (k0 : ℝ) + 1 with hk
  have hk1 : (1:ℝ) ≤ k := by
    rw [hk]
    have : (0:ℝ) ≤ (k0 : ℝ) := Nat.cast_nonneg k0
    linarith
  have hkpos : 0 < k := by linarith
  -- `√k ≥ A/24`
  have hsq : (A / 24) ^ 2 ≤ k := by rw [hk]; linarith
  have hsqrt : A / 24 ≤ Real.sqrt k := by
    have := Real.sqrt_le_sqrt hsq
    rwa [Real.sqrt_sq (by positivity : (0:ℝ) ≤ A / 24)] at this
  have hsk0 : 0 ≤ Real.sqrt k := Real.sqrt_nonneg _
  -- `A √k ≤ 24 k`
  have hAk : A * Real.sqrt k ≤ 24 * k := by
    have h1 : A * Real.sqrt k ≤ (24 * Real.sqrt k) * Real.sqrt k := by
      refine mul_le_mul_of_nonneg_right ?_ hsk0
      linarith
    have h2 : (24 * Real.sqrt k) * Real.sqrt k = 24 * k := by
      rw [mul_assoc, Real.mul_self_sqrt hkpos.le]
    linarith
  -- assemble
  have hCk : ((Lsched m : ℝ) * ((24 * (k0 + 1) : ℕ) : ℝ)) = C * k := by
    rw [hC, hk]; push_cast; ring
  unfold SchedOK
  rw [hCk]
  have hlogb : Real.logb 2 (C * k) ≤ 2 * Real.sqrt (C * k) / Real.log 2 := by
    rw [Real.logb, div_le_div_iff_of_pos_right hlog2]
    exact log_le_two_mul_sqrt (by positivity)
  have hsplit : Real.sqrt (C * k) = Real.sqrt C * Real.sqrt k :=
    Real.sqrt_mul hCpos.le k
  have hfinal : (39 / 4 : ℝ) * Real.logb 2 (C * k) ≤ A * Real.sqrt k := by
    calc (39 / 4 : ℝ) * Real.logb 2 (C * k)
        ≤ (39 / 4 : ℝ) * (2 * Real.sqrt (C * k) / Real.log 2) := by
          exact mul_le_mul_of_nonneg_left hlogb (by norm_num)
      _ = A * Real.sqrt k := by rw [hsplit, hA]; field_simp; ring
  have hcast : ((24 * (k0 + 1) : ℕ) : ℝ) = 24 * k := by rw [hk]; push_cast; ring
  rw [hcast]
  linarith

/-- Block length `b_m`: the schedule's positive multiple of 24. -/
noncomputable def bsched (m : ℕ) : ℕ :=
  24 * (Nat.find (exists_sched m) + 1)

lemma bsched_spec (m : ℕ) : SchedOK m (bsched m) := Nat.find_spec (exists_sched m)

lemma bsched_pos (m : ℕ) : 0 < bsched m := by
  unfold bsched
  omega

lemma bsched_ge_one (m : ℕ) : 1 ≤ bsched m := by
  have := bsched_pos m
  omega

lemma twentyfour_dvd_bsched (m : ℕ) : 24 ∣ bsched m := ⟨_, rfl⟩

/-- Direct length `N_m = L_m b_m`. -/
noncomputable def Nsched (m : ℕ) : ℕ := Lsched m * bsched m

lemma Nsched_eq (m : ℕ) : Nsched m = Lsched m * bsched m := rfl

lemma Lsched_le_Nsched (m : ℕ) : (Lsched m : ℝ) ≤ (Nsched m : ℝ) := by
  have h : Lsched m ≤ Nsched m := by
    rw [Nsched_eq]
    exact Nat.le_mul_of_pos_right _ (bsched_pos m)
  exact_mod_cast h

/-- `b_m ≥ (39/4) log₂ L_m`: the constraint survives dropping `b ≥ 1` inside
the logarithm, which is all that is needed for divergence. -/
lemma bsched_ge_logb (m : ℕ) :
    (39 / 4 : ℝ) * Real.logb 2 (Lsched m : ℝ) ≤ (bsched m : ℝ) := by
  refine le_trans ?_ (bsched_spec m)
  refine mul_le_mul_of_nonneg_left ?_ (by norm_num)
  refine Real.logb_le_logb_of_le (by norm_num) ?_ ?_
  · exact_mod_cast Lsched_pos m
  · have h1 : (1:ℝ) ≤ (bsched m : ℝ) := by exact_mod_cast bsched_ge_one m
    nlinarith [(by exact_mod_cast (Lsched_pos m) : (0:ℝ) < (Lsched m : ℝ))]

lemma Lsched_tendsto : Tendsto (fun m => (Lsched m : ℝ)) atTop atTop := by
  have h : Tendsto Lsched atTop atTop := by
    refine Filter.tendsto_atTop_atTop.mpr fun c => ⟨c, fun a ha => ?_⟩
    unfold Lsched
    omega
  exact tendsto_natCast_atTop_atTop.comp h

lemma bsched_tendsto : Tendsto bsched atTop atTop := by
  have hlb : Tendsto (fun m => (39 / 4 : ℝ) * Real.logb 2 (Lsched m : ℝ)) atTop atTop := by
    have h1 : Tendsto (fun m => Real.logb 2 (Lsched m : ℝ)) atTop atTop :=
      (Real.tendsto_logb_atTop (by norm_num)).comp Lsched_tendsto
    exact Tendsto.const_mul_atTop (by norm_num) h1
  refine tendsto_atTop.mpr fun C => ?_
  have hev := (tendsto_atTop.mp hlb) (C : ℝ)
  filter_upwards [hev] with m hm
  have : (C : ℝ) ≤ (bsched m : ℝ) := le_trans hm (bsched_ge_logb m)
  exact_mod_cast this

lemma Nsched_tendsto : Tendsto (fun m => (Nsched m : ℝ)) atTop atTop :=
  tendsto_atTop_mono (fun m => Lsched_le_Nsched m) Lsched_tendsto

/-! ## Consecutive lengths have relative gap `o(1)`

The paper derives this from `b_m = (39/4) log₂ N_m + O(1)`.  A cheaper route
works and needs no asymptotic characterization of `b_m` at all: `b` is
monotone (a longer outer makes the constraint strictly harder) and each step
raises it by at most one multiple of 24.  Since `b_m → ∞`, the ratio is
squeezed between `1` and `1 + 24/b_m`.
-/

/-- A longer outer makes the schedule constraint harder, so satisfying it at
`m+1` implies satisfying it at `m`. -/
lemma schedOK_of_succ {m b : ℕ} (hb : 0 < b) (h : SchedOK (m + 1) b) : SchedOK m b := by
  unfold SchedOK at h ⊢
  refine le_trans (mul_le_mul_of_nonneg_left ?_ (by norm_num)) h
  refine Real.logb_le_logb_of_le (by norm_num) ?_ ?_
  · have h1 : (0:ℝ) < (Lsched m : ℝ) := by exact_mod_cast Lsched_pos m
    have h2 : (0:ℝ) < (b : ℝ) := by exact_mod_cast hb
    positivity
  · have hL : (Lsched m : ℝ) ≤ (Lsched (m+1) : ℝ) := by
      have : Lsched m ≤ Lsched (m+1) := by unfold Lsched; omega
      exact_mod_cast this
    have h2 : (0:ℝ) ≤ (b : ℝ) := by positivity
    nlinarith

lemma bsched_le_succ (m : ℕ) : bsched m ≤ bsched (m + 1) := by
  unfold bsched
  have hmono : Nat.find (exists_sched m) ≤ Nat.find (exists_sched (m + 1)) :=
    Nat.find_mono (fun n hn => schedOK_of_succ (by positivity) hn)
  omega

lemma bsched_mono : Monotone bsched := monotone_nat_of_le_succ bsched_le_succ

lemma twentyfour_le_bsched (m : ℕ) : 24 ≤ bsched m := by
  unfold bsched; omega

lemma logb_two_four : Real.logb 2 4 = 2 := by
  have h2 : Real.log 2 ≠ 0 := ne_of_gt (Real.log_pos (by norm_num))
  have h4 : (4:ℝ) = 2 ^ (2:ℕ) := by norm_num
  rw [Real.logb, h4, Real.log_pow]
  field_simp
  norm_num

/-- One step of the schedule can always be absorbed by a single extra multiple
of 24: both `L` and `b` at most double, costing `log₂ 4 = 2`, and
`(39/4)·2 = 19.5 ≤ 24`. -/
lemma schedOK_succ_bsched_add (m : ℕ) : SchedOK (m + 1) (bsched m + 24) := by
  have hb := bsched_spec m
  unfold SchedOK at hb ⊢
  have hB24 : (24:ℝ) ≤ (bsched m : ℝ) := by exact_mod_cast twentyfour_le_bsched m
  have hM : (0:ℝ) ≤ (m : ℝ) := by positivity
  have hLm : ((Lsched m : ℕ) : ℝ) = 128 * ((m:ℝ) + 1) := by
    unfold Lsched; push_cast; ring
  have hLm1 : ((Lsched (m+1) : ℕ) : ℝ) = 128 * ((m:ℝ) + 2) := by
    unfold Lsched; push_cast; ring
  have hcast : ((bsched m + 24 : ℕ) : ℝ) = (bsched m : ℝ) + 24 := by push_cast; ring
  have hkey : (Lsched (m+1) : ℝ) * ((bsched m : ℝ) + 24)
      ≤ 4 * ((Lsched m : ℝ) * (bsched m : ℝ)) := by
    rw [hLm, hLm1]; nlinarith [hM, hB24]
  have hprodpos : (0:ℝ) < (Lsched m : ℝ) * (bsched m : ℝ) := by
    rw [hLm]; nlinarith [hM, hB24]
  have hstep : Real.logb 2 ((Lsched (m+1) : ℝ) * ((bsched m : ℝ) + 24))
      ≤ 2 + Real.logb 2 ((Lsched m : ℝ) * (bsched m : ℝ)) := by
    refine le_trans (Real.logb_le_logb_of_le (by norm_num) (by nlinarith) hkey) ?_
    rw [Real.logb_mul (by norm_num) (ne_of_gt hprodpos), logb_two_four]
  rw [hcast]
  linarith

lemma bsched_succ_le (m : ℕ) : bsched (m + 1) ≤ bsched m + 24 := by
  have hwit : SchedOK (m + 1) (24 * ((Nat.find (exists_sched m) + 1) + 1)) := by
    have h := schedOK_succ_bsched_add m
    unfold bsched at h
    have he : 24 * (Nat.find (exists_sched m) + 1) + 24
        = 24 * ((Nat.find (exists_sched m) + 1) + 1) := by omega
    rwa [he] at h
  have hle : Nat.find (exists_sched (m + 1)) ≤ Nat.find (exists_sched m) + 1 :=
    Nat.find_le hwit
  unfold bsched
  omega

lemma bsched_ratio_tendsto_one :
    Tendsto (fun m => (bsched (m + 1) : ℝ) / (bsched m : ℝ)) atTop (nhds 1) := by
  have hpos : ∀ m, (0:ℝ) < (bsched m : ℝ) := fun m => by exact_mod_cast bsched_pos m
  have hlow : ∀ m, (1:ℝ) ≤ (bsched (m + 1) : ℝ) / (bsched m : ℝ) := fun m => by
    rw [le_div_iff₀ (hpos m), one_mul]
    exact_mod_cast bsched_le_succ m
  have hhigh : ∀ m, (bsched (m + 1) : ℝ) / (bsched m : ℝ)
      ≤ 1 + 24 / (bsched m : ℝ) := fun m => by
    have hne : (bsched m : ℝ) ≠ 0 := ne_of_gt (hpos m)
    have h1 : ((bsched (m+1) : ℕ) : ℝ) ≤ (bsched m : ℝ) + 24 := by
      exact_mod_cast bsched_succ_le m
    rw [div_le_iff₀ (hpos m)]
    have h2 : (1 + 24 / (bsched m : ℝ)) * (bsched m : ℝ) = (bsched m : ℝ) + 24 := by
      field_simp
    linarith [h2.le, h2.ge]
  have hmaj : Tendsto (fun m => 1 + 24 / (bsched m : ℝ)) atTop (nhds 1) := by
    have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp bsched_tendsto
    have hinv := hb.inv_tendsto_atTop
    have h24 : Tendsto (fun m => 24 / (bsched m : ℝ)) atTop (nhds 0) := by
      have := hinv.const_mul (24:ℝ)
      rw [mul_zero] at this
      exact this.congr fun m => (div_eq_mul_inv 24 _).symm
    simpa using h24.const_add 1
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hmaj hlow hhigh

lemma Nsched_pos (m : ℕ) : 0 < Nsched m :=
  Nat.mul_pos (Lsched_pos m) (bsched_pos m)

lemma Nsched_mono : Monotone Nsched := by
  refine monotone_nat_of_le_succ fun m => ?_
  unfold Nsched
  exact Nat.mul_le_mul (by unfold Lsched; omega) (bsched_le_succ m)

/-- **Consecutive direct lengths have relative gap `o(1)`.** -/
theorem Nsched_ratio_tendsto_one :
    Tendsto (fun m => (Nsched (m + 1) : ℝ) / (Nsched m : ℝ)) atTop (nhds 1) := by
  have hL : Tendsto (fun m => (Lsched (m + 1) : ℝ) / (Lsched m : ℝ)) atTop (nhds 1) := by
    have hm1 : Tendsto (fun m : ℕ => (m : ℝ) + 1) atTop atTop :=
      tendsto_atTop_add_const_right atTop (1:ℝ) tendsto_natCast_atTop_atTop
    have h1 : Tendsto (fun m : ℕ => 1 + ((m : ℝ) + 1)⁻¹) atTop (nhds 1) := by
      simpa using hm1.inv_tendsto_atTop.const_add 1
    refine h1.congr fun m => ?_
    have hne : ((m:ℝ) + 1) ≠ 0 := by positivity
    have hA : ((Lsched (m+1) : ℕ) : ℝ) = 128 * ((m:ℝ) + 2) := by
      unfold Lsched; push_cast; ring
    have hB : ((Lsched m : ℕ) : ℝ) = 128 * ((m:ℝ) + 1) := by
      unfold Lsched; push_cast; ring
    rw [hA, hB]
    field_simp
    ring
  have hprod := hL.mul bsched_ratio_tendsto_one
  rw [one_mul] at hprod
  refine hprod.congr fun m => ?_
  have hL0 : (0:ℝ) < (Lsched m : ℝ) := by exact_mod_cast Lsched_pos m
  have hb0 : (0:ℝ) < (bsched m : ℝ) := by exact_mod_cast bsched_pos m
  unfold Nsched
  push_cast
  field_simp

/-! ## The requested-length wrapper's rate

`m(n)` is the largest index with `N_{m(n)} ≤ n`, and `N_{m(n)}/n = 1 - o(1)`
follows by squeezing between `N_m/N_{m+1}` and `1`.  With
`Spin.zeroExtend_preserves` (injectivity and minimum distance are preserved),
this is the rate half of `eq:structured-requested-length-wrapper`.
-/

lemma lt_Nsched (m : ℕ) : m < Nsched m := by
  have hb := bsched_pos m
  have h : 128 * (m + 1) ≤ 128 * (m + 1) * bsched m :=
    Nat.le_mul_of_pos_right _ hb
  unfold Nsched Lsched
  omega

/-- The largest schedule index whose direct length fits in `n`. -/
noncomputable def mOf (n : ℕ) : ℕ := Nat.findGreatest (fun m => Nsched m ≤ n) n

lemma Nsched_mOf_le {n : ℕ} (hn : Nsched 0 ≤ n) : Nsched (mOf n) ≤ n :=
  Nat.findGreatest_spec (P := fun m => Nsched m ≤ n) (Nat.zero_le n) hn

lemma lt_Nsched_mOf_succ (n : ℕ) : n < Nsched (mOf n + 1) := by
  by_contra hcon
  push_neg at hcon
  have hlt : mOf n + 1 ≤ n := by
    have := lt_Nsched (mOf n + 1)
    omega
  have hng := Nat.findGreatest_is_greatest (P := fun m => Nsched m ≤ n) (n := n)
    (k := mOf n + 1) (by simp [mOf]) hlt
  exact hng hcon

lemma mOf_tendsto : Tendsto mOf atTop atTop := by
  refine tendsto_atTop.mpr fun M => ?_
  filter_upwards [eventually_ge_atTop (Nsched M)] with n hn
  have hMn : M ≤ n := le_trans (lt_Nsched M).le hn
  exact Nat.le_findGreatest hMn hn

/-- **`N_{m(n)}/n = 1 - o(1)`.** -/
theorem wrapper_rate :
    Tendsto (fun n => (Nsched (mOf n) : ℝ) / (n : ℝ)) atTop (nhds 1) := by
  have hrecip : Tendsto (fun m => (Nsched m : ℝ) / (Nsched (m + 1) : ℝ))
      atTop (nhds 1) := by
    have h := Nsched_ratio_tendsto_one.inv₀ one_ne_zero
    rw [inv_one] at h
    refine h.congr fun m => ?_
    simp [inv_div]
  have hcomp := hrecip.comp mOf_tendsto
  have hNpos : ∀ m, (0:ℝ) < (Nsched m : ℝ) := fun m => by exact_mod_cast Nsched_pos m
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hcomp tendsto_const_nhds ?_ ?_
  · filter_upwards [eventually_ge_atTop (Nsched 0), eventually_ge_atTop 1] with n hn hn1
    have hnpos : (0:ℝ) < (n:ℝ) := by exact_mod_cast hn1
    have hupper : (n:ℝ) ≤ (Nsched (mOf n + 1) : ℝ) := by
      have := (lt_Nsched_mOf_succ n).le
      exact_mod_cast this
    show (Nsched (mOf n) : ℝ) / (Nsched (mOf n + 1) : ℝ) ≤ (Nsched (mOf n) : ℝ) / (n:ℝ)
    exact div_le_div_of_nonneg_left (hNpos _).le hnpos hupper
  · filter_upwards [eventually_ge_atTop (Nsched 0), eventually_ge_atTop 1] with n hn hn1
    have hnpos : (0:ℝ) < (n:ℝ) := by exact_mod_cast hn1
    rw [div_le_one hnpos]
    exact_mod_cast Nsched_mOf_le hn

end Spin.Structured