import SpinCodes.Structured.Instantiation

noncomputable section
namespace Spin.Structured
open Finset Filter

/-- Finite prefixes do not need the asymptotic sparse envelope. -/
theorem sparse_regime_eventually (EZ : ℕ → ℕ → ℝ) (ρ : ℕ → ℝ) (cut : ℕ → ℕ)
    (hnn : ∀ m Q, 0 ≤ EZ m Q) (hρ0 : ∀ m, 0 ≤ ρ m) (hρ1 : ∀ m, ρ m < 1)
    (hρ : Tendsto ρ atTop (nhds 0))
    (hbd : ∀ᶠ m in atTop, ∀ Q ∈ Ico 4096 (cut m+1), EZ m Q ≤ ρ m^Q) :
    Tendsto (fun m => ∑ Q ∈ Ico 4096 (cut m+1), EZ m Q) atTop (nhds 0) := by
  have h := sparse_regime (fun m Q => min (EZ m Q) (ρ m^Q)) ρ cut
    (fun m Q => le_min (hnn m Q) (pow_nonneg (hρ0 m) Q)) hρ0 hρ1 hρ
    (fun _ _ _ => min_le_right _ _)
  apply h.congr'
  filter_upwards [hbd] with m hm
  exact sum_congr rfl fun Q hQ => min_eq_left (hm Q hQ)

/-- The dense envelope, too, is only an eventual requirement. -/
theorem dense_regime_eventually (S : ℕ → ℝ) (L N : ℕ → ℕ) (η : ℝ) (ε : ℕ → ℝ)
    (hη : 0 < η) (hnn : ∀ m, 0 ≤ S m) (hL : ∀ m, (L m : ℝ) ≤ N m)
    (hbd : ∀ᶠ m in atTop, S m ≤ (L m : ℝ)*Real.exp (-η*N m+ε m*N m))
    (hε : Tendsto ε atTop (nhds 0)) (hN : Tendsto (fun m => (N m : ℝ)) atTop atTop) :
    Tendsto S atTop (nhds 0) := by
  have h := dense_regime (fun m => min (S m) ((L m : ℝ)*Real.exp (-η*N m+ε m*N m)))
    L N η ε hη (fun m => le_min (hnn m) (by positivity)) hL
    (fun _ => min_le_right _ _) hε hN
  apply h.congr'
  filter_upwards [hbd] with m hm
  exact min_eq_left hm

/-- Native-family assembly with eventual regime envelopes and the actual bad event. -/
theorem Family.distance_whp_native_eventually (F : Family)
    (hL : ∀ m, F.L m = Lsched m) (hN : ∀ m, F.N m = Nsched m)
    (hgood : Tendsto F.probNotGood atTop (nhds 0))
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => F.EZ m Q) atTop (nhds 0))
    (cut : ℕ → ℕ) (hcut : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ F.L m)
    (hsparse : ∀ᶠ m in atTop, ∀ Q ∈ Ico 4096 (cut m+1),
      F.EZ m Q ≤ Real.exp (-(0.006 : ℝ)*Q*bsched m))
    (η : ℝ) (hη : 0 < η) (ε : ℕ → ℝ) (hε : Tendsto ε atTop (nhds 0))
    (hdense : ∀ᶠ m in atTop, ∑ Q ∈ Ico (cut m+1) (F.L m+1), F.EZ m Q ≤
      (F.L m : ℝ)*Real.exp (-η*F.N m+ε m*F.N m)) :
    Tendsto F.probBad atTop (nhds 0) := by
  let ρ : ℕ → ℝ := fun m => Real.exp (-(0.006 : ℝ)*(bsched m : ℝ))
  have hρ0 (m : ℕ) : 0 ≤ ρ m := (Real.exp_pos _).le
  have hρ1 (m : ℕ) : ρ m < 1 := by
    apply Real.exp_lt_one_iff.mpr
    have hb : (0 : ℝ) < bsched m := by exact_mod_cast bsched_pos m
    nlinarith
  have hρ : Tendsto ρ atTop (nhds 0) := by
    have hb := (tendsto_natCast_atTop_atTop (R := ℝ)).comp bsched_tendsto
    have hm := Tendsto.const_mul_atTop (by norm_num : (0 : ℝ) < 0.006) hb
    have h := Real.tendsto_exp_atBot.comp (tendsto_neg_atTop_atBot.comp hm)
    apply h.congr
    intro m
    change Real.exp (-(0.006*(bsched m : ℝ))) = Real.exp (-(0.006 : ℝ)*(bsched m : ℝ))
    rw [neg_mul]
  have hsp : Tendsto (fun m => ∑ Q ∈ Ico 4096 (cut m+1), F.EZ m Q) atTop (nhds 0) := by
    apply sparse_regime_eventually F.EZ ρ cut F.EZ_nonneg hρ0 hρ1 hρ
    filter_upwards [hsparse] with m hm
    intro Q hQ
    refine (hm Q hQ).trans_eq ?_
    rw [show ρ m = Real.exp (-(0.006 : ℝ)*(bsched m : ℝ)) from rfl, ← Real.exp_nat_mul]
    congr 1
    ring
  have hde : Tendsto (fun m => ∑ Q ∈ Ico (cut m+1) (F.L m+1), F.EZ m Q) atTop (nhds 0) := by
    apply dense_regime_eventually _ F.L F.N η ε hη
      (fun m => sum_nonneg fun Q _ => F.EZ_nonneg m Q)
      (fun m => by rw [hL, hN]; exact Lsched_le_Nsched m) hdense hε
    simpa only [hN] using Nsched_tendsto
  have ht := total_regime_sum F.EZ F.L cut hcut (fixed_regime F.EZ hfixed) hsp hde
  have hmaj := hgood.add ht
  rw [add_zero] at hmaj
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hmaj F.probBad_nonneg F.probBad_le

end Spin.Structured
