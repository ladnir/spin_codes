import SpinCodes.Structured.ConcreteOuterTailMonotone

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset

/-- The first-stage moment of every sparse Golay message is controlled by its
number of active blocks, with the paper's exact constant. -/
theorem sparse_accumulator_moment_le {b ℓ j : ℕ} (hj : 0<j)
    (hjℓ : 4*j ≤ ℓ) (hℓj : ℓ ≤ 12*j) (hjb : 1000*j ≤ b) :
    (∑ h ∈ range (b+1), Pt b (2*ℓ) h*zeta^h) ≤
      (4*kappa*(j:ℝ)/((b:ℝ)*(1-8*eta)))^(4*j) := by
  have hb : 0<b := by omega
  have hℓ : 0<ℓ := by omega
  have hℓb : 2*ℓ ≤ b := by omega
  have hbR : (0:ℝ)<b := by exact_mod_cast hb
  have hjR : (0:ℝ)<j := by exact_mod_cast hj
  have hℓR : (0:ℝ)<ℓ := by exact_mod_cast hℓ
  have hjℓR : 4*(j:ℝ)≤ℓ := by exact_mod_cast hjℓ
  have hℓjR : (ℓ:ℝ)≤12*j := by exact_mod_cast hℓj
  have hjbR : 1000*(j:ℝ)≤b := by exact_mod_cast hjb
  have hden : 0<(b:ℝ)-2*ℓ := by linarith
  have hk := kappa_pos
  let x : ℝ := (ℓ:ℝ)/b
  let y : ℝ := 4*(j:ℝ)/b
  have hy : 0<y := by dsimp [y]; positivity
  have hyx : y≤x := (div_le_div_iff_of_pos_right hbR).mpr hjℓR
  have hx : x≤3/250 := (div_le_iff₀ hbR).mpr (by linarith)
  have hdx : 0<1-2*x := by linarith
  have hdy : 0<1-2*y := by linarith
  have hAx : 0<kappa*x/(1-2*x) := by have :=hy.trans_le hyx; positivity
  have hAy : 0<kappa*y/(1-2*y) := by positivity
  have hbase : kappa*(ℓ:ℝ)/(b+1-2*ℓ:ℕ) ≤ kappa*x/(1-2*x) := by
    rw [Nat.cast_sub (by omega : 2*ℓ≤b+1)]
    push_cast
    have hid : kappa*x/(1-2*x)=kappa*(ℓ:ℝ)/((b:ℝ)-2*ℓ) := by dsimp [x]; field_simp
    rw [hid]
    exact div_le_div_of_nonneg_left (by positivity) hden (by linarith)
  have hmono : (kappa*x/(1-2*x))^ℓ ≤ (kappa*y/(1-2*y))^(4*j) := by
    have hlog := mul_le_mul_of_nonneg_left (sparse_log_moment_antitone hy hyx hx) hbR.le
    have he₁ : (b:ℝ)*x=ℓ := by dsimp [x]; field_simp
    have he₂ : (b:ℝ)*y=4*j := by dsimp [y]; field_simp
    have hh : (ℓ:ℝ)*Real.log (kappa*x/(1-2*x)) ≤
        (4*j:ℕ)*Real.log (kappa*y/(1-2*y)) := by
      push_cast
      simpa only [←mul_assoc, he₁, he₂] using hlog
    rw [←Real.exp_log (pow_pos hAx ℓ), ←Real.exp_log (pow_pos hAy (4*j))]
    apply Real.exp_le_exp.mpr
    simpa only [Real.log_pow] using hh
  have hlast : kappa*y/(1-2*y) ≤ 4*kappa*(j:ℝ)/((b:ℝ)*(1-8*eta)) := by
    have hyeta : y≤4*eta := by unfold eta; exact (div_le_iff₀ hbR).mpr (by linarith)
    have hconst : 0<1-8*eta := by norm_num [eta]
    calc kappa*y/(1-2*y) ≤ kappa*y/(1-8*eta) :=
        div_le_div_of_nonneg_left (by positivity) hconst (by linarith)
      _ = _ := by dsimp [y]; field_simp
  exact (even_accumulator_moment_le hℓ hℓb).trans
    ((pow_le_pow_left₀ (by positivity) hbase ℓ).trans
      (hmono.trans (pow_le_pow_left₀ hAy.le hlast (4*j))))

end Spin.Structured.ConcreteOuter

