import SpinCodes.Structured.ConcreteOuterTailSumTerm

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset Filter

/-- The summation uses a summable geometric majorant, without choosing an
asymptotic cutoff depending on the message. -/
def sparseScaledTerm (b n : ℕ) : ℝ :=
  if 1000*(n+1) ≤ b then (b:ℝ)*(Cstar*(((n+1:ℕ):ℝ)/b)^3)^(n+1) else 0

def sparseBoundSeries (b : ℕ) : ℝ := ∑ j ∈ Ico 1 (b/1000+1), (Cstar*((j:ℝ)/b)^3)^j

lemma Cstar_nonneg : 0≤Cstar := by unfold Cstar; positivity

lemma sparseScaledTerm_nonneg (b n : ℕ) : 0 ≤ sparseScaledTerm b n := by
  unfold sparseScaledTerm
  split_ifs
  · exact mul_nonneg (by positivity) (pow_nonneg (mul_nonneg Cstar_nonneg (by positivity)) _)
  · exact le_refl _

lemma sparseScaledTerm_bound (b n : ℕ) :
    sparseScaledTerm b n ≤ Cstar*eta^2*((n:ℝ)+1)*(657/1000:ℝ)^n := by
  have hC := Cstar_nonneg
  have hη : 0≤eta := by norm_num [eta]
  unfold sparseScaledTerm
  split_ifs with h
  · have hb : (0:ℝ)<b := by exact_mod_cast (by omega : 0<b)
    have hj : (0:ℝ)≤((n+1:ℕ):ℝ)/b := by positivity
    have he : ((n+1:ℕ):ℝ)/b ≤ eta := by
      unfold eta
      apply (div_le_iff₀ hb).mpr
      have hh : (1000:ℝ)*((n+1:ℕ):ℝ)≤b := by exact_mod_cast h
      linarith
    have hbase : Cstar*(((n+1:ℕ):ℝ)/b)^3 ≤ (657/1000:ℝ) :=
      (mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hj he 3) hC).trans Cstar_eta_cubed_lt.le
    have hpow := pow_le_pow_left₀ (mul_nonneg hC (pow_nonneg hj _)) hbase n
    have hsq := pow_le_pow_left₀ hj he 2
    have hid : (b:ℝ)*(Cstar*(((n+1:ℕ):ℝ)/b)^3)^(n+1) =
        (Cstar*((n+1:ℕ):ℝ)*(((n+1:ℕ):ℝ)/b)^2)*(Cstar*(((n+1:ℕ):ℝ)/b)^3)^n := by
      rw [pow_succ]
      field_simp
    rw [hid]
    have hh := mul_le_mul (mul_le_mul_of_nonneg_left hsq (by positivity : 0≤Cstar*((n+1:ℕ):ℝ)))
      hpow (by positivity) (by positivity)
    convert hh using 1 <;> push_cast <;> ring
  · positivity

lemma sparseScaledTerm_tendsto (n : ℕ) :
    Tendsto (fun b : ℕ => sparseScaledTerm b n) atTop (nhds 0) := by
  have hb : Tendsto (fun b : ℕ => (b:ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hdiv := hb.const_div_atTop (((n+1:ℕ):ℝ))
  have ht := (((hdiv.pow 2).const_mul (Cstar*((n+1:ℕ):ℝ))).mul (((hdiv.pow 3).const_mul Cstar).pow n))
  simp only [zero_pow (by decide : 2≠0), zero_pow (by decide : 3≠0), mul_zero, zero_mul] at ht
  apply ht.congr'
  filter_upwards [eventually_ge_atTop (1000*(n+1))] with b h
  rw [sparseScaledTerm, if_pos h]
  have hb0 : (b:ℝ)≠0 := by exact_mod_cast (by omega : b≠0)
  have hid : (b:ℝ)*(Cstar*(((n+1:ℕ):ℝ)/b)^3) =
      Cstar*((n+1:ℕ):ℝ)*(((n+1:ℕ):ℝ)/b)^2 := by field_simp
  rw [←hid, pow_succ _ n]
  ring

/-- The sparse union, even multiplied by b for the upper-tail prefactor, vanishes. -/
theorem sparseScaledSeries_tendsto :
    Tendsto (fun b : ℕ => ∑' n, sparseScaledTerm b n) atTop (nhds 0) := by
  have hq : ‖(657/1000:ℝ)‖<1 := by norm_num
  have hs₁ := summable_pow_mul_geometric_of_norm_lt_one 1 hq
  have hs₀ := summable_geometric_of_norm_lt_one hq
  have hs : Summable (fun n : ℕ => Cstar*eta^2*((n:ℝ)+1)*(657/1000:ℝ)^n) := by
    convert (hs₁.add hs₀).mul_left (Cstar*eta^2) using 1
    funext n
    simp only [pow_one]
    ring
  have ht := tendsto_tsum_of_dominated_convergence hs sparseScaledTerm_tendsto
    (Eventually.of_forall (fun b n => by
      rw [Real.norm_eq_abs, abs_of_nonneg (sparseScaledTerm_nonneg b n)]
      exact sparseScaledTerm_bound b n))
  simpa using ht

end Spin.Structured.ConcreteOuter


