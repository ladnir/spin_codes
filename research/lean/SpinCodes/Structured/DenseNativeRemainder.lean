import SpinCodes.Structured.DenseOccupationOuterAffine
import SpinCodes.Structured.Schedule

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Filter ConcreteOuter

def rowLogCost (b : ℕ) : ℝ := Real.log 2+20*Real.log ((b : ℝ)+1)+13

/-- Profile count, active-set choice, selection cost, route cost, and witness constant. -/
def densePrefactor (L b Q : ℕ) (C : ℝ) : ℝ :=
  (L.choose Q : ℝ)*((b : ℝ)+1)^(2*Q)*outerCost b^Q*((L : ℝ)+1)^b*C

theorem outerCost_le_exp (b : ℕ) : outerCost b ≤ Real.exp (18*Real.log ((b : ℝ)+1)+13) := by
  have he : Real.exp (18*Real.log ((b : ℝ)+1)+13) =
      ((b : ℝ)+1)^2*Real.exp (spectrumLogError b) := by
    rw [show 18*Real.log ((b : ℝ)+1)+13 =
      (2 : ℕ)*Real.log ((b : ℝ)+1)+spectrumLogError b by unfold spectrumLogError; norm_num; ring,
      Real.exp_add, Real.exp_nat_mul, Real.exp_log (by positivity)]
  rw [he, outerCost]
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  nlinarith [Nat.cast_nonneg (α := ℝ) b]

theorem densePrefactor_le_exp {L b Q : ℕ} (hQL : Q ≤ L) {C : ℝ} (hC : 0<C) :
    densePrefactor L b Q C ≤
      Real.exp ((L : ℝ)*rowLogCost b+(b : ℝ)*Real.log ((L : ℝ)+1)+Real.log C) := by
  have hc : (L.choose Q : ℝ) ≤ Real.exp ((L : ℝ)*Real.log 2) := by
    rw [Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ)<2)]
    exact_mod_cast Nat.choose_le_two_pow L Q
  have hp : ((b : ℝ)+1)^(2*Q) ≤ Real.exp ((L : ℝ)*(2*Real.log ((b : ℝ)+1))) := by
    have he : Real.exp ((L : ℝ)*(2*Real.log ((b : ℝ)+1))) = ((b : ℝ)+1)^(2*L) := by
      rw [show (L : ℝ)*(2*Real.log ((b : ℝ)+1)) = (2*L : ℕ)*Real.log ((b : ℝ)+1) by push_cast; ring,
        Real.exp_nat_mul, Real.exp_log (by positivity)]
    rw [he]
    exact pow_le_pow_right₀ (by linarith [Nat.cast_nonneg (α := ℝ) b]) (by omega)
  have ho : outerCost b^Q ≤ Real.exp ((L : ℝ)*(18*Real.log ((b : ℝ)+1)+13)) := by
    have hn : 0 ≤ outerCost b := by unfold outerCost; positivity
    have hl : 0 ≤ Real.log ((b : ℝ)+1) := Real.log_nonneg (by linarith [Nat.cast_nonneg (α := ℝ) b])
    calc
      _ ≤ (Real.exp (18*Real.log ((b : ℝ)+1)+13))^Q := pow_le_pow_left₀ hn (outerCost_le_exp b) Q
      _ = Real.exp ((Q : ℝ)*(18*Real.log ((b : ℝ)+1)+13)) := (Real.exp_nat_mul _ _).symm
      _ ≤ _ := Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_right (by exact_mod_cast hQL) (by positivity))
  have hL : ((L : ℝ)+1)^b = Real.exp ((b : ℝ)*Real.log ((L : ℝ)+1)) := by
    rw [Real.exp_nat_mul, Real.exp_log (by positivity)]
  unfold densePrefactor
  have hoc0 : 0 ≤ outerCost b := by unfold outerCost; positivity
  calc
    _ ≤ Real.exp ((L : ℝ)*Real.log 2)*Real.exp ((L : ℝ)*(2*Real.log ((b : ℝ)+1)))*
        Real.exp ((L : ℝ)*(18*Real.log ((b : ℝ)+1)+13))*((L : ℝ)+1)^b*C := by
      gcongr
    _ = _ := by
      conv_lhs => rw [hL, ← Real.exp_log hC]
      simp only [← Real.exp_add]
      congr 1
      unfold rowLogCost
      ring

theorem log_succ_div_tendsto : Tendsto (fun n : ℕ => Real.log ((n : ℝ)+1)/(n : ℝ)) atTop (nhds 0) := by
  have hn : Tendsto (fun n : ℕ => (n : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hp := tendsto_atTop_add_const_right atTop 1 hn
  have hl := Real.isLittleO_log_id_atTop.tendsto_div_nhds_zero.comp hp
  have hi := hn.const_div_atTop (1 : ℝ)
  have hr : Tendsto (fun n : ℕ => 1+1/(n : ℝ)) atTop (nhds 1) := by
    simpa using tendsto_const_nhds.add hi
  have h := hl.mul hr
  rw [zero_mul] at h
  apply h.congr'
  filter_upwards [eventually_ge_atTop 1] with n hn
  have hne : (n : ℝ) ≠ 0 := by exact_mod_cast (show n≠0 by omega)
  change Real.log ((n : ℝ)+1)/((n : ℝ)+1)*(1+1/(n : ℝ)) = _
  field_simp

def denseRemainder (C : ℝ) (m : ℕ) : ℝ :=
  Real.log C/(Nsched m : ℝ)+rowLogCost (bsched m)/(bsched m : ℝ)+
    Real.log ((Lsched m : ℝ)+1)/(Lsched m : ℝ)

theorem denseRemainder_tendsto (C : ℝ) : Tendsto (denseRemainder C) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have hLnat : Tendsto Lsched atTop atTop := by
    apply tendsto_atTop_mono (fun m => ?_) tendsto_id
    change m ≤ 128*(m+1)
    omega
  have hlog := log_succ_div_tendsto.comp bsched_tendsto
  have hc := hb.const_div_atTop (Real.log 2+13)
  have hr := hc.add (hlog.const_mul 20)
  simp only [mul_zero, add_zero] at hr
  have he := ((Nsched_tendsto.const_div_atTop (Real.log C)).add hr).add
    (log_succ_div_tendsto.comp hLnat)
  simp only [add_zero] at he
  apply he.congr
  intro m
  unfold denseRemainder rowLogCost
  dsimp only [Function.comp_apply]
  ring

theorem native_densePrefactor_le_exp (m : ℕ) {Q : ℕ} (hQL : Q ≤ Lsched m) {C : ℝ} (hC : 0<C) :
    densePrefactor (Lsched m) (bsched m) Q C ≤ Real.exp (denseRemainder C m*(Nsched m : ℝ)) := by
  have h := densePrefactor_le_exp (b := bsched m) hQL hC
  convert h using 1
  congr 1
  have hb : (bsched m : ℝ) ≠ 0 := by exact_mod_cast (bsched_pos m).ne'
  have hL : (Lsched m : ℝ) ≠ 0 := by exact_mod_cast (Lsched_pos m).ne'
  unfold denseRemainder
  rw [Nsched_eq, Nat.cast_mul]
  field_simp
  <;> ring

#print axioms densePrefactor_le_exp
#print axioms denseRemainder_tendsto
#print axioms native_densePrefactor_le_exp

end Spin.Structured.DenseOccupationFixed
