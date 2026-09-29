import SpinCodes.Structured.ConcreteNativeSparseRate
import SpinCodes.Structured.ConcreteNativeTotal

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter

theorem native_good_spectrum_bound (m : ℕ) (seed : NativeSeed m) (hg : nativeGood m seed)
    {B : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ∈ weightWindow (bsched m),
      Abar (nativeSeedLaw m) spectrum w ≤ B * ((bsched m).choose w : ℝ)) :
    ∀ w ≤ bsched m, (spectrum seed w : ℝ) ≤
      ((bsched m : ℝ)^2 * B) * ((bsched m).choose w : ℝ) := by
  have h := good_spectrum_bound seed (weightWindow (nativeBlocks m*24)) hB0 hg
    (by simpa only [nativeSeedLaw, native_width] using hB)
  simpa only [native_width] using h

/-- The b² selection penalty is part of the explicit row cost. -/
theorem selected_EZ_sparse_rate (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m))
    (hQ : 4096 ≤ Q) (hQL : Q ≤ Lsched m) (hb : 1000 ≤ bsched m)
    {B ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ∈ weightWindow (bsched m),
      Abar (nativeSeedLaw m) spectrum w ≤ B * ((bsched m).choose w : ℝ))
    (hrow : (2 : ℝ)^(bsched m)*((bsched m : ℝ)^2*B) ≤
      Real.exp ((bsched m : ℝ)*(Real.log 2/2+1281/100000+ε)))
    (hsparse : (Q : ℝ)/(Lsched m : ℝ) ≤ 1/10000) (hε : ε ≤ 1/500) :
    concreteFamily.EZ m Q ≤ Real.exp (-(3/500)*(Q : ℝ)*bsched m) := by
  apply EZ_sparse_rate selectedGood selectedGood_pos m Q hQ hQL hb
    (by positivity) _ hrow hsparse hε
  intro seed hg
  rw [selectedGood_eq m hpos] at hg
  exact native_good_spectrum_bound m seed hg hB0 hB

def spectrumRatioBound (b : ℕ) (r : ℝ) : ℝ :=
  Real.exp ((b : ℝ)*(-Real.log 2/2+1281/100000+r))

def selectedRemainder (m : ℕ) (r : ℝ) : ℝ := 2*Real.log (bsched m : ℝ)/bsched m+r

theorem selected_row_cost (m : ℕ) (r : ℝ) :
    (2 : ℝ)^(bsched m)*((bsched m : ℝ)^2*spectrumRatioBound (bsched m) r) =
      Real.exp ((bsched m : ℝ)*(Real.log 2/2+1281/100000+selectedRemainder m r)) := by
  have hb : (0 : ℝ) < bsched m := by exact_mod_cast bsched_pos m
  have htwo : (2 : ℝ)^(bsched m) = Real.exp ((bsched m : ℝ)*Real.log 2) := by
    rw [Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ) < 2)]
  have hsq : (bsched m : ℝ)^2 = Real.exp (2*Real.log (bsched m : ℝ)) := by
    rw [show (2 : ℝ) = ((2 : ℕ) : ℝ) by norm_num,
      Real.exp_nat_mul, Real.exp_log hb]
  rw [htwo, hsq, spectrumRatioBound, ← Real.exp_add, ← Real.exp_add]
  congr 1
  unfold selectedRemainder
  field_simp
  ring

theorem selectedRemainder_tendsto {r : ℕ → ℝ} (hr : Tendsto r atTop (nhds 0)) :
    Tendsto (fun m => selectedRemainder m (r m)) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have hlog : Tendsto (fun m => Real.log (bsched m : ℝ)/(bsched m : ℝ)) atTop (nhds 0) :=
    Real.isLittleO_log_id_atTop.tendsto_div_nhds_zero.comp hb
  simpa only [selectedRemainder, mul_div_assoc, mul_zero, zero_add] using
    (hlog.const_mul 2).add hr

end Spin.Structured.ConcreteNativeFamily
