import SpinCodes.Structured.ConcreteOuterMajorantSpectrum
import SpinCodes.Structured.ConcreteOuterDefect
import SpinCodes.Structured.ConcreteNativeAsymptotic

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Filter

/-- The paper's numerical shell bound, now with the exact finite remainder. -/
theorem expected_spectrum_uniform {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ)) ≤
      Real.exp (((k*24:ℕ):ℝ)*(-Real.log 2/2+1281/100000)+shellLogError (k*24)) *
        ((k*24).choose w:ℝ) := by
  refine (expected_spectrum_shell_refined hk hw).trans ?_
  apply mul_le_mul_of_nonneg_right _ (by positivity)
  apply Real.exp_le_exp.mpr
  have h := mul_le_mul_of_nonneg_left (refined_defect hw) (show (0:ℝ)≤(k*24:ℕ) by positivity)
  linarith

theorem shellLogError_div_tendsto :
    Tendsto (fun b : ℕ => shellLogError b/(b:ℝ)) atTop (nhds 0) := by
  have hb : Tendsto (fun b : ℕ => (b:ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hb' : Tendsto (fun b : ℕ => (b:ℝ)+1) atTop atTop :=
    tendsto_atTop_add_const_right atTop 1 hb
  have hlog := Real.isLittleO_log_id_atTop.tendsto_div_nhds_zero.comp hb'
  have hinv := hb.const_div_atTop (1:ℝ)
  have hratio : Tendsto (fun b : ℕ => 1+1/(b:ℝ)) atTop (nhds 1) := by
    simpa using tendsto_const_nhds.add hinv
  have hprod := hlog.mul hratio
  have hquot : Tendsto (fun b : ℕ => Real.log ((b:ℝ)+1)/(b:ℝ)) atTop (nhds 0) := by
    simp only [zero_mul] at hprod
    apply hprod.congr'
    filter_upwards [eventually_ge_atTop 1] with b hb
    have hb0 : (b:ℝ)≠0 := by exact_mod_cast (show b≠0 by omega)
    change Real.log ((b:ℝ)+1)/((b:ℝ)+1)*(1+1/(b:ℝ)) = Real.log ((b:ℝ)+1)/(b:ℝ)
    field_simp
  have hconst := hb.const_div_atTop (13:ℝ)
  have hsum := (hquot.const_mul 17).add hconst
  simp only [mul_zero, zero_add] at hsum
  convert hsum using 1
  funext b
  unfold shellLogError
  ring

end Spin.Structured.ConcreteOuter

namespace Spin.Structured.ConcreteNativeFamily
open Spin.Numeric Filter ConcreteOuter Finset

def nativeShellRemainder (m : ℕ) : ℝ := shellLogError (bsched m)/(bsched m:ℝ)

theorem nativeShellRemainder_tendsto : Tendsto nativeShellRemainder atTop (nhds 0) :=
  shellLogError_div_tendsto.comp bsched_tendsto

/-- The concrete native outer satisfies the expected spectral envelope at every
schedule index, with a uniform explicitly vanishing remainder. -/
theorem native_expected_spectrum_uniform (m : ℕ) (w : ℕ)
    (hw : w ∈ weightWindow (bsched m)) :
    Abar (nativeSeedLaw m) spectrum w ≤
      spectrumRatioBound (bsched m) (nativeShellRemainder m) * ((bsched m).choose w:ℝ) := by
  have hk : 0 < nativeBlocks m := by
    have hh := bsched_pos m
    rw [←native_width m] at hh
    omega
  have hb : (0:ℝ)<bsched m := by exact_mod_cast bsched_pos m
  obtain ⟨_,hl,hu⟩ := mem_filter.mp hw
  have hlR : (13:ℝ)*bsched m ≤ 125*(w:ℝ) := by exact_mod_cast hl
  have huR : (125:ℝ)*w ≤ 112*(bsched m:ℝ) := by exact_mod_cast hu
  have hwR : (w:ℝ)/(nativeBlocks m*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125) := by
    rw [native_width]
    exact ⟨(le_div_iff₀ hb).mpr (by linarith), (div_le_iff₀ hb).mpr (by linarith)⟩
  have hh := expected_spectrum_uniform hk hwR
  simp only [native_width] at hh
  have he : (bsched m:ℝ)*(-Real.log 2/2+1281/100000+nativeShellRemainder m) =
      (bsched m:ℝ)*(-Real.log 2/2+1281/100000)+shellLogError (bsched m) := by
    unfold nativeShellRemainder
    field_simp
  simpa only [Abar, nativeSeedLaw, spectrumRatioBound, he] using hh

/-- After the envelope closure the native sparse regime has only the outside-window
tail limit as an outstanding outer premise. -/
theorem concrete_sparse_eventually_of_tail (htail : Tendsto nativeTail atTop (nhds 0)) :
    ∀ᶠ m in atTop, ∀ Q ∈ Ico 4096 (nativeCut m+1),
      concreteFamily.EZ m Q ≤ Real.exp (-(0.006:ℝ)*(Q:ℝ)*bsched m) :=
  concrete_sparse_eventually htail nativeShellRemainder_tendsto
    (Eventually.of_forall fun m w hw => native_expected_spectrum_uniform m w hw)

end Spin.Structured.ConcreteNativeFamily

