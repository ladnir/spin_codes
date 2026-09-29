import SpinCodes.Structured.ConcreteOuterTailDenseNear
import SpinCodes.Structured.ConcreteOuterEnvelope

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset

def denseLogError (b : ℕ) : ℝ := 48*(Real.log ((b:ℝ)+1)+1)+250

/-- Every supported dense-message path outside the retained window has a uniform gap. -/
theorem dense_pathEntropy_le {k w : ℕ} (hk : 0<k) {p : ℕ×ℕ}
    (hp : p∈supportedPaths k w) (ha : ((k*24:ℕ):ℝ)/125 ≤ p.1)
    (hw : 125*w ≤ 13*(k*24) ∨ 112*(k*24) ≤ 125*w) :
    pathEntropy k w p ≤ -(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24) := by
  obtain ⟨hprod,hcoef,hac,hcw⟩ := mem_filter.mp hp
  obtain ⟨hap,hcp⟩ := mem_product.mp hprod
  obtain ⟨ha0,hab⟩ := mem_Icc.mp hap
  have hcb : p.2 ≤ k*24 := by have := mem_range.mp hcp; omega
  have hb : 0<k*24 := by omega
  have hBR : (1:ℝ) ≤ (k*24:ℕ) := by exact_mod_cast hb
  have hBR0 : (0:ℝ)<(k*24:ℕ) := by positivity
  obtain ⟨hc0,_,_,_⟩ := accT_supported (by omega) hac
  obtain ⟨hf₁,hf₂⟩ := accT_even_feasible ha0
    (golayCoefficient_even hcoef) hac
  let t := clippedOutput (k*24) p.2 w
  obtain ⟨htlo,hthi,htw,hwt⟩ := clippedOutput_spec hc0 hcb hcw
  have h₁ := transitionEntropy_le_shift hb ha0 hab hac hf₁ hf₂ le_rfl (by norm_num)
  have h₂ := transitionEntropy_le_shift hb hc0 hcb hcw htlo hthi htw hwt
  have haB : (p.1:ℝ) ≤ (k*24:ℕ) := by exact_mod_cast hab
  have hg : ((k*24:ℕ):ℝ)*(gBA ((p.1:ℝ)/(k*24:ℕ))+
      piBA ((p.1:ℝ)/(k*24:ℕ)) ((p.2:ℝ)/(k*24:ℕ))+
      piBA ((p.2:ℝ)/(k*24:ℕ)) (t/(k*24:ℕ))) ≤
      -(768/10^10)*((k*24:ℕ):ℝ)+36*(Real.log (k*24:ℕ)+1)+250 := by
    rcases hw with hw|hw
    · apply dense_near_gap hBR ha haB hf₁ hf₂ htlo hthi
      have hwR : 125*(w:ℝ) ≤ 13*((k*24:ℕ):ℝ) := by exact_mod_cast hw
      change t ≤ (w:ℝ) at htw
      linarith
    · have hu₁ : (p.2:ℝ)/2 ≤ ((k*24:ℕ):ℝ)-t := by linarith
      have hu₂ : ((k*24:ℕ):ℝ)-t ≤ ((k*24:ℕ):ℝ)-(p.2:ℝ)/2 := by linarith
      have hwR : 112*((k*24:ℕ):ℝ) ≤ 125*(w:ℝ) := by exact_mod_cast hw
      have hg := dense_near_gap hBR ha haB hf₁ hf₂ hu₁ hu₂
        (show ((k*24:ℕ):ℝ)-t ≤ 13*((k*24:ℕ):ℝ)/125+1/2 by linarith)
      have he : (((k*24:ℕ):ℝ)-t)/(k*24:ℕ)=1-t/(k*24:ℕ) := by field_simp
      rw [he,piBA_reflect] at hg
      exact hg
  have hlog : Real.log ((k*24:ℕ):ℝ) ≤ Real.log (((k*24:ℕ):ℝ)+1) :=
    Real.log_le_log hBR0 (by linarith)
  unfold pathEntropy denseLogError
  change transitionEntropy (k*24) p.2 w ≤ _*piBA _ (t/_)+_ at h₂
  linarith

end Spin.Structured.ConcreteOuter

