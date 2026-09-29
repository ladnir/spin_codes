import SpinCodes.Structured.SparseRateCounted
import SpinCodes.Structured.Schedule

noncomputable section
namespace Spin.Structured.SparseRate
open Filter

theorem native_log_schedule (m : ℕ) :
    Real.log (Lsched m:ℝ) ≤ (4*Real.log 2/39)*(bsched m:ℝ) := by
  have hh := bsched_ge_logb m
  unfold Real.logb at hh
  have h2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hm := mul_le_mul_of_nonneg_right hh h2.le
  have he : ((39/4:ℝ)*(Real.log (Lsched m:ℝ)/Real.log 2))*Real.log 2 =
      (39/4)*Real.log (Lsched m:ℝ) := by field_simp
  rw [he] at hm
  nlinarith

theorem native_thresholds {ε : ℕ → ℝ} (hε : Tendsto ε atTop (nhds 0)) :
    ∀ᶠ m in atTop, 1000 ≤ bsched m ∧ ε m ≤ 1/500 := by
  have hb : ∀ᶠ m in atTop, 1000 ≤ bsched m := (tendsto_atTop.1 bsched_tendsto) 1000
  have he : ∀ᶠ m in atTop, ε m < 1/500 := hε.eventually (gt_mem_nhds (by norm_num))
  filter_upwards [hb, he] with m hm he
  exact ⟨hm, he.le⟩

/-- A vanishing certified spectral remainder yields a uniform native sparse bound. -/
theorem native_countedMassBound_eventually {ε : ℕ → ℝ}
    (hε : Tendsto ε atTop (nhds 0)) :
    ∀ᶠ m in atTop, ∀ Q R d : ℕ, ∀ α rowCost : ℝ,
      4096 ≤ Q → Q ≤ Lsched m →
      0 < α → α ≤ 1/10000 → (Lsched m:ℝ)*α = Q →
      128*R = Lsched m*bsched m →
      (d:ℝ) ≤ (11/100)*(Lsched m:ℝ)*bsched m →
      0 ≤ rowCost →
      rowCost ≤ Real.exp ((Q:ℝ)*bsched m*(Real.log 2/2+1281/100000+ε m)) →
      (Lsched m|>.choose Q : ℝ) * (rowCost *
        (((1/ConcreteMarked.markMass (Lsched m) Q ((8/5)*α))^(bsched m) *
          (2048*(1-96*α)^R))/((1-(8/5)*α)^d))) ≤
        Real.exp (-(3/500)*(Q:ℝ)*bsched m) := by
  filter_upwards [native_thresholds hε] with m hm
  intro Q R d α rowCost hQ hQL hα0 hα1 hα hrounds hd hrow0 hrow
  exact countedMassBound_le hQ hQL hm.1 hα0 hα1 hα hrounds hd
    (native_log_schedule m) hm.2 hrow0 hrow

end Spin.Structured.SparseRate
