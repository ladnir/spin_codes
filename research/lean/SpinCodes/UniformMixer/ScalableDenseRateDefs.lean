import SpinCodes.UniformMixer.ScalableDensePointRate
import SpinCodes.Structured.DenseOccupationFamilyRateDefs

/-! Interface between a checked dense witness and selected-outer profile counts. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.DenseOccupationFixed

def ScalableDensePointRate (α x slope intercept η C : ℝ) : Prop :=
  ∀ (m d : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))),
    0 < totalWeight rows →
    totalWeight rows < Lsched m * bsched m →
    density rows = α * x →
    (d : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m) →
    ∀ K : ℝ, 0 ≤ K →
      K ≤ Real.exp (((Lsched m : ℝ) * bsched m) *
        α * (slope * x + intercept)) →
      K * ((law rows).prod (scalableStateMapLaw m)).prob
        (fun ω => sampledOutputWeight
          (scalableFlatInput m ω.1) (scalableFlatMaps m ω.2) 0 ≤ d) ≤
        (((bsched m : ℝ) + 1) ^ activeRows rows *
          ((Lsched m : ℝ) + 1) ^ bsched m) * C *
          Real.exp (-η * ((Lsched m : ℝ) * bsched m))

def ScalableDenseRates (η C : ℝ) : Prop :=
  ∀ α x : ℝ,
    α ∈ Set.Icc (1 / 10000) 1 →
    x ∈ Set.Icc (13 / 125) (112 / 125) →
    ∃ p : ℚ × ℚ,
      p ∈ Spin.Majorant.refined.supports ∧
      ScalableDensePointRate α x p.1 p.2 η C

theorem scalable_dense_witness_pointRate
    {α x p y ρ z slope intercept η scale : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hx0 : 0 < x) (hx1 : x ≤ 1)
    (hp0 : 0 < p) (hp1 : p < 1)
    (hy0 : 0 < y) (hy1 : y < 1)
    (hρ : 0 < ρ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (hscale : 0 ≤ scale)
    (hones : ∀ i, (1 : ℝ) ≤ scale * v i)
    (hv : ∀ i,
      (paperSpectrumTransfer (p * y) z).mulVec v i ≤ ρ * v i)
    (hbox : scalableBoxExponent slope intercept p y ρ z α x ≤ -η) :
    ScalableDensePointRate α x slope intercept η (scale * v 0) := by
  intro m d rows hw0 hw1 hden hd K hK0 hK
  exact scalable_dense_counted_rate m d rows hw0 hw1
    hα0 hα1 hx0 hx1 hp0 hp1 hy0 hy1 hρ hz0 hz1
    v hscale hones hv hden hd hK0 hK hbox

theorem ScalableDensePointRate.scaled
    {α x slope intercept η C : ℝ}
    (h : ScalableDensePointRate α x slope intercept η C)
    (m d : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m)))
    (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < Lsched m * bsched m)
    (hden : density rows = α * x)
    (hd : (d : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m))
    {K D : ℝ} (hD : 0 < D) (hK0 : 0 ≤ K)
    (hK : K ≤ D * Real.exp (((Lsched m : ℝ) * bsched m) *
      α * (slope * x + intercept))) :
    K * ((law rows).prod (scalableStateMapLaw m)).prob
      (fun ω => sampledOutputWeight
        (scalableFlatInput m ω.1) (scalableFlatMaps m ω.2) 0 ≤ d) ≤
      D * (((bsched m : ℝ) + 1) ^ activeRows rows *
        ((Lsched m : ℝ) + 1) ^ bsched m) * C *
        Real.exp (-η * ((Lsched m : ℝ) * bsched m)) := by
  have hc : K / D ≤ Real.exp (((Lsched m : ℝ) * bsched m) *
      α * (slope * x + intercept)) :=
    (div_le_iff₀ hD).mpr (by simpa only [mul_comm] using hK)
  have hh := h m d rows hw0 hw1 hden hd (K / D)
    (div_nonneg hK0 hD.le) hc
  have hh' := mul_le_mul_of_nonneg_left hh hD.le
  simpa only [← mul_assoc, mul_div_cancel₀ _ hD.ne'] using hh'

#print axioms scalable_dense_witness_pointRate
#print axioms ScalableDensePointRate.scaled

end Spin.UniformMixer
