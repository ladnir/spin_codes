import SpinCodes.UniformMixer.ScalableDenseProfiles
import SpinCodes.Structured.DenseNativeRemainder

/-! Full dense occupation layer, conditional on a checked dense rate cover. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily
  Spin.Structured.DenseOccupationFixed Finset

theorem ScalableDenseRates.layer_bound {η C : ℝ}
    (h : ScalableDenseRates η C) (hC : 0 ≤ C)
    (m Q d : ℕ) (hQ : 0 < Q)
    (seed : NativeSeed m) (hg : nativeGood m seed)
    (hα : (Q : ℝ) / Lsched m ∈ Set.Icc (1 / 10000) 1)
    (hd : (d : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m)) :
    (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      (Spin.Structured.ConcreteRoute.law
        (fun i => (rowSupports seed x i).map
          (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d)) ≤
      densePrefactor (Lsched m) (bsched m) Q C *
        Real.exp (-η * ((Lsched m : ℝ) * bsched m)) := by
  let B : ℝ := outerCost (bsched m) ^ Q *
    (((bsched m : ℝ) + 1) ^ Q *
      ((Lsched m : ℝ) + 1) ^ bsched m) * C *
      Real.exp (-η * ((Lsched m : ℝ) * bsched m))
  have hh := sum_occupation_le_profile_bound seed Q
    (fun x => (Spin.Structured.ConcreteRoute.law
      (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d)) (B := B) (by
      intro S hS w
      have hSn : S.Nonempty := Finset.card_pos.mp (by omega)
      have hαS : (S.card : ℝ) / Lsched m ∈
          Set.Icc (1 / 10000) 1 := by simpa only [hS] using hα
      simpa only [hS, B] using
        h.profile_bound hC m seed hg S hSn hαS w d hd)
  simp only [native_width] at hh
  convert hh using 1
  dsimp [B, densePrefactor]
  rw [two_mul, pow_add]
  ring

#print axioms ScalableDenseRates.layer_bound

end Spin.UniformMixer
