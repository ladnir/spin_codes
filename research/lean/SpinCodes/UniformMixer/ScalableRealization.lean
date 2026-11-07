import SpinCodes.UniformMixer.ScalableEncoderKernel

/-!
# The complete sampled Scalable SPIN encoder weight

The outer word is routed into `bsched m` regions of `Lsched m` bits.
Each region is read in chronological 64-bit groups and receives its
own sampled state maps. The state continues across region boundaries,
as in the paper. This definition uses the same joint setup as
`ScalableSetup`.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily
open Spin.Structured.ConcreteOuter

def scalableRegionSupport (m : ℕ)
    (T : Finset (Fin (Lsched m))) :
    Finset (Fin (64 * scalableStepsPerRegion m)) :=
  T.map (finCongr (scalable_region_length m)).toEmbedding

def scalableRegionInput (m : ℕ)
    (T : Finset (Fin (Lsched m))) :
    Fin (scalableStepsPerRegion m) → Input64 :=
  fun i => regionInput (scalableRegionSupport m T) i

def scalableFlatInput (m : ℕ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    Fin (bsched m * scalableStepsPerRegion m) → Input64 :=
  fun i =>
    let jk := finProdFinEquiv.symm i
    scalableRegionInput m (regions jk.1) jk.2

def scalableFlatMaps (m : ℕ)
    (maps : ScalableStateMaps m) :
    Fin (bsched m * scalableStepsPerRegion m) → StateMap :=
  fun i =>
    let jk := finProdFinEquiv.symm i
    maps jk.1 jk.2

def scalableOutputWeight (m : ℕ)
    (message : Fin (Nsched m / 2) → Bool)
    (seed : ScalableSetupSeed m) : ℕ :=
  let routed := Spin.Structured.ConcreteRoute.routeEval
    (nativeRows m seed.1 message) seed.2.1
  sampledOutputWeight (scalableFlatInput m routed)
    (scalableFlatMaps m seed.2.2) 0

theorem scalable_single_region_moment (m : ℕ) (z : ℝ)
    (T : Finset (Fin (Lsched m))) (q : State22) :
    (piPMF (fun _ : Fin (scalableStepsPerRegion m) =>
      FinPMF.uniform StateMap)).expect
      (fun gs => z ^ sampledOutputWeight (scalableRegionInput m T) gs q) =
    ∑ r : State22,
      chronologicalRegionKernel z (scalableRegionSupport m T) q r := by
  exact sampled_region_moment z (scalableRegionSupport m T) q

end Spin.UniformMixer
