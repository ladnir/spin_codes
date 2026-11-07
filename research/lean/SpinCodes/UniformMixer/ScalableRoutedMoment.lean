import SpinCodes.UniformMixer.ScalableCastShuffle

/-!
# First moment after averaging the actual region shuffles

The route first shuffles outer rows, transposes them, and independently
shuffles the resulting regions. The exact encoder moment is the row-law
expectation of a product of averaged full-state regional kernels.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily
  Spin.Structured.ConcreteOuter Spin.Structured.ConcreteRoute
  Spin.Structured.Routing

theorem scalable_route_region_moment (m : ℕ) (z : ℝ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))) :
    (Spin.Structured.ConcreteRoute.law rows).expect (fun regions =>
      ∑ r : State22,
        inputListKernel z (List.ofFn (scalableFlatInput m regions)) 0 r) =
    ((rowLaw rows).map transpose).expect (fun S =>
      ∑ r : State22,
        (List.ofFn (fun j : Fin (bsched m) =>
          Matrix.of (averagedRegionKernel
            (scalableRegionSupport m (S j)) z))).prod 0 r) := by
  rw [Spin.Structured.ConcreteRoute.law, FinPMF.expect_bind]
  congr 1
  funext S
  rw [FinPMF.expect_sum]
  apply Finset.sum_congr rfl
  intro r _
  exact congrFun (congrFun
    (scalable_region_shuffle_product_actual m z S) 0) r

theorem scalable_setup_region_moment (m : ℕ) (z : ℝ)
    (message : Fin (Nsched m / 2) → Bool) :
    (scalableSetupLaw m).expect (fun seed =>
      z ^ scalableOutputWeight m message seed) =
    (nativeSeedLaw m).expect (fun out =>
      ((rowLaw (nativeRows m out message)).map transpose).expect
        (fun S =>
          ∑ r : State22,
            (List.ofFn (fun j : Fin (bsched m) =>
              Matrix.of (averagedRegionKernel
                (scalableRegionSupport m (S j)) z))).prod 0 r)) := by
  rw [scalable_setup_moment]
  congr 1
  funext out
  exact scalable_route_region_moment m z (nativeRows m out message)

end Spin.UniformMixer
