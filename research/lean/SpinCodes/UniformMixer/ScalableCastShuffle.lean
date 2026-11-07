import SpinCodes.UniformMixer.ScalableRouteAverage

/-!
# Regional shuffle after the length cast

The Scalable schedule writes a region length as `128(m+1)`, while the
regional encoder writes it as `64 · 2(m+1)`. Reindexing by this exact
equality preserves the uniform fixed-weight shuffle law.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Routing

attribute [local instance] Classical.propDecidable

theorem shuffleLaw_map_finsetCongr {n m : ℕ}
    (e : Fin n ≃ Fin m) (S : Finset (Fin n)) :
    (shuffleLaw S).map (Equiv.finsetCongr e) =
      shuffleLaw (Equiv.finsetCongr e S) := by
  have hnm : n = m := by
    simpa only [Fintype.card_fin] using Fintype.card_congr e
  subst m
  apply FinPMF.ext
  funext T
  rw [FinPMF.map_equiv_apply]
  have hT : ((Equiv.finsetCongr e).symm T).card = T.card := by
    simp [Equiv.finsetCongr_symm, Equiv.finsetCongr_apply]
  have hS : (Equiv.finsetCongr e S).card = S.card := by
    simp [Equiv.finsetCongr_apply]
  simp only [shuffleLaw_apply, hT, hS]

theorem scalableRegionSupport_eq_congr (m : ℕ)
    (T : Finset (Fin (Spin.Structured.Lsched m))) :
    scalableRegionSupport m T =
      Equiv.finsetCongr (finCongr (scalable_region_length m)) T := by
  simp [scalableRegionSupport, Equiv.finsetCongr_apply]

theorem scalableRegionSupport_shuffleLaw (m : ℕ)
    (S : Finset (Fin (Spin.Structured.Lsched m))) :
    (shuffleLaw S).map (scalableRegionSupport m) =
      shuffleLaw (scalableRegionSupport m S) := by
  have he : scalableRegionSupport m =
      Equiv.finsetCongr (finCongr (scalable_region_length m)) := by
    funext T
    exact scalableRegionSupport_eq_congr m T
  rw [he]
  exact shuffleLaw_map_finsetCongr _ S

theorem scalableAveragedRegion_eq (m : ℕ) (z : ℝ)
    (S : Finset (Fin (Spin.Structured.Lsched m))) :
    scalableAveragedRegion m z S =
      Matrix.of (averagedRegionKernel (scalableRegionSupport m S) z) := by
  ext q r
  change (shuffleLaw S).expect (fun T =>
      chronologicalRegionKernel z (scalableRegionSupport m T) q r) =
    (shuffleLaw (scalableRegionSupport m S)).expect
      (fun T => chronologicalRegionKernel z T q r)
  calc
    _ = ((shuffleLaw S).map (scalableRegionSupport m)).expect
          (fun T => chronologicalRegionKernel z T q r) := by
        rw [FinPMF.expect_map]
    _ = _ := by rw [scalableRegionSupport_shuffleLaw]

theorem scalable_region_shuffle_product_actual (m : ℕ) (z : ℝ)
    (S : Fin (Spin.Structured.bsched m) →
      Finset (Fin (Spin.Structured.Lsched m))) :
    matrixExpect22 (Spin.Structured.ConcreteRoute.regionKernel S)
      (fun regions => Matrix.of (inputListKernel z
        (List.ofFn (scalableFlatInput m regions)))) =
    (List.ofFn (fun j : Fin (Spin.Structured.bsched m) =>
      Matrix.of (averagedRegionKernel
        (scalableRegionSupport m (S j)) z))).prod := by
  simpa only [scalableAveragedRegion_eq] using
    scalable_region_shuffle_product m z S

end Spin.UniformMixer
