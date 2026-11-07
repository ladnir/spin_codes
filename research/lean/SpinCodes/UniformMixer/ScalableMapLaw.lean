import SpinCodes.UniformMixer.ScalableRealization

/-!
# Flattening the independently sampled state maps

The paper samples a map for every region and step. Serializing the
regions gives the same independent uniform law on the full step stream.
-/

noncomputable section

namespace Spin.UniformMixer

attribute [local instance] Classical.propDecidable

theorem pi_uniform_eq_uniform {ι α : Type*}
    [Fintype ι] [DecidableEq ι] [Fintype α] [Nonempty α] :
    (piPMF (fun _ : ι => FinPMF.uniform α)) =
      FinPMF.uniform (ι → α) := by
  apply FinPMF.ext
  funext x
  simp [piPMF_apply, FinPMF.uniform]

theorem scalableStateMapLaw_eq_uniform (m : ℕ) :
    scalableStateMapLaw m = FinPMF.uniform (ScalableStateMaps m) := by
  simp only [scalableStateMapLaw, pi_uniform_eq_uniform]

def scalableFlatMapsEquiv (m : ℕ) :
    ScalableStateMaps m ≃
      (Fin (Spin.Structured.bsched m * scalableStepsPerRegion m) → StateMap) where
  toFun := scalableFlatMaps m
  invFun := fun flat j k => flat (finProdFinEquiv (j, k))
  left_inv := by
    intro maps
    funext j k
    simp [scalableFlatMaps]
  right_inv := by
    intro flat
    funext i
    change flat (finProdFinEquiv (finProdFinEquiv.symm i)) = flat i
    exact congrArg flat (Equiv.apply_symm_apply finProdFinEquiv i)

theorem scalableFlatMaps_law (m : ℕ) :
    (scalableStateMapLaw m).map (scalableFlatMaps m) =
      piPMF (fun _ : Fin (Spin.Structured.bsched m *
        scalableStepsPerRegion m) => FinPMF.uniform StateMap) := by
  rw [scalableStateMapLaw_eq_uniform, pi_uniform_eq_uniform]
  apply FinPMF.ext
  funext flat
  rw [show scalableFlatMaps m = (scalableFlatMapsEquiv m :
    ScalableStateMaps m ≃
      (Fin (Spin.Structured.bsched m * scalableStepsPerRegion m) → StateMap))
      from rfl]
  rw [FinPMF.map_equiv_apply]
  simp only [FinPMF.uniform]
  rw [Fintype.card_congr (scalableFlatMapsEquiv m)]

theorem scalableFlatMaps_expect (m : ℕ)
    (f : (Fin (Spin.Structured.bsched m *
      scalableStepsPerRegion m) → StateMap) → ℝ) :
    (scalableStateMapLaw m).expect (fun maps => f (scalableFlatMaps m maps)) =
      (piPMF (fun _ : Fin (Spin.Structured.bsched m *
        scalableStepsPerRegion m) => FinPMF.uniform StateMap)).expect f := by
  rw [← FinPMF.expect_map, scalableFlatMaps_law]

end Spin.UniformMixer
