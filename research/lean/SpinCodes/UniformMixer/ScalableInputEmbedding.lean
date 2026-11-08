import SpinCodes.UniformMixer.ScalableRealization
import SpinCodes.Structured.ConcreteRoutedOutput

/-!
# Injectivity of the Scalable outer and input wiring

The old `ConcreteRoutedOutput.routeEval_injective` applies unchanged:
the revised encoder uses the same one-bit regional permutation. Only the
64-bit chronological reshape is new. These lemmas concern that reshape
and the shared outer, without invoking the old IMT inner.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
open scoped symmDiff

theorem regionInput_injective {R : ℕ} :
    Function.Injective (regionInput (R := R)) := by
  intro S T h
  ext x
  have hx := congrFun (congrFun h (regionStep x)) (regionLocal x)
  simp only [regionInput, regionBit_reconstruct] at hx
  by_cases hs : x ∈ S <;> by_cases ht : x ∈ T <;> simp_all

theorem scalableRegionInput_injective (m : ℕ) :
    Function.Injective (scalableRegionInput m) := by
  intro S T h
  have hh := regionInput_injective h
  exact Finset.map_injective (finCongr (scalable_region_length m)).toEmbedding hh

theorem scalableFlatInput_injective (m : ℕ) :
    Function.Injective (scalableFlatInput m) := by
  intro S T h
  funext j
  apply scalableRegionInput_injective m
  funext i
  have hh := congrFun h (finProdFinEquiv (j, i))
  simpa only [scalableFlatInput, Equiv.symm_apply_apply] using hh

theorem scalableRoutedInput_injective (m : ℕ)
    (out : NativeSeed m) (route : ScalableRouteSeed m) :
    Function.Injective (fun message =>
      scalableFlatInput m
        (Spin.Structured.ConcreteRoute.routeEval
          (nativeRows m out message) route)) :=
  (scalableFlatInput_injective m).comp
    ((Spin.Structured.ConcreteRoute.routeEval_injective route).comp
      (native_rows_injective m out))

theorem regionInput_xor {R : ℕ} (S T : Finset (Fin (64 * R))) :
    regionInput (S ∆ T) = regionInput S + regionInput T := by
  funext i j
  by_cases hs : regionBit i j ∈ S <;>
    by_cases ht : regionBit i j ∈ T <;>
    norm_num [regionInput, Finset.mem_symmDiff, hs, ht] <;> decide

theorem scalableRegionSupport_xor (m : ℕ)
    (S T : Finset (Fin (Lsched m))) :
    scalableRegionSupport m (S ∆ T) =
      scalableRegionSupport m S ∆ scalableRegionSupport m T := by
  ext i
  simp [scalableRegionSupport, Finset.mem_symmDiff]

theorem scalableRegionInput_xor (m : ℕ)
    (S T : Finset (Fin (Lsched m))) :
    scalableRegionInput m (S ∆ T) =
      scalableRegionInput m S + scalableRegionInput m T := by
  unfold scalableRegionInput
  rw [scalableRegionSupport_xor, regionInput_xor]

theorem scalableFlatInput_xor (m : ℕ)
    (S T : Fin (bsched m) → Finset (Fin (Lsched m))) :
    scalableFlatInput m (fun j => S j ∆ T j) =
      scalableFlatInput m S + scalableFlatInput m T := by
  funext i
  exact congrFun (scalableRegionInput_xor m
    (S (finProdFinEquiv.symm i).1) (T (finProdFinEquiv.symm i).1))
    (finProdFinEquiv.symm i).2

#print axioms scalableRoutedInput_injective

end Spin.UniformMixer
