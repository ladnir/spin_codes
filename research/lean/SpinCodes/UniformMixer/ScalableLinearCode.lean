import SpinCodes.UniformMixer.ScalableCodeProperties
import SpinCodes.UniformMixer.ScalableInputEmbedding
import SpinCodes.Structured.ConcreteOuterNativeLinearity

/-!
# The realized Scalable binary linear code

This follows the old `ConcreteNativeLinearCode` assembly, replacing its
IMT output by the actual sampled 64-bit-step GL inner. Output coordinates
remain grouped into 64-bit words; the coordinate count is exactly `Nsched m`.
-/

noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding Spin.Structured.ConcreteRoute
open scoped symmDiff

abbrev ScalableOutput (m : ℕ) :=
  Fin (bsched m * scalableStepsPerRegion m) → Input64

def scalableCodeword (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) : ScalableOutput m :=
  sampledOutput (scalableFlatInput m
    (routeEval (nativeRows m seed.1 ((wordEquiv (Nsched m / 2)).symm x)) seed.2.1))
    (scalableFlatMaps m seed.2.2) 0

theorem scalableCodeword_weight (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) :
    sampledWordWeight (scalableCodeword m seed x) =
      scalableOutputWeight m ((wordEquiv (Nsched m / 2)).symm x) seed :=
  sampledOutput_weight _ _ _

theorem scalableCodeword_injective (m : ℕ) (seed : ScalableSetupSeed m) :
    Function.Injective (scalableCodeword m seed) :=
  (sampledOutput_injective (scalableFlatMaps m seed.2.2) 0).comp
    ((scalableRoutedInput_injective m seed.1 seed.2.1).comp
      (wordEquiv (Nsched m / 2)).symm.injective)

theorem scalableCodeword_add (m : ℕ) (seed : ScalableSetupSeed m)
    (x y : Fin (Nsched m / 2) → ZMod 2) :
    scalableCodeword m seed (x + y) =
      scalableCodeword m seed x + scalableCodeword m seed y := by
  unfold scalableCodeword
  rw [wordEquiv_symm_add, nativeRows_xor, routeEval_xor, scalableFlatInput_xor]
  simpa only [add_zero] using sampledOutput_add _ _ (scalableFlatMaps m seed.2.2) 0 0

theorem scalableCodeword_zero (m : ℕ) (seed : ScalableSetupSeed m) :
    scalableCodeword m seed 0 = 0 := by
  have h := scalableCodeword_add m seed 0 0
  simp only [add_zero] at h
  exact (add_left_cancel (h.symm.trans (add_zero _).symm) : _)

def scalableEncoderLinearMap (m : ℕ) (seed : ScalableSetupSeed m) :
    (Fin (Nsched m / 2) → ZMod 2) →ₗ[ZMod 2] ScalableOutput m where
  toFun := scalableCodeword m seed
  map_add' := scalableCodeword_add m seed
  map_smul' c x := by
    have hc : c = 0 ∨ c = 1 := (by decide : ∀ c : ZMod 2, c = 0 ∨ c = 1) c
    rcases hc with rfl | rfl <;> simp [scalableCodeword_zero]

def scalableRealizedCode (m : ℕ) (seed : ScalableSetupSeed m) :
    Submodule (ZMod 2) (ScalableOutput m) :=
  (scalableEncoderLinearMap m seed).range

theorem scalableRealizedCode_finrank (m : ℕ) (seed : ScalableSetupSeed m) :
    Module.finrank (ZMod 2) (scalableRealizedCode m seed) = Nsched m / 2 := by
  rw [scalableRealizedCode, LinearMap.finrank_range_of_inj (scalableCodeword_injective m seed)]
  simp

theorem scalable_output_length (m : ℕ) :
    (bsched m * scalableStepsPerRegion m) * 64 = Nsched m := by
  rw [Nsched_eq, scalable_region_length m]
  ring

theorem scalableOutput_finrank (m : ℕ) :
    Module.finrank (ZMod 2) (ScalableOutput m) = Nsched m := by
  simpa only [ScalableOutput, Input64, Module.finrank_pi_fintype, Module.finrank_pi, Module.finrank_self,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, smul_eq_mul, mul_one] using
    scalable_output_length m

theorem scalableRealizedCode_rate (m : ℕ) (seed : ScalableSetupSeed m) :
    (Module.finrank (ZMod 2) (scalableRealizedCode m seed) : ℝ) /
      (Nsched m : ℝ) = 1 / 2 := by
  rw [scalableRealizedCode_finrank]
  have hN : (0 : ℝ) < Nsched m := by exact_mod_cast Nsched_pos m
  have he : Nsched m / 2 * 2 = Nsched m := by
    apply Nat.div_mul_cancel
    rw [Nsched_eq]
    apply dvd_mul_of_dvd_left
    unfold Lsched
    omega
  have he' : ((Nsched m / 2 : ℕ) : ℝ) * 2 = (Nsched m : ℝ) := by exact_mod_cast he
  apply (div_eq_iff hN.ne').mpr
  linarith

#print axioms scalableCodeword_weight
#print axioms scalableCodeword_injective
#print axioms scalableEncoderLinearMap
#print axioms scalableRealizedCode_rate

end Spin.UniformMixer
