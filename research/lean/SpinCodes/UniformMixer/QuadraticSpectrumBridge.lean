import SpinCodes.UniformMixer.QuadraticWalsh
import SpinCodes.UniformMixer.Map22PackedBridge

/-!
# Connect six-variable Walsh weights to the production expansion

The quadratic normal-form proof uses output indices `Fin 64`; Walsh
analysis uses six-bit vectors. This module identifies the two index sets
by their actual binary coordinates and transfers the production weight.
-/

noncomputable section

namespace Spin.UniformMixer

def outputIndexState (j : Fin 64) : State6 :=
  fun i => bitScalar (j.val.testBit i.val)

theorem outputIndexState_injective : Function.Injective outputIndexState := by
  intro j k h
  apply Fin.ext
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < 6
  · exact bitScalar_injective (congrFun h ⟨i, hi⟩)
  · have hpow : 2 ^ 6 ≤ 2 ^ i := Nat.pow_le_pow_right (by decide) (by omega)
    have hj0 : j.val < 2 ^ 6 := j.isLt
    have hk0 : k.val < 2 ^ 6 := k.isLt
    have hj : j.val < 2 ^ i := hj0.trans_le hpow
    have hk : k.val < 2 ^ i := hk0.trans_le hpow
    rw [Nat.testBit_lt_two_pow hj, Nat.testBit_lt_two_pow hk]

theorem outputIndexState_bijective : Function.Bijective outputIndexState := by
  have hcard : Fintype.card (Fin 64) = Fintype.card State6 := by
    simp [State6, ZMod.card]
  exact (Fintype.bijective_iff_injective_and_card _).mpr
    ⟨outputIndexState_injective, hcard⟩

def outputIndexEquiv : Fin 64 ≃ State6 :=
  Equiv.ofBijective outputIndexState outputIndexState_bijective

theorem productionWeight_as_quadraticWeight (q : State22) :
    coordinateWeight coordMap q =
      boolWeight6 (fun x => quadraticEvaluation
        (productionQuadraticCoefficients q) (outputIndexEquiv.symm x)) := by
  unfold coordinateWeight boolWeight6
  apply Fintype.sum_equiv outputIndexEquiv
  intro j
  change (if coordMap j q = 1 then 1 else 0) =
    if quadraticEvaluation (productionQuadraticCoefficients q)
      (outputIndexEquiv.symm (outputIndexEquiv j)) = 1 then 1 else 0
  rw [outputIndexEquiv.symm_apply_apply, production_quadratic_evaluation]

end Spin.UniformMixer
