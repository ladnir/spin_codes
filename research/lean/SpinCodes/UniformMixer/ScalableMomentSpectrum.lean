import SpinCodes.UniformMixer.ScalableEnvelopeRows
import SpinCodes.Structured.ConcreteScalarBits

/-!
# The Bernoulli row moment as an expansion-weight polynomial

The iid input bits XOR the fixed expansion of an entering state. This
identifies the weighted row sum with a polynomial in the Hamming weight
of that expansion. The numerical certificate can therefore use the
actual finite image spectrum of `coordMap`.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.ConcreteScalar
open scoped symmDiff

def expansionSupport (q : State22) : Finset (Fin 64) :=
  Finset.univ.filter (fun j => coordMap j q = 1)

theorem expansionSupport_card (q : State22) :
    (expansionSupport q).card = coordinateWeight coordMap q := by
  rw [expansionSupport, Finset.card_filter]
  rfl

theorem stepWeight_supportInput (S : Finset (Fin 64)) (q : State22) :
    stepWeight (supportInput S) q = (S ∆ expansionSupport q).card := by
  classical
  have hset : S ∆ expansionSupport q =
      Finset.univ.filter
        (fun j => coordMap j q + supportInput S j = 1) := by
    ext j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    by_cases hj : j ∈ S <;> by_cases hq : coordMap j q = 1
    · simp [Finset.mem_symmDiff, expansionSupport, supportInput, hj, hq]
    · have hzero : coordMap j q = 0 :=
        (zmod2_zero_or_one (coordMap j q)).resolve_right hq
      simp [Finset.mem_symmDiff, expansionSupport, supportInput, hj, hzero]
    · simp [Finset.mem_symmDiff, expansionSupport, supportInput, hj, hq]
    · have hzero : coordMap j q = 0 :=
        (zmod2_zero_or_one (coordMap j q)).resolve_right hq
      simp [Finset.mem_symmDiff, expansionSupport, supportInput, hj, hzero]
  rw [hset, Finset.card_filter]
  rfl

theorem scalableLiveMoment_weight_formula (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q : State22) :
    scalableLiveMoment β z hβ0 hβ1 q =
      (1 - β + β * z) ^ (64 - coordinateWeight coordMap q) *
        (β + (1 - β) * z) ^ coordinateWeight coordMap q := by
  unfold scalableLiveMoment
  rw [bernoulliInputLaw, FinPMF.expect_map]
  simp_rw [stepWeight_supportInput]
  rw [bernoulli_xor_moment β z hβ0 hβ1 (expansionSupport q)]
  rw [expansionSupport_card]

theorem stepWeight_supportInput_zero (S : Finset (Fin 64)) :
    stepWeight (supportInput S) 0 = S.card := by
  rw [stepWeight_supportInput]
  have h : expansionSupport (0 : State22) = ∅ := by
    ext j
    simp [expansionSupport]
  simp [h, Finset.symmDiff_def]

/-- The zero-to-zero row is the weight polynomial of the actual
feedback kernel. This identity does not assume that feedback is uniform. -/
theorem scalableZeroHold_kernel_sum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableZeroHold β z hβ0 hβ1 =
      ∑ S : Finset (Fin 64),
        if stepFeedback (supportInput S) = 0 then
          (β * z) ^ S.card * (1 - β) ^ (64 - S.card)
        else 0 := by
  unfold scalableZeroHold
  rw [bernoulliInputLaw, FinPMF.expect_map]
  simp only [FinPMF.expect, poissonBinom_const_apply]
  apply Finset.sum_congr rfl
  intro S _
  rw [stepWeight_supportInput_zero]
  split_ifs <;> ring

end Spin.UniformMixer
