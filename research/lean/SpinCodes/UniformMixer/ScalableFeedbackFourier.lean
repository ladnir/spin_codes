import SpinCodes.UniformMixer.ScalableMomentSpectrumExact
import SpinCodes.Structured.ConcreteFourierProduct

/-!
# Fourier characters for the actual feedback map

The zero-feedback event can be averaged over all entering-state
characters. The adjoint pairing identifies each such character with a
production expansion word, whose exact spectrum is already checked.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.ConcreteFourier

def statePairing (q r : State22) : ZMod 2 :=
  ∑ i : Fin 22, q i * r i

def statePairingLinear (r : State22) : State22 →ₗ[ZMod 2] ZMod 2 where
  toFun q := statePairing q r
  map_add' q t := by
    simp [statePairing, Pi.add_apply, add_mul, Finset.sum_add_distrib]
  map_smul' a q := by
    simp only [statePairing, Pi.smul_apply, smul_eq_mul,
      Finset.mul_sum, RingHom.id_apply]
    apply Finset.sum_congr rfl
    intro i _
    ring

theorem statePairingLinear_basis (r : State22) (i : Fin 22) :
    statePairingLinear r (basisVector i) = r i := by
  simp [statePairingLinear, statePairing, basisVector,
    Finset.sum_ite_eq']

theorem stateCharacter_orthogonality (r : State22) :
    (∑ q : State22, binarySign (statePairing q r)) =
      if r = 0 then (2 ^ 22 : ℤ) else 0 := by
  by_cases hr : r = 0
  · subst r
    simp [statePairing, binarySign_zero]
  · have hex : ∃ i : Fin 22, r i ≠ 0 := by
      by_contra h
      simp only [not_exists, not_not] at h
      apply hr
      funext i
      exact h i
    obtain ⟨i, hi⟩ := hex
    have hr1 : r i = 1 :=
      (zmod2_zero_or_one (r i)).resolve_left hi
    have hlin : statePairingLinear r (basisVector i) = 1 := by
      rw [statePairingLinear_basis, hr1]
    have hsum := binarySign_sum_linear (statePairingLinear r)
      (basisVector i) hlin
    change (∑ q : State22, binarySign (statePairing q r)) = 0 at hsum
    simpa only [hr, ite_false] using hsum

theorem statePairing_feedback (q : State22) (x : Input64) :
    statePairing q (stepFeedback x) =
      ∑ j : Fin 64, x j * coordMap j q := by
  unfold statePairing stepFeedback
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro j _
  change (∑ i : Fin 22, q i * (x j * coordCoeff j i)) =
    x j * (∑ i : Fin 22, coordCoeff j i * q i)
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem binarySign_finset_sum {α : Type*} [DecidableEq α]
    (s : Finset α) (f : α → ZMod 2) :
    binarySign (∑ i ∈ s, f i) =
      ∏ i ∈ s, binarySign (f i) := by
  induction s using Finset.induction_on with
  | empty => simp [binarySign_zero]
  | @insert i s hi ih =>
      simp [hi, binarySign_add, ih]

theorem feedback_indicator_character (x : Input64) :
    (2 ^ 22 : ℝ) * (if stepFeedback x = 0 then 1 else 0) =
      ∑ q : State22,
        (binarySign (statePairing q (stepFeedback x)) : ℝ) := by
  have h := stateCharacter_orthogonality (stepFeedback x)
  have hc := congrArg (fun n : ℤ => (n : ℝ)) h
  push_cast at hc
  by_cases hx : stepFeedback x = 0
  · simpa only [hx, ite_true, mul_one,
      show (2 ^ 22 : ℝ) = 4194304 by norm_num] using hc.symm
  · simpa only [hx, ite_false, mul_zero] using hc.symm

theorem feedback_character_support (q : State22)
    (S : Finset (Fin 64)) :
    binarySign (statePairing q (stepFeedback (supportInput S))) =
      ∏ j ∈ S, binarySign (coordMap j q) := by
  rw [statePairing_feedback, binarySign_finset_sum]
  have hfactor (j : Fin 64) :
      binarySign (supportInput S j * coordMap j q) =
        if j ∈ S then binarySign (coordMap j q) else 1 := by
    by_cases hj : j ∈ S <;>
      simp [supportInput, hj, binarySign_zero]
  simp_rw [hfactor]
  rw [Finset.prod_ite]
  have hyes : (Finset.univ : Finset (Fin 64)).filter (fun j => j ∈ S) = S := by
    ext j
    simp
  rw [hyes]
  simp

def feedbackWeightMoment (β z : ℝ) (w : Nat) : ℝ :=
  (1 - β + β * z) ^ (64 - w) * (1 - β - β * z) ^ w

theorem feedbackCharacterWeighted_moment (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q : State22) :
    (bernoulliInputLaw β hβ0 hβ1).expect (fun x =>
      z ^ stepWeight x 0 *
        (binarySign (statePairing q (stepFeedback x)) : ℝ)) =
      feedbackWeightMoment β z (coordinateWeight coordMap q) := by
  rw [bernoulliInputLaw, FinPMF.expect_map]
  have hfactor (S : Finset (Fin 64)) :
      z ^ stepWeight (supportInput S) 0 *
        (binarySign (statePairing q (stepFeedback (supportInput S))) : ℝ) =
      ∏ j ∈ S, z * (binarySign (coordMap j q) : ℝ) := by
    rw [stepWeight_supportInput_zero, feedback_character_support]
    push_cast
    rw [Finset.prod_mul_distrib, Finset.prod_const]
  simp_rw [hfactor]
  rw [poissonBinom_product_moment (fun _ : Fin 64 => β)
    (fun _ => hβ0) (fun _ => hβ1)
    (fun j => z * (binarySign (coordMap j q) : ℝ))]
  have hform :
      (∏ j : Fin 64,
        (1 - β + β * (z * (binarySign (coordMap j q) : ℝ)))) =
      ∏ j : Fin 64,
        if coordMap j q = 1 then 1 - β - β * z
        else 1 - β + β * z := by
    apply Finset.prod_congr rfl
    intro j _
    by_cases hj : coordMap j q = 1
    · simp only [hj, ite_true, binarySign_one, Int.cast_neg,
        Int.cast_one]
      ring
    · have hz : coordMap j q = 0 :=
        (zmod2_zero_or_one (coordMap j q)).resolve_right hj
      simp [hz, binarySign_zero]
  rw [hform, Finset.prod_ite]
  have hyes : (Finset.univ : Finset (Fin 64)).filter
      (fun j => coordMap j q = 1) = expansionSupport q := rfl
  have hno : (Finset.univ : Finset (Fin 64)).filter
      (fun j => ¬coordMap j q = 1) =
        Finset.univ \ expansionSupport q := by
    ext j
    simp [expansionSupport]
  rw [hyes, hno]
  simp only [Finset.prod_const,
    Finset.card_sdiff_of_subset (Finset.subset_univ (expansionSupport q)),
    Finset.card_univ, Fintype.card_fin, expansionSupport_card]
  unfold feedbackWeightMoment
  ring

end Spin.UniformMixer
