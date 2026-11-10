import SpinCodes.Turbo.PhysicalInputPolynomial22
import SpinCodes.Turbo.PhysicalFeedbackFourier22
import SpinCodes.Turbo.PhysicalProfile22

/-!
# Occupation polynomials for actual feedback characters

The adjoint identity from the old Scalable proof identifies a feedback
character with the expansion word of its character state. Exact packet
label sums, rather than independent-bit averages, then give the factors.
-/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial Spin.UniformMixer

theorem physicalStepCharacter22_product (q : State22) (x : PacketStep16) :
    (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ) =
      ∏ slot : Fin 16, packetCharacter22 (expandedPacket22 q slot) (x slot) := by
  rw [statePairing_feedback,
    ← Equiv.sum_comp (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64),
    Fintype.sum_prod_type, binarySign_finset_sum]
  push_cast
  apply prod_congr rfl
  intro slot _
  unfold packetCharacter22 expandedPacket22 packetStepInput22
  simp only [Equiv.symm_apply_apply]

def packetCharacterLabelPolynomial22 (v p : Packet4) (z : ℝ) : Polynomial ℝ :=
  Polynomial.monomial (if p = 0 then 0 else 1)
    (z ^ packetWeight22 p * packetCharacter22 v p)

theorem packetCharacterLabelPolynomial22_sum (v : Packet4) (z : ℝ) :
    (∑ p : Packet4, packetCharacterLabelPolynomial22 v p z) =
      Polynomial.monomial 0 1 + Polynomial.monomial 1
        ((1+z)^(4-packetWeight22 v) * (1-z)^packetWeight22 v - 1) := by
  classical
  have hs := Fintype.sum_subtype_add_sum_subtype (fun p : Packet4 => p ≠ 0)
    (fun p => packetCharacterLabelPolynomial22 v p z)
  have hn : (∑ p : {p : Packet4 // p ≠ 0},
      packetCharacterLabelPolynomial22 v p.val z) =
      Polynomial.monomial 1
        ((1+z)^(4-packetWeight22 v) * (1-z)^packetWeight22 v - 1) := by
    have he (p : {p : Packet4 // p ≠ 0}) :
        packetCharacterLabelPolynomial22 v p.val z =
          Polynomial.monomial 1 (z ^ packetWeight22 p.val * packetCharacter22 v p.val) := by
      unfold packetCharacterLabelPolynomial22
      rw [if_neg p.property]
    simp only [he]
    rw [← map_sum, packetCharacter22_nonzero_moment]
  have hz : (∑ p : {p : Packet4 // ¬ p ≠ 0},
      packetCharacterLabelPolynomial22 v p.val z) = Polynomial.monomial 0 1 := by
    have he (p : {p : Packet4 // ¬ p ≠ 0}) : p.val = 0 := not_not.mp p.property
    simp only [he, packetCharacterLabelPolynomial22, ite_true, packetWeight22_zero,
      pow_zero, packetCharacter22_zero, mul_one, sum_const, card_univ]
    have hc : Fintype.card {p : Packet4 // ¬ p ≠ 0} = 1 := by
      simp only [not_not]
      exact Fintype.card_subtype_eq 0
    rw [hc, one_nsmul]
  rw [hn, hz] at hs
  exact hs.symm.trans (add_comm _ _)

def physicalCharacterOccupationPolynomial22 (q : State22) (z : ℝ) : Polynomial ℝ :=
  degreePolynomial22 (Ω := PacketStep16) (fun x => (packetStepSupport22 x).card)
    (fun x : PacketStep16 => z ^ stepWeight (packetStepInput22 x) 0 *
      (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ))

theorem physicalCharacterOccupationPolynomial22_product (q : State22) (z : ℝ) :
    physicalCharacterOccupationPolynomial22 q z =
      ∏ slot : Fin 16, (Polynomial.monomial 0 1 + Polynomial.monomial 1
        ((1+z)^(4-packetWeight22 (expandedPacket22 q slot)) *
          (1-z)^packetWeight22 (expandedPacket22 q slot) - 1)) := by
  unfold physicalCharacterOccupationPolynomial22 degreePolynomial22
  have he (x : PacketStep16) :
      Polynomial.monomial (packetStepSupport22 x).card
        (z ^ stepWeight (packetStepInput22 x) 0 *
          (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) =
      ∏ slot : Fin 16, packetCharacterLabelPolynomial22 (expandedPacket22 q slot) (x slot) z := by
    unfold packetCharacterLabelPolynomial22
    rw [packetMonomial_product22]
    apply congrArg₂ (fun n c => Polynomial.monomial n c)
    · simp only [packetStepSupport22, card_filter]
      apply sum_congr rfl
      intro slot _
      split_ifs <;> simp_all
    · rw [prod_mul_distrib, prod_pow_eq_pow_sum, physicalStepCharacter22_product]
      apply congrArg (fun e : ℕ => z^e * ∏ slot : Fin 16,
        packetCharacter22 (expandedPacket22 q slot) (x slot))
      rw [stepWeight_eq_packet_weights22]
      apply sum_congr rfl
      intro slot _
      have hz : expandedPacket22 0 slot = 0 := by
        funext bit
        simp [expandedPacket22]
      rw [hz, zero_add]
  simp_rw [he]
  rw [← Fintype.prod_sum (fun (slot : Fin 16) (p : Packet4) =>
    packetCharacterLabelPolynomial22 (expandedPacket22 q slot) p z)]
  simp only [packetCharacterLabelPolynomial22_sum]

theorem physicalCharacterOccupationPolynomial22_eq_canonical (q : State22) (z : ℝ) :
    physicalCharacterOccupationPolynomial22 q z =
      degreePolynomial22 (Ω := PacketStep16) (fun x => (packetStepSupport22 x).card)
        (fun x : PacketStep16 => z ^ stepWeight (packetStepInput22 x) 0 *
          (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) := by
  rfl

theorem physicalCharacterMoment22_eq_coefficient (j : Fin 17) (q : State22) (z : ℝ) :
    (packetStepLaw22 j).expect (fun x => z ^ stepWeight (packetStepInput22 x) 0 *
      (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) =
      (physicalCharacterOccupationPolynomial22 q z).coeff j.val /
        ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val) := by
  conv_rhs => rw [physicalCharacterOccupationPolynomial22_eq_canonical]
  exact packetStepLaw22_expect_coefficient j
    (fun x : PacketStep16 => z ^ stepWeight (packetStepInput22 x) 0 *
      (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ))

#print axioms physicalStepCharacter22_product
#print axioms physicalCharacterOccupationPolynomial22_product
#print axioms physicalCharacterMoment22_eq_coefficient
end Spin.Turbo
