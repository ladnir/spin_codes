import SpinCodes.Turbo.PhysicalPacketMoment22
import SpinCodes.Turbo.PhysicalLayerPolynomial22
import Mathlib.Algebra.Polynomial.Coeff

/-!
# Polynomial of exact physical-input occupations

The checked exact layer law retains complete packet values. Its generating
polynomial is factored by packet, with one marker for each nonzero packet.
The coefficients therefore describe the actual law before binomial thinning.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer Polynomial

def expandedPacket22 (q : State22) (slot : Fin 16) : Packet4 :=
  fun bit => coordMap (finProdFinEquiv (slot,bit)) q

theorem stepWeight_eq_packet_weights22 (x : PacketStep16) (q : State22) :
    stepWeight (packetStepInput22 x) q =
      ∑ slot : Fin 16, packetWeight22 (expandedPacket22 q slot + x slot) := by
  unfold stepWeight
  rw [← Equiv.sum_comp (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64),
    Fintype.sum_prod_type]
  apply sum_congr rfl
  intro slot _
  unfold packetWeight22
  apply sum_congr rfl
  intro bit _
  simp only [packetStepInput22, Equiv.symm_apply_apply, expandedPacket22, Pi.add_apply]
  split_ifs <;> simp_all

def packetLabelPolynomial22 (v p : Packet4) (z : ℝ) : Polynomial ℝ :=
  Polynomial.monomial (if p = 0 then 0 else 1) (z ^ packetWeight22 (v+p))

theorem packetLabelPolynomial22_sum (v : Packet4) (z : ℝ) :
    (∑ p : Packet4, packetLabelPolynomial22 v p z) =
      Polynomial.monomial 0 (z ^ packetWeight22 v) +
        Polynomial.monomial 1 ((1+z)^4 - z ^ packetWeight22 v) := by
  classical
  have h := Fintype.sum_subtype_add_sum_subtype (fun p : Packet4 => p ≠ 0)
    (fun p => packetLabelPolynomial22 v p z)
  have hn : (∑ p : {p : Packet4 // p ≠ 0}, packetLabelPolynomial22 v p.val z) =
      Polynomial.monomial 1 ((1+z)^4 - z ^ packetWeight22 v) := by
    have he (p : {p : Packet4 // p ≠ 0}) : packetLabelPolynomial22 v p.val z =
        Polynomial.monomial 1 (z ^ packetWeight22 (v+p.val)) := by
      unfold packetLabelPolynomial22
      rw [if_neg p.property]
    simp only [he]
    rw [← map_sum, packetWeight22_nonzero_shift_sum]
  have hz : (∑ p : {p : Packet4 // ¬ p ≠ 0}, packetLabelPolynomial22 v p.val z) =
      Polynomial.monomial 0 (z ^ packetWeight22 v) := by
    have he (p : {p : Packet4 // ¬ p ≠ 0}) : p.val = 0 := not_not.mp p.property
    simp only [he, packetLabelPolynomial22, ite_true, add_zero, Finset.sum_const,
      Finset.card_univ]
    have hc : Fintype.card {p : Packet4 // ¬ p ≠ 0} = 1 := by
      simp only [not_not]
      exact Fintype.card_subtype_eq 0
    rw [hc, one_nsmul]
  rw [hn,hz] at h
  exact h.symm.trans (add_comm _ _)

theorem packetMonomial_product22 {I : Type*} (s : Finset I)
    (d : I → ℕ) (c : I → ℝ) :
    (∏ i ∈ s, Polynomial.monomial (d i) (c i)) =
      Polynomial.monomial (∑ i ∈ s, d i) (∏ i ∈ s, c i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert i s hi ih =>
    simp only [Finset.prod_insert hi, Finset.sum_insert hi, ih,
      Polynomial.monomial_mul_monomial]

def physicalOccupationPolynomial22 (q : State22) (z : ℝ) : Polynomial ℝ :=
  degreePolynomial22 (fun x => (packetStepSupport22 x).card)
    (fun x => z ^ stepWeight (packetStepInput22 x) q)

theorem physicalOccupationPolynomial22_product (q : State22) (z : ℝ) :
    physicalOccupationPolynomial22 q z =
      ∏ slot : Fin 16,
        (Polynomial.monomial 0 (z ^ packetWeight22 (expandedPacket22 q slot)) +
          Polynomial.monomial 1 ((1+z)^4-z ^ packetWeight22 (expandedPacket22 q slot))) := by
  rw [physicalOccupationPolynomial22, degreePolynomial22]
  have he (x : PacketStep16) :
      Polynomial.monomial (packetStepSupport22 x).card
          (z ^ stepWeight (packetStepInput22 x) q) =
        ∏ slot : Fin 16, packetLabelPolynomial22 (expandedPacket22 q slot) (x slot) z := by
    unfold packetLabelPolynomial22
    rw [packetMonomial_product22]
    apply congrArg₂ (fun n c => Polynomial.monomial n c)
    · simp only [packetStepSupport22, Finset.card_filter]
      apply sum_congr rfl
      intro i _
      split_ifs <;> simp_all
    · rw [stepWeight_eq_packet_weights22,
        Finset.prod_pow_eq_pow_sum]
  simp_rw [he]
  rw [← Fintype.prod_sum (fun (slot : Fin 16) (p : Packet4) =>
    packetLabelPolynomial22 (expandedPacket22 q slot) p z)]
  simp only [packetLabelPolynomial22_sum]

def packetStepMoment22 (j : Fin 17) (w : PacketStep16 → ℝ) : ℝ :=
  (packetStepLaw22 j).expect w

/-- The finite polynomial sum is independent of the synthesized enumeration
instance. Transport it explicitly instead of reducing a function-space
enumerator to compare separately hoisted proof arguments. -/
theorem degreePolynomial22_fintype_eq {Ω : Type*} (fi fj : Fintype Ω)
    (degree : Ω → ℕ) (w : Ω → ℝ) :
    @degreePolynomial22 Ω fi degree w = @degreePolynomial22 Ω fj degree w :=
  congrArg (fun ft : Fintype Ω => @degreePolynomial22 Ω ft degree w) (Subsingleton.elim fi fj)

theorem packetStepLaw22_expect_coefficient_fintype (fi : Fintype PacketStep16)
    (j : Fin 17) (w : PacketStep16 → ℝ) :
    packetStepMoment22 j w =
      (@degreePolynomial22 PacketStep16 fi (fun x => (packetStepSupport22 x).card) w).coeff j.val /
        ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val) := by
  conv_lhs => unfold packetStepMoment22
  have hcoeff := packetStepLaw22_expect_coefficient j w
  refine hcoeff.trans ?_
  refine congrArg (fun p : Polynomial ℝ => p.coeff j.val /
    ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val)) ?_
  exact degreePolynomial22_fintype_eq _ fi (fun x => (packetStepSupport22 x).card) w

private theorem physicalInputMoment22_expect_generic (P : FinPMF PacketStep16)
    (q : State22) (z : ℝ) :
    physicalInputMoment22 P z q = P.expect (fun x => z ^ stepWeight (packetStepInput22 x) q) := by
  rfl

/-- Canonicalize the polynomial's enumeration by equality before taking a
coefficient; reducing the two instantiated sums definitionally is expensive. -/
theorem physicalOccupationPolynomial22_eq_canonical (q : State22) (z : ℝ) :
    physicalOccupationPolynomial22 q z =
      degreePolynomial22 (fun x : PacketStep16 => (packetStepSupport22 x).card)
        (fun x : PacketStep16 => z ^ stepWeight (packetStepInput22 x) q) := by
  conv_lhs => unfold physicalOccupationPolynomial22

/-- Coefficients are normalized by the exact complete-word layer size. -/
theorem physicalInputMoment22_eq_coefficient (j : Fin 17) (q : State22) (z : ℝ) :
    physicalInputMoment22 (packetStepLaw22 j) z q =
      (physicalOccupationPolynomial22 q z).coeff j.val /
        ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val) := by
  -- Unrestricted rewriting also probes the polynomial coefficient on the
  -- right, potentially unfolding its 2^64-input enumerator during matching.
  conv_lhs => rw [physicalInputMoment22_expect_generic]
  conv_rhs => rw [physicalOccupationPolynomial22_eq_canonical]
  have hcoeff := packetStepLaw22_expect_coefficient j
    (fun x : PacketStep16 => z ^ stepWeight (packetStepInput22 x) q)
  exact hcoeff

#print axioms physicalOccupationPolynomial22_product
#print axioms packetStepLaw22_expect_coefficient_fintype
#print axioms physicalOccupationPolynomial22_eq_canonical
#print axioms physicalInputMoment22_eq_coefficient
end Spin.Turbo
