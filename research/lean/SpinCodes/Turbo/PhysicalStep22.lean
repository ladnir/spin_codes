import SpinCodes.Turbo.RemainingRoute22
import SpinCodes.UniformMixer.ScalableStepLaw

/-!
# Exact packet-valued physical step for the 22-state Turbo inner

A physical 64-bit step contains sixteen four-bit packet slots. At exact
occupancy `j`, its input is uniform among complete packet arrays with
exactly `j` nonzero slots. Equivalently the active support is a uniform
`j`-subset and its labels are independent uniform nonzero four-bit
values. The state update is an independent uniform GL(22,2) map.

This exact full-state kernel is the semantic target for the retained
11-coordinate outward physical operators. No quotient or numerical
bound is asserted in this file.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

abbrev PacketStep16 := Fin 16 → Packet4

def packetStepSupport22 (x : PacketStep16) : Finset (Fin 16) :=
  Finset.univ.filter (fun i => x i ≠ 0)

abbrev PacketStepLayer22 (j : Fin 17) :=
  {x : PacketStep16 // (packetStepSupport22 x).card = j.val}

noncomputable instance (j : Fin 17) : Fintype (PacketStepLayer22 j) :=
  Fintype.ofFinite _

noncomputable instance (j : Fin 17) : Nonempty (PacketStepLayer22 j) := by
  classical
  have hj : j.val ≤ (Finset.univ : Finset (Fin 16)).card := by
    simp only [Finset.card_univ, Fintype.card_fin]
    omega
  obtain ⟨S, -, hcard⟩ :=
    (Finset.univ : Finset (Fin 16)).exists_subset_card_eq hj
  let x : PacketStep16 := fun i => if i ∈ S then 1 else 0
  have hsupp : packetStepSupport22 x = S := by
    ext i
    simp [packetStepSupport22, x]
  exact ⟨⟨x, by rw [hsupp, hcard]⟩⟩

/-- Uniform complete packet input at a specified exact occupancy. -/
def packetStepLaw22 (j : Fin 17) : FinPMF PacketStep16 :=
  (FinPMF.uniform (PacketStepLayer22 j)).map Subtype.val

/-- A fixed support has one independent nonzero four-bit label per slot. -/
def packetStepSupportEquiv22 (S : Finset (Fin 16)) :
    {x : PacketStep16 // packetStepSupport22 x = S} ≃
      (S → {p : Packet4 // p ≠ 0}) := by
  classical
  refine ⟨(fun x i => ⟨x.val i.val, ?_⟩),
    (fun f => ⟨(fun i => if h : i ∈ S then (f ⟨i, h⟩).val else 0), ?_⟩),
    ?_, ?_⟩
  · have hi : i.val ∈ packetStepSupport22 x.val := by
      rw [x.property]
      exact i.property
    exact (Finset.mem_filter.mp hi).2
  · ext i
    simp only [packetStepSupport22, Finset.mem_filter,
      Finset.mem_univ, true_and]
    by_cases hi : i ∈ S
    · simp [hi, (f ⟨i, hi⟩).property]
    · simp [hi]
  · intro x
    apply Subtype.ext
    funext i
    by_cases hi : i ∈ S
    · simp [hi]
    · have hz : x.val i = 0 := by
        by_contra hn
        have hm : i ∈ packetStepSupport22 x.val :=
          Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
        rw [x.property] at hm
        exact hi hm
      simp [hi, hz]
  · intro f
    funext i
    apply Subtype.ext
    simp

theorem packetStepSupport_card22 (S : Finset (Fin 16)) :
    Fintype.card {x : PacketStep16 // packetStepSupport22 x = S} =
      15 ^ S.card := by
  rw [Fintype.card_congr (packetStepSupportEquiv22 S)]
  simp

/-- The complete exact-occupancy layer has `choose(16,j)*15^j`
elements. This is the count assumed by the physical operator producer. -/
theorem packetStepLayer_card22 (j : Fin 17) :
    Fintype.card (PacketStepLayer22 j) =
      Nat.choose 16 j.val * 15 ^ j.val := by
  classical
  let supports := (Finset.univ : Finset (Fin 16)).powersetCard j.val
  have hsum := Finset.sum_card_fiberwise_eq_card_filter
    (Finset.univ : Finset PacketStep16) supports packetStepSupport22
  have hterm (S : Finset (Fin 16)) :
      (Finset.univ.filter (fun x : PacketStep16 =>
        packetStepSupport22 x = S)).card = 15 ^ S.card := by
    rw [← Fintype.card_subtype]
    exact packetStepSupport_card22 S
  calc
    Fintype.card (PacketStepLayer22 j) =
        (Finset.univ.filter (fun x : PacketStep16 =>
          packetStepSupport22 x ∈ supports)).card := by
      rw [Fintype.card_subtype]
      apply congrArg Finset.card
      apply Finset.filter_congr
      intro x _
      simp [supports, Finset.mem_powersetCard]
    _ = ∑ S ∈ supports,
          (Finset.univ.filter (fun x : PacketStep16 =>
            packetStepSupport22 x = S)).card := hsum.symm
    _ = ∑ _S ∈ supports, 15 ^ j.val := by
      apply Finset.sum_congr rfl
      intro S hS
      rw [hterm, (Finset.mem_powersetCard.mp hS).2]
    _ = Nat.choose 16 j.val * 15 ^ j.val := by
      simp [supports, Finset.card_powersetCard]

/-- Every complete word in the exact-occupancy layer has the same
probability. The denominator counts both supports and nonzero labels. -/
theorem packetStepLaw22_apply (j : Fin 17) (x : PacketStep16) :
    (packetStepLaw22 j).p x =
      if (packetStepSupport22 x).card = j.val then
        1 / ((Nat.choose 16 j.val : ℝ) * (15 : ℝ) ^ j.val)
      else 0 := by
  classical
  rw [packetStepLaw22, FinPMF.map_p]
  by_cases hx : (packetStepSupport22 x).card = j.val
  · rw [if_pos hx]
    let xx : PacketStepLayer22 j := ⟨x, hx⟩
    have h (p : PacketStepLayer22 j) : (p.val = x) ↔ p = xx := by
      constructor
      · intro hp; exact Subtype.ext hp
      · intro hp; exact congrArg Subtype.val hp
    simp only [h, Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card (PacketStepLayer22 j) : ℝ) = _
    rw [packetStepLayer_card22]
    push_cast
    rfl
  · rw [if_neg hx]
    apply Finset.sum_eq_zero
    intro p _
    have hp : p.val ≠ x := by
      intro he
      apply hx
      rw [← he]
      exact p.property
    simp [hp]

def packetStepInput22 (x : PacketStep16) : Input64 :=
  fun bit =>
    let slot : Fin 16 × Fin 4 :=
      (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64).symm bit
    x slot.1 slot.2

/-- The physical-step coordinates in this local law are exactly the
sixteen packet slots used by the actual region-major Turbo encoder. -/
theorem stepInput_eq_packetStepInput22 {R : ℕ}
    (packets : Fin 64 → Fin (R * 16) → Packet4)
    (region : Fin 64) (step : Fin R) :
    stepInput packets region step =
      packetStepInput22 (fun slot =>
        packets region (finProdFinEquiv (step, slot))) := by
  rfl

/-- Averaged weighted transition on the *full* 22-bit state space.
The input distribution is exact, not Bernoulli or iid bitwise. -/
def turboPhysicalKernel22 (j : Fin 17) (z : ℝ) :
    State22 → State22 → ℝ :=
  fun q r => (packetStepLaw22 j).expect (fun x =>
    stepKernel z (packetStepInput22 x) q r)

theorem turboPhysicalKernel22_nonneg (j : Fin 17) (z : ℝ)
    (q r : State22) (hz : 0 ≤ z) :
    0 ≤ turboPhysicalKernel22 j z q r := by
  unfold turboPhysicalKernel22
  apply FinPMF.expect_nonneg
  intro x
  exact mul_nonneg (pow_nonneg hz _)
    (update_nonneg (stepFeedback (packetStepInput22 x)) q r)

theorem turboPhysicalKernel22_sum_one (j : Fin 17)
    (q : State22) :
    (∑ r : State22, turboPhysicalKernel22 j 1 q r) = 1 := by
  simp only [turboPhysicalKernel22]
  rw [← FinPMF.expect_sum]
  have h (x : PacketStep16) :
      (∑ r : State22,
        stepKernel 1 (packetStepInput22 x) q r) = 1 := by
    simp only [stepKernel, one_pow, one_mul]
    exact update_sum (stepFeedback (packetStepInput22 x)) q
  simp only [h, FinPMF.expect_const]

/-- The full-state physical kernel is precisely the expectation over
the packet input and one independently sampled invertible state map.
This reuses the checked Scalable GL transition for each fixed input. -/
theorem turboPhysicalKernel22_sampled
    (j : Fin 17) (z : ℝ) (q r : State22) :
    turboPhysicalKernel22 j z q r =
      (packetStepLaw22 j).expect (fun x =>
        (FinPMF.uniform StateMap).expect (fun g =>
          if g q + stepFeedback (packetStepInput22 x) = r then
            z ^ stepWeight (packetStepInput22 x) q else 0)) := by
  unfold turboPhysicalKernel22
  apply congrArg
  funext x
  let w : ℝ := z ^ stepWeight (packetStepInput22 x) q
  let c : State22 := stepFeedback (packetStepInput22 x)
  have h : (FinPMF.uniform StateMap).expect (fun g =>
      if g q + c = r then w else 0) =
      w * (FinPMF.uniform StateMap).prob (fun g => g q + c = r) := by
    rw [FinPMF.prob_eq_expect_indicator]
    unfold FinPMF.expect
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro g _
    by_cases hg : g q + c = r <;> simp [hg]; ring
  rw [stepKernel, ← uniform_state_map_update c q r]
  exact h.symm

#print axioms packetStepLayer_card22
#print axioms packetStepLaw22_apply
#print axioms stepInput_eq_packetStepInput22
#print axioms turboPhysicalKernel22_nonneg
#print axioms turboPhysicalKernel22_sum_one
#print axioms turboPhysicalKernel22_sampled

end Spin.Turbo
