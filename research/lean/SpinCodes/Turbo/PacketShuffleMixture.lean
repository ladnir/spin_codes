import SpinCodes.Turbo.SampledRS8Row

/-!
# Disintegrating a packet row by its support

The sampled RS8 row is uniform on each complete-word support fiber.
This file turns that pointwise symmetry into an exact mixture of uniform
active labels, then applies the independently sampled row permutation.
The support mixture is the value-carrying counterpart of the old
`Structured/ConcreteShuffleMixture.lean`.
-/

noncomputable section

set_option maxRecDepth 1000
set_option linter.constructorNameAsVariable false
attribute [local instance] Classical.propDecidable

namespace Spin.Turbo

open Finset

def uniformLabeledLaw (S : Finset (Fin 64)) : FinPMF PacketWord64 :=
  (FinPMF.uniform (LabeledSupport S)).map Subtype.val

theorem uniformLabeledLaw_apply (S : Finset (Fin 64))
    (y : PacketWord64) :
    (uniformLabeledLaw S).p y =
      if packetWordSupport y = S then
        1 / (15 : ℝ) ^ S.card else 0 := by
  classical
  rw [uniformLabeledLaw, FinPMF.map_p]
  by_cases hy : packetWordSupport y = S
  · rw [if_pos hy]
    let yy : LabeledSupport S := ⟨y, hy⟩
    have h (p : LabeledSupport S) : (p.val = y) ↔ p = yy := by
      constructor
      · intro hp; exact Subtype.ext hp
      · intro hp; exact congrArg Subtype.val hp
    simp only [h, Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card (LabeledSupport S) : ℝ) = _
    rw [packetWordSupport_card]
    push_cast
    rfl
  · rw [if_neg hy]
    apply Finset.sum_eq_zero
    intro p _
    have hp : p.val ≠ y := by
      intro he
      apply hy
      rw [← he]
      exact p.property
    simp [hp]

/-- The point mass at `y` equals its support-fiber mass divided by the
number of possible nonzero labelings, provided `P` is constant on that
support fiber. -/
theorem support_fiber_mass
    (P : FinPMF PacketWord64)
    (hP : ∀ y z, packetWordSupport y = packetWordSupport z →
      P.p y = P.p z)
    (y : PacketWord64) :
    (P.map packetWordSupport).p (packetWordSupport y) =
      (15 : ℝ) ^ (packetWordSupport y).card * P.p y := by
  classical
  have hp : (P.map packetWordSupport).p (packetWordSupport y) =
      (P.map packetWordSupport).prob
        (fun S => S = packetWordSupport y) := by
    simp [FinPMF.prob, Finset.sum_filter]
  rw [hp, FinPMF.map_prob]
  unfold FinPMF.prob
  have hterm (x : PacketWord64)
      (hx : x ∈ Finset.univ.filter
        (fun x => packetWordSupport x = packetWordSupport y)) :
      P.p x = P.p y :=
    hP x y (Finset.mem_filter.mp hx).2
  rw [Finset.sum_congr rfl hterm, Finset.sum_const,
    nsmul_eq_mul]
  have hc :
      (Finset.univ.filter
        (fun x : PacketWord64 =>
          packetWordSupport x = packetWordSupport y)).card =
        Fintype.card (LabeledSupport (packetWordSupport y)) := by
    rw [Fintype.card_subtype]
  rw [hc, packetWordSupport_card]
  push_cast
  rfl

/-- Exact disintegration of a finite packet-word law that is uniform on
each support fiber. The mass of the support is retained; its labels are
then sampled uniformly. -/
theorem support_disintegration
    (P : FinPMF PacketWord64)
    (hP : ∀ y z, packetWordSupport y = packetWordSupport z →
      P.p y = P.p z) :
    P = (P.map packetWordSupport).bind uniformLabeledLaw := by
  classical
  apply FinPMF.ext
  funext y
  rw [FinPMF.bind_p]
  have hvalue (S : Finset (Fin 64)) :
      (uniformLabeledLaw S).p y =
        if S = packetWordSupport y then
          1 / (15 : ℝ) ^ S.card else 0 := by
    rw [uniformLabeledLaw_apply]
    by_cases h : S = packetWordSupport y
    · simp [h]
    · have hr : packetWordSupport y ≠ S := Ne.symm h
      simp [h, hr]
  simp_rw [hvalue, mul_ite, mul_zero]
  simp only [Finset.sum_ite_eq', Finset.mem_univ, if_pos]
  rw [support_fiber_mass P hP y]
  have hp : (15 : ℝ) ^ (packetWordSupport y).card ≠ 0 := by positivity
  field_simp

/-- The actual sampled RS8 row is a support mixture with uniform labels.
This uses one sampled coefficient tuple per realized code, not fresh
coefficients for each message in the distance count. -/
theorem sampledRS8RowLaw_disintegration
    (points : Fin 8 → Field16) (rows : RS8Group Field16) :
    sampledRS8RowLaw points rows =
      ((sampledRS8RowLaw points rows).map packetWordSupport).bind
        uniformLabeledLaw :=
  support_disintegration _
    (sampledRS8RowLaw_fiber_uniform points rows)

def rowShuffleKernel (x : PacketWord64) : FinPMF PacketWord64 :=
  (FinPMF.uniform (Equiv.Perm (Fin 64))).map
    (fun σ => shufflePacketWord σ x)

/-- The actual independent row-permutation experiment. -/
def shuffledPacketLaw (P : FinPMF PacketWord64) : FinPMF PacketWord64 :=
  (P.prod (FinPMF.uniform (Equiv.Perm (Fin 64)))).map
    (fun p => shufflePacketWord p.2 p.1)

theorem shuffledPacketLaw_eq_bind (P : FinPMF PacketWord64) :
    shuffledPacketLaw P = P.bind rowShuffleKernel := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [shuffledPacketLaw, FinPMF.expect_map, FinPMF.expect_prod,
    FinPMF.expect_bind]
  simp only [rowShuffleKernel, FinPMF.expect_map]

theorem uniformLabeledLaw_shuffle (S : Finset (Fin 64)) :
    (uniformLabeledLaw S).bind rowShuffleKernel =
      packetShuffleLaw S := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [FinPMF.expect_bind, uniformLabeledLaw, FinPMF.expect_map,
    packetShuffleLaw, FinPMF.expect_map, FinPMF.expect_prod]
  simp only [rowShuffleKernel, FinPMF.expect_map]

/-- Uniform labels within each source support turn the independently
sampled row shuffle into an exact mixture of fixed-support shuffle laws. -/
theorem shuffledPacketLaw_disintegration
    (P : FinPMF PacketWord64)
    (hP : ∀ y z, packetWordSupport y = packetWordSupport z →
      P.p y = P.p z) :
    shuffledPacketLaw P =
      (P.map packetWordSupport).bind packetShuffleLaw := by
  calc
    shuffledPacketLaw P = P.bind rowShuffleKernel :=
      shuffledPacketLaw_eq_bind P
    _ = (((P.map packetWordSupport).bind uniformLabeledLaw).bind
          rowShuffleKernel) := by rw [← support_disintegration P hP]
    _ = (P.map packetWordSupport).bind
          (fun S => (uniformLabeledLaw S).bind rowShuffleKernel) := by
      rw [FinPMF.bind_assoc]
    _ = _ := by simp_rw [uniformLabeledLaw_shuffle]

theorem sampledRS8Row_shuffledLaw
    (points : Fin 8 → Field16) (rows : RS8Group Field16) :
    shuffledPacketLaw (sampledRS8RowLaw points rows) =
      ((sampledRS8RowLaw points rows).map packetWordSupport).bind
        packetShuffleLaw :=
  shuffledPacketLaw_disintegration _
    (sampledRS8RowLaw_fiber_uniform points rows)

/-- After row shuffling, every complete target packet word of weight `v`
has the source's total weight-`v` probability divided by the size of that
value-carrying weight layer. -/
theorem shuffledPacketLaw_apply
    (P : FinPMF PacketWord64)
    (hP : ∀ y z, packetWordSupport y = packetWordSupport z →
      P.p y = P.p z)
    (y : PacketWord64) :
    (shuffledPacketLaw P).p y =
      P.prob (fun x => (packetWordSupport x).card =
        (packetWordSupport y).card) /
        ((Nat.choose 64 (packetWordSupport y).card : ℝ) *
          (15 : ℝ) ^ (packetWordSupport y).card) := by
  classical
  let Q := P.map packetWordSupport
  let v := (packetWordSupport y).card
  let c : ℝ := 1 / ((Nat.choose 64 v : ℝ) * (15 : ℝ) ^ v)
  have hmass (E : Finset (Fin 64) → Prop)
      [DecidablePred E] :
      Q.expect (fun S => if E S then c else 0) = Q.prob E * c := by
    rw [FinPMF.prob_eq_expect_indicator]
    unfold FinPMF.expect
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro S _
    by_cases h : E S <;> simp [h]
  rw [shuffledPacketLaw_disintegration P hP, FinPMF.bind_p]
  have hterm (S : Finset (Fin 64)) :
      (packetShuffleLaw S).p y =
        if S.card = v then c else 0 := by
    rw [packetShuffleLaw_apply]
    change (if v = S.card then
        1 / ((Nat.choose 64 S.card : ℝ) * (15 : ℝ) ^ S.card)
      else 0) = (if S.card = v then c else 0)
    by_cases h : S.card = v
    · have hr : v = S.card := h.symm
      simp only [if_pos hr, if_pos h]
      rw [h]
    · have hr : v ≠ S.card := Ne.symm h
      simp only [if_neg hr, if_neg h]
  simp_rw [hterm]
  change Q.expect (fun S => if S.card = v then c else 0) = _
  rw [hmass (fun S => S.card = v)]
  rw [FinPMF.map_prob]
  simp only [c, v]
  ring

#print axioms uniformLabeledLaw_apply
#print axioms support_fiber_mass
#print axioms support_disintegration
#print axioms sampledRS8RowLaw_disintegration
#print axioms sampledRS8Row_shuffledLaw
#print axioms shuffledPacketLaw_apply

end Spin.Turbo
