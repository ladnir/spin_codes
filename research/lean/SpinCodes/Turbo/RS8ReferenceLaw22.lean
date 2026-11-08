import SpinCodes.Turbo.RS8ReferenceMoment22

/-!
# Normalized complete-word reference for the Turbo outer

An active group is uniform on all 64 four-bit packet values, including
zero. An inactive group retains its exact sampled zero-message law. The
independent product is the reference whose point density was used in
`RS8UniformReference22`.
-/

noncomputable section

namespace Spin.Turbo

theorem packetWord64_card : Fintype.card PacketWord64 = 2 ^ 256 := by
  change Fintype.card (Fin 64 → Fin 4 → ZMod 2) = _
  simp only [Fintype.card_fun, Fintype.card_fin, ZMod.card]
  norm_num [← pow_mul]

def turboReferenceRowLaw22 (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] : FinPMF PacketWord64 :=
  if active then FinPMF.uniform PacketWord64
  else shuffledPacketLaw (sampledRS8RowLaw points 0)

theorem turboReferenceRowLaw22_apply (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) :
    (turboReferenceRowLaw22 points active).p word =
      turboReferenceRowDensity22 points active word := by
  by_cases h : active
  · simp only [turboReferenceRowLaw22, turboReferenceRowDensity22, h,
      ite_true]
    change 1 / (Fintype.card PacketWord64 : ℝ) = _
    rw [packetWord64_card]
    norm_num
  · simp only [turboReferenceRowLaw22, turboReferenceRowDensity22, h,
      ite_false]

def turboReferenceArrayLaw22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) : FinPMF (TurboOuterWord22 R) :=
  Spin.piPMF (fun i => turboReferenceRowLaw22 points (i ∈ S))

theorem turboReferenceArrayLaw22_apply (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) :
    (turboReferenceArrayLaw22 R points S).p words =
      ∏ i : Fin (R * 16),
        turboReferenceRowDensity22 points (i ∈ S) (words i) := by
  simp only [turboReferenceArrayLaw22, Spin.piPMF_apply,
    turboReferenceRowLaw22_apply]

/-- The reference is a genuine probability law, and the original
active-set first moment is bounded by its expectation of the same
regional-route and inner failure event. -/
theorem turboSetup22_activeSet_moment_le_reference_expect
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (d : ℕ) :
    (turboOuterLaw22 R).expect (fun outerSeed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d
          (turboOuterWord22 points outerSeed msg)) ≤
    turboOldUniformFactor22 ^ S.card *
      (turboReferenceArrayLaw22 R points S).expect
        (fun words => turboRemainingFailure22 words d) := by
  calc
    _ ≤ turboOldUniformFactor22 ^ S.card *
          ∑ words : TurboOuterWord22 R,
            (∏ i : Fin (R * 16),
              turboReferenceRowDensity22 points (i ∈ S) (words i)) *
              turboRemainingFailure22 words d :=
        turboSetup22_activeSet_moment_le_reference R points hpoints S d
    _ = _ := by
      simp only [FinPMF.expect, turboReferenceArrayLaw22_apply]

#print axioms turboReferenceArrayLaw22_apply
#print axioms turboSetup22_activeSet_moment_le_reference_expect

end Spin.Turbo
