import SpinCodes.Turbo.RS8ReferenceLaw22

/-!
# The inactive Turbo reference row

An inactive RS8 group outputs the all-zero packet word for every sampled
column multiplier and row permutation. This is the first local fact needed
to identify the reference's regional selected-set and label law.
-/

noncomputable section

namespace Spin.Turbo

theorem rs8Bits_zero (points : Fin 8 → Field16) (j : Fin 8) :
    rs8Bits points (0 : RS8Group Field16) j = 0 := by
  have hc : rs8Column points (0 : RS8Group Field16) j = 0 := by
    funext r
    change (0 : Polynomial Field16).eval (points j) = 0
    simp
  simp [rs8Bits, hc]

theorem sampledRS8RowWord_zero (points : Fin 8 → Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    sampledRS8RowWord points 0 coeff = 0 := by
  funext k
  change sampledPacketColumn
    (rs8Bits points 0 ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm k).1)
    (coeff ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm k).1)
    ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm k).2 = 0
  rw [rs8Bits_zero, sampledPacketColumn_zero]
  rfl

theorem sampledRS8PacketColumns_zero (points : Fin 8 → Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    sampledRS8PacketColumns points 0 coeff = 0 := by
  funext j
  change sampledPacketColumn (rs8Bits points 0 j) (coeff j) = 0
  rw [rs8Bits_zero, sampledPacketColumn_zero]

theorem finite_map_const_apply
    {α β : Type*} [Fintype α] [Fintype β] [DecidableEq β]
    (P : FinPMF α) (value target : β) :
    (P.map (fun _ => value)).p target =
      if target = value then 1 else 0 := by
  classical
  rw [FinPMF.map_p]
  by_cases h : target = value
  · subst target
    simp only [ite_true]
    exact P.total
  · have hne : value ≠ target := Ne.symm h
    simp [hne]
    exact h

theorem sampledRS8RowLaw_zero_apply (points : Fin 8 → Field16)
    (word : PacketWord64) :
    (sampledRS8RowLaw points 0).p word =
      if word = 0 then 1 else 0 := by
  let P := Spin.piPMF (fun _ : Fin 8 =>
    FinPMF.uniform {c : Field32 // c ≠ 0})
  have he : sampledRS8RowLaw points 0 =
      P.map (fun _ => (0 : PacketWord64)) := by
    change (P.map (sampledRS8PacketColumns points 0)).map
      packetColumnsWordEquiv = _
    rw [FinPMF.map_comp]
    have hf : (packetColumnsWordEquiv ∘
        sampledRS8PacketColumns points 0) =
        (fun _ => (0 : PacketWord64)) := by
      funext coeff
      exact sampledRS8RowWord_zero points coeff
    exact congrArg (fun f => P.map f) hf
  rw [he]
  exact finite_map_const_apply P 0 word

theorem shufflePacketWord_zero (σ : Equiv.Perm (Fin 64)) :
    shufflePacketWord σ (0 : PacketWord64) = 0 := by
  rfl

theorem shuffledRS8RowLaw_zero_apply (points : Fin 8 → Field16)
    (word : PacketWord64) :
    (shuffledPacketLaw (sampledRS8RowLaw points 0)).p word =
      if word = 0 then 1 else 0 := by
  rw [shuffledPacketLaw_apply_inverse]
  have hz (σ : Equiv.Perm (Fin 64)) :
      shufflePacketWord σ.symm word = 0 ↔ word = 0 := by
    constructor
    · intro h
      have hh := congrArg (shufflePacketWord σ) h
      have he : shufflePacketWord σ
          (shufflePacketWord σ.symm word) = word :=
        (packetWordShuffleEquiv σ).right_inv word
      simpa only [he, shufflePacketWord_zero] using hh
    · intro h
      rw [h]
      exact shufflePacketWord_zero σ.symm
  have hf :
      (fun σ : Equiv.Perm (Fin 64) =>
        (sampledRS8RowLaw points 0).p
          (shufflePacketWord σ.symm word)) =
      (fun _ => if word = 0 then 1 else 0) := by
    funext σ
    rw [sampledRS8RowLaw_zero_apply]
    simp only [hz σ]
  rw [hf, FinPMF.expect_const]

#print axioms sampledRS8RowWord_zero
#print axioms sampledRS8RowLaw_zero_apply
#print axioms shuffledRS8RowLaw_zero_apply

end Spin.Turbo
