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

theorem shufflePacketWord_zero (σ : Equiv.Perm (Fin 64)) :
    shufflePacketWord σ (0 : PacketWord64) = 0 := by
  rfl

#print axioms sampledRS8RowWord_zero

end Spin.Turbo
