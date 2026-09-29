import SpinCodes.Structured.ConcreteReshape

namespace Spin.Structured.ConcreteRoute

/-- Serialize regions in increasing order and split the stream into 128-bit blocks. -/
def regionMajor {L b R : ℕ} (h : b * L = R * 128) :
    (Fin b × Fin L) ≃ (Fin R × Fin 128) :=
  finProdFinEquiv.trans ((finCongr h).trans finProdFinEquiv.symm)

theorem regionMajor_position {L b R : ℕ} (h : b * L = R * 128)
    (j : Fin b) (k : Fin L) :
    (regionMajor h (j, k)).2.val + 128 * (regionMajor h (j, k)).1.val =
      k.val + L * j.val := by
  have he : finProdFinEquiv (regionMajor h (j, k)) =
      finCongr h (finProdFinEquiv (j, k)) :=
    finProdFinEquiv.apply_symm_apply _
  exact congrArg Fin.val he

def streamWiring {L b : ℕ} (h : 128 ∣ b * L) :
    (Fin b × Fin L) ≃ (Fin (b * L / 128) × Fin 128) :=
  regionMajor (Nat.div_mul_cancel h).symm

end Spin.Structured.ConcreteRoute
