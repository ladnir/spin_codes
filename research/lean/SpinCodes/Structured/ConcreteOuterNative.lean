import SpinCodes.Structured.ConcreteOuter
import SpinCodes.Structured.Schedule

noncomputable section
namespace Spin.Structured.ConcreteOuter

/-- Number of Golay blocks in one native constituent. -/
def nativeBlocks (m : ℕ) : ℕ := bsched m / 24

theorem native_width (m : ℕ) : nativeBlocks m * 24 = bsched m :=
  Nat.div_mul_cancel (twentyfour_dvd_bsched m)

theorem native_dimension (m : ℕ) : Lsched m * (nativeBlocks m * 12) = Nsched m / 2 := by
  have hw := native_width m
  have h : Lsched m * (nativeBlocks m * 12) * 2 = Nsched m := by
    rw [Nsched_eq, ← hw]
    ring
  omega

theorem native_round_divisibility (m : ℕ) : 128 ∣ bsched m * Lsched m := by
  unfold Lsched
  exact ⟨bsched m * (m + 1), by ring⟩

/-- Flattening the row, Golay-block, and bit indices is an exact message bijection. -/
def messageEquiv (L k : ℕ) :
    (Fin (L * (k * 12)) → Bool) ≃ (Fin L → LocalMessage k) where
  toFun x i j t := x (finProdFinEquiv (i, finProdFinEquiv (j,t)))
  invFun x n := x (finProdFinEquiv.symm n).1
    (finProdFinEquiv.symm (finProdFinEquiv.symm n).2).1
    (finProdFinEquiv.symm (finProdFinEquiv.symm n).2).2
  left_inv x := by
    funext n
    simp only [Prod.mk.eta, Equiv.apply_symm_apply]
  right_inv x := by
    funext i j t
    simp

abbrev NativeMessage (m : ℕ) := Fin (Lsched m) → LocalMessage (nativeBlocks m)
abbrev NativeSeed (m : ℕ) := Seed (nativeBlocks m)

def nativeMessageEquiv (m : ℕ) : (Fin (Nsched m / 2) → Bool) ≃ NativeMessage m :=
  (Equiv.arrowCongr (finCongr (native_dimension m).symm) (Equiv.refl Bool)).trans
    (messageEquiv (Lsched m) (nativeBlocks m))

def nativeSeedLaw (m : ℕ) : FinPMF (NativeSeed m) := seedLaw (nativeBlocks m)

/-- Actual shared Golay-BAA rows, indexed by the native paper width. -/
def nativeRows (m : ℕ) (seed : NativeSeed m) (message : Fin (Nsched m / 2) → Bool) :
    Fin (Lsched m) → Finset (Fin (bsched m)) := fun i =>
  (rowSupports seed (nativeMessageEquiv m message) i).map (finCongr (native_width m)).toEmbedding

theorem native_rows_injective (m : ℕ) (seed : NativeSeed m) :
    Function.Injective (nativeRows m seed) := by
  intro x y h
  apply (nativeMessageEquiv m).injective
  apply rowSupports_injective seed
  funext i
  exact Finset.map_injective (finCongr (native_width m)).toEmbedding (congrFun h i)

@[simp] theorem native_message_zero (m : ℕ) : nativeMessageEquiv m 0 = 0 := by
  rfl

theorem native_activeRows (m : ℕ) (seed : NativeSeed m)
    (message : Fin (Nsched m / 2) → Bool) :
    ConcreteRoute.activeRows (nativeRows m seed message) =
      occupation (nativeMessageEquiv m message) := by
  unfold ConcreteRoute.activeRows nativeRows
  simp only [ne_eq, Finset.map_eq_empty, rowSupports_empty_iff]
  rfl

theorem native_totalWeight (m : ℕ) (seed : NativeSeed m)
    (message : Fin (Nsched m / 2) → Bool) :
    ConcreteRoute.totalWeight (nativeRows m seed message) =
      ∑ i, wtF (encode seed (nativeMessageEquiv m message i)) := by
  simp only [ConcreteRoute.totalWeight, nativeRows, Finset.card_map, rowSupports, support_card]

theorem native_occupation_mem (m : ℕ) (seed : NativeSeed m)
    (message : Fin (Nsched m / 2) → Bool) (hne : message ≠ 0) :
    ConcreteRoute.activeRows (nativeRows m seed message) ∈ Finset.Ico 1 (Lsched m + 1) := by
  rw [native_activeRows, Finset.mem_Ico]
  have hm : nativeMessageEquiv m message ≠ 0 := by
    intro h
    exact hne ((nativeMessageEquiv m).injective (h.trans (native_message_zero m).symm))
  exact ⟨occupation_pos _ hm, Nat.lt_succ_of_le (occupation_le _)⟩

end Spin.Structured.ConcreteOuter

