import SpinCodes.Turbo.RS8Spectrum
import Mathlib.Data.Finset.Sort

/-!
# Ordering arbitrary RS8 zero-column sets

The shortening checker accepts an ordered prefix of at most four distinct
columns. An arbitrary set of at most four columns can be put in that
prefix by a permutation of the eight positions. This is a positional
fact only; it does not permute the code or change its setup law.
-/

noncomputable section

namespace Spin.Turbo

private def complSubtypeEquiv (z : Finset (Fin 8)) :
    ↥(zᶜ) ≃ {j : Fin 8 // j ∉ z} where
  toFun x := ⟨x.val, Finset.mem_compl.mp x.property⟩
  invFun x := ⟨x.val, Finset.mem_compl.mpr x.property⟩
  left_inv _x := Subtype.ext rfl
  right_inv _x := Subtype.ext rfl

/-- Order the members of `z` first, then the other positions. -/
noncomputable def rs8PositionOrder (z : Finset (Fin 8)) :
    Fin 8 ≃ Fin 8 := by
  classical
  let m := z.card
  have hm : m ≤ 8 := by simpa [m] using Finset.card_le_card (Finset.subset_univ z)
  have hsum : m + (8-m) = 8 := Nat.add_sub_of_le hm
  have hcomp : zᶜ.card = 8-m := by
    simp [Finset.card_compl, m]
  let ez : Fin m ≃ z := (z.orderIsoOfFin rfl).toEquiv
  let ec : Fin (8-m) ≃ {j : Fin 8 // j ∉ z} :=
    (zᶜ.orderIsoOfFin hcomp).toEquiv.trans (complSubtypeEquiv z)
  let esum : Fin m ⊕ Fin (8-m) ≃ Fin 8 :=
    (ez.sumCongr ec).trans (Equiv.sumCompl (fun j : Fin 8 => j ∈ z))
  exact (finCongr hsum.symm).trans (finSumFinEquiv.symm.trans esum)

/-- The first `|z|` positions of the order enumerate exactly `z`. -/
theorem rs8PositionOrder_prefix (z : Finset (Fin 8)) :
    z = Finset.univ.image (fun i : Fin z.card =>
      rs8PositionOrder z ⟨i.val, by
        have := Finset.card_le_card (Finset.subset_univ z)
        simp at this
        omega⟩) := by
  classical
  symm
  apply Finset.eq_of_subset_of_card_le
  · intro j hj
    obtain ⟨i, _, rfl⟩ := Finset.mem_image.mp hj
    let m := z.card
    have hm : m ≤ 8 := by simpa [m] using Finset.card_le_card (Finset.subset_univ z)
    have hsum : m + (8-m) = 8 := Nat.add_sub_of_le hm
    have hcomp : zᶜ.card = 8-m := by simp [Finset.card_compl, m]
    let ez : Fin m ≃ z := (z.orderIsoOfFin rfl).toEquiv
    let ec : Fin (8-m) ≃ {j : Fin 8 // j ∉ z} :=
      (zᶜ.orderIsoOfFin hcomp).toEquiv.trans (complSubtypeEquiv z)
    change ((ez.sumCongr ec).trans (Equiv.sumCompl (fun j : Fin 8 => j ∈ z)))
      ((finSumFinEquiv).symm ((finCongr hsum.symm) ⟨i.val, by omega⟩)) ∈ z
    have hi : (finCongr hsum.symm) ⟨i.val, by omega⟩ =
        Fin.castAdd (8-m) i := Fin.ext rfl
    rw [hi, finSumFinEquiv_symm_apply_castAdd]
    exact (ez i).property
  · rw [Finset.card_image_of_injective]
    · simp
    · exact (rs8PositionOrder z).injective.comp (fun i j hij => Fin.ext (Fin.mk.inj hij))

def rs8ZeroOnSet {F : Type*} [Field F] (points : Fin 8 → F)
    (z : Finset (Fin 8)) (rows : RS8Group F) : Prop :=
  ∀ j ∈ z, rs8Column points rows j = 0

noncomputable instance {F : Type*} [Field F] [Fintype F]
    (points : Fin 8 → F) (z : Finset (Fin 8)) :
    Fintype {rows : RS8Group F // rs8ZeroOnSet points z rows} :=
  Fintype.ofFinite _

/-- Shortening on any set of at most four columns depends only on its size. -/
theorem rs8ZeroOnSet_card {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (z : Finset (Fin 8)) (hz : z.card ≤ 4) :
    Fintype.card {rows : RS8Group F // rs8ZeroOnSet points z rows} =
      (Fintype.card F) ^ (8 * (4-z.card)) := by
  classical
  let positions : Fin 4 → Fin 8 :=
    fun i => rs8PositionOrder z ⟨i.val, by omega⟩
  have hpositions : Function.Injective positions := by
    intro i j hij
    exact Fin.ext (Fin.mk.inj ((rs8PositionOrder z).injective hij))
  have hiff (rows : RS8Group F) :
      rs8ZeroOnSet points z rows ↔
      rs8ZeroOnPrefix z.card hz points positions rows := by
    constructor
    · intro h i
      apply h
      have hmem : positions ⟨i.val, by omega⟩ ∈
          Finset.univ.image (fun k : Fin z.card =>
            rs8PositionOrder z ⟨k.val, by
              have := Finset.card_le_card (Finset.subset_univ z)
              simp at this
              omega⟩) :=
        Finset.mem_image.mpr ⟨i, Finset.mem_univ _, rfl⟩
      exact Eq.mp (congrArg
        (fun s : Finset (Fin 8) => positions ⟨i.val, by omega⟩ ∈ s)
        (rs8PositionOrder_prefix z).symm) hmem
    · intro h j hj
      rw [rs8PositionOrder_prefix z] at hj
      obtain ⟨i, _, rfl⟩ := Finset.mem_image.mp hj
      exact h i
  let e : {rows : RS8Group F // rs8ZeroOnSet points z rows} ≃
      {rows : RS8Group F // rs8ZeroOnPrefix z.card hz points positions rows} :=
    { toFun := fun x => ⟨x.val, (hiff x.val).mp x.property⟩
      invFun := fun x => ⟨x.val, (hiff x.val).mpr x.property⟩
      left_inv := fun _ => Subtype.ext rfl
      right_inv := fun _ => Subtype.ext rfl }
  rw [Fintype.card_congr e]
  exact rs8Shortening_card z.card hz points hpoints positions hpositions

#print axioms rs8PositionOrder_prefix
#print axioms rs8ZeroOnSet_card

end Spin.Turbo
