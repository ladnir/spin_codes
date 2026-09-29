import SpinCodes.Structured.ConcreteFixedRoutedProfile

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteRoute
attribute [local instance] Classical.propDecidable

/-- Embed active rows and fill all other positions by the empty support. -/
def embedRows {Q L b : ℕ} (e : Fin Q ↪ Fin L) (rows : Fin Q → Finset (Fin b)) :
    Fin L → Finset (Fin b) := transpose (fun j => (transpose rows j).map e)

@[simp] theorem embedRows_apply {Q L b : ℕ} (e : Fin Q ↪ Fin L)
    (rows : Fin Q → Finset (Fin b)) (i : Fin Q) : embedRows e rows (e i)=rows i := by
  ext j
  simp [embedRows]

theorem embedRows_off {Q L b : ℕ} (e : Fin Q ↪ Fin L)
    (rows : Fin Q → Finset (Fin b)) (j : Fin L) (hj : j∉Set.range e) : embedRows e rows j=∅ := by
  ext k
  simp only [embedRows,mem_transpose,mem_map,Finset.notMem_empty,iff_false,not_exists,not_and]
  intro i _ hi
  exact hj ⟨i,hi⟩

theorem embedRows_injective {Q L b : ℕ} (e : Fin Q ↪ Fin L) : Function.Injective (embedRows (b:=b) e) := by
  intro rows ys h
  funext i
  simpa only [embedRows_apply] using congrFun h (e i)

theorem embedRows_extract {Q L b : ℕ} (e : Fin Q ↪ Fin L) (rows : Fin L → Finset (Fin b))
    (h : ∀ j, j∉Set.range e → rows j=∅) : embedRows e (fun i => rows (e i))=rows := by
  funext j
  by_cases hj : j∈Set.range e
  · obtain ⟨i,rfl⟩ := hj
    exact embedRows_apply e _ i
  · rw [embedRows_off e _ j hj,h j hj]

theorem rowLaw_embed_mass {Q L b : ℕ} (e : Fin Q ↪ Fin L)
    (rows ys : Fin Q → Finset (Fin b)) :
    (rowLaw (embedRows e rows)).p (embedRows e ys)=(rowLaw rows).p ys := by
  change (∏ j, (shuffleLaw (embedRows e rows j)).p (embedRows e ys j)) = ∏ i, (shuffleLaw (rows i)).p (ys i)
  have hr : (∏ j, (shuffleLaw (embedRows e rows j)).p (embedRows e ys j)) =
      ∏ j ∈ univ.map e, (shuffleLaw (embedRows e rows j)).p (embedRows e ys j) := by
    symm
    apply prod_subset (subset_univ _)
    intro j _ hj
    have he : j∉Set.range e := by simpa only [mem_map,mem_univ,true_and,Set.mem_range] using hj
    rw [embedRows_off e _ j he,embedRows_off e _ j he]
    simp [shuffleLaw_apply]
  rw [hr,prod_map]
  simp only [embedRows_apply]

/-- Independent shuffling commutes with insertion of inactive zero rows. -/
theorem rowLaw_embed {Q L b : ℕ} (e : Fin Q ↪ Fin L) (rows : Fin Q → Finset (Fin b)) :
    rowLaw (embedRows e rows)=(rowLaw rows).map (embedRows e) := by
  ext ys
  by_cases h : embedRows e (fun i => ys (e i))=ys
  · rw [←h,rowLaw_embed_mass,FinPMF.map_p]
    have he (x : Fin Q → Finset (Fin b)) :
        embedRows e x=embedRows e (fun i => ys (e i)) ↔ x=(fun i => ys (e i)) :=
      (embedRows_injective e).eq_iff
    simp only [he,sum_ite_eq',mem_univ,ite_true]
  · have hex : ∃ j, j∉Set.range e ∧ ys j≠∅ := by
      by_contra hh
      apply h
      apply embedRows_extract
      simpa only [not_exists,not_and,not_not] using hh
    obtain ⟨j,hj,hy⟩ := hex
    have hl : (rowLaw (embedRows e rows)).p ys=0 := by
      change (∏ i, (shuffleLaw (embedRows e rows i)).p (ys i))=0
      apply prod_eq_zero (mem_univ j)
      rw [embedRows_off e _ j hj,shuffleLaw_apply]
      simp [hy]
    rw [hl,FinPMF.map_p]
    symm
    apply sum_eq_zero
    intro x _
    apply if_neg
    intro hx
    apply h
    have he : (fun i => ys (e i))=x := by
      funext i
      rw [←hx,embedRows_apply]
    rw [he,hx]

#print axioms rowLaw_embed
end Spin.Structured.Placement
