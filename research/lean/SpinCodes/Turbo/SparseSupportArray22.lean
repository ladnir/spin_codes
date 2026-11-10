import SpinCodes.Turbo.SparseExactCell22

/-!
# Exact support-class arrays under the shared Turbo regional route

Audit: `RS8ReferenceRegion22.turboReferenceArrayLaw22_transpose` reindexes
independent complete cells, and `turboReferenceRoutedLaw22_product` pairs
one regional word with one shared permutation. The same two arguments
apply to exact nonzero label laws with region-dependent active supports.
-/
noncomputable section
namespace Spin.Turbo
open Finset

theorem regionalNonzeroGiven22_eq_uniformLabeledLaw (T : Finset (Fin 64)) :
    regionalNonzeroGiven22 T = uniformLabeledLaw T := by
  apply FinPMF.ext
  funext x
  rw [regionalNonzeroGiven22_apply, uniformLabeledLaw_apply]
  rfl

def sparseSupportArrayLaw22 (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) : FinPMF (TurboOuterWord22 R) :=
  Spin.piPMF (fun i => regionalNonzeroGiven22 (supports i))

def sparseRegionalSupport22 {R : ℕ}
    (supports : Fin (R * 16) → Finset (Fin 64)) (region : Fin 64) :
    Finset (Fin (R * 16)) :=
  Finset.univ.filter (fun i => region ∈ supports i)

theorem sparseSupportArrayLaw22_transpose (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) :
    (sparseSupportArrayLaw22 R supports).map (packetTransposeEquiv22 R) =
      Spin.piPMF (fun region => regionalNonzeroGiven22 (sparseRegionalSupport22 supports region)) := by
  apply FinPMF.ext
  funext regions
  rw [FinPMF.map_equiv_apply]
  simp only [sparseSupportArrayLaw22, regionalNonzeroGiven22_eq_cell_product, Spin.piPMF_apply]
  have he (i : Fin (R * 16)) (j : Fin 64) :
      sparseExactCellLaw22 (supports i) j =
        sparseExactCellLaw22 (sparseRegionalSupport22 supports j) i := by
    simp only [sparseExactCellLaw22, sparseRegionalSupport22, Finset.mem_filter,
      Finset.mem_univ, true_and]
  simp_rw [he]
  change (∏ i : Fin (R * 16), ∏ j : Fin 64,
      (sparseExactCellLaw22 (sparseRegionalSupport22 supports j) i).p (regions j i)) = _
  exact Finset.prod_comm

private theorem finite_map_prod_swap_bind {A B C : Type*}
    [Fintype A] [Fintype B] [Fintype C]
    (P : FinPMF A) (Q : FinPMF B) (f : A → B → C) :
    (P.prod Q).map (fun pair => f pair.1 pair.2) =
      Q.bind (fun b => P.map (fun a => f a b)) := by
  apply FinPMF.eq_of_expect_eq
  intro g
  simp only [FinPMF.expect_map, FinPMF.expect_prod, FinPMF.expect_bind]
  unfold FinPMF.expect
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro b _
  apply Finset.sum_congr rfl
  intro a _
  ring

def sparseSupportRoutedLaw22 (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) :
    FinPMF (Fin 64 → Fin (R * 16) → Packet4) :=
  ((sparseSupportArrayLaw22 R supports).prod (turboRegionPermLaw22 R)).map
    (fun seed => turboRouteRegions22 seed.1 seed.2)

/-- Complete exact values are independent between regions. Within each
region the same one permutation acts on all selected group positions. -/
theorem sparseSupportRoutedLaw22_product (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) :
    sparseSupportRoutedLaw22 R supports =
      Spin.piPMF (fun region => sparseExactRegionLaw22 R (sparseRegionalSupport22 supports region)) := by
  let P : Fin 64 → FinPMF (Fin (R * 16) → Packet4) :=
    fun region => regionalNonzeroGiven22 (sparseRegionalSupport22 supports region)
  let Q : Fin 64 → FinPMF (Equiv.Perm (Fin (R * 16))) :=
    fun _ => FinPMF.uniform (Equiv.Perm (Fin (R * 16)))
  let f : ∀ _ : Fin 64,
      ((Fin (R * 16) → Packet4) × Equiv.Perm (Fin (R * 16))) → (Fin (R * 16) → Packet4) :=
    fun _ pair i => pair.1 (pair.2 i)
  have hfirst : (sparseSupportArrayLaw22 R supports).map (packetTransposeEquiv22 R) =
      Spin.piPMF P := sparseSupportArrayLaw22_transpose R supports
  calc
    _ = (((sparseSupportArrayLaw22 R supports).map (packetTransposeEquiv22 R)).prod
        (turboRegionPermLaw22 R)).map (fun seed j i => seed.1 j (seed.2 j i)) := by
      apply FinPMF.eq_of_expect_eq
      intro g
      simp only [sparseSupportRoutedLaw22, FinPMF.expect_map, FinPMF.expect_prod]
      rfl
    _ = (((Spin.piPMF P).prod (Spin.piPMF Q)).map piPairEquiv).map
        (fun pairs j => f j (pairs j)) := by
      rw [hfirst]
      change ((Spin.piPMF P).prod (Spin.piPMF Q)).map _ = _
      rw [FinPMF.map_comp]
      rfl
    _ = (Spin.piPMF (fun j => (P j).prod (Q j))).map (fun pairs j => f j (pairs j)) := by
      rw [piPMF_prod_pair]
    _ = Spin.piPMF (fun j => ((P j).prod (Q j)).map (f j)) := Spin.piPMF_map _ _
    _ = _ := by
      apply congrArg Spin.piPMF
      funext j
      exact finite_map_prod_swap_bind (P j) (Q j)
        (fun word σ => regionalValueShuffleEquiv22 σ word)

#print axioms sparseSupportRoutedLaw22_product
end Spin.Turbo
