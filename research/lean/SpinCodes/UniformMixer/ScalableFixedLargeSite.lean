import SpinCodes.UniformMixer.ScalableFixedTwoComparison
import SpinCodes.Structured.ConcreteFixedInsertionNorm

/-!
The Scalable continuum uses the old insertion geometry with a different
two-state impulse. These identities are valid for every number of marks.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement
open MeasureTheory

attribute [local instance] Classical.propDecidable

def scalableImpulseFin : Matrix (Fin 2) (Fin 2) ℝ :=
  toFinTwo coarseP

def scalableSiteProduct (γ : ℝ) {a : ℕ} (x : Fin a → ℝ) :
    Matrix (Fin 2) (Fin 2) ℝ :=
  sitePath γ scalableImpulseFin 0 (List.ofFn x) 1

theorem toFinTwo_timeRegionProduct {a : ℕ} (γ : ℝ)
    (gaps : Fin a → ℝ) (last : ℝ) :
    toFinTwo (timeRegionProduct γ gaps last) =
      (List.ofFn (fun i =>
        timeEmpty γ (gaps i) * scalableImpulseFin)).prod *
        timeEmpty γ last := by
  induction a with
  | zero =>
      simp [timeRegionProduct, toFinTwo_gap]
  | succ a ih =>
      simp only [timeRegionProduct, toFinTwo_mul,
        toFinTwo_gap, toFinTwo_impulse,
        List.ofFn_succ, List.prod_cons, ih,
        scalableImpulseFin, Fin.tail, Matrix.mul_assoc]

theorem scalableSiteProduct_eq (γ : ℝ) {a : ℕ}
    (x : Fin a → ℝ) :
    toFinTwo (continuumSiteProduct γ x) =
      scalableSiteProduct γ x := by
  rw [continuumSiteProduct, scalableSiteProduct,
    toFinTwo_timeRegionProduct]
  rw [← factorSitePath_constant, List.map_ofFn,
    factorSitePath_ofFn]
  rfl

theorem scalableImpulseFin_le_plus (i j : Fin 2) :
    scalableImpulseFin i j ≤ Spin.Pplus (1 / 524287) i j :=
  scalable_impulse_le_plus i j

theorem scalableContinuumTwo_site (θ : ℝ) (a : ℕ)
    (i j : Fin 2) :
    toFinTwo (scalableContinuumTwo θ a) i j =
      (a.factorial : ℝ) *
        ∫ x in orderedSiteDomain a,
          scalableSiteProduct (regionGamma θ) x i j := by
  rw [toFinTwo_apply]
  unfold scalableContinuumTwo
  congr 1
  apply setIntegral_congr_fun (orderedSiteDomain_measurable a)
  intro x hx
  change continuumSiteProduct (regionGamma θ) x
    (finTwoEquiv i) (finTwoEquiv j) =
      scalableSiteProduct (regionGamma θ) x i j
  rw [← scalableSiteProduct_eq, toFinTwo_apply]

theorem scalable_selectedGapProduct_eq_subtuple {Q : ℕ}
    (γ : ℝ) (x : Fin Q → ℝ) (A : Finset (Fin Q)) :
    selectedGapProduct γ scalableImpulseFin
      (fun i => positionGaps x i.castSucc)
      (positionGaps x (Fin.last Q)) A =
    scalableSiteProduct γ
      (fun i : Fin A.card => x (A.orderEmbOfFin rfl i)) := by
  rw [selectedGapProduct_eq_selectedPath,
    selectedSites_ofFn]
  rfl

theorem scalable_selectedGapProduct_normalized {R Q a : ℕ}
    (γ : ℝ) (B : BlockSubset R Q)
    (A : Finset (Fin Q)) (hA : A.card = a) :
    selectedGapProduct γ scalableImpulseFin
      (fun i => positionGaps (normalizedSites B) i.castSucc)
      (positionGaps (normalizedSites B) (Fin.last Q)) A =
    scalableSiteProduct γ (normalizedSites
      (selectedBlocks B A hA)) := by
  subst a
  rw [scalable_selectedGapProduct_eq_subtuple,
    normalizedSites_selectedBlocks]

theorem scalable_fullSiteProduct_normalized_permuted {R Q : ℕ}
    (γ : ℝ) (B : BlockSubset R Q)
    (u : Fin Q → ℝ) (π : Equiv.Perm (Fin Q)) :
    fullSiteProduct γ scalableImpulseFin (normalizedSites B)
      (fun i => u (π.symm i)) =
    ∑ A : Finset (Fin Q), markCoefficient u A •
      scalableSiteProduct γ (normalizedSites
        (picked A B π)) := by
  rw [fullSiteProduct_permuted]
  apply Finset.sum_congr rfl
  intro A hA
  rw [scalable_selectedGapProduct_normalized γ B
    (A.map π.toEmbedding) (Finset.card_map _)]
  rw [selectedBlocks_permuted_eq_picked]

def scalableFiniteInsertionAverage (R : ℕ) {Q : ℕ}
    (γ : ℝ) (u : Fin Q → ℝ) (i j : Fin 2) : ℝ :=
  (∑ B : BlockSubset R Q, ∑ π : Equiv.Perm (Fin Q),
    fullSiteProduct γ scalableImpulseFin (normalizedSites B)
      (fun k => u (π.symm k)) i j) /
        ((R.choose Q : ℝ) * (Q.factorial : ℝ))

theorem scalableFiniteInsertionAverage_eq {R Q : ℕ}
    (hQR : Q ≤ R) (γ : ℝ) (u : Fin Q → ℝ)
    (i j : Fin 2) :
    scalableFiniteInsertionAverage R γ u i j =
      ∑ A : Finset (Fin Q), markCoefficient u A *
        ((∑ C : BlockSubset R A.card,
          scalableSiteProduct γ (normalizedSites C) i j) /
          (R.choose A.card : ℝ)) := by
  unfold scalableFiniteInsertionAverage
  simp only [scalable_fullSiteProduct_normalized_permuted,
    Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]
  have hs :
      (∑ B : BlockSubset R Q,
        ∑ π : Equiv.Perm (Fin Q),
          ∑ A : Finset (Fin Q),
            markCoefficient u A *
              scalableSiteProduct γ
                (normalizedSites (picked A B π)) i j) =
      ∑ A : Finset (Fin Q),
        ∑ B : BlockSubset R Q,
          ∑ π : Equiv.Perm (Fin Q),
            markCoefficient u A *
              scalableSiteProduct γ
                (normalizedSites (picked A B π)) i j := by
    calc
      _ = ∑ B : BlockSubset R Q,
          ∑ A : Finset (Fin Q),
            ∑ π : Equiv.Perm (Fin Q),
              markCoefficient u A *
                scalableSiteProduct γ
                  (normalizedSites (picked A B π)) i j := by
        apply Finset.sum_congr rfl
        intro B hB
        exact Finset.sum_comm
      _ = _ := Finset.sum_comm
  rw [hs, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro A hA
  simp only [← Finset.mul_sum, mul_div_assoc]
  exact congrArg (fun z => markCoefficient u A * z)
    (picked_blockSubset_average hQR A
      (fun C => scalableSiteProduct γ (normalizedSites C) i j))

end Spin.UniformMixer
