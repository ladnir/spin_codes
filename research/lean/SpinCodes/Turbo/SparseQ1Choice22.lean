import SpinCodes.Turbo.SparseOccupation22

/-! The retained q1 witnesses choose probability1 on some shells and a
Chernoff bound on the others. The choice is made inside the exact original
support measure, before its unconditional sum. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

def sparseOneChosenShellBound22 (points : Fin 8 → Field16)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (R d : ℕ)
    (z : ℕ → ℝ) (useChernoff : ℕ → Bool) : ℝ :=
  ∑ v ∈ Finset.range 65, if useChernoff v then
    ((rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ)) / (z v) ^ d) *
      ∑ k : Fin 11, orderedMatrixCoefficient sparseBoolDegree22
        (sparseOneRegionChoice22 (B v) R) 64 v 0 k
    else rs8ExpectedPacketCount points v

private theorem sparseSupportFailure22_le_one (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) (d : ℕ) :
    (sparseSupportArrayLaw22 R supports).expect
      (fun words => turboRemainingFailure22 words d) ≤ 1 := by
  calc
    _ ≤ (sparseSupportArrayLaw22 R supports).expect (fun _ => 1) := by
      apply FinPMF.expect_mono
      intro words
      exact FinPMF.prob_le_one _ _
    _ = _ := FinPMF.expect_const _ _

theorem turboSetup22_one_activeSet_le_chosen_shells
    {R : ℕ} (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin 1 ≃ {i : Fin (R * 16) // i ∈ S}) (d : ℕ)
    (z : ℕ → ℝ) (hz0 : ∀ v, 0 < z v) (hz1 : ∀ v, z v ≤ 1)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ v j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j (z v))) (B v j))
    (useChernoff : ℕ → Bool) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      sparseOneChosenShellBound22 points B R d z useChernoff := by
  let F := fun supports : Fin (R * 16) → Finset (Fin 64) =>
    (sparseSupportArrayLaw22 R supports).expect (fun words => turboRemainingFailure22 words d)
  have heq : (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) =
      ∑ T : Finset (Fin 64),
        (rs8ExpectedPacketCount points T.card / (Nat.choose 64 T.card : ℝ)) *
          F (sparseEnumeratedSupport22 S e (fun _ => T)) := by
    rw [turboSetup22_activeSet_expect_eq_support_measure]
    unfold sparseSupportProfileCount22
    rw [sparseSupportProfileCount22_reindex points S e F]
    rw [← Equiv.sum_comp sparseOneSupportEquiv22]
    simp only [sparseOneSupportEquiv22, Equiv.coe_fn_mk, Fin.prod_univ_one]
  rw [heq]
  calc
    _ ≤ ∑ T : Finset (Fin 64),
        (rs8ExpectedPacketCount points T.card / (Nat.choose 64 T.card : ℝ)) *
        (if useChernoff T.card then
          (∑ k : Fin 11, (List.ofFn (fun region : Fin 64 =>
            sparseOneRegionChoice22 (B T.card) R (decide (region ∈ T)))).prod 0 k) /
              (z T.card) ^ d else 1) := by
      apply Finset.sum_le_sum
      intro T _
      apply mul_le_mul_of_nonneg_left
      · split_ifs
        · apply (sparseSupportFailure22_le_kernel R
            (sparseEnumeratedSupport22 S e (fun _ => T)) d
            (z T.card) (hz0 T.card) (hz1 T.card)).trans
          apply div_le_div_of_nonneg_right _ (pow_nonneg (hz0 T.card).le d)
          have hm := sparseExactRegional_product_mass_le R (z T.card)
            (fun region => sparseRegionalSupport22
              (sparseEnumeratedSupport22 S e (fun _ => T)) region) (B T.card) (hB T.card)
          simpa only [sparsePhysicalRegionMatrix22_eq_occupancy,
            sparseEnumeratedSupport22_one_DP] using hm
        · exact sparseSupportFailure22_le_one R _ d
      · apply div_nonneg
        · exact Finset.sum_nonneg fun _ _ => FinPMF.prob_nonneg _ _
        · positivity
    _ = _ := by
      -- Explicit instantiation: `?w #T` is not a higher-order pattern, so a
      -- bare `rw` cannot infer the shell weights.
      refine (sparseOneWeightedSupportFold22_choice_mass (0 : Fin 11)
        (fun v => rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ))
        (fun v => (z v) ^ d) useChernoff
        (fun v => sparseOneRegionChoice22 (B v) R)).trans ?_
      unfold sparseOneChosenShellBound22
      apply Finset.sum_congr rfl
      intro v hv
      have hv64 : v ≤ 64 := by have := Finset.mem_range.mp hv; omega
      have hc : (Nat.choose 64 v : ℝ) ≠ 0 :=
        ne_of_gt (by exact_mod_cast Nat.choose_pos hv64)
      by_cases hu : useChernoff v = true
      · rw [if_pos hu, if_pos hu]
      · rw [if_neg hu, if_neg hu]
        exact div_mul_cancel₀ _ hc

theorem turboSetup22_occupation_one_le_chosen_shells
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (z : ℕ → ℝ) (hz0 : ∀ v, 0 < z v) (hz1 : ∀ v, z v ≤ 1)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ v j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j (z v))) (B v j))
    (useChernoff : ℕ → Bool) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) 1 ω : ℝ)) ≤
          ((R * 16 : ℕ) : ℝ) * sparseOneChosenShellBound22 points B R d z useChernoff := by
  simpa only [Nat.choose_one_right] using
    turboSetup22_occupation_le_of_activeSet_bound R points d 1 (by decide)
      (sparseOneChosenShellBound22 points B R d z useChernoff)
      (fun S hS => turboSetup22_one_activeSet_le_chosen_shells points S
        (S.orderIsoOfFin hS).toEquiv d z hz0 hz1 B hB useChernoff)

#print axioms turboSetup22_one_activeSet_le_chosen_shells
#print axioms turboSetup22_occupation_one_le_chosen_shells
end Spin.Turbo
