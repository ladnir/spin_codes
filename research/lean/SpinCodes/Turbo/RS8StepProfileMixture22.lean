import SpinCodes.Turbo.RS8StepCounts22

/-!
# Hypergeometric count-profile mixture for a Turbo region

The regional selected set is uniform at fixed total size. Conditional on
the vector of 16-slot counts, the step sets are independent uniform
subsets. The profile weight is the product of local binomial counts
divided by the global binomial count. We first prove the pointwise
mixture; its normalization follows from the two probability laws.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalStepProfileWeight22 (R : ℕ)
    (S : Finset (Fin (R * 16)))
    (counts : Fin R → Fin 17) : ℝ :=
  if (∑ step : Fin R, (counts step).val) = S.card then
    (1 / (Nat.choose (R * 16) S.card : ℝ)) *
      ∏ step : Fin R, (Nat.choose 16 (counts step).val : ℝ)
  else 0

theorem regionalStepProfileWeight22_nonneg (R : ℕ)
    (S : Finset (Fin (R * 16)))
    (counts : Fin R → Fin 17) :
    0 ≤ regionalStepProfileWeight22 R S counts := by
  unfold regionalStepProfileWeight22
  split_ifs <;> positivity

theorem regionalStepChoose_inverse_product (R : ℕ)
    (counts : Fin R → Fin 17) :
    (∏ step : Fin R, (Nat.choose 16 (counts step).val : ℝ)) *
      (∏ step : Fin R,
        1 / (Nat.choose 16 (counts step).val : ℝ)) = 1 := by
  rw [← Finset.prod_mul_distrib]
  have hterm (step : Fin R) :
      (Nat.choose 16 (counts step).val : ℝ) *
        (1 / (Nat.choose 16 (counts step).val : ℝ)) = 1 := by
    have hle : (counts step).val ≤ 16 := by omega
    have hne : (Nat.choose 16 (counts step).val : ℝ) ≠ 0 := by
      exact_mod_cast Nat.choose_ne_zero hle
    field_simp
  simp only [hterm, Finset.prod_const_one]

/-- Exact pointwise disintegration of a uniform `|S|`-subset of the
regional positions by its ordered vector of physical-step counts. -/
theorem regionalSelectedStepSetsLaw22_mixture (R : ℕ)
    (S : Finset (Fin (R * 16)))
    (sets : Fin R → Finset (Fin 16)) :
    (regionalSelectedStepSetsLaw22 R S).p sets =
      ∑ counts : Fin R → Fin 17,
        regionalStepProfileWeight22 R S counts *
          (uniformStepSelectedSets22 counts).p sets := by
  let actual := regionalStepCounts22 sets
  symm
  rw [Finset.sum_eq_single actual]
  · have hmatch : ∀ step : Fin R,
        (sets step).card = (actual step).val := by
      intro step
      rfl
    rw [uniformStepSelectedSets22_apply_profile, if_pos hmatch,
      regionalSelectedStepSetsLaw22_apply]
    by_cases ht : (∑ step : Fin R, (sets step).card) = S.card
    · rw [if_pos ht]
      have htActual : (∑ step : Fin R, (actual step).val) = S.card := by
        simpa only [actual, regionalStepCounts22] using ht
      rw [regionalStepProfileWeight22, if_pos htActual]
      rw [mul_assoc, regionalStepChoose_inverse_product]
      ring
    · rw [if_neg ht]
      have htActual : (∑ step : Fin R, (actual step).val) ≠ S.card := by
        simpa only [actual, regionalStepCounts22] using ht
      rw [regionalStepProfileWeight22, if_neg htActual]
      simp
  · intro counts _ hne
    have hnot : ¬∀ step : Fin R,
        (sets step).card = (counts step).val := by
      intro h
      apply hne
      funext step
      apply Fin.ext
      exact (h step).symm
    rw [uniformStepSelectedSets22_apply_profile, if_neg hnot]
    simp
  · intro h
    exact (h (Finset.mem_univ actual)).elim

theorem regionalStepProfileWeight22_total (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (∑ counts : Fin R → Fin 17,
      regionalStepProfileWeight22 R S counts) = 1 := by
  let P := regionalSelectedStepSetsLaw22 R S
  calc
    (∑ counts : Fin R → Fin 17,
        regionalStepProfileWeight22 R S counts) =
        ∑ counts : Fin R → Fin 17,
          ∑ sets : Fin R → Finset (Fin 16),
            regionalStepProfileWeight22 R S counts *
              (uniformStepSelectedSets22 counts).p sets := by
      apply Finset.sum_congr rfl
      intro counts _
      rw [← Finset.mul_sum,
        (uniformStepSelectedSets22 counts).total, mul_one]
    _ = ∑ sets : Fin R → Finset (Fin 16),
          ∑ counts : Fin R → Fin 17,
            regionalStepProfileWeight22 R S counts *
              (uniformStepSelectedSets22 counts).p sets :=
      Finset.sum_comm
    _ = ∑ sets : Fin R → Finset (Fin 16), P.p sets := by
      apply Finset.sum_congr rfl
      intro sets _
      exact (regionalSelectedStepSetsLaw22_mixture R S sets).symm
    _ = 1 := P.total

def regionalStepProfileLaw22 (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin R → Fin 17) where
  p := regionalStepProfileWeight22 R S
  nonneg := regionalStepProfileWeight22_nonneg R S
  total := regionalStepProfileWeight22_total R S

/-- The actual uniform regional selected-set law disintegrates by its
ordered count profile. Conditional on that profile, all physical-step
selected subsets are independent and uniform at their specified counts. -/
theorem regionalSelectedStepSetsLaw22_eq_profile_mixture (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    regionalSelectedStepSetsLaw22 R S =
      (regionalStepProfileLaw22 R S).bind
        (fun counts => uniformStepSelectedSets22 counts) := by
  apply FinPMF.ext
  funext sets
  rw [FinPMF.bind_p]
  exact regionalSelectedStepSetsLaw22_mixture R S sets

#print axioms regionalSelectedStepSetsLaw22_mixture
#print axioms regionalStepProfileWeight22_total
#print axioms regionalSelectedStepSetsLaw22_eq_profile_mixture

end Spin.Turbo
