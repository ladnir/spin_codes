import SpinCodes.Structured.ConcreteMarkedFairRows

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset

lemma fairOnMarks_mass {n : ℕ} (S T : Finset (Fin n)) :
    (fairOnMarks S).p T = if T⊆S then (1/2:ℝ)^S.card else 0 := by
  rw [fairOnMarks_eq_poisson,ConcreteRoute.poissonBinom_eq_prod]
  by_cases hsub : T⊆S
  · rw [if_pos hsub]
    have he : (∏ i:Fin n,if i∈T then activeProbability S i else 1-activeProbability S i)=
        ∏ i∈S,(if i∈T then activeProbability S i else 1-activeProbability S i) := by
      symm
      apply prod_subset (subset_univ S)
      intro i hi hiS
      have hiT : i∉T := fun h => hiS (hsub h)
      simp [activeProbability,hiS,hiT]
    rw [he]
    have he' : (∏ i∈S,(if i∈T then activeProbability S i else 1-activeProbability S i))=
        ∏ _i∈S,(1/2:ℝ) := by
      apply prod_congr rfl
      intro i hi
      by_cases hiT : i∈T <;> norm_num [activeProbability,hi,hiT]
    rw [he']
    simp
  · rw [if_neg hsub]
    obtain ⟨i,hiT,hiS⟩ := not_subset.mp hsub
    exact prod_eq_zero (mem_univ i) (by simp [hiT,activeProbability,hiS])

lemma subset_map_univ_iff {n Q : ℕ} (marks : Fin Q ↪ Fin n) (T : Finset (Fin n)) :
    T⊆univ.map marks ↔ ∃ A:Finset (Fin Q),A.map marks=T := by
  constructor
  · intro h
    refine ⟨univ.filter (fun i => marks i∈T),?_⟩
    ext j
    constructor
    · intro hj
      obtain ⟨i,hi,rfl⟩ := mem_map.mp hj
      exact (mem_filter.mp hi).2
    · intro hj
      obtain ⟨i,hi,rfl⟩ := mem_map.mp (h hj)
      exact mem_map.mpr ⟨i,mem_filter.mpr ⟨mem_univ _,hj⟩,rfl⟩
  · rintro ⟨A,rfl⟩
    exact map_subset_map.mpr (subset_univ _)

/-- Fair bits on Q embedded marked positions are exactly a uniformly sampled
subset of the Q labels, pushed through that embedding. -/
theorem fairOnMarks_uniform_map {n Q : ℕ} (marks : Fin Q ↪ Fin n) :
    fairOnMarks (univ.map marks)=
      (FinPMF.uniform (Finset (Fin Q))).map (fun A => A.map marks) := by
  ext T
  rw [fairOnMarks_mass,FinPMF.map_p]
  by_cases hT : ∃ A:Finset (Fin Q),A.map marks=T
  · obtain ⟨A,rfl⟩ := hT
    rw [if_pos (map_subset_map.mpr (subset_univ _))]
    have he (B:Finset (Fin Q)) : B.map marks=A.map marks ↔ B=A := (map_injective marks).eq_iff
    simp only [he,FinPMF.uniform,card_map,card_univ,Fintype.card_fin]
    simp [Fintype.card_finset,one_div_pow]
  · rw [if_neg ((subset_map_univ_iff marks T).not.mpr hT)]
    symm
    apply sum_eq_zero
    intro A hA
    exact if_neg (fun he => hT ⟨A,he⟩)

/-- The concrete uniform-subset expectation formula, including Q=0. -/
theorem fairOnMarks_expect_uniform_sum {n Q : ℕ} (marks : Fin Q ↪ Fin n)
    (F : Finset (Fin n)→ℝ) :
    (fairOnMarks (univ.map marks)).expect F = (∑ A:Finset (Fin Q),F (A.map marks))/(2:ℝ)^Q := by
  rw [fairOnMarks_uniform_map,FinPMF.expect_map]
  simp only [FinPMF.expect,FinPMF.uniform,Fintype.card_finset,Fintype.card_fin,Nat.cast_pow,Nat.cast_ofNat]
  rw [←mul_sum]
  ring

end Spin.Structured.ConcreteMarked
