import SpinCodes.Turbo.SparseRare22

/-! The actual setup pays the rare outer event once, and retains the
remaining count under its original product law. -/
noncomputable section
namespace Spin.Turbo

/-- The rare-complement count is an unconditional random variable. -/
def turboSparseTruncatedCount22 (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (ω : TurboOuterSeed22 R × RouteInnerSeed22 R) : ℕ :=
  if turboSparseRare22 R points ω.1 then 0 else (turboSetup22 R points).Z d ω

/-- No conditioned-law independence is used: the rare setup event is
charged once, followed by the unconditional truncated first moment. -/
theorem turboSetup22_bad_prob_le_rare_truncated
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      (R * 16 : ℕ) * ∑ v ∈ Finset.Icc 5 8, rs8ExpectedPacketCount points v +
        (turboSetup22 R points).joint.expect
          (fun ω => (turboSparseTruncatedCount22 R points d ω : ℝ)) := by
  let P := (turboSetup22 R points).joint
  have hsplit : P.prob (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      P.prob (fun ω => turboSparseRare22 R points ω.1) +
        P.prob (fun ω => 1 ≤ turboSparseTruncatedCount22 R points d ω) := by
    apply le_trans (P.prob_mono
      (s := fun ω => 1 ≤ (turboSetup22 R points).Z d ω)
      (t := fun ω => turboSparseRare22 R points ω.1 ∨
        1 ≤ turboSparseTruncatedCount22 R points d ω) (fun ω h => ?_)) (P.prob_or_le _ _)
    by_cases hr : turboSparseRare22 R points ω.1
    · exact Or.inl hr
    · exact Or.inr (by simpa [turboSparseTruncatedCount22, hr] using h)
  have hrare : P.prob (fun ω => turboSparseRare22 R points ω.1) ≤
      (R * 16 : ℕ) * ∑ v ∈ Finset.Icc 5 8, rs8ExpectedPacketCount points v := by
    change ((turboOuterLaw22 R).prod (routeInnerLaw22 R)).prob
      (fun ω => turboSparseRare22 R points ω.1) ≤ _
    rw [FinPMF.prob_prod_left]
    exact turboSparseRare22_prob_le R points
  exact hsplit.trans (add_le_add hrare (P.markov_one _))

/-- Fubini retains the indicator on the outer seed, before summing the
route/inner failure probabilities under the original law. -/
theorem turboSparseTruncatedCount22_expect_eq
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.expect
      (fun ω => (turboSparseTruncatedCount22 R points d ω : ℝ)) =
      (turboOuterLaw22 R).expect (fun seed =>
        if turboSparseRare22 R points seed then 0 else
          ∑ msg ∈ Spin.nonzeroMsgs (TurboMessage22 R),
            (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) := by
  change ((turboOuterLaw22 R).prod (routeInnerLaw22 R)).expect _ = _
  rw [FinPMF.expect_prod]
  apply congrArg ((turboOuterLaw22 R).expect)
  funext seed
  by_cases hr : turboSparseRare22 R points seed
  · simp only [turboSparseTruncatedCount22, hr, ite_true, ite_false, Nat.cast_zero,
      FinPMF.expect_const]
  · simp only [turboSparseTruncatedCount22, hr, ite_true, ite_false]
    have hcard : (fun inner : RouteInnerSeed22 R =>
        ((turboSetup22 R points).Z d (seed, inner) : ℝ)) =
        fun inner => ∑ msg ∈ Spin.nonzeroMsgs (TurboMessage22 R),
          if (turboSetup22 R points).innerWt inner
            (turboOuterWord22 points seed msg) ≤ d then (1 : ℝ) else 0 := by
      funext inner
      simp only [Spin.Setup.Z, Finset.card_filter, Nat.cast_sum]
      apply Finset.sum_congr rfl
      intro msg _
      split_ifs <;> simp_all [turboSetup22] <;> omega
    rw [hcard, FinPMF.expect_sum]
    apply Finset.sum_congr rfl
    intro msg _
    exact (FinPMF.prob_eq_expect_indicator (routeInnerLaw22 R)
      (fun inner => (turboSetup22 R points).innerWt inner
        (turboOuterWord22 points seed msg) ≤ d)).symm
/-- Heavy means every active group of this message has at least nine
packets under the same sampled outer seed. -/
def turboSparseHeavyMessage22 {R : ℕ} (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R) (msg : TurboMessage22 R) : Prop :=
  ∀ i, msg i ≠ 0 → 9 ≤ rs8SampledPacketWeight points (msg i) (seed i)

noncomputable instance {R : ℕ} (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R) :
    DecidablePred (turboSparseHeavyMessage22 points seed) :=
  fun _ => Classical.propDecidable _

/-- In particular the q=2 contribution can be restricted to heavy words
by discarding the rare-complement indicator, under the original law. -/
theorem turboSparseTruncatedOccupation22_le_heavy
    (R : ℕ) (points : Fin 8 → Field16) (hp : Function.Injective points)
    (d q : ℕ) (seed : TurboOuterSeed22 R) :
    (if turboSparseRare22 R points seed then 0 else
      ∑ msg ∈ (Spin.nonzeroMsgs (TurboMessage22 R)).filter
        (fun msg => turboOccupation22 msg = q),
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ msg ∈ (Spin.nonzeroMsgs (TurboMessage22 R)).filter
        (fun msg => turboOccupation22 msg = q ∧
          turboSparseHeavyMessage22 points seed msg),
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg) := by
  by_cases hr : turboSparseRare22 R points seed
  · simp only [if_pos hr]
    exact Finset.sum_nonneg fun _ _ => (turboSetup22 R points).qd_nonneg _ _
  · rw [if_neg hr]
    apply le_of_eq
    apply Finset.sum_congr
    · ext msg
      simp only [Finset.mem_filter]
      have hheavy : turboSparseHeavyMessage22 points seed msg := by
        intro i hi
        exact turboSparseRare22_compl_weight_ge_nine R points hp seed hr i (msg i) hi
      simp only [hheavy, and_true]
    · intro _ _
      rfl
#print axioms turboSparseTruncatedCount22_expect_eq
#print axioms turboSparseTruncatedOccupation22_le_heavy
#print axioms turboSetup22_bad_prob_le_rare_truncated
end Spin.Turbo







