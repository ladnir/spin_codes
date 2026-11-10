import SpinCodes.Turbo.Occupation22

/-! The realized code has a nonzero message of weight at most `d` only if
an occupation-1 or occupation-2 message is bad, or a message of occupation
3..16R is bad. Markov bounds the first two by their exact first moments; the
third is the event handled by the occupation-union tail theorem. -/
noncomputable section
namespace Spin.Turbo
open Finset

theorem turboSetup22_bad_prob_split (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.prob (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      (turboSetup22 R points).joint.expect
          (fun ω => ((turboSetup22 R points).ZQ d (fun _ msg => turboOccupation22 msg) 1 ω : ℝ)) +
        (turboSetup22 R points).joint.expect
          (fun ω => ((turboSetup22 R points).ZQ d (fun _ msg => turboOccupation22 msg) 2 ω : ℝ)) +
        (turboSetup22 R points).joint.prob (fun ω =>
          ∃ msg : TurboMessage22 R, turboOccupation22 msg ∈ Finset.Ico 3 (R * 16 + 1) ∧
            (turboSetup22 R points).innerWt ω.2 ((turboSetup22 R points).outer ω.1 msg) ≤ d) := by
  set S := turboSetup22 R points
  set occ : TurboOuterSeed22 R → TurboMessage22 R → ℕ := fun _ msg => turboOccupation22 msg
  set T : TurboOuterSeed22 R × RouteInnerSeed22 R → Prop := fun ω =>
    ∃ msg : TurboMessage22 R, turboOccupation22 msg ∈ Finset.Ico 3 (R * 16 + 1) ∧
      S.innerWt ω.2 (S.outer ω.1 msg) ≤ d
  have hsub : ∀ ω, 1 ≤ S.Z d ω →
      (1 ≤ S.ZQ d occ 1 ω ∨ 1 ≤ S.ZQ d occ 2 ω) ∨ T ω := by
    intro ω h
    unfold Spin.Setup.Z at h
    obtain ⟨msg, hmsg⟩ := Finset.card_pos.mp (Nat.lt_of_lt_of_le Nat.zero_lt_one h)
    obtain ⟨hnz, hw⟩ := Finset.mem_filter.mp hmsg
    have hne : msg ≠ 0 := (Finset.mem_filter.mp hnz).2
    have hocc := Finset.mem_Ico.mp (turboOccupation22_mem_Ico msg hne)
    have hq : ∀ q, turboOccupation22 msg = q → 1 ≤ S.ZQ d occ q ω := by
      intro q hq
      unfold Spin.Setup.ZQ
      exact Finset.card_pos.mpr ⟨msg, Finset.mem_filter.mpr ⟨hnz, hq, hw⟩⟩
    rcases (by omega : turboOccupation22 msg = 1 ∨ turboOccupation22 msg = 2 ∨
        3 ≤ turboOccupation22 msg) with h1 | h2 | h3
    · exact Or.inl (Or.inl (hq 1 h1))
    · exact Or.inl (Or.inr (hq 2 h2))
    · exact Or.inr ⟨msg, Finset.mem_Ico.mpr ⟨h3, hocc.2⟩, hw⟩
  calc S.joint.prob (fun ω => 1 ≤ S.Z d ω)
      ≤ S.joint.prob (fun ω => (1 ≤ S.ZQ d occ 1 ω ∨ 1 ≤ S.ZQ d occ 2 ω) ∨ T ω) :=
        FinPMF.prob_mono _ hsub
    _ ≤ S.joint.prob (fun ω => 1 ≤ S.ZQ d occ 1 ω ∨ 1 ≤ S.ZQ d occ 2 ω) + S.joint.prob T :=
        FinPMF.prob_or_le _ _ _
    _ ≤ (S.joint.prob (fun ω => 1 ≤ S.ZQ d occ 1 ω) + S.joint.prob (fun ω => 1 ≤ S.ZQ d occ 2 ω)) +
          S.joint.prob T := by
        gcongr
        exact FinPMF.prob_or_le _ _ _
    _ ≤ _ := by
        gcongr
        · exact FinPMF.markov_one _ _
        · exact FinPMF.markov_one _ _

#print axioms turboSetup22_bad_prob_split
end Spin.Turbo
