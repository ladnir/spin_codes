import SpinCodes.UniformMixer.ScalableSetup

/-!
# One sampled state map has the production transition kernel

The setup samples a uniform invertible linear map at each step. For any
fixed input and entering state, its next-state law is exactly the
`update` kernel used by the regional analysis.
-/

noncomputable section

namespace Spin.UniformMixer

attribute [local instance] Classical.propDecidable

theorem uniform_state_map_update (c q r : State22) :
    (FinPMF.uniform StateMap).prob
      (fun g => g q + c = r) = update c q r := by
  change (FinPMF.uniform StateMap).prob
    (fun g : StateMap => g • q + c = r) = update c q r
  rw [FinPMF.uniform_prob]
  have hcard : ((Finset.univ : Finset StateMap).filter
      (fun g => g • q + c = r)).card =
      Fintype.card (updateEvent (G := StateMap) c q r) := by
    calc
      _ = @Fintype.card {g : StateMap // g ∈
          (Finset.univ : Finset StateMap).filter (fun g => g • q + c = r)}
          (Finset.Subtype.fintype _) :=
        (Fintype.card_coe _).symm
      _ = Fintype.card (updateEvent (G := StateMap) c q r) := by
        exact @Fintype.card_congr _ _ (Finset.Subtype.fintype _) (inferInstance)
          { toFun := fun g => ⟨g.1, (Finset.mem_filter.mp g.2).2⟩
            invFun := fun g => ⟨g.1,
              Finset.mem_filter.mpr ⟨Finset.mem_univ _, g.2⟩⟩
            left_inv := by intro g; apply Subtype.ext; rfl
            right_inv := by intro g; apply Subtype.ext; rfl }
  rw [hcard]
  exact (groupUpdate_eq_event_count (G := StateMap) c q r).symm.trans
    (state22_uniform_gl_kernel c q r)

/-- An arbitrary observable of the next state has the checked update law. -/
theorem uniform_state_map_expect (c q : State22) (f : State22 → ℝ) :
    (FinPMF.uniform StateMap).expect (fun g => f (g q + c)) =
      ∑ r : State22, update c q r * f r := by
  rw [FinPMF.expect_comp (FinPMF.uniform StateMap)
    (fun g : StateMap => g q + c) Finset.univ
    (by intro _; exact Finset.mem_univ _) f]
  simp only [uniform_state_map_update]

end Spin.UniformMixer
