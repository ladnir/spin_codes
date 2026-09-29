import SpinCodes.Prob

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.FinPMF
open Filter
variable {Ω : Type*} [Fintype Ω]

/-- A total conditioning event: use the requested event whenever it has positive mass. -/
def positiveSelection (P : FinPMF Ω) (G : Ω → Prop) : Ω → Prop :=
  if 0 < P.prob G then G else fun _ => True

theorem positiveSelection_eq (P : FinPMF Ω) (G : Ω → Prop) (hG : 0 < P.prob G) :
    P.positiveSelection G = G := by simp only [positiveSelection, if_pos hG]

theorem positiveSelection_pos (P : FinPMF Ω) (G : Ω → Prop) :
    0 < P.prob (P.positiveSelection G) := by
  by_cases hG : 0 < P.prob G
  · simpa only [positiveSelection_eq P G hG] using hG
  · rw [positiveSelection, if_neg hG]
    simp only [prob, Finset.filter_true, P.total]
    norm_num

theorem positiveSelection_compl_le (P : FinPMF Ω) (G : Ω → Prop) :
    P.prob (fun ω => ¬ P.positiveSelection G ω) ≤ P.prob (fun ω => ¬ G ω) := by
  apply P.prob_mono
  intro ω h
  by_cases hG : 0 < P.prob G
  · simpa only [positiveSelection_eq P G hG] using h
  · simp [positiveSelection, hG] at h

theorem positiveSelection_eq_eventually {Ωm : ℕ → Type*} [∀ m, Fintype (Ωm m)]
    (P : ∀ m, FinPMF (Ωm m)) (G : ∀ m, Ωm m → Prop)
    (h : Tendsto (fun m => (P m).prob (fun ω => ¬ G m ω)) atTop (nhds 0)) :
    ∀ᶠ m in atTop, (P m).positiveSelection (G m) = G m := by
  have he := h.eventually (gt_mem_nhds (show (0 : ℝ) < 1 by norm_num))
  filter_upwards [he] with m hm
  apply positiveSelection_eq
  rw [prob_compl] at hm
  linarith

theorem positiveSelection_compl_tendsto {Ωm : ℕ → Type*} [∀ m, Fintype (Ωm m)]
    (P : ∀ m, FinPMF (Ωm m)) (G : ∀ m, Ωm m → Prop)
    (h : Tendsto (fun m => (P m).prob (fun ω => ¬ G m ω)) atTop (nhds 0)) :
    Tendsto (fun m => (P m).prob (fun ω => ¬ (P m).positiveSelection (G m) ω))
      atTop (nhds 0) := by
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds h
    (fun _ => prob_nonneg _ _) (fun m => positiveSelection_compl_le (P m) (G m))

end Spin.FinPMF

