import SpinCodes.Structured.ConcreteOuterEnvelope
import SpinCodes.Structured.ConcreteOuterEnvelopeShiftTransition

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset

lemma golay_weight_even (m : Fin 12 → Bool) : Even (Golay.W m) := by
  have hp : 0 < Golay.A (Golay.W m) := by
    rw [← Golay.card_weight]
    exact card_pos.mpr ⟨m, by simp⟩
  unfold Golay.A at hp
  split_ifs at hp with h₀ h₈ h₁₂ h₁₆ h₂₄
  all_goals first | omega | (rw [‹Golay.W m = _›]; decide)

lemma golayCoefficient_even {k a : ℕ} (ha : golayCoefficient k a ≠ 0) : Even a := by
  have he : (univ.filter (fun f : Fin k → (Fin 12 → Bool) => ∑ i, Golay.W (f i) = a)).card =
      golayCoefficient k a := card_tuples_weight Golay.W k a
  have hp : (univ.filter (fun f : Fin k → (Fin 12 → Bool) => ∑ i, Golay.W (f i) = a)).Nonempty :=
    card_pos.mp (by omega)
  obtain ⟨f,hf⟩ := hp
  have hf' := (mem_filter.mp hf).2
  rw [← hf']
  rw [even_iff_two_dvd]
  exact Finset.dvd_sum (fun i _ => (even_iff_two_dvd.mp (golay_weight_even (f i))))

/-- Clip only the final real count; at most half a bit is changed. -/
def clippedOutput (b a c : ℕ) : ℝ := min (c:ℝ) ((b:ℝ)-(a:ℝ)/2)

lemma clippedOutput_spec {b a c : ℕ} (ha : 0 < a) (hab : a ≤ b)
    (ht : accT b a c ≠ 0) :
    (a:ℝ)/2 ≤ clippedOutput b a c ∧
    clippedOutput b a c ≤ (b:ℝ)-(a:ℝ)/2 ∧
    clippedOutput b a c ≤ c ∧ (c:ℝ)-clippedOutput b a c ≤ 1/2 := by
  obtain ⟨hc,hcb,hr,hs⟩ := accT_supported (by omega) ht
  have hlo : a ≤ 2*c := by omega
  have hhi : 2*c ≤ 2*b-a+1 := by omega
  have hloR : (a:ℝ) ≤ 2*(c:ℝ) := by exact_mod_cast hlo
  have hhiR : 2*(c:ℝ) ≤ 2*(b:ℝ)-(a:ℝ)+1 := by
    have : 2*c+a ≤ 2*b+1 := by omega
    have hh : 2*(c:ℝ)+(a:ℝ) ≤ 2*(b:ℝ)+1 := by exact_mod_cast this
    linarith
  have habR : (a:ℝ) ≤ b := by exact_mod_cast hab
  unfold clippedOutput
  refine ⟨le_min (by linarith) (by linarith), min_le_right _ _, min_le_left _ _, ?_⟩
  rcases le_total (c:ℝ) ((b:ℝ)-(a:ℝ)/2) with h | h
  · rw [min_eq_left h]; norm_num
  · rw [min_eq_right h]; linarith

lemma accT_even_feasible {b a c : ℕ} (ha : 0 < a) (hpar : Even a)
    (ht : accT b a c ≠ 0) : (a:ℝ)/2 ≤ c ∧ (c:ℝ) ≤ (b:ℝ)-(a:ℝ)/2 := by
  obtain ⟨hc,hcb,hr,hs⟩ := accT_supported (by omega) ht
  obtain ⟨j,hj⟩ := hpar
  have h₁ : a ≤ 2*c := by omega
  have h₂ : 2*c+a ≤ 2*b := by omega
  have h₁R : (a:ℝ) ≤ 2*(c:ℝ) := by exact_mod_cast h₁
  have h₂R : 2*(c:ℝ)+(a:ℝ) ≤ 2*(b:ℝ) := by exact_mod_cast h₂
  constructor <;> linarith

end Spin.Structured.ConcreteOuter


