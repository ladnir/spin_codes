import SpinCodes.Structured.DenseOccupationRowParameters
import SpinCodes.Structured.DenseOccupationProfilePartition
import SpinCodes.Structured.ConcreteNativeSelection

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter ConcreteRoute ConcreteNativeFamily

theorem weightWindow_relative {b w : ℕ} (hb : 0 < b) (hw : w ∈ weightWindow b) :
    (w:ℝ)/b ∈ Set.Icc (13/125) (112/125) := by
  have hbR : (0:ℝ) < b := by exact_mod_cast hb
  have h := (Finset.mem_filter.mp hw).2
  have hl : (13:ℝ)*b ≤ 125*w := by exact_mod_cast h.1
  have hu : (125:ℝ)*w ≤ 112*b := by exact_mod_cast h.2
  constructor
  · apply (le_div_iff₀ hbR).mpr
    linarith
  · apply (div_le_iff₀ hbR).mpr
    linarith

/-- The selected common seed supplies the row window for every actual active message tuple. -/
theorem selected_active_rows {L k : ℕ} (hk : 0 < k) (seed : Seed k)
    (hg : constituentGood k seed) (S : Finset (Fin L))
    {messages : Fin L → LocalMessage k} (hm : messages ∈ activeMessages S) :
    (∀ i, i ∉ S → rowSupports seed messages i = ∅) ∧
    (∀ i ∈ S, (rowSupports seed messages i).card ∈ weightWindow (k*24)) ∧
    (∀ i ∈ S, ((rowSupports seed messages i).card:ℝ)/(k*24:ℕ) ∈
      Set.Icc (13/125) (112/125)) := by
  have ha := (Finset.mem_filter.mp hm).2
  have hz : ∀ i, i ∉ S → rowSupports seed messages i = ∅ := by
    intro i hi
    apply (rowSupports_empty_iff seed messages i).mpr
    exact not_not.mp (fun hn => hi ((ha i).mp hn))
  have hw : ∀ i ∈ S, (rowSupports seed messages i).card ∈ weightWindow (k*24) := by
    intro i hi
    simpa only [rowSupports,support_card] using good_weight_mem hg ((ha i).mpr hi)
  exact ⟨hz,hw,fun i hi => weightWindow_relative (by omega) (hw i hi)⟩

/-- All dense-box parameters follow from the actual selected experiment and its active positions. -/
theorem selected_message_parameters {L k : ℕ} (hL : 0 < L) (hk : 0 < k)
    (seed : Seed k) (hg : constituentGood k seed) (S : Finset (Fin L)) (hS : S.Nonempty)
    {messages : Fin L → LocalMessage k} (hm : messages ∈ activeMessages S) :
    ∃ x : ℝ, x ∈ Set.Icc (13/125) (112/125) ∧
      (∑ i ∈ S, ((rowSupports seed messages i).card:ℝ)/(k*24:ℕ)) = (S.card:ℝ)*x ∧
      density (rowSupports seed messages) = ((S.card:ℝ)/L)*x ∧
      activeRows (rowSupports seed messages) = S.card ∧
      0 < totalWeight (rowSupports seed messages) ∧
      totalWeight (rowSupports seed messages) < L*(k*24) := by
  have h := selected_active_rows hk seed hg S hm
  exact row_parameters hL (by omega) S hS _ h.1 h.2.2

#print axioms weightWindow_relative
#print axioms selected_active_rows
#print axioms selected_message_parameters
end Spin.Structured.DenseOccupationFixed
