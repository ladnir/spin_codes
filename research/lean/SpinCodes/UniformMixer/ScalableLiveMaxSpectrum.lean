import SpinCodes.UniformMixer.ScalableFeedbackMomentExact

/-!
# Live-row maximum over the eight realized nonzero weights

The production image contains every weight in the nonzero spectrum.
Consequently its maximum Bernoulli row moment can be evaluated over eight
weight values instead of over all 2^22 state vectors.
-/

noncomputable section

namespace Spin.UniformMixer

def liveSpectrumWeights : Finset Nat := spectrumWeights.erase 0

theorem liveSpectrumWeights_nonempty : liveSpectrumWeights.Nonempty := by
  refine ⟨16, ?_⟩
  decide

theorem productionWeight_mem_liveSpectrumWeights (q : State22) (hq : q ≠ 0) :
    coordinateWeight coordMap q ∈ liveSpectrumWeights := by
  apply Finset.mem_erase.mpr
  constructor
  · intro hw
    exact hq (productionWeight_zero_only q hw)
  · exact productionWeight_mem_spectrumWeights q

def scalableSpectrumMaxMoment (β z : ℝ) : ℝ :=
  liveSpectrumWeights.sup' liveSpectrumWeights_nonempty
    (productionWeightMoment β z)

theorem scalableMaxLiveMoment_le_spectrum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableMaxLiveMoment β z hβ0 hβ1 ≤
      scalableSpectrumMaxMoment β z := by
  unfold scalableMaxLiveMoment
  apply Finset.sup'_le
  intro q hq
  rw [scalableLiveMoment_weight_formula]
  have hw := productionWeight_mem_liveSpectrumWeights q
    (Finset.ne_of_mem_erase hq)
  change productionWeightMoment β z (coordinateWeight coordMap q) ≤
    liveSpectrumWeights.sup' liveSpectrumWeights_nonempty
      (productionWeightMoment β z)
  exact Finset.le_sup' (productionWeightMoment β z) hw

theorem productionWeightCount_pos_of_mem_live (w : Nat)
    (hw : w ∈ liveSpectrumWeights) : 0 < productionWeightCount w := by
  have hmem : w ∈ spectrumWeights := (Finset.mem_erase.mp hw).2
  have hne : w ≠ 0 := (Finset.mem_erase.mp hw).1
  simp only [spectrumWeights, Finset.mem_insert,
    Finset.mem_singleton] at hmem
  rcases hmem with h | h | h | h | h | h | h | h | h
  · exact (hne h).elim
  all_goals subst w
  all_goals simp [productionWeightCount_16, productionWeightCount_24,
    productionWeightCount_28, productionWeightCount_32,
    productionWeightCount_36, productionWeightCount_40,
    productionWeightCount_48, productionWeightCount_64]

theorem productionWeight_exists_of_mem_live (w : Nat)
    (hw : w ∈ liveSpectrumWeights) :
    ∃ q : State22, q ≠ 0 ∧ coordinateWeight coordMap q = w := by
  have hcard : 0 <
      (Finset.univ.filter (fun q : State22 =>
        coordinateWeight coordMap q = w)).card := by
    rw [← productionWeightCount_eq_card_filter]
    exact productionWeightCount_pos_of_mem_live w hw
  obtain ⟨q, hq⟩ := Finset.card_pos.mp hcard
  have hweight : coordinateWeight coordMap q = w :=
    (Finset.mem_filter.mp hq).2
  refine ⟨q, ?_, hweight⟩
  intro hzero
  subst q
  have hwzero : w ≠ 0 := (Finset.mem_erase.mp hw).1
  exact hwzero (by simpa [coordinateWeight_zero] using hweight.symm)

theorem scalableMaxLiveMoment_eq_spectrum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableMaxLiveMoment β z hβ0 hβ1 =
      scalableSpectrumMaxMoment β z := by
  apply le_antisymm
  · exact scalableMaxLiveMoment_le_spectrum β z hβ0 hβ1
  · unfold scalableSpectrumMaxMoment
    apply Finset.sup'_le
    intro w hw
    obtain ⟨q, hq0, hweight⟩ := productionWeight_exists_of_mem_live w hw
    calc
      productionWeightMoment β z w =
          scalableLiveMoment β z hβ0 hβ1 q := by
        rw [scalableLiveMoment_weight_formula, hweight]
        rfl
      _ ≤ scalableMaxLiveMoment β z hβ0 hβ1 :=
        Finset.le_sup' _ (Finset.mem_erase.mpr
          ⟨hq0, Finset.mem_univ q⟩)

end Spin.UniformMixer
