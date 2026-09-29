import SpinCodes.Structured.DenseOccupationFamilyRateProfiles

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter ConcreteNativeFamily

/-- Complete finite occupation-layer first moment for the actual selected shared outer code. -/
theorem DenseRates.layer_bound {η C : ℝ} (h : DenseRates η C) (hC : 0 ≤ C)
    {L k R d Q : ℕ} (hL : 0 < L) (hk : 0 < k) (hQ : 0 < Q)
    (seed : Seed k) (hg : constituentGood k seed)
    (hα : (Q:ℝ)/L ∈ Set.Icc (1/10000) 1)
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128))
    (hround : 128*R = L*(k*24)) (hd : (d:ℝ) ≤ (11/100)*((L:ℝ)*(k*24:ℕ))) :
    (∑ x ∈ univ.filter (fun x => occupation x = Q),
      failureProbability e (rowSupports seed x) d) ≤
      densePrefactor L (k*24) Q C * Real.exp (-η*((L:ℝ)*(k*24:ℕ))) := by
  let B : ℝ := outerCost (k*24)^Q *
    (((((k*24:ℕ):ℝ)+1)^Q*((L:ℝ)+1)^(k*24))*C*Real.exp (-η*((L:ℝ)*(k*24:ℕ))))
  have hh := sum_occupation_le_profile_bound seed Q
    (fun x => failureProbability e (rowSupports seed x) d) (B := B) (by
      intro S hS w
      have hSn : S.Nonempty := Finset.card_pos.mp (by omega)
      have hαS : (S.card:ℝ)/L ∈ Set.Icc (1/10000) 1 := by simpa only [hS] using hα
      simpa only [hS,B] using h.profile_bound hC hL hk seed hg S hSn hαS w e hround hd)
  convert hh using 1
  dsimp [B,densePrefactor]
  rw [two_mul,pow_add]
  ring

#print axioms DenseRates.layer_bound
end Spin.Structured.DenseOccupationFixed
