import SpinCodes.Structured.DenseOccupationFixedStrip
import SpinCodes.Structured.DenseOccupationOuterAffine

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open ConcreteOuter ConcreteRoute ConcreteRoutedEncoder

/-- The certified strip with its finite outer selection cost retained. -/
theorem selected_profile_routed_rate {L k R d : ℕ} (hL : 0 < L) (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (S : Finset (Fin L)) (rows : Fin L → Finset (Fin (k*24)))
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128))
    (hw : ∀ i ∈ S, (rows i).card ∈ W)
    (hr : ∀ i ∈ S, ((rows i).card:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125))
    {α x : ℝ} (hmean : ∑ i ∈ S, ((rows i).card:ℝ)/(k*24:ℕ) = (S.card:ℝ)*x)
    (hαQ : (L:ℝ)*α = S.card)
    (hα : α ∈ Set.Icc (1/10000) (10031/320000))
    (hx : x ∈ Set.Icc (13/125) (18296026121/120000000000))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*(k*24))
    (hden : density rows = α*x) (hround : 128*R = L*(k*24))
    (hd : (d:ℝ) ≤ (11/100)*((L:ℝ)*(k*24:ℕ))) :
    (∏ i ∈ S, (spectrum seed (rows i).card:ℝ)) *
      (experimentLaw L (k*24) R).prob (fun ω => weight e rows ω ≤ d) ≤
    outerCost (k*24)^S.card *
      (((((k*24:ℕ):ℝ)+1)^activeRows rows*((L:ℝ)+1)^(k*24))*568*
        Real.exp (-(4/10000000)*((L:ℝ)*(k*24:ℕ)))) := by
  let C : ℝ := ∏ i ∈ S, (spectrum seed (rows i).card:ℝ)
  let D : ℝ := outerCost (k*24)^S.card
  have hD : 0 < D := pow_pos (outerCost_pos (by omega)) _
  have hC : 0 ≤ C := Finset.prod_nonneg (fun _ _ => Nat.cast_nonneg _)
  have hb : 0 < k*24 := by omega
  have hp := selected_profile_affine hk seed W hg S (fun i => (rows i).card) hw hr hmean
  have he : ((k*24:ℕ):ℝ)*(S.card:ℝ)*((833/500)*x-43311/250000) =
      ((L:ℝ)*(k*24:ℕ))*α*((833/500)*x-43311/250000) := by
    rw [← hαQ]
    ring
  rw [he] at hp
  have hc : C/D ≤ Real.exp (((L:ℝ)*(k*24:ℕ))*α*((833/500)*x-43311/250000)) :=
    (div_le_iff₀ hD).mpr (by simpa only [mul_comm] using hp)
  have h := Strip.routed_rate hL hb e rows hw0 hw1 hα hx hden hround hd
    (div_nonneg hC hD.le) hc
  have h' := mul_le_mul_of_nonneg_left h hD.le
  simpa only [← mul_assoc, mul_div_cancel₀ _ hD.ne', C, D] using h'

#print axioms selected_profile_routed_rate
end Spin.Structured.DenseOccupationFixed
