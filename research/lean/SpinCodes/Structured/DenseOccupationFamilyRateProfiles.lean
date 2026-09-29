import SpinCodes.Structured.DenseOccupationFamilyRateDefs

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter ConcreteRoute ConcreteRoutedEncoder ConcreteNativeFamily

theorem profileMessages_congr {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    {w w' : Fin L → ℕ} (h : ∀ i ∈ S, w i = w' i) :
    profileMessages seed S w = profileMessages seed S w' := by
  classical
  ext x
  simp only [profileMessages, Finset.mem_filter, Finset.mem_univ, true_and]
  apply forall_congr'
  intro i
  by_cases hi : i ∈ S
  · simp only [hi, ite_true, h i hi]
  · simp only [hi, ite_false]

theorem PointRate.scaled {α x m c η C : ℝ} (h : PointRate α x m c η C)
    {L b R d : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    (hden : density rows = α*x) (hround : 128*R = L*b)
    (hd : (d:ℝ) ≤ (11/100)*((L:ℝ)*b)) {K D : ℝ} (hD : 0 < D) (hK0 : 0 ≤ K)
    (hK : K ≤ D*Real.exp (((L:ℝ)*b)*α*(m*x+c))) :
    K*(experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      D*((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*C*Real.exp (-η*((L:ℝ)*b))) := by
  have hc : K/D ≤ Real.exp (((L:ℝ)*b)*α*(m*x+c)) :=
    (div_le_iff₀ hD).mpr (by simpa only [mul_comm] using hK)
  have hh := h hL hb e rows hw0 hw1 hden hround hd (K/D) (div_nonneg hK0 hD.le) hc
  have hh' := mul_le_mul_of_nonneg_left hh hD.le
  simpa only [← mul_assoc, mul_div_cancel₀ _ hD.ne'] using hh'

/-- Full selected-seed bound for every actual fixed weight profile, including empty profiles. -/
theorem DenseRates.profile_bound {η C : ℝ} (h : DenseRates η C) (hC : 0 ≤ C)
    {L k R d : ℕ} (hL : 0 < L) (hk : 0 < k)
    (seed : Seed k) (hg : constituentGood k seed) (S : Finset (Fin L))
    (hS : S.Nonempty) (hα : (S.card:ℝ)/L ∈ Set.Icc (1/10000) 1)
    (w : WeightProfile S (k*24))
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128))
    (hround : 128*R = L*(k*24)) (hd : (d:ℝ) ≤ (11/100)*((L:ℝ)*(k*24:ℕ))) :
    (∑ x ∈ profileMessages seed S (profileWeights S w),
      failureProbability e (rowSupports seed x) d) ≤
      outerCost (k*24)^S.card *
        (((((k*24:ℕ):ℝ)+1)^S.card*((L:ℝ)+1)^(k*24))*C*
          Real.exp (-η*((L:ℝ)*(k*24:ℕ)))) := by
  classical
  by_cases hempty : (profileMessages seed S (profileWeights S w)).Nonempty
  · obtain ⟨message,hm⟩ := hempty
    have hact : message ∈ activeMessages S := by
      have hf : message ∈ (activeMessages S).filter (fun x => messageProfile seed S x = w) := by
        rw [profile_fiber_eq]
        exact hm
      exact (Finset.mem_filter.mp hf).1
    have hs := selected_active_rows hk seed hg S hact
    obtain ⟨x,hx,hmean,hden,hactive,hw0,hw1⟩ :=
      selected_message_parameters hL hk seed hg S hS hact
    obtain ⟨p,hp,hrate⟩ := h ((S.card:ℝ)/L) x hα hx
    have hweights : ∀ i ∈ S,
        (rowSupports seed message i).card = profileWeights S w i := by
      intro i hi
      have hh := (Finset.mem_filter.mp hm).2 i
      simp only [hi,ite_true] at hh
      simpa only [rowSupports,support_card] using hh.2
    rw [← profileMessages_congr seed S hweights, profile_failure_eq seed S _ hs.1 e d]
    have hpcount := selected_profile_support hk seed (weightWindow (k*24)) hg S
      (fun i => (rowSupports seed message i).card) hs.2.1 hs.2.2 hmean hp
    have hLR : (L:ℝ) ≠ 0 := by exact_mod_cast hL.ne'
    have he : ((k*24:ℕ):ℝ)*(S.card:ℝ)*((p.1:ℝ)*x+(p.2:ℝ)) =
        ((L:ℝ)*(k*24:ℕ))*((S.card:ℝ)/L)*((p.1:ℝ)*x+(p.2:ℝ)) := by field_simp
    rw [he] at hpcount
    have hh := hrate.scaled hL (by omega) e (rowSupports seed message) hw0 hw1 hden hround hd
      (pow_pos (outerCost_pos (by omega)) _) (Finset.prod_nonneg (fun _ _ => Nat.cast_nonneg _)) hpcount
    simpa only [hactive,failureProbability] using hh
  · rw [Finset.not_nonempty_iff_eq_empty.mp hempty, Finset.sum_empty]
    unfold outerCost
    positivity

#print axioms profileMessages_congr
#print axioms PointRate.scaled
#print axioms DenseRates.profile_bound
end Spin.Structured.DenseOccupationFixed
