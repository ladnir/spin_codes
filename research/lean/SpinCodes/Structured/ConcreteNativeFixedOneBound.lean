import SpinCodes.Structured.ConcreteFixedRoutedOne
import SpinCodes.Structured.ConcreteNativeDenseIntegration

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset ConcreteMarked

theorem occupation_failure_le_fair {L k T : ℕ} (seed : Seed k) (Q : ℕ)
    (e : (Fin (k*24) × Fin L) ≃ (Fin T × Fin 128)) (d : ℕ) {B P : ℝ}
    (hB : ∀ w≤k*24, (spectrum seed w:ℝ)≤B*((k*24).choose w:ℝ))
    (hP : ∀ S : Finset (Fin L), S.card=Q →
      (fairExperimentLaw S (k*24) T).prob (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2≤d)≤P) :
    (∑ x ∈ univ.filter (fun x => occupation x=Q), failureProbability e (rowSupports seed x) d) ≤
      (L.choose Q:ℝ)*(((2:ℝ)^(k*24)*B)^Q*P) := by
  have hB0 : 0≤B := by
    have hh := hB 0 (by omega)
    simp only [Nat.choose_zero_right,Nat.cast_one,mul_one] at hh
    exact (Nat.cast_nonneg _).trans hh
  rw [sum_occupation_eq_activeSets]
  calc
    _ ≤ ∑ _S ∈ (univ:Finset (Fin L)).powersetCard Q, ((2:ℝ)^(k*24)*B)^Q*P := by
      apply sum_le_sum
      intro S hS
      have hc := (mem_powersetCard.mp hS).2
      have hh := (activeMessages_failure_le_fair seed S e d hB).trans
        (mul_le_mul_of_nonneg_left (hP S hc) (by positivity))
      simpa only [hc] using hh
    _ = _ := by simp [card_powersetCard]

end Spin.Structured.ConcreteOuter
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteRoute ConcreteMarked Placement

theorem tupleWiring_eq_regionMajor (m : ℕ) :
    tupleWiring m = regionMajor (L:=Lsched m) (b:=nativeBlocks m*24) (R:=rounds m)
      (by rw [native_width]; simpa only [Nat.mul_comm] using (rounds_length m).symm) := by
  have key {b b' L T : ℕ} (hb : b=b') (h : b*L=T*128) (h' : b'*L=T*128) :
      ((Equiv.prodCongr (finCongr hb) (Equiv.refl (Fin L))).trans (regionMajor h')) = regionMajor h := by
    subst b'
    rfl
  exact key (native_width m) _ _

def fixedOneBound (m : ℕ) : ℝ :=
  (Lsched m:ℝ)*Real.exp ((bsched m:ℝ)*(Real.log 2/2+1281/100000+
    selectedRemainder m (nativeShellRemainder m)))*
    (1000*(1011/2000:ℝ)^(bsched m))/(Real.exp (-(2/(Lsched m:ℝ))))^(threshold m)

/-- Actual selected native one-row message layer, with no kernel inequality left as a premise. -/
theorem native_one_tuple_bound :
    ∀ᶠ m in atTop, ∀ seed, nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x=1),
        failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤ fixedOneBound m := by
  have hlim : Tendsto (fun m : ℕ => m+1) atTop atTop := tendsto_add_atTop_nat 1
  filter_upwards [hlim.eventually actual_one_fair_probability_regionMajor] with m hm
  intro seed hg
  have hB := native_good_spectrum_bound m seed hg
    (B:=spectrumRatioBound (bsched m) (nativeShellRemainder m)) (by unfold spectrumRatioBound; positivity)
    (native_expected_spectrum_uniform m)
  have hp : ∀ S : Finset (Fin (Lsched m)), S.card=1 →
      (fairExperimentLaw S (nativeBlocks m*24) (rounds m)).prob
        (fun ω => ConcreteRoutedEncoder.weight (tupleWiring m) ω.1 ω.2≤threshold m) ≤
          (1000*(1011/2000:ℝ)^(bsched m))/(Real.exp (-(2/(Lsched m:ℝ))))^(threshold m) := by
    intro S hS
    rw [tupleWiring_eq_regionMajor]
    have hh := hm (nativeBlocks m*24) (rounds m)
      (by change (nativeBlocks m*24)*Lsched m=rounds m*128
          rw [native_width]; simpa only [Nat.mul_comm] using (rounds_length m).symm) S hS (threshold m)
    simp only [native_width,Lsched,Nat.cast_mul,Nat.cast_ofNat,Nat.cast_add,Nat.cast_one] at hh ⊢
    convert hh using 1
    apply FinPMF.prob_congr
    intro ω
    rfl
  have hh := occupation_failure_le_fair seed 1 (tupleWiring m) (threshold m)
    (by simpa only [native_width] using hB) hp
  simpa only [native_width,Nat.choose_one_right,pow_one,selected_row_cost,fixedOneBound,mul_div_assoc,mul_assoc] using hh

theorem native_one_EZ_bound : ∀ᶠ m in atTop, concreteFamily.EZ m 1 ≤ fixedOneBound m := by
  filter_upwards [native_good_positive_eventually,native_one_tuple_bound] with m hm hb
  exact EZ_le_tuple_layer_bound m 1 hm hb

#print axioms tupleWiring_eq_regionMajor
#print axioms native_one_tuple_bound
#print axioms native_one_EZ_bound
end Spin.Structured.ConcreteNativeFamily
