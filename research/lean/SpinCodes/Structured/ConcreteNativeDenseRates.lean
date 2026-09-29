import SpinCodes.Structured.ConcreteNativeDenseIntegration
import SpinCodes.Structured.DenseOccupationFamilyRateLayer
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteRoute DenseOccupationFixed

theorem native_dense_density (m Q : ℕ)
    (hQ : Q ∈ Ico (nativeCut m+1) (Lsched m+1)) :
    (Q:ℝ)/(Lsched m:ℝ) ∈ Set.Icc (1/10000) 1 := by
  have hL : (0:ℝ)<Lsched m := by exact_mod_cast Lsched_pos m
  have hn := mem_Ico.mp hQ
  have hlo : Lsched m < 10000*Q := by unfold nativeCut at hn; omega
  have hhi : Q ≤ Lsched m := by omega
  constructor
  · apply (le_div_iff₀ hL).mpr
    have : (Lsched m:ℝ)<10000*(Q:ℝ) := by exact_mod_cast hlo
    linarith
  · apply (div_le_iff₀ hL).mpr
    simpa only [one_mul] using (show (Q:ℝ)≤(Lsched m:ℝ) by exact_mod_cast hhi)

/-- Mixed numerical rates bound every dense layer of the actual native code. -/
theorem dense_layers_of_rates {η C : ℝ} (h : DenseRates η C) (hC : 0≤C)
    (m Q : ℕ) (hQ : Q ∈ Ico (nativeCut m+1) (Lsched m+1))
    (seed : NativeSeed m) (hg : nativeGood m seed) :
    (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤
      densePrefactor (Lsched m) (bsched m) Q C * Real.exp (-η*(Nsched m:ℝ)) := by
  have hk : 0<nativeBlocks m := by
    have := native_width m
    have := bsched_pos m
    omega
  have hQ0 : 0<Q := by have := (mem_Ico.mp hQ).1; omega
  have hr : 128*rounds m=Lsched m*(nativeBlocks m*24) := by
    rw [native_width]
    exact rounds_length m
  have hd : (threshold m:ℝ)≤(11/100)*((Lsched m:ℝ)*(nativeBlocks m*24:ℕ)) := by
    rw [native_width]
    simpa only [mul_assoc] using threshold_le m
  have hh := h.layer_bound hC (Lsched_pos m) hk hQ0 seed hg
    (native_dense_density m Q hQ) (tupleWiring m) hr hd
  convert hh using 1
  rw [native_width,Nsched_eq,Nat.cast_mul]

/-- The exact minimum-distance theorem now requires fixed moments and mixed box certificates. -/
theorem minimum_distance_of_fixed_and_dense_rates
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0))
    {η C : ℝ} (hη : 0<η) (hC : 0<C) (hdense : DenseRates η C) :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
        (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊)) atTop (nhds 0) := by
  apply minimum_distance_of_fixed_and_dense_layers hfixed hC hη
  exact Filter.Eventually.of_forall (fun m Q hQ seed hg => dense_layers_of_rates hdense hC.le m Q hQ seed hg)

#print axioms native_dense_density
#print axioms dense_layers_of_rates
#print axioms minimum_distance_of_fixed_and_dense_rates
end Spin.Structured.ConcreteNativeFamily
