import SpinCodes.UniformMixer.ScalableFixedSmallExponent
import SpinCodes.Structured.ConcreteOuterCountingProbability
import SpinCodes.Structured.ConcreteNativeSparseSelected

/-! Reuse the selected shared outer's counting argument for a new inner functional. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute Spin.Structured.ConcreteMarked
  Spin.Structured.ConcreteNativeFamily Finset Filter

private def castRegions {L b b' : ℕ} (h : b = b')
    (regions : Fin b → Finset (Fin L)) :
    Fin b' → Finset (Fin L) :=
  fun j => regions ((finCongr h).symm j)

private theorem routeLaw_width_cast {L b b' : ℕ} (h : b = b')
    (rows : Fin L → Finset (Fin b)) :
    law (fun i => (rows i).map (finCongr h).toEmbedding) =
      (law rows).map (castRegions h) := by
  cases h
  have hc : (castRegions (L := L) (rfl : b = b)) = id := by
    funext regions j
    simp [castRegions, finCongr_refl]
  rw [hc, FinPMF.map_id]
  congr 1
  funext i
  simp [finCongr_refl]

private theorem regionsLaw_width_cast {L b b' : ℕ} (h : b = b')
    (S : Finset (Fin L)) :
    (regionsLaw S b).map (castRegions h) = regionsLaw S b' := by
  cases h
  have hc : (castRegions (L := L) (rfl : b = b)) = id := by
    funext regions j
    simp [castRegions, finCongr_refl]
  rw [hc, FinPMF.map_id]

theorem scalable_activeMessages_route_expect_le (m : ℕ)
    (seed : NativeSeed m) (S : Finset (Fin (Lsched m)))
    {B : ℝ}
    (hB : ∀ w ≤ bsched m, (spectrum seed w : ℝ) ≤
      B * ((bsched m).choose w : ℝ))
    (F : (Fin (bsched m) → Finset (Fin (Lsched m))) → ℝ)
    (hF : ∀ regions, 0 ≤ F regions) :
    (∑ x ∈ activeMessages S,
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect F) ≤
      ((2 : ℝ) ^ bsched m * B) ^ S.card *
        (regionsLaw S (bsched m)).expect F := by
  let G : (Fin (nativeBlocks m * 24) →
      Finset (Fin (Lsched m))) → ℝ :=
    fun regions => F (castRegions (native_width m) regions)
  have hG : ∀ regions, 0 ≤ G regions := fun regions => hF _
  have hh := activeMessages_route_expect_le seed S
    (by simpa only [native_width] using hB) G hG
  have he (x : NativeMessage m) :
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect F =
      (law (rowSupports seed x)).expect G := by
    rw [routeLaw_width_cast (native_width m), FinPMF.expect_map]
  have hf : (regionsLaw S (nativeBlocks m * 24)).expect G =
      (regionsLaw S (bsched m)).expect F := by
    rw [← regionsLaw_width_cast (native_width m) S,
      FinPMF.expect_map]
  simpa only [he, hf, native_width] using hh

def scalableRegionFailure (m d : ℕ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) : ℝ :=
  (scalableStateMapLaw m).prob (fun maps =>
    sampledOutputWeight (scalableFlatInput m regions)
      (scalableFlatMaps m maps) 0 ≤ d)

theorem scalableRegionFailure_nonneg (m d : ℕ) (regions) :
    0 ≤ scalableRegionFailure m d regions := FinPMF.prob_nonneg _ _

theorem scalableRegionFailure_fair (m d : ℕ)
    (S : Finset (Fin (Lsched m))) :
    (regionsLaw S (bsched m)).expect (scalableRegionFailure m d) =
      ((regionsLaw S (bsched m)).prod (scalableStateMapLaw m)).prob
        (fun ω => sampledOutputWeight (scalableFlatInput m ω.1)
          (scalableFlatMaps m ω.2) 0 ≤ d) := by
  exact (FinPMF.prob_prod_eq_expect _ _ _).symm

theorem scalable_occupation_failure_le_fair (m : ℕ)
    (seed : NativeSeed m) (Q d : ℕ) {B P : ℝ}
    (hB : ∀ w ≤ bsched m, (spectrum seed w : ℝ) ≤
      B * ((bsched m).choose w : ℝ))
    (hP : ∀ S : Finset (Fin (Lsched m)), S.card = Q →
      (regionsLaw S (bsched m)).expect (scalableRegionFailure m d) ≤ P) :
    (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d)) ≤
      ((Lsched m).choose Q : ℝ) *
        (((2 : ℝ) ^ bsched m * B) ^ Q * P) := by
  have hB0 : 0 ≤ B := by
    have hh := hB 0 (by omega)
    simp only [Nat.choose_zero_right, Nat.cast_one, mul_one] at hh
    exact (Nat.cast_nonneg _).trans hh
  rw [sum_occupation_eq_activeSets]
  calc
    _ ≤ ∑ _S ∈ (univ : Finset (Fin (Lsched m))).powersetCard Q,
        ((2 : ℝ) ^ bsched m * B) ^ Q * P := by
      apply sum_le_sum
      intro S hS
      have hc := (mem_powersetCard.mp hS).2
      have hh := (scalable_activeMessages_route_expect_le m seed S hB
        (scalableRegionFailure m d) (scalableRegionFailure_nonneg m d)).trans
        (mul_le_mul_of_nonneg_left (hP S hc) (by positivity))
      simpa only [hc] using hh
    _ = _ := by simp [card_powersetCard]

theorem scalable_native_one_tuple_bound :
    ∀ᶠ m : ℕ in atTop, ∀ seed : NativeSeed m,
      nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = 1),
        (law (fun i => (rowSupports seed x i).map
          (finCongr (native_width m)).toEmbedding)).expect
            (scalableRegionFailure m (threshold m))) ≤
        scalableFixedOneBound m := by
  filter_upwards [scalable_fair_one_probability] with m hm
  intro seed hg
  have hB := native_good_spectrum_bound m seed hg
    (B := spectrumRatioBound (bsched m) (nativeShellRemainder m))
    (by unfold spectrumRatioBound; positivity)
    (native_expected_spectrum_uniform m)
  have hp : ∀ S : Finset (Fin (Lsched m)), S.card = 1 →
      (regionsLaw S (bsched m)).expect
        (scalableRegionFailure m (threshold m)) ≤
      (1000 * (1011 / 2000 : ℝ) ^ bsched m) /
        Real.exp (-(scalableThetaOne / (Lsched m : ℝ))) ^
          threshold m := by
    intro S hS
    rw [scalableRegionFailure_fair]
    exact hm S hS (threshold m)
  have hh := scalable_occupation_failure_le_fair m seed 1
    (threshold m) (by simpa only [native_width] using hB) hp
  simpa only [Nat.choose_one_right, pow_one,
    selected_row_cost, scalableFixedOneBound,
    mul_div_assoc, mul_assoc] using hh

theorem scalable_native_two_tuple_bound :
    ∀ᶠ m : ℕ in atTop, ∀ seed : NativeSeed m,
      nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = 2),
        (law (fun i => (rowSupports seed x i).map
          (finCongr (native_width m)).toEmbedding)).expect
            (scalableRegionFailure m (threshold m))) ≤
        scalableFixedTwoBound m := by
  filter_upwards [scalable_fair_two_probability] with m hm
  intro seed hg
  have hB := native_good_spectrum_bound m seed hg
    (B := spectrumRatioBound (bsched m) (nativeShellRemainder m))
    (by unfold spectrumRatioBound; positivity)
    (native_expected_spectrum_uniform m)
  have hp : ∀ S : Finset (Fin (Lsched m)), S.card = 2 →
      (regionsLaw S (bsched m)).expect
        (scalableRegionFailure m (threshold m)) ≤
      (1000 * (1011 / 4000 : ℝ) ^ bsched m) /
        Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))) ^
          threshold m := by
    intro S hS
    rw [scalableRegionFailure_fair]
    exact hm S hS (threshold m)
  have hh := scalable_occupation_failure_le_fair m seed 2
    (threshold m) (by simpa only [native_width] using hB) hp
  have hc : ((Lsched m).choose 2 : ℝ) ≤ (Lsched m : ℝ) ^ 2 := by
    exact_mod_cast Nat.choose_le_pow (Lsched m) 2
  have hh' := hh.trans (mul_le_mul_of_nonneg_right hc (by positivity))
  simpa only [selected_row_cost, scalableFixedTwoBound,
    ← Real.exp_nat_mul, Nat.cast_ofNat,
    mul_div_assoc, mul_assoc] using hh'

end Spin.UniformMixer
