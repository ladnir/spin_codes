import SpinCodes.UniformMixer.ScalableLinearCode
import SpinCodes.UniformMixer.ScalableOriginalSetupFinal

/-! Minimum distance of the actual sampled Scalable binary linear code. -/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding Finset Filter
attribute [local instance] Classical.propDecidable

/-- All words emitted by this one shared sampled setup. -/
def scalableCodewords (m : ℕ) (seed : ScalableSetupSeed m) :
    Finset (ScalableOutput m) := univ.image (scalableCodeword m seed)

@[simp] theorem mem_scalableRealizedCode (m : ℕ) (seed : ScalableSetupSeed m)
    (c : ScalableOutput m) :
    c ∈ scalableRealizedCode m seed ↔ c ∈ scalableCodewords m seed := by
  simp only [scalableRealizedCode, LinearMap.mem_range, scalableCodewords,
    mem_image, mem_univ, true_and]
  rfl

lemma scalableCodewords_zero_mem (m : ℕ) (seed : ScalableSetupSeed m) :
    (0 : ScalableOutput m) ∈ scalableCodewords m seed :=
  (mem_scalableRealizedCode m seed _).mp (Submodule.zero_mem _)

lemma scalableCodewords_sub_mem {m : ℕ} {seed : ScalableSetupSeed m}
    {x y : ScalableOutput m} (hx : x ∈ scalableCodewords m seed)
    (hy : y ∈ scalableCodewords m seed) :
    x - y ∈ scalableCodewords m seed :=
  (mem_scalableRealizedCode m seed _).mp ((scalableRealizedCode m seed).sub_mem
    ((mem_scalableRealizedCode m seed _).mpr hx)
    ((mem_scalableRealizedCode m seed _).mpr hy))

def scalableHammingDist {m : ℕ} (x y : ScalableOutput m) : ℕ :=
  ∑ i, ∑ j, if x i j ≠ y i j then 1 else 0

lemma sampledWordWeight_sub {m : ℕ} (x y : ScalableOutput m) :
    sampledWordWeight (x - y) = scalableHammingDist x y := by
  have hb : ∀ a b : ZMod 2,
      (if a - b = 1 then (1 : ℕ) else 0) =
        (if a ≠ b then 1 else 0) := by decide
  unfold sampledWordWeight scalableHammingDist
  apply sum_congr rfl
  intro i hi
  apply sum_congr rfl
  intro j hj
  exact hb (x i j) (y i j)

/-- Low nonzero weight equals low distance between distinct emitted words. -/
theorem scalable_low_weight_iff_close_pair (m : ℕ) (seed : ScalableSetupSeed m)
    (d : ℕ) :
    (∃ c ∈ scalableCodewords m seed, c ≠ 0 ∧ sampledWordWeight c ≤ d) ↔
      ∃ x ∈ scalableCodewords m seed, ∃ y ∈ scalableCodewords m seed,
        x ≠ y ∧ scalableHammingDist x y ≤ d := by
  constructor
  · rintro ⟨c, hc, hn, hw⟩
    refine ⟨c, hc, 0, scalableCodewords_zero_mem m seed, hn, ?_⟩
    simpa only [← sampledWordWeight_sub, sub_zero] using hw
  · rintro ⟨x, hx, y, hy, hne, hd⟩
    exact ⟨x - y, scalableCodewords_sub_mem hx hy, sub_ne_zero.mpr hne,
      (sampledWordWeight_sub x y).symm ▸ hd⟩

def scalableNonzeroCodewords (m : ℕ) (seed : ScalableSetupSeed m) :
    Finset (ScalableOutput m) := (scalableCodewords m seed).filter (fun c => c ≠ 0)

lemma scalableNonzeroCodewords_nonempty (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalableNonzeroCodewords m seed).Nonempty := by
  have htwice : Nsched m / 2 * 2 = Nsched m := by
    apply Nat.div_mul_cancel
    rw [Nsched_eq]
    apply dvd_mul_of_dvd_left
    unfold Lsched
    omega
  have hK : 0 < Nsched m / 2 := by have := Nsched_pos m; omega
  let x : Fin (Nsched m / 2) → ZMod 2 := fun _ => 1
  have hx : x ≠ 0 := by
    intro h
    have hh := congrFun h ⟨0, hK⟩
    norm_num [x] at hh
  refine ⟨scalableCodeword m seed x,
    mem_filter.mpr ⟨mem_image.mpr ⟨x, mem_univ _, rfl⟩, ?_⟩⟩
  intro hz
  exact hx (scalableCodeword_injective m seed
    (hz.trans (scalableCodeword_zero m seed).symm))

/-- Minimum binary Hamming distance of the actual realized linear code. -/
def scalableMinimumDistance (m : ℕ) (seed : ScalableSetupSeed m) : ℕ :=
  (scalableNonzeroCodewords m seed).inf'
    (scalableNonzeroCodewords_nonempty m seed) sampledWordWeight

theorem scalableMinimumDistance_le_iff (m : ℕ) (seed : ScalableSetupSeed m)
    (d : ℕ) :
    scalableMinimumDistance m seed ≤ d ↔
      ∃ c ∈ scalableCodewords m seed, c ≠ 0 ∧ sampledWordWeight c ≤ d := by
  simp only [scalableMinimumDistance, inf'_le_iff, scalableNonzeroCodewords, mem_filter]
  aesop

theorem scalableMinimumDistance_le_iff_pair (m : ℕ) (seed : ScalableSetupSeed m)
    (d : ℕ) :
    scalableMinimumDistance m seed ≤ d ↔
      ∃ x ∈ scalableCodewords m seed, ∃ y ∈ scalableCodewords m seed,
        x ≠ y ∧ scalableHammingDist x y ≤ d :=
  (scalableMinimumDistance_le_iff m seed d).trans
    (scalable_low_weight_iff_close_pair m seed d)

theorem scalableMinimumDistance_gt_iff (m : ℕ) (seed : ScalableSetupSeed m)
    (d : ℕ) :
    d < scalableMinimumDistance m seed ↔
      ∀ x ∈ scalableCodewords m seed, ∀ y ∈ scalableCodewords m seed,
        x ≠ y → d < scalableHammingDist x y := by
  rw [← not_le, scalableMinimumDistance_le_iff_pair]
  push_neg
  rfl

theorem scalableBad_iff_outputWeight (m : ℕ) (seed : ScalableSetupSeed m) :
    scalableBad m seed ↔
      ∃ x : Fin (Nsched m / 2) → Bool, x ≠ 0 ∧
        scalableOutputWeight m x seed ≤
          Spin.Structured.ConcreteNativeFamily.threshold m := by
  constructor
  · rintro ⟨x, hx, hw⟩
    refine ⟨(nativeMessageEquiv m).symm x, ?_, hw⟩
    intro hz
    have hh := congrArg (nativeMessageEquiv m) hz
    simp only [Equiv.apply_symm_apply, native_message_zero] at hh
    exact hx hh
  · rintro ⟨x, hx, hw⟩
    refine ⟨nativeMessageEquiv m x, ?_, ?_⟩
    · intro hz
      exact hx ((nativeMessageEquiv m).injective
        (hz.trans (native_message_zero m).symm))
    · simpa only [scalableOutputWeight, Equiv.symm_apply_apply] using hw

/-- The setup bad event concerns the actual nonzero emitted codewords. -/
theorem scalableBad_iff_codewords (m : ℕ) (seed : ScalableSetupSeed m) :
    scalableBad m seed ↔
      ∃ c ∈ scalableCodewords m seed, c ≠ 0 ∧
        sampledWordWeight c ≤ Spin.Structured.ConcreteNativeFamily.threshold m := by
  rw [scalableBad_iff_outputWeight]
  constructor
  · rintro ⟨x, hx, hw⟩
    let y := wordEquiv (Nsched m / 2) x
    refine ⟨scalableCodeword m seed y,
      mem_image.mpr ⟨y, mem_univ _, rfl⟩, ?_, ?_⟩
    · intro hz
      have hy := scalableCodeword_injective m seed
        (hz.trans (scalableCodeword_zero m seed).symm)
      exact hx ((wordEquiv (Nsched m / 2)).injective
        (hy.trans (wordEquiv_zero _).symm))
    · simpa only [scalableCodeword_weight, y, Equiv.symm_apply_apply] using hw
  · rintro ⟨c, hc, hn, hw⟩
    obtain ⟨x, _, rfl⟩ := mem_image.mp hc
    refine ⟨(wordEquiv (Nsched m / 2)).symm x, ?_, ?_⟩
    · intro hz
      have hx := congrArg (wordEquiv (Nsched m / 2)) hz
      simp only [Equiv.apply_symm_apply, wordEquiv_zero] at hx
      exact hn (hx ▸ scalableCodeword_zero m seed)
    · rwa [scalableCodeword_weight] at hw

theorem scalableBad_iff_minimumDistance (m : ℕ) (seed : ScalableSetupSeed m) :
    scalableBad m seed ↔ scalableMinimumDistance m seed ≤
      Spin.Structured.ConcreteNativeFamily.threshold m :=
  (scalableBad_iff_codewords m seed).trans
    (scalableMinimumDistance_le_iff m seed _).symm

theorem scalable_minimum_distance_failure_tendsto_of_dense_rates {η C : ℝ}
    (h : ScalableDenseRates η C) (hη : 0 < η) (hC : 0 < C) :
    Tendsto (fun m => (scalableSetupLaw m).prob (fun seed =>
      scalableMinimumDistance m seed ≤
        Spin.Structured.ConcreteNativeFamily.threshold m)) atTop (nhds 0) := by
  have he : (fun m => (scalableSetupLaw m).prob (fun seed =>
      scalableMinimumDistance m seed ≤
        Spin.Structured.ConcreteNativeFamily.threshold m)) =
      (fun m => (scalableSetupLaw m).prob (scalableBad m)) := by
    funext m
    exact FinPMF.prob_congr _ (fun seed => (scalableBad_iff_minimumDistance m seed).symm)
  rw [he]
  exact scalable_failure_tendsto_of_dense_rates h hη hC

#print axioms scalableBad_iff_minimumDistance
#print axioms scalable_minimum_distance_failure_tendsto_of_dense_rates

end Spin.UniformMixer
