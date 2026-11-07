import SpinCodes.UniformMixer.ScalableWeightedRows

/-!
# Three-coordinate domination for a weighted state law

The coordinates record mass at zero, arbitrary mass on live states,
and a uniform pointwise allowance on live states. The free live mass
has an explicit witness, which is needed to account for its return
probability without assuming it is uniformly distributed.
-/

noncomputable section

namespace Spin.UniformMixer

def liveStateSet : Finset State22 :=
  (Finset.univ : Finset State22).erase 0

def WeightedEnvelope (μ : State22 → ℝ) (a f u : ℝ) : Prop :=
  ∃ ν : State22 → ℝ,
    (∀ q, 0 ≤ ν q) ∧
    (∑ q ∈ liveStateSet, ν q ≤ f) ∧
    μ 0 ≤ a ∧
    (∀ q : State22, q ≠ 0 →
      μ q ≤ ν q + u / (liveCount State22 : ℝ))

def weightedKernelStep (K : Matrix State22 State22 ℝ)
    (μ : State22 → ℝ) (r : State22) : ℝ :=
  ∑ q : State22, μ q * K q r

theorem state22_liveCount_real_pos :
    (0 : ℝ) < (liveCount State22 : ℝ) := by
  exact_mod_cast liveCount_pos (S := State22)

structure WeightedRowBounds (K : Matrix State22 State22 ℝ)
    (H G F U V W : ℝ) : Prop where
  kernel_nonneg : ∀ q r, 0 ≤ K q r
  H_nonneg : 0 ≤ H
  G_nonneg : 0 ≤ G
  F_nonneg : 0 ≤ F
  U_nonneg : 0 ≤ U
  V_nonneg : 0 ≤ V
  W_nonneg : 0 ≤ W
  zero_hold : K 0 0 ≤ H
  zero_live : ∑ r ∈ liveStateSet, K 0 r ≤ G
  live_atom : ∀ q r, q ≠ 0 → r ≠ 0 →
    K q r ≤ F / (liveCount State22 : ℝ)
  live_return : ∀ q, q ≠ 0 →
    K q 0 ≤ V / (liveCount State22 : ℝ)
  mean_atom : ∀ r, r ≠ 0 →
    ∑ q ∈ liveStateSet, K q r ≤ U
  mean_return : ∑ q ∈ liveStateSet, K q 0 ≤ W

theorem weighted_live_sum_le (K : Matrix State22 State22 ℝ)
    (ν : State22 → ℝ) (f c : ℝ) (r : State22)
    (hν : ∀ q, 0 ≤ ν q)
    (hmass : ∑ q ∈ liveStateSet, ν q ≤ f)
    (hc : 0 ≤ c)
    (hcap : ∀ q, q ∈ liveStateSet → K q r ≤ c) :
    ∑ q ∈ liveStateSet, ν q * K q r ≤ f * c := by
  calc
    _ ≤ ∑ q ∈ liveStateSet, ν q * c := by
      apply Finset.sum_le_sum
      intro q hq
      exact mul_le_mul_of_nonneg_left (hcap q hq) (hν q)
    _ = (∑ q ∈ liveStateSet, ν q) * c := by
      rw [Finset.sum_mul]
    _ ≤ f * c := mul_le_mul_of_nonneg_right hmass hc

theorem weightedKernelStep_le_decomp
    (K : Matrix State22 State22 ℝ) (μ ν : State22 → ℝ)
    (a u : ℝ)
    (hK : ∀ q r, 0 ≤ K q r)
    (hzero : μ 0 ≤ a)
    (hlive : ∀ q : State22, q ≠ 0 →
      μ q ≤ ν q + u / (liveCount State22 : ℝ))
    (r : State22) :
    weightedKernelStep K μ r ≤
      a * K 0 r +
      ∑ q ∈ liveStateSet, ν q * K q r +
      (u / (liveCount State22 : ℝ)) *
        ∑ q ∈ liveStateSet, K q r := by
  have hsplit := Finset.add_sum_erase (Finset.univ : Finset State22)
    (fun q => μ q * K q r) (Finset.mem_univ (0 : State22))
  change μ 0 * K 0 r + ∑ q ∈ liveStateSet, μ q * K q r =
    weightedKernelStep K μ r at hsplit
  rw [← hsplit]
  calc
    _ ≤ a * K 0 r +
        ∑ q ∈ liveStateSet,
          (ν q + u / (liveCount State22 : ℝ)) * K q r := by
      apply add_le_add (mul_le_mul_of_nonneg_right hzero (hK 0 r))
      apply Finset.sum_le_sum
      intro q hq
      exact mul_le_mul_of_nonneg_right
        (hlive q (Finset.ne_of_mem_erase hq)) (hK q r)
    _ = _ := by
      simp only [add_mul, Finset.sum_add_distrib, Finset.mul_sum]
      ring

theorem weightedEnvelope_step
    (K : Matrix State22 State22 ℝ) (H G F U V W : ℝ)
    (B : WeightedRowBounds K H G F U V W)
    (μ : State22 → ℝ) (a f u : ℝ)
    (ha : 0 ≤ a) (_hf : 0 ≤ f) (hu : 0 ≤ u)
    (hμ : WeightedEnvelope μ a f u) :
    WeightedEnvelope (weightedKernelStep K μ)
      (a * H + f * V / (liveCount State22 : ℝ) +
        u * W / (liveCount State22 : ℝ))
      (a * G) (f * F + u * U) := by
  obtain ⟨ν, hν, hmass, hzero, hlive⟩ := hμ
  let D : ℝ := (liveCount State22 : ℝ)
  have hD : 0 < D := state22_liveCount_real_pos
  have huD : 0 ≤ u / D := div_nonneg hu hD.le
  have hFdiv : 0 ≤ F / D := div_nonneg B.F_nonneg hD.le
  have hVdiv : 0 ≤ V / D := div_nonneg B.V_nonneg hD.le
  let ν' : State22 → ℝ := fun r =>
    if r = 0 then 0 else a * K 0 r
  refine ⟨ν', ?_, ?_, ?_, ?_⟩
  · intro r
    by_cases hr : r = 0
    · simp only [ν', hr, ite_true, le_refl]
    · simp only [ν', hr, ite_false]
      exact mul_nonneg ha (B.kernel_nonneg 0 r)
  · calc
      ∑ r ∈ liveStateSet, ν' r =
          a * ∑ r ∈ liveStateSet, K 0 r := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro r hr
        simp only [ν', Finset.ne_of_mem_erase hr, ite_false]
      _ ≤ a * G := mul_le_mul_of_nonneg_left B.zero_live ha
  · have hret :
        ∑ q ∈ liveStateSet, ν q * K q 0 ≤ f * (V / D) :=
      weighted_live_sum_le K ν f (V / D) 0 hν hmass hVdiv
        (fun q hq => B.live_return q (Finset.ne_of_mem_erase hq))
    have hdecomp := weightedKernelStep_le_decomp K μ ν a u
      B.kernel_nonneg hzero hlive 0
    change weightedKernelStep K μ 0 ≤
      a * K 0 0 + ∑ q ∈ liveStateSet, ν q * K q 0 +
      (u / D) * ∑ q ∈ liveStateSet, K q 0 at hdecomp
    calc
      _ ≤ a * K 0 0 + ∑ q ∈ liveStateSet, ν q * K q 0 +
          (u / D) * ∑ q ∈ liveStateSet, K q 0 := hdecomp
      _ ≤ a * H + f * (V / D) + (u / D) * W := by
        exact add_le_add
          (add_le_add (mul_le_mul_of_nonneg_left B.zero_hold ha) hret)
          (mul_le_mul_of_nonneg_left B.mean_return huD)
      _ = _ := by ring
  · intro r hr
    have hatom :
        ∑ q ∈ liveStateSet, ν q * K q r ≤ f * (F / D) :=
      weighted_live_sum_le K ν f (F / D) r hν hmass hFdiv
        (fun q hq => B.live_atom q r
          (Finset.ne_of_mem_erase hq) hr)
    have hdecomp := weightedKernelStep_le_decomp K μ ν a u
      B.kernel_nonneg hzero hlive r
    change weightedKernelStep K μ r ≤
      a * K 0 r + ∑ q ∈ liveStateSet, ν q * K q r +
      (u / D) * ∑ q ∈ liveStateSet, K q r at hdecomp
    calc
      _ ≤ a * K 0 r + ∑ q ∈ liveStateSet, ν q * K q r +
          (u / D) * ∑ q ∈ liveStateSet, K q r := hdecomp
      _ ≤ a * K 0 r + f * (F / D) + (u / D) * U := by
        exact add_le_add (add_le_add_right hatom _)
          (mul_le_mul_of_nonneg_left (B.mean_atom r hr) huD)
      _ = ν' r + (f * F + u * U) / D := by
        simp only [ν', hr, ite_false]
        ring

end Spin.UniformMixer
