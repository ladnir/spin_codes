import SpinCodes.UniformMixer.ScalableSparseWitness
import SpinCodes.UniformMixer.ScalableDualKernelDistance

/-!
# Elementary bounds for the sparse transfer

The certified dual distance controls the zero-to-zero entry. We use the
simple count of all 64-bit supports for the tail; its constant is loose but
the sparse interval is small enough that it may suffice.
-/

noncomputable section

namespace Spin.UniformMixer

def zeroHoldTerm (β z : ℝ) (S : Finset (Fin 64)) : ℝ :=
  if stepFeedback (supportInput S) = 0 then
    (β * z) ^ S.card * (1 - β) ^ (64 - S.card) else 0

theorem zeroHoldTerm_empty (β z : ℝ) :
    zeroHoldTerm β z ∅ = (1 - β) ^ 64 := by
  have hzero : stepFeedback (supportInput ∅) = 0 := by
    funext i
    simp [stepFeedback, supportInput]
  simp [zeroHoldTerm, hzero]

theorem zeroHoldTerm_le_beta8 (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (S : Finset (Fin 64)) (hS : S ≠ ∅) :
    zeroHoldTerm β z S ≤ β ^ 8 := by
  unfold zeroHoldTerm
  split_ifs with hkernel
  · have hdist : 8 ≤ S.card :=
      feedbackKernel_distance S hkernel (Finset.nonempty_iff_ne_empty.mpr hS)
    have hβz0 : 0 ≤ β * z := mul_nonneg hβ0 hz0
    have hβz1 : β * z ≤ β := by nlinarith
    have hpow : (β * z) ^ S.card ≤ β ^ S.card :=
      pow_le_pow_left₀ hβz0 hβz1 _
    have htail : (1 - β) ^ (64 - S.card) ≤ 1 :=
      pow_le_one₀ (by linarith) (by linarith)
    calc
      (β * z) ^ S.card * (1 - β) ^ (64 - S.card) ≤
          β ^ S.card * 1 := by
        exact mul_le_mul hpow htail (by positivity) (pow_nonneg hβ0 _)
      _ = β ^ S.card := by ring
      _ ≤ β ^ 8 := pow_le_pow_of_le_one hβ0 hβ1 hdist
  · exact pow_nonneg hβ0 _

theorem scalableZeroHold_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    scalableZeroHold β z hβ0 hβ1 ≤
      (1 - β) ^ 64 + (2 ^ 64 : ℝ) * β ^ 8 := by
  have htail :
      (∑ S ∈ (Finset.univ : Finset (Finset (Fin 64))).erase ∅,
        zeroHoldTerm β z S) ≤ (2 ^ 64 : ℝ) * β ^ 8 := by
    calc
      _ ≤ ∑ S ∈ (Finset.univ : Finset (Finset (Fin 64))).erase ∅,
            β ^ 8 := by
        apply Finset.sum_le_sum
        intro S hS
        exact zeroHoldTerm_le_beta8 β z hβ0 hβ1 hz0 hz1 S
          (Finset.ne_of_mem_erase hS)
      _ = (((Finset.univ : Finset (Finset (Fin 64))).erase ∅).card : ℝ) *
            β ^ 8 := by simp
      _ ≤ ((Finset.univ : Finset (Finset (Fin 64))).card : ℝ) *
            β ^ 8 := by
        apply mul_le_mul_of_nonneg_right _ (pow_nonneg hβ0 _)
        exact_mod_cast Finset.card_le_card (Finset.erase_subset _ _)
      _ = (2 ^ 64 : ℝ) * β ^ 8 := by
        simp [Fintype.card_finset]
        norm_num
  rw [scalableZeroHold_kernel_sum]
  change (∑ S : Finset (Fin 64), zeroHoldTerm β z S) ≤ _
  have hsplit := Finset.add_sum_erase
    (Finset.univ : Finset (Finset (Fin 64)))
    (zeroHoldTerm β z) (Finset.mem_univ (∅ : Finset (Fin 64)))
  rw [← hsplit, zeroHoldTerm_empty]
  linarith

theorem one_sub_pow_le_quadratic (n : ℕ) (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    (1 - β) ^ n ≤
      1 - (n : ℝ) * β + (Nat.choose n 2 : ℝ) * β ^ 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
      have hchoose : Nat.choose (n + 1) 2 = Nat.choose n 2 + n := by
        simp [Nat.choose_succ_succ, Nat.add_comm]
      have hmul := mul_le_mul_of_nonneg_right ih (by linarith : 0 ≤ 1 - β)
      have htail : 0 ≤ (Nat.choose n 2 : ℝ) * β ^ 3 := by positivity
      have hpoly :
          (1 - (n : ℝ) * β + (Nat.choose n 2 : ℝ) * β ^ 2) *
              (1 - β) =
            1 - ((n + 1 : ℕ) : ℝ) * β +
              (Nat.choose (n + 1) 2 : ℝ) * β ^ 2 -
              (Nat.choose n 2 : ℝ) * β ^ 3 := by
        rw [hchoose]
        push_cast
        ring
      rw [pow_succ]
      linarith

set_option maxRecDepth 8192 in
theorem scalableZeroHold_sparse_quadratic (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    scalableZeroHold β z hβ0 hβ1 ≤
      1 - 64 * β + 2016 * β ^ 2 + (2 ^ 64 : ℝ) * β ^ 8 := by
  have h := scalableZeroHold_sparse_upper β z hβ0 hβ1 hz0 hz1
  have hquad := one_sub_pow_le_quadratic 64 β hβ0 hβ1
  have hchoose64 : (Nat.choose 64 2 : ℝ) = 2016 := by
    exact_mod_cast (show Nat.choose 64 2 = 2016 by decide)
  have hquad' : (1 - β) ^ 64 ≤ 1 - 64 * β + 2016 * β ^ 2 := by
    simpa only [hchoose64, Nat.cast_ofNat] using hquad
  exact h.trans (add_le_add_left hquad' ((2 ^ 64 : ℝ) * β ^ 8))

theorem zeroHoldTerm_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz0 : 0 ≤ z)
    (S : Finset (Fin 64)) : 0 ≤ zeroHoldTerm β z S := by
  unfold zeroHoldTerm
  split_ifs <;> positivity

theorem scalableZeroHold_sparse_lower (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz0 : 0 ≤ z) :
    (1 - β) ^ 64 ≤ scalableZeroHold β z hβ0 hβ1 := by
  rw [scalableZeroHold_kernel_sum]
  change (1 - β) ^ 64 ≤
    ∑ S : Finset (Fin 64), zeroHoldTerm β z S
  rw [← zeroHoldTerm_empty β z]
  exact Finset.single_le_sum
    (fun S _ => zeroHoldTerm_nonneg β z hβ0 hβ1 hz0 S)
    (Finset.mem_univ ∅)

theorem scalableZeroLive_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    (1 - β + β * z) ^ 64 - scalableZeroHold β z hβ0 hβ1 ≤
      64 * β := by
  have hH := scalableZeroHold_sparse_lower β z hβ0 hβ1 hz0
  have hg0 : (1 - β + β * z) ^ 64 ≤ 1 := by
    apply pow_le_one₀
    · nlinarith
    · nlinarith
  have hbern : 1 - 64 * β ≤ (1 - β) ^ 64 := by
    have h := one_add_mul_le_pow (a := -β) (by linarith : -2 ≤ -β) 64
    simpa only [Nat.cast_ofNat, mul_neg, sub_eq_add_neg] using h
  calc
    (1 - β + β * z) ^ 64 - scalableZeroHold β z hβ0 hβ1 ≤
        1 - scalableZeroHold β z hβ0 hβ1 :=
      sub_le_sub_right hg0 _
    _ ≤ 1 - (1 - β) ^ 64 := sub_le_sub_left hH _
    _ ≤ 64 * β := by linarith only [hbern]

theorem spectrumG_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    spectrumG β z ≤ 64 * β := by
  have hspec : spectrumH β z = scalableZeroHold β z hβ0 hβ1 :=
    (scalableZeroHold_exact_spectrum β z hβ0 hβ1).symm
  unfold spectrumG
  rw [hspec]
  exact scalableZeroLive_sparse_upper β z hβ0 hβ1 hz0 hz1

end Spin.UniformMixer
