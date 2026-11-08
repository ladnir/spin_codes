import SpinCodes.Turbo.PhysicalLive22

/-!
# Exact-input-law live rows from explicit moment caps

This reuses the existing law-independent cone helpers. The input law remains
one complete packet-word distribution; no independence between feedback and
output weight is introduced. Numeric first/second-moment, feedback and Cauchy
certificates remain separate obligations for the supplied caps.
-/
noncomputable section
namespace Spin.Turbo
open Spin.UniformMixer Finset
attribute [local instance] Classical.propDecidable

def physicalInputKernel22 (P : FinPMF PacketStep16) (z : ℝ) :
    Matrix State22 State22 ℝ :=
  Matrix.of fun q r => P.expect (fun x => stepKernel z (packetStepInput22 x) q r)

def physicalInputMoment22 (P : FinPMF PacketStep16) (z : ℝ) (q : State22) : ℝ :=
  P.expect (fun x => z ^ stepWeight (packetStepInput22 x) q)

def physicalInputReturn22 (P : FinPMF PacketStep16) (z : ℝ) (q : State22) : ℝ :=
  P.expect (fun x => if stepFeedback (packetStepInput22 x) = 0 then 0
    else z ^ stepWeight (packetStepInput22 x) q)

theorem physicalInputKernel22_nonneg (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (q r : State22) : 0 ≤ physicalInputKernel22 P z q r := by
  simp only [physicalInputKernel22, Matrix.of_apply]
  apply FinPMF.expect_nonneg
  intro x
  exact mul_nonneg (pow_nonneg hz _) (update_nonneg _ _ _)

theorem physicalInputKernel22_live_entry (P : FinPMF PacketStep16) (z : ℝ)
    (q r : State22) (hq : q ≠ 0) :
    physicalInputKernel22 P z q r = P.expect (fun x =>
      if r = stepFeedback (packetStepInput22 x) then 0
      else z ^ stepWeight (packetStepInput22 x) q / (liveCount State22 : ℝ)) := by
  simp only [physicalInputKernel22, Matrix.of_apply]
  apply congrArg
  funext x
  change z ^ stepWeight (packetStepInput22 x) q * update _ q r = _
  rw [update_live _ q r hq]
  split_ifs <;> ring

theorem physicalInputKernel22_live_zero (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) (hq : q ≠ 0) :
    physicalInputKernel22 P z q 0 =
      physicalInputReturn22 P z q / (liveCount State22 : ℝ) := by
  rw [physicalInputKernel22_live_entry P z q 0 hq]
  unfold physicalInputReturn22 FinPMF.expect
  simp only [div_eq_mul_inv, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro x _
  by_cases hx : stepFeedback (packetStepInput22 x) = 0
  · simp only [hx, ite_true]; ring
  · simp only [hx, eq_comm, ite_false]; ring

theorem physicalInputKernel22_live_atom_le (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (q r : State22) (hq : q ≠ 0) :
    physicalInputKernel22 P z q r ≤
      physicalInputMoment22 P z q / (liveCount State22 : ℝ) := by
  rw [physicalInputKernel22_live_entry P z q r hq]
  calc
    _ ≤ P.expect (fun x => z ^ stepWeight (packetStepInput22 x) q /
        (liveCount State22 : ℝ)) := by
      apply FinPMF.expect_mono
      intro x
      split_ifs
      · exact div_nonneg (pow_nonneg hz _) (Nat.cast_nonneg _)
      · exact le_rfl
    _ = _ := by
      unfold physicalInputMoment22 FinPMF.expect
      simp only [div_eq_mul_inv, Finset.sum_mul, mul_assoc]

theorem physicalInputKernel22_row_sum (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) : (∑ r, physicalInputKernel22 P z q r) = physicalInputMoment22 P z q := by
  simp only [physicalInputKernel22, physicalInputMoment22, Matrix.of_apply]
  rw [← FinPMF.expect_sum]
  apply congrArg
  funext x
  simp only [stepKernel, ← Finset.mul_sum, update_sum, mul_one]

private theorem liveMeasure22_split (ν : State22 → ℝ) (Z U : ℝ)
    (hn : ∀ r, 0 ≤ ν r) (hm : (∑ r, ν r) ≤ U)
    (ha : ∀ r, ν r ≤ U / ((2 : ℝ)^22 - 1)) (hz : ν 0 ≤ Z) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, ν r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k (if k = 0 then Z else if k = 2 then U else 0) (parts k) := by
  let p0 : State22 → ℝ := fun r => if r = 0 then ν r else 0
  let p2 : State22 → ℝ := fun r => if r = 0 then 0 else ν r
  have hU : 0 ≤ U := (Finset.sum_nonneg (fun r _ => hn r)).trans hm
  have hp0 : physicalPiece22 0 Z p0 := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro r; dsimp [p0]; split_ifs; exact hn r; exact le_rfl
    · intro r hr
      have hr0 : r = 0 := by
        by_contra h; exact hr (by simp [p0, h])
      simpa [physicalClassAllows22] using hr0
    · simpa [p0] using hz
    · intro h; omega
  have hp2 : physicalPiece22 2 U p2 := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro r; dsimp [p2]; split_ifs; exact le_rfl; exact hn r
    · intro r hr
      have hr0 : r ≠ 0 := by intro h; subst r; exact hr (by simp [p2])
      simpa [physicalClassAllows22] using hr0
    · apply le_trans (Finset.sum_le_sum (fun r _ => ?_)) hm
      dsimp [p2]
      split_ifs
      · exact hn r
      · exact le_rfl
    · intro _ r
      dsimp [p2]
      split_ifs
      · exact div_nonneg hU (by norm_num)
      · exact ha r
  let parts : Fin 11 → State22 → ℝ := fun k r => if k = 0 then p0 r else if k = 2 then p2 r else 0
  refine ⟨parts, ?_, ?_⟩
  · intro r
    by_cases hr : r = 0
    · subst r; simp [parts, p0, p2]
    · simp only [parts, p0, p2, hr, ite_false]
      symm
      apply Finset.sum_eq_single 2
      · intro k _ hk2
        by_cases hk0 : k = 0 <;> simp [hk0, hk2]
      · intro h; exact (h (Finset.mem_univ 2)).elim
  · intro k
    by_cases hk0 : k = 0
    · subst k; simpa only [parts, ite_true] using hp0
    · by_cases hk2 : k = 2
      · subst k; simpa only [parts, hk0, ite_false, ite_true] using hp2
      · simpa only [parts, hk0, hk2, ite_false] using physicalPiece22_zero k 0 le_rfl

def physicalLiveCapRow22 (M R : ℝ) (k : Fin 11) : ℝ :=
  if k = 0 then R / (liveCount State22 : ℝ) else if k = 2 then M else 0

/-- One complete live-source row for any packet law, from weighted moment caps.
The caps are on the actual entering measure, permitting both class maxima
and uniformly dominated mean-state bounds. -/
theorem physicalInputKernel22_live_row_of_caps (P : FinPMF PacketStep16)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ)
    (M R : ℝ)
    (hM : (∑ q, μ q * physicalInputMoment22 P z q) ≤ b * M)
    (hR : (∑ q, μ q * physicalInputReturn22 P z q) ≤ b * R) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (physicalInputKernel22 P z) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k (b * physicalLiveCapRow22 M R k) (parts k) := by
  let ν := physicalKernelAction22 (physicalInputKernel22 P z) μ
  have hn (r : State22) : 0 ≤ ν r := by
    exact Finset.sum_nonneg (fun q _ => mul_nonneg (hμ.1 q)
      (physicalInputKernel22_nonneg P z hz q r))
  have hm : (∑ r, ν r) ≤ b * M := by
    dsimp [ν, physicalKernelAction22]
    rw [Finset.sum_comm]
    simpa only [← Finset.mul_sum, physicalInputKernel22_row_sum] using hM
  have ha (r : State22) : ν r ≤ (b * M) / (liveCount State22 : ℝ) := by
    calc
      ν r ≤ ∑ q, μ q * (physicalInputMoment22 P z q / (liveCount State22 : ℝ)) := by
        apply Finset.sum_le_sum
        intro q _
        by_cases hq : q = 0
        · subst q
          rw [physicalPiece22_live_zero hi hμ]
          simp only [zero_mul]
          exact le_rfl
        · exact mul_le_mul_of_nonneg_left (physicalInputKernel22_live_atom_le P z hz q r hq) (hμ.1 q)
      _ = (∑ q, μ q * physicalInputMoment22 P z q) / (liveCount State22 : ℝ) := by
        rw [Finset.sum_div]
        apply Finset.sum_congr rfl
        intro q _; ring
      _ ≤ _ := div_le_div_of_nonneg_right hM (Nat.cast_nonneg _)
  have hzero : ν 0 ≤ b * (R / (liveCount State22 : ℝ)) := by
    calc
      ν 0 = (∑ q, μ q * physicalInputReturn22 P z q) / (liveCount State22 : ℝ) := by
        dsimp [ν, physicalKernelAction22]
        rw [Finset.sum_div]
        apply Finset.sum_congr rfl
        intro q _
        by_cases hq : q = 0
        · subst q
          rw [physicalPiece22_live_zero hi hμ]
          simp only [zero_mul, zero_div]
        · rw [physicalInputKernel22_live_zero P z q hq]; ring
      _ ≤ (b * R) / (liveCount State22 : ℝ) := div_le_div_of_nonneg_right hR (Nat.cast_nonneg _)
      _ = _ := by ring
  obtain ⟨parts, hsum, hp⟩ := liveMeasure22_split ν
    (b * (R / (liveCount State22 : ℝ))) (b * M) hn hm
    (fun r => by simpa only [physicalLiveDenom22] using ha r) hzero
  refine ⟨parts, hsum, fun k => ?_⟩
  simpa only [physicalLiveCapRow22, mul_ite, mul_zero] using hp k

/-- Specialization to the actual exact-nonzero physical step, before thinning. -/
theorem turboPhysicalKernel22_live_row_of_caps (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ)
    (M R : ℝ)
    (hM : (∑ q, μ q * physicalInputMoment22 (packetStepLaw22 j) z q) ≤ b * M)
    (hR : (∑ q, μ q * physicalInputReturn22 (packetStepLaw22 j) z q) ≤ b * R) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k (b * physicalLiveCapRow22 M R k) (parts k) :=
  physicalInputKernel22_live_row_of_caps (packetStepLaw22 j) z hz i hi b μ hμ M R hM hR

#print axioms physicalInputKernel22_live_row_of_caps
#print axioms turboPhysicalKernel22_live_row_of_caps
end Spin.Turbo
