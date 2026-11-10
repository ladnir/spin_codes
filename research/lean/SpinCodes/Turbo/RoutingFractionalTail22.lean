import SpinCodes.Turbo.RoutingFractionalEndpoint22

/-!
# Occupation union with independent witnesses per occupation

Audit: the historical fixed-coefficient proof averages placement weights
linearly; the fractional cap applies to each active set before any union.
Accordingly the binomial cardinality is outside the fractional power.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem turboActiveBadCount22_prob_union_le_paired_dp
    (n : ℕ) (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (d : ℕ) (qs : Finset ℕ)
    (z α : ℕ → ℝ) (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hz0 : ∀ q ∈ qs, 0 < z q) (hz1 : ∀ q ∈ qs, z q ≤ 1)
    (ha0 : ∀ q ∈ qs, 0 < α q) (ha1 : ∀ q ∈ qs, α q ≤ 1)
    (hs : ∀ q ∈ qs, ∀ j,
      physicalRowSound22 (Matrix.of (turboSelectedKernel22 j (z q))) (B q j))
    (hB : ∀ q ∈ qs, ∀ j a b, 0 ≤ B q j a b) :
    (turboSetup22 (n * 2) points).joint.prob (fun ω =>
      ∃ q ∈ qs, ∃ S ∈ (univ : Finset (Fin ((n * 2) * 16))).powersetCard q,
        1 ≤ turboActiveBadCount22 (n * 2) points d S ω) ≤
      ∑ q ∈ qs, (Nat.choose ((n * 2) * 16) q : ℝ) *
        ((turboOldUniformFactor22 ^ q / z q ^ d) ^ α q *
          ∑ k : Fin 11, ((pairedRegionDP22 (B q) (α q) n q) ^ 64) 0 k) := by
  apply (FinPMF.prob_biUnion_le _ qs _).trans
  apply sum_le_sum
  intro q hq
  apply (FinPMF.prob_biUnion_le _
    ((univ : Finset (Fin ((n * 2) * 16))).powersetCard q) _).trans
  calc
    _ ≤ ∑ S ∈ (univ : Finset (Fin ((n * 2) * 16))).powersetCard q,
        ((turboOldUniformFactor22 ^ q / z q ^ d) ^ α q *
          ∑ k : Fin 11, ((pairedRegionDP22 (B q) (α q) n q) ^ 64) 0 k) := by
      apply sum_le_sum
      intro S hS
      have hc : S.card = q := (mem_powersetCard.mp hS).2
      have h := turboActiveBadCount22_prob_le_paired_dp n points hpoints S d
        (z q) (hz0 q hq) (hz1 q hq) (B q) (hs q hq) (hB q hq) (ha0 q hq) (ha1 q hq)
      simpa only [hc] using h
    _ = _ := by simp only [sum_const, card_powersetCard, card_univ, Fintype.card_fin,
        nsmul_eq_mul, Nat.cast_id]

/-- The union event covers the actual encoder tail, with every message using
one shared sampled setup. Taking qs = Ico 3 (32*n+1) selects q ≥ 3. -/
theorem turboOccupationTail22_prob_le_active_union
    (n : ℕ) (points : Fin 8 → Field16) (d : ℕ) (qs : Finset ℕ) :
    (turboSetup22 (n * 2) points).joint.prob (fun ω =>
      ∃ msg : TurboMessage22 (n * 2), turboOccupation22 msg ∈ qs ∧
        (turboSetup22 (n * 2) points).innerWt ω.2
          ((turboSetup22 (n * 2) points).outer ω.1 msg) ≤ d) ≤
    (turboSetup22 (n * 2) points).joint.prob (fun ω =>
      ∃ q ∈ qs, ∃ S ∈ (univ : Finset (Fin ((n * 2) * 16))).powersetCard q,
        1 ≤ turboActiveBadCount22 (n * 2) points d S ω) := by
  apply FinPMF.prob_mono
  intro ω h
  obtain ⟨msg, hq, hw⟩ := h
  refine ⟨turboOccupation22 msg, hq, turboActivePositions22 msg, ?_, ?_⟩
  · exact mem_powersetCard.mpr ⟨subset_univ _, rfl⟩
  · apply card_pos.mpr
    exact ⟨msg, mem_filter.mpr ⟨mem_filter.mpr ⟨mem_univ _, rfl⟩, hw⟩⟩
#print axioms turboActiveBadCount22_prob_union_le_paired_dp
end Spin.Turbo

