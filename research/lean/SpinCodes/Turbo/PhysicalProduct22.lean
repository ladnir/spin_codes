import SpinCodes.Turbo.PhysicalThinningRows22
import SpinCodes.Turbo.ConditionalFractionalMatrix

/-!
# Chronological products of complete physical envelopes

The historical `ScalableEnvelopeIterate` propagates a three-coordinate
envelope through a homogeneous matrix power. Turbo needs eleven coordinates
and a different complete row at each physical step. We reuse its induction
pattern with `physicalRowSound22_transfer`, preserving the continuous state.
This connects checked local rows to the scalar endpoint used by the retained
fractional calculation; it does not authenticate any proposed numeric row.
-/
noncomputable section
namespace Spin.Turbo
open Spin.UniformMixer

abbrev PhysicalPair22 := Matrix State22 State22 ℝ × Matrix (Fin 11) (Fin 11) ℝ

theorem physicalKernelAction22_one (μ : State22 → ℝ) :
    physicalKernelAction22 1 μ = μ := by
  exact Matrix.vecMul_one μ

theorem physicalEnvelopeApply22_one (v : Fin 11 → ℝ) :
    physicalEnvelopeApply22 v 1 = v := by
  exact Matrix.vecMul_one v

theorem physicalCone22_list_product (steps : List PhysicalPair22)
    (hs : ∀ p ∈ steps, physicalRowSound22 p.1 p.2)
    (v : Fin 11 → ℝ) (μ : State22 → ℝ) (hμ : physicalCone22 v μ) :
    physicalCone22
      (physicalEnvelopeApply22 v (steps.map Prod.snd).prod)
      (physicalKernelAction22 (steps.map Prod.fst).prod μ) := by
  induction steps generalizing v μ with
  | nil =>
    simpa only [List.map_nil, List.prod_nil,
      physicalKernelAction22_one, physicalEnvelopeApply22_one] using hμ
  | cons p steps ih =>
    simp only [List.map_cons, List.prod_cons,
      physicalKernelAction22_mul, physicalEnvelopeApply22_mul]
    exact ih (fun q hq => hs q (List.mem_cons_of_mem p hq)) _ _
      (physicalRowSound22_transfer (hs p List.mem_cons_self) hμ)

theorem physicalKernelAction22_initial (K : Matrix State22 State22 ℝ) :
    physicalKernelAction22 K physicalInitialMass22 = fun r => K 0 r := by
  funext r
  simp only [physicalKernelAction22, physicalInitialMass22,
    ite_mul, one_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

theorem physicalEnvelopeApply22_initial (B : Matrix (Fin 11) (Fin 11) ℝ) :
    physicalEnvelopeApply22 physicalInitialEnvelope22 B = fun k => B 0 k := by
  have hinit (i : Fin 11) : physicalInitialEnvelope22 i =
      if i = 0 then 1 else 0 := by
    have hi : i.val = 0 ↔ i = 0 := ⟨fun h => Fin.ext h, fun h => by simp [h]⟩
    simp only [physicalInitialEnvelope22, hi]
  funext k
  simp only [physicalEnvelopeApply22, hinit, ite_mul, one_mul, zero_mul,
    Finset.sum_ite_eq', Finset.mem_univ, ite_true]

/-- Sum every terminal state only after propagating the entire stream. -/
theorem physical_list_product_mass_le (steps : List PhysicalPair22)
    (hs : ∀ p ∈ steps, physicalRowSound22 p.1 p.2) :
    (∑ r : State22, (steps.map Prod.fst).prod 0 r) ≤
      ∑ k : Fin 11, (steps.map Prod.snd).prod 0 k := by
  have h := physicalCone22_list_product steps hs _ _ physicalInitialCone22
  rw [physicalKernelAction22_initial, physicalEnvelopeApply22_initial] at h
  exact physicalCone22_mass_le h

theorem physical_list_product_mass_nonneg (steps : List PhysicalPair22)
    (hs : ∀ p ∈ steps, physicalRowSound22 p.1 p.2) :
    0 ≤ ∑ r : State22, (steps.map Prod.fst).prod 0 r := by
  have h := physicalCone22_list_product steps hs _ _ physicalInitialCone22
  rw [physicalKernelAction22_initial, physicalEnvelopeApply22_initial] at h
  exact Finset.sum_nonneg (fun r _ => physicalCone22_nonneg h r)

/-- The exact physical endpoint is bounded before the fractional matrix
path inequality is applied to the complete envelope product. -/
theorem physical_list_product_fractional_le (steps : List PhysicalPair22)
    (hs : ∀ p ∈ steps, physicalRowSound22 p.1 p.2)
    (hB : ∀ p ∈ steps, ∀ i j, 0 ≤ p.2 i j)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) :
    (∑ r : State22, (steps.map Prod.fst).prod 0 r) ^ α ≤
      ∑ k : Fin 11,
        (steps.map (fun p => matrixEntrywiseRpow p.2 α)).prod 0 k := by
  apply (Real.rpow_le_rpow (physical_list_product_mass_nonneg steps hs)
    (physical_list_product_mass_le steps hs) hα0.le).trans
  have h := matrix_list_row_sum_fractional_le (steps.map Prod.snd)
    (by
      intro B hBmem i j
      obtain ⟨p, hp, rfl⟩ := List.mem_map.mp hBmem
      exact hB p hp i j) hα0 hα1 (0 : Fin 11)
  simpa only [List.map_map, Function.comp_def] using h

#print axioms physicalCone22_list_product
#print axioms physical_list_product_mass_le
#print axioms physical_list_product_fractional_le

/-- The retained two-step calculation: compare the exact stream with its
complete envelopes, multiply each chronological envelope pair, then apply
the fractional power. In particular, the two steps are not powered separately. -/
theorem physical_pairs_product_fractional_le
    (pairs : List (PhysicalPair22 × PhysicalPair22))
    (hs : ∀ p ∈ pairs, physicalRowSound22 p.1.1 p.1.2 ∧
      physicalRowSound22 p.2.1 p.2.2)
    (hB : ∀ p ∈ pairs, (∀ i j, 0 ≤ p.1.2 i j) ∧ (∀ i j, 0 ≤ p.2.2 i j))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) :
    (∑ r : State22,
      ((pairs.flatMap (fun p => [p.1, p.2])).map Prod.fst).prod 0 r) ^ α ≤
      ∑ k : Fin 11,
        (pairs.map (fun p => matrixEntrywiseRpow (p.1.2 * p.2.2) α)).prod 0 k := by
  let steps := pairs.flatMap (fun p => [p.1, p.2])
  have hsteps : ∀ p ∈ steps, physicalRowSound22 p.1 p.2 := by
    intro p hp
    obtain ⟨pair, hpair, hp⟩ := List.mem_flatMap.mp hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl
    · exact (hs pair hpair).1
    · exact (hs pair hpair).2
  apply (Real.rpow_le_rpow (physical_list_product_mass_nonneg steps hsteps)
    (physical_list_product_mass_le steps hsteps) hα0.le).trans
  have h := chronological_pairs_row_fractional_le
    (pairs.map (fun p => (p.1.2, p.2.2))) (by
      intro p hp
      obtain ⟨pair, hpair, rfl⟩ := List.mem_map.mp hp
      exact hB pair hpair) hα0 hα1 (0 : Fin 11)
  simpa only [steps, List.map_flatMap, List.flatMap_map, List.map_cons,
    List.map_nil, List.map_map, Function.comp_def] using h

#print axioms physical_pairs_product_fractional_le
end Spin.Turbo
