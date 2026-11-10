import SpinCodes.Turbo.PhysicalPacketMoment22
import SpinCodes.UniformMixer.Map22PacketProfileTable

/-!
Profile grouping and state averages follow the historical Scalable spectrum
mean/max route: first authenticate the support, then aggregate or bound the
observable. Here a profile retains all sixteen four-bit packet weights.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer
set_option maxHeartbeats 1000000

local instance profilePhysicalClassDecidable22 (i : Fin 11) (q : State22) :
    Decidable (physicalClassAllows22 i q) := by
  unfold physicalClassAllows22
  infer_instance

def physicalProfile22 (q : State22) : Fin 5 → Nat :=
  packetProfile64 (fun j => coordMap j q)

def physicalSlotWeight22 (q : State22) (v : Fin 16) : Nat :=
  profileNibbleWordWeight (fun vu => coordMap (finProdFinEquiv vu) q) v

theorem physicalSlotWeight22_bound (q : State22) (v : Fin 16) :
    physicalSlotWeight22 q v ≤ 4 := by
  unfold physicalSlotWeight22 profileNibbleWordWeight
  calc
    _ ≤ ∑ _u : Fin 4, 1 := sum_le_sum (fun u _ => by split_ifs <;> omega)
    _ = 4 := by simp

def physicalSlotWeightFin22 (q : State22) (v : Fin 16) : Fin 5 :=
  ⟨physicalSlotWeight22 q v, by have := physicalSlotWeight22_bound q v; omega⟩

theorem physicalProfile22_fiber (q : State22) (r : Fin 5) :
    (univ.filter (fun v => physicalSlotWeightFin22 q v = r)).card = physicalProfile22 q r := by
  unfold physicalProfile22 packetProfile64 profileOfNibbles
  apply congrArg Finset.card
  apply Finset.filter_congr
  intro v _
  exact Fin.ext_iff

def physicalProfileProduct22 {M : Type*} [CommMonoid M]
    (a : Fin 5 → M) (p : Fin 5 → Nat) : M := ∏ r, a r ^ p r

theorem physicalProfile22_product {M : Type*} [CommMonoid M]
    (q : State22) (a : Fin 5 → M) :
    (∏ v : Fin 16, a (physicalSlotWeightFin22 q v)) =
      physicalProfileProduct22 a (physicalProfile22 q) := by
  unfold physicalProfileProduct22
  rw [← Finset.prod_fiberwise' univ (physicalSlotWeightFin22 q) a]
  apply prod_congr rfl
  intro r _
  rw [prod_const, physicalProfile22_fiber]

theorem physicalProfile22_product_nat {M : Type*} [CommMonoid M]
    (q : State22) (a : Nat → M) :
    (∏ v : Fin 16, a (physicalSlotWeight22 q v)) =
      physicalProfileProduct22 (fun r => a r.val) (physicalProfile22 q) :=
  physicalProfile22_product q (fun r => a r.val)

theorem physicalProfile22_product_sum {R : Type*} [CommSemiring R]
    (a : Fin 5 → R) :
    (∑ q : State22, ∏ v : Fin 16, a (physicalSlotWeightFin22 q v)) =
      (productionPacketProfileTable.map (fun e =>
        (e.1 : R) * physicalProfileProduct22 a e.2)).sum := by
  simp_rw [physicalProfile22_product]
  simpa only [physicalProfile22, nsmul_eq_mul] using
    productionPacketProfile_table_sum (physicalProfileProduct22 a)

def physicalProfileWeight22 (p : Fin 5 → Nat) : Nat := ∑ r, r.val * p r

theorem physicalProfile22_weight (q : State22) :
    physicalProfileWeight22 (physicalProfile22 q) = coordinateWeight coordMap q := by
  have hg := Finset.sum_fiberwise' univ (physicalSlotWeightFin22 q) (fun r : Fin 5 => r.val)
  change (∑ r : Fin 5, ∑ v : Fin 16 with physicalSlotWeightFin22 q v = r, r.val) =
    ∑ v : Fin 16, physicalSlotWeight22 q v at hg
  have hh : (∑ v : Fin 16, physicalSlotWeight22 q v) = coordinateWeight coordMap q := by
    unfold physicalSlotWeight22 profileNibbleWordWeight coordinateWeight
    rw [← Equiv.sum_comp (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64), Fintype.sum_prod_type]
  rw [← hh, ← hg]
  unfold physicalProfileWeight22
  apply sum_congr rfl
  intro r _
  rw [sum_const, physicalProfile22_fiber, Nat.nsmul_eq_mul, Nat.mul_comm]

theorem physicalProfile22_zero : physicalProfile22 0 = ![16,0,0,0,0] := by
  funext r
  fin_cases r <;> simp [physicalProfile22, packetProfile64, profileOfNibbles,
    profileNibbleWordWeight]

theorem physicalProfile22_mem_table (q : State22) :
    ∃ e ∈ productionPacketProfileTable, e.2 = physicalProfile22 q := by
  classical
  have hp : 0 < ∑ s : State22, if physicalProfile22 s = physicalProfile22 q then (1 : Nat) else 0 := by
    apply sum_pos'
    · intro s _; split_ifs <;> omega
    · exact ⟨q, mem_univ q, by simp⟩
  have hs := productionPacketProfile_table_sum
    (fun p => if p = physicalProfile22 q then (1 : Nat) else 0)
  change (∑ s : State22, if physicalProfile22 s = physicalProfile22 q then (1 : Nat) else 0) = _ at hs
  rw [hs, List.sum_pos_iff_exists_pos_nat] at hp
  obtain ⟨x, hx, hpos⟩ := hp
  obtain ⟨e, he, rfl⟩ := List.mem_map.mp hx
  refine ⟨e, he, ?_⟩
  by_contra hn
  simp only [hn, ite_false, smul_zero] at hpos
  omega

def physicalProfileClassAllows22 (i : Fin 11) (p : Fin 5 → Nat) : Prop :=
  if i.val = 0 then physicalProfileWeight22 p = 0 else
  if i.val = 1 ∨ i.val = 2 then physicalProfileWeight22 p ≠ 0 else
    physicalProfileWeight22 p = physicalClassWeight22 i

instance physicalProfileClassAllows22_decidable (i : Fin 11) (p : Fin 5 → Nat) :
    Decidable (physicalProfileClassAllows22 i p) := by
  unfold physicalProfileClassAllows22
  infer_instance

theorem physicalProfile22_class (i : Fin 11) (q : State22) :
    physicalProfileClassAllows22 i (physicalProfile22 q) ↔ physicalClassAllows22 i q := by
  have hz : coordinateWeight coordMap q = 0 ↔ q = 0 :=
    ⟨productionWeight_zero_only q, fun h => by subst q; exact coordinateWeight_zero _⟩
  unfold physicalProfileClassAllows22 physicalClassAllows22
  rw [physicalProfile22_weight]
  split_ifs <;> simp only [hz, ne_eq]

theorem physicalProfile22_full_sum (f : (Fin 5 → Nat) → ℝ) :
    (∑ q : State22, f (physicalProfile22 q)) =
      (productionPacketProfileTable.map (fun e => (e.1 : ℝ) * f e.2)).sum := by
  simpa only [physicalProfile22, nsmul_eq_mul] using productionPacketProfile_table_sum f

theorem physicalProfile22_live_sum (f : (Fin 5 → Nat) → ℝ) :
    (∑ q : {q : State22 // q ≠ 0}, f (physicalProfile22 q.val)) =
      (productionPacketProfileTable.map (fun e => (e.1 : ℝ) * f e.2)).sum -
        f ![16,0,0,0,0] := by
  classical
  have hs := Fintype.sum_subtype_add_sum_subtype (fun q : State22 => q ≠ 0)
    (fun q => f (physicalProfile22 q))
  have hz : (∑ q : {q : State22 // ¬ q ≠ 0}, f (physicalProfile22 q.val)) =
      f ![16,0,0,0,0] := by
    have he (q : {q : State22 // ¬ q ≠ 0}) : q.val = 0 := not_not.mp q.property
    simp only [he, physicalProfile22_zero, sum_const, card_univ]
    have hc : Fintype.card {q : State22 // ¬ q ≠ 0} = 1 := by
      simp only [not_not]
      exact Fintype.card_subtype_eq 0
    rw [hc, one_nsmul]
  rw [hz, physicalProfile22_full_sum] at hs
  linarith

def physicalProfileLiveMean22 (f : (Fin 5 → Nat) → ℝ) : ℝ :=
  ((productionPacketProfileTable.map (fun e => (e.1 : ℝ) * f e.2)).sum -
    f ![16,0,0,0,0]) / (liveCount State22 : ℝ)

theorem physicalProfile22_live_mean (f : (Fin 5 → Nat) → ℝ) :
    (∑ q : {q : State22 // q ≠ 0}, f (physicalProfile22 q.val)) /
      (liveCount State22 : ℝ) = physicalProfileLiveMean22 f := by
  rw [physicalProfile22_live_sum]
  rfl

theorem physicalProfile22_class_sum (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) :
    (∑ q : State22, if physicalClassAllows22 i q then f (physicalProfile22 q) else 0) =
      (productionPacketProfileTable.map (fun e => if physicalProfileClassAllows22 i e.2 then
        (e.1 : ℝ) * f e.2 else 0)).sum := by
  classical
  have hs := physicalProfile22_full_sum
    (fun p => if physicalProfileClassAllows22 i p then f p else 0)
  simpa only [physicalProfile22_class, mul_ite, mul_zero] using hs

def physicalProfileClassCount22 (i : Fin 11) : Nat :=
  (productionPacketProfileTable.map (fun e => if physicalProfileClassAllows22 i e.2 then e.1 else 0)).sum

theorem physicalProfileClassCount22_checked : ∀ i : Fin 11,
    physicalProfileClassCount22 i =
      ([1,4194303,4194303,2604,291648,888832,1828134,888832,291648,2604,1] : List Nat).getD i.val 0 := by
  decide

theorem physicalProfile22_class_count (i : Fin 11) :
    (univ.filter (fun q : State22 => physicalClassAllows22 i q)).card = physicalProfileClassCount22 i := by
  have hs := physicalProfile22_class_sum i (fun _ => 1)
  simp only [mul_one] at hs
  apply Nat.cast_injective (R := ℝ)
  unfold physicalProfileClassCount22
  rw [card_filter, Nat.cast_sum, Nat.cast_list_sum, List.map_map]
  simpa only [Function.comp_def, apply_ite, Nat.cast_zero, Nat.cast_one] using hs

def physicalProfileClassMean22 (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) : ℝ :=
  (productionPacketProfileTable.map (fun e => if physicalProfileClassAllows22 i e.2 then
    (e.1 : ℝ) * f e.2 else 0)).sum / (physicalProfileClassCount22 i : ℝ)

theorem physicalProfile22_class_mean (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) :
    (∑ q : State22, if physicalClassAllows22 i q then f (physicalProfile22 q) else 0) /
      ((univ.filter (fun q : State22 => physicalClassAllows22 i q)).card : ℝ) =
        physicalProfileClassMean22 i f := by
  rw [physicalProfile22_class_sum, physicalProfile22_class_count]
  rfl

theorem physicalProfile22_class_cap (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) (B : ℝ)
    (h : ∀ e ∈ productionPacketProfileTable, physicalProfileClassAllows22 i e.2 → f e.2 ≤ B)
    (q : State22) (hq : physicalClassAllows22 i q) : f (physicalProfile22 q) ≤ B := by
  obtain ⟨e, he, heq⟩ := physicalProfile22_mem_table q
  rw [← heq]
  exact h e he (by rw [heq, physicalProfile22_class]; exact hq)

def physicalProfileClassSet22 (i : Fin 11) : Finset (Nat × (Fin 5 → Nat)) := by
  classical
  exact productionPacketProfileTable.toFinset.filter (fun e => physicalProfileClassAllows22 i e.2)

theorem physicalProfileClassSet22_nonempty (i : Fin 11) : (physicalProfileClassSet22 i).Nonempty := by
  have hh : ∀ k : Fin 11, ∃ e ∈ productionPacketProfileTable, physicalProfileClassAllows22 k e.2 := by
    decide
  obtain ⟨e, he, ha⟩ := hh i
  refine ⟨e, ?_⟩
  classical
  simpa only [physicalProfileClassSet22, mem_filter, List.mem_toFinset] using And.intro he ha

def physicalProfileClassMax22 (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) : ℝ :=
  (physicalProfileClassSet22 i).sup' (physicalProfileClassSet22_nonempty i) (fun e => f e.2)

theorem physicalProfile22_le_class_max (i : Fin 11) (f : (Fin 5 → Nat) → ℝ)
    (q : State22) (hq : physicalClassAllows22 i q) :
    f (physicalProfile22 q) ≤ physicalProfileClassMax22 i f := by
  classical
  obtain ⟨e, he, heq⟩ := physicalProfile22_mem_table q
  have hm : e ∈ physicalProfileClassSet22 i := by
    simp only [physicalProfileClassSet22, mem_filter, List.mem_toFinset]
    exact ⟨he, by rw [heq, physicalProfile22_class]; exact hq⟩
  rw [← heq]
  exact le_sup' (fun e => f e.2) hm

theorem physicalProfileClassMax22_le (i : Fin 11) (f : (Fin 5 → Nat) → ℝ) (B : ℝ)
    (h : ∀ e ∈ productionPacketProfileTable, physicalProfileClassAllows22 i e.2 → f e.2 ≤ B) :
    physicalProfileClassMax22 i f ≤ B := by
  classical
  apply sup'_le
  intro e he
  have hh : e ∈ productionPacketProfileTable ∧ physicalProfileClassAllows22 i e.2 := by
    simpa only [physicalProfileClassSet22, mem_filter, List.mem_toFinset] using he
  exact h e hh.1 hh.2

#print axioms physicalProfile22_product
#print axioms physicalProfile22_class_sum
#print axioms physicalProfile22_class_cap
end Spin.Turbo
