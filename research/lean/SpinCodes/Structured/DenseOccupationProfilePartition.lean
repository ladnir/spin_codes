import SpinCodes.Structured.DenseOccupationProfileFailure

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter

/-- Only occupied positions are indexed, so there are exactly (b+1)^Q profiles. -/
abbrev WeightProfile {L : ℕ} (S : Finset (Fin L)) (b : ℕ) := S → Fin (b+1)

def messageProfile {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (x : Fin L → LocalMessage k) : WeightProfile S (k*24) :=
  fun i => ⟨wtF (encode seed (x i)), by have := wtF_le (encode seed (x i)); omega⟩

def profileWeights {L b : ℕ} (S : Finset (Fin L)) (w : WeightProfile S b) : Fin L → ℕ :=
  fun i => if h : i ∈ S then (w ⟨i,h⟩).val else 0

theorem profile_fiber_eq {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (w : WeightProfile S (k*24)) :
    (activeMessages S).filter (fun x => messageProfile seed S x = w) =
      profileMessages seed S (profileWeights S w) := by
  classical
  ext x
  simp only [Finset.mem_filter, activeMessages, Finset.mem_univ, true_and, profileMessages]
  constructor
  · rintro ⟨ha,hp⟩ i
    by_cases hi : i ∈ S
    · simp only [hi, ite_true, profileWeights, dif_pos hi]
      refine ⟨(ha i).mpr hi, ?_⟩
      exact congrArg Fin.val (congrFun hp ⟨i,hi⟩)
    · simp only [hi, ite_false]
      exact not_not.mp (fun h => hi ((ha i).mp h))
  · intro h
    constructor
    · intro i
      by_cases hi : i ∈ S
      · have hh := h i
        simp only [hi, ite_true] at hh
        exact ⟨fun _ => hi, fun _ => hh.1⟩
      · have hh := h i
        simp only [hi, ite_false] at hh
        simp [hi, hh]
    · funext i
      apply Fin.ext
      have hh := h i.val
      simp only [i.property, ite_true] at hh
      simpa only [profileWeights, dif_pos i.property, messageProfile] using hh.2

theorem sum_active_eq_profiles {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (f : (Fin L → LocalMessage k) → ℝ) :
    (∑ x ∈ activeMessages S, f x) =
      ∑ w : WeightProfile S (k*24), ∑ x ∈ profileMessages seed S (profileWeights S w), f x := by
  rw [← Finset.sum_fiberwise (activeMessages S) (messageProfile seed S) f]
  simp only [profile_fiber_eq]

/-- Uniform profile bounds sum with the exact finite weight-profile cost. -/
theorem sum_active_le_profile_bound {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (f : (Fin L → LocalMessage k) → ℝ) {B : ℝ}
    (h : ∀ w : WeightProfile S (k*24),
      (∑ x ∈ profileMessages seed S (profileWeights S w), f x) ≤ B) :
    (∑ x ∈ activeMessages S, f x) ≤ (((k*24:ℕ):ℝ)+1)^S.card * B := by
  rw [sum_active_eq_profiles seed S f]
  calc
    _ ≤ ∑ _w : WeightProfile S (k*24), B := Finset.sum_le_sum (fun w _ => h w)
    _ = _ := by simp [WeightProfile, Fintype.card_fun]

#print axioms profile_fiber_eq
#print axioms sum_active_eq_profiles
#print axioms sum_active_le_profile_bound

/-- Exact active-position and profile costs for a whole occupation layer. -/
theorem sum_occupation_le_profile_bound {L k : ℕ} (seed : Seed k) (Q : ℕ)
    (f : (Fin L → LocalMessage k) → ℝ) {B : ℝ}
    (h : ∀ S : Finset (Fin L), S.card = Q → ∀ w : WeightProfile S (k*24),
      (∑ x ∈ profileMessages seed S (profileWeights S w), f x) ≤ B) :
    (∑ x ∈ univ.filter (fun x => occupation x = Q), f x) ≤
      (L.choose Q:ℝ)*((((k*24:ℕ):ℝ)+1)^Q*B) := by
  rw [sum_occupation_eq_activeSets]
  calc
    _ ≤ ∑ _S ∈ (univ : Finset (Fin L)).powersetCard Q,
        ((((k*24:ℕ):ℝ)+1)^Q*B) := by
      apply Finset.sum_le_sum
      intro S hS
      have hc : S.card = Q := (Finset.mem_powersetCard.mp hS).2
      simpa only [hc] using sum_active_le_profile_bound seed S f (h S hc)
    _ = _ := by simp only [Finset.sum_const, nsmul_eq_mul, Finset.card_powersetCard,
        Finset.card_univ, Fintype.card_fin]

#print axioms sum_occupation_le_profile_bound
end Spin.Structured.DenseOccupationFixed
