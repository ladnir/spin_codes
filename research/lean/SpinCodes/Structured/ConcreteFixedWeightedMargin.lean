import SpinCodes.Structured.ConcretePlacementLimitActual
import SpinCodes.Structured.WeightedNorm

/-! Strict continuum weighted margins transfer to actual finite region kernels.
A fixed per-region margin suffices, even for a growing number of regions. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder ConcreteMaps FiniteKernel Filter
open scoped Topology
attribute [local instance] Classical.propDecidable

def stateWeight (v : ℝ) (q : State) : ℝ := if q = ∅ then 1 else v

theorem stateWeight_pos {v : ℝ} (hv : 0 < v) (q : State) : 0 < stateWeight v q := by
  unfold stateWeight
  split_ifs <;> simp_all

theorem stateWeight_sum (v : ℝ) : (∑ q : State, stateWeight v q) = 1 + 524287 * v := by
  rw [Spin.sum_split]
  have he (q : State) (hq : q ∈ Spin.nonzeroStates 19) : stateWeight v q = v := by
    simp only [stateWeight, Spin.ne_empty_of_mem_nonzeroStates hq, ite_false]
  rw [Finset.sum_congr rfl he, sum_const, nsmul_eq_mul, nonzeroStates_card_actual]
  simp [stateWeight]

theorem lifted_weight_action (K : Matrix (Fin 2) (Fin 2) ℝ) (v : ℝ) (q : State) :
    (∑ r, (liveLift * K * liveProjection) q r * stateWeight v r) =
      K (stateClass q) 0 + K (stateClass q) 1 * v := by
  rw [Spin.sum_split]
  have he (r : State) (hr : r ∈ Spin.nonzeroStates 19) :
      (liveLift * K * liveProjection) q r * stateWeight v r = K (stateClass q) 1 / 524287 * v := by
    simp only [lifted_entry, stateWeight, Spin.ne_empty_of_mem_nonzeroStates hr, ite_false]
  rw [Finset.sum_congr rfl he, sum_const, nsmul_eq_mul, nonzeroStates_card_actual]
  simp only [lifted_entry, stateWeight, ite_true, mul_one]
  push_cast
  ring

theorem lifted_weight_bound {v c : ℝ} (hv : 0 ≤ v) {K : Matrix (Fin 2) (Fin 2) ℝ}
    (hK : Spin.RowNormLe v c K) (q : State) :
    (∑ r, (liveLift * K * liveProjection) q r * stateWeight v r) ≤ c * stateWeight v q := by
  rw [lifted_weight_action]
  have h0 := add_le_add (le_abs_self (K 0 0)) (mul_le_mul_of_nonneg_right (le_abs_self (K 0 1)) hv)
  have h1 := add_le_add (le_abs_self (K 1 0)) (mul_le_mul_of_nonneg_right (le_abs_self (K 1 1)) hv)
  by_cases hq : q = ∅
  · simpa only [stateClass, stateWeight, hq, ite_true, mul_one] using h0.trans hK.1
  · simpa only [stateClass, stateWeight, hq, ite_false] using h1.trans hK.2

theorem actual_region_weighted_margin (a : Nat) {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hδ : 0 < δ) (hK : Spin.RowNormLe v c (continuumRegionKernel θ a)) :
    ∀ᶠ R : Nat in atTop, ∀ S : Finset (Fin (128 * R)), S.card = a → ∀ q : State,
      (∑ r, shuffledRegionKernel θ S q r * stateWeight v r) ≤ (c + δ) * stateWeight v q := by
  let m := min 1 v
  let W := 1 + 524287 * v
  have hm : 0 < m := lt_min (by norm_num) hv
  have hW : 0 < W := by dsimp [W]; positivity
  have hε : 0 < δ * m / W := div_pos (mul_pos hδ hm) hW
  filter_upwards [shuffledRegionKernel_uniform_limit a hθ hε] with R hR
  intro S hS q
  have hh : (∑ r, shuffledRegionKernel θ S q r * stateWeight v r) ≤
      (∑ r, (liveLift * continuumRegionKernel θ a * liveProjection) q r * stateWeight v r) + δ * m := by
    calc
      _ ≤ ∑ r, ((liveLift * continuumRegionKernel θ a * liveProjection) q r + δ * m / W) * stateWeight v r := by
        apply sum_le_sum
        intro r _
        apply mul_le_mul_of_nonneg_right _ (stateWeight_pos hv r).le
        have he := (abs_lt.mp (hR S hS q r)).2
        linarith
      _ = _ := by
        simp only [add_mul, sum_add_distrib, ← mul_sum, stateWeight_sum]
        change _ + δ * m / W * W = _
        rw [div_mul_cancel₀ _ hW.ne']
  have hmq : m ≤ stateWeight v q := by
    by_cases hq : q = ∅
    · simpa only [stateWeight, hq, ite_true] using min_le_left (1 : ℝ) v
    · simpa only [stateWeight, hq, ite_false] using min_le_right (1 : ℝ) v
  have hb := lifted_weight_bound hv.le hK q
  have hmδ := mul_le_mul_of_nonneg_left hmq hδ.le
  nlinarith

end Spin.Structured.Placement

