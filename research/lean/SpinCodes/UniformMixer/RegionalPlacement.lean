import SpinCodes.UniformMixer.RegionalMatrix
import SpinCodes.Structured.ConcretePlacementGaps

/-!
# Bad placements in a shuffled region

The actual regional shuffle chooses a fixed-cardinality support uniformly.
We discard placements touching the first or last 64-bit step or containing
two active bits at distance at most 128. This is stronger than the paper's
minimum separation condition, but its probability remains `O(a²/L)` for
fixed occupation `a`.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement
open Spin.Structured.Routing
attribute [local instance] Classical.propDecidable

def badRegionalSites {n : ℕ} (T : Finset (Fin n)) : Prop :=
  (∃ x ∈ T, x.val < 64 ∨ n ≤ x.val + 64) ∨
  (∃ x ∈ T, ∃ y ∈ T,
    x ≠ y ∧ Nat.dist x.val y.val ≤ 128)

theorem shuffle_prob_badRegionalSites {n : ℕ} (hn : 2 ≤ n)
    (S : Finset (Fin n)) :
    (shuffleLaw S).prob badRegionalSites ≤
      128 * (S.card : ℝ) / n +
        257 * (S.card : ℝ) ^ 2 / ((n : ℝ) - 1) := by
  have hb := shuffle_prob_boundary S 64
  have hc := shuffle_prob_close hn S 128
  calc
    (shuffleLaw S).prob badRegionalSites ≤
      (shuffleLaw S).prob
        (fun T => ∃ x ∈ T, x.val < 64 ∨ n ≤ x.val + 64) +
      (shuffleLaw S).prob
        (fun T => ∃ x ∈ T, ∃ y ∈ T,
          x ≠ y ∧ Nat.dist x.val y.val ≤ 128) := by
            have heq : (shuffleLaw S).prob badRegionalSites =
                (shuffleLaw S).prob
                  (fun T => (∃ x ∈ T, x.val < 64 ∨ n ≤ x.val + 64) ∨
                    (∃ x ∈ T, ∃ y ∈ T,
                      x ≠ y ∧ Nat.dist x.val y.val ≤ 128)) := by
              apply Spin.FinPMF.prob_congr
              intro T
              rfl
            rw [heq]
            exact Spin.FinPMF.prob_or_le _ _ _
    _ ≤ 2 * (64 : ℝ) * S.card / n +
          (2 * (128 : ℝ) + 1) * (S.card : ℝ) ^ 2 /
            ((n : ℝ) - 1) := add_le_add hb hc
    _ = 128 * (S.card : ℝ) / n +
          257 * (S.card : ℝ) ^ 2 / ((n : ℝ) - 1) := by ring

theorem bit_distance_of_step_distance (x y : ℕ)
    (h : Nat.dist (x / 64) (y / 64) ≤ 1) :
    Nat.dist x y ≤ 127 := by
  have hx := Nat.mod_lt x (by omega : 0 < 64)
  have hy := Nat.mod_lt y (by omega : 0 < 64)
  have hxx := Nat.div_add_mod x 64
  have hyy := Nat.div_add_mod y 64
  unfold Nat.dist at *
  omega

theorem goodRegionalSites_step_separation {R : ℕ}
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (x y : Fin (64 * R)) (hx : x ∈ T) (hy : y ∈ T) (hxy : x ≠ y) :
    2 ≤ Nat.dist (x.val / 64) (y.val / 64) := by
  by_contra h
  have hsteps : Nat.dist (x.val / 64) (y.val / 64) ≤ 1 := by omega
  have hbits := bit_distance_of_step_distance x.val y.val hsteps
  exact hgood (Or.inr ⟨x, hx, y, hy, hxy, by omega⟩)

theorem goodRegionalSites_interior {R : ℕ}
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (x : Fin (64 * R)) (hx : x ∈ T) :
    1 ≤ x.val / 64 ∧ x.val / 64 + 1 < R := by
  have hleft : 64 ≤ x.val := by
    by_contra h
    exact hgood (Or.inl ⟨x, hx, Or.inl (by omega)⟩)
  have hright : x.val + 64 < 64 * R := by
    by_contra h
    exact hgood (Or.inl ⟨x, hx, Or.inr (by omega)⟩)
  have hxmod := Nat.mod_lt x.val (by omega : 0 < 64)
  have hxdiv := Nat.div_add_mod x.val 64
  constructor <;> omega

end Spin.UniformMixer
