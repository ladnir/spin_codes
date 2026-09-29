import SpinCodes.Structured.ConcreteTransfer

/-! The binomial mixture is the transition with independent Bernoulli input
bits, the concrete A/C maps, and the stated random transvection law. -/
noncomputable section
namespace Spin.Imt.Occupation

theorem matrix_apply {k : Nat} (n : Nat) (count : Fin k → ℝ) (M β : ℝ)
    (r : Fin (n + 1) → RowData k) (c : Coords k) :
    (matrix n count M β r).apply c = Coords.sum fun j : Fin (n + 1) =>
      Coords.smul (probability n j β) ((fixed count M (r j)).apply c) := by
  apply Coords.ext <;>
    simp only [matrix, Transfer.apply, Coords.sum, Coords.smul,
      Finset.mul_sum, Finset.sum_add_distrib, mul_add]
  all_goals try funext h
  all_goals
    try simp only [Finset.mul_sum, Finset.sum_add_distrib, mul_add]
    rw [Finset.sum_comm (f := fun i j => _)]
    simp only [mul_assoc, mul_left_comm, mul_comm]

end Spin.Imt.Occupation

namespace Spin.Structured.ConcreteMaps
open Spin.Imt
open scoped symmDiff

def mixtureStep (β z : ℝ) (μ : Finset (Fin 19) → ℝ) (r : Finset (Fin 19)) : ℝ :=
  ∑ j : Fin 129, Occupation.probability 128 j β * fixedStep j z μ r

def bernoulliRow (β z : ℝ) (q r : Finset (Fin 19)) : ℝ :=
  ∑ x : Finset (Fin 128),
    (β ^ x.card * (1 - β) ^ (128 - x.card)) * emitted z q x * actLaw q (r ∆ Cset x)

def bernoulliStep (β z : ℝ) := kernelApply (bernoulliRow β z)

theorem sum_layers (f : Finset (Fin 128) → ℝ) :
    (∑ x, f x) = ∑ j : Fin 129, ∑ x ∈ Spin.layer 128 j, f x := by
  let w (x : Finset (Fin 128)) : Fin 129 :=
    ⟨x.card, Nat.lt_succ_of_le (by simpa using Finset.card_le_univ x)⟩
  have he (j : Fin 129) : Finset.univ.filter (fun x => w x = j) = Spin.layer 128 j := by
    ext x
    simp [Spin.layer, w, Fin.ext_iff]
  rw [← Finset.sum_fiberwise Finset.univ w f]
  simp only [he]

theorem bernoulliRow_eq_mixture (β z : ℝ) (q r : Finset (Fin 19)) :
    bernoulliRow β z q r = ∑ j : Fin 129, Occupation.probability 128 j β * stepRow j z q r := by
  rw [bernoulliRow, sum_layers]
  apply Finset.sum_congr rfl
  intro j _
  have hn : (Nat.choose 128 j : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (Nat.le_of_lt_succ j.isLt)).ne'
  have he : (∑ x ∈ Spin.layer 128 j,
      (β ^ x.card * (1 - β) ^ (128 - x.card)) * emitted z q x * actLaw q (r ∆ Cset x)) =
      (β ^ j.val * (1 - β) ^ (128 - j.val)) *
        (∑ x ∈ Spin.layer 128 j, emitted z q x * actLaw q (r ∆ Cset x)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro x hx
    rw [(Finset.mem_powersetCard.mp hx).2]
    ring
  rw [he, Occupation.probability, stepRow]
  field_simp

theorem bernoulliStep_eq_mixture (β z : ℝ) (μ : Finset (Fin 19) → ℝ) :
    bernoulliStep β z μ = mixtureStep β z μ := by
  funext r
  simp only [bernoulliStep, mixtureStep, fixedStep, kernelApply,
    bernoulliRow_eq_mixture, Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro j _
  apply Finset.sum_congr rfl
  intro q _
  ring

end Spin.Structured.ConcreteMaps
