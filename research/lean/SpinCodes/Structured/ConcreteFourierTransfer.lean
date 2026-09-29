import SpinCodes.Structured.ConcreteFourierRow
import SpinCodes.Structured.ConcreteOccupation

noncomputable section
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteMaps ConcreteEncoder ConcreteScalar Spin.Imt

def atomBudget (c : ℝ) : Coords 5 :=
  ⟨c, 0, fun i => Occupation.Sparse.count i * c⟩

def diffuseCap (β z : ℝ) : ℝ :=
  univ.sup' univ_nonempty (fun i : Fin 5 => targetCap β z (shellWeight i))

/-- The Fourier alternative retains the exact zero row and bounds each live target atom. -/
def matrix (β z : ℝ) : Transfer 5 :=
  ⟨(Occupation.Sparse.numericalMatrix β z).rowZ,
    atomBudget (diffuseCap β z), fun i => atomBudget (targetCap β z (shellWeight i))⟩

theorem atomBudget_nonneg {c : ℝ} (hc : 0 ≤ c) : (atomBudget c).Nonneg :=
  ⟨hc, le_rfl, fun i => mul_nonneg (Occupation.Sparse.count_pos i).le hc⟩

theorem atomBudget_dominates {μ : State → ℝ} {c : ℝ} (h : ∀ q, μ q ≤ c) :
    LiveDominates actualShellSystem (atomBudget c) μ := by
  classical
  rw [liveDominates_iff_witness]
  refine ⟨fun _ => 0, fun _ => le_rfl, rfl, by simp [atomBudget], ?_⟩
  intro q
  change μ q ≤ (if q = ∅ then c else 0) + 0 +
    ∑ i, if q ∈ actualShellSystem.shell i then
      (Occupation.Sparse.count i*c)/(actualShellSystem.shell i).card else 0
  have he : (∑ i, if q ∈ actualShellSystem.shell i then
      (Occupation.Sparse.count i*c)/(actualShellSystem.shell i).card else 0) =
      if q = ∅ then 0 else c := by
    rw [← actualShellSystem.sum_mem q c]
    apply sum_congr rfl
    intro i _
    rw [actualShellSystem_card]
    simp only [mul_div_cancel_left₀ c (Occupation.Sparse.count_pos i).ne']
  rw [he]
  split_ifs <;> simpa only [add_zero, zero_add] using h q

theorem targetCap_nonneg {β z : ℝ} (hb0 : 0 ≤ β) (hb1 : β ≤ 1) (hz : 0 ≤ z) (d : ℕ) :
    0 ≤ targetCap β z d := by
  have hs : 0 ≤ spectrumCap d (ratio0 β z) (ratio1 β z) := by
    unfold spectrumCap
    apply div_nonneg _ (by norm_num)
    apply sum_nonneg
    intro w _
    apply mul_nonneg (Nat.cast_nonneg _) _
    unfold overlapCap
    exact le_trans (mul_nonneg (pow_nonneg (abs_nonneg _) _) (pow_nonneg (abs_nonneg _) _))
      (le_max_left _ _)
  have hb : 0 ≤ 1-β := sub_nonneg.mpr hb1
  unfold targetCap scalarF
  positivity

theorem diffuseCap_nonneg {β z : ℝ} (hb0 : 0 ≤ β) (hb1 : β ≤ 1) (hz : 0 ≤ z) :
    0 ≤ diffuseCap β z :=
  (targetCap_nonneg hb0 hb1 hz (shellWeight 0)).trans
    (le_sup' (fun i : Fin 5 => targetCap β z (shellWeight i)) (mem_univ 0))

theorem matrix_nonneg {β z : ℝ} (hb0 : 0 ≤ β) (hb1 : β ≤ 1) (hz : 0 ≤ z) :
    (matrix β z).Nonneg :=
  ⟨(Occupation.Sparse.numericalMatrix_nonneg hb0 hb1 hz).1,
    atomBudget_nonneg (diffuseCap_nonneg hb0 hb1 hz),
    fun i => atomBudget_nonneg (targetCap_nonneg hb0 hb1 hz _)⟩

theorem rowZ_dominates {β z : ℝ} (hb0 : 0 ≤ β) (hb1 : β ≤ 1) (hz : 0 ≤ z) :
    LiveDominates actualShellSystem (matrix β z).rowZ (bernoulliRow β z ∅) := by
  have h := liveDominates_sum (fun j : Fin 129 => liveDominates_smul
    (stepRow_dominates_Z hz j) (Occupation.probability_nonneg 128 j hb0 hb1))
  have he : bernoulliRow β z ∅ = fun r =>
      ∑ j : Fin 129, Occupation.probability 128 j β * stepRow j z ∅ r :=
    funext (bernoulliRow_eq_mixture β z ∅)
  rw [he]
  exact h

theorem rowD_dominates {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z)
    {q : State} (hq : q ≠ ∅) :
    LiveDominates actualShellSystem (matrix β z).rowD (bernoulliRow β z q) := by
  apply atomBudget_dominates
  intro r
  refine (bernoulliRow_target_le hb0 hb1 hz hq r).trans ?_
  have hmem : q ∈ (univ : Finset (Fin 5)).biUnion weightShell := by
    rw [weightShell_covers actual_shell_counts]
    exact Spin.mem_nonzeroStates hq
  obtain ⟨i, _, hi⟩ := mem_biUnion.mp hmem
  rw [weightShell_weight hi]
  exact le_sup' (fun i : Fin 5 => targetCap β z (shellWeight i)) (mem_univ i)

theorem rowS_dominates {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z) (i : Fin 5) :
    LiveDominates actualShellSystem ((matrix β z).rowS i)
      (fun r => (∑ q ∈ actualShellSystem.shell i, bernoulliRow β z q r) /
        (actualShellSystem.shell i).card) := by
  apply atomBudget_dominates
  intro r
  apply (div_le_iff₀ (actualShellSystem.card_pos i)).mpr
  calc
    _ ≤ ∑ q ∈ actualShellSystem.shell i, targetCap β z (shellWeight i) := by
      apply sum_le_sum
      intro q hq
      have h := bernoulliRow_target_le hb0 hb1 hz (actualShellSystem.ne_empty_of_mem hq) r
      simpa only [weightShell_weight hq] using h
    _ = _ := by simp only [sum_const, nsmul_eq_mul]; ring

theorem step_liveDominates {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z)
    (c : Coords 5) (μ : State → ℝ) (hc : c.Nonneg) (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem ((matrix β z).apply c) (bernoulliStep β z μ) := by
  apply kernelApply_liveDominates actualShellSystem _ (bernoulliRow β z) _
    (matrix_nonneg hb0.le hb1.le hz.le).2.1 (rowZ_dominates hb0.le hb1.le hz.le)
    (fun _ hq => rowD_dominates hb0 hb1 hz hq) (rowS_dominates hb0 hb1 hz) c μ hc hμ
  intro q r
  rw [bernoulliRow_eq_mixture]
  exact sum_nonneg fun j _ => mul_nonneg (Occupation.probability_nonneg 128 j hb0.le hb1.le)
    (stepRow_nonneg hz.le j q r)

theorem moment_bound {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z) (R : ℕ) :
    ∑ q, ((bernoulliStep β z)^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      (((matrix β z).apply)^[R] (Coords.eZ 5)).total :=
  live_moment_eZ _ (matrix_nonneg hb0.le hb1.le hz.le) _ (step_liveDominates hb0 hb1 hz) R

end Spin.Structured.ConcreteFourier
