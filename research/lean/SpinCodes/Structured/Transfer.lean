/-
The weighted-state abstraction of `app:imt-finite-transfers`.

    "We represent a weighted state measure by seven nonnegative coordinates
     (Z, D, S_1, ..., S_5).  The Z coordinate is mass at zero.  The D
     coordinate bounds arbitrary mass on nonzero states.  Coordinate S_i
     dominates S_i times the uniform probability measure on {q : wt(Aq)=d_i}.
     Terminal mass is bounded by the sum of coordinates."

The `D` coordinate is a *total* budget spread arbitrarily, not a pointwise
bound, so domination is an existential over that spread.  Everything else is
pointwise.

The shells are kept abstract: any family partitioning the nonzero states will
do, and the paper's five weight classes of `A` are one instance.  Nothing here
depends on there being five of them, or on `A` at all.

What is proved here is the architecture — the terminal bound, monotonicity in
the coordinates, and the exact shell decomposition of one transvection step.
The last of these is where `Transvection.lean` is consumed: the row of a
transfer matrix leaving a nonzero state is *exactly* `½` at that state plus
`c_h/(2M)` on each shell.
-/
import SpinCodes.Structured.Transvection

namespace Spin.Imt

open Finset
open scoped symmDiff

/-! ## Shells -/

/-- A partition of the nonzero states into shells. -/
structure ShellSystem (s k : ℕ) where
  shell : Fin k → Finset (Finset (Fin s))
  nonempty : ∀ i, (shell i).Nonempty
  disjoint : ∀ i j, i ≠ j → Disjoint (shell i) (shell j)
  covers : (univ : Finset (Fin k)).biUnion shell = nonzeroStates s

variable {s k : ℕ}

lemma ShellSystem.card_pos (sys : ShellSystem s k) (i : Fin k) :
    (0 : ℝ) < ((sys.shell i).card : ℝ) := by
  have : 0 < (sys.shell i).card := Finset.card_pos.mpr (sys.nonempty i)
  exact_mod_cast this

lemma ShellSystem.ne_empty_of_mem (sys : ShellSystem s k) {i : Fin k}
    {w : Finset (Fin s)} (hw : w ∈ sys.shell i) : w ≠ ∅ := by
  have : w ∈ nonzeroStates s := by
    rw [← sys.covers]
    exact Finset.mem_biUnion.mpr ⟨i, Finset.mem_univ i, hw⟩
  exact ne_empty_of_mem_nonzeroStates this

/-- Each nonzero state lies in exactly one shell; zero lies in none. -/
theorem ShellSystem.sum_mem (sys : ShellSystem s k) (w : Finset (Fin s)) (c : ℝ) :
    ∑ i : Fin k, (if w ∈ sys.shell i then c else 0) = if w = ∅ then 0 else c := by
  classical
  by_cases hw : w = ∅
  · rw [if_pos hw]
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [if_neg]
    intro hmem
    exact sys.ne_empty_of_mem hmem hw
  · rw [if_neg hw]
    have hmem : w ∈ (univ : Finset (Fin k)).biUnion sys.shell := by
      rw [sys.covers]
      exact mem_nonzeroStates hw
    obtain ⟨i₀, -, hi₀⟩ := Finset.mem_biUnion.mp hmem
    rw [Finset.sum_eq_single i₀]
    · rw [if_pos hi₀]
    · intro i _ hne
      rw [if_neg]
      intro hmem'
      exact (Finset.disjoint_left.mp (sys.disjoint i i₀ hne) hmem') hi₀
    · intro hcon
      exact absurd (Finset.mem_univ i₀) hcon

/-! ## The coordinates -/

/-- The weighted-state coordinates: mass at zero, a free budget on nonzero
states, and one coordinate per shell. -/
structure Coords (k : ℕ) where
  Z : ℝ
  D : ℝ
  S : Fin k → ℝ

/-- Terminal mass is bounded by the sum of coordinates. -/
def Coords.total (c : Coords k) : ℝ := c.Z + c.D + ∑ i, c.S i

/-- `c` dominates the weighted state measure `μ`.  The `D` budget is spread by
an arbitrary nonnegative `ν`; the rest is pointwise. -/
def Dominates (sys : ShellSystem s k) (c : Coords k)
    (μ : Finset (Fin s) → ℝ) : Prop :=
  ∃ ν : Finset (Fin s) → ℝ,
    (∀ q, 0 ≤ ν q) ∧ (∑ q, ν q ≤ c.D) ∧
    ∀ q, μ q ≤ (if q = ∅ then c.Z else 0) + ν q
        + ∑ i, (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0)

/-- **Terminal mass is bounded by the sum of coordinates.** -/
theorem Dominates.total_le {sys : ShellSystem s k} {c : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : Dominates sys c μ) :
    ∑ q, μ q ≤ c.total := by
  classical
  obtain ⟨ν, hν0, hνD, hbd⟩ := h
  have hstep : ∑ q : Finset (Fin s), μ q
      ≤ ∑ q : Finset (Fin s),
          ((if q = ∅ then c.Z else 0) + ν q
            + ∑ i, (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0)) :=
    Finset.sum_le_sum fun q _ => hbd q
  refine hstep.trans ?_
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
  have h1 : ∑ q : Finset (Fin s), (if q = ∅ then c.Z else 0) = c.Z := by
    rw [Finset.sum_ite_eq' univ (∅ : Finset (Fin s)) (fun _ => c.Z)]
    simp
  have h3 : ∑ q : Finset (Fin s),
        ∑ i, (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0)
      = ∑ i, c.S i := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hne := (sys.card_pos i).ne'
    rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.sum_const, nsmul_eq_mul]
    field_simp
  rw [h1, h3, Coords.total]
  linarith

/-- Domination is monotone in the coordinates. -/
theorem Dominates.mono {sys : ShellSystem s k} {c c' : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : Dominates sys c μ)
    (hZ : c.Z ≤ c'.Z) (hD : c.D ≤ c'.D) (hS : ∀ i, c.S i ≤ c'.S i) :
    Dominates sys c' μ := by
  classical
  obtain ⟨ν, hν0, hνD, hbd⟩ := h
  refine ⟨ν, hν0, hνD.trans hD, fun q => (hbd q).trans ?_⟩
  have h1 : (if q = ∅ then c.Z else 0) ≤ (if q = ∅ then c'.Z else 0) := by
    by_cases hq : q = ∅ <;> simp [hq, hZ]
  have h3 : ∀ i : Fin k,
      (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0)
        ≤ (if q ∈ sys.shell i then c'.S i / ((sys.shell i).card : ℝ) else 0) := by
    intro i
    by_cases hq : q ∈ sys.shell i
    · simp only [hq, if_true]
      exact (div_le_div_iff_of_pos_right (sys.card_pos i)).mpr (hS i)
    · simp [hq]
  have h3' : ∑ i, (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0)
      ≤ ∑ i, (if q ∈ sys.shell i then c'.S i / ((sys.shell i).card : ℝ) else 0) :=
    Finset.sum_le_sum fun i _ => h3 i
  linarith

/-! ## One transvection step out of a nonzero state

The row is *exact*: `½` stays at `q`, and the other half spreads flat over the
nonzero states, which in shell coordinates is `c_h / (2M)` on shell `h`. -/

/-- The transvection image of a point mass at a nonzero state. -/
noncomputable def actLaw (q : Finset (Fin s)) (w : Finset (Fin s)) : ℝ :=
  (((transPairs s).filter (fun p => act p.1 p.2 q = w)).card : ℝ)
    / ((transPairs s).card : ℝ)

/-- **The row leaving a nonzero state**, in shell coordinates. -/
theorem dominates_actLaw (sys : ShellSystem s k) {q : Finset (Fin s)} (hq : q ≠ ∅) :
    Dominates sys
      ⟨0, 1 / 2,
        fun h => ((sys.shell h).card : ℝ) / (2 * ((nonzeroStates s).card : ℝ))⟩
      (actLaw q) := by
  classical
  refine ⟨fun w => if w = q then 1 / 2 else 0, ?_, ?_, ?_⟩
  · intro w
    by_cases hw : w = q <;> simp [hw]
  · rw [Finset.sum_ite_eq' univ q (fun _ => (1 : ℝ) / 2)]
    simp
  · intro w
    have hrow := prob_act hq w
    have hshell : ∑ i : Fin k,
        (if w ∈ sys.shell i then
          (((sys.shell i).card : ℝ) / (2 * ((nonzeroStates s).card : ℝ)))
            / ((sys.shell i).card : ℝ) else 0)
        = if w = ∅ then 0 else 1 / (2 * ((nonzeroStates s).card : ℝ)) := by
      rw [← sys.sum_mem w (1 / (2 * ((nonzeroStates s).card : ℝ)))]
      refine Finset.sum_congr rfl fun i _ => ?_
      by_cases hw : w ∈ sys.shell i
      · simp only [hw, if_true]
        have := (sys.card_pos i).ne'
        field_simp
      · simp [hw]
    rw [actLaw, hrow, hshell]
    by_cases hwq : w = q
    · simp [hwq]
    · simp [hwq]

end Spin.Imt
