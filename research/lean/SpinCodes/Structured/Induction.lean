/-
The transfer induction, `eq:imt-moment`:

    E[z^{wt(Y)}]  ≤  e_Z T(β,z)^{N/128} 1.

"Rows of a transfer matrix index the entering coordinate, so weighted measures
propagate as row vectors.  Write `e_Z` for the zero-state unit row vector and
`1` for the all-ones column vector."

The induction itself is short once the one-step statement is in place; what it
needs is only that the step carries a measure dominated by `c` to one
dominated by `c · T`.  That one-step statement is the row-domination claim,
and is a hypothesis here — it is where the verifier's numerical envelopes
enter.  Everything below is the architecture around it.

Also here:

*  the initial condition, `δ_0` is dominated by `e_Z`;
*  mixtures, since `T_occ(β,z) = ∑_j C(128,j) β^j (1-β)^{128-j} T_j(z)` and a
   mixture of dominating matrices dominates the mixture of the measures;
*  the scalar envelope `T_sc`: if every row total is at most `λ`, then
   `e_Z T^R 1 ≤ λ^R`.
-/
import SpinCodes.Structured.Step

namespace Spin.Imt

open Finset

variable {s k : ℕ}

/-! ## Transfer matrices -/

/-- Rows of a transfer matrix, indexed by the entering coordinate. -/
structure Transfer (k : ℕ) where
  rowZ : Coords k
  rowD : Coords k
  rowS : Fin k → Coords k

/-- Weighted measures propagate as row vectors. -/
def Transfer.apply (T : Transfer k) (c : Coords k) : Coords k where
  Z := c.Z * T.rowZ.Z + c.D * T.rowD.Z + ∑ i, c.S i * (T.rowS i).Z
  D := c.Z * T.rowZ.D + c.D * T.rowD.D + ∑ i, c.S i * (T.rowS i).D
  S := fun h => c.Z * T.rowZ.S h + c.D * T.rowD.S h + ∑ i, c.S i * (T.rowS i).S h

/-- The zero-state unit row vector. -/
def Coords.eZ (k : ℕ) : Coords k where
  Z := 1
  D := 0
  S := fun _ => 0

@[simp] lemma Coords.eZ_total : (Coords.eZ k).total = 1 := by
  simp [Coords.total, Coords.eZ]

/-- The initial condition: all mass at the zero state. -/
theorem dominates_eZ (sys : ShellSystem s k) :
    Dominates sys (Coords.eZ k) (fun q => if q = ∅ then (1 : ℝ) else 0) := by
  classical
  refine ⟨fun _ => 0, fun _ => le_refl _, by simp [Coords.eZ], fun q => ?_⟩
  have hS : ∑ i : Fin k,
      (if q ∈ sys.shell i then (Coords.eZ k).S i / ((sys.shell i).card : ℝ) else 0) = 0 := by
    refine Finset.sum_eq_zero fun i _ => ?_
    by_cases hq : q ∈ sys.shell i <;> simp [hq, Coords.eZ]
  rw [hS]
  by_cases hq : q = ∅ <;> simp [hq, Coords.eZ]

/-! ## The induction -/

/-- **`eq:imt-moment`.**  If one step carries a measure dominated by `c` to one
dominated by `c · T`, then after `R` steps the terminal mass is at most
`c₀ T^R 1`. -/
theorem imt_moment {sys : ShellSystem s k} (T : Transfer k)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ c μ, Dominates sys c μ → Dominates sys (T.apply c) (step μ))
    {c₀ : Coords k} {μ₀ : Finset (Fin s) → ℝ} (h₀ : Dominates sys c₀ μ₀) (R : ℕ) :
    ∑ q, (step^[R] μ₀) q ≤ ((T.apply)^[R] c₀).total := by
  have key : ∀ n, Dominates sys ((T.apply)^[n] c₀) (step^[n] μ₀) := by
    intro n
    induction n with
    | zero => simpa using h₀
    | succ m ih =>
        rw [Function.iterate_succ_apply', Function.iterate_succ_apply']
        exact hstep _ _ ih
  exact (key R).total_le

/-- Started from the zero state, which is the paper's `e_Z T^R 1`. -/
theorem imt_moment_eZ {sys : ShellSystem s k} (T : Transfer k)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ c μ, Dominates sys c μ → Dominates sys (T.apply c) (step μ)) (R : ℕ) :
    ∑ q, (step^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q
      ≤ ((T.apply)^[R] (Coords.eZ k)).total :=
  imt_moment T step hstep (dominates_eZ sys) R

/-! ## Mixtures

`T_occ` is a binomial mixture of the fixed-weight transfers, so domination has
to survive nonnegative combinations. -/

def Coords.smul (a : ℝ) (c : Coords k) : Coords k where
  Z := a * c.Z
  D := a * c.D
  S := fun i => a * c.S i

def Coords.sum {ι : Type*} [Fintype ι] (f : ι → Coords k) : Coords k where
  Z := ∑ i, (f i).Z
  D := ∑ i, (f i).D
  S := fun h => ∑ i, (f i).S h

theorem dominates_smul {sys : ShellSystem s k} {c : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : Dominates sys c μ) {a : ℝ} (ha : 0 ≤ a) :
    Dominates sys (Coords.smul a c) (fun q => a * μ q) := by
  classical
  obtain ⟨ν, hν0, hνD, hbd⟩ := h
  refine ⟨fun q => a * ν q, fun q => mul_nonneg ha (hν0 q), ?_, fun q => ?_⟩
  · rw [← Finset.mul_sum]
    exact mul_le_mul_of_nonneg_left hνD ha
  · have hle := mul_le_mul_of_nonneg_left (hbd q) ha
    refine hle.trans_eq ?_
    simp only [Coords.smul, mul_add, Finset.mul_sum]
    congr 1
    · congr 1
      by_cases hq : q = ∅ <;> simp [hq]
    · refine Finset.sum_congr rfl fun i _ => ?_
      by_cases hq : q ∈ sys.shell i <;> simp [hq, mul_div_assoc]

theorem dominates_sum {sys : ShellSystem s k} {ι : Type*} [Fintype ι]
    {cs : ι → Coords k} {μs : ι → Finset (Fin s) → ℝ}
    (h : ∀ i, Dominates sys (cs i) (μs i)) :
    Dominates sys (Coords.sum cs) (fun q => ∑ i, μs i q) := by
  classical
  choose ν hν0 hνD hbd using h
  refine ⟨fun q => ∑ i, ν i q, fun q => Finset.sum_nonneg fun i _ => hν0 i q, ?_, fun q => ?_⟩
  · rw [Finset.sum_comm]
    exact Finset.sum_le_sum fun i _ => hνD i
  · refine (Finset.sum_le_sum fun i _ => hbd i q).trans_eq ?_
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
    congr 1
    · congr 1
      by_cases hq : q = ∅ <;> simp [hq, Coords.sum]
    · rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun i _ => ?_
      by_cases hq : q ∈ sys.shell i
      · simp only [hq, if_true, Coords.sum, Finset.sum_div]
      · simp [hq]

/-! ## The scalar envelope

`T_sc` bounds every row total by a single `λ`, which collapses the matrix
bound to `λ^R`. -/

/-- Nonnegative coordinates. -/
def Coords.Nonneg (c : Coords k) : Prop := 0 ≤ c.Z ∧ 0 ≤ c.D ∧ ∀ i, 0 ≤ c.S i

/-- All rows nonnegative. -/
def Transfer.Nonneg (T : Transfer k) : Prop :=
  T.rowZ.Nonneg ∧ T.rowD.Nonneg ∧ ∀ i, (T.rowS i).Nonneg

theorem Transfer.apply_nonneg {T : Transfer k} (hT : T.Nonneg) {c : Coords k}
    (hc : c.Nonneg) : (T.apply c).Nonneg := by
  obtain ⟨hZ, hD, hS⟩ := hc
  obtain ⟨hrZ, hrD, hrS⟩ := hT
  refine ⟨?_, ?_, fun h => ?_⟩ <;>
    exact add_nonneg (add_nonneg (mul_nonneg hZ (by first
        | exact hrZ.1 | exact hrZ.2.1 | exact hrZ.2.2 _))
      (mul_nonneg hD (by first
        | exact hrD.1 | exact hrD.2.1 | exact hrD.2.2 _)))
      (Finset.sum_nonneg fun i _ => mul_nonneg (hS i) (by first
        | exact (hrS i).1 | exact (hrS i).2.1 | exact (hrS i).2.2 _))

/-- `total` of a row-vector product, expanded. -/
lemma total_apply (T : Transfer k) (c : Coords k) :
    (T.apply c).total = c.Z * T.rowZ.total + c.D * T.rowD.total
      + ∑ i, c.S i * (T.rowS i).total := by
  classical
  have hswap : ∑ h : Fin k, ∑ i : Fin k, c.S i * (T.rowS i).S h
      = ∑ i : Fin k, c.S i * (∑ h : Fin k, (T.rowS i).S h) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun i _ => (Finset.mul_sum _ _ _).symm
  have e1 : ∑ h : Fin k, (T.apply c).S h
      = c.Z * (∑ h, T.rowZ.S h) + c.D * (∑ h, T.rowD.S h)
        + ∑ i, c.S i * (∑ h, (T.rowS i).S h) := by
    simp only [Transfer.apply]
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum, hswap]
  have e2 : ∑ i, c.S i * (T.rowS i).total
      = ∑ i, c.S i * (T.rowS i).Z + ∑ i, c.S i * (T.rowS i).D
        + ∑ i, c.S i * (∑ h, (T.rowS i).S h) := by
    rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by rw [Coords.total]; ring
  rw [Coords.total, e1, e2]
  simp only [Transfer.apply, Coords.total]
  ring

/-- One step of the scalar envelope. -/
theorem total_apply_le {T : Transfer k} {lam : ℝ}
    (hZ : T.rowZ.total ≤ lam) (hD : T.rowD.total ≤ lam)
    (hS : ∀ i, (T.rowS i).total ≤ lam) {c : Coords k} (hc : c.Nonneg) :
    (T.apply c).total ≤ lam * c.total := by
  classical
  obtain ⟨hcZ, hcD, hcS⟩ := hc
  rw [total_apply]
  have hct : lam * c.total = lam * c.Z + lam * c.D + ∑ i, lam * c.S i := by
    rw [Coords.total, mul_add, mul_add, Finset.mul_sum]
  rw [hct]
  have h1 : c.Z * T.rowZ.total ≤ lam * c.Z := by
    rw [mul_comm lam]
    exact mul_le_mul_of_nonneg_left hZ hcZ
  have h2 : c.D * T.rowD.total ≤ lam * c.D := by
    rw [mul_comm lam]
    exact mul_le_mul_of_nonneg_left hD hcD
  have h3 : ∑ i, c.S i * (T.rowS i).total ≤ ∑ i, lam * c.S i :=
    Finset.sum_le_sum fun i _ => by
      rw [mul_comm lam]
      exact mul_le_mul_of_nonneg_left (hS i) (hcS i)
  linarith

/-- **The scalar envelope**: `e_Z T^R 1 ≤ λ^R`. -/
theorem total_iterate_le {T : Transfer k} {lam : ℝ} (hlam : 0 ≤ lam)
    (hT : T.Nonneg)
    (hZ : T.rowZ.total ≤ lam) (hD : T.rowD.total ≤ lam)
    (hS : ∀ i, (T.rowS i).total ≤ lam) (R : ℕ) :
    ((T.apply)^[R] (Coords.eZ k)).total ≤ lam ^ R := by
  have hn : ∀ n, ((T.apply)^[n] (Coords.eZ k)).Nonneg := by
    intro n
    induction n with
    | zero => exact ⟨by norm_num [Coords.eZ], by norm_num [Coords.eZ], by
        intro i; norm_num [Coords.eZ]⟩
    | succ m ih =>
        rw [Function.iterate_succ_apply']
        exact T.apply_nonneg hT ih
  induction R with
  | zero => simp
  | succ m ih =>
      rw [Function.iterate_succ_apply', pow_succ]
      refine (total_apply_le hZ hD hS (hn m)).trans ?_
      rw [mul_comm (lam ^ m) lam]
      exact mul_le_mul_of_nonneg_left ih hlam

end Spin.Imt
