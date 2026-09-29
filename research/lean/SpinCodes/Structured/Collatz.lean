/-
The Collatz bound.

`eq:imt-dense-exponent` allows `λ` to be either a scalar row bound (done in
`Induction.lean`) or "a positive Collatz bound `T(py,z) w ≤ λ w`".  The second
is the useful one: it is a *right* eigen-inequality against a positive column
vector, and unlike the scalar bound it is not forced to dominate every row
total separately.

The argument is the adjointness `(c T) · w = c · (T w)`, which turns `R`
applications of the row action into `R` factors of `λ`:

    (e_Z T^R) · w  ≤  λ^R · (e_Z · w)  =  λ^R · w_Z,

and then `c.total · w_min ≤ c · w` converts the pairing back to the all-ones
column.  The conclusion is stated as

    (e_Z T^R 1) · w_min  ≤  λ^R · w_Z,

which avoids a division and is sharper than the `λ^R · (max w / min w)` of the
paper, since `e_Z` sees only `w_Z`.

"No numerical eigenvalue is trusted": `λ` and `w` are inputs, and the only
thing checked is the inequality `T w ≤ λ w`.
-/
import SpinCodes.Structured.Convexity

namespace Spin.Imt

open Finset

variable {k : ℕ}

/-! ## The pairing and the column action -/

/-- The pairing of a row vector with a column vector. -/
def Coords.dot (c d : Coords k) : ℝ := c.Z * d.Z + c.D * d.D + ∑ i, c.S i * d.S i

/-- Right multiplication by a column vector: `(T w)_i = ⟨row_i, w⟩`. -/
def Transfer.applyCol (T : Transfer k) (w : Coords k) : Coords k where
  Z := T.rowZ.dot w
  D := T.rowD.dot w
  S := fun h => (T.rowS h).dot w

/-- **Adjointness**: `(c T) · w = c · (T w)`. -/
theorem dot_apply (T : Transfer k) (c w : Coords k) :
    (T.apply c).dot w = c.dot (T.applyCol w) := by
  classical
  have hswap : ∑ h : Fin k, ∑ i : Fin k, c.S i * ((T.rowS i).S h * w.S h)
      = ∑ i : Fin k, ∑ h : Fin k, c.S i * ((T.rowS i).S h * w.S h) :=
    Finset.sum_comm
  simp only [Coords.dot, Transfer.apply, Transfer.applyCol, add_mul, mul_add,
    Finset.sum_mul, Finset.mul_sum, Finset.sum_add_distrib, mul_assoc]
  rw [hswap]
  ring

/-! ## Componentwise order -/

/-- Componentwise `≤` on coordinates. -/
def Coords.le (c d : Coords k) : Prop := c.Z ≤ d.Z ∧ c.D ≤ d.D ∧ ∀ i, c.S i ≤ d.S i

lemma Coords.dot_smul (c w : Coords k) (lam : ℝ) :
    c.dot (Coords.smul lam w) = lam * c.dot w := by
  simp only [Coords.dot, Coords.smul]
  have hs : ∑ i, c.S i * (lam * w.S i) = lam * ∑ i, c.S i * w.S i := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [hs]
  ring

lemma dot_mono {c d d' : Coords k} (hc : c.Nonneg) (h : d.le d') :
    c.dot d ≤ c.dot d' := by
  obtain ⟨hZ, hD, hS⟩ := hc
  obtain ⟨hdZ, hdD, hdS⟩ := h
  have h3 : ∑ i, c.S i * d.S i ≤ ∑ i, c.S i * d'.S i :=
    Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hdS i) (hS i)
  have h1 : c.Z * d.Z ≤ c.Z * d'.Z := mul_le_mul_of_nonneg_left hdZ hZ
  have h2 : c.D * d.D ≤ c.D * d'.D := mul_le_mul_of_nonneg_left hdD hD
  simp only [Coords.dot]
  linarith

/-! ## The bound -/

/-- One step against the eigen-inequality. -/
theorem collatz_step {T : Transfer k} {w : Coords k} {lam : ℝ}
    (hw : (T.applyCol w).le (Coords.smul lam w)) {c : Coords k} (hc : c.Nonneg) :
    (T.apply c).dot w ≤ lam * c.dot w := by
  rw [dot_apply]
  have := dot_mono hc hw
  rwa [Coords.dot_smul] at this

/-- `R` steps give `λ^R`. -/
theorem collatz_iterate {T : Transfer k} (hT : T.Nonneg) {w : Coords k} {lam : ℝ}
    (hlam : 0 ≤ lam) (hw : (T.applyCol w).le (Coords.smul lam w))
    {c : Coords k} (hc : c.Nonneg) (R : ℕ) :
    ((T.apply)^[R] c).dot w ≤ lam ^ R * c.dot w := by
  have hn : ∀ n, ((T.apply)^[n] c).Nonneg := by
    intro n
    induction n with
    | zero => simpa using hc
    | succ m ih =>
        rw [Function.iterate_succ_apply']
        exact T.apply_nonneg hT ih
  induction R with
  | zero => simp
  | succ m ih =>
      rw [Function.iterate_succ_apply', pow_succ]
      refine (collatz_step hw (hn m)).trans ?_
      rw [mul_comm (lam ^ m) lam, mul_assoc]
      exact mul_le_mul_of_nonneg_left ih hlam

/-- Pairing against a column bounded below converts to the all-ones column. -/
lemma total_mul_le_dot {c w : Coords k} (hc : c.Nonneg) {wmin : ℝ}
    (hZ : wmin ≤ w.Z) (hD : wmin ≤ w.D) (hS : ∀ i, wmin ≤ w.S i) :
    c.total * wmin ≤ c.dot w := by
  obtain ⟨hcZ, hcD, hcS⟩ := hc
  have h3 : ∑ i, c.S i * wmin ≤ ∑ i, c.S i * w.S i :=
    Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hS i) (hcS i)
  have h1 : c.Z * wmin ≤ c.Z * w.Z := mul_le_mul_of_nonneg_left hZ hcZ
  have h2 : c.D * wmin ≤ c.D * w.D := mul_le_mul_of_nonneg_left hD hcD
  have hsum : c.total * wmin = c.Z * wmin + c.D * wmin + ∑ i, c.S i * wmin := by
    simp only [Coords.total, add_mul, Finset.sum_mul]
  simp only [Coords.dot]
  rw [hsum]
  linarith

@[simp] lemma Coords.eZ_dot (w : Coords k) : (Coords.eZ k).dot w = w.Z := by
  simp [Coords.dot, Coords.eZ]

/-- **The Collatz bound.**  A positive column `w` with `T w ≤ λ w` bounds
`e_Z T^R 1`.  Stated multiplicatively, so no division and no eigenvalue. -/
theorem collatz_total {T : Transfer k} (hT : T.Nonneg) {w : Coords k} {lam : ℝ}
    (hw : (T.applyCol w).le (Coords.smul lam w)) {wmin : ℝ} (hwmin : 0 < wmin)
    (hZ : wmin ≤ w.Z) (hD : wmin ≤ w.D) (hS : ∀ i, wmin ≤ w.S i) (R : ℕ) :
    ((T.apply)^[R] (Coords.eZ k)).total * wmin ≤ lam ^ R * w.Z := by
  -- `λ ≥ 0` is forced: `T w` is nonnegative and `w.Z > 0`
  have hwnn : w.Nonneg := ⟨by linarith, by linarith, fun i => by linarith [hS i]⟩
  have hlam : 0 ≤ lam := by
    have hrow : 0 ≤ (T.applyCol w).Z := by
      obtain ⟨hrZ, -, -⟩ := hT
      obtain ⟨hwZ, hwD, hwS⟩ := hwnn
      obtain ⟨hz, hd, hs⟩ := hrZ
      have : 0 ≤ ∑ i, T.rowZ.S i * w.S i :=
        Finset.sum_nonneg fun i _ => mul_nonneg (hs i) (hwS i)
      simp only [Transfer.applyCol, Coords.dot]
      have h1 : 0 ≤ T.rowZ.Z * w.Z := mul_nonneg hz hwZ
      have h2 : 0 ≤ T.rowZ.D * w.D := mul_nonneg hd hwD
      linarith
    have hle : (T.applyCol w).Z ≤ lam * w.Z := hw.1
    have hwZpos : 0 < w.Z := by linarith
    nlinarith
  have heZ : (Coords.eZ k).Nonneg :=
    ⟨by norm_num [Coords.eZ], by norm_num [Coords.eZ], by intro i; norm_num [Coords.eZ]⟩
  have hn : ((T.apply)^[R] (Coords.eZ k)).Nonneg := by
    induction R with
    | zero => simpa using heZ
    | succ m ih =>
        rw [Function.iterate_succ_apply']
        exact T.apply_nonneg hT ih
  refine (total_mul_le_dot hn hZ hD hS).trans ?_
  have h := collatz_iterate hT hlam hw heZ R
  rwa [Coords.eZ_dot] at h

end Spin.Imt
