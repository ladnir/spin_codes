/-
The two-state coarse-graining of `app:imt-fixed`.

    "The zero state stays zero.  Define the lift `J` by `J(q,0)=1_{q=0}`,
     `J(q,1)=1_{q≠0}`, and the projection `R` by `R(0,q)=1_{q=0}`,
     `R(1,q)=1_{q≠0}/M`.  Thus `RJ = I₂`.  Put `D_γ(u) := diag(1, e^{-γu})`.
     ... equation (eq:imt-empty-limit) approximates `W₀^g` by `J D_γ(tg/L) R`
     with entry error `O(L^{-1/2})`."

Two things are proved here.

*  `R J = I₂`: the projection averages the live block and the lift spreads it
   back, so the round trip is the identity on the coarse space.  This is what
   makes `J D R` a legitimate stand-in rather than a lossy summary.

*  `J D R` has the live block *constant* at `d/M`, the zero-to-zero entry `1`,
   and both cross entries `0`.  So `eq:imt-empty-limit`, which bounds
   `|W₀^g(q,q') − e^{-θμg/L}/M|` uniformly over live `q,q'`, is literally the
   entrywise distance to `J D R`; no further estimate is needed to get from
   one to the other.

Also here: the Lipschitz step `|e^{-a} − e^{-b}| ≤ |a − b|` used to turn the
mean bound of `EmptyEpoch.lean` into `eq:imt-empty-limit`.
-/
import SpinCodes.Structured.EmptyEpoch

namespace Spin

open Finset

variable {s : ℕ}

/-! ## The Lipschitz step -/

/-- `x ↦ e^{-x}` is 1-Lipschitz on `[0,∞)`. -/
theorem abs_expNeg_sub_le {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    |Real.exp (-a) - Real.exp (-b)| ≤ |a - b| := by
  have key : ∀ x y : ℝ, 0 ≤ x → x ≤ y →
      Real.exp (-x) - Real.exp (-y) ≤ y - x := by
    intro x y hx hxy
    have h1 : 1 - (y - x) ≤ Real.exp (-(y - x)) := by
      have := Real.add_one_le_exp (-(y - x))
      linarith
    have hneg : -x + -(y - x) = -y := by ring
    have hfac : Real.exp (-x) - Real.exp (-y)
        = Real.exp (-x) * (1 - Real.exp (-(y - x))) := by
      rw [mul_sub, mul_one, ← Real.exp_add, hneg]
    have hex : Real.exp (-x) ≤ 1 := by
      rw [Real.exp_le_one_iff]
      linarith
    have hnn : 0 ≤ 1 - Real.exp (-(y - x)) := by
      have hle : Real.exp (-(y - x)) ≤ 1 := by
        rw [Real.exp_le_one_iff]
        linarith
      linarith
    rw [hfac]
    calc Real.exp (-x) * (1 - Real.exp (-(y - x)))
        ≤ 1 * (1 - Real.exp (-(y - x))) := mul_le_mul_of_nonneg_right hex hnn
      _ = 1 - Real.exp (-(y - x)) := one_mul _
      _ ≤ y - x := by linarith
  rcases le_total a b with hab | hab
  · have hmono : Real.exp (-b) ≤ Real.exp (-a) := Real.exp_le_exp.mpr (by linarith)
    rw [abs_of_nonneg (by linarith), abs_of_nonpos (by linarith)]
    linarith [key a b ha hab]
  · have hmono : Real.exp (-a) ≤ Real.exp (-b) := Real.exp_le_exp.mpr (by linarith)
    rw [abs_of_nonpos (by linarith), abs_of_nonneg (by linarith)]
    linarith [key b a hb hab]

/-! ## The lift and the projection -/

/-- Splitting a sum over states into the zero state and the live states. -/
lemma sum_split (f : Finset (Fin s) → ℝ) :
    ∑ q : Finset (Fin s), f q = f ∅ + ∑ q ∈ nonzeroStates s, f q := by
  rw [nonzeroStates_eq_erase]
  exact (Finset.add_sum_erase _ f (Finset.mem_univ _)).symm

/-- The lift `J`: coarse coordinate `0` is the zero state, `1` the live ones. -/
def lift (q : Finset (Fin s)) (i : Fin 2) : ℝ :=
  if i = 0 then (if q = ∅ then 1 else 0) else (if q = ∅ then 0 else 1)

/-- The projection `R`: reads off the zero state, and averages the live ones. -/
noncomputable def proj (i : Fin 2) (q : Finset (Fin s)) : ℝ :=
  if i = 0 then (if q = ∅ then 1 else 0)
  else (if q = ∅ then 0 else 1 / ((nonzeroStates s).card : ℝ))

/-- `diag(1, d)` on the coarse space. -/
def coarseD (d : ℝ) (i j : Fin 2) : ℝ :=
  if i = j then (if i = 0 then 1 else d) else 0

/-- **`R J = I₂`.** -/
theorem proj_lift (hs : 0 < (nonzeroStates s).card) (i j : Fin 2) :
    ∑ q : Finset (Fin s), proj i q * lift q j = if i = j then 1 else 0 := by
  have hcard : ((nonzeroStates s).card : ℝ) ≠ 0 := by positivity
  have hlive : ∀ q ∈ nonzeroStates s, q ≠ ∅ := fun q hq =>
    ne_empty_of_mem_nonzeroStates hq
  rw [sum_split]
  by_cases hi : i = 0 <;> by_cases hj : j = 0
  · have hz : ∀ q ∈ nonzeroStates s, proj i q * lift q j = 0 := by
      intro q hq; simp [proj, lift, hi, hj, hlive q hq]
    rw [Finset.sum_congr rfl hz]
    simp [proj, lift, hi, hj]
  · have hz : ∀ q ∈ nonzeroStates s, proj i q * lift q j = 0 := by
      intro q hq; simp [proj, lift, hi, hj, hlive q hq]
    rw [Finset.sum_congr rfl hz]
    simp [proj, lift, hi, hj]
    exact fun h => hj h.symm
  · have hz : ∀ q ∈ nonzeroStates s, proj i q * lift q j = 0 := by
      intro q hq; simp [proj, lift, hi, hj, hlive q hq]
    rw [Finset.sum_congr rfl hz]
    simp [proj, lift, hi, hj]
  · have hij : i = j := by
      apply Fin.ext
      have h1 : i.val < 2 := i.isLt
      have h2 : j.val < 2 := j.isLt
      have h3 : i.val ≠ 0 := fun h => hi (Fin.ext h)
      have h4 : j.val ≠ 0 := fun h => hj (Fin.ext h)
      omega
    have hz : ∀ q ∈ nonzeroStates s, proj i q * lift q j
        = 1 / ((nonzeroStates s).card : ℝ) := by
      intro q hq; simp [proj, lift, hi, hj, hlive q hq]
    rw [Finset.sum_congr rfl hz, Finset.sum_const, nsmul_eq_mul]
    simp [proj, lift, hi, hj, hij]
    field_simp

/-- **The coarse kernel `J D R`**, entry by entry: the live block is constant
at `d/M`, the zero state is fixed, and the cross entries vanish. -/
theorem lift_coarseD_proj (d : ℝ) (q q' : Finset (Fin s)) :
    ∑ i : Fin 2, ∑ j : Fin 2, lift q i * coarseD d i j * proj j q'
      = if q = ∅ then (if q' = ∅ then 1 else 0)
        else (if q' = ∅ then 0 else d / ((nonzeroStates s).card : ℝ)) := by
  simp only [Fin.sum_univ_two, lift, proj, coarseD]
  by_cases hq : q = ∅ <;> by_cases hq' : q' = ∅ <;> simp [hq, hq', div_eq_mul_inv]

/-- The live block of `J D R` is exactly `d/M`, so `eq:imt-empty-limit` is the
entrywise error against it with nothing further to prove. -/
theorem lift_coarseD_proj_live {d : ℝ} {q q' : Finset (Fin s)}
    (hq : q ≠ ∅) (hq' : q' ≠ ∅) :
    ∑ i : Fin 2, ∑ j : Fin 2, lift q i * coarseD d i j * proj j q'
      = d / ((nonzeroStates s).card : ℝ) := by
  rw [lift_coarseD_proj, if_neg hq, if_neg hq']

/-- **The approximation transfers verbatim.**  A uniform bound on the live
entries of `W` against `d/M` is a bound against `J D R`. -/
theorem entry_error_transfer {W : Finset (Fin s) → Finset (Fin s) → ℝ} {d eps : ℝ}
    (hW : ∀ q q', q ≠ ∅ → q' ≠ ∅ →
      |W q q' - d / ((nonzeroStates s).card : ℝ)| ≤ eps)
    {q q' : Finset (Fin s)} (hq : q ≠ ∅) (hq' : q' ≠ ∅) :
    |W q q' - ∑ i : Fin 2, ∑ j : Fin 2, lift q i * coarseD d i j * proj j q'|
      ≤ eps := by
  rw [lift_coarseD_proj_live hq hq']
  exact hW q q' hq hq'

end Spin
