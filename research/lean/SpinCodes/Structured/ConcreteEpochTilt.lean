import SpinCodes.Structured.ConcreteEpochVarianceBounds
import SpinCodes.Structured.ConcreteEpochKernel

/-! A finite tilted endpoint estimate for the actual empty-epoch kernel. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder

theorem exp_neg_sub_le {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) :
    Real.exp (-x) - Real.exp (-y) ≤ y - x := by
  have hfac : Real.exp (-y) = Real.exp (-x) * Real.exp (x - y) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have hbase := Real.add_one_le_exp (x - y)
  calc
    Real.exp (-x) - Real.exp (-y) = Real.exp (-x) * (1 - Real.exp (x - y)) := by
      rw [hfac]
      ring
    _ ≤ Real.exp (-x) * (y - x) :=
      mul_le_mul_of_nonneg_left (by linarith) (Real.exp_pos _).le
    _ ≤ 1 * (y - x) := mul_le_mul_of_nonneg_right
      (Real.exp_le_one_iff.mpr (by linarith)) (sub_nonneg.mpr hxy)
    _ = _ := one_mul _

theorem exp_neg_lipschitz {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) :
    |Real.exp (-x) - Real.exp (-y)| ≤ |x - y| := by
  rcases le_total x y with hxy | hyx
  · rw [abs_of_nonneg (sub_nonneg.mpr (Real.exp_le_exp.mpr (by linarith))),
      abs_of_nonpos (sub_nonpos.mpr hxy)]
    linarith [exp_neg_sub_le hx hxy]
  · rw [abs_sub_comm, abs_sub_comm x y]
    rw [abs_of_nonneg (sub_nonneg.mpr (Real.exp_le_exp.mpr (by linarith))),
      abs_of_nonpos (sub_nonpos.mpr hyx)]
    linarith [exp_neg_sub_le hy hyx]

theorem expect_abs_sub_le {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f h : α → ℝ) :
    |P.expect f - P.expect h| ≤ P.expect (fun x => |f x - h x|) := by
  unfold Spin.FinPMF.expect
  rw [← Finset.sum_sub_distrib]
  simp only [← mul_sub]
  calc |∑ x, P.p x * (f x - h x)| ≤ ∑ x, |P.p x * (f x - h x)| :=
      Finset.abs_sum_le_sum_abs _ _
    _ = _ := by
      apply Finset.sum_congr rfl
      intro x _
      rw [abs_mul, abs_of_nonneg (P.nonneg x)]

theorem exp_neg_pow (c : ℝ) (n : Nat) : Real.exp (-c) ^ n = Real.exp (-c * n) := by
  rw [← Real.exp_nat_mul]
  congr 1
  ring

theorem emptyKernel_exp (c : ℝ) (g : Nat) (q r : State) :
    emptyKernel (Real.exp (-c)) g q r =
      (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect (fun ts =>
        Real.exp (-c * (outputWeight (fun _ => ∅) ts q : ℝ)) *
          (if emptyState ts q = r then 1 else 0)) := by
  simp only [emptyKernel, exp_neg_pow]

theorem emptyKernel_tilt_deviation {c : ℝ} (hc : 0 ≤ c) (g : Nat) (q r : State) :
    |emptyKernel (Real.exp (-c)) g q r -
      Real.exp (-c * (epochMean * g)) * emptyKernel 1 g q r| ≤
        c * emptyAbsDeviation g q := by
  rw [emptyKernel_exp, emptyKernel_one]
  unfold emptyExpect emptyAbsDeviation
  rw [← expect_const_mul, ← expect_const_mul]
  refine (expect_abs_sub_le _ _ _).trans ?_
  apply Spin.FinPMF.expect_mono
  intro ts
  have hm : 0 ≤ epochMean * (g : ℝ) := by unfold epochMean; positivity
  have hw : (0 : ℝ) ≤ outputWeight (fun _ => ∅) ts q := Nat.cast_nonneg _
  have he := exp_neg_lipschitz (mul_nonneg hc hw) (mul_nonneg hc hm)
  rw [show c * (outputWeight (fun _ => ∅) ts q : ℝ) - c * (epochMean * g) =
    c * ((outputWeight (fun _ => ∅) ts q : ℝ) - epochMean * g) by ring,
    abs_mul, abs_of_nonneg hc] at he
  by_cases hr : emptyState ts q = r
  · simpa only [hr, ite_true, mul_one, neg_mul] using he
  · simp only [hr, ite_false, mul_zero, sub_self, abs_zero]
    exact mul_nonneg hc (abs_nonneg _)

/-- Finite form of the empty-interval approximation, with every law tied
to the concrete encoder. The scale c=θ/L gives the paper's endpoint bound. -/
theorem emptyKernel_tilt_live {c : ℝ} (hc : 0 ≤ c) (g : Nat)
    {q r : State} (hq : q ≠ ∅) (hr : r ≠ ∅) :
    |emptyKernel (Real.exp (-c)) g q r - Real.exp (-c * (epochMean * g)) / 524287| ≤
      c * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g := by
  let e := Real.exp (-c * (epochMean * g))
  have he0 : 0 ≤ e := (Real.exp_pos _).le
  have he1 : e ≤ 1 := by
    apply Real.exp_le_one_iff.mpr
    have hm : 0 ≤ epochMean * (g : ℝ) := by unfold epochMean; positivity
    nlinarith
  have ht : |emptyKernel (Real.exp (-c)) g q r - e * emptyKernel 1 g q r| ≤
      c * (128 * Real.sqrt (3 * g)) :=
    (emptyKernel_tilt_deviation hc g q r).trans
      (mul_le_mul_of_nonneg_left (emptyAbsDeviation_bound g hq) hc)
  have hp : 0 ≤ (1 / 2 : ℝ) ^ g := by positivity
  have hmix : |e * emptyKernel 1 g q r - e / 524287| ≤ (1 / 2 : ℝ) ^ g := by
    rw [show e * emptyKernel 1 g q r - e / 524287 =
      e * (emptyKernel 1 g q r - 1 / 524287) by ring,
      abs_mul, abs_of_nonneg he0]
    calc e * |emptyKernel 1 g q r - 1 / 524287| ≤ e * (1 / 2 : ℝ) ^ g :=
        mul_le_mul_of_nonneg_left (emptyKernel_one_error g hq hr) he0
      _ ≤ 1 * (1 / 2 : ℝ) ^ g := mul_le_mul_of_nonneg_right he1 hp
      _ = _ := one_mul _
  exact (abs_sub_le _ (e * emptyKernel 1 g q r) _).trans (add_le_add ht hmix)

theorem emptyKernel_tilt_scaled {θ L : ℝ} (hθ : 0 ≤ θ) (hL : 0 < L) (g : Nat)
    {q r : State} (hq : q ≠ ∅) (hr : r ≠ ∅) :
    |emptyKernel (Real.exp (-(θ / L))) g q r -
      Real.exp (-(θ / L) * (epochMean * g)) / 524287| ≤
      θ / L * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g :=
  emptyKernel_tilt_live (div_nonneg hθ hL.le) g hq hr

end Spin.Structured.ConcreteEncoder
