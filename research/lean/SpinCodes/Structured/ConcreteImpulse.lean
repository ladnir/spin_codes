import SpinCodes.Structured.ConcreteEpoch

/-! Exact zero/live impulse law for the concrete encoder, including the
uniformly chosen bit position used by the fixed-occupation argument. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps
open scoped symmDiff

theorem liveMean_transvection (f : State → ℝ) :
    liveMean (fun q => transvectionLaw.expect (fun t => f (Spin.act t.val.1 t.val.2 q))) =
      liveMean f := by
  conv_lhs => unfold liveMean
  rw [Finset.sum_congr rfl (fun q hq =>
    transvection_expect_live (Spin.ne_empty_of_mem_nonzeroStates hq) f)]
  simp only [add_div, Finset.sum_add_distrib, ← Finset.sum_div,
    Finset.sum_const, nsmul_eq_mul, nonzeroStates_card_actual]
  unfold liveMean
  push_cast
  ring

theorem liveMean_impulse (x : Input) (f : State → ℝ) :
    liveMean (fun q => transvectionLaw.expect (fun t => f (nextState q x t))) =
      liveMean (fun q => f (q ∆ Cset x)) := by
  exact liveMean_transvection (fun q => f (q ∆ Cset x))

theorem Cset_singleton_nonzero (i : Fin 128) : Cset {i} ≠ ∅ := by
  intro h
  have hc := Cset_singleton_card i
  rw [h, Finset.card_empty] at hc
  omega

def coarseValue (a b : ℝ) (q : State) : ℝ := if q = ∅ then a else b

theorem liveMean_coarse_shift {c : State} (hc : c ≠ ∅) (a b : ℝ) :
    liveMean (fun q => coarseValue a b (q ∆ c)) =
      (1 / 524287 : ℝ) * a + (1 - 1 / 524287 : ℝ) * b := by
  have hmem := Spin.mem_nonzeroStates hc
  have hf : (fun q : State => coarseValue a b (q ∆ c)) =
      fun q => b + (if q = c then a - b else 0) := by
    funext q
    simp only [coarseValue, Finset.symmDiff_eq_empty]
    split_ifs <;> ring
  rw [hf]
  unfold liveMean
  rw [Finset.sum_add_distrib, Finset.sum_const, nsmul_eq_mul,
    Finset.sum_ite_eq', if_pos hmem, nonzeroStates_card_actual]
  push_cast
  ring

theorem impulse_from_zero (i : Fin 128) (a b : ℝ) :
    transvectionLaw.expect (fun t => coarseValue a b (nextState ∅ {i} t)) = b := by
  have he (t : Transvection) : nextState ∅ {i} t = Cset {i} := by
    simp only [nextState, Spin.act_empty]
    exact bot_symmDiff _
  simp only [he, coarseValue, if_neg (Cset_singleton_nonzero i)]
  exact Spin.FinPMF.expect_const _ _

theorem impulse_from_live (i : Fin 128) (a b : ℝ) :
    liveMean (fun q => transvectionLaw.expect (fun t => coarseValue a b (nextState q {i} t))) =
      (1 / 524287 : ℝ) * a + (1 - 1 / 524287 : ℝ) * b := by
  rw [liveMean_impulse]
  exact liveMean_coarse_shift (Cset_singleton_nonzero i) a b

/-- Projection of one actual impulse onto zero/live states; the input live
state is uniform, as in the paper's projection R. -/
def impulseCoarse (i : Fin 128) : Matrix (Fin 2) (Fin 2) ℝ := fun a b =>
  if a = 0 then
    transvectionLaw.expect (fun t => coarseValue (if b = 0 then 1 else 0)
      (if b = 1 then 1 else 0) (nextState ∅ {i} t))
  else
    liveMean (fun q => transvectionLaw.expect (fun t =>
      coarseValue (if b = 0 then 1 else 0) (if b = 1 then 1 else 0)
        (nextState q {i} t)))

def impulseMatrix : Matrix (Fin 2) (Fin 2) ℝ :=
  !![0, 1; 1 / 524287, 1 - 1 / 524287]

theorem impulseCoarse_eq (i : Fin 128) : impulseCoarse i = impulseMatrix := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    norm_num [impulseCoarse, impulseMatrix, impulse_from_zero, impulse_from_live]

def uniformImpulseCoarse : Matrix (Fin 2) (Fin 2) ℝ :=
  fun a b => (∑ i : Fin 128, impulseCoarse i a b) / 128

/-- Exact RW₁(1)J = P₀ for uniform local singleton input. -/
theorem uniformImpulseCoarse_eq : uniformImpulseCoarse = impulseMatrix := by
  ext a b
  simp only [uniformImpulseCoarse, impulseCoarse_eq, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul]
  ring

end Spin.Structured.ConcreteEncoder
