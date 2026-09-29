import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.ConcreteScalarRouted

noncomputable section
namespace Spin.Structured.DenseScalarExact
open ConcreteScalar ConcreteMaps Finset

theorem scalarF_rational (qn qd zn zd : Int) (hqd : 0 < qd) (hzd : 0 < zd)
    (d : ℕ) (hd : d ≤ 128) :
    scalarF ((qn : ℝ)/qd) ((zn : ℝ)/zd) d =
      ((g0 qn qd zn zd : ℝ)^(128-d)*(g1 qn qd zn zd : ℝ)^d)/((qd : ℝ)*zd)^128 := by
  have hq : (qd : ℝ) ≠ 0 := by exact_mod_cast hqd.ne'
  have hz : (zd : ℝ) ≠ 0 := by exact_mod_cast hzd.ne'
  have h0 : 1-(qn : ℝ)/qd+((qn : ℝ)/qd)*((zn : ℝ)/zd) =
      (g0 qn qd zn zd : ℝ)/((qd : ℝ)*zd) := by
    unfold g0
    push_cast
    field_simp
  have h1 : (qn : ℝ)/qd+(1-(qn : ℝ)/qd)*((zn : ℝ)/zd) =
      (g1 qn qd zn zd : ℝ)/((qd : ℝ)*zd) := by
    unfold g1
    push_cast
    field_simp
  rw [scalarF, h0, h1, div_pow, div_pow, div_mul_div_comm, ← pow_add, Nat.sub_add_cancel hd]

theorem check_sound {qn qd zn zd rn rd : Int} (hqd : 0 < qd) (hzd : 0 < zd)
    (hrd : 0 < rd) {d : ℕ} (hd : d ≤ 128) (h : check qn qd zn zd rn rd d = true) :
    scalarF ((qn : ℝ)/qd) ((zn : ℝ)/zd) d ≤ (rn : ℝ)/rd := by
  have hh : g0 qn qd zn zd^(128-d)*g1 qn qd zn zd^d*rd ≤ rn*(qd*zd)^128 := by
    simpa only [check, decide_eq_true_eq] using h
  rw [scalarF_rational qn qd zn zd hqd hzd d hd]
  have hq : (0 : ℝ) < qd := by exact_mod_cast hqd
  have hz : (0 : ℝ) < zd := by exact_mod_cast hzd
  have hr : (0 : ℝ) < rd := by exact_mod_cast hrd
  apply (div_le_div_iff₀ (pow_pos (mul_pos hq hz) _) hr).mpr
  exact_mod_cast hh

/-- Exact integer comparisons certify all six rows of the actual scalar transfer. -/
theorem checked_scalar {qn qd zn zd rn rd : Int} (hqd : 0 < qd) (hzd : 0 < zd)
    (hrd : 0 < rd) (h0 : check qn qd zn zd rn rd 0 = true)
    (hs : ∀ i : Fin 5, check qn qd zn zd rn rd (shellWeight i) = true) :
    scalarBound ((qn : ℝ)/qd) ((zn : ℝ)/zd) ≤ (rn : ℝ)/rd := by
  apply max_le (check_sound hqd hzd hrd (by omega) h0)
  apply (Finset.sup'_le_iff _ _).mpr
  intro i _
  apply check_sound hqd hzd hrd _ (hs i)
  fin_cases i <;> decide

#print axioms check_sound
#print axioms checked_scalar

end Spin.Structured.DenseScalarExact
