import SpinCodes.Structured.DenseOccupationFixedMixture

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric Spin.Imt

lemma le_of_interval_separation {a b : Fix} {x y : ℝ}
    (hx : Fix.Mem a x) (hy : Fix.Mem b y) (h : a.hi ≤ b.lo) : x≤y :=
  hx.2.trans ((div_le_div_of_nonneg_right (by exact_mod_cast h) scaleR_pos.le).trans hy.1)

/-- Integer comparisons against an enclosed right-hand side certify the actual real matrix. -/
theorem checked_collatz (qn qd : Int) (hd : 0<qd) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) {v : FCoords} {w : Coords 5} (hv : v.Mem w)
    {radius : Fix} {lam : ℝ} (hr : Fix.Mem radius lam)
    (hZ : (column qn qd zp v).Z.hi ≤ (Fix.mul radius v.Z).lo)
    (hD : (column qn qd zp v).D.hi ≤ (Fix.mul radius v.D).lo)
    (hS : ∀ i, ((column qn qd zp v).S i).hi ≤ (Fix.mul radius (v.S i)).lo) :
    ((Occupation.Sparse.numericalMatrix ((qn:ℝ)/qd) z).applyCol w).le (Coords.smul lam w) := by
  have hm := column_mem qn qd hd zp z hz hv
  exact ⟨le_of_interval_separation hm.1 (Fix.mul_mem hr hv.1) hZ,
    le_of_interval_separation hm.2.1 (Fix.mul_mem hr hv.2.1) hD,
    fun i => le_of_interval_separation (hm.2.2 i) (Fix.mul_mem hr (hv.2.2 i)) (hS i)⟩

end Spin.Structured.DenseOccupationFixed
