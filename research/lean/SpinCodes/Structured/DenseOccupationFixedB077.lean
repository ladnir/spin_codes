import SpinCodes.Structured.DenseOccupationFixedB077Data
import SpinCodes.Structured.DenseOccupationFixedVertex
import SpinCodes.Structured.DenseOccupationFixedW032

noncomputable section
namespace Spin.Structured.DenseOccupationFixed.B077
open Spin.Numeric Set

lemma vertex00 : boxExponent m.real c.real p.real y.real radius.real z.real a0.real x0.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a0 x0 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_00
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex01 : boxExponent m.real c.real p.real y.real radius.real z.real a0.real x1.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a0 x1 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_01
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex10 : boxExponent m.real c.real p.real y.real radius.real z.real a1.real x0.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a1 x0 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_10
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex11 : boxExponent m.real c.real p.real y.real radius.real z.real a1.real x1.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a1 x1 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_11
  norm_num [upper,scale] at h ⊢
  exact h

/-- Every point in the occupation box 77 has the stated negative exponent. -/
theorem exponent_bound {α x : ℝ} (hα : α∈Icc a0.real a1.real) (hx : x∈Icc x0.real x1.real) :
    boxExponent m.real c.real p.real y.real radius.real z.real α x ≤ -(4/10000000:ℝ) := by
  apply boxExponent_le_of_vertices m.real c.real p.real y.real radius.real z.real _
    (by norm_num [p,QInput.real]) (by norm_num [p,QInput.real])
    (by norm_num [y,QInput.real]) (by norm_num [y,QInput.real])
    (by norm_num [a0,QInput.real]) (by norm_num [a1,QInput.real])
    (by norm_num [x0,QInput.real]) (by norm_num [x1,QInput.real]) hα hx
    vertex00 vertex01 vertex10 vertex11

theorem parameters_match : p.real*y.real = W032.qReal ∧ radius.real = W032.radiusReal ∧
    z.real = W032.zReal := by
  norm_num [p,y,radius,z,QInput.real,W032.qReal,W032.qn,W032.qd,
    W032.radiusReal,W032.radius,W032.zReal,W032.z,Fix.sc,scale]

/-- Matrix contraction and the exponent enclosure use exactly the same reference parameters. -/
theorem collatz : ((Spin.Imt.Occupation.Sparse.numericalMatrix (p.real*y.real) z.real).applyCol W032.w).le
    (Spin.Imt.Coords.smul radius.real W032.w) := by
  rw [parameters_match.1, parameters_match.2.1, parameters_match.2.2]
  exact W032.collatz

#print axioms exponent_bound
#print axioms collatz
end Spin.Structured.DenseOccupationFixed.B077
