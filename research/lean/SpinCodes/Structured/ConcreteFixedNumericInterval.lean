import SpinCodes.Structured.ConcreteFixedNumericSound
import SpinCodes.Majorant.RefinedData

noncomputable section
namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric DenseOccupationFixed Set

def rate (x u : ℝ) : ℝ := Spin.Majorant.refined.toFun x-hEnt x-x*Real.log u+
  Real.log (Spin.mv (127/250) (3/1600) (1/524287) u)+gThree+
  (11/100)*(133/125)/(262144/524287)+(4/39)*Real.log 2

theorem rowRate_interval {x lo hi s i u v d : ℝ}
    (hl : 0 ≤ lo) (hh : hi ≤ 1) (hx : x∈Icc lo hi)
    (hlo : rowRate lo s i u v ≤ d) (hhi : rowRate hi s i u v ≤ d) :
    rowRate x s i u v ≤ d := by
  have hb := ConcreteOuter.entropy_defect_interval (s:=s-Real.log u)
    (i:=i+Real.log v+gThree+(11/100)*(133/125)/(262144/524287)+(4/39)*Real.log 2)
    hl hh hx (d:=d) ?_ ?_
  · convert hb using 1 <;> unfold rowRate <;> ring
  · convert hlo using 1 <;> unfold rowRate <;> ring
  · convert hhi using 1 <;> unfold rowRate <;> ring

theorem rate_le_rowRate {s i : ℚ} (h : (s,i)∈Spin.Majorant.refined.supports)
    (x u v : ℝ) (hv : v=Spin.mv (127/250) (3/1600) (1/524287) u) :
    rate x u ≤ rowRate x s i u v := by
  have hh := Spin.Majorant.refined.toFun_le_support h x
  unfold rate rowRate
  rw [hv]
  linarith

#print axioms rowRate_interval
#print axioms rate_le_rowRate
end Spin.Structured.ConcreteFixedNumeric
