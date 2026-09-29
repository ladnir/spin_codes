import SpinCodes.Structured.ConcreteFixedNumericDefs
import SpinCodes.Structured.DenseOccupationFixedVertex
import SpinCodes.Structured.ConcreteOuterDefectNumeric
import SpinCodes.Structured.ConcreteFixedLargeExponent
import SpinCodes.Structured.WeightedNorm

noncomputable section
namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric Spin.Numeric.Fix DenseOccupationFixed

def gThree : ℝ := -(1203/5125)+(4/3)*Real.log (20500/16891)

theorem gThree_eq : gThree = Placement.paperLargeExponent := by
  norm_num [gThree,Placement.paperLargeExponent,Placement.largeExponent]

def rowRate (x s i u v : ℝ) : ℝ :=
  s*x+i-hEnt x-x*Real.log u+Real.log v+gThree+
    (11/100)*(133/125)/(262144/524287)+(4/39)*Real.log 2

theorem rowCheck_sound {x s i u v : QInput}
    (hx : 0 < x.den) (hs : 0 < s.den) (hi : 0 < i.den)
    (h : rowCheck x s i u v=true) :
    rowRate x.real s.real i.real u.real v.real ≤ -(8679/10000000) := by
  unfold rowCheck rowEval at h
  cases hh : fhEnt x.interval 30 with
  | none => simp [hh] at h
  | some eh =>
    cases hu : flogQ u.num u.den 30 with
    | none => simp [hh,hu] at h
    | some eu =>
      cases hv : flogQ v.num v.den 30 with
      | none => simp [hh,hu,hv] at h
      | some ev =>
        cases hg : flogQ 20500 16891 30 with
        | none => simp [hh,hu,hv,hg] at h
        | some eg =>
          simp only [hh,hu,hv,hg,Option.bind_eq_bind,Option.bind_some,Option.pure_def,decide_eq_true_eq] at h
          have mx : Mem x.interval x.real := ofFrac_mem hx
          have ms : Mem s.interval s.real := ofFrac_mem hs
          have mi : Mem i.interval i.real := ofFrac_mem hi
          have mh := fhEnt_mem mx hh
          have mu := flogQ_mem hu
          have mv := flogQ_mem hv
          have mg := flogQ_mem hg
          have mm := add_mem (add_mem (add_mem (add_mem
            (add_mem (sub_mem (sub_mem (add_mem (mul_mem ms mx) mi) mh) (mul_mem mx mu)) mv)
            (add_mem (ofFrac_mem (p:= -1203) (q:=5125) (by norm_num))
              (mul_mem (ofFrac_mem (p:=4) (q:=3) (by norm_num)) mg)))
            (ofFrac_mem (p:=767031881) (q:=3276800000) (by norm_num)))
            (mul_mem (ofFrac_mem (p:=4) (q:=39) (by norm_num)) log2_mem)) Fix.zero_mem
          have bound := div_le_div_of_nonneg_right (show ((_:Int):ℝ)≤((ofFrac (-8679) 10000000).lo:ℝ) by exact_mod_cast h) scaleR_pos.le
          have target := (ofFrac_mem (p:= -8679) (q:=10000000) (by norm_num)).1
          have result := mm.2.trans (bound.trans target)
          norm_num only [Int.cast_ofNat,Int.cast_neg,add_zero] at result
          convert result using 1 <;> unfold rowRate gThree QInput.real <;> ring

#print axioms gThree_eq
#print axioms rowCheck_sound
end Spin.Structured.ConcreteFixedNumeric
