import SpinCodes.Structured.DenseFourierExactZero

noncomputable section
namespace Spin.Structured.DenseFourierExact
open Spin.Imt ConcreteFourier ConcreteMaps Finset

lemma weightedTotal_cast (v : ICoords) :
    (weightedTotal v:ℝ)=v.real.Z+∑ i:Fin 5,Occupation.Sparse.count i*v.real.S i := by
  have hc : ∀ i:Fin 5,(SparsePolynomial.counts.getD i 1:ℝ)=Occupation.Sparse.count i := by
    intro i; fin_cases i <;> rfl
  simp only [weightedTotal,Int.cast_add,Int.cast_list_sum,List.map_map]
  change (v.Z:ℝ)+((List.finRange 5).map (fun (i:Fin 5) => (((SparsePolynomial.counts.getD i 1:Int)*v.S i:ℤ):ℝ))).sum=_
  simp only [Int.cast_mul,Int.cast_natCast,hc,←List.ofFn_eq_map,List.sum_ofFn,ICoords.real]

lemma atomBudget_dot (c : ℝ) (v : ICoords) :
    (atomBudget c).dot v.real=c*(weightedTotal v:ℝ) := by
  rw [weightedTotal_cast]
  simp only [atomBudget,Coords.dot,zero_mul,add_zero,mul_add,mul_sum]
  congr 1
  apply sum_congr rfl
  intro i hi
  ring

lemma target_le_targetD (qn qd zn zd : Int) (i : Fin 5) :
    target qn qd zn zd (shellWeight i) ≤ targetD qn qd zn zd := by
  fin_cases i <;> simp [shellWeight,SparsePolynomial.levels,targetD] <;> omega

lemma denominator_pos {qd zd : Int} (hqd : 0<qd) (hzd : 0<zd) : 0<denominator qd zd := by
  unfold denominator
  positivity

lemma diffuseCap_le (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<DenseScalarExact.g0 qn qd zn zd) (hB : 0<DenseScalarExact.g1 qn qd zn zd) :
    diffuseCap ((qn:ℝ)/qd) ((zn:ℝ)/zd) ≤ (targetD qn qd zn zd:ℝ)/(denominator qd zd:ℝ) := by
  apply (sup'_le_iff _ _).mpr
  intro i hi
  rw [targetCap_rational qn qd zn zd hqd hzd hA hB (show shellWeight i≤128 by fin_cases i <;> decide)]
  apply div_le_div_of_nonneg_right _ (by exact_mod_cast (denominator_pos hqd hzd).le)
  exact_mod_cast target_le_targetD qn qd zn zd i

lemma le_of_cross {a d r s v : Int} (hd : 0<d) (hs : 0<s) (h : a*s≤r*d*v) :
    (a:ℝ)/d≤(r:ℝ)/s*(v:ℝ) := by
  rw [div_mul_eq_mul_div]
  apply (div_le_div_iff₀ (by exact_mod_cast hd : (0:ℝ)<d) (by exact_mod_cast hs : (0:ℝ)<s)).mpr
  have hh : (a:ℝ)*s≤(r:ℝ)*d*v := by exact_mod_cast h
  nlinarith

/-- Exact integer Collatz checks certify the concrete Fourier transfer. -/
theorem checked_collatz {qn qd zn zd rn rd : Int} (hqd : 0<qd) (hzd : 0<zd) (hrd : 0<rd)
    (hA : 0<DenseScalarExact.g0 qn qd zn zd) (hB : 0<DenseScalarExact.g1 qn qd zn zd)
    (v : ICoords) (hv : 0≤weightedTotal v)
    (hZ : checkZ qn qd zn zd rn rd v=true)
    (hD : checkD qn qd zn zd rn rd v=true)
    (hS : ∀ i,checkS qn qd zn zd rn rd v i=true) :
    ((matrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).le (Coords.smul ((rn:ℝ)/rd) v.real) := by
  refine ⟨?_,?_,fun i => ?_⟩
  · rw [zeroColumn_rational qn qd zn zd hqd hzd]
    have hh : zeroColumn qn qd zn zd v*rd≤rn*(qd*zd)^128*v.Z := by simpa only [checkZ,decide_eq_true_eq] using hZ
    have hl := le_of_cross (show 0<(qd*zd)^128 by positivity) hrd hh
    simpa only [Int.cast_pow,Int.cast_mul,Coords.smul,ICoords.real] using hl
  · change (atomBudget (diffuseCap ((qn:ℝ)/qd) ((zn:ℝ)/zd))).dot v.real≤_
    rw [atomBudget_dot]
    have hd := mul_le_mul_of_nonneg_right (diffuseCap_le qn qd zn zd hqd hzd hA hB)
      (show (0:ℝ)≤weightedTotal v by exact_mod_cast hv)
    have hh : targetD qn qd zn zd*weightedTotal v*rd≤rn*denominator qd zd*v.D := by
      simpa only [checkD,decide_eq_true_eq] using hD
    have hl := le_of_cross (denominator_pos hqd hzd) hrd hh
    refine hd.trans ?_
    simpa only [Int.cast_mul,div_mul_eq_mul_div,Coords.smul,ICoords.real] using hl
  · change (atomBudget (targetCap ((qn:ℝ)/qd) ((zn:ℝ)/zd) (shellWeight i))).dot v.real≤_
    rw [atomBudget_dot,targetCap_rational qn qd zn zd hqd hzd hA hB (show shellWeight i≤128 by fin_cases i <;> decide)]
    have he : SparsePolynomial.levels.getD i 0=shellWeight i := by fin_cases i <;> rfl
    have hh : target qn qd zn zd (shellWeight i)*weightedTotal v*rd≤rn*denominator qd zd*v.S i := by
      simpa only [checkS,decide_eq_true_eq,he] using hS i
    have hl := le_of_cross (denominator_pos hqd hzd) hrd hh
    simpa only [Int.cast_mul,div_mul_eq_mul_div,Coords.smul,ICoords.real] using hl

end Spin.Structured.DenseFourierExact

