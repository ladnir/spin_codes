import SpinCodes.Structured.DenseFourierExactMonomial

noncomputable section
namespace Spin.Structured.DenseFourierExact
open DenseScalarExact ConcreteFourier ConcreteScalar Finset

lemma weighted_overlap_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<g0 qn qd zn zd) (hB : 0<g1 qn qd zn zd)
    {d w : ℕ} (hd : d≤128) (hw : w≤128) :
    scalarF ((qn:ℝ)/qd) ((zn:ℝ)/zd) d*
      overlapCap 128 d w (ratio0 ((qn:ℝ)/qd) ((zn:ℝ)/zd)) (ratio1 ((qn:ℝ)/qd) ((zn:ℝ)/zd))=
      (overlap qn qd zn zd d w:ℝ)/((qd:ℝ)*zd)^128 := by
  have hA' : (0:ℝ)≤g0 qn qd zn zd := by exact_mod_cast hA.le
  have hB' : (0:ℝ)≤g1 qn qd zn zd := by exact_mod_cast hB.le
  have hqd' : (0:ℝ)<qd := by exact_mod_cast hqd
  have hzd' : (0:ℝ)<zd := by exact_mod_cast hzd
  have hF : 0 ≤ scalarF ((qn:ℝ)/qd) ((zn:ℝ)/zd) d := by
    rw [scalarF_rational qn qd zn zd hqd hzd d hd]
    positivity
  unfold overlapCap
  rw [mul_max_of_nonneg _ _ hF,
    weighted_monomial_rational qn qd zn zd hqd hzd hA hB hd (by omega) (by omega) (by omega),
    weighted_monomial_rational qn qd zn zd hqd hzd hA hB hd (by omega) (by omega) (by omega),
    max_div_div_right (by positivity)]
  simp only [overlap,Int.cast_max]

lemma spectrum_cast (qn qd zn zd : Int) (d : ℕ) :
    (spectrum qn qd zn zd d:ℝ)=∑ w:Fin 129,
      (FiberNumerics.Data.spectrum.getD w 0:ℝ)*(overlap qn qd zn zd d w:ℝ) := by
  rw [spectrum,Int.cast_list_sum,List.map_map]
  rw [SparsePolynomial.sum_fin_range (fun w => (FiberNumerics.Data.spectrum.getD w 0:ℝ)*(overlap qn qd zn zd d w:ℝ)) 129]
  congr 2
  funext w
  simp

lemma weighted_spectrum_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<g0 qn qd zn zd) (hB : 0<g1 qn qd zn zd)
    {d : ℕ} (hd : d≤128) :
    scalarF ((qn:ℝ)/qd) ((zn:ℝ)/zd) d*
      spectrumCap d (ratio0 ((qn:ℝ)/qd) ((zn:ℝ)/zd)) (ratio1 ((qn:ℝ)/qd) ((zn:ℝ)/zd))=
      (spectrum qn qd zn zd d:ℝ)/(524288*((qd:ℝ)*zd)^128) := by
  unfold spectrumCap
  rw [←mul_div_assoc,mul_sum,spectrum_cast]
  simp only [←mul_assoc]
  have hterm (w:Fin 129) :
      scalarF ((qn:ℝ)/qd) ((zn:ℝ)/zd) d*(FiberNumerics.Data.spectrum.getD w 0:ℝ)*
        overlapCap 128 d w (ratio0 ((qn:ℝ)/qd) ((zn:ℝ)/zd)) (ratio1 ((qn:ℝ)/qd) ((zn:ℝ)/zd))=
      ((FiberNumerics.Data.spectrum.getD w 0:ℝ)*(overlap qn qd zn zd d w:ℝ))/((qd:ℝ)*zd)^128 := by
    rw [mul_right_comm,weighted_overlap_rational qn qd zn zd hqd hzd hA hB hd (by omega)]
    ring
  simp_rw [hterm]
  rw [←sum_div]
  ring

/-- Each actual Fourier target cap is an exact integer numerator over one
common positive denominator; no tiny absolute interval is introduced. -/
theorem targetCap_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<g0 qn qd zn zd) (hB : 0<g1 qn qd zn zd)
    {d : ℕ} (hd : d≤128) :
    targetCap ((qn:ℝ)/qd) ((zn:ℝ)/zd) d=(target qn qd zn zd d:ℝ)/(denominator qd zd:ℝ) := by
  unfold targetCap
  rw [mul_add,←mul_div_assoc,weighted_spectrum_rational qn qd zn zd hqd hzd hA hB hd,
    scalarF_rational qn qd zn zd hqd hzd d hd]
  simp only [target,denominator,Int.cast_add,Int.cast_mul,Int.cast_pow,Int.cast_ofNat]
  ring

end Spin.Structured.DenseFourierExact

