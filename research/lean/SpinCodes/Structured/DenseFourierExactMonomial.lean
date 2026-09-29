import SpinCodes.Structured.DenseFourierExactDefs
import SpinCodes.Structured.DenseScalarExactSound
import SpinCodes.Structured.ConcreteFourierRouted

noncomputable section
namespace Spin.Structured.DenseFourierExact
open DenseScalarExact ConcreteFourier ConcreteScalar Finset

lemma bases_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd) :
    1-(qn:ℝ)/qd+((qn:ℝ)/qd)*((zn:ℝ)/zd)=(g0 qn qd zn zd:ℝ)/((qd:ℝ)*zd) ∧
    (qn:ℝ)/qd+(1-(qn:ℝ)/qd)*((zn:ℝ)/zd)=(g1 qn qd zn zd:ℝ)/((qd:ℝ)*zd) := by
  have hq : (qd:ℝ)≠0 := by exact_mod_cast hqd.ne'
  have hz : (zd:ℝ)≠0 := by exact_mod_cast hzd.ne'
  constructor <;> simp only [g0,g1] <;> push_cast <;> field_simp <;> ring

lemma ratios_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<g0 qn qd zn zd) (hB : 0<g1 qn qd zn zd) :
    ratio0 ((qn:ℝ)/qd) ((zn:ℝ)/zd)=(u0 qn qd zn zd:ℝ)/(g0 qn qd zn zd:ℝ) ∧
    ratio1 ((qn:ℝ)/qd) ((zn:ℝ)/zd)=(u1 qn qd zn zd:ℝ)/(g1 qn qd zn zd:ℝ) := by
  have hq : (qd:ℝ)≠0 := by exact_mod_cast hqd.ne'
  have hz : (zd:ℝ)≠0 := by exact_mod_cast hzd.ne'
  have hA' : (0:ℝ)<g0 qn qd zn zd := by exact_mod_cast hA
  have hB' : (0:ℝ)<g1 qn qd zn zd := by exact_mod_cast hB
  obtain ⟨h₀,h₁⟩ := bases_rational qn qd zn zd hqd hzd
  constructor
  · rw [ratio0,h₀]
    have he : 1-2*(((qn:ℝ)/qd)*((zn:ℝ)/zd)/((g0 qn qd zn zd:ℝ)/((qd:ℝ)*zd)))=
        (((qd-qn)*zd-qn*zn:ℤ):ℝ)/(g0 qn qd zn zd:ℝ) := by
      field_simp
      unfold g0
      push_cast
      ring
    rw [he,abs_div,abs_of_pos hA']
    congr 1
    simp [u0,Int.natCast_natAbs,Int.cast_abs]
  · rw [ratio1,h₁]
    have he : 1-2*((qn:ℝ)/qd/((g1 qn qd zn zd:ℝ)/((qd:ℝ)*zd)))=
        (-((qn*zd-(qd-qn)*zn:ℤ):ℝ))/(g1 qn qd zn zd:ℝ) := by
      field_simp
      unfold g1
      push_cast
      ring
    rw [he,abs_div,abs_of_pos hB',abs_neg]
    congr 1
    simp [u1,Int.natCast_natAbs,Int.cast_abs]

lemma monomial_cancel {A B U V : ℝ} (hA : A≠0) (hB : B≠0)
    {d w h : ℕ} (hd : d≤128) (hh : h≤d) (hw : h≤w) (hl : w-h≤128-d) :
    (A^(128-d)*B^d)*((U/A)^(w-h)*(V/B)^h)=
      A^(128-d-(w-h))*B^(d-h)*U^(w-h)*V^h := by
  have heA : 128-d=(128-d-(w-h))+(w-h) := by omega
  have heB : d=(d-h)+h := by omega
  have hpowA : A^(128-d)=A^(128-d-(w-h))*A^(w-h) := by rw [←pow_add,←heA]
  have hpowB : B^d=B^(d-h)*B^h := by rw [←pow_add,←heB]
  rw [hpowA,hpowB]
  simp only [div_pow]
  field_simp

lemma weighted_monomial_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (hA : 0<g0 qn qd zn zd) (hB : 0<g1 qn qd zn zd)
    {d w h : ℕ} (hd : d≤128) (hh : h≤d) (hw : h≤w) (hl : w-h≤128-d) :
    scalarF ((qn:ℝ)/qd) ((zn:ℝ)/zd) d*
      ((ratio0 ((qn:ℝ)/qd) ((zn:ℝ)/zd))^(w-h)*(ratio1 ((qn:ℝ)/qd) ((zn:ℝ)/zd))^h)=
      (monomial qn qd zn zd d w h:ℝ)/((qd:ℝ)*zd)^128 := by
  have hA' : (g0 qn qd zn zd:ℝ)≠0 := by exact_mod_cast hA.ne'
  have hB' : (g1 qn qd zn zd:ℝ)≠0 := by exact_mod_cast hB.ne'
  obtain ⟨h₀,h₁⟩ := ratios_rational qn qd zn zd hqd hzd hA hB
  rw [scalarF_rational qn qd zn zd hqd hzd d hd,h₀,h₁]
  rw [div_mul_eq_mul_div,monomial_cancel hA' hB' hd hh hw hl]
  simp only [monomial,Int.cast_mul,Int.cast_pow]

end Spin.Structured.DenseFourierExact


