import SpinCodes.Structured.DenseFourierExactColumnDefs
import SpinCodes.Structured.DenseFourierExactTarget

noncomputable section
namespace Spin.Structured.DenseFourierExact
open Spin.Imt ConcreteFourier Finset

def ICoords.real (v : ICoords) : Coords 5 := ⟨v.Z,v.D,fun i => v.S i⟩

lemma zeroTerm_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd)
    (v : ICoords) {j : ℕ} (hj : j≤128) :
    Occupation.probability 128 j ((qn:ℝ)/qd)*
      (((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.Data.weights.getD j SparsePolynomial.Data.weight0) ((zn:ℝ)/zd))).applyCol v.real).Z)=
      (zeroTerm qn qd zn zd v j:ℝ)/((qd:ℝ)*zd)^128 := by
  have hq : (qd:ℝ)≠0 := by exact_mod_cast hqd.ne'
  have hz : (zd:ℝ)≠0 := by exact_mod_cast hzd.ne'
  have hc : (Nat.choose 128 j:ℝ)≠0 := by exact_mod_cast (Nat.choose_pos hj).ne'
  rw [Occupation.fixed_col_Z]
  simp only [Occupation.probability,Occupation.Sparse.row,Occupation.Sparse.live,ICoords.real,zeroTerm,
    Int.cast_add,Int.cast_mul,Int.cast_sub,Int.cast_pow,Int.cast_natCast,polyChoose_eq]
  have heq : 1-(qn:ℝ)/qd=((qd:ℝ)-qn)/qd := by field_simp
  have heD : ((qd:ℝ)*zd)^128=((qd:ℝ)*zd)^j*((qd:ℝ)*zd)^(128-j) := by
    rw [←pow_add,Nat.add_sub_of_le hj]
  rw [heq,heD]
  simp only [div_pow,mul_pow]
  field_simp
  ring

/-- Exact common-denominator polynomial for the zero row of the actual Fourier transfer. -/
theorem zeroColumn_rational (qn qd zn zd : Int) (hqd : 0<qd) (hzd : 0<zd) (v : ICoords) :
    ((matrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).Z=
      (zeroColumn qn qd zn zd v:ℝ)/((qd:ℝ)*zd)^128 := by
  change ((Occupation.Sparse.numericalMatrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).Z=_
  rw [Occupation.Sparse.numericalMatrix,Occupation.matrix_col]
  change (∑ j:Fin 129,Occupation.probability 128 j ((qn:ℝ)/qd)*
    (((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.Data.weights.getD j SparsePolynomial.Data.weight0) ((zn:ℝ)/zd))).applyCol v.real).Z))=_
  have ht (j:Fin 129) := zeroTerm_rational qn qd zn zd hqd hzd v (show (j:ℕ)≤128 by omega)
  simp_rw [ht]
  rw [←sum_div]
  congr 1
  rw [zeroColumn,Int.cast_list_sum,List.map_map,SparsePolynomial.sum_fin_range (fun j => (zeroTerm qn qd zn zd v j:ℝ)) 129]
  rfl

end Spin.Structured.DenseFourierExact


