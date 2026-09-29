import SpinCodes.Structured.DenseOccupationFixedFirstBox
import SpinCodes.Structured.DenseOccupationFixedW001Box
import SpinCodes.Structured.DenseOccupationFixedW002Box
import SpinCodes.Structured.DenseOccupationFixedW003Box
import SpinCodes.Structured.DenseOccupationFixedW004Box
import SpinCodes.Structured.DenseOccupationFixedW005Box
import SpinCodes.Structured.DenseOccupationFixedW006Box
import SpinCodes.Structured.DenseOccupationFixedW007Box
import SpinCodes.Structured.DenseOccupationFixedW008Box

import SpinCodes.Structured.DenseOccupationFixedRate

noncomputable section
namespace Spin.Structured.DenseOccupationFixed.Strip
open Spin.Numeric Spin.Imt Set

/-- A concrete contraction witness for this affine outer segment, with a common finite prefactor. -/
structure PointCertificate (α x : ℝ) where
  p : ℝ
  y : ℝ
  z : ℝ
  radius : ℝ
  w : Coords 5
  p_pos : 0<p
  p_lt_one : p<1
  y_pos : 0<y
  y_lt_one : y<1
  z_pos : 0<z
  z_lt_one : z<1
  radius_pos : 0<radius
  radius_lt_one : radius<1
  zero_eq_one : w.Z=1
  diffuse_floor : (1:ℝ)/568≤w.D
  shell_floor : ∀ i, (1:ℝ)/568≤w.S i
  collatz : ((Occupation.Sparse.numericalMatrix (p*y) z).applyCol w).le (Coords.smul radius w)
  exponent_le : boxExponent (833/500) (-43311/250000) p y radius z α x≤-(4/10000000)

def witnessFirst {α x : ℝ}
    (hα : α∈Icc FirstBox.a0.real FirstBox.a1.real)
    (hx : x∈Icc FirstBox.x0.real FirstBox.x1.real) : PointCertificate α x where
  p := FirstBox.p.real
  y := FirstBox.y.real
  z := FirstBox.z.real
  radius := FirstBox.radius.real
  w := First.w
  p_pos := by norm_num [FirstBox.p,QInput.real]
  p_lt_one := by norm_num [FirstBox.p,QInput.real]
  y_pos := by norm_num [FirstBox.y,QInput.real]
  y_lt_one := by norm_num [FirstBox.y,QInput.real]
  z_pos := by rw [FirstBox.parameters_match.2.2]; exact First.parameters.2.2.1
  z_lt_one := by rw [FirstBox.parameters_match.2.2]; exact First.parameters.2.2.2.1
  radius_pos := by rw [FirstBox.parameters_match.2.1]; exact First.parameters.2.2.2.2.1
  radius_lt_one := by rw [FirstBox.parameters_match.2.1]; exact First.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [First.w,First.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/372 by norm_num).trans First.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/372 by norm_num).trans (First.witness_floor.2.2 i)
  collatz := FirstBox.collatz
  exponent_le := by
    simpa [FirstBox.m,FirstBox.c,QInput.real] using FirstBox.exponent_bound hα hx

def witnessW001 {α x : ℝ}
    (hα : α∈Icc W001Box.a0.real W001Box.a1.real)
    (hx : x∈Icc W001Box.x0.real W001Box.x1.real) : PointCertificate α x where
  p := W001Box.p.real
  y := W001Box.y.real
  z := W001Box.z.real
  radius := W001Box.radius.real
  w := W001.w
  p_pos := by norm_num [W001Box.p,QInput.real]
  p_lt_one := by norm_num [W001Box.p,QInput.real]
  y_pos := by norm_num [W001Box.y,QInput.real]
  y_lt_one := by norm_num [W001Box.y,QInput.real]
  z_pos := by rw [W001Box.parameters_match.2.2]; exact W001.parameters.2.2.1
  z_lt_one := by rw [W001Box.parameters_match.2.2]; exact W001.parameters.2.2.2.1
  radius_pos := by rw [W001Box.parameters_match.2.1]; exact W001.parameters.2.2.2.2.1
  radius_lt_one := by rw [W001Box.parameters_match.2.1]; exact W001.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W001.w,W001.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/176 by norm_num).trans W001.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/176 by norm_num).trans (W001.witness_floor.2.2 i)
  collatz := W001Box.collatz
  exponent_le := by
    simpa [W001Box.m,W001Box.c,QInput.real] using W001Box.exponent_bound hα hx

def witnessW002 {α x : ℝ}
    (hα : α∈Icc W002Box.a0.real W002Box.a1.real)
    (hx : x∈Icc W002Box.x0.real W002Box.x1.real) : PointCertificate α x where
  p := W002Box.p.real
  y := W002Box.y.real
  z := W002Box.z.real
  radius := W002Box.radius.real
  w := W002.w
  p_pos := by norm_num [W002Box.p,QInput.real]
  p_lt_one := by norm_num [W002Box.p,QInput.real]
  y_pos := by norm_num [W002Box.y,QInput.real]
  y_lt_one := by norm_num [W002Box.y,QInput.real]
  z_pos := by rw [W002Box.parameters_match.2.2]; exact W002.parameters.2.2.1
  z_lt_one := by rw [W002Box.parameters_match.2.2]; exact W002.parameters.2.2.2.1
  radius_pos := by rw [W002Box.parameters_match.2.1]; exact W002.parameters.2.2.2.2.1
  radius_lt_one := by rw [W002Box.parameters_match.2.1]; exact W002.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W002.w,W002.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/100 by norm_num).trans W002.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/100 by norm_num).trans (W002.witness_floor.2.2 i)
  collatz := W002Box.collatz
  exponent_le := by
    simpa [W002Box.m,W002Box.c,QInput.real] using W002Box.exponent_bound hα hx

def witnessW003 {α x : ℝ}
    (hα : α∈Icc W003Box.a0.real W003Box.a1.real)
    (hx : x∈Icc W003Box.x0.real W003Box.x1.real) : PointCertificate α x where
  p := W003Box.p.real
  y := W003Box.y.real
  z := W003Box.z.real
  radius := W003Box.radius.real
  w := W003.w
  p_pos := by norm_num [W003Box.p,QInput.real]
  p_lt_one := by norm_num [W003Box.p,QInput.real]
  y_pos := by norm_num [W003Box.y,QInput.real]
  y_lt_one := by norm_num [W003Box.y,QInput.real]
  z_pos := by rw [W003Box.parameters_match.2.2]; exact W003.parameters.2.2.1
  z_lt_one := by rw [W003Box.parameters_match.2.2]; exact W003.parameters.2.2.2.1
  radius_pos := by rw [W003Box.parameters_match.2.1]; exact W003.parameters.2.2.2.2.1
  radius_lt_one := by rw [W003Box.parameters_match.2.1]; exact W003.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W003.w,W003.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/59 by norm_num).trans W003.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/59 by norm_num).trans (W003.witness_floor.2.2 i)
  collatz := W003Box.collatz
  exponent_le := by
    simpa [W003Box.m,W003Box.c,QInput.real] using W003Box.exponent_bound hα hx

def witnessW004 {α x : ℝ}
    (hα : α∈Icc W004Box.a0.real W004Box.a1.real)
    (hx : x∈Icc W004Box.x0.real W004Box.x1.real) : PointCertificate α x where
  p := W004Box.p.real
  y := W004Box.y.real
  z := W004Box.z.real
  radius := W004Box.radius.real
  w := W004.w
  p_pos := by norm_num [W004Box.p,QInput.real]
  p_lt_one := by norm_num [W004Box.p,QInput.real]
  y_pos := by norm_num [W004Box.y,QInput.real]
  y_lt_one := by norm_num [W004Box.y,QInput.real]
  z_pos := by rw [W004Box.parameters_match.2.2]; exact W004.parameters.2.2.1
  z_lt_one := by rw [W004Box.parameters_match.2.2]; exact W004.parameters.2.2.2.1
  radius_pos := by rw [W004Box.parameters_match.2.1]; exact W004.parameters.2.2.2.2.1
  radius_lt_one := by rw [W004Box.parameters_match.2.1]; exact W004.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W004.w,W004.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/37 by norm_num).trans W004.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/37 by norm_num).trans (W004.witness_floor.2.2 i)
  collatz := W004Box.collatz
  exponent_le := by
    simpa [W004Box.m,W004Box.c,QInput.real] using W004Box.exponent_bound hα hx

def witnessW005 {α x : ℝ}
    (hα : α∈Icc W005Box.a0.real W005Box.a1.real)
    (hx : x∈Icc W005Box.x0.real W005Box.x1.real) : PointCertificate α x where
  p := W005Box.p.real
  y := W005Box.y.real
  z := W005Box.z.real
  radius := W005Box.radius.real
  w := W005.w
  p_pos := by norm_num [W005Box.p,QInput.real]
  p_lt_one := by norm_num [W005Box.p,QInput.real]
  y_pos := by norm_num [W005Box.y,QInput.real]
  y_lt_one := by norm_num [W005Box.y,QInput.real]
  z_pos := by rw [W005Box.parameters_match.2.2]; exact W005.parameters.2.2.1
  z_lt_one := by rw [W005Box.parameters_match.2.2]; exact W005.parameters.2.2.2.1
  radius_pos := by rw [W005Box.parameters_match.2.1]; exact W005.parameters.2.2.2.2.1
  radius_lt_one := by rw [W005Box.parameters_match.2.1]; exact W005.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W005.w,W005.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/26 by norm_num).trans W005.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/26 by norm_num).trans (W005.witness_floor.2.2 i)
  collatz := W005Box.collatz
  exponent_le := by
    simpa [W005Box.m,W005Box.c,QInput.real] using W005Box.exponent_bound hα hx

def witnessW006 {α x : ℝ}
    (hα : α∈Icc W006Box.a0.real W006Box.a1.real)
    (hx : x∈Icc W006Box.x0.real W006Box.x1.real) : PointCertificate α x where
  p := W006Box.p.real
  y := W006Box.y.real
  z := W006Box.z.real
  radius := W006Box.radius.real
  w := W006.w
  p_pos := by norm_num [W006Box.p,QInput.real]
  p_lt_one := by norm_num [W006Box.p,QInput.real]
  y_pos := by norm_num [W006Box.y,QInput.real]
  y_lt_one := by norm_num [W006Box.y,QInput.real]
  z_pos := by rw [W006Box.parameters_match.2.2]; exact W006.parameters.2.2.1
  z_lt_one := by rw [W006Box.parameters_match.2.2]; exact W006.parameters.2.2.2.1
  radius_pos := by rw [W006Box.parameters_match.2.1]; exact W006.parameters.2.2.2.2.1
  radius_lt_one := by rw [W006Box.parameters_match.2.1]; exact W006.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W006.w,W006.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/97 by norm_num).trans W006.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/97 by norm_num).trans (W006.witness_floor.2.2 i)
  collatz := W006Box.collatz
  exponent_le := by
    simpa [W006Box.m,W006Box.c,QInput.real] using W006Box.exponent_bound hα hx

def witnessW007 {α x : ℝ}
    (hα : α∈Icc W007Box.a0.real W007Box.a1.real)
    (hx : x∈Icc W007Box.x0.real W007Box.x1.real) : PointCertificate α x where
  p := W007Box.p.real
  y := W007Box.y.real
  z := W007Box.z.real
  radius := W007Box.radius.real
  w := W007.w
  p_pos := by norm_num [W007Box.p,QInput.real]
  p_lt_one := by norm_num [W007Box.p,QInput.real]
  y_pos := by norm_num [W007Box.y,QInput.real]
  y_lt_one := by norm_num [W007Box.y,QInput.real]
  z_pos := by rw [W007Box.parameters_match.2.2]; exact W007.parameters.2.2.1
  z_lt_one := by rw [W007Box.parameters_match.2.2]; exact W007.parameters.2.2.2.1
  radius_pos := by rw [W007Box.parameters_match.2.1]; exact W007.parameters.2.2.2.2.1
  radius_lt_one := by rw [W007Box.parameters_match.2.1]; exact W007.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W007.w,W007.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/124 by norm_num).trans W007.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/124 by norm_num).trans (W007.witness_floor.2.2 i)
  collatz := W007Box.collatz
  exponent_le := by
    simpa [W007Box.m,W007Box.c,QInput.real] using W007Box.exponent_bound hα hx

def witnessW008 {α x : ℝ}
    (hα : α∈Icc W008Box.a0.real W008Box.a1.real)
    (hx : x∈Icc W008Box.x0.real W008Box.x1.real) : PointCertificate α x where
  p := W008Box.p.real
  y := W008Box.y.real
  z := W008Box.z.real
  radius := W008Box.radius.real
  w := W008.w
  p_pos := by norm_num [W008Box.p,QInput.real]
  p_lt_one := by norm_num [W008Box.p,QInput.real]
  y_pos := by norm_num [W008Box.y,QInput.real]
  y_lt_one := by norm_num [W008Box.y,QInput.real]
  z_pos := by rw [W008Box.parameters_match.2.2]; exact W008.parameters.2.2.1
  z_lt_one := by rw [W008Box.parameters_match.2.2]; exact W008.parameters.2.2.2.1
  radius_pos := by rw [W008Box.parameters_match.2.1]; exact W008.parameters.2.2.2.2.1
  radius_lt_one := by rw [W008Box.parameters_match.2.1]; exact W008.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [W008.w,W008.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/568 by norm_num).trans W008.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/568 by norm_num).trans (W008.witness_floor.2.2 i)
  collatz := W008Box.collatz
  exponent_le := by
    simpa [W008Box.m,W008Box.c,QInput.real] using W008Box.exponent_bound hα hx

/-- Exact coverage of the first nine adjacent boxes, with one of nine checked witnesses at every point. -/
theorem covered {α x : ℝ} (hα : α∈Icc (1/10000) (10031/320000))
    (hx : x∈Icc (13/125) (18296026121/120000000000)) : Nonempty (PointCertificate α x) := by
  by_cases hFirst : α≤FirstBox.a1.real
  · apply Nonempty.intro (witnessFirst (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [FirstBox.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hFirst
    · simpa [FirstBox.x0,FirstBox.x1,QInput.real] using hx
  by_cases hW001 : α≤W001Box.a1.real
  · apply Nonempty.intro (witnessW001 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W001Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW001
    · simpa [W001Box.x0,W001Box.x1,QInput.real] using hx
  by_cases hW002 : α≤W002Box.a1.real
  · apply Nonempty.intro (witnessW002 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W002Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW002
    · simpa [W002Box.x0,W002Box.x1,QInput.real] using hx
  by_cases hW003 : α≤W003Box.a1.real
  · apply Nonempty.intro (witnessW003 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W003Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW003
    · simpa [W003Box.x0,W003Box.x1,QInput.real] using hx
  by_cases hW004 : α≤W004Box.a1.real
  · apply Nonempty.intro (witnessW004 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W004Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW004
    · simpa [W004Box.x0,W004Box.x1,QInput.real] using hx
  by_cases hW005 : α≤W005Box.a1.real
  · apply Nonempty.intro (witnessW005 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W005Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW005
    · simpa [W005Box.x0,W005Box.x1,QInput.real] using hx
  by_cases hW006 : α≤W006Box.a1.real
  · apply Nonempty.intro (witnessW006 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W006Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW006
    · simpa [W006Box.x0,W006Box.x1,QInput.real] using hx
  by_cases hW007 : α≤W007Box.a1.real
  · apply Nonempty.intro (witnessW007 (α := α) (x := x) ?_ ?_)
    · constructor
      · norm_num [W007Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
        linarith
      · exact hW007
    · simpa [W007Box.x0,W007Box.x1,QInput.real] using hx
  apply Nonempty.intro (witnessW008 (α := α) (x := x) ?_ ?_)
  · constructor
    · norm_num [W008Box.a0,FirstBox.a1,W001Box.a1,W002Box.a1,W003Box.a1,W004Box.a1,W005Box.a1,W006Box.a1,W007Box.a1,QInput.real] at *
      linarith
    · simpa [W008Box.a1,QInput.real] using hα.2
  · simpa [W008Box.x0,W008Box.x1,QInput.real] using hx
/-- On this covered strip the full finite routing factor is retained and the matrix prefactor is 568. -/
theorem routed_rate {L b R d : Nat} (hL : 0<L) (hb : 0<b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0<ConcreteRoute.totalWeight rows) (hw1 : ConcreteRoute.totalWeight rows<L*b)
    {α x C : ℝ} (hα : α∈Icc (1/10000) (10031/320000))
    (hx : x∈Icc (13/125) (18296026121/120000000000))
    (hden : ConcreteRoute.density rows=α*x) (hround : 128*R=L*b)
    (hd : (d:ℝ)≤(11/100)*((L:ℝ)*b)) (hC0 : 0≤C)
    (hC : C≤Real.exp (((L:ℝ)*b)*α*((833/500)*x-43311/250000))) :
    C*(ConcreteRoutedEncoder.experimentLaw L b R).prob
      (fun ω => ConcreteRoutedEncoder.weight e rows ω≤d) ≤
      (((b:ℝ)+1)^ConcreteRoute.activeRows rows*((L:ℝ)+1)^b)*568*
        Real.exp (-(4/10000000)*((L:ℝ)*b)) := by
  obtain ⟨v⟩ := covered hα hx
  have h := routed_counted_rate hL hb e rows hw0 hw1
    (by linarith [hα.1]) (by linarith [hα.2])
    (by linarith [hx.1]) (by linarith [hx.2]) hden
    v.p_pos v.p_lt_one v.y_pos v.y_lt_one v.radius_pos v.z_pos v.z_lt_one.le
    (by norm_num : (0:ℝ)<1/568)
    (by rw [v.zero_eq_one]; norm_num : (1:ℝ)/568≤v.w.Z)
    v.diffuse_floor v.shell_floor v.collatz hround hd hC0 hC (by simpa only [neg_div] using v.exponent_le)
  simpa only [v.zero_eq_one,one_div,inv_inv] using h

#print axioms covered
#print axioms routed_rate
end Spin.Structured.DenseOccupationFixed.Strip
