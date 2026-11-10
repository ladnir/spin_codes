import SpinCodes.Turbo.PhysicalProfileNumericDefs22
import SpinCodes.Numeric.Fixed
import Mathlib.Algebra.Polynomial.Coeff

/-!
Soundness of the stored linear-factor recurrence. The kernel executes only
the integer operations in the companion definitions file. Real polynomials
appear here solely in the once-proved semantic bridge.
-/
noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix Polynomial

private theorem profileFixZeroMem22 : Mem Fix.zero 0 := by
  simp [Mem, Fix.zero]

def profileRealLinearAux22 (a b prev : ℝ) : List ℝ → List ℝ
  | [] => [b*prev]
  | c :: cs => (a*c+b*prev) :: profileRealLinearAux22 a b c cs

def profileRealLinear22 (a b : ℝ) (cs : List ℝ) : List ℝ :=
  profileRealLinearAux22 a b 0 cs

def profileRealCoefficients22 : List (ℝ × ℝ) → List ℝ
  | [] => [1]
  | (a,b) :: factors => profileRealLinear22 a b (profileRealCoefficients22 factors)

def profileRealPolynomial22 : List ℝ → Polynomial ℝ
  | [] => 0
  | c :: cs => C c + X * profileRealPolynomial22 cs

theorem profileFixLinearAux22_mem {a b : Fix} {ar br : ℝ}
    (ha : Mem a ar) (hb : Mem b br) {cs : List Fix} {ds : List ℝ}
    (hc : List.Forall₂ Mem cs ds) {prev : Fix} {pr : ℝ} (hp : Mem prev pr) :
    List.Forall₂ Mem (profileFixLinearAux22 a b prev cs) (profileRealLinearAux22 ar br pr ds) := by
  induction hc generalizing prev pr with
  | nil => exact .cons (mul_mem hb hp) .nil
  | cons hcd htail ih =>
    exact .cons (add_mem (mul_mem ha hcd) (mul_mem hb hp)) (ih hcd)

theorem profileFixLinear22_mem {a b : Fix} {ar br : ℝ}
    (ha : Mem a ar) (hb : Mem b br) {cs : List Fix} {ds : List ℝ}
    (hc : List.Forall₂ Mem cs ds) :
    List.Forall₂ Mem (profileFixLinear22 a b cs) (profileRealLinear22 ar br ds) :=
  profileFixLinearAux22_mem ha hb hc profileFixZeroMem22

theorem profileFixCoefficients22_mem {fs : List (Fix × Fix)} {rs : List (ℝ × ℝ)}
    (hf : List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2) fs rs) :
    List.Forall₂ Mem (profileFixCoefficients22 fs) (profileRealCoefficients22 rs) := by
  induction hf with
  | nil => exact .cons (by simpa only [Int.cast_one] using ofInt_mem 1) .nil
  | cons hfr htail ih => exact profileFixLinear22_mem hfr.1 hfr.2 ih

theorem profileRealLinearAux22_polynomial (a b prev : ℝ) (cs : List ℝ) :
    profileRealPolynomial22 (profileRealLinearAux22 a b prev cs) =
      (C a + X*C b) * profileRealPolynomial22 cs + C (b*prev) := by
  induction cs generalizing prev with
  | nil => simp [profileRealLinearAux22, profileRealPolynomial22]
  | cons c cs ih =>
    simp only [profileRealLinearAux22, profileRealPolynomial22, ih, map_add, map_mul]
    ring

theorem profileRealLinear22_polynomial (a b : ℝ) (cs : List ℝ) :
    profileRealPolynomial22 (profileRealLinear22 a b cs) =
      (C a + X*C b) * profileRealPolynomial22 cs := by
  simpa only [profileRealLinear22, mul_zero, map_zero, add_zero] using
    profileRealLinearAux22_polynomial a b 0 cs

theorem profileRealCoefficients22_polynomial (factors : List (ℝ × ℝ)) :
    profileRealPolynomial22 (profileRealCoefficients22 factors) =
      (factors.map (fun p : ℝ × ℝ => C p.1 + X*C p.2)).prod := by
  induction factors with
  | nil => simp [profileRealCoefficients22, profileRealPolynomial22]
  | cons factor factors ih =>
    rcases factor with ⟨a,b⟩
    simp only [profileRealCoefficients22, profileRealLinear22_polynomial, ih, List.map_cons,
      List.prod_cons]

theorem profileRealPolynomial22_coeff (cs : List ℝ) (j : Nat) :
    (profileRealPolynomial22 cs).coeff j = cs.getD j 0 := by
  induction cs generalizing j with
  | nil => simp [profileRealPolynomial22]
  | cons c cs ih =>
    cases j with
    | zero => simp [profileRealPolynomial22]
    | succ j => simp [profileRealPolynomial22, ih]

theorem profileForall₂_mem_getD {cs : List Fix} {ds : List ℝ}
    (h : List.Forall₂ Mem cs ds) (j : Nat) : Mem (cs.getD j Fix.zero) (ds.getD j 0) := by
  induction h generalizing j with
  | nil => exact profileFixZeroMem22
  | cons hcd htail ih =>
    cases j with
    | zero => exact hcd
    | succ j => exact ih j

theorem profileFixCoefficients22_sound {fs : List (Fix × Fix)} {rs : List (ℝ × ℝ)}
    (hf : List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2) fs rs) (j : Nat) :
    Mem ((profileFixCoefficients22 fs).getD j Fix.zero)
      (((rs.map (fun p : ℝ × ℝ => C p.1 + X*C p.2)).prod).coeff j) := by
  rw [← profileRealCoefficients22_polynomial, profileRealPolynomial22_coeff]
  exact profileForall₂_mem_getD (profileFixCoefficients22_mem hf) j

theorem profileFixUpperCheck22_sound {box : Fix} {x : ℝ} (hx : Mem box x)
    {num den : Int} (hd : 0 < den) (h : profileFixUpperCheck22 box num den = true) :
    x ≤ (num : ℝ) / (den : ℝ) := by
  have hb : box.hi * den ≤ num * Spin.Numeric.scale := of_decide_eq_true h
  have hr : (box.hi : ℝ) * (den : ℝ) ≤ (num : ℝ) * (Spin.Numeric.scale : ℝ) := by
    exact_mod_cast hb
  have hdR : (0 : ℝ) < (den : ℝ) := by exact_mod_cast hd
  exact hx.2.trans ((div_le_div_iff₀ scaleR_pos hdR).mpr hr)

theorem profileFixLowerCheck22_sound {box : Fix} {x : ℝ} (hx : Mem box x)
    {num den : Int} (hd : 0 < den) (h : profileFixLowerCheck22 box num den = true) :
    (num : ℝ) / (den : ℝ) ≤ x := by
  have hb : num * Spin.Numeric.scale ≤ box.lo * den := of_decide_eq_true h
  have hr : (num : ℝ) * (Spin.Numeric.scale : ℝ) ≤ (box.lo : ℝ) * (den : ℝ) := by
    exact_mod_cast hb
  have hdR : (0 : ℝ) < (den : ℝ) := by exact_mod_cast hd
  exact ((div_le_div_iff₀ hdR scaleR_pos).mpr hr).trans hx.1

#print axioms profileFixCoefficients22_sound
end Spin.Turbo
