import SpinCodes.Turbo.PhysicalProfile22
import SpinCodes.Turbo.Hi.PhysicalProfileNumericData22
import SpinCodes.Turbo.Hi.PhysicalProfileNumericFactors22
import SpinCodes.Turbo.Hi.PhysicalProfileNumericAggregate22

noncomputable section
namespace Spin.Turbo.Hi
open Spin.UniformMixer Polynomial
open Spin.NumericHi Spin.NumericHi.Fix

theorem profileNumericEntries22_authenticated :
    productionPacketProfileTable.map (fun e => (e.1,List.ofFn e.2)) = profileNumericEntries22 := by
  decide

theorem profileNumericEntries22_wellformed :
    profileNumericEntries22.all (fun e => decide (e.2.length = 5 ∧ e.2.sum = 16 ∧ 0 < e.1)) = true := by
  decide

theorem profileNumericClassCount22_checked : ∀ i : Fin 11,
    profileListClassCount22 profileNumericEntries22 i.val =
      ([1,4194303,4194303,2604,291648,888832,1828134,888832,291648,2604,1] : List Nat).getD i.val 0 := by
  decide

theorem profileNumericClassCount22_authenticated (i : Fin 11) :
    profileListClassCount22 profileNumericEntries22 i.val = physicalProfileClassCount22 i := by
  rw [profileNumericClassCount22_checked, physicalProfileClassCount22_checked]

theorem profileNumericClassCount22_pos (i : Fin 11) :
    0 < profileListClassCount22 profileNumericEntries22 i.val := by
  exact (show ∀ k : Fin 11, 0 < profileListClassCount22 profileNumericEntries22 k.val by decide) i

theorem profileListWeight22_ofFn (p : Fin 5 → Nat) :
    profileListWeight22 (List.ofFn p) = physicalProfileWeight22 p := by
  simp [profileListWeight22, physicalProfileWeight22, List.range_succ,
    List.ofFn_succ, Fin.sum_univ_succ, Nat.add_assoc]

theorem profileListClassWeight22_checked : ∀ i : Fin 11,
    ([0,0,0,16,24,28,32,36,40,48,64] : List Nat).getD i.val 0 = physicalClassWeight22 i := by
  decide

theorem profileListClassAllows22_ofFn (i : Fin 11) (p : Fin 5 → Nat) :
    profileListClassAllows22 i.val (List.ofFn p) = decide (physicalProfileClassAllows22 i p) := by
  unfold profileListClassAllows22 physicalProfileClassAllows22
  rw [profileListWeight22_ofFn, profileListClassWeight22_checked]
  split_ifs <;> rfl

theorem profileListGetD_ofFn (p : Fin 5 → Nat) (r : Fin 5) :
    (List.ofFn p).getD r.val 0 = p r := by
  rw [List.getD_eq_getElem _ _ (by simp only [List.length_ofFn]; exact r.isLt)]
  exact List.getElem_ofFn _

theorem profileRealMomentBank22_getD (z : ℝ) (r : Fin 5) :
    (profileRealMomentBank22 z).getD r.val (0,0) = (z^r.val,(1+z)^4-z^r.val) := by
  fin_cases r <;> rfl

theorem profileRealFourierBank22_getD (z : ℝ) (r : Fin 5) :
    (profileRealFourierBank22 z).getD r.val (0,0) = (1,(1+z)^(4-r.val)*(1-z)^r.val-1) := by
  fin_cases r <;> rfl

theorem profileBankPolynomial22_moment (p : Fin 5 → Nat) (z : ℝ) :
    profileBankPolynomial22 (List.ofFn p) (profileRealMomentBank22 z) =
      physicalProfileProduct22 (fun r : Fin 5 => C (z^r.val) + X*C ((1+z)^4-z^r.val)) p := by
  unfold profileBankPolynomial22 physicalProfileProduct22
  apply Finset.prod_congr rfl
  intro r _
  rw [profileListGetD_ofFn, profileRealMomentBank22_getD]

theorem profileBankPolynomial22_fourier (p : Fin 5 → Nat) (z : ℝ) :
    profileBankPolynomial22 (List.ofFn p) (profileRealFourierBank22 z) =
      physicalProfileProduct22 (fun r : Fin 5 => C (1 : ℝ) +
        X*C ((1+z)^(4-r.val)*(1-z)^r.val-1)) p := by
  unfold profileBankPolynomial22 physicalProfileProduct22
  apply Finset.prod_congr rfl
  intro r _
  rw [profileListGetD_ofFn, profileRealFourierBank22_getD]

theorem profileNumericClassUpperCheck22_actual (i : Fin 11)
    (boxValue : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (boxValue p) (realValue p)) {num den : Int} (hd : 0 < den)
    (h : profileFixClassUpperCheck22 profileNumericEntries22 i.val boxValue num den = true)
    (q : State22) (hq : physicalClassAllows22 i q) :
    realValue (List.ofFn (physicalProfile22 q)) ≤ (num : ℝ)/(den : ℝ) := by
  apply physicalProfile22_class_cap i (fun p => realValue (List.ofFn p)) _ _ q hq
  intro e he ha
  have he' : (e.1,List.ofFn e.2) ∈ profileNumericEntries22 := by
    rw [← profileNumericEntries22_authenticated]
    exact List.mem_map.mpr ⟨e,he,rfl⟩
  have ha' : profileListClassAllows22 i.val (List.ofFn e.2) = true := by
    rw [profileListClassAllows22_ofFn]
    simpa only [decide_eq_true_eq] using ha
  exact profileFixClassUpperCheck22_sound profileNumericEntries22 i.val boxValue realValue hv hd h
    (e.1,List.ofFn e.2) he' ha'

theorem profileNumericClassMean22_actual_mem (i : Fin 11)
    (boxValue : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (boxValue p) (realValue p)) :
    Mem (profileFixClassMean22 profileNumericEntries22 i.val boxValue)
      (physicalProfileClassMean22 i (fun p => realValue (List.ofFn p))) := by
  have hh := profileFixClassMean22_mem profileNumericEntries22 i.val boxValue realValue hv
    (profileNumericClassCount22_pos i)
  have hn : (profileNumericEntries22.map (fun e => if profileListClassAllows22 i.val e.2 then
      (e.1 : ℝ)*realValue e.2 else 0)).sum =
      (productionPacketProfileTable.map (fun e => if physicalProfileClassAllows22 i e.2 then
        (e.1 : ℝ)*realValue (List.ofFn e.2) else 0)).sum := by
    rw [← profileNumericEntries22_authenticated, List.map_map]
    simp only [Function.comp_def, profileListClassAllows22_ofFn, decide_eq_true_eq]
  rw [hn, profileNumericClassCount22_authenticated i] at hh
  exact hh

theorem profileNumericFullSum22_authenticated (realValue : List Nat → ℝ) :
    (profileNumericEntries22.map (fun e => (e.1 : ℝ)*realValue e.2)).sum =
      ∑ q : State22, realValue (List.ofFn (physicalProfile22 q)) := by
  rw [← profileNumericEntries22_authenticated, List.map_map]
  simpa only [Function.comp_def] using
    (physicalProfile22_full_sum (fun p => realValue (List.ofFn p))).symm

theorem profileNumericAllFourier22_actual_mem {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Fin 17) :
    Mem (profileFixAllFourier22 profileNumericEntries22 box j.val)
      ((∑ q : State22, (profileBankPolynomial22 (List.ofFn (physicalProfile22 q))
        (profileRealFourierBank22 z)).coeff j.val) /
          ((4194304 : ℝ) * ((Nat.choose 16 j.val : ℝ)*15^j.val))) := by
  have hh := profileFixAllFourier22_mem profileNumericEntries22 hz j
  rw [profileNumericFullSum22_authenticated
    (fun p => (profileBankPolynomial22 p (profileRealFourierBank22 z)).coeff j.val)] at hh
  exact hh

#print axioms profileNumericEntries22_authenticated
#print axioms profileBankPolynomial22_moment
#print axioms profileNumericClassUpperCheck22_actual
#print axioms profileNumericAllFourier22_actual_mem
end Spin.Turbo.Hi
