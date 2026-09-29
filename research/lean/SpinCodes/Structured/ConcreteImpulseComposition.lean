import SpinCodes.Structured.ConcreteEpochComposition
import SpinCodes.Structured.ConcreteImpulse

/-! Exact products for isolated singleton impulses separated by empty gaps.
The local coordinates can be fixed or independently uniform. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder

def impulseInputs : {a : Nat} → (Fin a → Nat) → (Fin a → Fin 128) → Nat → List Input
  | 0, _, _, last => List.replicate last ∅
  | _ + 1, gaps, coords, last => List.replicate (gaps 0) ∅ ++
      {coords 0} :: impulseInputs (Fin.tail gaps) (Fin.tail coords) last

def impulseProduct (z : ℝ) :
    {a : Nat} → (Fin a → Nat) → (Fin a → Fin 128) → Nat → Matrix State State ℝ
  | 0, _, _, last => emptyKernel z last
  | _ + 1, gaps, coords, last => Matrix.of (emptyKernel z (gaps 0)) * roundKernel z {coords 0} *
      impulseProduct z (Fin.tail gaps) (Fin.tail coords) last

theorem pathKernel_impulseInputs {a : Nat} (z : ℝ) (gaps : Fin a → Nat)
    (coords : Fin a → Fin 128) (last : Nat) :
    pathKernel z (impulseInputs gaps coords last) = impulseProduct z gaps coords last := by
  induction a with
  | zero => exact pathKernel_empty z last
  | succ a ih =>
    simp only [impulseInputs, pathKernel_append, pathKernel_empty, pathKernel,
      impulseProduct, ih, Matrix.mul_assoc]
    rfl

/-- Actual iid transvections on the deterministic impulse/gap input schedule. -/
theorem endpointKernel_impulseInputs {a : Nat} (z : ℝ) (gaps : Fin a → Nat)
    (coords : Fin a → Fin 128) (last : Nat) :
    endpointKernel z (impulseInputs gaps coords last).get = impulseProduct z gaps coords last := by
  rw [endpointKernel_eq_path, List.ofFn_get, pathKernel_impulseInputs]

def matrixExpect {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → Matrix State State ℝ) : Matrix State State ℝ := fun q r =>
  P.expect (fun x => f x q r)

theorem matrixExpect_const {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (K : Matrix State State ℝ) : matrixExpect P (fun _ => K) = K := by
  ext q r
  exact Spin.FinPMF.expect_const _ _

theorem matrixExpect_mul_left {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (K : Matrix State State ℝ) (f : α → Matrix State State ℝ) :
    matrixExpect P (fun x => K * f x) = K * matrixExpect P f := by
  ext q r
  simp only [matrixExpect, Matrix.mul_apply, Spin.FinPMF.expect_sum, expect_const_mul]

theorem matrixExpect_mul_right {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → Matrix State State ℝ) (K : Matrix State State ℝ) :
    matrixExpect P (fun x => f x * K) = matrixExpect P f * K := by
  ext q r
  simp only [matrixExpect, Matrix.mul_apply, Spin.FinPMF.expect_sum, expect_mul_const]

theorem matrixExpect_iid_succ {α : Type*} [Fintype α] (P : Spin.FinPMF α) {a : Nat}
    (f : (Fin (a + 1) → α) → Matrix State State ℝ) :
    matrixExpect (Spin.piPMF (fun _ : Fin (a + 1) => P)) f =
      matrixExpect P (fun x => matrixExpect (Spin.piPMF (fun _ : Fin a => P))
        (fun xs => f (Fin.cons x xs))) := by
  ext q r
  exact expect_iid_succ _ _

def coordinateLaw : Spin.FinPMF (Fin 128) := Spin.FinPMF.uniform _

def uniformImpulseKernel (z : ℝ) : Matrix State State ℝ :=
  matrixExpect coordinateLaw (fun i => roundKernel z {i})

theorem uniformImpulseKernel_apply (z : ℝ) (q r : State) :
    uniformImpulseKernel z q r = (∑ i : Fin 128, roundKernel z {i} q r) / 128 := by
  simp only [uniformImpulseKernel, matrixExpect, coordinateLaw, Spin.FinPMF.expect,
    Spin.FinPMF.uniform, Fintype.card_fin, div_mul_eq_mul_div, one_mul, ← Finset.sum_div]
  norm_num

def uniformImpulseProduct (z : ℝ) : {a : Nat} → (Fin a → Nat) → Nat → Matrix State State ℝ
  | 0, _, last => emptyKernel z last
  | _ + 1, gaps, last => Matrix.of (emptyKernel z (gaps 0)) * uniformImpulseKernel z *
      uniformImpulseProduct z (Fin.tail gaps) last

theorem impulseProduct_average {a : Nat} (z : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
      (fun coords => impulseProduct z gaps coords last) = uniformImpulseProduct z gaps last := by
  induction a with
  | zero => exact matrixExpect_const _ _
  | succ a ih =>
    rw [matrixExpect_iid_succ]
    simp only [impulseProduct, Fin.cons_zero, Fin.tail_cons]
    simp_rw [matrixExpect_mul_left, ih]
    rw [matrixExpect_mul_right, matrixExpect_mul_left]
    rfl

/-- Full finite region product for independently uniform local bit coordinates.
Both the input-coordinate law and the actual transvection law are explicit. -/
theorem endpointKernel_uniform_impulses {a : Nat} (z : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
      (fun coords => endpointKernel z (impulseInputs gaps coords last).get) =
        uniformImpulseProduct z gaps last := by
  simp only [endpointKernel_impulseInputs, impulseProduct_average]

end Spin.Structured.ConcreteEncoder
