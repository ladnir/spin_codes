import SpinCodes.UniformMixer.ScalableRegionComposition
import SpinCodes.Structured.ConcreteFixedRegionComposition

/-!
# Averaging independent region shuffles

After the independently shuffled outer rows are transposed, the route
shuffles each region independently. The full-state regional kernels
therefore compose after expectation, in their encoder order.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute Spin.Structured.Routing
  Spin.Structured.ConcreteEncoder

def matrixExpect22 {α : Type*} [Fintype α]
    (P : Spin.FinPMF α) (f : α → Matrix State22 State22 ℝ) :
    Matrix State22 State22 ℝ :=
  fun q r => P.expect (fun x => f x q r)

theorem matrixExpect22_const {α : Type*} [Fintype α]
    (P : Spin.FinPMF α) (K : Matrix State22 State22 ℝ) :
    matrixExpect22 P (fun _ => K) = K := by
  ext q r
  exact Spin.FinPMF.expect_const _ _

theorem matrixExpect22_mul_left {α : Type*} [Fintype α]
    (P : Spin.FinPMF α) (K : Matrix State22 State22 ℝ)
    (f : α → Matrix State22 State22 ℝ) :
    matrixExpect22 P (fun x => K * f x) = K * matrixExpect22 P f := by
  ext q r
  simp only [matrixExpect22, Matrix.mul_apply, Spin.FinPMF.expect_sum,
    expect_const_mul]

theorem matrixExpect22_mul_right {α : Type*} [Fintype α]
    (P : Spin.FinPMF α) (f : α → Matrix State22 State22 ℝ)
    (K : Matrix State22 State22 ℝ) :
    matrixExpect22 P (fun x => f x * K) = matrixExpect22 P f * K := by
  ext q r
  simp only [matrixExpect22, Matrix.mul_apply, Spin.FinPMF.expect_sum,
    expect_mul_const]

theorem matrixExpect22_pi_succ {α : Type*} [Fintype α] {b : ℕ}
    (P : Fin (b + 1) → Spin.FinPMF α)
    (K : (Fin (b + 1) → α) → Matrix State22 State22 ℝ) :
    matrixExpect22 (Spin.piPMF P) K =
      matrixExpect22 (P 0) (fun x =>
        matrixExpect22 (Spin.piPMF (Fin.tail P))
          (fun xs => K (Fin.cons x xs))) := by
  ext q r
  exact Spin.Structured.Placement.expect_pi_succ P _

theorem matrixExpect22_pi_product {α : Type*} [Fintype α] {b : ℕ}
    (P : Fin b → Spin.FinPMF α)
    (K : Fin b → α → Matrix State22 State22 ℝ) :
    matrixExpect22 (Spin.piPMF P)
      (fun xs => (List.ofFn (fun i => K i (xs i))).prod) =
    (List.ofFn (fun i => matrixExpect22 (P i) (K i))).prod := by
  induction b with
  | zero => simp only [List.ofFn_zero, List.prod_nil, matrixExpect22_const]
  | succ b ih =>
      rw [matrixExpect22_pi_succ]
      simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
      have he (x : α) :
          matrixExpect22 (Spin.piPMF (Fin.tail P)) (fun xs : Fin b → α =>
            K 0 x * (List.ofFn (fun i : Fin b => K i.succ (xs i))).prod) =
          K 0 x * (List.ofFn (fun i : Fin b =>
            matrixExpect22 (P i.succ) (K i.succ))).prod := by
        rw [matrixExpect22_mul_left, ih]
        rfl
      simp only [he, matrixExpect22_mul_right]

theorem kernelListProduct_matrix
    (Ks : List (State22 → State22 → ℝ)) :
    kernelListProduct Ks = (Ks.map Matrix.of).prod := by
  induction Ks with
  | nil =>
      ext q r
      simp [kernelListProduct, identity, Matrix.one_apply]
  | cons K tail ih =>
      change compose K (kernelListProduct tail) =
        Matrix.of K * (tail.map Matrix.of).prod
      rw [ih]
      exact compose_eq_matrixMul (Matrix.of K) ((tail.map Matrix.of).prod)

def scalableAveragedRegion (m : ℕ) (z : ℝ)
    (S : Finset (Fin (Lsched m))) : Matrix State22 State22 ℝ :=
  matrixExpect22 (shuffleLaw S) (fun T =>
    Matrix.of (chronologicalRegionKernel z (scalableRegionSupport m T)))

theorem scalable_region_shuffle_product (m : ℕ) (z : ℝ)
    (S : Fin (bsched m) → Finset (Fin (Lsched m))) :
    matrixExpect22 (regionKernel S) (fun regions =>
      Matrix.of (inputListKernel z
        (List.ofFn (scalableFlatInput m regions)))) =
    (List.ofFn (fun j : Fin (bsched m) =>
      scalableAveragedRegion m z (S j))).prod := by
  have hpoint (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
      Matrix.of (inputListKernel z
        (List.ofFn (scalableFlatInput m regions))) =
      (List.ofFn (fun j : Fin (bsched m) =>
        Matrix.of (chronologicalRegionKernel z
          (scalableRegionSupport m (regions j))))).prod := by
    rw [scalable_kernel_region_product, kernelListProduct_matrix]
    simp only [List.map_ofFn]
    rfl
  simp only [regionKernel, hpoint, scalableAveragedRegion]
  exact matrixExpect22_pi_product
    (fun j => shuffleLaw (S j))
    (fun _ T => Matrix.of (chronologicalRegionKernel z
      (scalableRegionSupport m T)))

end Spin.UniformMixer
