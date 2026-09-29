import SpinCodes.Structured.ConcreteFixedFugacity

/-! Exact region-wise composition under independent actual region shuffles. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder ConcreteMaps FiniteKernel Filter
open scoped Topology
attribute [local instance] Classical.propDecidable

theorem expect_pi_succ {α : Type*} [Fintype α] {b : Nat}
    (P : Fin (b + 1) → Spin.FinPMF α) (f : (Fin (b + 1) → α) → ℝ) :
    (Spin.piPMF P).expect f = (P 0).expect (fun x =>
      (Spin.piPMF (Fin.tail P)).expect (fun xs => f (Fin.cons x xs))) := by
  simp only [Spin.FinPMF.expect, Spin.piPMF_apply, sum_tuple_succ,
    Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, Fin.tail, Finset.mul_sum, mul_assoc]

theorem matrixExpect_pi_succ {α : Type*} [Fintype α] {b : Nat}
    (P : Fin (b + 1) → Spin.FinPMF α) (K : (Fin (b + 1) → α) → Matrix State State ℝ) :
    matrixExpect (Spin.piPMF P) K = matrixExpect (P 0) (fun x =>
      matrixExpect (Spin.piPMF (Fin.tail P)) (fun xs => K (Fin.cons x xs))) := by
  ext q r
  exact expect_pi_succ P _

theorem matrixExpect_pi_product {α : Type*} [Fintype α] {b : Nat}
    (P : Fin b → Spin.FinPMF α) (K : Fin b → α → Matrix State State ℝ) :
    matrixExpect (Spin.piPMF P) (fun xs => (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => matrixExpect (P i) (K i))).prod := by
  induction b with
  | zero => simp only [List.ofFn_zero, List.prod_nil, matrixExpect_const]
  | succ b ih =>
    rw [matrixExpect_pi_succ]
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    have he (x : α) :
        matrixExpect (Spin.piPMF (Fin.tail P)) (fun xs : Fin b → α =>
          K 0 x * (List.ofFn (fun i : Fin b => K i.succ (xs i))).prod) =
          K 0 x * (List.ofFn (fun i : Fin b => matrixExpect (P i.succ) (K i.succ))).prod := by
      rw [matrixExpect_mul_left, ih]
      rfl
    simp only [he, matrixExpect_mul_right]

def regionStream {R b : Nat} (regions : Fin b → Finset (Fin (128 * R))) : List Input :=
  (List.ofFn (fun j => List.ofFn (serializedInput (regions j)))).flatten

theorem pathKernel_regionStream {R b : Nat} (z : ℝ) (regions : Fin b → Finset (Fin (128 * R))) :
    pathKernel z (regionStream regions) =
      (List.ofFn (fun j => endpointKernel z (serializedInput (regions j)))).prod := by
  induction b with
  | zero => simp [regionStream, pathKernel]
  | succ b ih =>
    simp only [regionStream, List.ofFn_succ, List.flatten_cons, pathKernel_append, List.prod_cons]
    rw [← endpointKernel_eq_path]
    change endpointKernel z (serializedInput (regions 0)) * pathKernel z (regionStream (Fin.tail regions)) = _
    rw [ih]
    rfl

theorem endpointKernel_regionStream {R b : Nat} (z : ℝ) (regions : Fin b → Finset (Fin (128 * R))) :
    endpointKernel z (regionStream regions).get =
      (List.ofFn (fun j => endpointKernel z (serializedInput (regions j)))).prod := by
  rw [endpointKernel_eq_path, List.ofFn_get, pathKernel_regionStream]

theorem shuffled_regionStream_endpoint {R b : Nat} (θ : ℝ) (S : Fin b → Finset (Fin (128 * R))) :
    matrixExpect (Spin.piPMF (fun j => shuffleLaw (S j)))
      (fun regions => endpointKernel (Real.exp (-(θ / (128 * R)))) (regionStream regions).get) =
        (List.ofFn (fun j => Matrix.of (shuffledRegionKernel θ (S j)))).prod := by
  simp only [endpointKernel_regionStream]
  exact matrixExpect_pi_product (fun j => shuffleLaw (S j)) (fun j T => endpointKernel (Real.exp (-(θ / (128 * R)))) (serializedInput T))

theorem shuffled_regionStream_moment {R b : Nat} (θ : ℝ) (S : Fin b → Finset (Fin (128 * R))) (q : State) :
    (Spin.piPMF (fun j => shuffleLaw (S j))).expect
      (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get q) =
        ∑ r, (List.ofFn (fun j => Matrix.of (shuffledRegionKernel θ (S j)))).prod q r := by
  simp only [← endpointKernel_total, Spin.FinPMF.expect_sum]
  apply sum_congr rfl
  intro r _
  exact congrFun (congrFun (shuffled_regionStream_endpoint θ S) q) r

end Spin.Structured.Placement

