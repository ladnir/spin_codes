import SpinCodes.Structured.ConcreteFixedInsertionPaths
import SpinCodes.Structured.ConcreteFixedInsertionEnlarged
import SpinCodes.Structured.ConcreteFixedLargeSpacings

/-! Gap exponentials in each matrix state path are the actual live-spacing duration. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder Matrix
attribute [local instance] Classical.propDecidable

def pathLive {Q : Nat} (i : Fin 2) (ys : Fin Q → Fin 2) : Finset (Fin (Q+1)) :=
  univ.filter (fun k => (Fin.cons i ys : Fin (Q+1) → Fin 2) k = 1)

def pathWeight {Q : Nat} (M : Fin Q → Matrix (Fin 2) (Fin 2) ℝ) (w : Fin 2 → ℝ)
    (i : Fin 2) (ys : Fin Q → Fin 2) : ℝ :=
  (∏ k : Fin Q, M k ((Fin.cons i ys : Fin (Q+1) → Fin 2) k.castSucc) (ys k)) *
    w ((Fin.cons i ys : Fin (Q+1) → Fin 2) (Fin.last Q))

def exponentialDiagonal {Q : Nat} (γ : ℝ) (gaps : Fin (Q+1) → ℝ) (k : Fin (Q+1)) (s : Fin 2) : ℝ :=
  Real.exp (-γ * gaps k * (if s=1 then 1 else 0))

theorem exponentialDiagonal_eq {Q : Nat} (γ : ℝ) (gaps : Fin (Q+1) → ℝ) (k : Fin (Q+1)) :
    Matrix.diagonal (exponentialDiagonal γ gaps k) = timeEmpty γ (gaps k) := by
  ext i j
  fin_cases i <;> fin_cases j <;> norm_num [Matrix.diagonal_apply, exponentialDiagonal, timeEmpty]

theorem fullSiteProduct_diagonalPath {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (x u : Fin Q → ℝ) :
    fullSiteProduct γ P x u = diagonalPathProduct (exponentialDiagonal γ (positionGaps x)) (fun k => 1+u k • P) := by
  rw [diagonalPathProduct_list]
  simp only [exponentialDiagonal_eq]
  rfl

theorem exponentialDiagonal_prod {Q : Nat} (γ : ℝ) (x : Fin Q → ℝ) (i : Fin 2) (ys : Fin Q → Fin 2) :
    (∏ k, exponentialDiagonal γ (positionGaps x) k ((Fin.cons i ys : Fin (Q+1) → Fin 2) k)) =
      Real.exp (-γ * liveDuration (pathLive i ys) x) := by
  unfold exponentialDiagonal
  rw [← Real.exp_sum]
  congr 1
  simp only [liveDuration, pathLive, sum_filter, Finset.mul_sum]
  apply sum_congr rfl
  intro k hk
  split_ifs <;> ring

theorem fullSiteProduct_path_expansion {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (x u : Fin Q → ℝ) (w : Fin 2 → ℝ) (i : Fin 2) :
    (fullSiteProduct γ P x u *ᵥ w) i = ∑ ys : Fin Q → Fin 2,
      Real.exp (-γ * liveDuration (pathLive i ys) x) * pathWeight (fun k => 1+u k • P) w i ys := by
  rw [fullSiteProduct_diagonalPath, diagonalPathProduct_expansion]
  simp only [exponentialDiagonal_prod, pathWeight, mul_assoc]

def constantDiagonal (σ : ℝ) (s : Fin 2) : ℝ := if s=1 then σ else 1

theorem constantDiagonal_eq (σ : ℝ) : Matrix.diagonal (constantDiagonal σ) = Spin.Dsig σ := by
  ext i j
  fin_cases i <;> fin_cases j <;> norm_num [Matrix.diagonal_apply, constantDiagonal, Spin.Dsig]

theorem constantDiagonal_prod {Q : Nat} (σ : ℝ) (i : Fin 2) (ys : Fin Q → Fin 2) :
    (∏ k : Fin (Q+1), constantDiagonal σ ((Fin.cons i ys : Fin (Q+1) → Fin 2) k)) = σ ^ (pathLive i ys).card := by
  simp [constantDiagonal, pathLive, ← prod_filter]

theorem constant_path_expansion {Q : Nat} (σ : ℝ) (M : Fin Q → Matrix (Fin 2) (Fin 2) ℝ)
    (w : Fin 2 → ℝ) (i : Fin 2) :
    (diagonalPathProduct (fun _ : Fin (Q+1) => constantDiagonal σ) M *ᵥ w) i =
      ∑ ys : Fin Q → Fin 2, σ ^ (pathLive i ys).card * pathWeight M w i ys := by
  rw [diagonalPathProduct_expansion]
  simp only [constantDiagonal_prod, pathWeight, mul_assoc]

theorem constant_pathProduct {Q : Nat} (σ r : ℝ) (u : Fin Q → ℝ) :
    diagonalPathProduct (fun _ : Fin (Q+1) => constantDiagonal σ) (fun k => 1+u k • Spin.Pplus r) =
      Spin.Dsig σ * Spin.pathProduct σ r (List.ofFn u) := by
  induction Q with
  | zero => simp [diagonalPathProduct, constantDiagonal_eq, Spin.pathProduct]
  | succ Q ih =>
    simp only [diagonalPathProduct, constantDiagonal_eq, List.ofFn_succ, Spin.pathProduct, Spin.IuP_eq]
    change Spin.Dsig σ * (1+u 0 • Spin.Pplus r) *
      diagonalPathProduct (fun _ : Fin (Q+1) => constantDiagonal σ) (fun k : Fin Q => 1+u k.succ • Spin.Pplus r) = _
    rw [ih]
    simp only [Matrix.mul_assoc]

end Spin.Structured.Placement

