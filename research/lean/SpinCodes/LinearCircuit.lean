import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.ToLin
import Mathlib.Data.Matrix.Mul

/-!
# Costed binary linear circuits

A dense `m` by `n` macro evaluates its displayed dot products. Its charge
`2*m*n` bounds one coefficient-bit AND and one accumulator XOR per entry.
Composition, disjoint parallel execution, and tiling retain their actual
coordinate semantics. Reindexing denotes wires and has zero gate cost.
Input/output movement can therefore be charged separately by the interface
width. No arbitrary linear map is admitted as a free primitive.
-/
noncomputable section
namespace Spin.LinearCircuit
abbrev Bit := ZMod 2

inductive Circuit : Type → Type → Type 1 where
  | dense {ι κ : Type} [Fintype ι] [Fintype κ] (M : Matrix κ ι Bit) : Circuit ι κ
  | reindex {ι κ : Type} [Fintype ι] [Fintype κ] (e : ι ≃ κ) : Circuit ι κ
  | comp {ι κ υ : Type} [Fintype ι] [Fintype κ] [Fintype υ] (g : Circuit κ υ) (f : Circuit ι κ) : Circuit ι υ
  | parallel {ι κ ι' κ' : Type} [Fintype ι] [Fintype κ] [Fintype ι'] [Fintype κ']
      (f : Circuit ι κ) (g : Circuit ι' κ') :
      Circuit (Sum ι ι') (Sum κ κ')
  | tiled {σ ι κ : Type} [Fintype σ] [Fintype ι] [Fintype κ] (f : Circuit ι κ) :
      Circuit (σ × ι) (σ × κ)

namespace Circuit

def eval : {ι κ : Type} → Circuit ι κ → (ι → Bit) → (κ → Bit)
  | _, _, @Circuit.dense _ _ _ _ M, x => M.mulVec x
  | _, _, @Circuit.reindex _ _ _ _ e, x => fun j => x (e.symm j)
  | _, _, @Circuit.comp _ _ _ _ _ _ g f, x => eval g (eval f x)
  | _, _, @Circuit.parallel _ _ _ _ _ _ _ _ f g, x => Sum.elim
      (eval f (fun i => x (.inl i))) (eval g (fun i => x (.inr i)))
  | _, _, @Circuit.tiled _ _ _ _ _ _ f, x => fun j => eval f (fun i => x (j.1, i)) j.2

def size : {ι κ : Type} → Circuit ι κ → ℕ
  | _, _, @Circuit.dense ι κ _ _ _ => 2 * Fintype.card ι * Fintype.card κ
  | _, _, @Circuit.reindex _ _ _ _ _ => 0
  | _, _, @Circuit.comp _ _ _ _ _ _ g f => size g + size f
  | _, _, @Circuit.parallel _ _ _ _ _ _ _ _ f g => size f + size g
  | _, _, @Circuit.tiled σ _ _ _ _ _ f => Fintype.card σ * size f

def transpose : {ι κ : Type} → Circuit ι κ → Circuit κ ι
  | _, _, @Circuit.dense _ _ _ _ M => .dense M.transpose
  | _, _, @Circuit.reindex _ _ _ _ e => .reindex e.symm
  | _, _, @Circuit.comp _ _ _ _ _ _ g f => .comp (transpose f) (transpose g)
  | _, _, @Circuit.parallel _ _ _ _ _ _ _ _ f g => .parallel (transpose f) (transpose g)
  | _, _, @Circuit.tiled _ _ _ _ _ _ f => .tiled (transpose f)

@[simp] theorem size_transpose {ι κ : Type} (c : Circuit ι κ) :
    c.transpose.size = c.size := by
  induction c with
  | dense M => simp [transpose, size, Nat.mul_comm, Nat.mul_left_comm]
  | reindex e => rfl
  | comp g f hg hf => simp only [transpose, size, hg, hf, Nat.add_comm]
  | parallel f g hf hg => simp only [transpose, size, hf, hg]
  | tiled f hf => simp only [transpose, size, hf]

theorem eval_add {ι κ : Type} (c : Circuit ι κ) (x y : ι → Bit) :
    c.eval (x + y) = c.eval x + c.eval y := by
  induction c with
  | dense M => exact Matrix.mulVec_add M x y
  | reindex e => rfl
  | comp g f hg hf => simp only [eval, hf, hg]
  | parallel f g hf hg =>
    funext j
    cases j with
    | inl j => exact congrFun (hf (fun i => x (.inl i)) (fun i => y (.inl i))) j
    | inr j => exact congrFun (hg (fun i => x (.inr i)) (fun i => y (.inr i))) j
  | tiled f hf =>
    funext j
    exact congrFun (hf (fun i => x (j.1, i)) (fun i => y (j.1, i))) j.2

@[simp] theorem eval_zero {ι κ : Type} (c : Circuit ι κ) : c.eval 0 = 0 := by
  have h := c.eval_add 0 0
  simp only [add_zero] at h
  exact add_left_cancel (h.symm.trans (add_zero _).symm)

def linearMap {ι κ : Type} (c : Circuit ι κ) : (ι → Bit) →ₗ[Bit] (κ → Bit) where
  toFun := c.eval
  map_add' := c.eval_add
  map_smul' a x := by
    have ha : a = 0 ∨ a = 1 := (by decide : ∀ a : Bit, a = 0 ∨ a = 1) a
    rcases ha with rfl | rfl <;> simp

/-- Dense implementation of a linear map, charged for its full matrix. -/
def ofLinearMap {ι κ : Type} [Fintype ι] [Fintype κ] [DecidableEq ι]
    (f : (ι → Bit) →ₗ[Bit] (κ → Bit)) : Circuit ι κ :=
  .dense (LinearMap.toMatrix' f)

@[simp] theorem eval_ofLinearMap {ι κ : Type} [Fintype ι] [Fintype κ] [DecidableEq ι]
    (f : (ι → Bit) →ₗ[Bit] (κ → Bit)) (x : ι → Bit) :
    (ofLinearMap f).eval x = f x :=
  LinearMap.toMatrix'_mulVec f x

@[simp] theorem size_ofLinearMap {ι κ : Type} [Fintype ι] [Fintype κ] [DecidableEq ι]
    (f : (ι → Bit) →ₗ[Bit] (κ → Bit)) :
    (ofLinearMap f).size = 2 * Fintype.card ι * Fintype.card κ := rfl
/-- Structural transposition computes the binary adjoint of the same circuit. -/
theorem transpose_pairing {ι κ : Type} [fi : Fintype ι] [fk : Fintype κ]
    (c : Circuit ι κ) (x : ι → Bit) (y : κ → Bit) :
    (c.transpose.eval y) ⬝ᵥ x = y ⬝ᵥ c.eval x := by
  induction c generalizing fi fk with
  | @dense a b ia ib M =>
    cases Subsingleton.elim fi ia
    cases Subsingleton.elim fk ib
    simp only [transpose, eval]
    rw [dotProduct_comm]
    exact Matrix.dotProduct_transpose_mulVec M x y
  | @reindex a b ia ib e =>
    cases Subsingleton.elim fi ia
    cases Subsingleton.elim fk ib
    change (∑ i, y (e i) * x i) = ∑ j, y j * x (e.symm j)
    simpa only [Equiv.symm_apply_apply] using
      (Equiv.sum_comp e (fun j => y j * x (e.symm j)))
  | @comp a b d ia ib id g f hg hf =>
    cases Subsingleton.elim fi ia
    cases Subsingleton.elim fk id
    change f.transpose.eval (g.transpose.eval y) ⬝ᵥ x = y ⬝ᵥ g.eval (f.eval x)
    rw [hf, hg]
  | @parallel a b a' b' ia ib ia' ib' f g hf hg =>
    cases Subsingleton.elim fi (@instFintypeSum a a' ia ia')
    cases Subsingleton.elim fk (@instFintypeSum b b' ib ib')
    simp only [transpose, eval, dotProduct, Fintype.sum_sum_type, Sum.elim_inl,
      Sum.elim_inr]
    exact congrArg₂ (· + ·) (hf _ _) (hg _ _)
  | @tiled s a b is ia ib f hf =>
    cases Subsingleton.elim fi (@instFintypeProd s a is ia)
    cases Subsingleton.elim fk (@instFintypeProd s b is ib)
    simp only [transpose, eval, dotProduct, Fintype.sum_prod_type]
    apply Finset.sum_congr rfl
    intro a _
    exact hf (fun i => x (a, i)) (fun j => y (a, j))
#print axioms size_transpose
#print axioms transpose_pairing
end Circuit
end Spin.LinearCircuit
