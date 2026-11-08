import SpinCodes.LinearCircuit

/-!
# Scalar expansion of costed binary linear circuits

Each leaf is an explicit expression in zero/input wires, coefficient ANDs,
and XOR gates. Composition connects one component's output wires to the next
component's input wires; it does not substitute expression trees or duplicate
an earlier component. Wires and zero constants are free in this gate model.
Theorems below establish the scalar expansion's semantics and gate count.
They make no claim about the wall-clock evaluation of these Lean definitions.
-/
noncomputable section
namespace Spin.LinearCircuit

inductive ScalarExpr (ι : Type) where
  | zero : ScalarExpr ι
  | input (i : ι) : ScalarExpr ι
  | coeffAnd (a : Bit) (e : ScalarExpr ι) : ScalarExpr ι
  | xor (e f : ScalarExpr ι) : ScalarExpr ι

namespace ScalarExpr

def eval {ι : Type} : ScalarExpr ι → (ι → Bit) → Bit
  | .zero, _ => 0
  | .input i, x => x i
  | .coeffAnd a e, x => a * e.eval x
  | .xor e f, x => e.eval x + f.eval x

def gates {ι : Type} : ScalarExpr ι → ℕ
  | .zero => 0
  | .input _ => 0
  | .coeffAnd _ e => 1 + e.gates
  | .xor e f => 1 + e.gates + f.gates

/-- Truth-table interpretation of the two scalar gates as bit operations. -/
def bitValue (b : Bool) : Bit := if b then 1 else 0

@[simp] theorem bitValue_and (a b : Bool) :
    bitValue (a && b) = bitValue a * bitValue b := by
  cases a <;> cases b <;> decide

@[simp] theorem bitValue_xor (a b : Bool) :
    bitValue (Bool.xor a b) = bitValue a + bitValue b := by
  cases a <;> cases b <;> decide

@[simp] theorem bitValue_decode (a : Bit) : bitValue (decide (a ≠ 0)) = a := by
  exact (by decide : ∀ a : Bit, bitValue (decide (a ≠ 0)) = a) a

def evalBool {ι : Type} : ScalarExpr ι → (ι → Bool) → Bool
  | .zero, _ => false
  | .input i, x => x i
  | .coeffAnd a e, x => decide (a ≠ 0) && e.evalBool x
  | .xor e f, x => Bool.xor (e.evalBool x) (f.evalBool x)

theorem evalBool_agrees {ι : Type} (e : ScalarExpr ι) (x : ι → Bool) :
    bitValue (e.evalBool x) = e.eval (fun i => bitValue (x i)) := by
  induction e with
  | zero => rfl
  | input i => rfl
  | coeffAnd a e ih => simp only [evalBool, bitValue_and, bitValue_decode, ih, eval]
  | xor e f he hf => simp only [evalBool, bitValue_xor, he, hf, eval]
/-- One coefficient AND and one accumulator XOR per input coordinate. -/
def row {ι : Type} (a : ι → Bit) : List ι → ScalarExpr ι
  | [] => .zero
  | i :: is => .xor (.coeffAnd (a i) (.input i)) (row a is)

theorem row_eval {ι : Type} (a x : ι → Bit) (is : List ι) :
    (row a is).eval x = (is.map (fun i => a i * x i)).sum := by
  induction is with
  | nil => rfl
  | cons i is ih => simp only [row, eval, List.map_cons, List.sum_cons, ih]

theorem row_gates {ι : Type} (a : ι → Bit) (is : List ι) :
    (row a is).gates = 2 * is.length := by
  induction is with
  | nil => rfl
  | cons i is ih => simp [row, gates, ih, Nat.mul_add, Nat.add_comm]

end ScalarExpr

/-- A hierarchical circuit of explicit scalar gate expressions. Intermediate
outputs in `comp` are shared wires; components are never syntactically inlined
into each individual output expression. -/
inductive PrimitiveCircuit : Type → Type → Type 1 where
  | rows {ι κ : Type} [Fintype ι] [Fintype κ]
      (es : κ → ScalarExpr ι) : PrimitiveCircuit ι κ
  | reindex {ι κ : Type} [Fintype ι] [Fintype κ] (e : ι ≃ κ) : PrimitiveCircuit ι κ
  | comp {ι κ υ : Type} [Fintype ι] [Fintype κ] [Fintype υ]
      (g : PrimitiveCircuit κ υ) (f : PrimitiveCircuit ι κ) : PrimitiveCircuit ι υ
  | parallel {ι κ ι' κ' : Type} [Fintype ι] [Fintype κ] [Fintype ι'] [Fintype κ']
      (f : PrimitiveCircuit ι κ) (g : PrimitiveCircuit ι' κ') :
      PrimitiveCircuit (Sum ι ι') (Sum κ κ')
  | tiled {σ ι κ : Type} [Fintype σ] [Fintype ι] [Fintype κ]
      (f : PrimitiveCircuit ι κ) : PrimitiveCircuit (σ × ι) (σ × κ)

namespace PrimitiveCircuit

def eval : {ι κ : Type} → PrimitiveCircuit ι κ → (ι → Bit) → (κ → Bit)
  | _, _, @PrimitiveCircuit.rows _ _ _ _ es, x => fun j => (es j).eval x
  | _, _, @PrimitiveCircuit.reindex _ _ _ _ e, x => fun j => x (e.symm j)
  | _, _, @PrimitiveCircuit.comp _ _ _ _ _ _ g f, x => eval g (eval f x)
  | _, _, @PrimitiveCircuit.parallel _ _ _ _ _ _ _ _ f g, x => Sum.elim
      (eval f (fun i => x (.inl i))) (eval g (fun i => x (.inr i)))
  | _, _, @PrimitiveCircuit.tiled _ _ _ _ _ _ f, x => fun j => eval f (fun i => x (j.1, i)) j.2

def gates : {ι κ : Type} → PrimitiveCircuit ι κ → ℕ
  | _, _, @PrimitiveCircuit.rows _ _ _ _ es => ∑ j, (es j).gates
  | _, _, @PrimitiveCircuit.reindex _ _ _ _ _ => 0
  | _, _, @PrimitiveCircuit.comp _ _ _ _ _ _ g f => gates g + gates f
  | _, _, @PrimitiveCircuit.parallel _ _ _ _ _ _ _ _ f g => gates f + gates g
  | _, _, @PrimitiveCircuit.tiled σ _ _ _ _ _ f => Fintype.card σ * gates f

end PrimitiveCircuit

namespace Circuit

def lower : {ι κ : Type} → Circuit ι κ → PrimitiveCircuit ι κ
  | _, _, @Circuit.dense ι _ _ _ M =>
      .rows (fun j => ScalarExpr.row (M j) (Finset.univ.toList : List ι))
  | _, _, @Circuit.reindex _ _ _ _ e => .reindex e
  | _, _, @Circuit.comp _ _ _ _ _ _ g f => .comp g.lower f.lower
  | _, _, @Circuit.parallel _ _ _ _ _ _ _ _ f g => .parallel f.lower g.lower
  | _, _, @Circuit.tiled _ _ _ _ _ _ f => .tiled f.lower

theorem lower_eval {ι κ : Type} (c : Circuit ι κ) (x : ι → Bit) :
    c.lower.eval x = c.eval x := by
  induction c with
  | dense M =>
    funext j
    simp only [lower, PrimitiveCircuit.eval, ScalarExpr.row_eval, Circuit.eval,
      Matrix.mulVec, dotProduct]
    exact Finset.sum_map_toList _ _
  | reindex e => rfl
  | comp g f hg hf => simp only [lower, PrimitiveCircuit.eval, Circuit.eval, hf, hg]
  | parallel f g hf hg =>
    funext j
    cases j with
    | inl j => exact congrFun (hf (fun i => x (.inl i))) j
    | inr j => exact congrFun (hg (fun i => x (.inr i))) j
  | tiled f hf =>
    funext j
    exact congrFun (hf (fun i => x (j.1, i))) j.2

theorem lower_gates {ι κ : Type} (c : Circuit ι κ) :
    c.lower.gates = c.size := by
  induction c with
  | dense M =>
    simp [lower, PrimitiveCircuit.gates, ScalarExpr.row_gates, Circuit.size,
      Nat.mul_comm, Nat.mul_left_comm]
  | reindex e => rfl
  | comp g f hg hf => simp only [lower, PrimitiveCircuit.gates, Circuit.size, hg, hf]
  | parallel f g hf hg => simp only [lower, PrimitiveCircuit.gates, Circuit.size, hf, hg]
  | tiled f hf => simp only [lower, PrimitiveCircuit.gates, Circuit.size, hf]

/-- The same primitive expansion certifies the actual adjoint circuit. -/
theorem lower_transpose_pairing {ι κ : Type} [Fintype ι] [Fintype κ]
    (c : Circuit ι κ) (x : ι → Bit) (y : κ → Bit) :
    c.transpose.lower.eval y ⬝ᵥ x = y ⬝ᵥ c.lower.eval x := by
  simp only [lower_eval]
  exact c.transpose_pairing x y

@[simp] theorem lower_transpose_gates {ι κ : Type} (c : Circuit ι κ) :
    c.transpose.lower.gates = c.lower.gates := by
  rw [lower_gates, lower_gates, size_transpose]

#print axioms lower_eval
#print axioms lower_gates
#print axioms lower_transpose_pairing
end Circuit
end Spin.LinearCircuit
