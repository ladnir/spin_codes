import Mathlib.Data.Rat.Cast.Order

namespace Spin.Structured.DenseGeometry

structure Rect where
  alo : ℚ
  ahi : ℚ
  xlo : ℚ
  xhi : ℚ
  deriving DecidableEq

def zeroRect : Rect := ⟨0,0,0,0⟩

inductive Tree where
  | leaf (index : ℕ)
  | alpha (cut : ℚ) (left right : Tree)
  | density (cut : ℚ) (left right : Tree)

def check (boxes : List Rect) : Tree → Rect → Bool
  | .leaf i, r => decide (i < boxes.length ∧
      (boxes.getD i zeroRect).alo ≤ r.alo ∧ r.ahi ≤ (boxes.getD i zeroRect).ahi ∧
      (boxes.getD i zeroRect).xlo ≤ r.xlo ∧ r.xhi ≤ (boxes.getD i zeroRect).xhi)
  | .alpha c l h, r => check boxes l {r with ahi := c} && check boxes h {r with alo := c}
  | .density c l h, r => check boxes l {r with xhi := c} && check boxes h {r with xlo := c}

end Spin.Structured.DenseGeometry
