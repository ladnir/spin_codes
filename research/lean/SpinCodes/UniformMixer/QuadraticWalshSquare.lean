import SpinCodes.UniformMixer.QuadraticWalsh

/-!
# Walsh squares through quadratic autocorrelation

The identity here is valid for every Boolean function on six variables.
For a quadratic function, its polar derivative is linear in the second
argument; character orthogonality will then reduce the square to the
size and character of the radical.
-/

namespace Spin.UniformMixer

def polar6 (f : State6 → ZMod 2) (z x : State6) : ZMod 2 :=
  f (x + z) + f x + f z

theorem binarySign_shift_identity (u v w s t : ZMod 2) :
    binarySign (u + v) * binarySign (w + (v + t)) =
      binarySign (s + t) * binarySign (w + u + s) := by
  fin_cases u <;> fin_cases v <;> fin_cases w <;>
    fin_cases s <;> fin_cases t <;> decide

theorem walsh6_square_autocorrelation
    (f : State6 → ZMod 2) (a : State6) :
    walsh6 f a * walsh6 f a =
      ∑ z : State6, binarySign (f z + dot6Linear z a) *
        (∑ x : State6, binarySign (polar6 f z x)) := by
  classical
  have hshift (x : State6) :
      (∑ y : State6, binarySign (f y + dot6Linear y a)) =
        ∑ z : State6,
          binarySign (f (x + z) + dot6Linear (x + z) a) := by
    symm
    apply Fintype.sum_equiv (Equiv.addRight x)
    intro z
    change binarySign (f (x + z) + dot6Linear (x + z) a) =
      binarySign (f (z + x) + dot6Linear (z + x) a)
    rw [add_comm x z]
  calc
    _ = ∑ x : State6, ∑ y : State6,
        binarySign (f x + dot6Linear x a) *
          binarySign (f y + dot6Linear y a) := by
      unfold walsh6
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro x _
      rw [Finset.mul_sum]
    _ = ∑ x : State6, ∑ z : State6,
        binarySign (f x + dot6Linear x a) *
          binarySign (f (x + z) + dot6Linear (x + z) a) := by
      apply Finset.sum_congr rfl
      intro x _
      rw [← Finset.mul_sum, hshift x, Finset.mul_sum]
    _ = ∑ z : State6, binarySign (f z + dot6Linear z a) *
          (∑ x : State6, binarySign (polar6 f z x)) := by
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro z _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro x _
      rw [dot6Linear_add]
      change binarySign (f x + dot6Linear x a) *
          binarySign (f (x + z) +
            (dot6Linear x a + dot6Linear z a)) =
        binarySign (f z + dot6Linear z a) *
          binarySign (f (x + z) + f x + f z)
      exact binarySign_shift_identity
        (f x) (dot6Linear x a) (f (x + z)) (f z) (dot6Linear z a)

end Spin.UniformMixer
