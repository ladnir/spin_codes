import SpinCodes.FiniteProductLaw

/-!
# Uniform scalar action on a nonzero outer column

This follows the same checked-action pattern as
`UniformMixer.ScalableStepLaw.uniform_state_map_update`: a uniform
invertible choice acts transitively on every fixed nonzero input.
The theorem here is field-generic. A later bridge must identify the
binary coordinates and the adjoint action used by Turbo SPIN.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F] [Fintype F] [DecidableEq F]

instance : Nonempty {c : F // c ≠ 0} := ⟨⟨1, one_ne_zero⟩⟩

/-- Multiplication by a fixed nonzero field element permutes nonzero scalars. -/
def nonzeroMulEquiv (x : F) (hx : x ≠ 0) :
    {c : F // c ≠ 0} ≃ {y : F // y ≠ 0} where
  toFun c := ⟨c.1 * x, mul_ne_zero c.2 hx⟩
  invFun y := ⟨y.1 / x, div_ne_zero y.2 hx⟩
  left_inv := by
    intro c
    apply Subtype.ext
    simp [hx]
  right_inv := by
    intro y
    apply Subtype.ext
    simp [hx]

/-- The output of a uniform nonzero scalar is uniform for each nonzero input. -/
theorem uniform_nonzero_mul (x : F) (hx : x ≠ 0) :
    (FinPMF.uniform {c : F // c ≠ 0}).map (nonzeroMulEquiv x hx) =
      FinPMF.uniform {y : F // y ≠ 0} := by
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply]
  rfl

/-- Independently sampled nonzero coefficients act on all columns. -/
def scaledColumns {n : ℕ} (x : Fin n → F)
    (coeff : Fin n → {c : F // c ≠ 0}) : Fin n → F :=
  fun i => (coeff i).1 * x i

def scaledColumnLaw {n : ℕ} (x : Fin n → F) : FinPMF (Fin n → F) :=
  (piPMF (fun _ : Fin n => FinPMF.uniform {c : F // c ≠ 0})).map
    (scaledColumns x)

/-- The complete column law is the product of the individual laws. -/
theorem scaledColumnLaw_product {n : ℕ} (x : Fin n → F) :
    scaledColumnLaw x =
      piPMF (fun i : Fin n =>
        (FinPMF.uniform {c : F // c ≠ 0}).map
          (fun c => c.1 * x i)) := by
  exact piPMF_map (fun _ : Fin n =>
    FinPMF.uniform {c : F // c ≠ 0})
    (fun i c => c.1 * x i)

omit [Fintype F] [DecidableEq F] in
/-- Nonzero scalar maps preserve the exact set of nonzero columns. -/
theorem scaledColumns_zero_iff {n : ℕ} (x : Fin n → F)
    (coeff : Fin n → {c : F // c ≠ 0}) (i : Fin n) :
    scaledColumns x coeff i = 0 ↔ x i = 0 := by
  simp [scaledColumns, (coeff i).2]

#print axioms nonzeroMulEquiv
#print axioms uniform_nonzero_mul
#print axioms scaledColumnLaw_product
#print axioms scaledColumns_zero_iff

end Spin.Turbo
