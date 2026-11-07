import SpinCodes.Turbo.ColumnMixer

/-!
# Coordinates of a uniformly randomized nonzero column

The old Scalable proof uses a checked transitive action on nonzero
states. Here a nonzero field scalar acts transitively on each nonzero
outer column. A zero-preserving coordinate equivalence transfers this
law to the 32 binary coordinates used to form four-bit packets. The
coordinate equivalence itself is a separate implementation bridge.
-/

noncomputable section

namespace Spin.Turbo

variable {F B : Type*} [Field F] [Fintype F] [DecidableEq F]
  [Zero B] [Fintype B] [DecidableEq B]

/-- Restrict a zero-preserving coordinate bijection to nonzero values. -/
def nonzeroCoordinateEquiv (e : F ≃ B) (hz : e 0 = 0) :
    {x : F // x ≠ 0} ≃ {y : B // y ≠ 0} where
  toFun x := ⟨e x.val, by
    intro he
    exact x.property (e.injective (he.trans hz.symm))⟩
  invFun y := ⟨e.symm y.val, by
    intro he
    apply y.property
    calc
      y.val = e (e.symm y.val) := (e.apply_symm_apply _).symm
      _ = e 0 := by rw [he]
      _ = 0 := hz⟩
  left_inv x := Subtype.ext (e.symm_apply_apply _)
  right_inv y := Subtype.ext (e.apply_symm_apply _)

/-- For any fixed nonzero column, ideal nonzero scalar multiplication
produces a uniform nonzero binary-coordinate column. -/
theorem uniform_nonzero_coordinates
    [Nonempty {y : B // y ≠ 0}]
    (e : F ≃ B) (hz : e 0 = 0)
    (x : F) (hx : x ≠ 0) :
    (FinPMF.uniform {c : F // c ≠ 0}).map
      (fun c => nonzeroCoordinateEquiv e hz (nonzeroMulEquiv x hx c)) =
      FinPMF.uniform {y : B // y ≠ 0} := by
  change (FinPMF.uniform {c : F // c ≠ 0}).map
    ((nonzeroMulEquiv x hx).trans (nonzeroCoordinateEquiv e hz)) = _
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply]
  change 1 / (Fintype.card {c : F // c ≠ 0} : ℝ) =
    1 / (Fintype.card {y : B // y ≠ 0} : ℝ)
  rw [Fintype.card_congr
    ((nonzeroMulEquiv x hx).trans (nonzeroCoordinateEquiv e hz))]

#print axioms nonzeroCoordinateEquiv
#print axioms uniform_nonzero_coordinates

end Spin.Turbo
