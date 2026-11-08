import SpinCodes.Turbo.ColumnCoordinates

/-!
# Uniform action of adjoint field multiplication

This algebraic bridge is conditional on a nondegenerate binary pairing and
the defining adjoint identity. The concrete field basis and transposed
kernel must be checked against those hypotheses separately.
-/

noncomputable section

namespace Spin.Turbo

variable {F B : Type*} [Field F] [Fintype F] [DecidableEq F]
  [AddCommGroup B] [Fintype B] [DecidableEq B]

omit [Fintype F] [DecidableEq F] [Fintype B] [DecidableEq B] in
private theorem primal_surjective (e : F ≃+ B) (c : F) (hc : c ≠ 0) :
    Function.Surjective (fun y : B => e (c * e.symm y)) := by
  intro z
  refine ⟨e (c⁻¹ * e.symm z), ?_⟩
  simp [hc]

/-- Adjoint multiplication by a nonzero scalar has no nonzero kernel. -/
theorem adjoint_ne_zero
    (e : F ≃+ B) (pair : B → B →+ ZMod 2)
    (nondeg : ∀ x : B, (∀ y : B, pair x y = 0) → x = 0)
    (zero_left : ∀ y : B, pair 0 y = 0)
    (adj : F → B → B)
    (hadj : ∀ (c : F) (x y : B), pair (adj c x) y =
      pair x (e (c * e.symm y)))
    {c : F} {x : B} (hc : c ≠ 0) (hx : x ≠ 0) : adj c x ≠ 0 := by
  intro hzero
  apply hx
  apply nondeg
  intro z
  obtain ⟨y, hy⟩ := primal_surjective e c hc z
  rw [← hy]
  change pair x (e (c * e.symm y)) = 0
  rw [← hadj, hzero, zero_left]

/-- Distinct coefficients give distinct adjoint images of a fixed nonzero
column. This is the key bijection behind the uniform-column law. -/
theorem adjoint_coeff_injective
    (e : F ≃+ B) (pair : B → B →+ ZMod 2)
    (nondeg : ∀ x : B, (∀ y : B, pair x y = 0) → x = 0)
    (adj : F → B → B)
    (hadj : ∀ (c : F) (x y : B), pair (adj c x) y =
      pair x (e (c * e.symm y)))
    {x : B} (hx : x ≠ 0) : Function.Injective (fun c : F => adj c x) := by
  intro c d hcd
  change adj c x = adj d x at hcd
  by_contra hne
  have hdiff : c - d ≠ 0 := sub_ne_zero.mpr hne
  have hzero : x = 0 := by
    apply nondeg
    intro z
    obtain ⟨y, hy⟩ := primal_surjective e (c - d) hdiff z
    rw [← hy]
    change pair x (e ((c - d) * e.symm y)) = 0
    have heq : pair x (e (c * e.symm y)) =
        pair x (e (d * e.symm y)) := by
      calc
        _ = pair (adj c x) y := (hadj c x y).symm
        _ = pair (adj d x) y := by rw [hcd]
        _ = pair x (e (d * e.symm y)) := hadj d x y
    rw [sub_mul, e.map_sub, map_sub, heq, sub_self]
  exact hx hzero

/-- Restrict the adjoint action of a fixed nonzero column to a permutation
of nonzero columns. -/
def nonzeroAdjointEquiv
    (e : F ≃+ B) (pair : B → B →+ ZMod 2)
    (nondeg : ∀ x : B, (∀ y : B, pair x y = 0) → x = 0)
    (zero_left : ∀ y : B, pair 0 y = 0)
    (adj : F → B → B)
    (hadj : ∀ (c : F) (x y : B), pair (adj c x) y =
      pair x (e (c * e.symm y)))
    (x : B) (hx : x ≠ 0) :
    {c : F // c ≠ 0} ≃ {y : B // y ≠ 0} := by
  let f : {c : F // c ≠ 0} → {y : B // y ≠ 0} :=
    fun c => ⟨adj c.1 x, adjoint_ne_zero e pair nondeg zero_left adj hadj c.2 hx⟩
  have hinj : Function.Injective f := by
    intro c d h
    apply Subtype.ext
    exact adjoint_coeff_injective e pair nondeg adj hadj hx
      (congrArg Subtype.val h)
  exact Equiv.ofBijective f ⟨hinj,
    hinj.surjective_of_finite (nonzeroCoordinateEquiv e.toEquiv e.map_zero)⟩

/-- Uniform nonzero coefficients produce a uniform nonzero adjoint image. -/
theorem uniform_nonzero_adjoint
    [Nonempty {y : B // y ≠ 0}]
    (e : F ≃+ B) (pair : B → B →+ ZMod 2)
    (nondeg : ∀ x : B, (∀ y : B, pair x y = 0) → x = 0)
    (zero_left : ∀ y : B, pair 0 y = 0)
    (adj : F → B → B)
    (hadj : ∀ (c : F) (x y : B), pair (adj c x) y =
      pair x (e (c * e.symm y)))
    (x : B) (hx : x ≠ 0) :
    (FinPMF.uniform {c : F // c ≠ 0}).map
      (nonzeroAdjointEquiv e pair nondeg zero_left adj hadj x hx) =
      FinPMF.uniform {y : B // y ≠ 0} := by
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply]
  change 1 / (Fintype.card {c : F // c ≠ 0} : ℝ) =
    1 / (Fintype.card {y : B // y ≠ 0} : ℝ)
  rw [Fintype.card_congr
    (nonzeroAdjointEquiv e pair nondeg zero_left adj hadj x hx)]

/-- Independent sampled coefficients give independent uniform nonzero
adjoint images for any fixed family of nonzero input columns. -/
theorem uniform_active_adjoint_columns
    [Nonempty {y : B // y ≠ 0}]
    {n : ℕ}
    (e : F ≃+ B) (pair : B → B →+ ZMod 2)
    (nondeg : ∀ x : B, (∀ y : B, pair x y = 0) → x = 0)
    (zero_left : ∀ y : B, pair 0 y = 0)
    (adj : F → B → B)
    (hadj : ∀ (c : F) (x y : B), pair (adj c x) y =
      pair x (e (c * e.symm y)))
    (x : Fin n → B) (hx : ∀ i, x i ≠ 0) :
    (piPMF (fun _ : Fin n => FinPMF.uniform {c : F // c ≠ 0})).map
      (fun coeff i =>
        nonzeroAdjointEquiv e pair nondeg zero_left adj hadj
          (x i) (hx i) (coeff i)) =
      piPMF (fun _ : Fin n =>
        FinPMF.uniform {y : B // y ≠ 0}) := by
  calc
    _ = piPMF (fun i : Fin n =>
        (FinPMF.uniform {c : F // c ≠ 0}).map
          (nonzeroAdjointEquiv e pair nondeg zero_left adj hadj
            (x i) (hx i))) :=
      piPMF_map (fun _ : Fin n => FinPMF.uniform {c : F // c ≠ 0})
        (fun i => nonzeroAdjointEquiv e pair nondeg zero_left adj hadj
          (x i) (hx i))
    _ = _ := by simp_rw [uniform_nonzero_adjoint]

#print axioms adjoint_ne_zero
#print axioms adjoint_coeff_injective
#print axioms uniform_nonzero_adjoint
#print axioms uniform_active_adjoint_columns

end Spin.Turbo
