import SpinCodes.UniformMixer.QuadraticWalshSquare

/-!
# Polar derivatives of six-variable quadratic forms

Each quadratic part has fifteen independent coefficients, one per pair of
variables. Its polar derivative is linear in each argument. This is the
algebraic step that makes the Walsh-square autocorrelation vanish away
from the radical.
-/

namespace Spin.UniformMixer

def pairLeft : Fin 15 → Fin 6 := ![
  0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4]

def pairRight : Fin 15 → Fin 6 := ![
  1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 4, 5, 5]

def quadraticPart (c : Fin 15 → ZMod 2) (x : State6) : ZMod 2 :=
  ∑ p : Fin 15, c p * x (pairLeft p) * x (pairRight p)

theorem quadraticPart_zero (c : Fin 15 → ZMod 2) :
    quadraticPart c 0 = 0 := by
  simp [quadraticPart]

def quadraticPolarLinear (c : Fin 15 → ZMod 2) (z : State6) :
    State6 →ₗ[ZMod 2] ZMod 2 where
  toFun x := ∑ p : Fin 15,
    c p * (x (pairLeft p) * z (pairRight p) +
      z (pairLeft p) * x (pairRight p))
  map_add' x y := by
    simp [Pi.add_apply, mul_add, add_mul, Finset.sum_add_distrib]
    ring
  map_smul' a x := by
    simp only [Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro p _
    simp only [RingHom.id_apply]
    ring

theorem quadraticPart_polar (c : Fin 15 → ZMod 2)
    (z x : State6) :
    polar6 (quadraticPart c) z x = quadraticPolarLinear c z x := by
  simp only [polar6, quadraticPart, quadraticPolarLinear,
    LinearMap.coe_mk, AddHom.coe_mk]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro p _
  simp only [Pi.add_apply]
  ring_nf
  simp [show (2 : ZMod 2) = 0 by decide]

theorem binarySign_sum_state6_linear_zero
    (ℓ : State6 →ₗ[ZMod 2] ZMod 2) (hℓ : ℓ = 0) :
    (∑ x : State6, binarySign (ℓ x)) = 64 := by
  classical
  subst ℓ
  simp [binarySign_zero]

theorem binarySign_sum_state6_linear_nonzero
    (ℓ : State6 →ₗ[ZMod 2] ZMod 2) (hℓ : ℓ ≠ 0) :
    (∑ x : State6, binarySign (ℓ x)) = 0 := by
  have hex : ∃ x : State6, ℓ x ≠ 0 := by
    by_contra h
    push Not at h
    apply hℓ
    apply LinearMap.ext
    intro x
    exact h x
  obtain ⟨x, hx⟩ := hex
  have hx1 : ℓ x = 1 :=
    (zmod2_zero_or_one (ℓ x)).resolve_left hx
  exact binarySign_sum_linear ℓ x hx1

noncomputable def quadraticRadical (c : Fin 15 → ZMod 2) : Finset State6 := by
  classical
  exact Finset.univ.filter (fun z => quadraticPolarLinear c z = 0)

theorem quadraticPart_walsh_square (c : Fin 15 → ZMod 2)
    (a : State6) :
    walsh6 (quadraticPart c) a * walsh6 (quadraticPart c) a =
      64 * ∑ z ∈ quadraticRadical c,
        binarySign (quadraticPart c z + dot6Linear z a) := by
  classical
  rw [walsh6_square_autocorrelation]
  simp_rw [quadraticPart_polar]
  calc
    _ = ∑ z : State6,
        if z ∈ quadraticRadical c then
          64 * binarySign (quadraticPart c z + dot6Linear z a)
        else 0 := by
      apply Finset.sum_congr rfl
      intro z _
      by_cases hz : quadraticPolarLinear c z = 0
      · simp [quadraticRadical, hz, binarySign_zero, mul_comm]
      · simp [quadraticRadical, hz,
          binarySign_sum_state6_linear_nonzero (quadraticPolarLinear c z) hz]
    _ = 64 * ∑ z ∈ quadraticRadical c,
        binarySign (quadraticPart c z + dot6Linear z a) := by
      simp [Finset.mul_sum]

end Spin.UniformMixer
