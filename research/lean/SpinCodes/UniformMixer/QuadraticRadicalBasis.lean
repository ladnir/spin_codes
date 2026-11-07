import SpinCodes.UniformMixer.QuadraticRadical

/-!
# Six basis tests for the quadratic radical

The polar map is linear in its second argument. Checking six basis
vectors therefore determines whether an input lies in its radical.
This is the bridge to a compact finite radical-size certificate.
-/

namespace Spin.UniformMixer

def basis6 (i : Fin 6) : State6 :=
  fun j => if j = i then 1 else 0

theorem basis6_eq_basisFun (i : Fin 6) :
    basis6 i = (Pi.basisFun (ZMod 2) (Fin 6)) i := by
  funext j
  simp [basis6, Pi.basisFun_apply, Pi.single_apply]

theorem linearMap_state6_eq_zero_iff
    (ℓ : State6 →ₗ[ZMod 2] ZMod 2) :
    ℓ = 0 ↔ ∀ i : Fin 6,
      ℓ (basis6 i) = 0 := by
  constructor
  · intro h i
    simp [h]
  · intro h
    apply (Pi.basisFun (ZMod 2) (Fin 6)).ext
    intro i
    simpa [basis6_eq_basisFun] using h i

theorem quadraticRadical_basis_mem (c : Fin 15 → ZMod 2)
    (z : State6) :
    z ∈ quadraticRadical c ↔
      ∀ i : Fin 6,
        quadraticPolarLinear c z
          (basis6 i) = 0 := by
  classical
  simp [quadraticRadical, linearMap_state6_eq_zero_iff]

def basisRadical (c : Fin 15 → ZMod 2) : Finset State6 :=
  Finset.univ.filter (fun z => decide
    (∀ i : Fin 6,
      quadraticPolarLinear c z
        (basis6 i) = 0) = true)

theorem basisRadical_eq_quadraticRadical (c : Fin 15 → ZMod 2) :
    basisRadical c = quadraticRadical c := by
  classical
  ext z
  simp [basisRadical, quadraticRadical_basis_mem]

end Spin.UniformMixer
