import SpinCodes.UniformMixer.QuadraticPolar

/-!
# Additive structure of the quadratic radical

The radical is an additive subgroup of the six-variable state space.
The quadratic phase restricts to an additive character there, which
allows the remaining character sum to be evaluated by translation.
-/

namespace Spin.UniformMixer

theorem quadraticPolarLinear_add (c : Fin 15 → ZMod 2)
    (z w : State6) :
    quadraticPolarLinear c (z + w) =
      quadraticPolarLinear c z + quadraticPolarLinear c w := by
  apply LinearMap.ext
  intro x
  simp [quadraticPolarLinear, Pi.add_apply,
    mul_add, add_mul, Finset.sum_add_distrib]
  ring

theorem state6_neg_eq_self (x : State6) : -x = x := by
  funext i
  rcases zmod2_zero_or_one (x i) with h0 | h1
  · simp [h0]
  · simp [h1]

def quadraticRadicalGroup (c : Fin 15 → ZMod 2) : AddSubgroup State6 where
  carrier := {z | quadraticPolarLinear c z = 0}
  zero_mem' := by
    change quadraticPolarLinear c 0 = 0
    apply LinearMap.ext
    intro x
    simp [quadraticPolarLinear]
  add_mem' := by
    intro z w hz hw
    change quadraticPolarLinear c (z + w) = 0
    rw [quadraticPolarLinear_add, hz, hw]
    simp
  neg_mem' := by
    intro z hz
    simpa only [state6_neg_eq_self] using hz

noncomputable instance (c : Fin 15 → ZMod 2) :
    Fintype (quadraticRadicalGroup c) := by
  classical
  exact Subtype.fintype (fun z : State6 => quadraticPolarLinear c z = 0)

theorem quadraticRadicalGroup_mem (c : Fin 15 → ZMod 2) (z : State6) :
    z ∈ quadraticRadicalGroup c ↔ quadraticPolarLinear c z = 0 := Iff.rfl

theorem zmod2_three_add_eq_zero (a b d : ZMod 2) :
    a + b + d = 0 ↔ a = b + d := by
  constructor
  · intro h
    have hself := bitScalar_add_self (b + d)
    calc
      a = a + 0 := by simp
      _ = a + ((b + d) + (b + d)) := by rw [hself]
      _ = (a + b + d) + (b + d) := by abel
      _ = b + d := by rw [h]; simp
  · intro h
    rw [h]
    calc
      (b + d) + b + d = (b + d) + (b + d) := by abel
      _ = 0 := bitScalar_add_self (b + d)

theorem quadraticPart_add_of_radical (c : Fin 15 → ZMod 2)
    {z : State6} (hz : z ∈ quadraticRadicalGroup c) (w : State6) :
    quadraticPart c (z + w) = quadraticPart c z + quadraticPart c w := by
  have hp := quadraticPart_polar c z w
  have hz0 : quadraticPolarLinear c z = 0 := hz
  rw [hz0] at hp
  have hzero : quadraticPart c (w + z) +
      quadraticPart c w + quadraticPart c z = 0 := by
    simpa [polar6] using hp
  have h := (zmod2_three_add_eq_zero _ _ _).mp hzero
  simpa only [add_comm z w, add_comm (quadraticPart c w)
    (quadraticPart c z)] using h

def quadraticRadicalPhase (c : Fin 15 → ZMod 2) (a : State6) :
    quadraticRadicalGroup c →+ ZMod 2 where
  toFun z := quadraticPart c z.val + dot6Linear z.val a
  map_zero' := by
    change quadraticPart c 0 + dot6Linear 0 a = 0
    simp [quadraticPart_zero, dot6Linear]
  map_add' z w := by
    change quadraticPart c (z.val + w.val) +
      dot6Linear (z.val + w.val) a =
        (quadraticPart c z.val + dot6Linear z.val a) +
          (quadraticPart c w.val + dot6Linear w.val a)
    rw [quadraticPart_add_of_radical c z.property w.val,
      dot6Linear_add]
    abel

theorem binarySign_sum_addHom {V : Type*} [AddCommGroup V] [Fintype V]
    (χ : V →+ ZMod 2) (u : V) (hu : χ u = 1) :
    (∑ v : V, binarySign (χ v)) = 0 := by
  classical
  have hshift :
      (∑ v : V, binarySign (χ (v + u))) =
        ∑ v : V, binarySign (χ v) := by
    apply Fintype.sum_equiv (Equiv.addRight u)
    intro v
    rfl
  have hneg (v : V) : binarySign (χ (v + u)) =
      -binarySign (χ v) := by
    rw [map_add, hu, binarySign_neg]
  simp_rw [hneg] at hshift
  simp only [Finset.sum_neg_distrib] at hshift
  omega

theorem binarySign_sum_addHom_nonzero {V : Type*}
    [AddCommGroup V] [Fintype V]
    (χ : V →+ ZMod 2) (hχ : χ ≠ 0) :
    (∑ v : V, binarySign (χ v)) = 0 := by
  have hex : ∃ v : V, χ v ≠ 0 := by
    by_contra h
    push Not at h
    apply hχ
    apply AddMonoidHom.ext
    intro v
    exact h v
  obtain ⟨v, hv⟩ := hex
  have hv1 : χ v = 1 :=
    (zmod2_zero_or_one (χ v)).resolve_left hv
  exact binarySign_sum_addHom χ v hv1

theorem quadraticRadical_sum_eq_subtype
    (c : Fin 15 → ZMod 2) (a : State6) :
    (∑ z ∈ quadraticRadical c,
      binarySign (quadraticPart c z + dot6Linear z a)) =
    ∑ z : quadraticRadicalGroup c,
      binarySign (quadraticRadicalPhase c a z) := by
  classical
  change (∑ z ∈ Finset.univ.filter
      (fun z : State6 => quadraticPolarLinear c z = 0),
        binarySign (quadraticPart c z + dot6Linear z a)) = _
  rw [← Finset.sum_subtype_eq_sum_filter]
  rw [Finset.subtype_univ]
  change (∑ z : {z : State6 // quadraticPolarLinear c z = 0},
      binarySign (quadraticPart c z.val + dot6Linear z.val a)) =
    ∑ z : quadraticRadicalGroup c,
      binarySign (quadraticRadicalPhase c a z)
  apply Fintype.sum_equiv (Equiv.refl _)
  intro z
  rfl

theorem quadraticRadical_card (c : Fin 15 → ZMod 2) :
    Fintype.card (quadraticRadicalGroup c) =
      (quadraticRadical c).card := by
  classical
  simpa [quadraticRadical, quadraticRadicalGroup] using
    (Fintype.card_subtype
      (fun z : State6 => quadraticPolarLinear c z = 0))

theorem quadraticPart_walsh_square_dichotomy
    (c : Fin 15 → ZMod 2) (a : State6) :
    walsh6 (quadraticPart c) a * walsh6 (quadraticPart c) a = 0 ∨
    walsh6 (quadraticPart c) a * walsh6 (quadraticPart c) a =
      64 * (quadraticRadical c).card := by
  classical
  rw [quadraticPart_walsh_square, quadraticRadical_sum_eq_subtype]
  by_cases hphase : quadraticRadicalPhase c a = 0
  · right
    have hsum : (∑ z : quadraticRadicalGroup c,
        binarySign (quadraticRadicalPhase c a z)) =
          (Fintype.card (quadraticRadicalGroup c) : ℤ) := by
      rw [hphase]
      simp [binarySign_zero]
    rw [hsum, quadraticRadical_card]
  · left
    rw [binarySign_sum_addHom_nonzero
      (quadraticRadicalPhase c a) hphase]
    simp

noncomputable def activeWalsh (c : Fin 15 → ZMod 2) : Finset State6 := by
  classical
  exact Finset.univ.filter (fun a => walsh6 (quadraticPart c) a ≠ 0)

theorem quadraticPart_walsh_square_nonzero
    (c : Fin 15 → ZMod 2) {a : State6}
    (ha : a ∈ activeWalsh c) :
    walsh6 (quadraticPart c) a * walsh6 (quadraticPart c) a =
      64 * (quadraticRadical c).card := by
  classical
  have hne : walsh6 (quadraticPart c) a ≠ 0 := by
    simpa [activeWalsh] using ha
  rcases quadraticPart_walsh_square_dichotomy c a with h0 | hR
  · rcases mul_eq_zero.mp h0 with h | h <;> exact False.elim (hne h)
  · exact hR

theorem sum_square_by_support {α : Type*} [Fintype α]
    (w : α → ℤ) (C : ℤ)
    (h : ∀ a, w a ≠ 0 → w a * w a = C) :
    (∑ a : α, w a * w a) =
      C * ((Finset.univ.filter fun a => w a ≠ 0).card : ℤ) := by
  classical
  calc
    _ = ∑ a ∈ Finset.univ.filter (fun a => w a ≠ 0), C := by
      rw [Finset.sum_filter]
      apply Finset.sum_congr rfl
      intro a _
      by_cases ha : w a = 0
      · simp [ha]
      · simp [ha, h a ha]
    _ = C * ((Finset.univ.filter fun a => w a ≠ 0).card : ℤ) := by
      simp [mul_comm]

theorem quadraticPart_activeWalsh_count (c : Fin 15 → ZMod 2) :
    (4096 : ℤ) = 64 * ((quadraticRadical c).card : ℤ) *
      ((activeWalsh c).card : ℤ) := by
  classical
  have h := sum_square_by_support
    (fun a : State6 => walsh6 (quadraticPart c) a)
    ((64 : ℤ) * (quadraticRadical c).card)
    (by
      intro a ha
      exact quadraticPart_walsh_square_nonzero c
        (by simpa [activeWalsh] using ha))
  rw [walsh6_parseval] at h
  have hfilter :
      (Finset.univ.filter fun a : State6 =>
        walsh6 (quadraticPart c) a ≠ 0) = activeWalsh c := rfl
  rw [hfilter] at h
  exact h

end Spin.UniformMixer
