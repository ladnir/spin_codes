import SpinCodes.UniformMixer.Map22QuadraticSpace

/-!
# Binary character sums for quadratic-spectrum certificates

The production image is the full degree-at-most-two Boolean code on six
variables. Its spectrum can be obtained by grouping words by quadratic
part and evaluating the Walsh transform over their affine parts. These
generic lemmas provide the additive character and the cancellation used
in that reduction; they assert no numerical spectrum on their own.
-/

namespace Spin.UniformMixer

def binarySign (a : ZMod 2) : ℤ := if a = 0 then 1 else -1

theorem binarySign_zero : binarySign 0 = 1 := by decide

theorem binarySign_one : binarySign 1 = -1 := by decide

theorem binarySign_add (a b : ZMod 2) :
    binarySign (a + b) = binarySign a * binarySign b := by
  fin_cases a <;> fin_cases b <;> decide

theorem binarySign_neg (a : ZMod 2) :
    binarySign (a + 1) = -binarySign a := by
  rw [binarySign_add, binarySign_one]
  ring

theorem binarySign_sum_linear {V : Type*}
    [AddCommGroup V] [Module (ZMod 2) V] [Fintype V]
    (ℓ : V →ₗ[ZMod 2] ZMod 2) (u : V) (hu : ℓ u = 1) :
    (∑ v : V, binarySign (ℓ v)) = 0 := by
  classical
  have hshift :
      (∑ v : V, binarySign (ℓ (v + u))) =
        ∑ v : V, binarySign (ℓ v) := by
    apply Fintype.sum_equiv (Equiv.addRight u)
    intro v
    rfl
  have hneg (v : V) : binarySign (ℓ (v + u)) =
      -binarySign (ℓ v) := by
    rw [map_add, hu, binarySign_neg]
  simp_rw [hneg] at hshift
  simp only [Finset.sum_neg_distrib] at hshift
  omega

theorem binarySign_sum_zero_linear {V : Type*}
    [AddCommGroup V] [Module (ZMod 2) V] [Fintype V]
    (ℓ : V →ₗ[ZMod 2] ZMod 2) (hℓ : ℓ = 0) :
    (∑ v : V, binarySign (ℓ v)) = Fintype.card V := by
  classical
  subst ℓ
  simp [binarySign_zero]

abbrev State6 := Fin 6 → ZMod 2

/-- The full family of six-variable binary linear characters. -/
def dot6Linear (x : State6) : State6 →ₗ[ZMod 2] ZMod 2 where
  toFun a := ∑ i : Fin 6, x i * a i
  map_add' a b := by
    simp [Pi.add_apply, mul_add, Finset.sum_add_distrib]
  map_smul' c a := by
    simp only [Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    simp only [RingHom.id_apply]
    ring

theorem dot6Linear_separates (x : State6) (hx : x ≠ 0) :
    ∃ a : State6, dot6Linear x a = 1 := by
  have hi : ∃ i : Fin 6, x i ≠ 0 := by
    by_contra h
    push Not at h
    apply hx
    funext i
    exact h i
  obtain ⟨i, hxi⟩ := hi
  have hx1 : x i = 1 := by
    rcases zmod2_zero_or_one (x i) with h0 | h1
    · exact False.elim (hxi h0)
    · exact h1
  refine ⟨fun j => if j = i then 1 else 0, ?_⟩
  change (∑ j : Fin 6, x j * (if j = i then 1 else 0)) = 1
  simp [hx1]

theorem dot6_character_sum (x : State6) :
    (∑ a : State6, binarySign (dot6Linear x a)) =
      if x = 0 then 64 else 0 := by
  by_cases hx : x = 0
  · subst x
    simp [dot6Linear, binarySign_zero]
  · obtain ⟨a, ha⟩ := dot6Linear_separates x hx
    simp [hx, binarySign_sum_linear (dot6Linear x) a ha]

theorem dot6Linear_add (x y a : State6) :
    dot6Linear (x + y) a = dot6Linear x a + dot6Linear y a := by
  simp [dot6Linear, Pi.add_apply, add_mul, Finset.sum_add_distrib]

theorem bitScalar_add_self (a : ZMod 2) : a + a = 0 := by
  fin_cases a <;> decide

theorem state6_add_self (x : State6) : x + x = 0 := by
  funext i
  exact bitScalar_add_self (x i)

theorem state6_add_eq_zero_iff (x y : State6) : x + y = 0 ↔ x = y := by
  constructor
  · intro h
    calc
      x = x + 0 := by simp
      _ = x + (x + y) := by rw [h]
      _ = (x + x) + y := by abel
      _ = y := by rw [state6_add_self]; simp
  · rintro rfl
    exact state6_add_self x

/-- Orthogonality of two Walsh characters, including their Boolean phases. -/
theorem dot6_pair_character_sum (f : State6 → ZMod 2)
    (x y : State6) :
    (∑ a : State6,
      binarySign (f x + dot6Linear x a) *
        binarySign (f y + dot6Linear y a)) =
      if x = y then 64 else 0 := by
  classical
  calc
    _ = binarySign (f x + f y) *
        (∑ a : State6, binarySign (dot6Linear (x + y) a)) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro a _
      rw [← binarySign_add, ← binarySign_add, dot6Linear_add]
      congr 1
      abel
    _ = binarySign (f x + f y) *
        (if x + y = 0 then 64 else 0) := by
      rw [dot6_character_sum]
    _ = if x = y then 64 else 0 := by
      simp only [state6_add_eq_zero_iff]
      by_cases hxy : x = y
      · subst y
        simp [bitScalar_add_self, binarySign_zero]
      · simp [hxy]

/-- Walsh transform over the six-variable affine characters. -/
def walsh6 (f : State6 → ZMod 2) (a : State6) : ℤ :=
  ∑ x : State6, binarySign (f x + dot6Linear x a)

theorem walsh6_parseval (f : State6 → ZMod 2) :
    (∑ a : State6, walsh6 f a * walsh6 f a) = 4096 := by
  classical
  unfold walsh6
  calc
    _ = ∑ a : State6, ∑ x : State6, ∑ y : State6,
        binarySign (f x + dot6Linear x a) *
          binarySign (f y + dot6Linear y a) := by
      apply Finset.sum_congr rfl
      intro a _
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro x _
      rw [Finset.mul_sum]
    _ = ∑ x : State6, ∑ y : State6, ∑ a : State6,
        binarySign (f x + dot6Linear x a) *
          binarySign (f y + dot6Linear y a) := by
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro x _
      rw [Finset.sum_comm]
    _ = ∑ x : State6, ∑ y : State6,
        if x = y then (64 : ℤ) else 0 := by
      simp_rw [dot6_pair_character_sum]
    _ = 4096 := by
      simp

def boolWeight6 (f : State6 → ZMod 2) : ℕ :=
  ∑ x : State6, if f x = 1 then 1 else 0

theorem binarySign_as_weight (a : ZMod 2) :
    binarySign a = 1 - 2 * (if a = 1 then (1 : ℤ) else 0) := by
  fin_cases a <;> decide

theorem binarySign_sum_weight (f : State6 → ZMod 2) :
    (∑ x : State6, binarySign (f x)) =
      64 - 2 * (boolWeight6 f : ℤ) := by
  classical
  simp_rw [binarySign_as_weight]
  simp [boolWeight6, Finset.sum_sub_distrib, Finset.sum_ite, mul_comm]

theorem walsh6_constant_shift (f : State6 → ZMod 2)
    (c : ZMod 2) (a : State6) :
    walsh6 (fun x => f x + c) a = binarySign c * walsh6 f a := by
  unfold walsh6
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro x _
  rw [← binarySign_add]
  congr 1
  change f x + c + dot6Linear x a = c + (f x + dot6Linear x a)
  abel

theorem walsh6_weight (f : State6 → ZMod 2) (a : State6) :
    walsh6 f a = 64 - 2 *
      (boolWeight6 (fun x => f x + dot6Linear x a) : ℤ) := by
  exact binarySign_sum_weight _

end Spin.UniformMixer
