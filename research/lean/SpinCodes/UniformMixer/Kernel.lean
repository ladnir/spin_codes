import Mathlib.Data.Fintype.Card
import Mathlib.Data.Fintype.EquivFin
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.GroupTheory.GroupAction.Basic
import Mathlib.LinearAlgebra.FiniteDimensional.Basic
import Mathlib.LinearAlgebra.LinearIndependent.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Field.ZMod
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp

/-!
# Uniform-mixer transition kernel

This file formalizes the finite-state law used by the new Scalable SPIN inner
and derives it from a uniformly sampled invertible linear map. Nothing here
uses the older IMT transvection law.
-/

noncomputable section

namespace Spin.UniformMixer

variable {S : Type*} [Fintype S] [DecidableEq S] [AddCommGroup S]

/-- Number of nonzero states. -/
def liveCount (S : Type*) [Fintype S] : ℕ := Fintype.card S - 1

omit [DecidableEq S] [AddCommGroup S] in
theorem liveCount_pos [Nontrivial S] : 0 < liveCount S := by
  simpa [liveCount] using (Fintype.one_lt_card : 1 < Fintype.card S)

/-- The state transition after a uniform refresh of a nonzero state, followed
by an additive feedback value `c`. A zero incoming state receives no refresh. -/
def update (c q r : S) : ℝ :=
  if q = 0 then (if r = c then 1 else 0)
  else if r = c then 0 else 1 / (liveCount S : ℝ)

@[simp] theorem update_zero (c r : S) :
    update c 0 r = if r = c then 1 else 0 := by
  simp [update]

@[simp] theorem update_live (c q r : S) (hq : q ≠ 0) :
    update c q r = if r = c then 0 else 1 / (liveCount S : ℝ) := by
  simp [update, hq]

@[simp] theorem update_zero_zero : update (0 : S) 0 0 = 1 := by
  simp

@[simp] theorem update_live_to_zero (c q : S) (hq : q ≠ 0) :
    update c q 0 = if c = 0 then 0 else 1 / (liveCount S : ℝ) := by
  simp [update_live c q 0 hq, eq_comm]

theorem update_nonneg [Nontrivial S] (c q r : S) : 0 ≤ update c q r := by
  by_cases hq : q = 0
  · subst q
    simp only [update_zero]
    split_ifs <;> simp
  · rw [update_live c q r hq]
    split_ifs
    · exact le_refl _
    · exact one_div_nonneg.mpr (Nat.cast_nonneg _)

theorem live_to_zero_probability [Nontrivial S] (c q : S)
    (hq : q ≠ 0) (hc : c ≠ 0) :
    update c q 0 = 1 / (liveCount S : ℝ) := by
  simp [update_live_to_zero c q hq, hc]

theorem update_sum [Nontrivial S] (c q : S) :
    ∑ r : S, update c q r = 1 := by
  by_cases hq : q = 0
  · subst q
    simp [update_zero]
  · simp_rw [update_live c q _ hq]
    let k : ℝ := 1 / (liveCount S : ℝ)
    have hsplit : (∑ r : S, if r = c then (0 : ℝ) else k) =
        (∑ r : S, k) - (∑ r : S, if r = c then k else 0) := by
      rw [← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro r _
      split_ifs <;> simp
    change (∑ r : S, if r = c then (0 : ℝ) else k) = 1
    rw [hsplit]
    simp only [Finset.sum_ite_eq', Finset.sum_const, Finset.card_univ,
      Finset.mem_univ, ite_true]
    rw [nsmul_eq_mul]
    have hle : 1 ≤ Fintype.card S := le_of_lt Fintype.one_lt_card
    have hd : (liveCount S : ℝ) = (Fintype.card S : ℝ) - 1 := by
      simp [liveCount, Nat.cast_sub hle]
    calc
      (Fintype.card S : ℝ) * k - k = ((Fintype.card S : ℝ) - 1) * k := by ring
      _ = (liveCount S : ℝ) * (1 / (liveCount S : ℝ)) := by rw [← hd]
      _ = 1 := by
        have hn : (liveCount S : ℝ) ≠ 0 := by
          exact_mod_cast (Nat.ne_of_gt (liveCount_pos (S := S)))
        simp [hn]

/-- Probability mass assigned to the nonzero outgoing states. -/
def liveMass (c q : S) : ℝ :=
  ∑ r ∈ (Finset.univ : Finset S).erase 0, update c q r

theorem liveMass_eq [Nontrivial S] (c q : S) :
    liveMass c q = 1 - update c q 0 := by
  have h := Finset.add_sum_erase (Finset.univ : Finset S)
    (fun r => update c q r) (Finset.mem_univ (0 : S))
  change update c q 0 + liveMass c q = ∑ r : S, update c q r at h
  rw [update_sum] at h
  linarith

theorem zero_to_live [Nontrivial S] (c : S) (hc : c ≠ 0) :
    liveMass c 0 = 1 := by
  rw [liveMass_eq]
  simp [update_zero, hc, eq_comm]

theorem live_to_live [Nontrivial S] (c q : S) (hq : q ≠ 0)
    (hc : c ≠ 0) :
    liveMass c q = 1 - 1 / (liveCount S : ℝ) := by
  rw [liveMass_eq, live_to_zero_probability c q hq hc]

/-- The one-bit impulse quotient. The first row enters at zero; the
second enters at any nonzero state. Output columns are zero and live. -/
theorem impulse_quotient [Nontrivial S] (c q : S)
    (hc : c ≠ 0) (hq : q ≠ 0) :
    ((update c 0 0, liveMass c 0),
      (update c q 0, liveMass c q)) =
      (((0 : ℝ), 1),
        (1 / (liveCount S : ℝ), 1 - 1 / (liveCount S : ℝ))) := by
  simp [update_zero, hc, eq_comm, zero_to_live,
    live_to_live, hq]

/-- An empty input emits the expansion of the entering state and then
refreshes that state. `w` is the Hamming weight of the expansion. -/
def emptyStep (z : ℝ) (w : S → ℕ) (q r : S) : ℝ :=
  z ^ w q * update 0 q r

theorem emptyStep_zero (z : ℝ) (w : S → ℕ) (hw : w 0 = 0) (r : S) :
    emptyStep z w 0 r = if r = 0 then 1 else 0 := by
  simp [emptyStep, hw]

theorem emptyStep_live_zero (z : ℝ) (w : S → ℕ) (q : S) (hq : q ≠ 0) :
    emptyStep z w q 0 = 0 := by
  simp [emptyStep, update_live, hq]

theorem emptyStep_live (z : ℝ) (w : S → ℕ) (q r : S)
    (hq : q ≠ 0) (hr : r ≠ 0) :
    emptyStep z w q r = z ^ w q / (liveCount S : ℝ) := by
  simp [emptyStep, update_live, hq, hr, div_eq_mul_inv]

/-- The endpoint kernel for `g` empty steps, without conditioning on a
zero or uniform incoming state. -/
def emptyGap (z : ℝ) (w : S → ℕ) : ℕ → S → S → ℝ
  | 0, q, r => if q = r then 1 else 0
  | g + 1, q, r => ∑ v : S, emptyGap z w g q v * emptyStep z w v r

theorem emptyGap_one (z : ℝ) (w : S → ℕ) (q r : S) :
    emptyGap z w 1 q r = emptyStep z w q r := by
  simp [emptyGap]

theorem emptyGap_zero (z : ℝ) (w : S → ℕ) (hw : w 0 = 0)
    (g : ℕ) (r : S) :
    emptyGap z w g 0 r = if r = 0 then 1 else 0 := by
  induction g generalizing r with
  | zero => simp [emptyGap, eq_comm]
  | succ g ih =>
      rw [emptyGap]
      simp_rw [ih]
      simp [emptyStep_zero z w hw]

theorem emptyGap_live_zero (z : ℝ) (w : S → ℕ) (g : ℕ)
    (q : S) (hq : q ≠ 0) :
    emptyGap z w (g + 1) q 0 = 0 := by
  induction g with
  | zero => rw [emptyGap_one, emptyStep_live_zero z w q hq]
  | succ g ih =>
      rw [emptyGap]
      apply Finset.sum_eq_zero
      intro v _
      by_cases hv : v = 0
      · subst v
        simp [ih]
      · simp [emptyStep_live_zero z w v hv]

/-- The mean weight generating function over nonzero states. -/
def liveMean (z : ℝ) (w : S → ℕ) : ℝ :=
  (∑ v ∈ (Finset.univ : Finset S).erase 0, z ^ w v) /
    (liveCount S : ℝ)

theorem emptyGap_live (z : ℝ) (w : S → ℕ) (g : ℕ)
    (q r : S) (hq : q ≠ 0) (hr : r ≠ 0) :
    emptyGap z w (g + 1) q r =
      z ^ w q / (liveCount S : ℝ) * (liveMean z w) ^ g := by
  induction g generalizing r with
  | zero => simp [emptyGap_one, emptyStep_live, hq, hr]
  | succ g ih =>
      rw [emptyGap]
      have hz := emptyGap_live_zero z w g q hq
      rw [← Finset.add_sum_erase (Finset.univ : Finset S)
        (fun v => emptyGap z w (g + 1) q v * emptyStep z w v r)
        (Finset.mem_univ (0 : S))]
      simp only [hz, zero_mul, zero_add]
      calc
        (∑ v ∈ (Finset.univ : Finset S).erase 0,
          emptyGap z w (g + 1) q v * emptyStep z w v r) =
            ∑ v ∈ (Finset.univ : Finset S).erase 0,
              (z ^ w q / (liveCount S : ℝ) * liveMean z w ^ g) *
                (z ^ w v / (liveCount S : ℝ)) := by
                  apply Finset.sum_congr rfl
                  intro v hv
                  have hv0 : v ≠ 0 := (Finset.mem_erase.mp hv).1
                  rw [ih v hv0, emptyStep_live z w v r hv0 hr]
        _ = z ^ w q / (liveCount S : ℝ) * liveMean z w ^ (g + 1) := by
          rw [pow_succ]
          simp only [liveMean]
          rw [← Finset.mul_sum]
          simp only [div_eq_mul_inv]
          rw [← Finset.sum_mul]
          ring

end Spin.UniformMixer

namespace Spin.UniformMixer

variable {F V : Type*} [Field F] [AddCommGroup V] [Module F V]
  [FiniteDimensional F V]

theorem linearEquiv_transitive_live (a b : V) (ha : a ≠ 0) (hb : b ≠ 0) :
    ∃ g : V ≃ₗ[F] V, g • a = b := by
  let va : Unit → V := fun _ => a
  let vb : Unit → V := fun _ => b
  have hia : LinearIndependent F va := linearIndependent_unique_iff.mpr ha
  have hib : LinearIndependent F vb := linearIndependent_unique_iff.mpr hb
  let f := hia.linearCombinationEquiv.symm ≪≫ₗ hib.linearCombinationEquiv
  have hfin : FiniteDimensional F (Submodule.span F (Set.range va)) :=
    inferInstance
  obtain ⟨g, hg⟩ := Submodule.exists_linearEquiv_restrict_eq f
  refine ⟨g, ?_⟩
  let e : Unit →₀ F := Finsupp.single () 1
  have hae : ((hia.linearCombinationEquiv e : Submodule.span F (Set.range va)) : V) = a := by
    simp [e, va]
  have hbe : ((hib.linearCombinationEquiv e : Submodule.span F (Set.range vb)) : V) = b := by
    simp [e, vb]
  calc
    g • a = g ((hia.linearCombinationEquiv e : Submodule.span F (Set.range va)) : V) := by
      rw [hae]
      rfl
    _ = (f (hia.linearCombinationEquiv e) : V) := (hg _).symm
    _ = b := by simpa [f] using hbe

end Spin.UniformMixer

namespace Spin.UniformMixer

variable {G S : Type*} [Group G] [Fintype G]
  [AddCommGroup S] [Fintype S] [MulAction G S]

/-- The group elements sending `q` to `r`. -/
def actionFiber (q r : S) := {g : G // g • q = r}

instance (q r : S) : Fintype (actionFiber (G := G) q r) := by
  classical
  unfold actionFiber
  infer_instance

/-- Left multiplication identifies any two fibers in the same orbit. -/
def actionFiberEquiv (q r s : S) (h : G) (hh : h • r = s) :
    actionFiber (G := G) q r ≃ actionFiber (G := G) q s where
  toFun g := ⟨h * g.1, by rw [mul_smul, g.2, hh]⟩
  invFun g := ⟨h⁻¹ * g.1, by
    rw [mul_smul, g.2]
    calc
      h⁻¹ • s = h⁻¹ • (h • r) := by rw [hh]
      _ = r := inv_smul_smul h r⟩
  left_inv g := by apply Subtype.ext; simp
  right_inv g := by apply Subtype.ext; simp

omit [AddCommGroup S] [Fintype S] in
theorem actionFiber_card_eq (q r s : S) (h : G) (hh : h • r = s) :
    Fintype.card (actionFiber (G := G) q r) =
      Fintype.card (actionFiber (G := G) q s) :=
  Fintype.card_congr (actionFiberEquiv q r s h hh)

/-- A group element is uniquely determined by its image of `q` together
with its position inside that image's fiber. -/
def actionSigma (q : S) : G ≃ Σ r : S, actionFiber (G := G) q r where
  toFun g := ⟨g • q, ⟨g, rfl⟩⟩
  invFun x := x.2.1
  left_inv g := rfl
  right_inv x := by
    rcases x with ⟨r, g, hg⟩
    dsimp
    cases hg
    rfl

omit [AddCommGroup S] in
theorem card_group_eq_sum_fibers (q : S) :
    Fintype.card G = ∑ r : S, Fintype.card (actionFiber (G := G) q r) := by
  rw [Fintype.card_congr (actionSigma q), Fintype.card_sigma]

omit [Fintype S] in
theorem actionFiber_zero_card (q : S) (hq : q ≠ 0)
    (hfix0 : ∀ g : G, g • (0 : S) = 0) :
    Fintype.card (actionFiber (G := G) q 0) = 0 := by
  have hempty : IsEmpty (actionFiber (G := G) q 0) :=
    ⟨fun g => hq (by
      calc
        q = g.1⁻¹ • (g.1 • q) := (inv_smul_smul g.1 q).symm
        _ = g.1⁻¹ • 0 := by rw [g.2]
        _ = 0 := hfix0 g.1⁻¹)⟩
  simp

theorem card_group_eq_liveCount_mul_fiber (q r : S)
    (hq : q ≠ 0) (hr : r ≠ 0)
    (hfix0 : ∀ g : G, g • (0 : S) = 0)
    (htrans : ∀ a b : S, a ≠ 0 → b ≠ 0 → ∃ g : G, g • a = b) :
    Fintype.card G = liveCount S * Fintype.card (actionFiber (G := G) q r) := by
  classical
  rw [card_group_eq_sum_fibers (G := G) q]
  have hsplit := Finset.add_sum_erase (Finset.univ : Finset S)
    (fun v => Fintype.card (actionFiber (G := G) q v))
    (Finset.mem_univ (0 : S))
  rw [actionFiber_zero_card q hq hfix0, zero_add] at hsplit
  rw [← hsplit]
  have hsum : (∑ v ∈ (Finset.univ : Finset S).erase 0,
      Fintype.card (actionFiber (G := G) q v)) =
      ∑ _v ∈ (Finset.univ : Finset S).erase 0,
        Fintype.card (actionFiber (G := G) q r) := by
    apply Finset.sum_congr rfl
    intro v hv
    obtain ⟨g, hg⟩ := htrans v r (Finset.mem_erase.mp hv).1 hr
    exact actionFiber_card_eq q v r g hg
  rw [hsum, Finset.sum_const, nsmul_eq_mul]
  simp [liveCount]

/-- The state update induced by choosing a uniformly random element of `G`.
The event is `g • q + c = r`, equivalently `g • q = r - c`. -/
def groupUpdate (c q r : S) : ℝ := by
  classical
  exact (Fintype.card (actionFiber (G := G) q (r - c)) : ℝ) / Fintype.card G

/-- The literal update event in the encoder. -/
def updateEvent (c q r : S) := {g : G // g • q + c = r}

instance (c q r : S) : Fintype (updateEvent (G := G) c q r) := by
  classical
  unfold updateEvent
  infer_instance

omit [Fintype S] in
theorem groupUpdate_eq_event_count (c q r : S) :
    groupUpdate (G := G) c q r =
      (Fintype.card (updateEvent (G := G) c q r) : ℝ) / Fintype.card G := by
  have hcard : Fintype.card (actionFiber (G := G) q (r - c)) =
      Fintype.card (updateEvent (G := G) c q r) := by
    apply Fintype.card_congr
    exact
      { toFun := fun g => ⟨g.1, by rw [g.2, sub_add_cancel]⟩
        invFun := fun g => ⟨g.1, by
          have hg := g.2
          exact eq_sub_iff_add_eq.mpr hg⟩
        left_inv := by intro g; rfl
        right_inv := by intro g; rfl }
  simp only [groupUpdate, hcard]

theorem groupUpdate_live [DecidableEq S] [Nontrivial S] (c q r : S) (hq : q ≠ 0)
    (hfix0 : ∀ g : G, g • (0 : S) = 0)
    (htrans : ∀ a b : S, a ≠ 0 → b ≠ 0 → ∃ g : G, g • a = b) :
    groupUpdate (G := G) c q r = update c q r := by
  classical
  by_cases hrc : r = c
  · subst r
    rw [update_live c q c hq]
    simp only [ite_true, sub_self, groupUpdate]
    rw [actionFiber_zero_card q hq hfix0]
    simp
  · rw [update_live c q r hq]
    simp only [hrc, ite_false]
    have hv : r - c ≠ 0 := sub_ne_zero.mpr hrc
    obtain ⟨g, hg⟩ := htrans q (r - c) hq hv
    have hFpos : 0 < Fintype.card (actionFiber (G := G) q (r - c)) :=
      Fintype.card_pos_iff.mpr ⟨⟨g, hg⟩⟩
    have hF : (Fintype.card (actionFiber (G := G) q (r - c)) : ℝ) ≠ 0 := by
      exact_mod_cast (Nat.ne_of_gt hFpos)
    have hD : (liveCount S : ℝ) ≠ 0 := by
      exact_mod_cast (Nat.ne_of_gt (liveCount_pos (S := S)))
    have hcard := card_group_eq_liveCount_mul_fiber q (r - c) hq hv hfix0 htrans
    have hcardReal : (Fintype.card G : ℝ) =
        (liveCount S : ℝ) *
          (Fintype.card (actionFiber (G := G) q (r - c)) : ℝ) := by
      exact_mod_cast hcard
    rw [groupUpdate, hcardReal]
    field_simp

omit [Fintype S] in
theorem actionFiber_zero_zero_card
    (hfix0 : ∀ g : G, g • (0 : S) = 0) :
    Fintype.card (actionFiber (G := G) (0 : S) (0 : S)) = Fintype.card G := by
  apply Fintype.card_congr
  exact
    { toFun := fun g => g.1
      invFun := fun g => ⟨g, hfix0 g⟩
      left_inv := by intro g; apply Subtype.ext; rfl
      right_inv := by intro g; rfl }

omit [Fintype S] in
theorem actionFiber_zero_ne_card (r : S) (hr : r ≠ 0)
    (hfix0 : ∀ g : G, g • (0 : S) = 0) :
    Fintype.card (actionFiber (G := G) (0 : S) r) = 0 := by
  have hempty : IsEmpty (actionFiber (G := G) 0 r) :=
    ⟨fun g => hr (by rw [← g.2, hfix0])⟩
  simp

theorem groupUpdate_zero [DecidableEq S] (c r : S)
    (hfix0 : ∀ g : G, g • (0 : S) = 0) :
    groupUpdate (G := G) c 0 r = update c 0 r := by
  classical
  rw [update_zero]
  by_cases hrc : r = c
  · subst r
    simp only [ite_true, sub_self, groupUpdate]
    rw [actionFiber_zero_zero_card hfix0]
    have hG : (Fintype.card G : ℝ) ≠ 0 := by
      exact_mod_cast (Nat.ne_of_gt Fintype.card_pos)
    simp [hG]
  · simp only [hrc, ite_false, groupUpdate]
    rw [actionFiber_zero_ne_card (r - c) (sub_ne_zero.mpr hrc) hfix0]
    simp

theorem groupUpdate_eq_update [DecidableEq S] [Nontrivial S] (c q r : S)
    (hfix0 : ∀ g : G, g • (0 : S) = 0)
    (htrans : ∀ a b : S, a ≠ 0 → b ≠ 0 → ∃ g : G, g • a = b) :
    groupUpdate (G := G) c q r = update c q r := by
  by_cases hq : q = 0
  · subst q
    exact groupUpdate_zero c r hfix0
  · exact groupUpdate_live c q r hq hfix0 htrans

end Spin.UniformMixer

namespace Spin.UniformMixer

variable {F V : Type*} [Field F] [AddCommGroup V] [Module F V]
  [FiniteDimensional F V] [Fintype V] [DecidableEq V] [Nontrivial V]
  [Fintype (V ≃ₗ[F] V)]

/-- The abstract refresh kernel is the actual law of a uniform invertible
linear map on any finite-dimensional vector space. -/
theorem linearEquiv_groupUpdate_eq_update (c q r : V) :
    groupUpdate (G := V ≃ₗ[F] V) c q r = update c q r := by
  apply groupUpdate_eq_update c q r
  · intro g
    exact map_zero g
  · intro a b ha hb
    exact linearEquiv_transitive_live a b ha hb

end Spin.UniformMixer

namespace Spin.UniformMixer

abbrev State22 := Fin 22 → ZMod 2

theorem card_binary_vector (m : ℕ) :
    Fintype.card (Fin m → ZMod 2) = 2 ^ m := by
  classical
  rw [Fintype.card_fun, ZMod.card]
  simp

theorem liveCount_binary_vector (m : ℕ) :
    liveCount (Fin m → ZMod 2) = 2 ^ m - 1 := by
  unfold liveCount
  rw [card_binary_vector]

theorem state22_liveCount : liveCount State22 = 2 ^ 22 - 1 :=
  liveCount_binary_vector 22

/-- Finiteness follows by embedding every automorphism in the finite set
of all functions on the state space. -/
theorem state22GL_finite : Finite (State22 ≃ₗ[ZMod 2] State22) :=
  Finite.of_injective
    (fun g : State22 ≃ₗ[ZMod 2] State22 => (g : State22 → State22))
    (fun _ _ h => LinearEquiv.ext (congrFun h))

theorem state22_uniform_gl_kernel
    [Fintype (State22 ≃ₗ[ZMod 2] State22)] (c q r : State22) :
    groupUpdate (G := State22 ≃ₗ[ZMod 2] State22) c q r = update c q r := by
  apply groupUpdate_eq_update c q r
  · intro g
    exact map_zero g
  · intro a b ha hb
    exact linearEquiv_transitive_live (F := ZMod 2) a b ha hb

end Spin.UniformMixer
