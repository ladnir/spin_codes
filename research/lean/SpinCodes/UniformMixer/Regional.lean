import SpinCodes.UniformMixer.Kernel

/-!
# Full-state regional kernels

This file develops the composition algebra for the uniform-mixer kernels.
-/

noncomputable section

namespace Spin.UniformMixer

variable {S : Type*} [Fintype S] [DecidableEq S] [AddCommGroup S]

/-- Composition of full-state kernels. -/
def compose (K L : S → S → ℝ) (q r : S) : ℝ :=
  ∑ v : S, K q v * L v r

/-- The identity full-state kernel. -/
def identity (q r : S) : ℝ := if q = r then 1 else 0

omit [AddCommGroup S] in
theorem compose_identity_left (K : S → S → ℝ) (q r : S) :
    compose identity K q r = K q r := by
  simp [compose, identity]

omit [AddCommGroup S] in
theorem compose_identity_right (K : S → S → ℝ) (q r : S) :
    compose K identity q r = K q r := by
  simp [compose, identity]

omit [DecidableEq S] [AddCommGroup S] in
theorem compose_assoc (K L M : S → S → ℝ) (q r : S) :
    compose (compose K L) M q r = compose K (compose L M) q r := by
  simp only [compose, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro v _
  apply Finset.sum_congr rfl
  intro u _
  ring

omit [DecidableEq S] [AddCommGroup S] in
theorem compose_nonneg (K L : S → S → ℝ)
    (hK : ∀ q r, 0 ≤ K q r) (hL : ∀ q r, 0 ≤ L q r) (q r : S) :
    0 ≤ compose K L q r := by
  unfold compose
  apply Finset.sum_nonneg
  intro v _
  exact mul_nonneg (hK q v) (hL v r)

omit [DecidableEq S] [AddCommGroup S] in
theorem compose_scaled_bound (K L K' L' : S → S → ℝ) (a b : ℝ)
    (ha : 0 ≤ a) (hL : ∀ q r, 0 ≤ L q r)
    (hK' : ∀ q r, 0 ≤ K' q r)
    (hK : ∀ q r, K q r ≤ a * K' q r)
    (hL' : ∀ q r, L q r ≤ b * L' q r) (q r : S) :
    compose K L q r ≤ a * b * compose K' L' q r := by
  calc
    compose K L q r ≤
        ∑ v : S, (a * K' q v) * (b * L' v r) := by
          unfold compose
          apply Finset.sum_le_sum
          intro v _
          exact mul_le_mul (hK q v) (hL' v r)
            (hL v r) (mul_nonneg ha (hK' q v))
    _ = a * b * compose K' L' q r := by
          rw [compose, Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro v _
          ring

theorem emptyGap_succ (z : ℝ) (w : S → ℕ) (g : ℕ) :
    emptyGap z w (g + 1) = compose (emptyGap z w g) (emptyStep z w) := by
  funext q r
  rfl

theorem emptyGap_add (z : ℝ) (w : S → ℕ) (g h : ℕ) :
    emptyGap z w (g + h) = compose (emptyGap z w g) (emptyGap z w h) := by
  induction h with
  | zero =>
      funext q r
      simp only [Nat.add_zero]
      change emptyGap z w g q r = compose (emptyGap z w g) identity q r
      exact (compose_identity_right (S := S) (emptyGap z w g) q r).symm
  | succ h ih =>
      calc
        emptyGap z w (g + (h + 1)) =
            compose (emptyGap z w (g + h)) (emptyStep z w) := by
              rw [Nat.add_succ]
              exact emptyGap_succ z w (g + h)
        _ = compose (compose (emptyGap z w g) (emptyGap z w h))
              (emptyStep z w) := by rw [ih]
        _ = compose (emptyGap z w g)
              (compose (emptyGap z w h) (emptyStep z w)) := by
                funext q r
                exact compose_assoc (emptyGap z w g) (emptyGap z w h)
                  (emptyStep z w) q r
        _ = compose (emptyGap z w g) (emptyGap z w (h + 1)) := by
              rw [emptyGap_succ]

/-- Lift the two classes (zero and live) to the state space. -/
def classLift (q : S) (b : Bool) : ℝ :=
  if b then (if q = 0 then 0 else 1) else (if q = 0 then 1 else 0)

/-- Project the live class uniformly onto the nonzero states. -/
def classProject (b : Bool) (q : S) : ℝ :=
  if b then (if q = 0 then 0 else 1 / (liveCount S : ℝ))
  else (if q = 0 then 1 else 0)

theorem project_lift [Nontrivial S] (b d : Bool) :
    (∑ q : S, classProject b q * classLift q d) =
      if b = d then 1 else 0 := by
  cases b <;> cases d
  · simp [classProject, classLift]
  · apply Finset.sum_eq_zero
    intro q _
    by_cases hq : q = 0 <;> simp [classProject, classLift, hq]
  · simp [classProject, classLift]
  · obtain ⟨q0, hq0⟩ := exists_ne (0 : S)
    calc
      (∑ q : S, classProject true q * classLift q true) =
          ∑ q : S, update 0 q0 q := by
            apply Finset.sum_congr rfl
            intro q _
            by_cases hq : q = 0
            · subst q
              simp [classProject, classLift, update_live 0 q0 0 hq0]
            · simp [classProject, classLift, update_live 0 q0 q hq0, hq]
      _ = 1 := update_sum 0 q0

theorem project_live_sum [Nontrivial S] :
    (∑ q : S, classProject true q) = 1 := by
  calc
    (∑ q : S, classProject true q) =
        ∑ q : S, classProject true q * classLift q true := by
          apply Finset.sum_congr rfl
          intro q _
          by_cases hq : q = 0 <;> simp [classProject, classLift, hq]
    _ = 1 := by simpa using project_lift (S := S) true true

theorem project_live_const [Nontrivial S] (f : S → ℝ) (k : ℝ)
    (hf : ∀ q : S, q ≠ 0 → f q = k) :
    (∑ q : S, classProject true q * f q) = k := by
  have hp (q : S) : classProject true q * f q =
      classProject true q * k := by
    by_cases hq : q = 0
    · subst q
      simp [classProject]
    · rw [hf q hq]
  simp_rw [hp]
  rw [← Finset.sum_mul, project_live_sum]
  simp

theorem project_zero_const (f : S → ℝ) :
    (∑ q : S, classProject false q * f q) = f 0 := by
  simp [classProject]

/-- The left factor of a nonempty gap. Its live entry retains the
weight of the actual entering state; the right factor is uniform live. -/
def gapClassFactor (z : ℝ) (w : S → ℕ) (g : ℕ) (q : S) (b : Bool) : ℝ :=
  if b then
    if q = 0 then 0 else z ^ w q * (liveMean z w) ^ g
  else
    if q = 0 then 1 else 0

theorem emptyGap_factor [Nontrivial S] (z : ℝ) (w : S → ℕ)
    (hw : w 0 = 0) (g : ℕ) (q r : S) :
    emptyGap z w (g + 1) q r =
      ∑ b : Bool, gapClassFactor z w g q b * classProject b r := by
  by_cases hq : q = 0
  · subst q
    rw [emptyGap_zero z w hw]
    by_cases hr : r = 0 <;>
      simp [gapClassFactor, classProject, hr]
  · by_cases hr : r = 0
    · subst r
      rw [emptyGap_live_zero z w g q hq]
      simp [gapClassFactor, classProject, hq]
    · rw [emptyGap_live z w g q r hq hr]
      simp [gapClassFactor, classProject, hq, hr]
      ring

theorem update_class_sum (c q : S) (b : Bool) :
    (∑ r : S, update c q r * classLift r b) =
      if b then liveMass c q else update c q 0 := by
  cases b
  · simp [classLift]
  · simp only [ite_true]
    simp [classLift, liveMass]
    have hp (r : S) : (if r = 0 then (0 : ℝ) else update c q r) =
        update c q r - (if r = 0 then update c q 0 else 0) := by
      by_cases hr : r = 0 <;> simp [hr]
    simp_rw [hp]
    rw [Finset.sum_sub_distrib]
    simp

/-- The literal matrix product `R V_c(1) J`: average the incoming state
within its class, apply the one-step state transition, and record the
outgoing class. -/
def classImpulse (c : S) (b d : Bool) : ℝ :=
  ∑ q : S, classProject b q *
    ∑ r : S, update c q r * classLift r d

theorem classImpulse_zero_zero [Nontrivial S] (c : S) (hc : c ≠ 0) :
    classImpulse c false false = 0 := by
  calc
    classImpulse c false false =
        ∑ q : S, classProject false q * update c q 0 := by
          simp [classImpulse, update_class_sum]
    _ = update c 0 0 := project_zero_const _
    _ = 0 := by simp [update_zero, hc, eq_comm]

theorem classImpulse_zero_live [Nontrivial S] (c : S) (hc : c ≠ 0) :
    classImpulse c false true = 1 := by
  calc
    classImpulse c false true =
        ∑ q : S, classProject false q * liveMass c q := by
          simp [classImpulse, update_class_sum]
    _ = liveMass c 0 := project_zero_const _
    _ = 1 := zero_to_live c hc

theorem classImpulse_live_zero [Nontrivial S] (c : S) (hc : c ≠ 0) :
    classImpulse c true false = 1 / (liveCount S : ℝ) := by
  calc
    classImpulse c true false =
        ∑ q : S, classProject true q * update c q 0 := by
          simp [classImpulse, update_class_sum]
    _ = 1 / (liveCount S : ℝ) :=
      project_live_const _ _ (fun q hq => live_to_zero_probability c q hq hc)

theorem classImpulse_live_live [Nontrivial S] (c : S) (hc : c ≠ 0) :
    classImpulse c true true = 1 - 1 / (liveCount S : ℝ) := by
  calc
    classImpulse c true true =
        ∑ q : S, classProject true q * liveMass c q := by
          simp [classImpulse, update_class_sum]
    _ = 1 - 1 / (liveCount S : ℝ) :=
      project_live_const _ _ (fun q hq => live_to_live c q hq hc)

theorem classImpulse_matrix [Nontrivial S] (c : S) (hc : c ≠ 0) :
    ((classImpulse c false false, classImpulse c false true),
      (classImpulse c true false, classImpulse c true true)) =
      (((0 : ℝ), 1),
        (1 / (liveCount S : ℝ), 1 - 1 / (liveCount S : ℝ))) := by
  simp [classImpulse_zero_zero c hc, classImpulse_zero_live c hc,
    classImpulse_live_zero c hc, classImpulse_live_live c hc]

/-- An impulse with state-dependent emitted weight. -/
def weightedImpulse (z : ℝ) (h : S → ℕ) (c q r : S) : ℝ :=
  z ^ h q * update c q r

theorem weightedImpulse_bounds [Nontrivial S] (z : ℝ) (h : S → ℕ)
    (c q r : S) (t : ℕ) (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (hwt : h q ≤ t) :
    z ^ t * update c q r ≤ weightedImpulse z h c q r ∧
      weightedImpulse z h c q r ≤ update c q r := by
  have hu := update_nonneg c q r
  constructor
  · exact mul_le_mul_of_nonneg_right
      (pow_le_pow_of_le_one hz0 hz1 hwt) hu
  · unfold weightedImpulse
    calc
      z ^ h q * update c q r ≤ 1 * update c q r :=
        mul_le_mul_of_nonneg_right (pow_le_one₀ hz0 hz1) hu
      _ = update c q r := one_mul _

/-- The two-state gap approximation lifted back to full states. Its live
mass is `α`, distributed uniformly over the nonzero exit states. -/
def liftedGap (α : ℝ) (q r : S) : ℝ :=
  if q = 0 then (if r = 0 then 1 else 0)
  else if r = 0 then 0 else α / (liveCount S : ℝ)

theorem liftedGap_nonneg (α : ℝ) (hα : 0 ≤ α) (q r : S) :
    0 ≤ liftedGap α q r := by
  by_cases hq : q = 0 <;> by_cases hr : r = 0 <;>
    simp [liftedGap, hq, hr, div_nonneg hα (Nat.cast_nonneg _)]

theorem liftedGap_as_JDR (α : ℝ) (q r : S) :
    liftedGap α q r =
      ∑ b : Bool, classLift q b * (if b then α else 1) * classProject b r := by
  by_cases hq : q = 0 <;> by_cases hr : r = 0 <;>
    simp [liftedGap, classLift, classProject, hq, hr, div_eq_mul_inv]

/-- The algebraic part of the nonempty-gap upper comparison. The analytic
work is isolated in `hfactor`, a bound on the actual entering-state weight
and the live-state mean. -/
theorem emptyGap_upper_from_factor [Nontrivial S]
    (z : ℝ) (w : S → ℕ) (hw : w 0 = 0) (g : ℕ)
    (α a : ℝ) (ha : 1 ≤ a)
    (hfactor : ∀ q : S, q ≠ 0 →
      z ^ w q * (liveMean z w) ^ g ≤ a * α)
    (q r : S) :
    emptyGap z w (g + 1) q r ≤ a * liftedGap α q r := by
  by_cases hq : q = 0
  · subst q
    rw [emptyGap_zero z w hw]
    by_cases hr : r = 0 <;> simp [liftedGap, hr, ha]
  · by_cases hr : r = 0
    · subst r
      rw [emptyGap_live_zero z w g q hq]
      simp [liftedGap, hq]
    · rw [emptyGap_live z w g q r hq hr]
      simp only [liftedGap, hq, hr, ite_false]
      have hd : 0 ≤ (1 / (liveCount S : ℝ)) :=
        one_div_nonneg.mpr (Nat.cast_nonneg _)
      have hh := mul_le_mul_of_nonneg_right (hfactor q hq) hd
      simp only [div_eq_mul_inv] at hh ⊢
      nlinarith [hh]

theorem emptyGap_lower_from_factor [Nontrivial S]
    (z : ℝ) (w : S → ℕ) (hw : w 0 = 0) (g : ℕ)
    (α b : ℝ) (hb : b ≤ 1)
    (hfactor : ∀ q : S, q ≠ 0 →
      b * α ≤ z ^ w q * (liveMean z w) ^ g)
    (q r : S) :
    b * liftedGap α q r ≤ emptyGap z w (g + 1) q r := by
  by_cases hq : q = 0
  · subst q
    rw [emptyGap_zero z w hw]
    by_cases hr : r = 0 <;> simp [liftedGap, hr, hb]
  · by_cases hr : r = 0
    · subst r
      rw [emptyGap_live_zero z w g q hq]
      simp [liftedGap, hq]
    · rw [emptyGap_live z w g q r hq hr]
      simp only [liftedGap, hq, hr, ite_false]
      have hd : 0 ≤ (1 / (liveCount S : ℝ)) :=
        one_div_nonneg.mpr (Nat.cast_nonneg _)
      have hh := mul_le_mul_of_nonneg_right (hfactor q hq) hd
      simp only [div_eq_mul_inv] at hh ⊢
      nlinarith [hh]

end Spin.UniformMixer
