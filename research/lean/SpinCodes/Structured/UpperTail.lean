/-
The upper-tail counting step of `lem:structured-ba-tails`.

Conditioning on the first bit of a uniform weight-`v` word gives

    Pr[wt(Acc V) ≥ (1-δ_o)b] ≤ (v/(b-v+1))·p_{v-1} + ((b-v)/(v+1))·p_{v+1}.

The mechanism is `Acc(V + e₁) = ¬Acc(V)` (`accWtL_flipHead`): flipping the
first bit turns the upper tail into a lower tail, and moves the input weight
by one in the direction set by that bit.  So the upper-tail words inject into
the union of the weight-`(v-1)` and weight-`(v+1)` lower-tail words, and the
count inequality follows with no probability and no division.

Flipping the first bit is an involution, hence injective — which is all the
counting argument needs.  Working over `Fin (b+1) → Bool` rather than lists
makes that a one-line `Function.update` fact.

This is a leaf module: it does not invalidate the expensive Golay computation.
-/
import SpinCodes.Structured.AccTuple
import SpinCodes.Structured.BATails

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-- Add `e₁` to a tuple: flip the first coordinate. -/
def flipHeadF {b : ℕ} (v : Fin (b + 1) → Bool) : Fin (b + 1) → Bool :=
  Function.update v 0 (!(v 0))

@[simp] lemma flipHeadF_zero {b : ℕ} (v : Fin (b + 1) → Bool) :
    flipHeadF v 0 = !(v 0) := by
  simp [flipHeadF]

lemma flipHeadF_succ {b : ℕ} (v : Fin (b + 1) → Bool) (i : Fin b) :
    flipHeadF v i.succ = v i.succ := by
  simp [flipHeadF, Function.update, Fin.succ_ne_zero]

/-- Flipping twice is the identity, so `flipHeadF` is injective. -/
lemma flipHeadF_involutive {b : ℕ} : Function.Involutive (flipHeadF (b := b)) := by
  intro v
  funext i
  by_cases h : i = 0
  · subst h; simp
  · have : flipHeadF v i = v i := by
      simp [flipHeadF, Function.update, h]
    simp [flipHeadF, Function.update, h, this]

lemma flipHeadF_injective {b : ℕ} : Function.Injective (flipHeadF (b := b)) :=
  flipHeadF_involutive.injective

lemma ofFn_flipHeadF {b : ℕ} (v : Fin (b + 1) → Bool) :
    List.ofFn (flipHeadF v) = flipHead (List.ofFn v) := by
  rw [List.ofFn_succ, List.ofFn_succ]
  simp [flipHead, flipHeadF_succ]

/-- The accumulated weight is complemented. -/
theorem accWtF_flipHeadF {b : ℕ} (v : Fin (b + 1) → Bool) :
    accWtF (flipHeadF v) = (b + 1) - accWtF v := by
  unfold accWtF
  rw [ofFn_flipHeadF, accWtL_flipHead]
  simp

/-- The input weight moves by one, in the direction set by the first bit. -/
theorem wtF_flipHeadF {b : ℕ} (v : Fin (b + 1) → Bool) :
    wtF (flipHeadF v) = if v 0 then wtF v - 1 else wtF v + 1 := by
  unfold wtF
  rw [ofFn_flipHeadF, List.ofFn_succ]
  cases h : v 0 <;>
    simp [flipHead, wtL, h, List.ofFn_succ] <;> omega

/-- **The upper-tail counting step.**

Upper-tail words of weight `v` inject, via `V ↦ V + e₁`, into the union of the
lower-tail words of weights `v-1` and `v+1`.  Stated as a count, with no
probability and no division. -/
theorem card_upper_tail_le {b D v : ℕ} (hv : 1 ≤ v) (hD : D ≤ b + 1) :
    (univ.filter (fun w : Fin (b + 1) → Bool =>
        wtF w = v ∧ (b + 1) - D ≤ accWtF w)).card
      ≤ (univ.filter (fun w : Fin (b + 1) → Bool =>
            wtF w = v - 1 ∧ accWtF w ≤ D)).card
        + (univ.filter (fun w : Fin (b + 1) → Bool =>
            wtF w = v + 1 ∧ accWtF w ≤ D)).card := by
  classical
  set A := univ.filter (fun w : Fin (b + 1) → Bool => wtF w = v - 1 ∧ accWtF w ≤ D) with hA
  set B := univ.filter (fun w : Fin (b + 1) → Bool => wtF w = v + 1 ∧ accWtF w ≤ D) with hB
  have hdisj : Disjoint A B := by
    rw [Finset.disjoint_left]
    intro w hwA hwB
    simp only [hA, hB, Finset.mem_filter] at hwA hwB
    omega
  have hcard : A.card + B.card = (A ∪ B).card := (Finset.card_union_of_disjoint hdisj).symm
  rw [hcard]
  refine Finset.card_le_card_of_injOn flipHeadF ?_ (fun x _ y _ h => flipHeadF_injective h)
  intro w hw
  obtain ⟨hwt, hacc⟩ := (Finset.mem_filter.mp hw).2
  have hle : accWtF w ≤ b + 1 := by
    rw [← wtF_accF]
    exact wtF_le _
  have haccF : accWtF (flipHeadF w) ≤ D := by
    rw [accWtF_flipHeadF]
    omega
  refine Finset.mem_union.mpr ?_
  have hflip := wtF_flipHeadF w
  cases h : w 0
  · right
    rw [hB]
    refine Finset.mem_filter.mpr ⟨Finset.mem_univ _, ?_, haccF⟩
    rw [hflip, h, hwt]
    simp
  · left
    rw [hA]
    refine Finset.mem_filter.mpr ⟨Finset.mem_univ _, ?_, haccF⟩
    rw [hflip, h, hwt]
    simp


/-! ## From counts to the paper's coefficients

The paper writes the upper-tail bound with the coefficients
`C(b,v-1)/C(b,v) = v/(b-v+1)` and `C(b,v+1)/C(b,v) = (b-v)/(v+1)`.  Both are
the *same* Mathlib identity, `Nat.choose_succ_right_eq`, read at `k = v-1` and
`k = v`.  Since each coefficient is at most `b`, the composite bound only needs
`C(b, v±1) ≤ b · C(b, v)`. -/

/-- `C(b,v)·v = C(b,v-1)·(b-(v-1))` — the paper's `v/(b-v+1)` coefficient. -/
lemma choose_pred_ratio {b v : ℕ} (hv : 1 ≤ v) :
    b.choose v * v = b.choose (v - 1) * (b - (v - 1)) := by
  have h := Nat.choose_succ_right_eq b (v - 1)
  rwa [show v - 1 + 1 = v from by omega] at h

/-- `C(b,v+1)·(v+1) = C(b,v)·(b-v)` — the paper's `(b-v)/(v+1)` coefficient. -/
lemma choose_succ_ratio (b v : ℕ) :
    b.choose (v + 1) * (v + 1) = b.choose v * (b - v) :=
  Nat.choose_succ_right_eq b v

/-- Each coefficient is at most `b`, so `C(b,v-1) ≤ b·C(b,v)`. -/
lemma choose_pred_le_mul {b v : ℕ} (hv : 1 ≤ v) (hvb : v ≤ b) :
    b.choose (v - 1) ≤ b * b.choose v := by
  have hr := choose_pred_ratio (b := b) hv
  have hge : 1 ≤ b - (v - 1) := by omega
  calc b.choose (v - 1) ≤ b.choose (v - 1) * (b - (v - 1)) :=
        Nat.le_mul_of_pos_right _ (by omega)
    _ = b.choose v * v := hr.symm
    _ ≤ b.choose v * b := Nat.mul_le_mul_left _ hvb
    _ = b * b.choose v := by ring

/-- And `C(b,v+1) ≤ b·C(b,v)`. -/
lemma choose_succ_le_mul (b v : ℕ) : b.choose (v + 1) ≤ b * b.choose v := by
  calc b.choose (v + 1) ≤ b.choose (v + 1) * (v + 1) :=
        Nat.le_mul_of_pos_right _ (by omega)
    _ = b.choose v * (b - v) := choose_succ_ratio b v
    _ ≤ b.choose v * b := Nat.mul_le_mul_left _ (by omega)
    _ = b * b.choose v := by ring

/-! ## The lower-tail count as a cardinality

`Σ_{h ≤ D} T_b(u,h)` counts the words of weight `u` whose accumulation has
weight at most `D`. -/

theorem sum_accT_eq_card (b u D : ℕ) :
    ∑ h ∈ range (D + 1), accT b u h
      = (univ.filter (fun w : Fin b → Bool => wtF w = u ∧ accWtF w ≤ D)).card := by
  classical
  have hmaps : ∀ w ∈ univ.filter (fun w : Fin b → Bool => wtF w = u ∧ accWtF w ≤ D),
      accWtF w ∈ range (D + 1) := by
    intro w hw
    obtain ⟨-, -, hle⟩ := Finset.mem_filter.mp hw
    rw [Finset.mem_range]
    omega
  rw [Finset.card_eq_sum_card_fiberwise hmaps]
  refine Finset.sum_congr rfl fun h hh => ?_
  rw [Finset.mem_range] at hh
  rw [← accCountF_eq_accT]
  unfold accCountF
  refine congrArg Finset.card ?_
  ext w
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  constructor
  · rintro ⟨hw, hacc⟩
    exact ⟨⟨hw, by omega⟩, hacc⟩
  · rintro ⟨⟨hw, -⟩, hacc⟩
    exact ⟨hw, hacc⟩


/-- **The upper tail against the two neighbouring lower tails.**

`card_upper_tail_le` phrased with the lower tails as `accT` sums, which is the
form `eq:structured-ba-one-acc-tail` bounds.  In the paper this appears divided
by `C(b,v)`, producing the coefficients `v/(b-v+1)` and `(b-v)/(v+1)`; those
are `choose_pred_ratio` and `choose_succ_ratio` above, and each is at most `b`. -/
theorem upper_tail_le_sums {b D v : ℕ} (hv : 1 ≤ v) (hD : D ≤ b + 1) :
    (univ.filter (fun w : Fin (b + 1) → Bool =>
        wtF w = v ∧ (b + 1) - D ≤ accWtF w)).card
      ≤ (∑ h ∈ range (D + 1), accT (b + 1) (v - 1) h)
        + (∑ h ∈ range (D + 1), accT (b + 1) (v + 1) h) := by
  rw [sum_accT_eq_card, sum_accT_eq_card]
  exact card_upper_tail_le hv hD


/-! ## The upper tail, numerically

Substituting `eq:structured-ba-one-acc-tail` into `upper_tail_le_sums`.  Both
neighbours `v∓1` share a parity, so this is a single case split rather than
four.  The odd-`v` case is stated here; it is edge-case free, because `v = 2m+1`
gives neighbours `2m` and `2m+2` with `ℓ = m` and `ℓ = m+1`, both at least one.
(The even case meets `ℓ = 0` at `v = 2`, where the tail lemma does not apply
and the bound is trivial instead.) -/

/-- **The upper tail at odd weight.**

`U_{2m+1} · 28^{m+1} ≤ 41·13^m · b·C(b, 2m+1)`, which is the paper's
`b(ζ⁻²+1)ζ^v` with `ζ² = 13/28`: the `41 = 28 + 13` collects the two
neighbouring contributions. -/
theorem upper_tail_odd_bound {b D m : ℕ} (hm : 1 ≤ m) (hmD : m + 1 ≤ D)
    (hDb : 125 * D ≤ 13 * (b + 1)) (hDle : D ≤ b + 1) :
    (univ.filter (fun w : Fin (b + 1) → Bool =>
        wtF w = 2 * m + 1 ∧ (b + 1) - D ≤ accWtF w)).card * 28 ^ (m + 1)
      ≤ 41 * 13 ^ m * ((b + 1) * (b + 1).choose (2 * m + 1)) := by
  set B := b + 1 with hB
  set U := (univ.filter (fun w : Fin (b + 1) → Bool =>
    wtF w = 2 * m + 1 ∧ (b + 1) - D ≤ accWtF w)).card with hU
  have hvB : 2 * m + 1 ≤ B := by omega
  -- the two neighbouring lower-tail sums
  have hsplit : U ≤ (∑ h ∈ range (D + 1), accT B (2 * m) h)
      + (∑ h ∈ range (D + 1), accT B (2 * m + 2) h) := by
    have h := upper_tail_le_sums (b := b) (D := D) (v := 2 * m + 1) (by omega) hDle
    rw [show 2 * m + 1 - 1 = 2 * m from by omega] at h
    exact h
  -- tail bounds at the two neighbours
  have hlow : (∑ h ∈ range (D + 1), accT B (2 * m) h) * 28 ^ m
      ≤ 13 ^ m * B.choose (2 * m) := sum_accT_even_le_pow hm (by omega) hDb
  have hhigh : (∑ h ∈ range (D + 1), accT B (2 * (m + 1)) h) * 28 ^ (m + 1)
      ≤ 13 ^ (m + 1) * B.choose (2 * (m + 1)) :=
    sum_accT_even_le_pow (by omega) (by omega) hDb
  rw [show 2 * (m + 1) = 2 * m + 2 from by ring] at hhigh
  -- both neighbouring binomials are at most `B · C(B, 2m+1)`
  have hpred : B.choose (2 * m) ≤ B * B.choose (2 * m + 1) := by
    have h := choose_pred_le_mul (b := B) (v := 2 * m + 1) (by omega) hvB
    rwa [show 2 * m + 1 - 1 = 2 * m from by omega] at h
  have hsucc : B.choose (2 * m + 2) ≤ B * B.choose (2 * m + 1) := by
    have h := choose_succ_le_mul B (2 * m + 1)
    rwa [show 2 * m + 1 + 1 = 2 * m + 2 from by omega] at h
  -- assemble
  calc U * 28 ^ (m + 1)
      ≤ ((∑ h ∈ range (D + 1), accT B (2 * m) h)
          + (∑ h ∈ range (D + 1), accT B (2 * m + 2) h)) * 28 ^ (m + 1) :=
        Nat.mul_le_mul_right _ hsplit
    _ = 28 * ((∑ h ∈ range (D + 1), accT B (2 * m) h) * 28 ^ m)
          + (∑ h ∈ range (D + 1), accT B (2 * m + 2) h) * 28 ^ (m + 1) := by ring
    _ ≤ 28 * (13 ^ m * B.choose (2 * m))
          + 13 ^ (m + 1) * B.choose (2 * m + 2) := by
        exact Nat.add_le_add (Nat.mul_le_mul_left _ hlow) hhigh
    _ ≤ 28 * (13 ^ m * (B * B.choose (2 * m + 1)))
          + 13 ^ (m + 1) * (B * B.choose (2 * m + 1)) := by
        exact Nat.add_le_add (Nat.mul_le_mul_left _ (Nat.mul_le_mul_left _ hpred))
          (Nat.mul_le_mul_left _ hsucc)
    _ = 41 * 13 ^ m * (B * B.choose (2 * m + 1)) := by ring


/-- **The upper tail at even weight.**

`U_{2(k+1)} · 28^{k+1} ≤ 41·13^k · b·C(b, 2(k+1))`.  Writing the weight as
`2(k+1)` rather than `2m` keeps both neighbours at `ℓ = k` and `ℓ = k+1` with
no truncated subtraction, and `sum_accT_odd_le_pow'` covers `ℓ = 0`, so the
`v = 2` case needs no separate treatment. -/
theorem upper_tail_even_bound {b D k : ℕ} (hkD : k + 1 ≤ D)
    (hDb : 125 * D ≤ 13 * (b + 1)) (hDle : D ≤ b + 1) :
    (univ.filter (fun w : Fin (b + 1) → Bool =>
        wtF w = 2 * (k + 1) ∧ (b + 1) - D ≤ accWtF w)).card * 28 ^ (k + 1)
      ≤ 41 * 13 ^ k * ((b + 1) * (b + 1).choose (2 * (k + 1))) := by
  set B := b + 1 with hB
  set U := (univ.filter (fun w : Fin (b + 1) → Bool =>
    wtF w = 2 * (k + 1) ∧ (b + 1) - D ≤ accWtF w)).card with hU
  have hvB : 2 * (k + 1) ≤ B := by omega
  have hsplit : U ≤ (∑ h ∈ range (D + 1), accT B (2 * k + 1) h)
      + (∑ h ∈ range (D + 1), accT B (2 * (k + 1) + 1) h) := by
    have h := upper_tail_le_sums (b := b) (D := D) (v := 2 * (k + 1)) (by omega) hDle
    rw [show 2 * (k + 1) - 1 = 2 * k + 1 from by omega] at h
    exact h
  have hlow : (∑ h ∈ range (D + 1), accT B (2 * k + 1) h) * 28 ^ k
      ≤ 13 ^ k * B.choose (2 * k + 1) := sum_accT_odd_le_pow' (by omega) hDb
  have hhigh : (∑ h ∈ range (D + 1), accT B (2 * (k + 1) + 1) h) * 28 ^ (k + 1)
      ≤ 13 ^ (k + 1) * B.choose (2 * (k + 1) + 1) :=
    sum_accT_odd_le_pow' (by omega) hDb
  have hpred : B.choose (2 * k + 1) ≤ B * B.choose (2 * (k + 1)) := by
    have h := choose_pred_le_mul (b := B) (v := 2 * (k + 1)) (by omega) hvB
    rwa [show 2 * (k + 1) - 1 = 2 * k + 1 from by omega] at h
  have hsucc : B.choose (2 * (k + 1) + 1) ≤ B * B.choose (2 * (k + 1)) :=
    choose_succ_le_mul B (2 * (k + 1))
  calc U * 28 ^ (k + 1)
      ≤ ((∑ h ∈ range (D + 1), accT B (2 * k + 1) h)
          + (∑ h ∈ range (D + 1), accT B (2 * (k + 1) + 1) h)) * 28 ^ (k + 1) :=
        Nat.mul_le_mul_right _ hsplit
    _ = 28 * ((∑ h ∈ range (D + 1), accT B (2 * k + 1) h) * 28 ^ k)
          + (∑ h ∈ range (D + 1), accT B (2 * (k + 1) + 1) h) * 28 ^ (k + 1) := by ring
    _ ≤ 28 * (13 ^ k * B.choose (2 * k + 1))
          + 13 ^ (k + 1) * B.choose (2 * (k + 1) + 1) :=
        Nat.add_le_add (Nat.mul_le_mul_left _ hlow) hhigh
    _ ≤ 28 * (13 ^ k * (B * B.choose (2 * (k + 1))))
          + 13 ^ (k + 1) * (B * B.choose (2 * (k + 1))) :=
        Nat.add_le_add (Nat.mul_le_mul_left _ (Nat.mul_le_mul_left _ hpred))
          (Nat.mul_le_mul_left _ hsucc)
    _ = 41 * 13 ^ k * (B * B.choose (2 * (k + 1))) := by ring

end Spin.Structured
