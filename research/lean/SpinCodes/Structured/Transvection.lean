/-
The IMT transvection law.

Setup samples `u ← 𝔽₂^s \ {0}` and then `v ← {v : ⟨u,v⟩ = 0}`, and sets
`M := I + u vᵀ`.  The paper records, for each fixed `q ≠ 0`,

    L(M q) = ½ δ_q + ½ Unif(𝔽₂^s \ {0}),

with the one-line reason: `u = q` always fixes `q`, and for every other `u`
the bit `⟨v,q⟩` is fair.  Every entry of the transfer matrices in
`app:imt-finite-transfers` is built on this law -- the `1/2` on the lazy
branch and the `1/(2M)` per target on the refresh branch are exactly its two
halves -- so it is proved here before the matrices are assembled.

As elsewhere in this development, vectors are subsets: `⟨a,b⟩` is the parity
of `|a ∩ b|` and addition is `∆`.

The law is stated as three exact counts over the sample space, which avoids
introducing `2^(s-1)` and so has no edge cases in `s`:

*  zero is unreachable from a nonzero `q`;
*  every other nonzero target is hit by half of one `perp`;
*  `q` itself is hit by all of `perp q` plus half of every other `perp`.

Together with `card_perp_double` these give `½ + 1/(2M)` and `1/(2M)`.
-/
import SpinCodes.Structured.VarianceBound

namespace Spin

open Finset
open scoped symmDiff

variable {s : ℕ}

/-! ## The action and the sample space -/

/-- The transvection `I + u vᵀ` applied to `q`: it fixes `q` unless
`⟨v,q⟩ = 1`, in which case it adds `u`. -/
def act (u v q : Finset (Fin s)) : Finset (Fin s) :=
  if Even (q ∩ v).card then q else q ∆ u

/-- `u^⊥ = {v : ⟨u,v⟩ = 0}`. -/
def perp (u : Finset (Fin s)) : Finset (Finset (Fin s)) :=
  univ.filter (fun v => Even (u ∩ v).card)

/-- The sample space: `u` nonzero, `v` orthogonal to `u`. -/
def transPairs (s : ℕ) : Finset (Finset (Fin s) × Finset (Fin s)) :=
  univ.filter (fun p => p.1 ≠ ∅ ∧ Even (p.1 ∩ p.2).card)

/-- The nonzero states. -/
def nonzeroStates (s : ℕ) : Finset (Finset (Fin s)) :=
  univ.filter (fun u : Finset (Fin s) => u ≠ ∅)

/-- Counting over the sample space row by row. -/
lemma card_filter_transPairs (P : Finset (Fin s) × Finset (Fin s) → Prop)
    [DecidablePred P] :
    ((transPairs s).filter P).card
      = ∑ u ∈ nonzeroStates s, ((perp u).filter (fun v => P (u, v))).card := by
  classical
  have key : ∀ u : Finset (Fin s),
      (if u ≠ ∅ then ((perp u).filter (fun v => P (u, v))).card else 0)
        = ∑ v : Finset (Fin s),
            (if (u ≠ ∅ ∧ Even (u ∩ v).card) ∧ P (u, v) then 1 else 0) := by
    intro u
    by_cases hu : u ≠ ∅
    · rw [if_pos hu, perp, Finset.filter_filter, Finset.card_filter]
      exact Finset.sum_congr rfl fun v _ => by simp [hu, and_assoc]
    · rw [if_neg hu]
      exact (Finset.sum_eq_zero fun v _ => by simp [hu]).symm
  rw [transPairs, Finset.filter_filter, Finset.card_filter, Fintype.sum_prod_type,
    nonzeroStates, Finset.sum_filter]
  exact Finset.sum_congr rfl fun u _ => (key u).symm

/-! ## Parity is additive under `∆` -/

lemma even_inter_symmDiff (t v w : Finset (Fin s)) :
    (Even (t ∩ (v ∆ w)).card ↔ (Even (t ∩ v).card ↔ Even (t ∩ w).card)) := by
  have hd : t ∩ (v ∆ w) = (t ∩ v) ∆ (t ∩ w) := inter_symmDiff_right t v w
  have hc := card_symmDiff_add_two_mul_card_inter (t ∩ v) (t ∩ w)
  rw [hd, Nat.even_iff, Nat.even_iff, Nat.even_iff]
  omega

/-! ## A parity functional splits a `∆`-closed family in half -/

/-- If `w` has odd overlap with `t` and translating by `w` preserves `S`, then
`t`'s parity splits `S` evenly. -/
theorem card_split_by_parity (S : Finset (Finset (Fin s))) (t w : Finset (Fin s))
    (hclosed : ∀ v ∈ S, v ∆ w ∈ S) (hw : ¬ Even (t ∩ w).card) :
    (S.filter (fun v => Even (t ∩ v).card)).card
      = (S.filter (fun v => ¬ Even (t ∩ v).card)).card := by
  classical
  refine Finset.card_nbij' (fun v => v ∆ w) (fun v => v ∆ w) ?_ ?_ ?_ ?_
  · intro v hv
    simp only [Finset.mem_coe, Finset.mem_filter] at hv ⊢
    refine ⟨hclosed v hv.1, ?_⟩
    rw [even_inter_symmDiff]
    intro h
    exact hw (h.mp hv.2)
  · intro v hv
    simp only [Finset.mem_coe, Finset.mem_filter] at hv ⊢
    refine ⟨hclosed v hv.1, ?_⟩
    rw [even_inter_symmDiff]
    exact ⟨fun h => absurd h hv.2, fun h => absurd h hw⟩
  · intro v _
    exact symmDiff_cancel_tail v w
  · intro v _
    exact symmDiff_cancel_tail v w

/-! ## Two applications of the split -/

/-- `|u^⊥| = 2^(s-1)`, stated as `2·|u^⊥| = 2^s`. -/
theorem card_perp_double {u : Finset (Fin s)} (hu : u ≠ ∅) :
    2 * (perp u).card = 2 ^ s := by
  classical
  obtain ⟨a, ha⟩ := Finset.nonempty_iff_ne_empty.mpr hu
  have hodd : ¬ Even (u ∩ ({a} : Finset (Fin s))).card := by
    have : u ∩ ({a} : Finset (Fin s)) = {a} := by
      ext b; simp only [Finset.mem_inter, Finset.mem_singleton]
      exact ⟨fun h => h.2, fun h => ⟨h ▸ ha, h⟩⟩
    rw [this, Finset.card_singleton]
    decide
  have hsplit := card_split_by_parity (univ : Finset (Finset (Fin s))) u {a}
    (fun v _ => Finset.mem_univ _) hodd
  have htotal :
      ((univ : Finset (Finset (Fin s))).filter (fun v => Even (u ∩ v).card)).card
        + ((univ : Finset (Finset (Fin s))).filter
            (fun v => ¬ Even (u ∩ v).card)).card
        = (univ : Finset (Finset (Fin s))).card :=
    Finset.card_filter_add_card_filter_not _
  have hcard : (univ : Finset (Finset (Fin s))).card = 2 ^ s := by
    rw [Finset.card_univ, Fintype.card_finset, Fintype.card_fin]
  rw [hcard] at htotal
  unfold perp
  omega

/-- A vector orthogonal to `u` but not to `t`, when `t ∉ {∅, u}`. -/
theorem exists_perp_odd {u t : Finset (Fin s)} (ht : t ≠ ∅) (hne : t ≠ u) :
    ∃ w, Even (u ∩ w).card ∧ ¬ Even (t ∩ w).card := by
  classical
  by_cases hsub : t ⊆ u
  · -- `t ⊊ u`: take one element of `t` and one of `u \ t`
    obtain ⟨a, ha⟩ := Finset.nonempty_iff_ne_empty.mpr ht
    have hlt : t ⊂ u := lt_of_le_of_ne hsub hne
    obtain ⟨b, hbu, hbt⟩ := Finset.exists_of_ssubset hlt
    have hab : a ≠ b := fun h => hbt (h ▸ ha)
    refine ⟨{a, b}, ?_, ?_⟩
    · have : u ∩ ({a, b} : Finset (Fin s)) = {a, b} := by
        ext c
        simp only [Finset.mem_inter, Finset.mem_insert, Finset.mem_singleton]
        exact ⟨fun h => h.2, fun h => ⟨h.elim (fun hc => hc ▸ hsub ha) (fun hc => hc ▸ hbu), h⟩⟩
      rw [this, Finset.card_insert_of_notMem (by simpa using hab), Finset.card_singleton]
      decide
    · have : t ∩ ({a, b} : Finset (Fin s)) = {a} := by
        ext c
        simp only [Finset.mem_inter, Finset.mem_insert, Finset.mem_singleton]
        constructor
        · rintro ⟨hct, hc | hc⟩
          · exact hc
          · exact absurd (hc ▸ hct) hbt
        · rintro rfl
          exact ⟨ha, Or.inl rfl⟩
      rw [this, Finset.card_singleton]
      decide
  · -- some element of `t` lies outside `u`
    obtain ⟨a, hat, hau⟩ := Finset.not_subset.mp hsub
    refine ⟨{a}, ?_, ?_⟩
    · have : u ∩ ({a} : Finset (Fin s)) = ∅ := by
        ext c
        simp only [Finset.mem_inter, Finset.mem_singleton, Finset.notMem_empty, iff_false]
        rintro ⟨hcu, rfl⟩
        exact hau hcu
      rw [this, Finset.card_empty]
      decide
    · have : t ∩ ({a} : Finset (Fin s)) = {a} := by
        ext c; simp only [Finset.mem_inter, Finset.mem_singleton]
        exact ⟨fun h => h.2, fun h => ⟨h ▸ hat, h⟩⟩
      rw [this, Finset.card_singleton]
      decide

/-- `perp u` is closed under translation by any of its members. -/
lemma perp_symmDiff_mem {u v w : Finset (Fin s)} (hv : v ∈ perp u) (hw : w ∈ perp u) :
    v ∆ w ∈ perp u := by
  simp only [perp, Finset.mem_filter, Finset.mem_univ, true_and] at hv hw ⊢
  rw [even_inter_symmDiff]
  exact ⟨fun _ => hw, fun _ => hv⟩

/-- **The fair-bit step.**  For `u ≠ q` with `q` nonzero, `⟨v,q⟩` is fair on
`u^⊥`. -/
theorem card_perp_fair {u q : Finset (Fin s)} (hq : q ≠ ∅) (hne : q ≠ u) :
    2 * ((perp u).filter (fun v => ¬ Even (q ∩ v).card)).card = (perp u).card := by
  classical
  obtain ⟨w, hwu, hwq⟩ := exists_perp_odd hq hne
  have hwmem : w ∈ perp u := by
    simp only [perp, Finset.mem_filter, Finset.mem_univ, true_and]
    exact hwu
  have hsplit := card_split_by_parity (perp u) q w
    (fun v hv => perp_symmDiff_mem hv hwmem) hwq
  have htotal :
      ((perp u).filter (fun v => Even (q ∩ v).card)).card
        + ((perp u).filter (fun v => ¬ Even (q ∩ v).card)).card = (perp u).card :=
    Finset.card_filter_add_card_filter_not _
  omega

/-! ## The law -/

lemma symmDiff_empty_right (a : Finset (Fin s)) : a ∆ (∅ : Finset (Fin s)) = a := by
  rw [show (∅ : Finset (Fin s)) = ⊥ from rfl, symmDiff_bot]

lemma symmDiff_eq_iff {q u q' : Finset (Fin s)} : q ∆ u = q' ↔ u = q ∆ q' := by
  constructor
  · rintro rfl
    rw [symmDiff_symmDiff_cancel_left]
  · rintro rfl
    rw [symmDiff_symmDiff_cancel_left]

/-- **Zero is unreachable** from a nonzero state. -/
theorem count_target_zero {q : Finset (Fin s)} (hq : q ≠ ∅) :
    ((transPairs s).filter (fun p => act p.1 p.2 q = ∅)).card = 0 := by
  classical
  rw [card_filter_transPairs]
  refine Finset.sum_eq_zero fun u hu => ?_
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro v hv
  simp only [perp, Finset.mem_filter, Finset.mem_univ, true_and] at hv
  unfold act
  by_cases hbit : Even (q ∩ v).card
  · simpa [hbit] using hq
  · simp only [hbit, if_false]
    rw [symmDiff_eq_iff]
    intro hcon
    -- `u = q ∆ ∅ = q`, but then `v ⊥ q` forces the bit even
    rw [symmDiff_empty_right] at hcon
    subst hcon
    exact hbit hv

/-- **Every other nonzero target** is hit by exactly half of one `u^⊥`. -/
theorem count_target_other {q q' : Finset (Fin s)} (hq : q ≠ ∅) (hq' : q' ≠ ∅)
    (hne : q' ≠ q) :
    2 * ((transPairs s).filter (fun p => act p.1 p.2 q = q')).card
      = (perp (q ∆ q')).card := by
  classical
  rw [card_filter_transPairs, Finset.mul_sum]
  have hu0 : q ∆ q' ≠ ∅ := by
    intro h
    exact hne (symmDiff_eq_bot.mp (by simpa using h)).symm
  have hmem : q ∆ q' ∈ nonzeroStates s := by
    simp only [nonzeroStates, Finset.mem_filter, Finset.mem_univ, true_and]
    exact hu0
  rw [← Finset.add_sum_erase _ _ hmem]
  have hzero : ∀ u ∈ (nonzeroStates s).erase (q ∆ q'),
      2 * ((perp u).filter (fun v => act u v q = q')).card = 0 := by
    intro u hu
    have hune : u ≠ q ∆ q' := (Finset.mem_erase.mp hu).1
    have hempty : ((perp u).filter (fun v => act u v q = q')) = ∅ := by
      rw [Finset.filter_eq_empty_iff]
      intro v _
      unfold act
      by_cases hbit : Even (q ∩ v).card
      · simp only [hbit, if_true]
        exact fun h => hne h.symm
      · simp only [hbit, if_false]
        rw [symmDiff_eq_iff]
        exact hune
    rw [hempty, Finset.card_empty, Nat.mul_zero]
  rw [Finset.sum_congr rfl hzero, Finset.sum_const_zero, Nat.add_zero]
  -- at `u = q ∆ q'` the condition is exactly the odd branch
  have hcond : ∀ v ∈ perp (q ∆ q'),
      (act (q ∆ q') v q = q') ↔ ¬ Even (q ∩ v).card := by
    intro v _
    unfold act
    by_cases hbit : Even (q ∩ v).card
    · simp only [hbit, if_true, not_true]
      exact ⟨fun h => hne h.symm, False.elim⟩
    · simp only [hbit, if_false, not_false_iff, iff_true]
      rw [symmDiff_eq_iff]
  rw [Finset.filter_congr hcond]
  have hqne : q ≠ q ∆ q' := by
    intro h
    have hc := symmDiff_symmDiff_cancel_left q q'
    rw [← h, symmDiff_self] at hc
    exact hq' (by simpa using hc.symm)
  exact card_perp_fair hq hqne

/-- **The lazy branch.**  `q` keeps all of `perp q` and half of every other
`perp`. -/
theorem count_target_self {q : Finset (Fin s)} (hq : q ≠ ∅) :
    2 * ((transPairs s).filter (fun p => act p.1 p.2 q = q)).card
      = 2 * (perp q).card
        + ∑ u ∈ (nonzeroStates s).erase q, (perp u).card := by
  classical
  rw [card_filter_transPairs, Finset.mul_sum]
  have hmem : q ∈ nonzeroStates s := by
    simp only [nonzeroStates, Finset.mem_filter, Finset.mem_univ, true_and]
    exact hq
  rw [← Finset.add_sum_erase _ _ hmem]
  congr 1
  · -- at `u = q` every `v ∈ perp q` has an even bit
    have hall : ∀ v ∈ perp q, act q v q = q := by
      intro v hv
      simp only [perp, Finset.mem_filter, Finset.mem_univ, true_and] at hv
      unfold act
      simp [hv]
    rw [Finset.filter_true_of_mem hall]
  · refine Finset.sum_congr rfl fun u hu => ?_
    have hune : u ≠ q := (Finset.mem_erase.mp hu).1
    have hu0 : u ≠ ∅ := by
      have hmem2 := (Finset.mem_erase.mp hu).2
      simpa [nonzeroStates] using hmem2
    have hcond : ∀ v ∈ perp u, (act u v q = q) ↔ Even (q ∩ v).card := by
      intro v _
      unfold act
      by_cases hbit : Even (q ∩ v).card
      · simp [hbit]
      · simp only [hbit, if_false, iff_false]
        rw [symmDiff_eq_iff]
        intro hcon
        rw [symmDiff_self] at hcon
        exact hu0 (by simpa using hcon)
    rw [Finset.filter_congr hcond]
    have htotal :
        ((perp u).filter (fun v => Even (q ∩ v).card)).card
          + ((perp u).filter (fun v => ¬ Even (q ∩ v).card)).card = (perp u).card :=
      Finset.card_filter_add_card_filter_not _
    have hfair := card_perp_fair (u := u) hq (Ne.symm hune)
    omega

/-! ## The law in the form the transfer matrices consume

Every row of `T_j(z)` is an expectation against the transvection law, so what
the matrix entries actually need is

    E[g(M q)] = ½ g(q) + (1/(2M)) ∑_{w ≠ 0} g(w),

an *exact* identity for every `g` -- no sign condition, and the `½` and the
`1/(2M)` are literally the two halves of `½ δ_q + ½ Unif`. -/

lemma card_perp_eq {u u' : Finset (Fin s)} (hu : u ≠ ∅) (hu' : u' ≠ ∅) :
    (perp u).card = (perp u').card := by
  have h1 := card_perp_double hu
  have h2 := card_perp_double hu'
  omega

lemma nonzeroStates_eq_erase :
    nonzeroStates s = (univ : Finset (Finset (Fin s))).erase ∅ := by
  ext w
  simp only [nonzeroStates, Finset.mem_filter, Finset.mem_univ, true_and,
    Finset.mem_erase, and_true]

lemma mem_nonzeroStates {w : Finset (Fin s)} (hw : w ≠ ∅) : w ∈ nonzeroStates s := by
  simp only [nonzeroStates, Finset.mem_filter, Finset.mem_univ, true_and]
  exact hw

lemma ne_empty_of_mem_nonzeroStates {w : Finset (Fin s)} (hw : w ∈ nonzeroStates s) :
    w ≠ ∅ := by
  simpa [nonzeroStates] using hw

lemma card_perp_pos {q : Finset (Fin s)} (hq : q ≠ ∅) : 0 < (perp q).card := by
  have h2 := card_perp_double hq
  have hpow : 0 < 2 ^ s := by positivity
  omega

/-- The sample space has `M · |u^⊥|` points. -/
theorem card_transPairs {q : Finset (Fin s)} (hq : q ≠ ∅) :
    (transPairs s).card = (nonzeroStates s).card * (perp q).card := by
  classical
  have h := card_filter_transPairs (s := s) (fun _ => True)
  simp only [Finset.filter_true] at h
  rw [h, Finset.sum_congr rfl (fun u hu =>
    card_perp_eq (ne_empty_of_mem_nonzeroStates hu) hq),
    Finset.sum_const, smul_eq_mul]

/-- **The transvection expectation**, in exact integer form. -/
theorem sum_act_eq {q : Finset (Fin s)} (hq : q ≠ ∅) (g : Finset (Fin s) → ℝ) :
    2 * ∑ p ∈ transPairs s, g (act p.1 p.2 q)
      = ((perp q).card : ℝ)
        * (((nonzeroStates s).card : ℝ) * g q + ∑ w ∈ nonzeroStates s, g w) := by
  classical
  set cnt : Finset (Fin s) → ℕ :=
    fun w => ((transPairs s).filter (fun p => act p.1 p.2 q = w)).card with hcnt
  have hqmem : q ∈ nonzeroStates s := mem_nonzeroStates hq
  -- group the sample space by target
  have hinner : ∀ w : Finset (Fin s), (cnt w : ℝ) * g w
      = ∑ p ∈ (transPairs s).filter (fun p => act p.1 p.2 q = w),
          g (act p.1 p.2 q) := by
    intro w
    rw [Finset.sum_congr rfl (fun p hp => by rw [(Finset.mem_filter.mp hp).2]),
      Finset.sum_const, nsmul_eq_mul]
  have hfib : ∑ w : Finset (Fin s), (cnt w : ℝ) * g w
      = ∑ p ∈ transPairs s, g (act p.1 p.2 q) := by
    rw [Finset.sum_congr rfl (fun w _ => hinner w)]
    exact Finset.sum_fiberwise_of_maps_to
      (fun p _ => Finset.mem_univ (act p.1 p.2 q))
      (fun p : Finset (Fin s) × Finset (Fin s) => g (act p.1 p.2 q))
  -- the zero target contributes nothing
  have hz : cnt ∅ = 0 := count_target_zero hq
  have hsplit : ∑ w : Finset (Fin s), (cnt w : ℝ) * g w
      = ∑ w ∈ nonzeroStates s, (cnt w : ℝ) * g w := by
    rw [nonzeroStates_eq_erase, ← Finset.add_sum_erase _ (fun w => (cnt w : ℝ) * g w)
      (Finset.mem_univ (∅ : Finset (Fin s))), hz]
    simp
  have hpeel : ∑ w ∈ nonzeroStates s, (cnt w : ℝ) * g w
      = (cnt q : ℝ) * g q + ∑ w ∈ (nonzeroStates s).erase q, (cnt w : ℝ) * g w :=
    (Finset.add_sum_erase _ (fun w => (cnt w : ℝ) * g w) hqmem).symm
  -- the two counting laws, cast to `ℝ`
  have hother : ∀ w ∈ (nonzeroStates s).erase q,
      2 * (cnt w : ℝ) = ((perp q).card : ℝ) := by
    intro w hw
    have hwq : w ≠ q := (Finset.mem_erase.mp hw).1
    have hw0 : w ≠ ∅ :=
      ne_empty_of_mem_nonzeroStates (Finset.mem_erase.mp hw).2
    have hc := count_target_other hq hw0 hwq
    have hperp : (perp (q ∆ w)).card = (perp q).card := by
      refine card_perp_eq ?_ hq
      intro hcon
      exact hwq (symmDiff_eq_bot.mp (by simpa using hcon)).symm
    rw [hperp] at hc
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) hc
  have hself : 2 * (cnt q : ℝ)
      = 2 * ((perp q).card : ℝ)
        + (((nonzeroStates s).erase q).card : ℝ) * ((perp q).card : ℝ) := by
    have hc := count_target_self hq
    have hsum : ∑ u ∈ (nonzeroStates s).erase q, (perp u).card
        = ((nonzeroStates s).erase q).card * (perp q).card := by
      rw [Finset.sum_congr rfl (fun u hu => card_perp_eq
        (ne_empty_of_mem_nonzeroStates (Finset.mem_erase.mp hu).2) hq),
        Finset.sum_const, smul_eq_mul]
    rw [hsum] at hc
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) hc
  have hcard : (((nonzeroStates s).erase q).card : ℝ)
      = ((nonzeroStates s).card : ℝ) - 1 := by
    have h := Finset.card_erase_add_one hqmem
    have : (((nonzeroStates s).erase q).card : ℝ) + 1 = ((nonzeroStates s).card : ℝ) := by
      exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) h
    linarith
  -- assemble
  have hA : 2 * ∑ p ∈ transPairs s, g (act p.1 p.2 q)
      = (2 * (cnt q : ℝ)) * g q
        + ∑ w ∈ (nonzeroStates s).erase q, 2 * ((cnt w : ℝ) * g w) := by
    rw [← hfib, hsplit, hpeel, mul_add, Finset.mul_sum]
    ring
  have hB : ∑ w ∈ (nonzeroStates s).erase q, 2 * ((cnt w : ℝ) * g w)
      = ((perp q).card : ℝ) * ∑ w ∈ (nonzeroStates s).erase q, g w := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun w hw => ?_
    rw [← mul_assoc, hother w hw]
  have hnz : ∑ w ∈ nonzeroStates s, g w
      = g q + ∑ w ∈ (nonzeroStates s).erase q, g w :=
    (Finset.add_sum_erase _ g hqmem).symm
  rw [hA, hB, hself, hnz, hcard]
  ring

/-- **The normalised form**: `E[g(M q)] = ½ g(q) + (1/(2M)) ∑_{w ≠ 0} g(w)`. -/
theorem expect_act {q : Finset (Fin s)} (hq : q ≠ ∅) (g : Finset (Fin s) → ℝ) :
    (∑ p ∈ transPairs s, g (act p.1 p.2 q)) / ((transPairs s).card : ℝ)
      = g q / 2
        + (∑ w ∈ nonzeroStates s, g w) / (2 * ((nonzeroStates s).card : ℝ)) := by
  classical
  have hqmem : q ∈ nonzeroStates s := mem_nonzeroStates hq
  have hMne : ((nonzeroStates s).card : ℝ) ≠ 0 := by
    have : 0 < (nonzeroStates s).card := Finset.card_pos.mpr ⟨q, hqmem⟩
    positivity
  have hPne : ((perp q).card : ℝ) ≠ 0 := by
    have := card_perp_pos hq
    positivity
  have hcard : ((transPairs s).card : ℝ)
      = ((nonzeroStates s).card : ℝ) * ((perp q).card : ℝ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (card_transPairs hq)
  have hmain := sum_act_eq hq g
  have hS : (∑ p ∈ transPairs s, g (act p.1 p.2 q))
      = ((perp q).card : ℝ) * (((nonzeroStates s).card : ℝ) * g q
          + ∑ w ∈ nonzeroStates s, g w) / 2 := by linarith
  rw [hcard, hS]
  field_simp
/-! ## The law as pointwise probabilities

This is the shape the transfer-matrix entries are written in: a `1/2` on the
diagonal and a flat `1/(2M)` on every nonzero target. -/

lemma two_mul_count_other {q w : Finset (Fin s)} (hq : q ≠ ∅) (hw0 : w ≠ ∅)
    (hwq : w ≠ q) :
    2 * ((((transPairs s).filter (fun p => act p.1 p.2 q = w)).card : ℝ))
      = ((perp q).card : ℝ) := by
  have hc := count_target_other hq hw0 hwq
  have hperp : (perp (q ∆ w)).card = (perp q).card := by
    refine card_perp_eq ?_ hq
    intro hcon
    exact hwq (symmDiff_eq_bot.mp (by simpa using hcon)).symm
  rw [hperp] at hc
  exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) hc

lemma two_mul_count_self {q : Finset (Fin s)} (hq : q ≠ ∅) :
    2 * ((((transPairs s).filter (fun p => act p.1 p.2 q = q)).card : ℝ))
      = (((nonzeroStates s).card : ℝ) + 1) * ((perp q).card : ℝ) := by
  have hqmem : q ∈ nonzeroStates s := mem_nonzeroStates hq
  have hc := count_target_self hq
  have hsum : ∑ u ∈ (nonzeroStates s).erase q, (perp u).card
      = ((nonzeroStates s).erase q).card * (perp q).card := by
    rw [Finset.sum_congr rfl (fun u hu => card_perp_eq
      (ne_empty_of_mem_nonzeroStates (Finset.mem_erase.mp hu).2) hq),
      Finset.sum_const, smul_eq_mul]
  rw [hsum] at hc
  have hcR : 2 * ((((transPairs s).filter (fun p => act p.1 p.2 q = q)).card : ℝ))
      = 2 * ((perp q).card : ℝ)
        + (((nonzeroStates s).erase q).card : ℝ) * ((perp q).card : ℝ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) hc
  have herase : (((nonzeroStates s).erase q).card : ℝ)
      = ((nonzeroStates s).card : ℝ) - 1 := by
    have h := Finset.card_erase_add_one hqmem
    have h' : (((nonzeroStates s).erase q).card : ℝ) + 1
        = ((nonzeroStates s).card : ℝ) := by
      exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) h
    linarith
  rw [hcR, herase]
  ring

/-- **The transvection law, pointwise.**  `½` on the diagonal, `1/(2M)` on
every nonzero target, `0` at zero. -/
theorem prob_act {q : Finset (Fin s)} (hq : q ≠ ∅) (w : Finset (Fin s)) :
    ((((transPairs s).filter (fun p => act p.1 p.2 q = w)).card : ℝ))
        / ((transPairs s).card : ℝ)
      = (if w = q then 1 / 2 else 0)
        + (if w = ∅ then 0 else 1 / (2 * ((nonzeroStates s).card : ℝ))) := by
  classical
  have hqmem : q ∈ nonzeroStates s := mem_nonzeroStates hq
  have hNne : ((nonzeroStates s).card : ℝ) ≠ 0 := by
    have : 0 < (nonzeroStates s).card := Finset.card_pos.mpr ⟨q, hqmem⟩
    positivity
  have hPne : ((perp q).card : ℝ) ≠ 0 := by
    have := card_perp_pos hq
    positivity
  have hcard : ((transPairs s).card : ℝ)
      = ((nonzeroStates s).card : ℝ) * ((perp q).card : ℝ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (card_transPairs hq)
  rw [hcard]
  by_cases hw0 : w = ∅
  · subst hw0
    have hz : (((transPairs s).filter (fun p => act p.1 p.2 q = ∅)).card : ℝ) = 0 := by
      exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (count_target_zero hq)
    rw [hz, if_neg (fun h => hq h.symm), if_pos rfl]
    simp
  · by_cases hwq : w = q
    · subst hwq
      have h := two_mul_count_self hq
      rw [if_pos rfl, if_neg hw0]
      field_simp
      linarith
    · have h := two_mul_count_other hq hw0 hwq
      rw [if_neg hwq, if_neg hw0]
      field_simp
      linarith

end Spin
