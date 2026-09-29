/-
Soundness of the fixed-point BA evaluators.

Each lemma has the same shape: *if the evaluator returned a value at all, that
value encloses the real quantity.*  The `Option` is what carries the shift
search's success, so none of these lemmas needs a hypothesis about the search —
which is the point of returning `Option` rather than assuming a good shift.
-/
import SpinCodes.Numeric.BAEvalDefs
import SpinCodes.Numeric.FixedLog
import SpinCodes.Numeric.BAExponent

set_option linter.unusedSectionVars false

namespace Spin.Numeric

namespace Fix

open Real

/-! ## The two-sided shift -/

private lemma log_shift {p q : Int} {k j : Nat} (hp : 0 < p) (hq : 0 < q) :
    Real.log (((p * ipow2 k : Int) : ℝ) / ((q * ipow2 j : Int) : ℝ))
      - (k : ℝ) * Real.log 2 + (j : ℝ) * Real.log 2
      = Real.log ((p : ℝ) / (q : ℝ)) := by
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hk : ((p * ipow2 k : Int) : ℝ) = (p : ℝ) * 2 ^ k := by
    push_cast [ipow2_cast]; ring
  have hj : ((q * ipow2 j : Int) : ℝ) = (q : ℝ) * 2 ^ j := by
    push_cast [ipow2_cast]; ring
  rw [hk, hj, Real.log_div (by positivity) (by positivity),
    Real.log_mul (ne_of_gt hpR) (by positivity),
    Real.log_mul (ne_of_gt hqR) (by positivity),
    Real.log_pow, Real.log_pow, Real.log_div (ne_of_gt hpR) (ne_of_gt hqR)]
  push_cast
  ring

/-- Positivity of a scaled numerator gives positivity of the numerator. -/
private lemma pos_of_scaled {p : Int} {k : Nat} (h : 0 < p * ipow2 k) : 0 < p := by
  rcases lt_trichotomy p 0 with hlt | heq | hgt
  · exact absurd h (not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg hlt.le (ipow2_pos k).le))
  · rw [heq] at h; simp at h
  · exact hgt

/-- **`flogKJ` encloses `Real.log (p/q)` whenever it answers.** -/
theorem flogKJ_mem {p q : Int} {k j n : Nat} {v : Fix} (h : flogKJ p q k j n = some v) :
    Mem v (Real.log ((p : ℝ) / (q : ℝ))) := by
  unfold flogKJ at h
  split_ifs at h with hcond
  · obtain ⟨hP, hQ, hw1, hw2⟩ := hcond
    have hp : 0 < p := pos_of_scaled hP
    have hq : 0 < q := pos_of_scaled hQ
    have hbase := flogA_mem (P := p * ipow2 k) (Q := q * ipow2 j) (m := n)
      hP hQ hw1 hw2
    have hk2 : Mem (mul (ofInt (k : Int)) log2) ((k : ℝ) * Real.log 2) := by
      have hki : Mem (ofInt (k : Int)) ((k : ℝ)) := by
        have := ofInt_mem ((k : Int)); simpa using this
      exact mul_mem hki log2_mem
    have hj2 : Mem (mul (ofInt (j : Int)) log2) ((j : ℝ) * Real.log 2) := by
      have hji : Mem (ofInt (j : Int)) ((j : ℝ)) := by
        have := ofInt_mem ((j : Int)); simpa using this
      exact mul_mem hji log2_mem
    have hcomb := add_mem (sub_mem hbase hk2) hj2
    rw [log_shift hp hq] at hcomb
    rw [← Option.some.inj h]
    exact hcomb

theorem flogQ_mem {p q : Int} {n : Nat} {v : Fix} (h : flogQ p q n = some v) :
    Mem v (Real.log ((p : ℝ) / (q : ℝ))) := flogKJ_mem h

/-- The guard `flogKJ` checks includes positivity of the numerator. -/
lemma flogQ_pos {p q : Int} {n : Nat} {v : Fix} (h : flogQ p q n = some v) : 0 < p := by
  unfold flogQ flogKJ at h
  split_ifs at h with hc
  exact pos_of_scaled hc.1

/-! ## The log oracle -/

/-- What an evaluator needs of its logarithm source: an answer is positive at
the numerator and encloses the logarithm. -/
def Oracle (lg : LogFn) : Prop :=
  ∀ p v, lg p = some v → 0 < p ∧ Mem v (Real.log ((p : ℝ) / ((scale : Int) : ℝ)))

lemma directLog_oracle (n : Nat) : Oracle (directLog n) :=
  fun _ _ h => ⟨flogQ_pos h, flogQ_mem h⟩

/-- A verified node is all a lookup needs: `find` answers only where it has
compared the key equal, so no ordering invariant enters the proof. -/
lemma find_check {n : Nat} : ∀ {t : LogTree} {k : Int} {v : Fix},
    t.check n = true → t.find k = some v → flogQ k scale n = some v
  | .leaf, _, _, _, h => by simp [LogTree.find] at h
  | .node key val l r, k, v, hc, h => by
      simp only [LogTree.check, Bool.and_eq_true, decide_eq_true_eq] at hc
      obtain ⟨⟨hnode, hl⟩, hr⟩ := hc
      simp only [LogTree.find] at h
      split_ifs at h with h1 h2
      · rw [h1, ← Option.some.inj h]
        exact hnode
      · exact find_check hl h
      · exact find_check hr h

lemma treeLog_oracle {t : LogTree} {n : Nat} (hc : t.check n = true) :
    Oracle (treeLog t) := by
  intro p v h
  have hq : flogQ p scale n = some v := find_check hc h
  exact ⟨flogQ_pos hq, flogQ_mem hq⟩

/-! ### A table split across modules -/

lemma treeLogN_find {n : Nat} : ∀ {ts : List LogTree} {p : Int} {v : Fix},
    checkAll n ts = true → treeLogN ts p = some v → flogQ p scale n = some v
  | [], _, _, _, h => by simp [treeLogN] at h
  | t :: ts, p, v, hc, h => by
      simp only [checkAll, Bool.and_eq_true] at hc
      simp only [treeLogN] at h
      rcases hf : t.find p with _ | u
      · rw [hf] at h
        simp only [Option.orElse_none] at h
        exact treeLogN_find hc.2 h
      · rw [hf] at h
        simp only [Option.orElse_some, Option.some.injEq] at h
        rw [← h]
        exact find_check hc.1 hf

lemma treeLogN_oracle {ts : List LogTree} {n : Nat} (hc : checkAll n ts = true) :
    Oracle (treeLogN ts) := by
  intro p v h
  have hq : flogQ p scale n = some v := treeLogN_find hc h
  exact ⟨flogQ_pos hq, flogQ_mem hq⟩

/-! ## The evaluators, over any oracle -/

theorem flogIWL_mem {lg : LogFn} (hlg : Oracle lg) {a : Fix} {x : ℝ} {v : Fix}
    (hx : Mem a x) (h : flogIWL lg a = some v) : Mem v (Real.log x) := by
  simp only [flogIWL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨L, hL, H, hH, rfl⟩ := h
  obtain ⟨halo, hLm⟩ := hlg _ _ hL
  obtain ⟨_, hHm⟩ := hlg _ _ hH
  have hlopos : (0 : ℝ) < (a.lo : ℝ) / ((scale : Int) : ℝ) :=
    div_pos (by exact_mod_cast halo) scaleR_pos
  have hxpos : 0 < x := lt_of_lt_of_le hlopos hx.1
  exact ⟨le_trans hLm.1 (Real.log_le_log hlopos hx.1),
         le_trans (Real.log_le_log hxpos hx.2) hHm.2⟩

/-! ## `x log x` -/

lemma one_mem : Mem one 1 := by
  have := ofInt_mem 1
  simpa [one] using this

lemma xlPointL_nonneg {lg : LogFn} (hlg : Oracle lg) {p : Int} {v : Fix}
    (h : xlPointL lg p = some v) : 0 ≤ p := by
  unfold xlPointL at h
  split_ifs at h with hp
  · exact le_of_eq hp.symm
  · rcases hq : lg p with _ | l
    · rw [hq] at h; simp at h
    · exact (hlg _ _ hq).1.le

theorem xlPointL_mem {lg : LogFn} (hlg : Oracle lg) {p : Int} {v : Fix}
    (h : xlPointL lg p = some v) :
    Mem v (xlogx ((p : ℝ) / ((scale : Int) : ℝ))) := by
  unfold xlPointL at h
  split_ifs at h with hp
  · subst hp
    rw [← Option.some.inj h]
    simpa using zero_mem
  · rcases hq : lg p with _ | l
    · rw [hq] at h; simp at h
    · rw [hq] at h
      simp only [Option.bind_some, Option.some.injEq] at h
      rw [← h]
      exact mul_mem (sc_mem p) (hlg _ _ hq).2

theorem fxlogxL_mem {lg : LogFn} (hlg : Oracle lg) {a : Fix} {x : ℝ} {v : Fix}
    (hx : Mem a x) (h : fxlogxL lg a = some v) : Mem v (xlogx x) := by
  simp only [fxlogxL, Option.bind_eq_some_iff] at h
  obtain ⟨L, hL, H, hH, h⟩ := h
  have hS := scaleR_pos
  have hlo0 : 0 ≤ a.lo := xlPointL_nonneg hlg hL
  have hx0 : 0 ≤ x := le_trans (div_nonneg (by exact_mod_cast hlo0) hS.le) hx.1
  have hLm := xlPointL_mem hlg hL
  have hHm := xlPointL_mem hlg hH
  split_ifs at h with hzero
  · have hxz : x = 0 := by
      have h1 : x ≤ (a.hi : ℝ) / ((scale : Int) : ℝ) := hx.2
      rw [hzero] at h1
      simp at h1
      linarith
    rw [← Option.some.inj h, hxz]
    simpa using zero_mem
  obtain ⟨lm, hlm, hv⟩ := Option.bind_eq_some_iff.mp h
  rw [← Option.some.inj hv]
  refine ⟨?_, ?_⟩
  · have hmpos : 0 < tangentPt a := (hlg _ _ hlm).1
    have hmR : (0 : ℝ) < ((tangentPt a : Int) : ℝ) / ((scale : Int) : ℝ) :=
      div_pos (by exact_mod_cast hmpos) hS
    have haff := sub_mem (mul_mem hx (add_mem (hlg _ _ hlm).2 one_mem))
      (sc_mem (tangentPt a))
    refine le_trans haff.1 ?_
    have ht := xlogx_tangent hx0 hmR
    nlinarith [ht]
  · have hmax := xlogx_le_max (lo := (a.lo : ℝ) / ((scale : Int) : ℝ))
      (hi := (a.hi : ℝ) / ((scale : Int) : ℝ)) (x := x)
      (div_nonneg (by exact_mod_cast hlo0) hS.le) ⟨hx.1, hx.2⟩
    have h1 : (L.hi : ℝ) / ((scale : Int) : ℝ)
        ≤ ((imax L.hi H.hi : Int) : ℝ) / ((scale : Int) : ℝ) :=
      div_le_div_right_of_le (by exact_mod_cast left_le_imax L.hi H.hi) hS
    have h2 : (H.hi : ℝ) / ((scale : Int) : ℝ)
        ≤ ((imax L.hi H.hi : Int) : ℝ) / ((scale : Int) : ℝ) :=
      div_le_div_right_of_le (by exact_mod_cast right_le_imax L.hi H.hi) hS
    exact le_trans hmax (max_le (le_trans hLm.2 h1) (le_trans hHm.2 h2))

/-! ## Entropy -/

theorem fhEntL_mem {lg : LogFn} (hlg : Oracle lg) {x : Fix} {t : ℝ} {v : Fix}
    (hx : Mem x t) (h : fhEntL lg x = some v) : Mem v (hEnt t) := by
  simp only [fhEntL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨t1, hA, t2, hB, rfl⟩ := h
  have h1x : Mem (sub one x) (1 - t) := sub_mem one_mem hx
  rw [hEnt_eq_xlogx]
  exact neg_mem (add_mem (fxlogxL_mem hlg hx hA) (fxlogxL_mem hlg h1x hB))

/-! ## The accumulator exponent -/

theorem fpiEvalL_mem {lg : LogFn} (hlg : Oracle lg) {a c : Fix} {al ga : ℝ} {v : Fix}
    (ha : Mem a al) (hc : Mem c ga) (h : fpiEvalL lg a c = some v) :
    Mem v (piEval al ga) := by
  simp only [fpiEvalL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨t1, h1, t2, h2, t3, h3, t4, h4, t5, h5, rfl⟩ := h
  have hah : Mem (divInt a 2) (al / 2) := by
    have := divInt_mem ha (k := 2) (by norm_num)
    simpa using this
  have hmc : Mem (sub one c) (1 - ga) := sub_mem one_mem hc
  have hmu : Mem (sub c (divInt a 2)) (ga - al / 2) := sub_mem hc hah
  have hmv : Mem (sub (sub one c) (divInt a 2)) (1 - ga - al / 2) := by
    have hh := sub_mem hmc hah
    rw [show (1 - ga) - al / 2 = 1 - ga - al / 2 by ring] at hh
    exact hh
  have hma : Mem (sub one a) (1 - al) := sub_mem one_mem ha
  rw [piEval_eq_xlogx]
  exact add_mem (sub_mem (sub_mem (add_mem (add_mem (mul_mem ha log2_mem)
    (fxlogxL_mem hlg hc h1)) (fxlogxL_mem hlg hmc h2)) (fxlogxL_mem hlg hmu h3))
    (fxlogxL_mem hlg hmv h4)) (fxlogxL_mem hlg hma h5)

/-! ## The centered bound -/

lemma le_of_mem {a : Fix} {x : ℝ} (h : Mem a x) : a.lo ≤ a.hi := by
  have := le_trans h.1 h.2
  have := (div_le_div_iff_of_pos_right scaleR_pos).mp this
  exact_mod_cast this

/-- The midpoint, as a value of the enclosing box. -/
lemma midPt_val_mem {A : Fix} {x : ℝ} (hx : Mem A x) :
    Mem A (((fdiv (A.lo + A.hi) 2 : Int) : ℝ) / ((scale : Int) : ℝ)) := by
  have hle : A.lo ≤ A.hi := le_of_mem hx
  have hb : A.lo ≤ fdiv (A.lo + A.hi) 2 ∧ fdiv (A.lo + A.hi) 2 ≤ A.hi := by
    unfold fdiv; omega
  exact ⟨div_le_div_right_of_le (by exact_mod_cast hb.1) scaleR_pos,
         div_le_div_right_of_le (by exact_mod_cast hb.2) scaleR_pos⟩

theorem fDpaL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {x y : ℝ} {v : Fix}
    (hx : Mem A x) (hy : Mem C y) (h : fDpaL lg A C = some v) : Mem v (Dpa x y) := by
  simp only [fDpaL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨l1, h1, l2, h2, l3, h3, rfl⟩ := h
  have hah : Mem (divInt A 2) (x / 2) := by
    have := divInt_mem hx (k := 2) (by norm_num); simpa using this
  have m1 : Mem (sub C (divInt A 2)) (y - x / 2) := sub_mem hy hah
  have m2 : Mem (sub (sub one C) (divInt A 2)) (1 - y - x / 2) := by
    have hh := sub_mem (sub_mem one_mem hy) hah
    rw [show (1 - y) - x / 2 = 1 - y - x / 2 by ring] at hh
    exact hh
  have m3 : Mem (sub one A) (1 - x) := sub_mem one_mem hx
  have e1 := flogIWL_mem hlg m1 h1
  have e2 := flogIWL_mem hlg m2 h2
  have e3 := flogIWL_mem hlg m3 h3
  have hhalf : Mem (divInt (add l1 l2) 2)
      ((Real.log (y - x / 2) + Real.log (1 - y - x / 2)) / 2) := by
    have := divInt_mem (add_mem e1 e2) (k := 2) (by norm_num)
    simpa using this
  unfold Dpa
  exact sub_mem (add_mem log2_mem hhalf) e3

theorem fDpcL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {x y : ℝ} {v : Fix}
    (hx : Mem A x) (hy : Mem C y) (h : fDpcL lg A C = some v) : Mem v (Dpc x y) := by
  simp only [fDpcL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨l1, h1, l2, h2, l3, h3, l4, h4, rfl⟩ := h
  have hah : Mem (divInt A 2) (x / 2) := by
    have := divInt_mem hx (k := 2) (by norm_num); simpa using this
  have m2 : Mem (sub one C) (1 - y) := sub_mem one_mem hy
  have m3 : Mem (sub C (divInt A 2)) (y - x / 2) := sub_mem hy hah
  have m4 : Mem (sub (sub one C) (divInt A 2)) (1 - y - x / 2) := by
    have hh := sub_mem m2 hah
    rw [show (1 - y) - x / 2 = 1 - y - x / 2 by ring] at hh
    exact hh
  unfold Dpc
  exact add_mem (sub_mem (sub_mem (flogIWL_mem hlg hy h1) (flogIWL_mem hlg m2 h2))
    (flogIWL_mem hlg m3 h3)) (flogIWL_mem hlg m4 h4)

/-- **The centered bound encloses `π` on the box.**

The five integer hypotheses are the verifier's feasibility guard, scaled: the
box must sit strictly inside the region where all five logarithms have
positive arguments. -/
theorem fpiCenteredL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c)
    (f1 : 0 < 2 * C.lo - A.hi) (f2 : 0 < 2 * scale - 2 * C.hi - A.hi)
    (f3 : A.hi < scale) (f4 : 0 < C.lo) (f5 : C.hi < scale)
    (h : fpiCenteredL lg A C = some v) : Mem v (piEval a c) := by
  have hS := scaleR_pos
  simp only [fpiCenteredL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨p0, hp0, da, hda, dc, hdc, rfl⟩ := h
  -- the real box
  set alo : ℝ := (A.lo : ℝ) / ((scale : Int) : ℝ) with halo
  set ahi : ℝ := (A.hi : ℝ) / ((scale : Int) : ℝ) with hahi
  set clo : ℝ := (C.lo : ℝ) / ((scale : Int) : ℝ) with hclo
  set chi : ℝ := (C.hi : ℝ) / ((scale : Int) : ℝ) with hchi
  have c1 : 0 < clo - ahi / 2 := by
    have hpos : (0 : ℝ) < 2 * (C.lo : ℝ) - (A.hi : ℝ) := by
      have hx : (0 : ℝ) < ((2 * C.lo - A.hi : Int) : ℝ) := by exact_mod_cast f1
      push_cast at hx; linarith
    have heq : clo - ahi / 2
        = (2 * (C.lo : ℝ) - (A.hi : ℝ)) / (2 * ((scale : Int) : ℝ)) := by
      rw [hclo, hahi]; field_simp; try ring
    rw [heq]
    exact div_pos hpos (by positivity)
  have c2 : 0 < 1 - chi - ahi / 2 := by
    have hpos : (0 : ℝ) < 2 * ((scale : Int) : ℝ) - 2 * (C.hi : ℝ) - (A.hi : ℝ) := by
      have hx : (0 : ℝ) < ((2 * scale - 2 * C.hi - A.hi : Int) : ℝ) := by exact_mod_cast f2
      push_cast at hx; linarith
    have heq : 1 - chi - ahi / 2
        = (2 * ((scale : Int) : ℝ) - 2 * (C.hi : ℝ) - (A.hi : ℝ))
          / (2 * ((scale : Int) : ℝ)) := by
      rw [hchi, hahi]; field_simp; try ring
    rw [heq]
    exact div_pos hpos (by positivity)
  have c3 : ahi < 1 := by
    rw [hahi, div_lt_one hS]; exact_mod_cast f3
  have c4 : 0 < clo := by
    rw [hclo]; exact div_pos (by exact_mod_cast f4) hS
  have c5 : chi < 1 := by
    rw [hchi, div_lt_one hS]; exact_mod_cast f5
  -- the midpoints
  set am : ℝ := ((fdiv (A.lo + A.hi) 2 : Int) : ℝ) / ((scale : Int) : ℝ) with ham0
  set cm : ℝ := ((fdiv (C.lo + C.hi) 2 : Int) : ℝ) / ((scale : Int) : ℝ) with hcm0
  have hamA : Mem A am := midPt_val_mem ha
  have hcmC : Mem C cm := midPt_val_mem hc
  obtain ⟨ξ, hξ, η, hη, hid⟩ :=
    piEval_centered (alo := alo) (ahi := ahi) (clo := clo) (chi := chi)
      c1 c2 c3 c4 c5 (a := a) (c := c) (am := am) (cm := cm)
      ⟨ha.1, ha.2⟩ ⟨hc.1, hc.2⟩ ⟨hamA.1, hamA.2⟩ ⟨hcmC.1, hcmC.2⟩
  have hmA : Mem (midPt A) am := sc_mem _
  have hmC : Mem (midPt C) cm := sc_mem _
  have e0 := fpiEvalL_mem hlg hmA hmC hp0
  have ea := fDpaL_mem hlg (A := A) (C := C) (x := ξ) (y := c) ⟨hξ.1, hξ.2⟩ hc hda
  have ec := fDpcL_mem hlg (A := A) (C := C) (x := am) (y := η) hamA ⟨hη.1, hη.2⟩ hdc
  have da1 : Mem (sub A (midPt A)) (a - am) := sub_mem ha hmA
  have dc1 : Mem (sub C (midPt C)) (c - cm) := sub_mem hc hmC
  rw [hid]
  exact add_mem (add_mem e0 (mul_mem ea da1)) (mul_mem ec dc1)

/-! ## Intersecting enclosures, and the guarded centered bound -/

lemma meet_mem {a b : Fix} {x : ℝ} (ha : Mem a x) (hb : Mem b x) : Mem (meet a b) x := by
  have hS := scaleR_pos
  refine ⟨?_, ?_⟩
  · show ((imax a.lo b.lo : Int) : ℝ) / ((scale : Int) : ℝ) ≤ x
    unfold imax
    split_ifs with h
    · exact hb.1
    · exact ha.1
  · show x ≤ ((imin a.hi b.hi : Int) : ℝ) / ((scale : Int) : ℝ)
    unfold imin
    split_ifs with h
    · exact ha.2
    · exact hb.2

theorem fpiCenteredG_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c) (h : fpiCenteredG lg A C = some v) :
    Mem v (piEval a c) := by
  unfold fpiCenteredG at h
  split_ifs at h with hf
  · simp only [feasOK, Bool.and_eq_true, decide_eq_true_eq] at hf
    obtain ⟨⟨⟨⟨f1, f2⟩, f3⟩, f4⟩, f5⟩ := hf
    exact fpiCenteredL_mem hlg ha hc f1 f2 f3 f4 f5 h

theorem fpiBestL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ} {v : Fix}
    (ha : Mem A a) (hc : Mem C c) (h : fpiBestL lg A C = some v) :
    Mem v (piEval a c) := by
  unfold fpiBestL at h
  rcases hp : fpiEvalL lg A C with _ | p <;> rcases hq : fpiCenteredG lg A C with _ | q <;>
    rw [hp, hq] at h <;> simp only [Option.some.injEq] at h
  · exact absurd h (by simp)
  · rw [← h]; exact fpiCenteredG_mem hlg ha hc hq
  · rw [← h]; exact fpiEvalL_mem hlg ha hc hp
  · rw [← h]
    exact meet_mem (fpiEvalL_mem hlg ha hc hp) (fpiCenteredG_mem hlg ha hc hq)

/-! ## The split tree covers the box

The whole covering argument.  At a split the two children agree with the
parent except in one coordinate, where one takes `[lo, m]` and the other
`[m, hi]`; any real in `[lo, hi]` is in one of those two.  No property of the
leaves enters, so `leafOK` and `P` stay parameters. -/

lemma mem_setHi_or_setLo {X : Fix} {x : ℝ} (hx : Mem X x) (m : Int) :
    Mem (setHi X m) x ∨ Mem (setLo X m) x := by
  rcases le_or_gt x ((m : ℝ) / ((scale : Int) : ℝ)) with h | h
  · exact Or.inl ⟨hx.1, h⟩
  · exact Or.inr ⟨h.le, hx.2⟩

theorem checkTree_sound {leafOK : Int → Fix → Fix → Fix → Bool}
    {P : ℝ → ℝ → ℝ → Prop}
    (hleaf : ∀ u A B W a b w, leafOK u A B W = true →
      Mem A a → Mem B b → Mem W w → P a b w) :
    ∀ (t : BoxTree) (A B W : Fix) (a b w : ℝ),
      checkTree leafOK t A B W = true →
      Mem A a → Mem B b → Mem W w → P a b w
  | .leaf u, A, B, W, a, b, w, h, ha, hb, hw => hleaf u A B W a b w h ha hb hw
  | .split 0 m l r, A, B, W, a, b, w, h, ha, hb, hw => by
      simp only [checkTree, Bool.and_eq_true] at h
      rcases mem_setHi_or_setLo ha m with hc | hc
      · exact checkTree_sound hleaf l _ B W a b w h.1 hc hb hw
      · exact checkTree_sound hleaf r _ B W a b w h.2 hc hb hw
  | .split 1 m l r, A, B, W, a, b, w, h, ha, hb, hw => by
      simp only [checkTree, Bool.and_eq_true] at h
      rcases mem_setHi_or_setLo hb m with hc | hc
      · exact checkTree_sound hleaf l A _ W a b w h.1 ha hc hw
      · exact checkTree_sound hleaf r A _ W a b w h.2 ha hc hw
  | .split (n + 2) m l r, A, B, W, a, b, w, h, ha, hb, hw => by
      simp only [checkTree, Bool.and_eq_true] at h
      rcases mem_setHi_or_setLo hw m with hc | hc
      · exact checkTree_sound hleaf l A B _ a b w h.1 ha hb hc
      · exact checkTree_sound hleaf r A B _ a b w h.2 ha hb hc

/-! ## Splitting a region at the claim level

`checkTree_sound` composes *checks*, which forces every subtree to share one
oracle — and therefore one table, loaded into every module that proves a
subtree.  Splitting at the level of the conclusion instead lets each module
carry its own table and prove its own piece independently; only the finished
claims are combined.

The argument is the same one `checkTree_sound` makes at each split, lifted out
so it can be used across module boundaries. -/

lemma claim_split_fst {P : ℝ → ℝ → ℝ → Prop} {A B W : Fix} {m : Int}
    (hl : ∀ a b w, Mem (setHi A m) a → Mem B b → Mem W w → P a b w)
    (hr : ∀ a b w, Mem (setLo A m) a → Mem B b → Mem W w → P a b w)
    (a b w : ℝ) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) : P a b w := by
  rcases mem_setHi_or_setLo ha m with h | h
  · exact hl a b w h hb hw
  · exact hr a b w h hb hw

lemma claim_split_snd {P : ℝ → ℝ → ℝ → Prop} {A B W : Fix} {m : Int}
    (hl : ∀ a b w, Mem A a → Mem (setHi B m) b → Mem W w → P a b w)
    (hr : ∀ a b w, Mem A a → Mem (setLo B m) b → Mem W w → P a b w)
    (a b w : ℝ) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) : P a b w := by
  rcases mem_setHi_or_setLo hb m with h | h
  · exact hl a b w ha h hw
  · exact hr a b w ha h hw

lemma claim_split_thd {P : ℝ → ℝ → ℝ → Prop} {A B W : Fix} {m : Int}
    (hl : ∀ a b w, Mem A a → Mem B b → Mem (setHi W m) w → P a b w)
    (hr : ∀ a b w, Mem A a → Mem B b → Mem (setLo W m) w → P a b w)
    (a b w : ℝ) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) : P a b w := by
  rcases mem_setHi_or_setLo hw m with h | h
  · exact hl a b w ha hb h
  · exact hr a b w ha hb h

/-! ## The Golay exponent -/

theorem fGolayG_mem {u : Fix} {w : ℝ} (hu : Mem u w) : Mem (fGolayG u) (GolayG w) := by
  have c759 : Mem (ofInt 759) (759 : ℝ) := by simpa using ofInt_mem 759
  have c2576 : Mem (ofInt 2576) (2576 : ℝ) := by simpa using ofInt_mem 2576
  have hu2 := mul_mem hu hu
  have hv := mul_mem hu2 hu2
  have hv2 := mul_mem hv hv
  have hres := add_mem one_mem (mul_mem hv2 (add_mem c759
    (mul_mem hv (add_mem c2576 (mul_mem hv (add_mem c759 hv2))))))
  unfold fGolayG GolayG
  convert hres using 1
  ring

theorem fgObjL_mem {lg : LogFn} (hlg : Oracle lg) {u a : Fix} {w al : ℝ} {v : Fix}
    (hu : Mem u w) (ha : Mem a al) (h : fgObjL lg u a = some v) :
    Mem v (gObj w al) := by
  simp only [fgObjL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨lgG, h1, lu, h2, rfl⟩ := h
  have e1 := flogIWL_mem hlg (fGolayG_mem hu) h1
  have e2 := flogIWL_mem hlg hu h2
  have hdiv : Mem (divInt lgG 24) (Real.log (GolayG w) / 24) := by
    have := divInt_mem e1 (k := 24) (by norm_num)
    simpa using this
  exact sub_mem hdiv (mul_mem ha e2)

/-! ## Clamping to the feasible part -/

lemma clampLo_mem {X : Fix} {x : ℝ} (hx : Mem X x) (h0 : 0 ≤ x) :
    Mem (clampLo X) x := by
  refine ⟨?_, hx.2⟩
  show ((imax 0 X.lo : Int) : ℝ) / ((scale : Int) : ℝ) ≤ x
  unfold imax
  split_ifs with h
  · exact hx.1
  · simpa using h0

/-- **The clamped plain bound.**  Sound wherever the point is feasible, which
is the only place the claim asserts anything. -/
theorem fpiEvalClampL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c)
    (hf1 : a / 2 < c) (hf2 : c < 1 - a / 2)
    (h : fpiEvalClampL lg A C = some v) : Mem v (piEval a c) := by
  simp only [fpiEvalClampL, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨t1, h1, t2, h2, t3, h3, t4, h4, t5, h5, rfl⟩ := h
  have hah : Mem (divInt A 2) (a / 2) := by
    have := divInt_mem ha (k := 2) (by norm_num)
    simpa using this
  have hmc : Mem (sub one C) (1 - c) := sub_mem one_mem hc
  have hmu : Mem (clampLo (sub C (divInt A 2))) (c - a / 2) :=
    clampLo_mem (sub_mem hc hah) (by linarith)
  have hmv : Mem (clampLo (sub (sub one C) (divInt A 2))) (1 - c - a / 2) := by
    have hh := sub_mem hmc hah
    rw [show (1 - c) - a / 2 = 1 - c - a / 2 by ring] at hh
    exact clampLo_mem hh (by linarith)
  have hma : Mem (sub one A) (1 - a) := sub_mem one_mem ha
  rw [piEval_eq_xlogx]
  exact add_mem (sub_mem (sub_mem (add_mem (add_mem (mul_mem ha log2_mem)
    (fxlogxL_mem hlg hc h1)) (fxlogxL_mem hlg hmc h2)) (fxlogxL_mem hlg hmu h3))
    (fxlogxL_mem hlg hmv h4)) (fxlogxL_mem hlg hma h5)

theorem fpiBestClampL_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c)
    (hf1 : a / 2 < c) (hf2 : c < 1 - a / 2)
    (h : fpiBestClampL lg A C = some v) : Mem v (piEval a c) := by
  unfold fpiBestClampL at h
  rcases hp : fpiEvalClampL lg A C with _ | p <;> rcases hq : fpiCenteredG lg A C with _ | q <;>
    rw [hp, hq] at h <;> simp only [Option.some.injEq] at h
  · exact absurd h (by simp)
  · rw [← h]; exact fpiCenteredG_mem hlg ha hc hq
  · rw [← h]; exact fpiEvalClampL_mem hlg ha hc hf1 hf2 hp
  · rw [← h]
    exact meet_mem (fpiEvalClampL_mem hlg ha hc hf1 hf2 hp) (fpiCenteredG_mem hlg ha hc hq)

/-! ## The dense-tail claim

What the cover proves, one point at a time.  The hypothesis is the paper's
feasibility region (strictly, which is where `piBA` and `piEval` agree), and
the conclusion is the strict negative gap. -/

/-- The claim of `eq:structured-ba-dense-tail` at a single point. -/
def DenseClaim (thr : Int) (a b w : ℝ) : Prop :=
  (a / 2 < b ∧ b < 1 - a / 2 ∧ b / 2 < w ∧ w < 1 - b / 2 ∧ 0 ≤ a ∧ a < 1) →
    gBA a + piBA a b + piBA b w < (thr : ℝ) / ((scale : Int) : ℝ)

/-- An infeasible box has nothing to prove. -/
lemma infeasible_vacuous {A B W : Fix} {a b w : ℝ} (ha : Mem A a) (hb : Mem B b)
    (hw : Mem W w) (h : infeasible A B W = true) (thr : Int) :
    DenseClaim thr a b w := by
  intro hfeas
  obtain ⟨f1, f2, f3, f4, _, _⟩ := hfeas
  exfalso
  have hS := scaleR_pos
  simp only [infeasible, Bool.or_eq_true, decide_eq_true_eq] at h
  rcases h with ((h | h) | h) | h
  · -- every `b` in the box is below `a/2`
    have h1 : (b : ℝ) ≤ (B.hi : ℝ) / ((scale : Int) : ℝ) := hb.2
    have h2 : (A.lo : ℝ) / ((scale : Int) : ℝ) ≤ a := ha.1
    have hR : (2 : ℝ) * (B.hi : ℝ) < (A.lo : ℝ) := by exact_mod_cast h
    have : b < a / 2 := by
      have : (2 : ℝ) * b ≤ 2 * ((B.hi : ℝ) / ((scale : Int) : ℝ)) := by linarith
      rw [div_le_iff₀ hS] at h2
      rw [le_div_iff₀ hS] at h1
      nlinarith [h1, h2, hR, hS]
    linarith
  · have h1 : (B.lo : ℝ) / ((scale : Int) : ℝ) ≤ b := hb.1
    have h2 : (A.lo : ℝ) / ((scale : Int) : ℝ) ≤ a := ha.1
    have hR : (2 : ℝ) * ((scale : Int) : ℝ) - (A.lo : ℝ) < 2 * (B.lo : ℝ) := by
      exact_mod_cast h
    have : 1 - a / 2 < b := by
      rw [div_le_iff₀ hS] at h1 h2
      nlinarith [h1, h2, hR, hS]
    linarith
  · have h1 : (w : ℝ) ≤ (W.hi : ℝ) / ((scale : Int) : ℝ) := hw.2
    have h2 : (B.lo : ℝ) / ((scale : Int) : ℝ) ≤ b := hb.1
    have hR : (2 : ℝ) * (W.hi : ℝ) < (B.lo : ℝ) := by exact_mod_cast h
    have : w < b / 2 := by
      rw [le_div_iff₀ hS] at h1
      rw [div_le_iff₀ hS] at h2
      nlinarith [h1, h2, hR, hS]
    linarith
  · have h1 : (W.lo : ℝ) / ((scale : Int) : ℝ) ≤ w := hw.1
    have h2 : (B.lo : ℝ) / ((scale : Int) : ℝ) ≤ b := hb.1
    have hR : (2 : ℝ) * ((scale : Int) : ℝ) - (B.lo : ℝ) < 2 * (W.lo : ℝ) := by
      exact_mod_cast h
    have : 1 - b / 2 < w := by
      rw [div_le_iff₀ hS] at h1 h2
      nlinarith [h1, h2, hR, hS]
    linarith

/-- **A leaf that passes proves the claim on its box.** -/
theorem leafOK_sound {lg : LogFn} (hlg : Oracle lg) {thr u : Int} {A B W : Fix}
    {a b w : ℝ} (ha : Mem A a) (hb : Mem B b) (hw : Mem W w)
    (h : leafOK lg thr u A B W = true) : DenseClaim thr a b w := by
  simp only [leafOK, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  rcases h with hinf | ⟨hu, hbody⟩
  · exact infeasible_vacuous ha hb hw hinf thr
  intro hfeas
  obtain ⟨f1, f2, f3, f4, ha0, ha1⟩ := hfeas
  have hS := scaleR_pos
  have hb0 : 0 < b := by linarith
  have hb1 : b < 1 := by linarith
  -- the three enclosures
  rcases hg : fgObjL lg (sc u) A with _ | g
  · rw [hg] at hbody; simp at hbody
  rcases hp1 : fpiBestClampL lg A B with _ | p1
  · rw [hg, hp1] at hbody; simp at hbody
  rcases hp2 : fpiBestClampL lg B W with _ | p2
  · rw [hg, hp1, hp2] at hbody; simp at hbody
  rw [hg, hp1, hp2] at hbody
  simp only [decide_eq_true_eq] at hbody
  have huR : (0 : ℝ) < (u : ℝ) / ((scale : Int) : ℝ) :=
    div_pos (by exact_mod_cast hu) hS
  have eg := fgObjL_mem hlg (sc_mem u) ha hg
  have e1 := fpiBestClampL_mem hlg ha hb f1 f2 hp1
  have e2 := fpiBestClampL_mem hlg hb hw f3 f4 hp2
  have esum := add_mem (add_mem eg e1) e2
  -- the sum is below the threshold
  have hlt : gObj ((u : ℝ) / ((scale : Int) : ℝ)) a + piEval a b + piEval b w
      < (thr : ℝ) / ((scale : Int) : ℝ) := by
    refine lt_of_le_of_lt esum.2 ?_
    have hnum : ((add (add g p1) p2).hi : ℝ) < (thr : ℝ) := by exact_mod_cast hbody
    gcongr
  -- and `g` and `π` are dominated by what was evaluated
  have hgle : gBA a ≤ gObj ((u : ℝ) / ((scale : Int) : ℝ)) a := gBA_le ha0 ha1.le huR
  have hpi1 : piBA a b = piEval a b := piBA_eq_piEval ha0 ha1 f1 f2
  have hpi2 : piBA b w = piEval b w := piBA_eq_piEval hb0.le hb1 f3 f4
  rw [hpi1, hpi2]
  linarith

/-- **The bridge.**  A tree that checks out proves the dense-tail claim at
every point of the root box.

Everything specific to this certificate is in `leafOK`; the covering is
`checkTree_sound`, which knows nothing about `π` or `g`. -/
theorem denseTail_of_checkTree {lg : LogFn} (hlg : Oracle lg) {thr : Int}
    (t : BoxTree) {A B W : Fix} (h : checkTree (leafOK lg thr) t A B W = true)
    {a b w : ℝ} (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) :
    DenseClaim thr a b w :=
  checkTree_sound (P := DenseClaim thr)
    (fun _ _ _ _ _ _ _ hl ha' hb' hw' => leafOK_sound hlg ha' hb' hw' hl)
    t A B W a b w h ha hb hw

/-! ## The direct instances

These are the shapes the pins name; each is the oracle version at
`directLog n`. -/

theorem flogIW_mem {a : Fix} {x : ℝ} {n : Nat} {v : Fix}
    (hx : Mem a x) (h : flogIW a n = some v) : Mem v (Real.log x) :=
  flogIWL_mem (directLog_oracle n) hx h

theorem xlPoint_mem {p : Int} {n : Nat} {v : Fix} (h : xlPoint p n = some v) :
    Mem v (xlogx ((p : ℝ) / ((scale : Int) : ℝ))) :=
  xlPointL_mem (directLog_oracle n) h

theorem fxlogx_mem {a : Fix} {x : ℝ} {n : Nat} {v : Fix}
    (hx : Mem a x) (h : fxlogx a n = some v) : Mem v (xlogx x) :=
  fxlogxL_mem (directLog_oracle n) hx h

theorem fhEnt_mem {x : Fix} {t : ℝ} {n : Nat} {v : Fix}
    (hx : Mem x t) (h : fhEnt x n = some v) : Mem v (hEnt t) :=
  fhEntL_mem (directLog_oracle n) hx h

theorem fpiEval_mem {a c : Fix} {al ga : ℝ} {n : Nat} {v : Fix}
    (ha : Mem a al) (hc : Mem c ga) (h : fpiEval a c n = some v) :
    Mem v (piEval al ga) :=
  fpiEvalL_mem (directLog_oracle n) ha hc h

theorem fgObj_mem {u a : Fix} {w al : ℝ} {n : Nat} {v : Fix}
    (hu : Mem u w) (ha : Mem a al) (h : fgObj u a n = some v) :
    Mem v (gObj w al) :=
  fgObjL_mem (directLog_oracle n) hu ha h

theorem fpiCentered_mem {A C : Fix} {a c : ℝ} {n : Nat} {v : Fix}
    (ha : Mem A a) (hc : Mem C c)
    (f1 : 0 < 2 * C.lo - A.hi) (f2 : 0 < 2 * scale - 2 * C.hi - A.hi)
    (f3 : A.hi < scale) (f4 : 0 < C.lo) (f5 : C.hi < scale)
    (h : fpiCentered A C n = some v) : Mem v (piEval a c) :=
  fpiCenteredL_mem (directLog_oracle n) ha hc f1 f2 f3 f4 f5 h

end Fix

end Spin.Numeric
