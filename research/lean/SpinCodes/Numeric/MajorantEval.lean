import SpinCodes.Numeric.BAEval
import SpinCodes.Numeric.MajorantDefs
import SpinCodes.Numeric.BAClosed

set_option linter.unusedSectionVars false
namespace Spin.Numeric.Fix

noncomputable def midVal (A : Fix) : ℝ :=
  ((fdiv (A.lo + A.hi) 2 : Int) : ℝ) / ((scale : Int) : ℝ)

lemma midVal_mem {A : Fix} {a : ℝ} (ha : Mem A a) : Mem A (midVal A) :=
  midPt_val_mem ha

lemma midVal_point (A : Fix) : Mem (midPt A) (midVal A) := sc_mem _

lemma pi_centered_of_guard {A C : Fix} {a c : ℝ}
    (ha : Mem A a) (hc : Mem C c) (hf : feasOK A C = true) :
    ∃ ξ, Mem A ξ ∧ ∃ η, Mem C η ∧
      piEval a c = piEval (midVal A) (midVal C) +
        Dpa ξ c * (a - midVal A) + Dpc (midVal A) η * (c - midVal C) := by
  simp only [feasOK, Bool.and_eq_true, decide_eq_true_eq] at hf
  obtain ⟨⟨⟨⟨f1, f2⟩, f3⟩, f4⟩, f5⟩ := hf
  have hS := scaleR_pos
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
  exact ⟨ξ, ⟨hξ.1, hξ.2⟩, η, ⟨hη.1, hη.2⟩, hid⟩

theorem fResidualCentered_mem {lg : LogFn} (hlg : Oracle lg)
    {U A B W S I v : Fix} {u a b w s i : ℝ}
    (hu : Mem U u) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w)
    (hs : Mem S s) (hi : Mem I i)
    (h : fResidualCentered lg U A B W S I = some v) :
    Mem v (gObj u a + piEval a b + piEval b w - (s * w + i)) := by
  unfold fResidualCentered at h
  split_ifs at h with hf
  simp only [Bool.and_eq_true] at hf
  simp only [Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨g, hg, p, hp, q, hq, lu, hlu, pa, hpa, pb, hpb,
    qb, hqb, qw, hqw, rfl⟩ := h
  obtain ⟨ξ, hξ, η, hη, ep⟩ := pi_centered_of_guard ha hb hf.1
  obtain ⟨ζ, hζ, θ, hθ, eq⟩ := pi_centered_of_guard hb hw hf.2
  have ma := midVal_point A
  have mb := midVal_point B
  have mw := midVal_point W
  have eg := fgObjL_mem hlg hu ma hg
  have e0 := sub_mem (add_mem (add_mem eg (fpiEvalL_mem hlg ma mb hp))
    (fpiEvalL_mem hlg mb mw hq)) (add_mem (mul_mem hs mw) hi)
  have ea := sub_mem (fDpaL_mem hlg hξ hb hpa) (flogIWL_mem hlg hu hlu)
  have eb := add_mem (fDpcL_mem hlg (midVal_mem ha) hη hpb)
    (fDpaL_mem hlg hζ hw hqb)
  have ew := sub_mem (fDpcL_mem hlg (midVal_mem hb) hθ hqw) hs
  have result := add_mem (add_mem (add_mem e0 (mul_mem ea (sub_mem ha ma)))
    (mul_mem eb (sub_mem hb mb))) (mul_mem ew (sub_mem hw mw))
  convert result using 1
  rw [ep, eq]
  unfold gObj
  ring

theorem fpiEvalClampL_closed_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c)
    (hf1 : a / 2 ≤ c) (hf2 : c ≤ 1 - a / 2)
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

theorem fpiBestClampL_closed_mem {lg : LogFn} (hlg : Oracle lg) {A C : Fix} {a c : ℝ}
    {v : Fix} (ha : Mem A a) (hc : Mem C c)
    (hf1 : a / 2 ≤ c) (hf2 : c ≤ 1 - a / 2)
    (h : fpiBestClampL lg A C = some v) : Mem v (piEval a c) := by
  unfold fpiBestClampL at h
  rcases hp : fpiEvalClampL lg A C with _ | p <;> rcases hq : fpiCenteredG lg A C with _ | q <;>
    rw [hp, hq] at h <;> simp only [Option.some.injEq] at h
  · exact absurd h (by simp)
  · rw [← h]; exact fpiCenteredG_mem hlg ha hc hq
  · rw [← h]; exact fpiEvalClampL_closed_mem hlg ha hc hf1 hf2 hp
  · rw [← h]
    exact meet_mem (fpiEvalClampL_closed_mem hlg ha hc hf1 hf2 hp) (fpiCenteredG_mem hlg ha hc hq)

theorem fResidual_mem {lg : LogFn} (hlg : Oracle lg)
    {U A B W S I v : Fix} {u a b w s i : ℝ}
    (hu : Mem U u) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w)
    (hs : Mem S s) (hi : Mem I i)
    (f1 : a / 2 ≤ b) (f2 : b ≤ 1 - a / 2)
    (f3 : b / 2 ≤ w) (f4 : w ≤ 1 - b / 2)
    (h : fResidual lg U A B W S I = some v) :
    Mem v (gObj u a + piEval a b + piEval b w - (s * w + i)) := by
  simp only [fResidual, Option.bind_eq_some_iff, Option.some.injEq] at h
  obtain ⟨g, hg, p, hp, q, hq, rfl⟩ := h
  exact sub_mem (add_mem (add_mem (fgObjL_mem hlg hu ha hg)
    (fpiBestClampL_closed_mem hlg ha hb f1 f2 hp))
    (fpiBestClampL_closed_mem hlg hb hw f3 f4 hq)) (add_mem (mul_mem hs hw) hi)

/-- The affine bound on the full closed accumulator feasibility region. -/
def MajorantClaim (s i a b w : ℝ) : Prop :=
  (a / 2 ≤ b ∧ b ≤ 1 - a / 2 ∧ b / 2 ≤ w ∧ w ≤ 1 - b / 2 ∧ 0 ≤ a ∧ a ≤ 1) →
    gBA a + piBA a b + piBA b w ≤ s * w + i

lemma majorant_infeasible {A B W : Fix} {a b w : ℝ} (ha : Mem A a) (hb : Mem B b)
    (hw : Mem W w) (h : infeasible A B W = true)  :
    ¬ (a / 2 ≤ b ∧ b ≤ 1 - a / 2 ∧ b / 2 ≤ w ∧ w ≤ 1 - b / 2 ∧ 0 ≤ a ∧ a ≤ 1) := by
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

theorem majorantLeaf_sound {lg : LogFn} (hlg : Oracle lg)
    {S I A B W : Fix} {s i a b w : ℝ} {u : Int}
    (hs : Mem S s) (hi : Mem I i) (ha : Mem A a) (hb : Mem B b) (hw : Mem W w)
    (h : majorantLeaf lg S I u A B W = true) : MajorantClaim s i a b w := by
  simp only [majorantLeaf, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  rcases h with hinf | ⟨hu, hbody⟩
  · exact fun hf => (majorant_infeasible ha hb hw hinf hf).elim
  intro hf
  obtain ⟨f1, f2, f3, f4, ha0, ha1⟩ := hf
  have huN : (0 : Int) < |u| := abs_pos.mpr hu
  have huR : 0 < ((|u| : Int) : ℝ) / ((scale : Int) : ℝ) :=
    div_pos (by exact_mod_cast huN) scaleR_pos
  have hres : ∃ v, Mem v
      (gObj (((|u| : Int) : ℝ) / ((scale : Int) : ℝ)) a +
        piEval a b + piEval b w - (s * w + i)) ∧ v.hi < 0 := by
    split_ifs at hbody with hpos
    · cases he : fResidualCentered lg (sc (|u|)) A B W S I with
      | none => simp [he] at hbody
      | some v =>
        exact ⟨v, fResidualCentered_mem hlg (sc_mem _) ha hb hw hs hi he,
          by simpa [he] using hbody⟩
    · cases he : fResidual lg (sc (|u|)) A B W S I with
      | none => simp [he] at hbody
      | some v =>
        exact ⟨v, fResidual_mem hlg (sc_mem _) ha hb hw hs hi f1 f2 f3 f4 he,
          by simpa [he] using hbody⟩
  obtain ⟨v, hv, hv0⟩ := hres
  have hneg : (v.hi : ℝ) / ((scale : Int) : ℝ) < 0 :=
    div_neg_of_neg_of_pos (by exact_mod_cast hv0) scaleR_pos
  have hg := gBA_le ha0 ha1 huR
  have hb0 : 0 ≤ b := by linarith
  have hb1 : b ≤ 1 := by linarith
  rw [piBA_eq_piEval_closed ha0 f1 f2, piBA_eq_piEval_closed hb0 f3 f4]
  linarith [hv.2]

theorem majorant_of_checkTree {lg : LogFn} (hlg : Oracle lg)
    {S I A B W : Fix} {s i : ℝ} (hs : Mem S s) (hi : Mem I i)
    (t : BoxTree) (h : checkTree (majorantLeaf lg S I) t A B W = true)
    {a b w : ℝ} (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) :
    MajorantClaim s i a b w :=
  checkTree_sound (P := MajorantClaim s i)
    (fun _ _ _ _ _ _ _ hl ha' hb' hw' => majorantLeaf_sound hlg hs hi ha' hb' hw' hl)
    t A B W a b w h ha hb hw

end Spin.Numeric.Fix
