import SpinCodes.Numeric.MajorantEval

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Spin.Numeric.Fix

def ClosedDenseClaim (thr : Int) (a b w : ℝ) : Prop :=
  (a/2≤b ∧ b≤1-a/2 ∧ b/2≤w ∧ w≤1-b/2 ∧ 0≤a ∧ a≤1) →
    gBA a+piBA a b+piBA b w < (thr:ℝ)/(scale:ℝ)

lemma infeasible_closed_vacuous {A B W : Fix} {a b w : ℝ} (ha : Mem A a) (hb : Mem B b)
    (hw : Mem W w) (h : infeasible A B W = true) (thr : Int) :
    ClosedDenseClaim thr a b w := by
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
theorem leafOK_closed_sound {lg : LogFn} (hlg : Oracle lg) {thr u : Int} {A B W : Fix}
    {a b w : ℝ} (ha : Mem A a) (hb : Mem B b) (hw : Mem W w)
    (h : leafOK lg thr u A B W = true) : ClosedDenseClaim thr a b w := by
  simp only [leafOK, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  rcases h with hinf | ⟨hu, hbody⟩
  · exact infeasible_closed_vacuous ha hb hw hinf thr
  intro hfeas
  obtain ⟨f1, f2, f3, f4, ha0, ha1⟩ := hfeas
  have hS := scaleR_pos
  have hb0 : 0 ≤ b := by linarith
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
  have e1 := fpiBestClampL_closed_mem hlg ha hb f1 f2 hp1
  have e2 := fpiBestClampL_closed_mem hlg hb hw f3 f4 hp2
  have esum := add_mem (add_mem eg e1) e2
  -- the sum is below the threshold
  have hlt : gObj ((u : ℝ) / ((scale : Int) : ℝ)) a + piEval a b + piEval b w
      < (thr : ℝ) / ((scale : Int) : ℝ) := by
    refine lt_of_le_of_lt esum.2 ?_
    have hnum : ((add (add g p1) p2).hi : ℝ) < (thr : ℝ) := by exact_mod_cast hbody
    gcongr
  -- and `g` and `π` are dominated by what was evaluated
  have hgle : gBA a ≤ gObj ((u : ℝ) / ((scale : Int) : ℝ)) a := gBA_le ha0 ha1 huR
  have hpi1 : piBA a b = piEval a b := piBA_eq_piEval_closed ha0 f1 f2
  have hpi2 : piBA b w = piEval b w := piBA_eq_piEval_closed hb0 f3 f4
  rw [hpi1, hpi2]
  linarith

/-- **The bridge.**  A tree that checks out proves the dense-tail claim at
every point of the root box.

Everything specific to this certificate is in `leafOK`; the covering is
`checkTree_sound`, which knows nothing about `π` or `g`. -/
theorem denseTail_closed_of_checkTree {lg : LogFn} (hlg : Oracle lg) {thr : Int}
    (t : BoxTree) {A B W : Fix} (h : checkTree (leafOK lg thr) t A B W = true)
    {a b w : ℝ} (ha : Mem A a) (hb : Mem B b) (hw : Mem W w) :
    ClosedDenseClaim thr a b w :=
  checkTree_sound (P := ClosedDenseClaim thr)
    (fun _ _ _ _ _ _ _ hl ha' hb' hw' => leafOK_closed_sound hlg ha' hb' hw' hl)
    t A B W a b w h ha hb hw


end Spin.Structured.ConcreteOuter
