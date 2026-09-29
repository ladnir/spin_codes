import SpinCodes.Numeric.BAEval

namespace Spin.Structured.ConcreteOuter
open Spin.Numeric
open Spin.Numeric.Fix

/-- A direct fixed-point certificate for an affine entropy defect. -/
def defectValue (x s i : Fix) : Option Fix :=
  (fhEnt x 24).map (fun h => add (sub (add (mul s x) i) h) (divInt log2 2))

def defectOK (x s i : Fix) : Bool :=
  match defectValue x s i with
  | none => false
  | some v => decide (v.hi ≤ (ofFrac 1281 100000).lo)

theorem defectOK_sound {X S I : Fix} {x s i : ℝ}
    (hx : Mem X x) (hs : Mem S s) (hi : Mem I i)
    (h : defectOK X S I = true) :
    s*x+i-hEnt x+Real.log 2/2 ≤ 1281/100000 := by
  unfold defectOK defectValue at h
  cases he : fhEnt X 24 with
  | none => simp [he] at h
  | some v =>
    simp only [he, Option.map_some, decide_eq_true_eq] at h
    have hm := add_mem (sub_mem (add_mem (mul_mem hs hx) hi) (fhEnt_mem hx he))
      (divInt_mem log2_mem (by norm_num : (0:Int)<2))
    norm_num only [Int.cast_ofNat] at hm
    have hbound : ((add (sub (add (mul S X) I) v) (divInt log2 2)).hi:ℝ) /
        (scale:ℝ) ≤ ((ofFrac 1281 100000).lo:ℝ)/(scale:ℝ) := by
      exact div_le_div_of_nonneg_right (by exact_mod_cast h) scaleR_pos.le
    have hfrac := (ofFrac_mem (p:=1281) (q:=100000) (by norm_num)).1
    norm_num only [Int.cast_ofNat] at hfrac
    exact hm.2.trans (hbound.trans hfrac)

lemma hEnt_eq_binEntropy (x : ℝ) : hEnt x = Real.binEntropy x := by
  simp only [hEnt, Real.binEntropy, Real.log_inv]
  ring

lemma entropy_defect_convex (s i : ℝ) :
    ConvexOn ℝ (Set.Icc (0:ℝ) 1) (fun x => s*x+i-hEnt x) := by
  have hc : ConcaveOn ℝ (Set.Icc (0:ℝ) 1) hEnt := by
    simpa only [funext hEnt_eq_binEntropy] using Real.strictConcave_binEntropy.concaveOn
  refine ⟨convex_Icc 0 1, ?_⟩
  intro x hx y hy a b ha hb hab
  have hh := hc.2 hx hy ha hb hab
  simp only [smul_eq_mul] at hh ⊢
  nlinarith [congrArg (fun t : ℝ => i*t) hab]

/-- One affine support controls the entire interval from its two endpoints. -/
theorem entropy_defect_interval {s i lo hi x d : ℝ}
    (hlo : 0 ≤ lo) (hhi : hi ≤ 1) (hx : x ∈ Set.Icc lo hi)
    (hl : s*lo+i-hEnt lo ≤ d) (hh : s*hi+i-hEnt hi ≤ d) :
    s*x+i-hEnt x ≤ d := by
  have hlo' : lo ∈ Set.Icc (0:ℝ) 1 := ⟨hlo, (hx.1.trans hx.2).trans hhi⟩
  have hhi' : hi ∈ Set.Icc (0:ℝ) 1 := ⟨hlo.trans (hx.1.trans hx.2), hhi⟩
  exact ((entropy_defect_convex s i).le_max_of_mem_Icc hlo' hhi' hx).trans (max_le hl hh)

end Spin.Structured.ConcreteOuter

