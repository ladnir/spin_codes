import SpinCodes.UniformMixer.ScalableDenseLeaf
import SpinCodes.UniformMixer.ScalableDenseCoverCore

/-! A compact checked leaf for the revised Scalable dense cover.

The certificate generator supplies rational data. `densePackedLeafCheck` is
decidable, and the theorem below connects its answer to the actual sampled
inner law through `scalable_dense_checked_leaf`. The generator is not trusted.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.DenseOccupationFixed
open Spin.Structured.DenseGeometry

structure DenseLine where
  slope : QInput
  intercept : QInput

structure DenseWitness where
  p : QInput
  y : QInput
  ρ : QInput
  z : QInput
  vector : Fin 3 → QInput

structure DenseBox where
  a₀ : QInput
  a₁ : QInput
  x₀ : QInput
  x₁ : QInput

def DenseLine.rat (l : DenseLine) : ℚ × ℚ :=
  (transferQInputRat l.slope, transferQInputRat l.intercept)

def DenseBox.rect (b : DenseBox) : Rect :=
  ⟨transferQInputRat b.a₀, transferQInputRat b.a₁,
   transferQInputRat b.x₀, transferQInputRat b.x₁⟩

def DenseLine.Valid (l : DenseLine) : Prop :=
  0 < l.slope.den ∧ 0 < l.intercept.den ∧
    l.rat ∈ Spin.Majorant.refined.supports

def DenseWitness.Valid (w : DenseWitness) : Prop :=
  0 < w.p.num ∧ w.p.num < w.p.den ∧
  0 < w.y.num ∧ w.y.num < w.y.den ∧
  0 < transferQInputRat w.p ∧ transferQInputRat w.p < 1 ∧
  0 < transferQInputRat w.y ∧ transferQInputRat w.y < 1 ∧
  0 < transferQInputRat w.ρ ∧
  0 < transferQInputRat w.z ∧ transferQInputRat w.z ≤ 1 ∧
  transferQInputRat (w.vector 0) = 1 ∧
  (∀ i : Fin 3, 1 ≤ 10000000 * transferQInputRat (w.vector i))

def DenseBox.Valid (b : DenseBox) : Prop :=
  0 < b.a₀.den ∧ 0 < b.a₁.den ∧
  0 < b.x₀.den ∧ 0 < b.x₁.den ∧
  0 < transferQInputRat b.a₀ ∧ transferQInputRat b.a₁ ≤ 1 ∧
  0 < transferQInputRat b.x₀ ∧ transferQInputRat b.x₁ ≤ 1

instance (l : DenseLine) : Decidable l.Valid := by
  unfold DenseLine.Valid
  infer_instance

instance (w : DenseWitness) : Decidable w.Valid := by
  unfold DenseWitness.Valid
  infer_instance

instance (b : DenseBox) : Decidable b.Valid := by
  unfold DenseBox.Valid
  infer_instance

def densePackedLeafCheck (l : DenseLine) (w : DenseWitness)
    (b : DenseBox) : Bool :=
  decide (l.Valid ∧ w.Valid ∧ b.Valid ∧
    transferRatCheck
      (transferQInputRat w.p * transferQInputRat w.y)
      (transferQInputRat w.z) (transferQInputRat w.ρ)
      (fun i => transferQInputRat (w.vector i)) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₀ b.x₀ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₀ b.x₁ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₁ b.x₀ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₁ b.x₁ 50 (-100000000000000000000000) = true)

private theorem qInput_real_pos (a : QInput)
    (h : 0 < transferQInputRat a) : 0 < a.real := by
  have hh : (0 : ℝ) < (transferQInputRat a : ℝ) := by exact_mod_cast h
  simpa only [transferQInputRat_cast] using hh

private theorem qInput_real_lt_one (a : QInput)
    (h : transferQInputRat a < 1) : a.real < 1 := by
  have hh : (transferQInputRat a : ℝ) < 1 := by exact_mod_cast h
  simpa only [transferQInputRat_cast] using hh

private theorem qInput_real_le_one (a : QInput)
    (h : transferQInputRat a ≤ 1) : a.real ≤ 1 := by
  have hh : (transferQInputRat a : ℝ) ≤ 1 := by exact_mod_cast h
  simpa only [transferQInputRat_cast] using hh

theorem densePackedLeafCheck_sound (l : DenseLine) (w : DenseWitness)
    (b : DenseBox) (hc : densePackedLeafCheck l w b = true) :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) b.rect := by
  have h := of_decide_eq_true hc
  rcases h with ⟨⟨hs, hi, hsupport⟩,
    ⟨hpnum, hpnum1, hynum, hynum1,
      hp0, hp1, hy0, hy1, hρ0, hz0, hz1, hv0, hones⟩,
    ⟨ha₀d, ha₁d, hx₀d, hx₁d, ha₀0, ha₁1, hx₀0, hx₁1⟩,
    ht, h00, h01, h10, h11⟩
  intro α x hr
  have hα : α ∈ Set.Icc b.a₀.real b.a₁.real := by
    simpa [DenseBox.rect, Rect.Contains, transferQInputRat_cast] using
      (show (b.rect.alo : ℝ) ≤ α ∧ α ≤ (b.rect.ahi : ℝ) from ⟨hr.1, hr.2.1⟩)
  have hx : x ∈ Set.Icc b.x₀.real b.x₁.real := by
    simpa [DenseBox.rect, Rect.Contains, transferQInputRat_cast] using
      (show (b.rect.xlo : ℝ) ≤ x ∧ x ≤ (b.rect.xhi : ℝ) from
        ⟨hr.2.2.1, hr.2.2.2⟩)
  have honesReal : ∀ i, (1 : ℝ) ≤ 10000000 * (w.vector i).real := by
    intro i
    have hi := hones i
    have hh : (1 : ℝ) ≤ 10000000 *
        (transferQInputRat (w.vector i) : ℝ) := by exact_mod_cast hi
    simpa only [transferQInputRat_cast] using hh
  have hrate := scalable_dense_checked_leaf
    l.slope l.intercept w.p w.y w.ρ w.z
    b.a₀ b.a₁ b.x₀ b.x₁ w.vector (10000000 : ℝ)
    hs hi ha₀d ha₁d hx₀d hx₁d
    hpnum hpnum1 hynum hynum1
    (qInput_real_pos w.p hp0) (qInput_real_lt_one w.p hp1)
    (qInput_real_pos w.y hy0) (qInput_real_lt_one w.y hy1)
    (qInput_real_pos w.ρ hρ0) (qInput_real_pos w.z hz0)
    (qInput_real_le_one w.z hz1) (qInput_real_pos b.a₀ ha₀0)
    (qInput_real_le_one b.a₁ ha₁1) (qInput_real_pos b.x₀ hx₀0)
    (qInput_real_le_one b.x₁ hx₁1)
    (by norm_num) honesReal ht h00 h01 h10 h11 hα hx
  refine ⟨l.rat, hsupport, ?_⟩
  have hv0Real : (w.vector 0).real = (1 : ℝ) := by
    have hh : (transferQInputRat (w.vector 0) : ℝ) = 1 := by
      exact_mod_cast hv0
    simpa only [transferQInputRat_cast] using hh
  simpa [DenseLine.rat, transferQInputRat_cast, hv0Real] using hrate

structure DenseLeafRef where
  witnessIndex : Nat
  box : DenseBox

/-- An absent witness is a failed check, never a default witness. -/
def denseLeafRefCheck (l : DenseLine) (ws : List DenseWitness)
    (r : DenseLeafRef) : Bool :=
  match ws[r.witnessIndex]? with
  | none => false
  | some w => densePackedLeafCheck l w r.box

theorem denseLeafRefCheck_sound (l : DenseLine) (ws : List DenseWitness)
    (r : DenseLeafRef) (hc : denseLeafRefCheck l ws r = true) :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) r.box.rect := by
  unfold denseLeafRefCheck at hc
  cases hw : ws[r.witnessIndex]? with
  | none => simp [hw] at hc
  | some w =>
    exact densePackedLeafCheck_sound l w r.box (by simpa [hw] using hc)

def denseChunkCheck (l : DenseLine) (ws : List DenseWitness)
    (refs : List DenseLeafRef) (tree : Tree) (root : Rect) : Bool :=
  check (refs.map (fun r => r.box.rect)) tree root &&
    refs.all (denseLeafRefCheck l ws)

theorem denseChunkCheck_sound (l : DenseLine) (ws : List DenseWitness)
    (refs : List DenseLeafRef) (tree : Tree) (root : Rect)
    (hc : denseChunkCheck l ws refs tree root = true) :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  simp only [denseChunkCheck, Bool.and_eq_true] at hc
  obtain ⟨hgeometry, hnumeric⟩ := hc
  apply denseClaim_of_checked_tree
    (DenseTarget (1 / 10000000) 10000000)
    (refs.map (fun r => r.box.rect)) tree root hgeometry
  intro i hi
  have hi' : i < refs.length := by simpa using hi
  have hmem : refs[i] ∈ refs := List.getElem_mem hi'
  have hall := List.all_eq_true.mp hnumeric refs[i] hmem
  have h := denseLeafRefCheck_sound l ws refs[i] hall
  simpa only [List.getD_eq_getElem _ _ hi, List.getElem_map] using h

def denseWitnessCheck (w : DenseWitness) : Bool :=
  decide (w.Valid ∧ transferRatCheck
    (transferQInputRat w.p * transferQInputRat w.y)
    (transferQInputRat w.z) (transferQInputRat w.ρ)
    (fun i => transferQInputRat (w.vector i)) = true)

def denseBoxCheck (l : DenseLine) (w : DenseWitness)
    (b : DenseBox) : Bool :=
  decide (b.Valid ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₀ b.x₀ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₀ b.x₁ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₁ b.x₀ 50 (-100000000000000000000000) = true ∧
    scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
      b.a₁ b.x₁ 50 (-100000000000000000000000) = true)

def denseRefBoxCheck (l : DenseLine) (ws : List DenseWitness)
    (r : DenseLeafRef) : Bool :=
  match ws[r.witnessIndex]? with
  | none => false
  | some w => denseBoxCheck l w r.box

def denseChunkCheckShared (l : DenseLine) (ws : List DenseWitness)
    (refs : List DenseLeafRef) (tree : Tree) (root : Rect) : Bool :=
  check (refs.map (fun r => r.box.rect)) tree root &&
    decide l.Valid && ws.all denseWitnessCheck &&
    refs.all (denseRefBoxCheck l ws)

theorem denseChunkCheckShared_sound (l : DenseLine) (ws : List DenseWitness)
    (refs : List DenseLeafRef) (tree : Tree) (root : Rect)
    (hc : denseChunkCheckShared l ws refs tree root = true) :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  simp only [denseChunkCheckShared, Bool.and_eq_true] at hc
  obtain ⟨⟨⟨hgeometry, hline⟩, hwitnesses⟩, hboxes⟩ := hc
  have hline' : l.Valid := of_decide_eq_true hline
  have hwitnesses' := List.all_eq_true.mp hwitnesses
  have hboxes' := List.all_eq_true.mp hboxes
  apply denseClaim_of_checked_tree
    (DenseTarget (1 / 10000000) 10000000)
    (refs.map (fun r => r.box.rect)) tree root hgeometry
  intro i hi
  have hi' : i < refs.length := by simpa using hi
  have hmem : refs[i] ∈ refs := List.getElem_mem hi'
  have hb := hboxes' refs[i] hmem
  unfold denseRefBoxCheck at hb
  cases hget : ws[refs[i].witnessIndex]? with
  | none => simp [hget] at hb
  | some w =>
    have hwmem : w ∈ ws := List.mem_of_getElem? hget
    have hw := hwitnesses' w hwmem
    have hw' : w.Valid ∧ transferRatCheck
        (transferQInputRat w.p * transferQInputRat w.y)
        (transferQInputRat w.z) (transferQInputRat w.ρ)
        (fun j => transferQInputRat (w.vector j)) = true :=
      of_decide_eq_true (by simpa [denseWitnessCheck] using hw)
    have hb' : refs[i].box.Valid ∧
        scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
          refs[i].box.a₀ refs[i].box.x₀ 50
            (-100000000000000000000000) = true ∧
        scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
          refs[i].box.a₀ refs[i].box.x₁ 50
            (-100000000000000000000000) = true ∧
        scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
          refs[i].box.a₁ refs[i].box.x₀ 50
            (-100000000000000000000000) = true ∧
        scalableVertexCheck l.slope l.intercept w.p w.y w.ρ w.z
          refs[i].box.a₁ refs[i].box.x₁ 50
            (-100000000000000000000000) = true :=
      of_decide_eq_true (by simpa [hget, denseBoxCheck] using hb)
    have hfull : densePackedLeafCheck l w refs[i].box = true := by
      apply decide_eq_true_eq.mpr
      exact ⟨hline', hw'.1, hb'.1, hw'.2,
        hb'.2.1, hb'.2.2.1, hb'.2.2.2.1, hb'.2.2.2.2⟩
    have h := densePackedLeafCheck_sound l w refs[i].box hfull
    simpa only [List.getD_eq_getElem _ _ hi, List.getElem_map] using h

#print axioms densePackedLeafCheck_sound
#print axioms denseChunkCheck_sound
#print axioms denseChunkCheckShared_sound

end Spin.UniformMixer
