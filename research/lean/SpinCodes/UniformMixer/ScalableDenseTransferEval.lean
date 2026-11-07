import SpinCodes.UniformMixer.ScalableDenseTransferEvalDefs
import SpinCodes.UniformMixer.ScalableTransferPaperOrientation
import SpinCodes.Numeric.Fixed

/-! Soundness of the fixed-point transfer evaluator against the production map. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Numeric

private theorem transferOne_mem : Fix.Mem transferOne 1 :=
  by simpa [transferOne] using (Fix.ofInt_mem 1)

theorem transferG0_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferG0 β ζ) (1 - b + b * z) := by
  exact Fix.add_mem (Fix.sub_mem transferOne_mem hb)
    (Fix.mul_mem hb hz)

theorem transferG1_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferG1 β ζ) (b + (1 - b) * z) := by
  exact Fix.add_mem hb (Fix.mul_mem (Fix.sub_mem transferOne_mem hb) hz)

theorem transferH1_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferH1 β ζ) (1 - b - b * z) := by
  exact Fix.sub_mem (Fix.sub_mem transferOne_mem hb)
    (Fix.mul_mem hb hz)

theorem transferWeightMoment_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) (w : Nat) :
    Fix.Mem (transferWeightMoment β ζ w) (productionWeightMoment b z w) := by
  exact Fix.mul_mem (Fix.pow_mem (transferG0_mem hb hz) _)
    (Fix.pow_mem (transferG1_mem hb hz) _)

theorem transferFeedbackMoment_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) (w : Nat) :
    Fix.Mem (transferFeedbackMoment β ζ w) (feedbackWeightMoment b z w) := by
  exact Fix.mul_mem (Fix.pow_mem (transferG0_mem hb hz) _)
    (Fix.pow_mem (transferH1_mem hb hz) _)

theorem transferNeutralMoment_mem {b : ℝ} {β : Fix}
    (hb : Fix.Mem β b) (w : Nat) :
    Fix.Mem (transferNeutralMoment β w) ((1 - 2 * b) ^ w) := by
  exact Fix.pow_mem (Fix.sub_mem transferOne_mem
    (Fix.mul_mem (Fix.ofInt_mem 2) hb)) _

private theorem transferFold_mem (terms : List (Nat × Int))
    (f : Nat → Fix) (g : Nat → ℝ)
    (hf : ∀ w, Fix.Mem (f w) (g w))
    (a : Fix) (r : ℝ) (ha : Fix.Mem a r) :
    Fix.Mem
      (terms.foldl
        (fun acc term => Fix.add acc
          (Fix.mul (Fix.ofInt term.2) (f term.1))) a)
      (terms.foldl
        (fun acc term => acc + (term.2 : ℝ) * g term.1) r) := by
  induction terms generalizing a r with
  | nil => simpa using ha
  | cons t ts ih =>
    simp only [List.foldl_cons]
    exact ih _ _ (Fix.add_mem ha
      (Fix.mul_mem (Fix.ofInt_mem t.2) (hf t.1)))

theorem transferSpectrumPoly_mem (f : Nat → Fix) (g : Nat → ℝ)
    (hf : ∀ w, Fix.Mem (f w) (g w)) :
    Fix.Mem (transferSpectrumPoly f) (productionSpectrumPolynomial g) := by
  have hzero : Fix.Mem Fix.zero (0 : ℝ) := by
    simpa [Fix.zero, Fix.sc] using (Fix.sc_mem 0)
  have h := transferFold_mem transferSpectrumTerms f g hf Fix.zero 0 hzero
  convert h using 1
  · rfl
  · simp [transferSpectrumTerms, productionSpectrumPolynomial]

theorem transferEvalH_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferEvalH β ζ) (spectrumH b z) := by
  rw [spectrumH_polynomial]
  have h := Fix.divInt_mem
    (transferSpectrumPoly_mem _ _ (transferFeedbackMoment_mem hb hz))
    (by norm_num : (0 : Int) < 4194304)
  convert h using 1 <;> norm_num [transferEvalH]

theorem transferEvalG_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferEvalG β ζ) (spectrumG b z) := by
  exact Fix.sub_mem (Fix.pow_mem (transferG0_mem hb hz) 64)
    (transferEvalH_mem hb hz)

theorem transferEvalU_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z) :
    Fix.Mem (transferEvalU β ζ) (spectrumU b z) := by
  rw [spectrumU_polynomial]
  rw [state22_liveCount]
  have h := Fix.divInt_mem
    (Fix.sub_mem
      (transferSpectrumPoly_mem _ _ (transferWeightMoment_mem hb hz))
      (transferWeightMoment_mem hb hz 0))
    (by norm_num : (0 : Int) < 4194303)
  convert h using 1; norm_num [transferEvalU]

theorem transferEvalu_mem {b : ℝ} {β : Fix}
    (hb : Fix.Mem β b) :
    Fix.Mem (transferEvalu β) (spectrumu b) := by
  rw [spectrumu_polynomial]
  have h := Fix.divInt_mem
    (transferSpectrumPoly_mem _ _ (transferNeutralMoment_mem hb))
    (by norm_num : (0 : Int) < 4194304)
  have hh := Fix.sub_mem transferOne_mem h
  convert hh using 1 <;> norm_num [transferEvalu]

private theorem initial_le_foldl_max (a : Int) :
    ∀ xs : List Int, a ≤ xs.foldl max a := by
  intro xs
  induction xs generalizing a with
  | nil => exact le_refl a
  | cons x xs ih =>
    simp only [List.foldl_cons]
    exact (le_max_left a x).trans (ih _)

private theorem member_le_foldl_max (a x : Int) (xs : List Int)
    (hx : x ∈ xs) : x ≤ xs.foldl max a := by
  induction xs generalizing a with
  | nil => simp at hx
  | cons y ys ih =>
    simp only [List.mem_cons] at hx
    simp only [List.foldl_cons]
    rcases hx with rfl | htail
    · exact (le_max_right a x).trans (initial_le_foldl_max _ _)
    · exact ih _ htail

private theorem live_weight_in_list {w : Nat}
    (hw : w ∈ liveSpectrumWeights) : w ∈ transferLiveWeights := by
  simp [liveSpectrumWeights, spectrumWeights, transferLiveWeights] at hw ⊢
  tauto

theorem transferEvalF_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z)
    (hb0 : 0 ≤ b) (hb1 : b ≤ 1) (hz0 : 0 ≤ z) :
    Fix.Mem (transferEvalF β ζ) (spectrumF b z) := by
  have hweight : 0 ≤ productionWeightMoment b z 16 := by
    unfold productionWeightMoment
    have h0 : 0 ≤ 1 - b + b * z := by
      nlinarith [mul_nonneg hb0 hz0]
    have h1 : 0 ≤ b + (1 - b) * z := by
      have : 0 ≤ 1 - b := by linarith
      positivity
    positivity
  have hw16 : (16 : Nat) ∈ liveSpectrumWeights := by decide
  have hF0 : 0 ≤ spectrumF b z :=
    hweight.trans (Finset.le_sup' (productionWeightMoment b z) hw16)
  have hFhi : spectrumF b z ≤
      ((transferEvalF β ζ).hi : ℝ) / (scale : ℝ) := by
    unfold spectrumF scalableSpectrumMaxMoment
    apply Finset.sup'_le
    intro w hw
    have hmem : (transferWeightMoment β ζ w).hi ∈
        transferLiveWeights.map
          (fun u => (transferWeightMoment β ζ u).hi) :=
      List.mem_map_of_mem (live_weight_in_list hw)
    have hle := member_le_foldl_max 0
      (transferWeightMoment β ζ w).hi _ hmem
    have hm := (transferWeightMoment_mem hb hz w).2
    exact hm.trans (by
      dsimp only [transferEvalF]
      exact div_le_div_of_nonneg_right (by exact_mod_cast hle) scaleR_pos.le)
  exact ⟨by simpa [Fix.Mem, transferEvalF, Spin.Numeric.scale] using hF0,
    hFhi⟩

theorem transferMin_mem {a b : Fix} {x y : ℝ}
    (ha : Fix.Mem a x) (hb : Fix.Mem b y) :
    Fix.Mem (transferMin a b) (min x y) := by
  unfold transferMin Fix.Mem
  simp only [Int.cast_min]
  constructor
  · rw [← min_div_div_right (le_of_lt scaleR_pos)]
    exact min_le_min ha.1 hb.1
  · rw [← min_div_div_right (le_of_lt scaleR_pos)]
    exact min_le_min ha.2 hb.2

theorem paperSpectrumTransfer_components (b z : ℝ)
    (w : Fin 3 → ℝ) (i : Fin 3) :
    (paperSpectrumTransfer b z).mulVec w i =
      if i = 0 then spectrumH b z * w 0 + spectrumG b z * w 1
      else if i = 1 then
        min (spectrumF b z) (spectrumu b) /
          (liveCount State22 : ℝ) * w 0 + spectrumF b z * w 2
      else
        min (spectrumU b z) (spectrumu b) /
          (liveCount State22 : ℝ) * w 0 + spectrumU b z * w 2 := by
  fin_cases i
  all_goals dsimp [paperSpectrumTransfer, spectrumTransfer,
    weightedTransfer, Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  all_goals norm_num only [Fin.sum_univ_three]
  all_goals simp only [show (0 : Fin 3) ≠ 1 by decide,
    show (1 : Fin 3) ≠ 0 by decide,
    show (2 : Fin 3) ≠ 0 by decide,
    show (2 : Fin 3) ≠ 1 by decide,
    ite_true, ite_false, zero_mul, add_zero]

theorem transferEvalRows_mem {b z : ℝ} {β ζ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z)
    (hb0 : 0 ≤ b) (hb1 : b ≤ 1) (hz0 : 0 ≤ z)
    {v : Fin 3 → Fix} {w : Fin 3 → ℝ}
    (hv : ∀ i, Fix.Mem (v i) (w i)) :
    ∀ i, Fix.Mem (transferEvalRows β ζ v i)
      ((paperSpectrumTransfer b z).mulVec w i) := by
  have hH := transferEvalH_mem hb hz
  have hG := transferEvalG_mem hb hz
  have hF := transferEvalF_mem hb hz hb0 hb1 hz0
  have hU := transferEvalU_mem hb hz
  have hu := transferEvalu_mem hb
  have hminF := transferMin_mem hF hu
  have hminU := transferMin_mem hU hu
  have hrow0 : Fix.Mem (transferEvalRows β ζ v 0)
      (spectrumH b z * w 0 + spectrumG b z * w 1) := by
    simpa [transferEvalRows] using
      (Fix.add_mem (Fix.mul_mem hH (hv 0)) (Fix.mul_mem hG (hv 1)))
  have hrow1 : Fix.Mem (transferEvalRows β ζ v 1)
      (min (spectrumF b z) (spectrumu b) / 4194303 * w 0 +
        spectrumF b z * w 2) := by
    simpa [transferEvalRows] using
      (Fix.add_mem
        (Fix.mul_mem (Fix.divInt_mem hminF
          (by norm_num : (0 : Int) < 4194303)) (hv 0))
        (Fix.mul_mem hF (hv 2)))
  have hrow2 : Fix.Mem (transferEvalRows β ζ v 2)
      (min (spectrumU b z) (spectrumu b) / 4194303 * w 0 +
        spectrumU b z * w 2) := by
    simpa [transferEvalRows] using
      (Fix.add_mem
        (Fix.mul_mem (Fix.divInt_mem hminU
          (by norm_num : (0 : Int) < 4194303)) (hv 0))
        (Fix.mul_mem hU (hv 2)))
  intro i
  fin_cases i
  · simpa [paperSpectrumTransfer_components] using hrow0
  · simpa [paperSpectrumTransfer_components, state22_liveCount] using hrow1
  · simpa [paperSpectrumTransfer_components, state22_liveCount] using hrow2

private theorem le_of_fix_gap {a c : Fix} {x y : ℝ}
    (ha : Fix.Mem a x) (hc : Fix.Mem c y)
    (hgap : a.hi ≤ c.lo) : x ≤ y := by
  exact ha.2.trans ((div_le_div_of_nonneg_right
    (by exact_mod_cast hgap) scaleR_pos.le).trans hc.1)

theorem transferCheck_sound {b z radius : ℝ} {β ζ ρ : Fix}
    (hb : Fix.Mem β b) (hz : Fix.Mem ζ z)
    (hr : Fix.Mem ρ radius)
    (hb0 : 0 ≤ b) (hb1 : b ≤ 1) (hz0 : 0 ≤ z)
    {v : Fin 3 → Fix} {w : Fin 3 → ℝ}
    (hv : ∀ i, Fix.Mem (v i) (w i))
    (hcheck : transferCheck β ζ ρ v = true) :
    ∀ i, (paperSpectrumTransfer b z).mulVec w i ≤ radius * w i := by
  have hrows := transferEvalRows_mem hb hz hb0 hb1 hz0 hv
  simp only [transferCheck, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨h0, h1⟩, h2⟩ := hcheck
  intro i
  fin_cases i
  · exact le_of_fix_gap (hrows 0) (Fix.mul_mem hr (hv 0)) h0
  · exact le_of_fix_gap (hrows 1) (Fix.mul_mem hr (hv 1)) h1
  · exact le_of_fix_gap (hrows 2) (Fix.mul_mem hr (hv 2)) h2

#print axioms transferCheck_sound

end Spin.UniformMixer
