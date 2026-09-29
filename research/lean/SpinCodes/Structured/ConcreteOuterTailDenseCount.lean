import SpinCodes.Structured.ConcreteOuterTailDensePath
import SpinCodes.Structured.ConcreteOuterTailSparse

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

def denseMessages (k : ℕ) : Finset (LocalMessage k) :=
  univ.filter (fun x => x≠0 ∧ k*24<1000*golayOccupation x)

def denseShell (k w : ℕ) : ℝ :=
  ∑ x ∈ denseMessages k, (seedLaw k).prob (fun seed => wtF (encode seed x)=w)

lemma dense_inputWeight {k : ℕ} {x : LocalMessage k} (hx : x∈denseMessages k) :
    k*24 ≤ 125*wtF (outerWord Golay.enc x) := by
  have hd := (mem_filter.mp hx).2.2
  have hl := (golay_inputWeight_bounds x).1
  omega

lemma golay_weight_fiber_count (k a : ℕ) :
    (univ.filter (fun x : LocalMessage k => wtF (outerWord Golay.enc x)=a)).card =
      golayCoefficient k a := by
  simp only [wtF_outerWord]
  exact card_tuples_weight Golay.W k a

lemma denseShell_coefficient_le (k w : ℕ) :
    denseShell k w ≤ ∑ a ∈ Icc 1 (k*24), if k*24 ≤ 125*a then
      (golayCoefficient k a:ℝ)*(∑ c ∈ range (k*24+1), Pt (k*24) a c*Pt (k*24) c w) else 0 := by
  have hmaps : ∀ x ∈ denseMessages k, wtF (outerWord Golay.enc x) ∈ Icc 1 (k*24) := by
    intro x hx
    have hp := golayOccupation_pos (mem_filter.mp hx).2.1
    have hl := (golay_inputWeight_bounds x).1
    exact mem_Icc.mpr ⟨by omega,wtF_le _⟩
  unfold denseShell
  simp only [weight_probability]
  rw [←sum_fiberwise_of_maps_to hmaps]
  apply sum_le_sum
  intro a ha
  by_cases hd : k*24 ≤ 125*a
  · rw [if_pos hd]
    have he : (∑ x ∈ (denseMessages k).filter (fun x => wtF (outerWord Golay.enc x)=a),
        ∑ c ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) c*Pt (k*24) c w) =
        (((denseMessages k).filter (fun x => wtF (outerWord Golay.enc x)=a)).card:ℝ)*
          (∑ c ∈ range (k*24+1), Pt (k*24) a c*Pt (k*24) c w) := by
      rw [←nsmul_eq_mul,←sum_const]
      apply sum_congr rfl
      intro x hx
      rw [(mem_filter.mp hx).2]
    rw [he]
    apply mul_le_mul_of_nonneg_right _ (by apply sum_nonneg; intro c hc; unfold Pt; positivity)
    exact_mod_cast (show ((denseMessages k).filter (fun x => wtF (outerWord Golay.enc x)=a)).card ≤ golayCoefficient k a from by
      rw [←golay_weight_fiber_count]
      apply card_le_card
      intro x hx
      exact mem_filter.mpr ⟨mem_univ _,(mem_filter.mp hx).2⟩)
  · rw [if_neg hd]
    have he : (denseMessages k).filter (fun x => wtF (outerWord Golay.enc x)=a)=∅ := by
      apply filter_eq_empty_iff.mpr
      intro x hx hxa
      exact hd (hxa ▸ dense_inputWeight hx)
    rw [he,sum_empty]

lemma dense_spectrum_term_le {k : ℕ} (hk : 0<k) (w a c : ℕ)
    (ha : a∈Icc 1 (k*24)) (hc : c∈range (k*24+1))
    (hd : k*24 ≤ 125*a) (hw : 125*w ≤ 13*(k*24) ∨ 112*(k*24) ≤ 125*w) :
    (golayCoefficient k a:ℝ)*Pt (k*24) a c*Pt (k*24) c w ≤
      (((k*24:ℕ):ℝ)+1)^2*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) := by
  have hpnon (v t : ℕ) : 0 ≤ Pt (k*24) v t := by unfold Pt; positivity
  by_cases hp : (a,c)∈supportedPaths k w
  · have h₁ := golayCoefficient_le_exp_gBA hk a
    have h₂ := transition_le_exp_entropy (mem_Icc.mp ha).2 c
    have h₃ := transition_le_exp_entropy (show c≤k*24 by have := mem_range.mp hc; omega) w
    have hE := dense_pathEntropy_le hk hp
      (show ((k*24:ℕ):ℝ)/125 ≤ a by
        have hh : ((k*24:ℕ):ℝ) ≤ 125*(a:ℝ) := by exact_mod_cast hd
        linarith) hw
    calc
      (golayCoefficient k a:ℝ)*Pt (k*24) a c*Pt (k*24) c w ≤
          Real.exp (((k*24:ℕ):ℝ)*gBA ((a:ℝ)/(k*24:ℕ)))*
            ((((k*24:ℕ):ℝ)+1)*Real.exp (transitionEntropy (k*24) a c))*
            ((((k*24:ℕ):ℝ)+1)*Real.exp (transitionEntropy (k*24) c w)) :=
        mul_le_mul (mul_le_mul h₁ h₂ (hpnon _ _) (by positivity)) h₃ (hpnon _ _) (by positivity)
      _ = (((k*24:ℕ):ℝ)+1)^2*Real.exp (pathEntropy k w (a,c)) := by
        rw [pathEntropy,Real.exp_add,Real.exp_add]; ring
      _ ≤ _ := mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr hE) (sq_nonneg _)
  · have hz : golayCoefficient k a=0 ∨ accT (k*24) a c=0 ∨ accT (k*24) c w=0 := by
      simpa only [supportedPaths,mem_filter,mem_product,ha,hc,true_and,not_and_or,not_not] using hp
    rcases hz with hz|hz|hz <;> simp only [Pt,hz,Nat.cast_zero,zero_mul,mul_zero,zero_div] <;> positivity

theorem denseShell_le {k : ℕ} (hk : 0<k) (w : ℕ)
    (hw : 125*w ≤ 13*(k*24) ∨ 112*(k*24) ≤ 125*w) :
    denseShell k w ≤ (((k*24:ℕ):ℝ)+1)^4*
      Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) := by
  refine (denseShell_coefficient_le k w).trans ?_
  have hh : (∑ a ∈ Icc 1 (k*24), if k*24 ≤ 125*a then
      (golayCoefficient k a:ℝ)*(∑ c ∈ range (k*24+1), Pt (k*24) a c*Pt (k*24) c w) else 0) ≤
      ∑ _a ∈ Icc 1 (k*24), ∑ _c ∈ range (k*24+1),
        (((k*24:ℕ):ℝ)+1)^2*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) := by
    apply sum_le_sum
    intro a ha
    split_ifs with hd
    · rw [mul_sum]
      apply sum_le_sum
      intro c hc
      simpa only [mul_assoc] using dense_spectrum_term_le hk w a c ha hc hd hw
    · positivity
  simp only [sum_const,card_range,Nat.card_Icc,Nat.add_sub_cancel,nsmul_eq_mul,Nat.cast_add,Nat.cast_one] at hh
  refine hh.trans ?_
  have h := mul_le_mul_of_nonneg_right (show ((k*24:ℕ):ℝ) ≤ ((k*24:ℕ):ℝ)+1 by linarith)
    (show 0 ≤ (((k*24:ℕ):ℝ)+1)^3*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) by positivity)
  convert h using 1 <;> ring

end Spin.Structured.ConcreteOuter


