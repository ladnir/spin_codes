import SpinCodes.Structured.ConcreteOuterTailUpper
import SpinCodes.Structured.ConcreteOuterTailSeriesFinite

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset Filter

def firstMoment {k : ℕ} (x : LocalMessage k) : ℝ :=
  ∑ h ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) h*zeta^h

def sparseMessages (k : ℕ) : Finset (LocalMessage k) :=
  univ.filter (fun x => x≠0 ∧ 1000*golayOccupation x ≤ k*24)

theorem sum_sparse_moments_le {k : ℕ} (hk : 0<k) :
    (∑ x ∈ sparseMessages k, firstMoment x) ≤ sparseBoundSeries (k*24) := by
  have hmaps : ∀ x ∈ sparseMessages k, golayOccupation x ∈ Ico 1 (k*24/1000+1) := by
    intro x hx
    obtain ⟨_,hx0,hsp⟩ := mem_filter.mp hx
    have hp := golayOccupation_pos hx0
    exact mem_Ico.mpr ⟨hp, by omega⟩
  rw [←sum_fiberwise_of_maps_to hmaps]
  unfold sparseBoundSeries
  apply sum_le_sum
  intro j hj
  have hj0 : 0<j := (mem_Ico.mp hj).1
  have hsp : 1000*j≤k*24 := by have := (mem_Ico.mp hj).2; omega
  have he : (sparseMessages k).filter (fun x => golayOccupation x=j) =
      univ.filter (fun x : LocalMessage k => golayOccupation x=j) := by
    ext x
    simp only [sparseMessages,mem_filter,mem_univ,true_and]
    constructor
    · exact fun h => h.2
    · intro hx
      have hx0 : x≠0 := by
        intro hz
        have hh : golayOccupation x=0 := by rw [hz]; simp [golayOccupation,golayActive]
        omega
      exact ⟨⟨hx0, by omega⟩,hx⟩
  rw [he]
  exact golay_occupation_moment_le hk hj0 hsp

lemma weight_probability {k : ℕ} (x : LocalMessage k) (w : ℕ) :
    (seedLaw k).prob (fun seed => wtF (encode seed x)=w) =
      ∑ c ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) c * Pt (k*24) c w := by
  simpa only [seedLaw,encode,wtF_accF] using prob_two_stage (w:=w) (outerWord Golay.enc x) rfl

lemma double_upperMass_le_moment {b a D : ℕ} (ha : 0<a) (hDb : 125*D≤13*b) :
    (∑ w ∈ Icc (b-D) b, ∑ c ∈ range (b+1), Pt b a c*Pt b c w) ≤
      ((b:ℝ)*(zeta⁻¹^2+1))*(∑ c ∈ range (b+1), Pt b a c*zeta^c) := by
  rw [sum_comm,mul_sum]
  apply sum_le_sum
  intro c hc
  rw [←mul_sum]
  by_cases hc0 : c=0
  · subst c
    simp [Pt,accT_zero_right b (by omega : a≠0)]
  · have hcb : c≤b := by have := mem_range.mp hc; omega
    have hh := mul_le_mul_of_nonneg_left (upperMass_le_zeta (by omega) hcb hDb)
      (show 0≤Pt b a c by unfold Pt; positivity)
    simpa only [upperMass,mul_assoc,mul_left_comm] using hh

def sparseLowerTail (k : ℕ) : ℝ :=
  ∑ x ∈ sparseMessages k, ∑ w ∈ range ((k*24)*13/125+1),
    (seedLaw k).prob (fun seed => wtF (encode seed x)=w)

def sparseUpperTail (k : ℕ) : ℝ :=
  ∑ x ∈ sparseMessages k, ∑ w ∈ Icc (k*24-(k*24)*13/125) (k*24),
    (seedLaw k).prob (fun seed => wtF (encode seed x)=w)

theorem sparseLowerTail_le {k : ℕ} (hk : 0<k) :
    sparseLowerTail k ≤ zeta⁻¹*sparseBoundSeries (k*24) := by
  calc
    sparseLowerTail k ≤ zeta⁻¹*(∑ x ∈ sparseMessages k,firstMoment x) := by
      unfold sparseLowerTail
      rw [mul_sum]
      apply sum_le_sum
      intro x hx
      simp only [weight_probability]
      exact double_lowerMass_le_moment (k*24) (wtF (outerWord Golay.enc x)) _ (by omega)
    _ ≤ _ := mul_le_mul_of_nonneg_left (sum_sparse_moments_le hk) (inv_nonneg.mpr zeta_pos.le)

theorem sparseUpperTail_le {k : ℕ} (hk : 0<k) :
    sparseUpperTail k ≤ (zeta⁻¹^2+1)*((k*24:ℕ):ℝ)*sparseBoundSeries (k*24) := by
  have hh : sparseUpperTail k ≤ ((k*24:ℕ):ℝ)*(zeta⁻¹^2+1)*(∑ x ∈ sparseMessages k,firstMoment x) := by
    unfold sparseUpperTail
    rw [mul_sum]
    apply sum_le_sum
    intro x hx
    have hx0 := (mem_filter.mp hx).2.1
    have ha : 0<wtF (outerWord Golay.enc x) := by
      have hp := golayOccupation_pos hx0
      have hl := (golay_inputWeight_bounds x).1
      omega
    simp only [weight_probability]
    exact double_upperMass_le_moment ha (by omega)
  have hb := hh.trans (mul_le_mul_of_nonneg_left (sum_sparse_moments_le hk) (by positivity))
  convert hb using 1 <;> ring

end Spin.Structured.ConcreteOuter

