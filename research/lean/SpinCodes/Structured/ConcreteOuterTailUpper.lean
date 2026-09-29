import SpinCodes.Structured.ConcreteOuterTailTransition
import SpinCodes.Structured.UpperTail

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

lemma lowerCount_le_zeta {b v D : ℕ} (hDb : 125*D≤13*b) :
    ((∑ h ∈ range (D+1), accT b v h:ℕ):ℝ) ≤ (b.choose v:ℝ)*(zeta⁻¹*zeta^v) := by
  have hh : ((∑ h ∈ range (D+1), accT b v h:ℕ):ℝ)*(28:ℝ)^(v/2) ≤
      (13:ℝ)^(v/2)*(b.choose v:ℝ) := by exact_mod_cast sum_accT_lower_uniform hDb
  have hp : ((∑ h ∈ range (D+1), accT b v h:ℕ):ℝ) ≤ (b.choose v:ℝ)*(13/28:ℝ)^(v/2) := by
    rw [div_pow, ←mul_div_assoc]
    apply (le_div_iff₀ (by positivity)).mpr
    nlinarith
  exact hp.trans (mul_le_mul_of_nonneg_left (half_power_le_zeta v) (by positivity))

lemma sum_accT_upper_eq_card (b v D : ℕ) :
    ∑ h ∈ Icc (b-D) b, accT b v h =
      (univ.filter (fun x : Fin b→Bool => wtF x=v ∧ b-D≤accWtF x)).card := by
  have hmaps : ∀ x ∈ univ.filter (fun x : Fin b→Bool => wtF x=v ∧ b-D≤accWtF x),
      accWtF x ∈ Icc (b-D) b := by
    intro x hx
    exact mem_Icc.mpr ⟨(mem_filter.mp hx).2.2, accWtF_le x⟩
  rw [Finset.card_eq_sum_card_fiberwise hmaps]
  apply sum_congr rfl
  intro h hh
  rw [←accCountF_eq_accT]
  unfold accCountF
  apply congrArg Finset.card
  ext x
  simp only [mem_filter, mem_univ, true_and]
  have hlo := (mem_Icc.mp hh).1
  constructor
  · rintro ⟨hv,hc⟩
    exact ⟨⟨hv, by omega⟩,hc⟩
  · rintro ⟨⟨hv,_⟩,hc⟩
    exact ⟨hv,hc⟩

def upperMass (b v D : ℕ) : ℝ := ∑ h ∈ Icc (b-D) b, Pt b v h

/-- The complementary tail, uniformly in every nonzero supported input weight. -/
theorem upperMass_le_zeta {b v D : ℕ} (hv : 0<v) (hvb : v≤b)
    (hDb : 125*D≤13*b) :
    upperMass b v D ≤ (b:ℝ)*(zeta⁻¹^2+1)*zeta^v := by
  have hb : 0<b := by omega
  obtain ⟨b',rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : b≠0)
  have hc : (0:ℝ)<(b'+1).choose v := by exact_mod_cast Nat.choose_pos hvb
  have hD : D≤b'+1 := by omega
  have hu := upper_tail_le_sums (b:=b') hv hD
  have huR : ((univ.filter (fun x : Fin (b'+1)→Bool => wtF x=v ∧ b'+1-D≤accWtF x)).card:ℝ) ≤
      ((∑ h ∈ range (D+1), accT (b'+1) (v-1) h:ℕ):ℝ)+
      ((∑ h ∈ range (D+1), accT (b'+1) (v+1) h:ℕ):ℝ) := by exact_mod_cast hu
  have hpre := lowerCount_le_zeta (v:=v-1) hDb
  have hsuc := lowerCount_le_zeta (v:=v+1) hDb
  have hcp : ((b'+1).choose (v-1):ℝ) ≤ ((b'+1:ℕ):ℝ)*((b'+1).choose v:ℝ) := by
    exact_mod_cast choose_pred_le_mul hv hvb
  have hcs : ((b'+1).choose (v+1):ℝ) ≤ ((b'+1:ℕ):ℝ)*((b'+1).choose v:ℝ) := by
    exact_mod_cast choose_succ_le_mul (b'+1) v
  have hz := zeta_pos
  have hp' := hpre.trans (mul_le_mul_of_nonneg_right hcp (by positivity))
  have hs' := hsuc.trans (mul_le_mul_of_nonneg_right hcs (by positivity))
  have he : zeta⁻¹*zeta^(v-1)+zeta⁻¹*zeta^(v+1)=(zeta⁻¹^2+1)*zeta^v := by
    have hv' : v-1+1=v := by omega
    have hh : zeta^(v-1)*zeta=zeta^v := by rw [←pow_succ, hv']
    rw [pow_succ]
    field_simp
    nlinarith [hh]
  unfold upperMass
  simp only [Pt, ←sum_div, ←Nat.cast_sum, sum_accT_upper_eq_card]
  apply (div_le_iff₀ hc).mpr
  have hh := huR.trans (add_le_add hp' hs')
  rw [←mul_add, he] at hh
  convert hh using 1 <;> ring

end Spin.Structured.ConcreteOuter

