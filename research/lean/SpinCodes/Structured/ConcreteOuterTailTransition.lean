import SpinCodes.Structured.BATails
import SpinCodes.Structured.Composition
import SpinCodes.Numeric.BAConstants

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

/-- A uniform integer lower-tail bound; empty support and both parities included. -/
theorem sum_accT_lower_uniform {b v D : ℕ} (hDb : 125*D ≤ 13*b) :
    (∑ h ∈ range (D+1), accT b v h)*28^(v/2) ≤ 13^(v/2)*b.choose v := by
  by_cases hv : v=0
  · subst v
    simp [accT_zero_left]
  by_cases hp : v%2=0
  · have he : v=2*(v/2) := by omega
    rw [he]
    simp only [Nat.mul_div_cancel_left _ (by omega : 0<2)]
    by_cases hD : v/2 ≤ D
    · exact sum_accT_even_le_pow (by omega) hD hDb
    · rw [sum_accT_even_eq_zero (by omega) (by omega), zero_mul]
      positivity
  · have he : v=2*(v/2)+1 := by omega
    rw [he]
    have he' : (2*(v/2)+1)/2=v/2 := by omega
    rw [he']
    by_cases hD : v/2 ≤ D
    · exact sum_accT_odd_le_pow' hD hDb
    · have hz : ∑ h ∈ range (D+1), accT b (2*(v/2)+1) h = 0 := by
        apply Nat.le_zero.mp
        have hh := sum_accT_odd_le (b:=b) (ℓ:=v/2) (D:=D)
        simpa only [Nat.choose_eq_zero_of_lt (by omega : D<v/2+1), mul_zero] using hh
      rw [hz, zero_mul]
      positivity

def lowerMass (b v D : ℕ) : ℝ := ∑ h ∈ range (D+1), Pt b v h

lemma lowerMass_nonneg (b v D : ℕ) : 0 ≤ lowerMass b v D := by
  apply Finset.sum_nonneg
  intro h _
  unfold Pt
  positivity

lemma lowerMass_eq (b v D : ℕ) :
    lowerMass b v D = ((∑ h ∈ range (D+1), accT b v h:ℕ):ℝ)/(b.choose v:ℝ) := by
  simp only [lowerMass, Pt, ←Finset.sum_div, Nat.cast_sum]

theorem lowerMass_le_pow {b v D : ℕ} (hv : v ≤ b) (hDb : 125*D ≤ 13*b) :
    lowerMass b v D ≤ (13/28:ℝ)^(v/2) := by
  have hc : (0:ℝ)<b.choose v := by exact_mod_cast Nat.choose_pos hv
  have hh : ((∑ h ∈ range (D+1), accT b v h:ℕ):ℝ)*(28:ℝ)^(v/2) ≤
      (13:ℝ)^(v/2)*(b.choose v:ℝ) := by exact_mod_cast sum_accT_lower_uniform hDb
  rw [lowerMass_eq, div_pow]
  apply (div_le_iff₀ hc).mpr
  rw [div_mul_eq_mul_div]
  exact (le_div_iff₀ (by positivity : (0:ℝ)<(28:ℝ)^(v/2))).mpr hh

lemma half_power_le_zeta (v : ℕ) : (13/28:ℝ)^(v/2) ≤ zeta⁻¹*zeta^v := by
  have hz := zeta_pos
  have hz1 := zeta_lt_one
  rw [inv_mul_eq_div, le_div_iff₀ hz, ←zeta_sq, ←pow_mul]
  have hv : v=2*(v/2) ∨ v=2*(v/2)+1 := by omega
  rcases hv with h | h
  · conv_rhs => rw [h]
    exact mul_le_of_le_one_right (by positivity) hz1.le
  · conv_rhs => rw [h, pow_succ]

/-- The one-step tail is bounded by a geometric statistic of the input weight. -/
theorem lowerMass_le_zeta {b v D : ℕ} (hv : v ≤ b) (hDb : 125*D ≤ 13*b) :
    lowerMass b v D ≤ zeta⁻¹*zeta^v :=
  (lowerMass_le_pow hv hDb).trans (half_power_le_zeta v)

/-- Conditioning on the first weight reduces the two-accumulator lower tail to
one geometric moment, with the actual finite transition probabilities. -/
theorem double_lowerMass_le_moment (b a D : ℕ) (hDb : 125*D ≤ 13*b) :
    (∑ w ∈ range (D+1), ∑ c ∈ range (b+1), Pt b a c*Pt b c w) ≤
      zeta⁻¹*(∑ c ∈ range (b+1), Pt b a c*zeta^c) := by
  rw [Finset.sum_comm, Finset.mul_sum]
  apply Finset.sum_le_sum
  intro c hc
  rw [←Finset.mul_sum]
  have hh := mul_le_mul_of_nonneg_left (lowerMass_le_zeta (by simpa using hc) hDb)
    (show 0 ≤ Pt b a c by unfold Pt; positivity)
  simpa only [lowerMass, mul_left_comm] using hh

end Spin.Structured.ConcreteOuter

