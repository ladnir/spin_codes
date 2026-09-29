import SpinCodes.Structured.ConcreteOuterTailBlocks

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

lemma choose_le_exp_ratio (k : ℕ) {j : ℕ} (hj : 0<j) :
    (k.choose j:ℝ) ≤ (Real.exp 1*(k:ℝ)/(j:ℝ))^j := by
  have hjR : (0:ℝ)<j := by exact_mod_cast hj
  have hf : (0:ℝ)<j.factorial := by exact_mod_cast Nat.factorial_pos j
  have he := Real.pow_div_factorial_le_exp (j:ℝ) (show (0:ℝ)≤j by positivity) j
  have hmul := mul_le_mul_of_nonneg_left he (show (0:ℝ)≤((k:ℝ)/j)^j by positivity)
  have heq : ((k:ℝ)/j)^j*((j:ℝ)^j/(j.factorial:ℝ)) = (k:ℝ)^j/(j.factorial:ℝ) := by
    rw [div_pow]
    field_simp
  rw [heq] at hmul
  have heq' : ((k:ℝ)/j)^j*Real.exp (j:ℝ) = (Real.exp 1*(k:ℝ)/(j:ℝ))^j := by
    rw [show (j:ℝ)=(j:ℝ)*1 by ring, Real.exp_nat_mul, ←mul_pow]
    congr 1
    ring
  exact (Nat.choose_le_pow_div j k).trans (hmul.trans_eq heq')

/-- The precise paper summand after combining the message count and moment. -/
theorem sparse_occupation_term {k j : ℕ} (hk : 0<k) (hj : 0<j) :
    ((k.choose j:ℝ)*4095^j) *
      (4*kappa*(j:ℝ)/((k*24:ℕ)*(1-8*eta)))^(4*j) ≤
        (Cstar*((j:ℝ)/(k*24:ℕ))^3)^j := by
  have hkR : (0:ℝ)<k := by exact_mod_cast hk
  have hjR : (0:ℝ)<j := by exact_mod_cast hj
  have hz := kappa_pos
  have hden : 0<1-8*eta := by norm_num [eta]
  have hh := mul_le_mul_of_nonneg_right
    (mul_le_mul_of_nonneg_right (choose_le_exp_ratio k hj) (show (0:ℝ)≤4095^j by positivity))
    (show 0≤(4*kappa*(j:ℝ)/((k*24:ℕ)*(1-8*eta)))^(4*j) by positivity)
  refine hh.trans_eq ?_
  rw [pow_mul, ←mul_pow, ←mul_pow]
  congr 1
  unfold Cstar
  push_cast
  field_simp

/-- Summing the actual first-stage moments at one active-block count. -/
theorem golay_occupation_moment_le {k j : ℕ} (hk : 0<k) (hj : 0<j)
    (hsp : 1000*j ≤ k*24) :
    (∑ x ∈ univ.filter (fun x : LocalMessage k => golayOccupation x=j),
      ∑ h ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) h*zeta^h) ≤
        (Cstar*((j:ℝ)/(k*24:ℕ))^3)^j := by
  calc
    (∑ x ∈ univ.filter (fun x : LocalMessage k => golayOccupation x=j),
      ∑ h ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) h*zeta^h) ≤
        ∑ _x ∈ univ.filter (fun x : LocalMessage k => golayOccupation x=j),
          (4*kappa*(j:ℝ)/((k*24:ℕ)*(1-8*eta)))^(4*j) := by
      apply sum_le_sum
      intro x hx
      have he := (mem_filter.mp hx).2
      have hx0 : x≠0 := by
        intro hz
        have : golayOccupation x=0 := by rw [hz]; simp [golayOccupation,golayActive]
        omega
      simpa only [he] using golay_sparse_moment x hx0 (by omega)
    _ = ((k.choose j:ℝ)*4095^j) *
        (4*kappa*(j:ℝ)/((k*24:ℕ)*(1-8*eta)))^(4*j) := by
      rw [sum_const, nsmul_eq_mul, golay_occupation_count]
      push_cast
      rfl
    _ ≤ _ := sparse_occupation_term hk hj

end Spin.Structured.ConcreteOuter

