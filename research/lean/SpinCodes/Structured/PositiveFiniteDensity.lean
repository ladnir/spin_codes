import SpinCodes.Structured.ConcreteRoutePermutation
import SpinCodes.Structured.SparseConditioningMassLikelihood

noncomputable section
namespace Spin.Structured.PositiveFinite
open Finset ConcreteRoute ConcreteMarked

/-- Point mass of an iid bit array depends only on total Hamming weight. -/
theorem iidLaw_apply_weight {L b : ℕ} {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (x : Fin b → Finset (Fin L)) :
    (iidLaw L b p hp0 hp1).p x =
      p ^ totalWeight x * (1-p) ^ (L*b-totalWeight x) := by
  simp only [iidLaw, piPMF_apply, poissonBinom_const_apply, Finset.prod_mul_distrib,
    Finset.prod_pow_eq_pow_sum, totalWeight]
  congr 2
  rw [Finset.sum_tsub_distrib _ (fun i _ => by simpa using Finset.card_le_univ (x i))]
  simp [Nat.mul_comm]

/-- Exact iid likelihood ratio on the fixed weight preserved by every route seed. -/
theorem iidLaw_weight_ratio {L b Q : ℕ} (hQ : 0 < Q) (hQN : Q < L*b)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1)
    (x : Fin b → Finset (Fin L)) (hx : totalWeight x = Q) :
    (iidLaw L b ((Q:ℝ)/(L*b)) (by positivity)
      (by apply (div_le_one (show (0:ℝ) < (L:ℝ)*b by exact_mod_cast (lt_trans hQ hQN))).mpr; exact_mod_cast hQN.le)).p x =
    Real.exp (((L:ℝ)*b)*binKL ((Q:ℝ)/(L*b)) p) *
      (iidLaw L b p hp0.le hp1.le).p x := by
  have hN : 0 < L*b := lt_trans hQ hQN
  have hn : (0:ℝ) < (L*b:ℕ) := by exact_mod_cast hN
  have ha : 0 < (Q:ℝ)/(L*b:ℕ) := by positivity
  have ha1 : (Q:ℝ)/(L*b:ℕ) < 1 := (div_lt_one hn).mpr (by exact_mod_cast hQN)
  have hc : (0:ℝ) < ((L*b).choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQN.le
  have hr := markMass_likelihood_ratio hQ hQN.le hp0 hp1
  have he : Real.exp ((L*b:ℕ)*binKL ((Q:ℝ)/(L*b:ℕ)) p) *
      Real.exp (-(L*b:ℕ)*binKL ((Q:ℝ)/(L*b:ℕ)) p) = 1 := by
    rw [← Real.exp_add]
    ring_nf
    exact Real.exp_zero
  have hr' : markMass (L*b) Q ((Q:ℝ)/(L*b:ℕ)) =
      Real.exp ((L*b:ℕ)*binKL ((Q:ℝ)/(L*b:ℕ)) p) * markMass (L*b) Q p := by
    rw [hr]
    nlinarith [congrArg (fun t => markMass (L*b) Q ((Q:ℝ)/(L*b:ℕ)) * t) he]
  rw [iidLaw_apply_weight, iidLaw_apply_weight, hx]
  unfold markMass at hr'
  have hh : ((L*b).choose Q:ℝ) *
      (((Q:ℝ)/(L*b:ℕ))^Q * (1-(Q:ℝ)/(L*b:ℕ))^(L*b-Q)) =
      ((L*b).choose Q:ℝ) * (Real.exp ((L*b:ℕ)*binKL ((Q:ℝ)/(L*b:ℕ)) p) *
      (p^Q*(1-p)^(L*b-Q))) := by nlinarith [hr']
  have hhh := (mul_left_cancel₀ hc.ne') hh
  simpa only [Nat.cast_mul] using hhh

end Spin.Structured.PositiveFinite
