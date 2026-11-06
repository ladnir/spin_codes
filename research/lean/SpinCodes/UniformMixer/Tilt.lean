import SpinCodes.UniformMixer.Map22
import SpinCodes.UniformMixer.Regional
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

/-!
# Elementary finite-state tilt estimates

The quadratic bound below controls the error from replacing the average
empty-step tilt by the exponential of its exact live-state mean. It avoids
any stationary-state assumption about the state entering a region.
-/

noncomputable section

namespace Spin.UniformMixer

/-- A quadratic sandwich for the exponential tilt. -/
theorem exp_neg_quadratic (x : ℝ) (hx : 0 ≤ x) :
    1 - x ≤ Real.exp (-x) ∧ Real.exp (-x) ≤ 1 - x + x ^ 2 := by
  constructor
  · have h := Real.add_one_le_exp (-x)
    linarith
  · have hpos : 0 < 1 + x := by linarith
    have hexp : 1 + x ≤ Real.exp x := by
      simpa [add_comm] using Real.add_one_le_exp x
    have hmul : (1 + x) * Real.exp (-x) ≤ 1 := by
      calc
        (1 + x) * Real.exp (-x) ≤ Real.exp x * Real.exp (-x) :=
          mul_le_mul_of_nonneg_right hexp (Real.exp_pos _).le
        _ = 1 := by rw [← Real.exp_add]; simp
    have hpoly : 1 ≤ (1 + x) * (1 - x + x ^ 2) := by
      nlinarith [mul_nonneg hx (sq_nonneg x)]
    nlinarith

/-- Finite-sum form: averaging an exponential tilt differs from its
first-order approximation by at most the square of the largest exponent. -/
theorem finite_tilt_quadratic {α : Type*} (s : Finset α)
    (w : α → ℕ) (t : ℕ) (c : ℝ)
    (hc : 0 ≤ c) (hw : ∀ v ∈ s, w v ≤ t) :
    (s.sum fun v => (1 : ℝ) - c * (w v : ℝ)) ≤
        s.sum (fun v => Real.exp (-(c * (w v : ℝ)))) ∧
    (s.sum fun v => Real.exp (-(c * (w v : ℝ)))) ≤
        (s.sum fun v => (1 : ℝ) - c * (w v : ℝ)) +
          (s.card : ℝ) * (c * (t : ℝ)) ^ 2 := by
  have hpoint (v : α) (hv : v ∈ s) :
      (1 : ℝ) - c * (w v : ℝ) ≤ Real.exp (-(c * (w v : ℝ))) ∧
      Real.exp (-(c * (w v : ℝ))) ≤
        (1 : ℝ) - c * (w v : ℝ) + (c * (t : ℝ)) ^ 2 := by
    have hcv : 0 ≤ c * (w v : ℝ) := mul_nonneg hc (Nat.cast_nonneg _)
    have hct : c * (w v : ℝ) ≤ c * (t : ℝ) :=
      mul_le_mul_of_nonneg_left (by exact_mod_cast hw v hv) hc
    obtain ⟨hlo, hhi⟩ := exp_neg_quadratic _ hcv
    constructor
    · exact hlo
    · have hsq : (c * (w v : ℝ)) ^ 2 ≤ (c * (t : ℝ)) ^ 2 := by
        nlinarith [sq_nonneg (c * (t : ℝ) - c * (w v : ℝ))]
      linarith
  constructor
  · apply Finset.sum_le_sum
    intro v hv
    exact (hpoint v hv).1
  · calc
      (s.sum fun v => Real.exp (-(c * (w v : ℝ)))) ≤
          s.sum (fun v => (1 : ℝ) - c * (w v : ℝ) + (c * (t : ℝ)) ^ 2) := by
            apply Finset.sum_le_sum
            intro v hv
            exact (hpoint v hv).2
      _ = (s.sum fun v => (1 : ℝ) - c * (w v : ℝ)) +
            (s.card : ℝ) * (c * (t : ℝ)) ^ 2 := by
          rw [Finset.sum_add_distrib]
          simp

/-- The empty-step tilt for a binary linear output map is within a quadratic
error of the exponential's first-order approximation at its exact mean. -/
theorem liveMean_quadratic {V : Type*} [AddCommGroup V]
    [Module (ZMod 2) V] [Fintype V] [DecidableEq V] [Nontrivial V]
    {t : ℕ} (A : Fin t → V →ₗ[ZMod 2] ZMod 2) (c : ℝ) (hc : 0 ≤ c) :
    1 - c * liveWeightAverage A ≤
      liveMean (Real.exp (-c)) (coordinateWeight A) ∧
    liveMean (Real.exp (-c)) (coordinateWeight A) ≤
      1 - c * liveWeightAverage A + (c * (t : ℝ)) ^ 2 := by
  let s : Finset V := (Finset.univ : Finset V).erase 0
  have hs : s.card = liveCount V := by
    simp [s, liveCount]
  have hd : (0 : ℝ) < (liveCount V : ℝ) := by
    exact_mod_cast liveCount_pos (S := V)
  have hexp (v : V) :
      Real.exp (-c) ^ coordinateWeight A v =
        Real.exp (-(c * (coordinateWeight A v : ℝ))) := by
    rw [← Real.exp_nat_mul]
    congr 1
    ring
  have hraw := finite_tilt_quadratic s (coordinateWeight A) t c hc
    (fun v _ => coordinateWeight_le A v)
  have hlinear :
      (s.sum fun v => (1 : ℝ) - c * (coordinateWeight A v : ℝ)) /
          (liveCount V : ℝ) = 1 - c * liveWeightAverage A := by
    rw [Finset.sum_sub_distrib, liveWeightAverage_erase]
    simp only [Finset.sum_const, nsmul_eq_mul, mul_one]
    rw [hs]
    field_simp
    rw [Finset.mul_sum]
  have htilt :
      (s.sum fun v => Real.exp (-(c * (coordinateWeight A v : ℝ)))) /
          (liveCount V : ℝ) =
        liveMean (Real.exp (-c)) (coordinateWeight A) := by
    simp only [liveMean, s, ← hexp]
  constructor
  · rw [← hlinear, ← htilt]
    exact div_le_div_of_nonneg_right hraw.1 hd.le
  · rw [← hlinear, ← htilt]
    have h := div_le_div_of_nonneg_right hraw.2 hd.le
    rw [add_div, hs] at h
    have hden : (liveCount V : ℝ) ≠ 0 := ne_of_gt hd
    simpa [hden] using h

/-- Exact production mean and a second-order bound on the scalar
empty-step tilt. This is the analytic input to the regional comparison. -/
theorem production_liveMean_close (c : ℝ) (hc : 0 ≤ c) :
    |liveMean (Real.exp (-c)) (coordinateWeight coordMap) -
      Real.exp (-(c * liveWeightAverage coordMap))| ≤
        (c * 64) ^ 2 := by
  have hmean := liveMean_quadratic coordMap c hc
  have hμ0 : 0 ≤ liveWeightAverage coordMap := by
    rw [production_liveWeightAverage]
    norm_num
  have hμt : liveWeightAverage coordMap ≤ 64 := by
    rw [production_liveWeightAverage]
    norm_num
  have hquad := exp_neg_quadratic
    (c * liveWeightAverage coordMap) (mul_nonneg hc hμ0)
  have hsq : (c * liveWeightAverage coordMap) ^ 2 ≤ (c * 64) ^ 2 := by
    nlinarith [sq_nonneg (c * (64 - liveWeightAverage coordMap))]
  rcases hmean with ⟨hmeanlo, hmeanhi⟩
  rcases hquad with ⟨hquadlo, hquadhi⟩
  norm_num only [Nat.cast_ofNat] at hmeanhi
  rw [abs_le]
  constructor <;> nlinarith [hmeanlo, hmeanhi, hquadlo, hquadhi, hsq]

theorem liveMean_bounds {S : Type*} [Fintype S] [DecidableEq S]
    [AddCommGroup S] [Nontrivial S] (z : ℝ) (w : S → ℕ)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    0 ≤ liveMean z w ∧ liveMean z w ≤ 1 := by
  let s : Finset S := (Finset.univ : Finset S).erase 0
  have hs : s.card = liveCount S := by simp [s, liveCount]
  have hd : (0 : ℝ) < (liveCount S : ℝ) := by
    exact_mod_cast liveCount_pos (S := S)
  have hsum0 : 0 ≤ s.sum (fun v => z ^ w v) := by
    apply Finset.sum_nonneg
    intro v _
    exact pow_nonneg hz0 _
  have hsum1 : s.sum (fun v => z ^ w v) ≤ (liveCount S : ℝ) := by
    calc
      s.sum (fun v => z ^ w v) ≤ s.sum (fun _ => (1 : ℝ)) := by
        apply Finset.sum_le_sum
        intro v _
        exact pow_le_one₀ hz0 hz1
      _ = (liveCount S : ℝ) := by simp [hs]
  constructor
  · exact div_nonneg hsum0 hd.le
  · exact (div_le_one hd).mpr hsum1

/-- Error through an empty gap, with a bound linear in its length. -/
theorem production_liveMean_power_close (c : ℝ) (hc : 0 ≤ c) (g : ℕ) :
    |(liveMean (Real.exp (-c)) (coordinateWeight coordMap)) ^ g -
      (Real.exp (-(c * liveWeightAverage coordMap))) ^ g| ≤
        (g : ℝ) * (c * 64) ^ 2 := by
  let a := liveMean (Real.exp (-c)) (coordinateWeight coordMap)
  let b := Real.exp (-(c * liveWeightAverage coordMap))
  have ha := liveMean_bounds (Real.exp (-c))
    (coordinateWeight coordMap) (Real.exp_pos _).le
    (Real.exp_le_one_iff.mpr (by linarith))
  have hb0 : 0 ≤ b := (Real.exp_pos _).le
  have hb1 : b ≤ 1 := by
    apply Real.exp_le_one_iff.mpr
    have hμ : 0 ≤ liveWeightAverage coordMap := by
      rw [production_liveWeightAverage]
      norm_num
    nlinarith
  have hmax : max |a| |b| ≤ 1 := by
    rw [abs_of_nonneg ha.1, abs_of_nonneg hb0]
    exact max_le ha.2 hb1
  have hpow : max |a| |b| ^ (g - 1) ≤ 1 :=
    pow_le_one₀ (le_max_of_le_left (abs_nonneg a)) hmax
  have hdiff : |a - b| ≤ (c * 64) ^ 2 :=
    production_liveMean_close c hc
  calc
    |a ^ g - b ^ g| ≤ |a - b| * (g : ℝ) * max |a| |b| ^ (g - 1) :=
      abs_pow_sub_pow_le a b g
    _ ≤ (c * 64) ^ 2 * (g : ℝ) * 1 := by
      gcongr
    _ = (g : ℝ) * (c * 64) ^ 2 := by ring

end Spin.UniformMixer
