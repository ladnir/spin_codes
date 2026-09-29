import SpinCodes.Structured.SparseConditioningMassMoments

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset

def centralSet (N Q r : ℕ) : Finset (Fin N) :=
  univ.filter (fun j => Q < j.val+r ∧ j.val < Q+r)

theorem centralSet_card_le (N Q r : ℕ) (hr : 0 < r) :
    (centralSet N Q r).card ≤ 2*r-1 := by
  have hi : (centralSet N Q r).card ≤ (Icc (Q+1-r) (Q+r-1)).card := by
    apply Finset.card_le_card_of_injOn Fin.val
    · intro j hj
      change j ∈ centralSet N Q r at hj
      simp only [centralSet, mem_filter, mem_univ, true_and] at hj
      change j.val ∈ (Icc (Q+1-r) (Q+r-1) : Finset ℕ)
      exact Finset.mem_Icc.mpr (by omega)
    · intro x _ y _ hxy
      exact Fin.ext hxy
  rw [Nat.card_Icc] at hi
  omega

theorem outside_central_square {N Q r : ℕ} (j : Fin N)
    (hj : j ∉ centralSet N Q r) : (r:ℝ)^2 ≤ ((j.val:ℝ)-Q)^2 := by
  have hh : j.val+r ≤ Q ∨ Q+r ≤ j.val := by
    simp only [centralSet, mem_filter, mem_univ, true_and] at hj
    omega
  have habs : (r:ℝ) ≤ |(j.val:ℝ)-Q| := by
    rcases hh with hh | hh
    · have hc : (j.val:ℝ)+(r:ℝ) ≤ Q := by exact_mod_cast hh
      exact le_trans (by linarith) (neg_le_abs _)
    · have hc : (Q:ℝ)+(r:ℝ) ≤ j.val := by exact_mod_cast hh
      exact le_trans (by linarith) (le_abs_self _)
  simpa only [sq_abs] using (sq_le_sq₀ (Nat.cast_nonneg r) (abs_nonneg _)).mpr habs

/-- A mode with second moment at most Q about integer Q has mass at least 1/(8 sqrt Q). -/
theorem mode_mass_of_variance {N Q : ℕ} (hQ : 0 < Q)
    (P : FinPMF (Fin N)) (m : Fin N) (hmode : ∀ j, P.p j ≤ P.p m)
    (hv : P.expect (fun j => ((j.val:ℝ)-Q)^2) ≤ Q) :
    1 / (8 * Real.sqrt Q) ≤ P.p m := by
  classical
  let r := Nat.ceil (2 * Real.sqrt (Q:ℝ))
  have hQp : (0:ℝ) < Q := by exact_mod_cast hQ
  have hQ1 : (1:ℝ) ≤ Q := by exact_mod_cast hQ
  have hs : 0 < Real.sqrt (Q:ℝ) := Real.sqrt_pos.2 hQp
  have hs1 : 1 ≤ Real.sqrt (Q:ℝ) := by
    nlinarith [Real.sq_sqrt hQp.le, Real.sqrt_nonneg (Q:ℝ)]
  have hr : 0 < r := Nat.ceil_pos.mpr (by positivity)
  have hrlo : 2 * Real.sqrt (Q:ℝ) ≤ (r:ℝ) := Nat.le_ceil _
  have hrhi : (r:ℝ) < 2 * Real.sqrt (Q:ℝ)+1 := Nat.ceil_lt_add_one (by positivity)
  have hrsq : 4*(Q:ℝ) ≤ (r:ℝ)^2 := by
    have hh := (sq_le_sq₀ (by positivity : 0 ≤ 2*Real.sqrt (Q:ℝ))
      (Nat.cast_nonneg r)).mpr hrlo
    nlinarith [Real.sq_sqrt hQp.le]
  have hcard : ((centralSet N Q r).card:ℝ) ≤ 5*Real.sqrt (Q:ℝ) := by
    have hh : ((centralSet N Q r).card:ℝ) ≤ ((2*r-1:ℕ):ℝ) := by
      exact_mod_cast centralSet_card_le N Q r hr
    rw [Nat.cast_sub (by omega), Nat.cast_mul, Nat.cast_ofNat, Nat.cast_one] at hh
    linarith
  have hout : P.prob (fun j => j ∉ centralSet N Q r) ≤ 1/4 := by
    have hm := P.markov (fun j => sq_nonneg ((j.val:ℝ)-Q))
      (by positivity : 0 < 4*(Q:ℝ))
    have hle := P.prob_mono (s := fun j => j ∉ centralSet N Q r)
      (t := fun j => 4*(Q:ℝ) ≤ ((j.val:ℝ)-Q)^2)
      (fun j hj => hrsq.trans (outside_central_square j hj))
    refine hle.trans (hm.trans ?_)
    apply (div_le_iff₀ (by positivity : 0 < 4*(Q:ℝ))).mpr
    linarith
  have hin : 3/4 ≤ P.prob (fun j => j ∈ centralSet N Q r) := by
    have hc := P.prob_compl (fun j => j ∈ centralSet N Q r)
    linarith
  have hmodeSum : P.prob (fun j => j ∈ centralSet N Q r) ≤
      ((centralSet N Q r).card:ℝ) * P.p m := by
    unfold FinPMF.prob
    have he : univ.filter (fun j => j ∈ centralSet N Q r) = centralSet N Q r := by ext j; simp
    rw [he]
    calc ∑ j ∈ centralSet N Q r, P.p j ≤ ∑ _j ∈ centralSet N Q r, P.p m :=
        sum_le_sum (fun j _ => hmode j)
      _ = _ := by simp
  have hh := hin.trans (hmodeSum.trans (mul_le_mul_of_nonneg_right hcard (P.nonneg m)))
  apply (div_le_iff₀ (by positivity : 0 < 8*Real.sqrt (Q:ℝ))).mpr
  nlinarith

theorem markMass_mode_lower {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L) :
    1 / (8 * Real.sqrt (Q:ℝ)) ≤ markMass L Q ((Q:ℝ)/L) := by
  have hL : 0 < L := lt_of_lt_of_le hQ hQL
  have h0 : (0:ℝ) ≤ (Q:ℝ)/L := by positivity
  have h1 : (Q:ℝ)/L ≤ 1 :=
    div_le_one_of_le₀ (by exact_mod_cast hQL) (by positivity)
  exact mode_mass_of_variance hQ (binPMF L ((Q:ℝ)/L) h0 h1) ⟨Q, by omega⟩
    (binPMF_mode hL hQL h0 h1) (binPMF_mode_variance hL hQL h0 h1)

end Spin.Structured.ConcreteMarked


