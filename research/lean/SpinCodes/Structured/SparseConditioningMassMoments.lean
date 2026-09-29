import SpinCodes.Structured.ConcreteMarkedConditioning
import SpinCodes.Structured.ConcreteEncoder

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset

def ones {n : ℕ} (x : Fin n → Bool) : ℕ := (support x).card

theorem ones_eq_sum {n : ℕ} (x : Fin n → Bool) :
    ones x = ∑ i, if x i then 1 else 0 := by
  exact Finset.card_filter _ _

theorem ones_cons {n : ℕ} (bit : Bool) (x : Fin n → Bool) :
    ones (Fin.cons bit x) = (if bit then 1 else 0) + ones x := by
  simp only [ones_eq_sum, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ]

theorem expect_affine_square {Ω : Type*} [Fintype Ω] (P : FinPMF Ω)
    (X : Ω → ℝ) (a b c : ℝ) :
    P.expect (fun x => a + b*X x + c*(X x)^2) =
      a + b*P.expect X + c*P.expect (fun x => (X x)^2) := by
  unfold FinPMF.expect
  simp only [mul_add, Finset.sum_add_distrib]
  rw [← Finset.sum_mul, P.total, one_mul]
  congr 1
  · congr 1
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro x _
    ring
  · congr 1
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro x _
    ring

theorem coin_ones_moments (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect (fun x => (ones x : ℝ)) = n*p ∧
    (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect (fun x => (ones x : ℝ)^2) =
      n*p*(1-p)+(n*p)^2 := by
  induction n with
  | zero => simp [ones, support, FinPMF.expect_const]
  | succ n ih =>
    have ha : (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
        (fun x => 1+(ones x : ℝ)) = 1+n*p := by
      have hh := expect_affine_square (piPMF (fun _ : Fin n => coin p hp0 hp1))
        (fun x => (ones x : ℝ)) 1 1 0
      simpa only [one_mul, zero_mul, add_zero, ih.1] using hh
    have hb : (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
        (fun x => (1+(ones x : ℝ))^2) = 1+2*(n*p)+(n*p*(1-p)+(n*p)^2) := by
      have hh := expect_affine_square (piPMF (fun _ : Fin n => coin p hp0 hp1))
        (fun x => (ones x : ℝ)) 1 2 1
      have he : (fun x : Fin n → Bool => (1+(ones x : ℝ))^2) =
          (fun x => 1+2*(ones x : ℝ)+1*(ones x : ℝ)^2) := by funext x; ring
      rw [he, hh, ih.1, ih.2]
      ring
    constructor
    · rw [ConcreteEncoder.expect_iid_succ]
      simp only [FinPMF.expect, coin, Fintype.sum_bool, ones_cons, Bool.false_eq_true,
        ite_false, Bool.true_eq, ite_true, Nat.cast_zero, zero_add, Nat.cast_add, Nat.cast_one]
      conv_lhs => rw [add_comm]
      change (1-p) * (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
        (fun x => (ones x : ℝ)) + p *
        (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
          (fun x => 1+(ones x : ℝ)) = _
      rw [ih.1, ha]
      ring
    · rw [ConcreteEncoder.expect_iid_succ]
      simp only [FinPMF.expect, coin, Fintype.sum_bool, ones_cons, Bool.false_eq_true,
        ite_false, Bool.true_eq, ite_true, Nat.cast_zero, zero_add, Nat.cast_add, Nat.cast_one]
      conv_lhs => rw [add_comm]
      change (1-p) * (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
        (fun x => (ones x : ℝ)^2) + p *
        (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
          (fun x => (1+(ones x : ℝ))^2) = _
      rw [ih.2, hb]
      ring

def cardIndex {n : ℕ} (S : Finset (Fin n)) : Fin (n+1) :=
  ⟨S.card, Nat.lt_succ_of_le (Routing.card_le_width S)⟩

theorem count_law (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (iidBits n p hp0 hp1).map cardIndex = binPMF n p hp0 hp1 := by
  classical
  ext k
  have he : ((iidBits n p hp0 hp1).map cardIndex).p k =
      ((iidBits n p hp0 hp1).map cardIndex).prob (fun j => j = k) := by
    simp [FinPMF.prob, Finset.sum_filter]
  rw [he, FinPMF.map_prob]
  calc (iidBits n p hp0 hp1).prob (fun S => cardIndex S = k)
      = (iidBits n p hp0 hp1).prob (fun S => S.card = k.val) :=
        FinPMF.prob_congr _ (fun S => by simp [cardIndex, Fin.ext_iff])
    _ = _ := iidBits_prob_card n k.val p hp0 hp1
theorem binPMF_expect (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (F : Fin (n+1) → ℝ) :
    (binPMF n p hp0 hp1).expect F =
      (piPMF (fun _ : Fin n => coin p hp0 hp1)).expect
        (fun x => F (cardIndex (support x))) := by
  rw [← count_law, FinPMF.expect_map, ← support_coin, FinPMF.expect_map]

theorem binPMF_moments (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (binPMF n p hp0 hp1).expect (fun x => (x.val : ℝ)) = n*p ∧
    (binPMF n p hp0 hp1).expect (fun x => (x.val : ℝ)^2) =
      n*p*(1-p)+(n*p)^2 := by
  simp only [binPMF_expect]
  exact coin_ones_moments n p hp0 hp1

theorem binPMF_variance (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (binPMF n p hp0 hp1).expect (fun x => ((x.val : ℝ)-n*p)^2) = n*p*(1-p) := by
  have hm := binPMF_moments n p hp0 hp1
  have hh := expect_affine_square (binPMF n p hp0 hp1)
    (fun x => (x.val : ℝ)) ((n*p)^2) (-(2*n*p)) 1
  have he : (fun x : Fin (n+1) => ((x.val : ℝ)-n*p)^2) =
      (fun x => (n*p)^2 + (-(2*n*p))*(x.val:ℝ) + 1*(x.val:ℝ)^2) := by
    funext x
    ring
  rw [he, hh, hm.1, hm.2]
  ring

theorem binPMF_mode_variance {L Q : ℕ} (hL : 0 < L) (hQL : Q ≤ L)
    (h0 : (0:ℝ) ≤ (Q:ℝ)/L) (h1 : (Q:ℝ)/L ≤ 1) :
    (binPMF L ((Q:ℝ)/L) h0 h1).expect (fun j => ((j.val:ℝ)-Q)^2) ≤ Q := by
  have hLp : (L:ℝ) ≠ 0 := by exact_mod_cast hL.ne'
  have hn : (L:ℝ)*((Q:ℝ)/L) = Q := by field_simp
  have hv := binPMF_variance L ((Q:ℝ)/L) h0 h1
  rw [hn] at hv
  rw [hv]
  nlinarith [show (0:ℝ) ≤ (Q:ℝ) from Nat.cast_nonneg Q]

end Spin.Structured.ConcreteMarked



