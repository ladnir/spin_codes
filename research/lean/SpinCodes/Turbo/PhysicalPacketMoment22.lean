import SpinCodes.Turbo.PhysicalGenericLive22
import Mathlib.Algebra.BigOperators.Ring.Finset

/-!
# Four-bit packet moment factors

The historical one-bit regional proof averages a fair bit explicitly.
Turbo's physical input uses a nonzero four-bit packet. Sum all sixteen
labels by the coordinate product identity, then remove the zero label.
This derives the local factor used by the compact profile calculation;
there is no assumption that the four bits of a nonzero packet are independent.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

instance packet4NonzeroNonempty22 : Nonempty {p : Packet4 // p ≠ 0} :=
  ⟨⟨1, by
    intro h
    have hb := congrFun h 0
    norm_num at hb⟩⟩

def packetWeight22 (p : Packet4) : ℕ :=
  ∑ i : Fin 4, if p i = 1 then 1 else 0

theorem packetWeight22_zero : packetWeight22 0 = 0 := by
  simp [packetWeight22]

theorem packetWeight22_all_moment (z : ℝ) :
    (∑ p : Packet4, z ^ packetWeight22 p) = (1 + z) ^ 4 := by
  have hbit : (∑ b : ZMod 2, z ^ (if b = 1 then 1 else 0 : ℕ)) = 1 + z := by
    change (∑ b : Fin 2, z ^ (if b = 1 then 1 else 0 : ℕ)) = 1 + z
    simp [Fin.sum_univ_two]
  calc
    _ = ∑ p : Packet4, ∏ i : Fin 4,
        z ^ (if p i = 1 then 1 else 0 : ℕ) := by
      apply sum_congr rfl
      intro p _
      exact (Finset.prod_pow_eq_pow_sum univ
        (fun i : Fin 4 => if p i = 1 then 1 else 0) z).symm
    _ = ∏ _i : Fin 4, ∑ b : ZMod 2,
        z ^ (if b = 1 then 1 else 0 : ℕ) :=
      (Fintype.prod_sum (fun (_i : Fin 4) (b : ZMod 2) =>
        z ^ (if b = 1 then 1 else 0 : ℕ))).symm
    _ = _ := by simp only [hbit, Finset.prod_const, Finset.card_univ, Fintype.card_fin]

theorem packetWeight22_shift_moment (v : Packet4) (z : ℝ) :
    (∑ p : Packet4, z ^ packetWeight22 (v + p)) = (1 + z) ^ 4 := by
  calc
    _ = ∑ p : Packet4, z ^ packetWeight22 p :=
      Equiv.sum_comp (Equiv.addLeft v) (fun p => z ^ packetWeight22 p)
    _ = _ := packetWeight22_all_moment z

theorem packetWeight22_nonzero_shift_sum (v : Packet4) (z : ℝ) :
    (∑ p : {p : Packet4 // p ≠ 0}, z ^ packetWeight22 (v + p.val)) =
      (1 + z) ^ 4 - z ^ packetWeight22 v := by
  classical
  have hsum := Fintype.sum_subtype_add_sum_subtype (fun p : Packet4 => p ≠ 0)
    (fun p : Packet4 => z ^ packetWeight22 (v + p))
  have hzero : (∑ p : {p : Packet4 // ¬ p ≠ 0},
      z ^ packetWeight22 (v + p.val)) = z ^ packetWeight22 v := by
    have he (p : {p : Packet4 // ¬ p ≠ 0}) : p.val = 0 := not_not.mp p.property
    simp only [he, add_zero, Finset.sum_const, Finset.card_univ]
    have hc : Fintype.card {p : Packet4 // ¬ p ≠ 0} = 1 := by
      simp only [not_not]
      exact Fintype.card_subtype_eq 0
    rw [hc, one_nsmul]
  rw [hzero, packetWeight22_shift_moment] at hsum
  linarith

/-- The exact expectation over fifteen nonzero labels, not four independent
conditioned bits. -/
theorem packetWeight22_nonzero_shift_mean (v : Packet4) (z : ℝ) :
    (FinPMF.uniform {p : Packet4 // p ≠ 0}).expect
      (fun p => z ^ packetWeight22 (v + p.val)) =
        ((1 + z) ^ 4 - z ^ packetWeight22 v) / 15 := by
  simp only [FinPMF.expect, FinPMF.uniform, packet4_nonzero_card]
  rw [← Finset.mul_sum, packetWeight22_nonzero_shift_sum]
  ring

#print axioms packetWeight22_all_moment
#print axioms packetWeight22_nonzero_shift_mean
end Spin.Turbo
