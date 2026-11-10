import SpinCodes.Turbo.PhysicalPacketMoment22
import SpinCodes.UniformMixer.ScalableFeedbackFourier

/-!
# Exact packet character factors

The historical feedback proof uses the adjoint identity and character
orthogonality for the actual 22-bit map. Those identities remain unchanged.
Here the local input is a nonzero four-bit packet rather than independent
bits. We sum all sixteen labels and remove zero before conditioning.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer Spin.Structured.ConcreteFourier

def packetCharacter22 (v p : Packet4) : ℝ :=
  (binarySign (∑ i : Fin 4, p i * v i) : ℝ)

theorem packetCharacter22_zero (v : Packet4) : packetCharacter22 v 0 = 1 := by
  simp [packetCharacter22, binarySign_zero]

theorem packetCharacter22_product (v p : Packet4) :
    packetCharacter22 v p = ∏ i : Fin 4, (binarySign (p i * v i) : ℝ) := by
  unfold packetCharacter22
  rw [binarySign_finset_sum]
  push_cast
  rfl

theorem packetCharacter22_all_moment (v : Packet4) (z : ℝ) :
    (∑ p : Packet4, z ^ packetWeight22 p * packetCharacter22 v p) =
      (1+z)^(4-packetWeight22 v) * (1-z)^packetWeight22 v := by
  have hb (b : ZMod 2) :
      (∑ a : ZMod 2, z^(if a = 1 then 1 else 0 : ℕ) *
        (binarySign (a*b) : ℝ)) = if b = 1 then 1-z else 1+z := by
    have hsum (f : ZMod 2 → ℝ) : (∑ a, f a) = f 0 + f 1 :=
      Fin.sum_univ_two f
    rw [hsum]
    rcases zmod2_zero_or_one b with h | h
    · norm_num [h, binarySign_zero, binarySign_one]
    · norm_num [h, binarySign_zero, binarySign_one, sub_eq_add_neg]
  calc
    _ = ∑ p : Packet4, ∏ i : Fin 4,
        z^(if p i = 1 then 1 else 0 : ℕ) * (binarySign (p i * v i) : ℝ) := by
      apply sum_congr rfl
      intro p _
      rw [packetCharacter22_product, prod_mul_distrib, prod_pow_eq_pow_sum]
      rfl
    _ = ∏ i : Fin 4, if v i = 1 then 1-z else 1+z := by
      rw [← Fintype.prod_sum (fun (i : Fin 4) (a : ZMod 2) =>
        z^(if a = 1 then 1 else 0 : ℕ) * (binarySign (a*v i) : ℝ))]
      simp only [hb]
    _ = _ := by
      rw [prod_ite]
      have hw : (univ.filter (fun i : Fin 4 => v i = 1)).card = packetWeight22 v := by
        simp only [card_filter, packetWeight22]
      have hc : (univ.filter (fun i : Fin 4 => ¬ v i = 1)) =
          univ \ (univ.filter (fun i : Fin 4 => v i = 1)) := by
        ext i
        simp
      rw [hc]
      simp only [prod_const, card_sdiff_of_subset (filter_subset _ _), card_univ,
        Fintype.card_fin, hw]
      exact mul_comm _ _

theorem packetCharacter22_nonzero_moment (v : Packet4) (z : ℝ) :
    (∑ p : {p : Packet4 // p ≠ 0},
      z ^ packetWeight22 p.val * packetCharacter22 v p.val) =
        (1+z)^(4-packetWeight22 v) * (1-z)^packetWeight22 v - 1 := by
  classical
  have h := Fintype.sum_subtype_add_sum_subtype (fun p : Packet4 => p ≠ 0)
    (fun p => z ^ packetWeight22 p * packetCharacter22 v p)
  have hz : (∑ p : {p : Packet4 // ¬ p ≠ 0},
      z ^ packetWeight22 p.val * packetCharacter22 v p.val) = 1 := by
    have he (p : {p : Packet4 // ¬ p ≠ 0}) : p.val = 0 := not_not.mp p.property
    simp only [he, packetWeight22_zero, pow_zero, packetCharacter22_zero, mul_one,
      sum_const, card_univ]
    have hc : Fintype.card {p : Packet4 // ¬ p ≠ 0} = 1 := by
      simp only [not_not]
      exact Fintype.card_subtype_eq 0
    rw [hc, one_nsmul]
  rw [hz, packetCharacter22_all_moment] at h
  linarith

#print axioms packetCharacter22_all_moment
#print axioms packetCharacter22_nonzero_moment
end Spin.Turbo
