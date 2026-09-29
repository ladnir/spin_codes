import SpinCodes.Structured.ConcreteMarkedBits

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteScalar
open Finset ConcreteRoute
open scoped symmDiff

def flipProb {n : ℕ} (p : Fin n → ℝ) (Y : Finset (Fin n)) (i : Fin n) : ℝ :=
  if i ∈ Y then 1-p i else p i

theorem flipProb_nonneg {n : ℕ} (p : Fin n → ℝ) (hp0 : ∀ i, 0 ≤ p i)
    (hp1 : ∀ i, p i ≤ 1) (Y : Finset (Fin n)) (i : Fin n) : 0 ≤ flipProb p Y i := by
  unfold flipProb
  split_ifs <;> linarith [hp0 i, hp1 i]

theorem flipProb_le_one {n : ℕ} (p : Fin n → ℝ) (hp0 : ∀ i, 0 ≤ p i)
    (hp1 : ∀ i, p i ≤ 1) (Y : Finset (Fin n)) (i : Fin n) : flipProb p Y i ≤ 1 := by
  unfold flipProb
  split_ifs <;> linarith [hp0 i, hp1 i]

/-- XOR with a fixed word flips the corresponding independent Bernoulli parameters. -/
theorem poissonBinom_xor {n : ℕ} (p : Fin n → ℝ) (hp0 : ∀ i, 0 ≤ p i)
    (hp1 : ∀ i, p i ≤ 1) (Y : Finset (Fin n)) :
    (poissonBinom p hp0 hp1).map (fun X => X ∆ Y) =
      poissonBinom (flipProb p Y) (flipProb_nonneg p hp0 hp1 Y) (flipProb_le_one p hp0 hp1 Y) := by
  let e : Equiv.Perm (Finset (Fin n)) := Function.Involutive.toPerm
    (fun X => X ∆ Y) (fun X => by simp)
  ext X
  change ((poissonBinom p hp0 hp1).map e).p X = _
  rw [FinPMF.map_equiv_apply]
  change (poissonBinom p hp0 hp1).p (X ∆ Y) = _
  simp only [poissonBinom_eq_prod]
  apply Finset.prod_congr rfl
  intro i _
  by_cases hiX : i ∈ X <;> by_cases hiY : i ∈ Y <;>
    simp [flipProb, hiX, hiY, Finset.mem_symmDiff] <;> ring

/-- Exact one-round scalar generating function, uniform over the entering state. -/
theorem bernoulli_xor_moment {n : ℕ} (β z : ℝ) (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (Y : Finset (Fin n)) :
    (poissonBinom (fun _ : Fin n => β) (fun _ => hβ0) (fun _ => hβ1)).expect
      (fun X => z ^ (X ∆ Y).card) =
      (1-β+β*z)^(n-Y.card) * (β+(1-β)*z)^Y.card := by
  rw [← FinPMF.expect_map _ (fun X => X ∆ Y) (fun X => z^X.card), poissonBinom_xor]
  change (∑ X : Finset (Fin n), _ * z^X.card) = _
  rw [poissonBinom_mgf]
  have he : (∏ i : Fin n, (1-flipProb (fun _ => β) Y i+flipProb (fun _ => β) Y i*z)) =
      ∏ i : Fin n, if i ∈ Y then β+(1-β)*z else 1-β+β*z := by
    apply prod_congr rfl
    intro i _
    by_cases hi : i ∈ Y <;> simp [flipProb, hi] <;> ring
  rw [he, Finset.prod_ite]
  have hyes : (univ : Finset (Fin n)).filter (fun i => i ∈ Y) = Y := by ext; simp
  have hno : (univ : Finset (Fin n)).filter (fun i => i ∉ Y) = univ \ Y := by ext; simp
  rw [hyes, hno]
  simp only [prod_const, card_sdiff, inter_univ, card_univ, Fintype.card_fin]
  ring

end Spin.Structured.ConcreteScalar



