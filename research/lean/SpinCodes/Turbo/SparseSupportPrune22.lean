import SpinCodes.Turbo.SparseSupportMeasure22

/-! Eliminate inactive support coordinates by their exact zero-support
atom. This is finite product disintegration, without probabilistic
conditioning or independent resampling of messages. -/
noncomputable section
namespace Spin.Turbo
open Finset

private theorem finite_sum_pi_delta_prune
    {A B : Type*} [Fintype A] [DecidableEq A] [Fintype B] [DecidableEq B]
    (p : A → Prop) [DecidablePred p] (b0 : B) (w : A → B → ℝ) (F : (A → B) → ℝ) :
    (∑ x : A → B, (∏ i, if p i then w i (x i) else if x i = b0 then 1 else 0) * F x) =
      ∑ a : {i // p i} → B, (∏ i, w i.val (a i)) *
        F ((Equiv.piEquivPiSubtypeProd p (fun _ => B)).symm (a, fun _ => b0)) := by
  let e := Equiv.piEquivPiSubtypeProd p (fun _ => B)
  rw [← Equiv.sum_comp e.symm, Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro a _
  have hprod (b : {i // ¬p i} → B) :
      (∏ i, if p i then w i (e.symm (a, b) i)
        else if e.symm (a, b) i = b0 then 1 else 0) =
      (∏ i : {i // p i}, w i.val (a i)) * (if b = (fun _ => b0) then 1 else 0) := by
    have hd : (∏ i : {i // ¬p i}, if b i = b0 then (1 : ℝ) else 0) =
        (if b = (fun _ => b0) then 1 else 0) := by
      rw [Fintype.prod_boole]
      simp only [funext_iff]
    rw [← hd, ← Fintype.prod_subtype_mul_prod_subtype p]
    congr 1
    · apply Finset.prod_congr rfl
      intro i _
      simp [e, Equiv.piEquivPiSubtypeProd, i.property]
    · apply Finset.prod_congr rfl
      intro i _
      simp [e, Equiv.piEquivPiSubtypeProd, i.property]
  simp_rw [hprod]
  simp only [mul_ite, mul_one, mul_zero, ite_mul, zero_mul,
    Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rfl

def sparseExtendSupport22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (a : {i : Fin (R * 16) // i ∈ S} → Finset (Fin 64)) : Fin (R * 16) → Finset (Fin 64) :=
  (Equiv.piEquivPiSubtypeProd (fun i => i ∈ S) (fun _ => Finset (Fin 64))).symm
    (a, fun _ => ∅)

theorem sparseSupportProfileCount22_prune
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (F : (Fin (R * 16) → Finset (Fin 64)) → ℝ) :
    letI : Fintype {i : Fin (R * 16) // i ∈ S} :=
      Subtype.fintype (fun i : Fin (R * 16) => i ∈ S)
    (∑ supports : Fin (R * 16) → Finset (Fin 64),
      (∏ i, sparseSupportCount22 points (i ∈ S) (supports i)) * F supports) =
      ∑ a : {i : Fin (R * 16) // i ∈ S} → Finset (Fin 64),
        (∏ i : {i : Fin (R * 16) // i ∈ S}, rs8ExpectedPacketCount points (a i).card /
          (Nat.choose 64 (a i).card : ℝ)) * F (sparseExtendSupport22 S a) := by
  have h := finite_sum_pi_delta_prune (fun i => i ∈ S) ∅
    (fun _ T => rs8ExpectedPacketCount points T.card / (Nat.choose 64 T.card : ℝ)) F
  simpa only [sparseSupportCount22, sparseExtendSupport22] using h

#print axioms sparseSupportProfileCount22_prune
end Spin.Turbo



