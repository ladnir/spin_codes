import SpinCodes.Turbo.RoutingConditionalKernel22
import SpinCodes.Turbo.MatrixCoefficient22

/-!
# Actual regional count profiles and ordered paired coefficients

The checked sixteen-slot profile law is reindexed into chronological pairs.
Its binomial weights remain outside the fractional power. The ordered
coefficient recurrence is exactly `MatrixCoefficient22`, with one regional
normalizer; no independent-slot approximation is introduced.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset

/-- Consecutive positions 2b,2b+1 form one chronological pair. -/
def pairedCountEquiv22 (n : ℕ) :
    (Fin (n * 2) → Fin 17) ≃ (Fin n → Fin 17 × Fin 17) :=
  (Equiv.arrowCongr (finProdFinEquiv : Fin n × Fin 2 ≃ Fin (n * 2)).symm
    (Equiv.refl (Fin 17))).trans
      ((Equiv.curry (Fin n) (Fin 2) (Fin 17)).trans
        (Equiv.piCongrRight (fun _ : Fin n => finTwoArrowEquiv (Fin 17))))

theorem pairedCountEquiv22_degree (n : ℕ) (counts : Fin (n * 2) → Fin 17) :
    (∑ step, (counts step).val) =
      ∑ b, pairedProfileDegree22 (pairedCountEquiv22 n counts b) := by
  rw [← Equiv.sum_comp (finProdFinEquiv : Fin n × Fin 2 ≃ Fin (n * 2)),
    Fintype.sum_prod_type]
  simp only [Fin.sum_univ_two, pairedProfileDegree22, pairedCountEquiv22,
    Equiv.trans_apply, Equiv.arrowCongr_apply, Equiv.symm_symm,
    Equiv.refl_apply, Equiv.curry_apply, Equiv.piCongrRight_apply,
    finTwoArrowEquiv, piFinTwoEquiv, Equiv.coe_fn_mk]
  rfl

theorem pairedCountEquiv22_choose_product (n : ℕ) (counts : Fin (n * 2) → Fin 17) :
    (∏ step, (Nat.choose 16 (counts step).val : ℝ)) =
      ∏ b, ((Nat.choose 16 (pairedCountEquiv22 n counts b).1.val : ℝ) *
        (Nat.choose 16 (pairedCountEquiv22 n counts b).2.val : ℝ)) := by
  rw [← Equiv.prod_comp (finProdFinEquiv : Fin n × Fin 2 ≃ Fin (n * 2)),
    Fintype.prod_prod_type]
  simp only [Fin.prod_univ_two, pairedCountEquiv22,
    Equiv.trans_apply, Equiv.arrowCongr_apply, Equiv.symm_symm,
    Equiv.refl_apply, Equiv.curry_apply, Equiv.piCongrRight_apply,
    finTwoArrowEquiv, piFinTwoEquiv, Equiv.coe_fn_mk]
  rfl

private theorem matrix_ofFn_smul_prod {I : Type*} [Fintype I] [DecidableEq I]
    {n : ℕ} (c : Fin n → ℝ) (K : Fin n → Matrix I I ℝ) :
    (List.ofFn (fun i => c i • K i)).prod = (∏ i, c i) • (List.ofFn K).prod := by
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [List.ofFn_succ, List.prod_cons, Fin.prod_univ_succ]
    rw [ih, smul_mul_smul_comm]

def regionalPairedMatrix22 (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (α : ℝ) (n : ℕ) (counts : Fin (n * 2) → Fin 17) : Matrix (Fin 11) (Fin 11) ℝ :=
  (List.ofFn (fun b : Fin n => matrixEntrywiseRpow
    (B (pairedCountEquiv22 n counts b).1 * B (pairedCountEquiv22 n counts b).2) α)).prod

/-- Exact actual regional profile average equals the paired coefficient DP. -/
theorem regionalPairedMatrix22_expect_eq_dp
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ) (n : ℕ)
    (S : Finset (Fin ((n * 2) * 16))) :
    (fun i j => (regionalStepProfileLaw22 (n * 2) S).expect
      (fun counts => regionalPairedMatrix22 B α n counts i j)) =
      pairedRegionDP22 B α n S.card := by
  rw [pairedRegionDP22_eq_profile_sum]
  have hn : (n * 2) * 16 = 32 * n := by omega
  ext i j
  simp only [FinPMF.expect, regionalStepProfileLaw22,
    regionalStepProfileWeight22, regionalPairedMatrix22]
  rw [← Equiv.sum_comp (pairedCountEquiv22 n).symm]
  simp only [pairedCountEquiv22_degree, pairedCountEquiv22_choose_product,
    Equiv.apply_symm_apply]
  simp only [orderedMatrixCoefficient, Matrix.smul_apply, Matrix.sum_apply, smul_eq_mul,
    Finset.mul_sum]
  apply sum_congr rfl
  intro xs _
  have hc : Nat.choose ((n * 2) * 16) S.card = Nat.choose (32 * n) S.card :=
    congrArg (fun len => Nat.choose len S.card) hn
  rw [hc]
  split_ifs with hq
  · have he := matrix_ofFn_smul_prod
      (fun b : Fin n => (Nat.choose 16 (xs b).1.val : ℝ) * (Nat.choose 16 (xs b).2.val : ℝ))
      (fun b : Fin n => matrixEntrywiseRpow (B (xs b).1 * B (xs b).2) α)
    have hv := congrFun (congrFun he i) j
    change _ = (1 / (Nat.choose (32 * n) S.card : ℝ)) *
      (List.ofFn (fun b : Fin n => pairedProfileMatrix22 B α (xs b))).prod i j
    rw [show (List.ofFn (fun b : Fin n => pairedProfileMatrix22 B α (xs b))).prod =
      (∏ b : Fin n, (Nat.choose 16 (xs b).1.val : ℝ) *
        (Nat.choose 16 (xs b).2.val : ℝ)) •
        (List.ofFn (fun b : Fin n => matrixEntrywiseRpow (B (xs b).1 * B (xs b).2) α)).prod
      from he]
    simp only [Matrix.smul_apply, smul_eq_mul]
    ring
  · simp

private theorem matrix_sum_ofFn_prod {A I : Type*} [Fintype A]
    [Fintype I] [DecidableEq I] {n : ℕ}
    (K : Fin n → A → Matrix I I ℝ) :
    (∑ xs : Fin n → A, (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => ∑ a, K i a)).prod := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (n + 1) => A)),
      Fintype.sum_prod_type]
    change (∑ a : A, ∑ xs : Fin n → A,
      (List.ofFn (fun i : Fin (n + 1) => K i ((Fin.cons a xs : Fin (n + 1) → A) i))).prod) = _
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    simp_rw [← Finset.mul_sum, ih]
    rw [← Finset.sum_mul]

private def pairedMatrixExpect22 {A : Type*} [Fintype A]
    (P : FinPMF A) (f : A → Matrix (Fin 11) (Fin 11) ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ := fun i j => P.expect (fun a => f a i j)

private theorem pairedMatrixExpect22_pi_product {A : Type*} [Fintype A] {n : ℕ}
    (P : Fin n → FinPMF A) (K : Fin n → A → Matrix (Fin 11) (Fin 11) ℝ) :
    pairedMatrixExpect22 (Spin.piPMF P)
      (fun xs => (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => pairedMatrixExpect22 (P i) (K i))).prod := by
  have hleft : pairedMatrixExpect22 (Spin.piPMF P)
      (fun xs => (List.ofFn (fun i => K i (xs i))).prod) =
      ∑ xs : Fin n → A, (∏ i, (P i).p (xs i)) •
        (List.ofFn (fun i => K i (xs i))).prod := by
    ext i j
    simp only [pairedMatrixExpect22, FinPMF.expect, Spin.piPMF_apply,
      Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]
  rw [hleft]
  simp_rw [← matrix_ofFn_smul_prod]
  rw [matrix_sum_ofFn_prod (fun i a => (P i).p a • K i a)]
  apply congrArg List.prod
  apply congrArg List.ofFn
  funext k
  ext i j
  simp only [Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul,
    pairedMatrixExpect22, FinPMF.expect]

/-- Independent regional counts are averaged outside every fractional power.
The final full matrix power retains the state across all 64 regions. -/
theorem turboPairedProfiles22_expect_eq_dp_power
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ) (n : ℕ)
    (S : Finset (Fin ((n * 2) * 16))) :
    (fun i j => (turboPotentialCountLaw22 (n * 2) S).expect
      (fun counts => (List.ofFn (fun region : Fin 64 =>
        regionalPairedMatrix22 B α n (counts region))).prod i j)) =
      (pairedRegionDP22 B α n S.card) ^ 64 := by
  change pairedMatrixExpect22
    (Spin.piPMF (fun _ : Fin 64 => regionalStepProfileLaw22 (n * 2) S))
      (fun counts => (List.ofFn (fun region : Fin 64 =>
        regionalPairedMatrix22 B α n (counts region))).prod) = _
  rw [pairedMatrixExpect22_pi_product
    (fun _ : Fin 64 => regionalStepProfileLaw22 (n * 2) S)
    (fun (_ : Fin 64) counts => regionalPairedMatrix22 B α n counts)]
  have he : pairedMatrixExpect22 (regionalStepProfileLaw22 (n * 2) S)
      (regionalPairedMatrix22 B α n) = pairedRegionDP22 B α n S.card :=
    regionalPairedMatrix22_expect_eq_dp B α n S
  simp only [he]
  rw [List.ofFn_const, List.prod_replicate]
  rfl

#print axioms turboPairedProfiles22_expect_eq_dp_power
#print axioms regionalPairedMatrix22_expect_eq_dp
end Spin.Turbo




