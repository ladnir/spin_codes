import SpinCodes.Turbo.SparseExactStep22
import SpinCodes.Turbo.SparseRegionKernel22
import SpinCodes.Turbo.SparseMatrixAverage22
import SpinCodes.Turbo.MatrixCoefficient22

/-! Exact sparse regional kernels use the existing hypergeometric profile
weights and actual exact-nonzero physical kernels. The coefficient DP is
unnormalized until its one regional binomial denominator is applied. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

private theorem sparse_matrix_ofFn_smul_prod22 {I : Type*} [Fintype I] [DecidableEq I]
    {n : ℕ} (c : Fin n → ℝ) (K : Fin n → Matrix I I ℝ) :
    (List.ofFn (fun i => c i • K i)).prod = (∏ i, c i) • (List.ofFn K).prod := by
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [List.ofFn_succ, List.prod_cons, Fin.prod_univ_succ]
    rw [ih, smul_mul_smul_comm]

def sparsePhysicalRegionMatrix22 {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) (T : Finset (Fin (R * 16))) : Matrix I I ℝ :=
  SparseMatrixAverage22.average (regionalStepProfileLaw22 R T)
    (fun counts => (List.ofFn (fun step => B (counts step))).prod)

/-- Transport only the finite summation population; never normalize a
large enumerator merely to compare hidden Fintype proof fields. -/
private theorem sparseExpectFintypeTransport22 {A : Type*} {ft1 ft2 : Fintype A}
    (P1 : @FinPMF A ft1) (P2 : @FinPMF A ft2)
    (hp : @FinPMF.p A ft1 P1 = @FinPMF.p A ft2 P2) (f : A → ℝ) :
    @FinPMF.expect A ft1 P1 f = @FinPMF.expect A ft2 P2 f := by
  unfold FinPMF.expect
  have hu : @Finset.univ A ft1 = @Finset.univ A ft2 :=
    congrArg (fun ft : Fintype A => @Finset.univ A ft) (Subsingleton.elim ft1 ft2)
  rw [hu, hp]
private def sparseLawEnumeration22 {A : Type*} {fi : Fintype A}
    (_P : @FinPMF A fi) : Fintype A := fi

private def sparsePhysicalProductLaw22 {R : ℕ} (counts : Fin R → Fin 17) :=
  @Spin.piPMF (Fin R) (inferInstance) (inferInstance) (fun _ => PacketStep16)
    (fun _ => sparseLawEnumeration22 (packetStepLaw22 0))
    (fun step => packetStepLaw22 (counts step))
/-- Infer the population instance from the fixed physical law itself. -/
theorem sparsePhysicalKernel22_matrix (j : Fin 17) (z : ℝ) :
    @matrixExpect22 PacketStep16 (sparseLawEnumeration22 (packetStepLaw22 j)) (packetStepLaw22 j)
      (fun word => Matrix.of (stepKernel z (packetStepInput22 word))) =
      Matrix.of (turboPhysicalKernel22 j z) := by
  rfl

theorem sparsePhysicalStepProduct22 {R : ℕ} (counts : Fin R → Fin 17) (z : ℝ) :
    @matrixExpect22 (Fin R → PacketStep16) (sparseLawEnumeration22 (sparsePhysicalProductLaw22 counts)) (sparsePhysicalProductLaw22 counts)
      (fun words => Matrix.of (inputListKernel z (turboRegionStepStream22 words))) =
      (List.ofFn (fun step => Matrix.of (turboPhysicalKernel22 (counts step) z))).prod := by
  simp_rw [inputListKernel_turbo_steps22]
  calc
    _ = (List.ofFn (fun step : Fin R => @matrixExpect22 PacketStep16 (sparseLawEnumeration22 (packetStepLaw22 (counts step))) (packetStepLaw22 (counts step))
        (fun word => Matrix.of (stepKernel z (packetStepInput22 word))))).prod :=
      @matrixExpect22_pi_product PacketStep16 (sparseLawEnumeration22 (packetStepLaw22 0)) R
        (fun step : Fin R => packetStepLaw22 (counts step))
        (fun _ word => Matrix.of (stepKernel z (packetStepInput22 word)))
    _ = _ := by simp only [sparsePhysicalKernel22_matrix]
theorem sparsePhysicalStepProduct22_public {R : ℕ} (counts : Fin R → Fin 17) (z : ℝ) :
    matrixExpect22 (Spin.piPMF (fun step : Fin R => packetStepLaw22 (counts step)))
      (fun words => Matrix.of (inputListKernel z (turboRegionStepStream22 words))) =
      (List.ofFn (fun step => Matrix.of (turboPhysicalKernel22 (counts step) z))).prod := by
  exact sparsePhysicalStepProduct22 counts z
theorem sparseExactRegionKernel22_eq_physical_profile (R : ℕ)
    (T : Finset (Fin (R * 16))) (z : ℝ) :
    sparseExactRegionKernel22 R T z =
      sparsePhysicalRegionMatrix22 (fun j => Matrix.of (turboPhysicalKernel22 j z)) R T := by
  unfold sparseExactRegionKernel22
  rw [← matrixExpect22_map_turbo (sparseExactRegionLaw22 R T) (regionalStepWordEquiv22 R)
    (fun words => Matrix.of (inputListKernel z (turboRegionStepStream22 words)))]
  rw [sparseExactRegionLaw22_step_mixture, matrixExpect22_bind_turbo]
  simp_rw [sparsePhysicalStepProduct22_public]
  rfl

theorem sparsePhysicalRegionMatrix22_eq_coefficient {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) (T : Finset (Fin (R * 16))) :
    sparsePhysicalRegionMatrix22 B R T =
      (1 / (Nat.choose (R * 16) T.card : ℝ)) •
        orderedMatrixCoefficient Fin.val
          (fun j : Fin 17 => (Nat.choose 16 j.val : ℝ) • B j) R T.card := by
  ext i j
  simp only [sparsePhysicalRegionMatrix22, SparseMatrixAverage22.average,
    FinPMF.expect, regionalStepProfileLaw22, regionalStepProfileWeight22,
    orderedMatrixCoefficient, Matrix.smul_apply, Matrix.sum_apply, smul_eq_mul,
    Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro counts _
  split_ifs with hq
  · have h := sparse_matrix_ofFn_smul_prod22
      (fun step : Fin R => (Nat.choose 16 (counts step).val : ℝ))
      (fun step : Fin R => B (counts step))
    have he := congrFun (congrFun h i) j
    simp only [Matrix.smul_apply, smul_eq_mul] at he
    rw [he]
    ring
  · simp

theorem sparsePhysicalRegionMatrix22_eq_dp {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) (T : Finset (Fin (R * 16))) :
    sparsePhysicalRegionMatrix22 B R T =
      (1 / (Nat.choose (R * 16) T.card : ℝ)) •
        orderedMatrixDP Fin.val (fun j : Fin 17 => (Nat.choose 16 j.val : ℝ) • B j) R T.card := by
  rw [sparsePhysicalRegionMatrix22_eq_coefficient, orderedMatrixCoefficient_eq_dp]

#print axioms sparsePhysicalKernel22_matrix
#print axioms sparsePhysicalStepProduct22
#print axioms sparseExactRegionKernel22_eq_physical_profile
#print axioms sparsePhysicalRegionMatrix22_eq_dp
end Spin.Turbo











