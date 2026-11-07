import SpinCodes.UniformMixer.ScalableSparseCounted
import SpinCodes.UniformMixer.ScalableFixedSelected
import SpinCodes.Structured.ConcreteNativeSparseSelected
import SpinCodes.Structured.ConcreteNativeAsymptotic
import SpinCodes.Structured.ConcreteOuterDefectEnvelope

/-! Uniform growing-sparse first moment for the sampled Scalable setup. -/

noncomputable section
set_option maxHeartbeats 1000000
attribute [local instance] Classical.propDecidable
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute
  Spin.Structured.SparseRate Finset Filter

private theorem nativeRows_symm_eq_sparse (m : ℕ)
    (seed : NativeSeed m) (x : NativeMessage m) :
    nativeRows m seed ((nativeMessageEquiv m).symm x) =
      fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding := by
  funext i
  simp only [nativeRows, Equiv.apply_symm_apply]

theorem scalable_native_sparse_tuple_bound (m Q : ℕ)
    (seed : NativeSeed m) (hg : nativeGood m seed)
    (hQ : 4096 ≤ Q) (hQL : Q ≤ Lsched m)
    (hb : 1000 ≤ bsched m)
    {B ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ∈ weightWindow (bsched m),
      Spin.Abar (nativeSeedLaw m) spectrum w ≤
        B * ((bsched m).choose w : ℝ))
    (hrow : (2 : ℝ) ^ bsched m *
      ((bsched m : ℝ) ^ 2 * B) ≤
        Real.exp ((bsched m : ℝ) *
          (Real.log 2 / 2 + 1281 / 100000 + ε)))
    (hsparse : (Q : ℝ) / (Lsched m : ℝ) ≤ 1 / 10000)
    (hε : ε ≤ 1 / 500) :
    (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m (threshold m))) ≤
      Real.exp (-(3 / 500) * (Q : ℝ) * bsched m) := by
  let α : ℝ := (Q : ℝ) / (Lsched m : ℝ)
  have hL : (0 : ℝ) < Lsched m := by
    exact_mod_cast Lsched_pos m
  have hQ0 : (0 : ℝ) < Q := by
    exact_mod_cast (by omega : 0 < Q)
  have hα0 : 0 < α := div_pos hQ0 hL
  have hα1 : α ≤ 1 / 10000 := hsparse
  have hα : (Lsched m : ℝ) * α = Q := by
    dsimp [α]
    field_simp
  have hsteps : 64 * (bsched m * scalableStepsPerRegion m) =
      Lsched m * bsched m := by
    rw [scalable_region_length]
    ring
  have hspectrum := native_good_spectrum_bound m seed hg hB0 hB
  have hraw := scalable_sparse_tuple_raw m Q (threshold m)
    seed hspectrum hα0 hα1
  have hrowQ :
      (((2 : ℝ) ^ bsched m *
        ((bsched m : ℝ) ^ 2 * B)) ^ Q) ≤
      Real.exp ((Q : ℝ) * (bsched m : ℝ) *
        (Real.log 2 / 2 + 1281 / 100000 + ε)) := by
    have hh := pow_le_pow_left₀
      (by positivity : 0 ≤ (2 : ℝ) ^ bsched m *
        ((bsched m : ℝ) ^ 2 * B)) hrow Q
    rw [← Real.exp_nat_mul] at hh
    convert hh using 1 <;> ring
  have hcount := scalable_sparse_counted_mass hQ hQL hb
    hα0 hα1 hα hsteps (threshold_le m)
    (native_log_schedule m) hε
    (by positivity : 0 ≤ (((2 : ℝ) ^ bsched m *
      ((bsched m : ℝ) ^ 2 * B)) ^ Q)) hrowQ
  exact hraw.trans (by
    simpa only [sparseContraction, sparseTilt] using hcount)

theorem scalableSelectedLayer_sparse_bound (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m))
    (hQ : 4096 ≤ Q) (hQL : Q ≤ Lsched m)
    (hb : 1000 ≤ bsched m)
    {B ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ∈ weightWindow (bsched m),
      Spin.Abar (nativeSeedLaw m) spectrum w ≤
        B * ((bsched m).choose w : ℝ))
    (hrow : (2 : ℝ) ^ bsched m *
      ((bsched m : ℝ) ^ 2 * B) ≤
        Real.exp ((bsched m : ℝ) *
          (Real.log 2 / 2 + 1281 / 100000 + ε)))
    (hsparse : (Q : ℝ) / (Lsched m : ℝ) ≤ 1 / 10000)
    (hε : ε ≤ 1 / 500) :
    scalableSelectedLayer m Q hpos ≤
      Real.exp (-(3 / 500) * (Q : ℝ) * bsched m) := by
  rw [scalableSelectedLayer_eq_route]
  calc
    _ ≤ ((nativeSeedLaw m).condition (nativeGood m) hpos).expect
        (fun _ => Real.exp (-(3 / 500) * (Q : ℝ) * bsched m)) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      simpa only [nativeRows_symm_eq_sparse] using
        scalable_native_sparse_tuple_bound m Q seed hg hQ hQL hb
          hB0 hB hrow hsparse hε
    _ = _ := FinPMF.expect_const _ _

theorem scalableGoodLayer_sparse_eventually :
    ∀ᶠ m in atTop, ∀ Q ∈ Ico 4096 (nativeCut m + 1),
      scalableGoodLayer m Q ≤
        Real.exp (-(0.006 : ℝ) * (Q : ℝ) * bsched m) := by
  filter_upwards [native_good_positive_eventually,
    native_thresholds
      (selectedRemainder_tendsto nativeShellRemainder_tendsto)]
      with m hp hth
  intro Q hQ
  have hQlo := (mem_Ico.mp hQ).1
  have hQcut : Q ≤ nativeCut m := by
    have := (mem_Ico.mp hQ).2
    omega
  have hQL : Q ≤ Lsched m :=
    le_trans hQcut (Nat.div_le_self _ _)
  have hh := scalableSelectedLayer_sparse_bound m Q hp
    hQlo hQL hth.1
    (Real.exp_pos _).le
    (native_expected_spectrum_uniform m)
    (selected_row_cost m (nativeShellRemainder m)).le
    (nativeCut_density m Q hQcut) hth.2
  simp only [scalableGoodLayer, dite_eq_left hp]
  simpa only [show (3 / 500 : ℝ) = 0.006 by norm_num] using hh

#print axioms scalable_native_sparse_tuple_bound
#print axioms scalableGoodLayer_sparse_eventually
end Spin.UniformMixer
