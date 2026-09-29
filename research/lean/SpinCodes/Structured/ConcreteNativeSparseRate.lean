import SpinCodes.Structured.ConcreteNativeSparse
import SpinCodes.Structured.SparseRateOuter

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteRoute

theorem rounds_length (m : ℕ) : 128 * rounds m = Lsched m * bsched m := by
  simpa only [rounds, Nat.mul_comm] using Nat.div_mul_cancel (native_round_divisibility m)

theorem threshold_le (m : ℕ) :
    (threshold m : ℝ) ≤ (11/100) * (Lsched m : ℝ) * bsched m := by
  rw [threshold_floor]
  have h := Nat.floor_le (show (0 : ℝ) ≤ (0.11 : ℝ) * Nsched m by positivity)
  convert h using 1
  rw [Nsched_eq, Nat.cast_mul]
  ring

/-- The actual binary native message sum has the sparse rate once the realized row envelope holds. -/
theorem qd_sum_sparse_rate (m Q : ℕ) (seed : NativeSeed m)
    (hQ : 4096 ≤ Q) (hQL : Q ≤ Lsched m) (hb : 1000 ≤ bsched m)
    {B ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ≤ bsched m, (spectrum seed w : ℝ) ≤ B * ((bsched m).choose w : ℝ))
    (hrow : (2 : ℝ)^(bsched m)*B ≤
      Real.exp ((bsched m : ℝ)*(Real.log 2/2+1281/100000+ε)))
    (hsparse : (Q : ℝ)/(Lsched m : ℝ) ≤ 1/10000) (hε : ε ≤ 1/500) :
    (∑ x ∈ (nonzeroMsgs (Fin (Nsched m / 2) → ZMod 2)).filter (fun x => occ m seed x = Q),
      (nativeSetup m).qd (threshold m) ((nativeSetup m).outer seed x)) ≤
      Real.exp (-(3/500)*(Q : ℝ)*bsched m) := by
  rw [qd_sum_eq_tuple]
  have hL : (0 : ℝ) < Lsched m := by exact_mod_cast Lsched_pos m
  have hQ0 : (0 : ℝ) < Q := by exact_mod_cast (show 0 < Q by omega)
  have ha : (Lsched m : ℝ) * ((Q : ℝ)/Lsched m) = Q := by field_simp
  have hb' : 1000 ≤ nativeBlocks m * 24 := by simpa only [native_width] using hb
  have hB' : ∀ w ≤ nativeBlocks m * 24,
      (spectrum seed w : ℝ) ≤ B * ((nativeBlocks m * 24).choose w : ℝ) := by
    simpa only [native_width] using hB
  have hr := SparseRate.outer_occupation_sparse_rate seed (tupleWiring m) hQ hQL hb' hB0 hB'
    (by simpa only [native_width] using hrow) (div_pos hQ0 hL) hsparse ha
    (by simpa only [native_width] using rounds_length m)
    (by simpa only [native_width] using threshold_le m)
    (by simpa only [native_width] using SparseRate.native_log_schedule m) hε
  simpa only [native_width] using hr

/-- Selection is over the shared outer seed; its conditioning costs no further factor here. -/
theorem EZ_sparse_rate (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m Q : ℕ)
    (hQ : 4096 ≤ Q) (hQL : Q ≤ Lsched m) (hb : 1000 ≤ bsched m)
    {B ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ seed, G m seed → ∀ w ≤ bsched m,
      (spectrum seed w : ℝ) ≤ B * ((bsched m).choose w : ℝ))
    (hrow : (2 : ℝ)^(bsched m)*B ≤
      Real.exp ((bsched m : ℝ)*(Real.log 2/2+1281/100000+ε)))
    (hsparse : (Q : ℝ)/(Lsched m : ℝ) ≤ 1/10000) (hε : ε ≤ 1/500) :
    (family G hG).EZ m Q ≤ Real.exp (-(3/500)*(Q : ℝ)*bsched m) := by
  change (((nativeSeedLaw m).condition (G m) (hG m)).prod (nativeSetup m).Pin).expect
    (fun ω => ((nativeSetup m).ZQ (threshold m) (occ m) Q ω : ℝ)) ≤ _
  rw [Setup.expect_ZQ_eq_outerLaw]
  calc
    _ ≤ ((nativeSeedLaw m).condition (G m) (hG m)).expect
        (fun _ => Real.exp (-(3/500)*(Q : ℝ)*bsched m)) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      exact qd_sum_sparse_rate m Q seed hQ hQL hb hB0 (hB seed hg) hrow hsparse hε
    _ = _ := FinPMF.expect_const _ _

end Spin.Structured.ConcreteNativeFamily
