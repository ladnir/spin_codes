import SpinCodes.Turbo.CodeProperties22

/-!
# Exact rate of every realized 22-state Turbo code

As in `ConcreteNativeLinearCode`, the code is the range of the actual
injective linear encoder. RS polynomial messages have the same dimension
as the systematic information coordinates (`rsRow_four_eval_equiv`).
-/
noncomputable section
namespace Spin.Turbo

open Spin.UniformMixer

theorem turboMessage22_finrank (R : ℕ) :
    Module.finrank (ZMod 2) (TurboMessage22 R) = 2048 * R := by
  have hrow : Module.finrank Field16 (RSRow Field16) = 4 := by
    rw [(Polynomial.degreeLTEquiv Field16 4).finrank_eq]
    exact Module.finrank_fin_fun Field16
  have hbinary : Module.finrank (ZMod 2) (RSRow Field16) = 16 := by
    rw [← Module.finrank_mul_finrank (ZMod 2) Field16 (RSRow Field16), hrow]
    rw [GaloisField.finrank 2 (n := 4) (by decide)]
  simp only [TurboMessage22, RS8Group, Module.finrank_pi_fintype,
    hbinary, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    smul_eq_mul]
  omega

def turboRealizedCode22 (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    Submodule (ZMod 2) (Fin (64 * R) → Input64) :=
  (turboEncoderLinearMap22 R points seed).range

theorem turboRealizedCode22_finrank (R : ℕ) (points : Fin 8 → Field16)
    (hp : Function.Injective points)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    Module.finrank (ZMod 2) (turboRealizedCode22 R points seed) = 2048 * R := by
  rw [turboRealizedCode22,
    LinearMap.finrank_range_of_inj (turboCodeword22_injective R points hp seed)]
  exact turboMessage22_finrank R

theorem turboOutput22_finrank (R : ℕ) :
    Module.finrank (ZMod 2) (Fin (64 * R) → Input64) = 4096 * R := by
  simp only [Input64, Module.finrank_pi_fintype, Module.finrank_self,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, smul_eq_mul, mul_one]
  omega

theorem turboRealizedCode22_rate (R : ℕ) (hR : 0 < R)
    (points : Fin 8 → Field16) (hp : Function.Injective points)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    (Module.finrank (ZMod 2) (turboRealizedCode22 R points seed) : ℝ) /
      (Module.finrank (ZMod 2) (Fin (64 * R) → Input64) : ℝ) = 1 / 2 := by
  rw [turboRealizedCode22_finrank R points hp, turboOutput22_finrank]
  have hr : (0 : ℝ) < R := by exact_mod_cast hR
  push_cast
  field_simp
  <;> ring

#print axioms turboRealizedCode22_finrank
#print axioms turboRealizedCode22_rate
end Spin.Turbo
