import SpinCodes.Structured.ConcretePairing
import SpinCodes.Structured.FiberBounds

/-! Fourier inversion and Parseval for the actual selected feedback map.
These identities identify concrete kernel/fiber counts before substituting
the independently checked spectrum of the transpose. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open scoped symmDiff

def kernelLayer (j : ℕ) : Finset (Finset (Fin 128)) :=
  (Spin.layer 128 j).filter (fun x => Cset x = ∅)

def syndromeFiber (j : ℕ) (q : Finset (Fin 19)) : Finset (Finset (Fin 128)) :=
  (Spin.layer 128 j).filter (fun x => Cset x = q)

def equalSyndromePairs (j : ℕ) : Finset (Finset (Fin 128) × Finset (Fin 128)) :=
  ((Spin.layer 128 j) ×ˢ (Spin.layer 128 j)).filter (fun p => Cset p.1 = Cset p.2)

theorem sum_dual {M : Type*} [AddCommMonoid M] (f : Finset (Fin 128) → M) :
    ∑ v ∈ dual, f v = ∑ q : Finset (Fin 19), f (CtransposeSet q) := by
  rw [dual, Finset.sum_image]
  intro q _ p _ h
  exact CtransposeSet_injective h

theorem kernelLayer_eq_orth (j : ℕ) :
    kernelLayer j = (Spin.layer 128 j).filter (fun x => Spin.OrthTo dual x) := by
  ext x
  simp only [kernelLayer, Finset.mem_filter, orthTo_iff_kernel]

theorem syndromeFiber_eq (j : ℕ) (x : Finset (Fin 128)) :
    syndromeFiber j (Cset x) = Spin.fiber dual j x := (actual_fiber j x).symm

theorem equalSyndromePairs_eq_orth (j : ℕ) :
    equalSyndromePairs j = ((Spin.layer 128 j) ×ˢ (Spin.layer 128 j)).filter
      (fun p => Spin.OrthTo dual (p.1 ∆ p.2)) := by
  ext p
  simp only [equalSyndromePairs, Finset.mem_filter, orthTo_iff_kernel,
    Cset_xor, Finset.symmDiff_eq_empty]

theorem kernel_fourier (j : ℕ) :
    (2 ^ 19 : ℤ) * ((kernelLayer j).card : ℤ) =
      ∑ q : Finset (Fin 19), Spin.krawtchouk 128 j (CtransposeSet q).card := by
  have h := Spin.card_orth_layer dual (fun _ hu _ hv => dual_xor_closed hu hv) j
  rw [dual_card, ← kernelLayer_eq_orth, sum_dual] at h
  simpa only [Nat.cast_pow, Nat.cast_ofNat] using h

theorem fiber_fourier (j : ℕ) (x : Finset (Fin 128)) :
    (2 ^ 19 : ℤ) * ((syndromeFiber j (Cset x)).card : ℤ) =
      ∑ q : Finset (Fin 19), Spin.chi q (Cset x) *
        Spin.krawtchouk 128 j (CtransposeSet q).card := by
  have h := Spin.card_fiber dual (fun _ hu _ hv => dual_xor_closed hu hv) j x
  rw [dual_card, ← syndromeFiber_eq, sum_dual] at h
  simpa only [Nat.cast_pow, Nat.cast_ofNat, C_character_pairing] using h

theorem feedback_parseval (j : ℕ) :
    (2 ^ 19 : ℤ) * ((equalSyndromePairs j).card : ℤ) =
      ∑ q : Finset (Fin 19), (Spin.krawtchouk 128 j (CtransposeSet q).card) ^ 2 := by
  have h := Spin.card_orth_pairs dual (fun _ hu _ hv => dual_xor_closed hu hv) j
  rw [dual_card, ← equalSyndromePairs_eq_orth, sum_dual] at h
  simpa only [Nat.cast_pow, Nat.cast_ofNat] using h

theorem nonzero_fiber_le_complement (j : ℕ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (syndromeFiber j q).card + (kernelLayer j).card ≤ Nat.choose 128 j := by
  obtain ⟨x, rfl⟩ := Cset_surjective q
  have h := Spin.card_fiber_add_card_kernel_le dual j
    (x := x) (by simpa only [orthTo_iff_kernel] using hq)
  rwa [← syndromeFiber_eq, ← kernelLayer_eq_orth] at h

theorem fiber_le_fourier_abs (j : ℕ) (q : Finset (Fin 19)) :
    (2 ^ 19 : ℤ) * ((syndromeFiber j q).card : ℤ) ≤
      ∑ p : Finset (Fin 19), |Spin.krawtchouk 128 j (CtransposeSet p).card| := by
  obtain ⟨x, rfl⟩ := Cset_surjective q
  have h := Spin.card_fiber_le_abs_sum dual (fun _ hu _ hv => dual_xor_closed hu hv) j x
  rw [dual_card, ← syndromeFiber_eq, sum_dual] at h
  simpa only [Nat.cast_pow, Nat.cast_ofNat] using h

end Spin.Structured.ConcreteMaps
